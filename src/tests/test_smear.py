"""
Unit tests for the readout-smear correction (calibration/smear.py) and the
on-frame measurement in calibration/smear_check.py.

Covers in particular:
  - the de-smear inverting the forward model exactly, along either axis and
    in either direction, and agreeing with the running-subtraction form
  - a single star's trail being one-sided and tau * F per pixel
  - NaN pixels neither spreading along the row nor being filled in
  - the shutter decision: override window beats CAMSHUTT beats the default,
    with CAMSHUTT accepted as a bool, 0/1 or a string because it is set by hand
  - correct_frame leaving shuttered frames and non-ANDOR instruments alone
  - smear_check recovering an injected t_read, and ~100% removal afterwards
  - the header record of clipped rows reproducing saturated_trail_mask exactly,
    and apertures_on_trail agreeing with a brute-force overlap test
  - SPlightcurve keeping stars on a smear trail out of the comparison ensemble
  - t_read from an override window, else Astra's READOUTM, else T6's READMODE; a shutterless
    frame in a mode with no known t_read left alone, flagged and warned about
  - clipped rows judged at the frame's own SATURATE (plus pedestal), falling
    back to the config's saturation_threshold
  - an unreadable CAMSHUTT taken as shutter used, never an exception
  - T6 (pipered.reduce_file) and T4 (pipeflat.reducer) end to end: a shuttered
    frame calibrated exactly as before, shutterless ones by readout mode
  - FLAGS keeping BW_FLAG, which the promotion check reads

Run:  pytest src/tests/test_smear.py -q
"""
import os
import sys
from datetime import datetime

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import smear                       # noqa: E402
from calibration import smear_check                 # noqa: E402

ANDOR = {
    'saturation_threshold': 62000,
    'smear': {'axis': 1, 'direction': 1, 'readout_time': {'NORMAL': 4.33, 'FAST': 1.5},
              'readout_rate_time': {'1MHz': 4.33, '3MHz': 1.5}, 'n_shifts': 2048},
}
TAU_20S = smear.smear_tau(4.33, 2048, 20.0)        # ~1.06e-4, Ganymede 20250408


def star_field(shape=(400, 600), n=40, seed=1, sky=800.0):
    """Gaussian stars on a flat sky, no noise."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    img = np.full(shape, sky)
    for _ in range(n):
        y0, x0 = rng.uniform(30, shape[0] - 30), rng.uniform(30, shape[1] - 30)
        amp = rng.uniform(2e3, 4e4)
        img += amp * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2 * 2.4 ** 2))
    return img


@pytest.mark.parametrize('axis', [0, 1])
@pytest.mark.parametrize('direction', [1, -1])
def test_desmear_inverts_forward_model(axis, direction):
    true = star_field()
    smeared = smear.add_smear(true, 0.01, axis=axis, direction=direction)
    assert not np.allclose(smeared, true)
    recovered = smear.desmear(smeared, 0.01, axis=axis, direction=direction)
    np.testing.assert_allclose(recovered, true, rtol=1e-10, atol=1e-8)


def test_iir_and_running_subtraction_agree():
    img = smear.add_smear(star_field(), TAU_20S)
    a = smear.desmear(img, TAU_20S)
    b = smear.desmear_loop(img, TAU_20S)
    assert np.max(np.abs(a - b)) < 1e-9


def test_single_star_trail_is_one_sided():
    true = np.zeros((5, 50))
    true[2, 10] = 1e5
    smeared = smear.add_smear(true, 1e-3)
    assert np.all(smeared[2, :10] == 0)                         # nothing upstream
    np.testing.assert_allclose(smeared[2, 11:], 100.0)          # tau * F downstream
    assert np.all(smeared[[0, 1, 3, 4]] == 0)                   # other rows untouched
    np.testing.assert_allclose(smear.desmear(smeared, 1e-3), true, atol=1e-9)


def test_nan_stays_put_and_does_not_spread():
    img = smear.add_smear(star_field(), TAU_20S)
    img[100, 200] = np.nan
    out = smear.desmear(img, TAU_20S)
    assert np.isnan(out[100, 200])
    assert np.isfinite(np.delete(out.ravel(), 100 * img.shape[1] + 200)).all()


def test_masked_input_keeps_its_mask():
    img = np.ma.array(smear.add_smear(star_field(), TAU_20S), mask=False)
    img.mask[5, 5] = True
    out = smear.desmear(img, TAU_20S)
    assert np.ma.isMaskedArray(out) and out.mask[5, 5] and out.mask.sum() == 1


def test_saturated_trail_mask_marks_downstream_only():
    raw = np.full((3, 10), 1000.0)
    raw[1, 4] = 65000
    mask = smear.saturated_trail_mask(raw, 62000)
    assert not mask[1, :5].any()
    assert mask[1, 5:].all()
    assert not mask[[0, 2]].any()
    flipped = smear.saturated_trail_mask(raw, 62000, direction=-1)
    assert flipped[1, :4].all() and not flipped[1, 4:].any()


def header(camshutt=None, date_obs='2026-08-01T03:00:00', telescope='SPECULOOS-EUROPA', exptime=20.0,
           readmode='NORMAL'):
    """A frame header as correct_frame sees it in T6: READMODE is already written."""
    h = fits.Header()
    h['TELESCOP'] = telescope
    h['DATE-OBS'] = date_obs
    h['EXPTIME'] = exptime
    if camshutt is not None:
        h['CAMSHUTT'] = camshutt
    if readmode is not None:
        h['READMODE'] = readmode
    return h


@pytest.mark.parametrize('value, used', [(True, True), (False, False), (0, False), (1, True),
                                         ('F', False), ('true', True)])
def test_camshutt_values(value, used):
    assert smear.shutter_state(header(value), overrides=[]) == (used, 'CAMSHUTT', None)


def test_no_header_means_shutter_used():
    assert smear.shutter_state(header(), overrides=[]) == (True, 'default', None)


def test_override_window_beats_header():
    rows = [{'telescope': 'SPECULOOS-EUROPA', 'start': datetime(2026, 7, 31, 12),
             'end': datetime(2026, 8, 1, 12), 'camshutt': False, 'readout_time': None, 'note': ''}]
    assert smear.shutter_state(header(True), rows) == (False, 'override', None)
    # outside the window, or another telescope: back to the header
    assert smear.shutter_state(header(True, date_obs='2026-08-01T12:00:00'), rows) == (True, 'CAMSHUTT', None)
    assert smear.shutter_state(header(True, telescope='SPECULOOS-IO'), rows) == (True, 'CAMSHUTT', None)


def test_override_file_parses(tmp_path):
    path = tmp_path / 'overrides.csv'
    path.write_text('# comment\ntelescope,start_utc,end_utc,camshutt,readout_time,note\n'
                    'SPECULOOS-GANYMEDE,2025-04-08T12:00:00,2025-04-09T12:00:00,F,,checked\n'
                    'SPECULOOS-GANYMEDE,2026-02-27T16:00:00,2026-03-02T16:00:00,F,1.5,fast readout\n')
    rows = smear.load_shutter_overrides(str(path))
    assert rows == [{'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2025, 4, 8, 12),
                     'end': datetime(2025, 4, 9, 12), 'camshutt': False, 'readout_time': None, 'note': 'checked'},
                    {'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2026, 2, 27, 16),
                     'end': datetime(2026, 3, 2, 16), 'camshutt': False, 'readout_time': 1.5,
                     'note': 'fast readout'}]


def test_shipped_override_file_parses():
    for row in smear.load_shutter_overrides():
        assert row['start'] < row['end']


def test_shipped_override_windows_do_not_overlap():
    rows = sorted(smear.load_shutter_overrides(), key=lambda r: (r['telescope'], r['start']))
    for a, b in zip(rows, rows[1:]):
        assert a['telescope'] != b['telescope'] or a['end'] <= b['start']


def test_override_readout_time_is_used():
    rows = [{'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2026, 2, 27, 16),
             'end': datetime(2026, 3, 2, 16), 'camshutt': False, 'readout_time': 1.5, 'note': ''}]
    true = star_field()
    h = header(None, date_obs='2026-02-28T03:00:00', telescope='SPECULOOS-GANYMEDE', exptime=41.0)
    tau = smear.smear_tau(1.5, 2048, 41.0)
    out = smear.correct_frame(smear.add_smear(true, tau), h, ANDOR, raw=true, overrides=rows)
    np.testing.assert_allclose(out, true, rtol=1e-10)
    assert h['SMEARTRD'] == 1.5 and h['SMEARSRC'] == 'override'


def test_correct_frame_shuttered_is_untouched():
    img = star_field()
    h = header(True)
    out = smear.correct_frame(img, h, ANDOR, overrides=[])
    assert out is img
    assert h['SMEARCOR'] is False


def test_correct_frame_shutterless_is_corrected():
    true = star_field()
    h = header(False)
    out = smear.correct_frame(smear.add_smear(true, TAU_20S), h, ANDOR, raw=true, overrides=[])
    np.testing.assert_allclose(out, true, rtol=1e-10)
    assert h['SMEARCOR'] is True
    assert h['SMEARTAU'] == pytest.approx(TAU_20S)
    assert h['SMEARSRC'] == 'CAMSHUTT'
    assert h['SMSATROW'] == 0


def test_correct_frame_ignores_instruments_without_smear_parameters():
    img = star_field()
    h = header(False)
    assert smear.correct_frame(img, h, {'saturation_threshold': 11000}, overrides=[]) is img
    assert 'SMEARCOR' not in h


def test_smear_check_recovers_injected_readout_time():
    rng = np.random.default_rng(7)
    shape = (1000, 1400)
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    true = np.full(shape, 600.0)
    # a sparse grid of bright isolated stars, plus fainter ones that must not matter
    for y0 in range(60, 960, 70):
        for x0 in range(200, 1200, 330):
            true += 2.5e4 * np.exp(-((yy - y0 - rng.uniform(-3, 3)) ** 2
                                     + (xx - x0 - rng.uniform(-3, 3)) ** 2) / (2 * 2.5 ** 2))
    exptime, t_read, n = 20.0, 4.33, 2048
    tau = smear.smear_tau(t_read, n, exptime)
    observed = smear.add_smear(true, tau) + rng.normal(0, 25, shape)

    stars = smear_check.find_stars(observed, observed, 62000)
    assert len(stars) >= 20
    tau_hat, err, _ = smear_check.measure(observed, stars)
    assert tau_hat * n * exptime == pytest.approx(t_read, rel=0.05)
    assert tau_hat > 5 * err

    after, _, _ = smear_check.measure(smear.desmear(observed, tau), stars)
    assert abs(after / tau_hat) < 0.03

    unsmeared = true + rng.normal(0, 25, shape)
    tau0, err0, _ = smear_check.measure(unsmeared, smear_check.find_stars(unsmeared, unsmeared, 62000))
    assert abs(tau0) < 3 * err0


@pytest.mark.parametrize('t_read', [0.0, 1.5, 4.33])
def test_flat_ramp_recovers_readout_time(t_read):
    rng = np.random.default_rng(3)
    yy, xx = np.mgrid[:600, :2046]
    vignetting = 1 - 0.02 * (((xx - 1023) / 1023) ** 2 + ((yy - 300) / 300) ** 2)
    flat = 20000.0 * vignetting
    exptime = 7.75
    observed = smear.add_smear(flat, smear.smear_tau(t_read, 2048, exptime)) + rng.normal(0, 30, flat.shape)
    assert smear_check.ramp_t_read(observed, 2048, exptime) == pytest.approx(t_read, abs=0.1)
    expected = 'SHUTTERLESS' if t_read else 'shuttered'
    assert smear_check.flat_verdict(smear_check.ramp_t_read(observed, 2048, exptime), exptime) == expected


def test_flat_verdict_on_real_values():
    # Ganymede, 2026-02 and -03: shuttered, 1 MHz shutterless, 3 MHz shutterless,
    # and a 120 s shuttered flat whose twilight gradient mimics t_read 2.4 s
    assert smear_check.flat_verdict(0.08, 7.6) == 'shuttered'
    assert smear_check.flat_verdict(4.30, 6.0) == 'SHUTTERLESS'
    assert smear_check.flat_verdict(1.50, 6.0) == 'SHUTTERLESS'
    assert smear_check.flat_verdict(2.43, 120.0) == 'inconclusive'


def clipped_frame(shape=(60, 80), n=12, seed=5):
    """Sky with a few clipped pixels, some sharing a row, one in the last column."""
    rng = np.random.default_rng(seed)
    raw = np.full(shape, 1000.0)
    raw[rng.integers(0, shape[0], n), rng.integers(0, shape[1], n)] = 65535
    raw[7, 20:24] = 65535                      # a saturated core spanning a row
    raw[8, shape[1] - 1] = 65535               # nothing downstream of this one
    raw[9, 0] = 65535                          # nor, for direction -1, of this one
    return raw


def mask_from_lines(shape, lines, first, axis, direction):
    mask = np.zeros(shape, dtype=bool)
    for line, pos in zip(lines, first):
        idx = np.arange(shape[axis])
        trail = idx > pos if direction > 0 else idx < pos
        if axis == 1:
            mask[line, trail] = True
        else:
            mask[trail, line] = True
    return mask


@pytest.mark.parametrize('axis', [0, 1])
@pytest.mark.parametrize('direction', [1, -1])
def test_clipped_lines_reproduce_trail_mask(axis, direction):
    raw = clipped_frame()
    lines, first = smear.clipped_lines(raw, 62000, axis=axis, direction=direction)
    expected = smear.saturated_trail_mask(raw, 62000, axis=axis, direction=direction)
    assert expected.any()
    np.testing.assert_array_equal(mask_from_lines(raw.shape, lines, first, axis, direction), expected)


def test_trail_header_round_trip():
    lines = np.concatenate([[3, 4, 5], [9], np.arange(100, 400)])             # a 300-row run spans cards
    first = np.concatenate([[812, 810, 811], [0], np.arange(1700, 2000)])
    h = fits.Header()
    h['SMEARCOR'] = True
    smear.write_trail_header(h, lines, first, axis=1, direction=1)
    assert h['SMSATROW'] == len(lines)
    cards = [k for k in h if k.startswith('SMTR') and k[4:].isdigit()]
    assert 10 < len(cards) < 40 and all(len(str(h[k])) <= 68 for k in cards)
    assert h['SMTR0001'] == '4:813,811,812' and h['SMTR0002'] == '10:1'     # stored 1-based
    got_lines, got_first, axis, direction = smear.read_trail_header(h)
    np.testing.assert_array_equal(got_lines, lines)
    np.testing.assert_array_equal(got_first, first)
    assert (axis, direction) == (1, 1)

    # rewriting replaces every earlier card; the header survives a FITS round trip
    smear.write_trail_header(h, np.array([50]), np.array([60]), axis=1, direction=1)
    h2 = fits.Header.fromstring(h.tostring())
    got_lines, got_first, _, _ = smear.read_trail_header(h2)
    assert list(got_lines) == [50] and list(got_first) == [60]
    assert [k for k in h2 if k.startswith('SMTR') and k[4:].isdigit()] == ['SMTR0001']


def test_read_trail_header_states():
    shuttered = header(True)
    shuttered['SMEARCOR'] = False
    lines, first, _, _ = smear.read_trail_header(shuttered)
    assert len(lines) == 0 and len(first) == 0
    assert len(smear.read_trail_header(header())[0]) == 0           # no correction recorded at all
    corrected_blind = header(False)
    corrected_blind['SMEARCOR'] = True                              # de-smeared, rows not recorded
    assert smear.read_trail_header(corrected_blind) is None
    radii = [2.0, 4.0]
    assert list(smear.apertures_on_trail(corrected_blind, [10.0, 20.0], [5.0, 5.0], radii)) == [-1, -1]
    assert list(smear.apertures_on_trail(shuttered, [10.0, 20.0], [5.0, 5.0], radii)) == [0, 0]


def test_trail_distance_geometry():
    lines, first = np.array([10]), np.array([50])                   # trail is row 10, x >= 51
    d = smear.trail_distance([40.0, 60.0, 40.0, 60.0, 50.0], [10.0, 10.0, 13.0, 13.0, 10.0], lines, first)
    np.testing.assert_allclose(d, [10.5, 0.0, np.hypot(10.5, 2.5), 2.5, 0.5])
    # direction -1: the trail is x <= 49
    d = smear.trail_distance([60.0, 40.0], [10.0, 10.0], lines, first, direction=-1)
    np.testing.assert_allclose(d, [10.5, 0.0])
    # axis 0: lines are columns, the trail runs down column 10 from y = 51
    d = smear.trail_distance([10.0, 13.0], [40.0, 60.0], lines, first, axis=0)
    np.testing.assert_allclose(d, [10.5, 2.5])
    assert np.isinf(smear.trail_distance([5.0], [5.0], np.zeros(0, int), np.zeros(0, int))).all()


@pytest.mark.parametrize('axis', [0, 1])
@pytest.mark.parametrize('direction', [1, -1])
def test_apertures_on_trail_match_brute_force(axis, direction):
    raw = clipped_frame()
    mask = smear.saturated_trail_mask(raw, 62000, axis=axis, direction=direction)
    h = fits.Header()
    h['SMEARCOR'] = True
    smear.write_trail_header(h, *smear.clipped_lines(raw, 62000, axis=axis, direction=direction), axis, direction)

    rng = np.random.default_rng(11)
    x = rng.uniform(0, raw.shape[1] - 1, 400)
    y = rng.uniform(0, raw.shape[0] - 1, 400)
    radii = 2.0 * np.array([0.5, 1 / np.sqrt(2), 1, np.sqrt(2), 2, 2 * np.sqrt(2), 4, 5, 6, 7, 8, 10, 12])
    got = smear.apertures_on_trail(h, x, y, radii)

    # an aperture reaches a pixel when the pixel's square comes closer to the centre than r
    ty, tx = np.nonzero(mask)
    for k in range(len(x)):
        dx = np.maximum(np.abs(tx - x[k]) - 0.5, 0)
        dy = np.maximum(np.abs(ty - y[k]) - 0.5, 0)
        nearest = np.min(np.hypot(dx, dy))
        reaching = np.flatnonzero(radii > nearest)
        assert got[k] == (reaching[0] + 1 if len(reaching) else 0), (k, x[k], y[k], nearest)
    assert len(np.unique(got)) >= 5                                  # many apertures exercised


def test_correct_frame_records_trails_for_the_photometry():
    true = star_field(n=0)
    yy, xx = np.mgrid[:true.shape[0], :true.shape[1]]
    true += 4e5 * np.exp(-((yy - 200) ** 2 + (xx - 150) ** 2) / (2 * 2.4 ** 2))   # saturates
    raw = np.minimum(smear.add_smear(true, TAU_20S), 65535)
    h = header(False)
    smear.correct_frame(raw, h, ANDOR, raw=raw, overrides=[])
    assert h['SMSATROW'] > 0 and 'SMTR0001' in h
    radii = [2.0, 4.0, 8.0]
    # downstream in the same rows, upstream in the same rows, downstream but well clear of the rows
    flags = smear.apertures_on_trail(fits.Header.fromstring(h.tostring()),
                                     [400.0, 60.0, 400.0], [200.0, 200.0, 230.0], radii)
    assert list(flags) == [1, 0, 0]


def test_splightcurve_keeps_trail_stars_out_of_the_ensemble():
    pytest.importorskip('fitsio')
    pytest.importorskip('matplotlib')
    from diff_photometry import SPlightcurve as spl

    smear_ap = np.zeros((4, 6), dtype=np.int16)
    smear_ap[1, 2] = 4          # star 1 reaches a trail at aperture 4 and above, in one frame
    smear_ap[2, :] = 7          # star 2 only from aperture 7
    smear_ap[3, :] = -1         # star 3 not assessed
    on = spl.smear_trail_frames(smear_ap, 5, 4, 6)
    assert on[1].tolist() == [False, False, True, False, False, False]
    assert not on[[0, 2, 3]].any()
    assert not spl.smear_trail_frames(None, 5, 4, 6).any()

    class Tel:
        rcore = np.full(6, 4.0)
        exposure = np.full(6, 20.0)
        gain = np.full(6, 1.0)

    jd = 2460774.5 + np.arange(6) / 1440
    flux = [np.full(6, 5e4) for _ in range(4)]
    nflux = [np.ones(6) for _ in range(4)]
    saturated = np.zeros(4, dtype=bool)          # saturated_stars() found none
    nan4 = [np.nan] * 4
    lcurves = spl.init_all_lcurves(Tel(), ['SP0', 'SP1', 'SP2', 'SP3'], nan4, nan4, nflux, flux, saturated, jd, jd, jd,
                                   nan4, nan4, nan4, nan4, nan4, [np.zeros(6)] * 4, [np.zeros(6)] * 4, nan4, nan4,
                                   nan4, nan4, 5, 'I+z', '20250408', np.full(6, 20.0), smear_ap)
    assert [lc.smeared for lc in lcurves] == [False, True, False, False]
    assert [lc.weight for lc in lcurves] == [1, 0, 1, 1]
    assert lcurves[1].smear.tolist() == on[1].tolist()
    lcurves[1].w3 = 0                    # and init_weights, as the iteration's reset calls it, keeps it out
    lcurves[1].init_weights()
    assert lcurves[1].weight == 0 and lcurves[1].w1 == 0


def test_best_aperture_avoids_smear_trails():
    pytest.importorskip('fitsio')
    pytest.importorskip('matplotlib')
    from diff_photometry import SPlightcurve as spl

    def choose(runs):
        """The pipeline's j = 3..8 loop: each run compares itself with what the file holds."""
        stored = None
        for ap, metric, on_trail in runs:
            if stored is None or spl.replaces_best_ap(stored, metric, on_trail):
                stored = '%d %r %d' % (ap, metric, on_trail)
        return stored.split()

    # TOI-1227 on Ganymede 20250408: apertures 6-8 reach the G 8.8 star's trail
    assert choose([(3, 5e-6, 0), (4, 4e-6, 0), (5, 6e-6, 0), (6, 1e-6, 1), (7, 9e-7, 1), (8, 8e-7, 1)])[0] == '4'
    # a clear aperture displaces a better one on a trail
    assert choose([(3, 1e-6, 1), (4, 9e-6, 0)])[0] == '4'
    # all on trails: the best metric still wins, and the file says so
    assert choose([(3, 3e-6, 1), (4, 2e-6, 1), (5, 4e-6, 1)]) == ['4', '2e-06', '1']
    # files written before the flag have two fields and count as clear
    assert spl.replaces_best_ap('6 4.2e-06', 1e-6, False)
    assert not spl.replaces_best_ap('6 4.2e-06', 1e-6, True)


# ------------------------------------------- readout mode, saturation, CAMSHUTT

@pytest.fixture
def fresh_warnings():
    smear._warned.clear()
    yield
    smear._warned.clear()


@pytest.mark.parametrize('mode, t_read', [('NORMAL', 4.33), ('FAST', 1.5)])
def test_readout_time_follows_readmode(mode, t_read):
    true = star_field()
    h = header(False, readmode=mode)
    out = smear.correct_frame(smear.add_smear(true, smear.smear_tau(t_read, 2048, 20.0)), h, ANDOR, raw=true,
                              overrides=[])
    np.testing.assert_allclose(out, true, rtol=1e-10)
    assert h['SMEARCOR'] is True and h['SMEARTRD'] == t_read
    assert ('READMODE ' + mode) in h.comments['SMEARTRD']


def test_override_readout_time_beats_readmode():
    rows = [{'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2026, 2, 27, 16),
             'end': datetime(2026, 3, 2, 16), 'camshutt': False, 'readout_time': 1.5, 'note': ''}]
    h = header(None, date_obs='2026-02-28T03:00:00', telescope='SPECULOOS-GANYMEDE', readmode='NORMAL')
    smear.correct_frame(star_field(), h, ANDOR, overrides=rows)
    assert h['SMEARTRD'] == 1.5 and '(override)' in h.comments['SMEARTRD']


def test_override_without_readout_time_uses_readmode():
    rows = [{'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2026, 2, 15, 12),
             'end': datetime(2026, 2, 27, 16), 'camshutt': False, 'readout_time': None, 'note': ''}]
    h = header(None, date_obs='2026-02-21T03:00:00', telescope='SPECULOOS-GANYMEDE', readmode='FAST')
    smear.correct_frame(star_field(), h, ANDOR, overrides=rows)
    assert h['SMEARTRD'] == 1.5 and h['SMEARSRC'] == 'override'


@pytest.mark.parametrize('readmode', ['PED7500', 'GAIN3.8', 'PED2900', 'UNKNOWN', 'DEFAULT', '', None])
def test_unknown_readout_mode_is_left_alone_and_flagged(readmode, capsys, fresh_warnings):
    img = star_field()
    h = header(False, readmode=readmode)
    out = smear.correct_frame(img, h, ANDOR, raw=img, overrides=[])
    assert out is img
    assert h['SMEARCOR'] is False and h['SMEARSRC'] == 'CAMSHUTT'
    assert h['SMEARERR'].startswith('no readout time for READMODE')
    assert not any(k in h for k in ('SMEARTRD', 'SMEARTAU', 'SMTRAXIS', 'SMSATROW'))
    assert 'readout smear NOT removed' in capsys.readouterr().out
    # once per process, not once per frame
    smear.correct_frame(img, header(False, readmode=readmode), ANDOR, overrides=[])
    assert capsys.readouterr().out == ''
    # and the photometry treats it as it did before the correction existed: no trail record
    assert len(smear.read_trail_header(h)[0]) == 0


@pytest.mark.parametrize('readoutm, t_read', [(2, 4.33), (1, 1.5)])
def test_readoutm_gives_the_readout_time_and_beats_readmode(readoutm, t_read):
    true = star_field()
    # the overscan class says otherwise: the camera's own card wins
    h = header(False, readmode='FAST' if readoutm == 2 else 'NORMAL')
    h['READOUTM'] = readoutm
    out = smear.correct_frame(smear.add_smear(true, smear.smear_tau(t_read, 2048, 20.0)), h, ANDOR, raw=true,
                              overrides=[])
    np.testing.assert_allclose(out, true, rtol=1e-10)
    assert h['SMEARTRD'] == t_read and 'READOUTM %d (' % readoutm in h.comments['SMEARTRD']


@pytest.mark.parametrize('readoutm', [0, 3, 4, -1, 'x', 2.5, True])
def test_readoutm_without_a_measured_time_is_left_alone(readoutm, capsys, fresh_warnings):
    img = star_field()
    h = header(False, readmode='NORMAL')            # no fallback to the overscan class
    h['READOUTM'] = readoutm
    assert smear.correct_frame(img, h, ANDOR, raw=img, overrides=[]) is img
    assert h['SMEARCOR'] is False and 'READOUTM' in h['SMEARERR']
    assert 'readout smear NOT removed' in capsys.readouterr().out


def test_override_readout_time_beats_readoutm():
    rows = [{'telescope': 'SPECULOOS-GANYMEDE', 'start': datetime(2026, 10, 18, 12),
             'end': datetime(2026, 10, 20, 12), 'camshutt': False, 'readout_time': 2.0, 'note': ''}]
    h = header(None, date_obs='2026-10-19T03:00:00', telescope='SPECULOOS-GANYMEDE')
    h['READOUTM'] = 2
    smear.correct_frame(star_field(), h, ANDOR, overrides=rows)
    assert h['SMEARTRD'] == 2.0 and '(override)' in h.comments['SMEARTRD']


def test_readoutm_rates():
    assert [smear.readoutm_rate(i) for i in range(4)] == ['5MHz', '3MHz', '1MHz', '0.05MHz']
    assert smear.readoutm_rate('2') == '1MHz' and smear.readoutm_rate(2.0) == '1MHz'
    assert smear.readoutm_rate(None) is None and smear.readoutm_rate(False) is None


def test_shipped_config_gives_the_measured_readout_times():
    path = os.path.join(SRC, 'calibration', 'instrument_config.json')
    import json
    with open(path) as f:
        config = json.load(f)['telescopes']
    for tel in ('SPECULOOS-IO', 'SPECULOOS-EUROPA', 'SPECULOOS-GANYMEDE', 'SPECULOOS-CALLISTO'):
        cfg = config[tel]['andor']['smear']
        assert smear.mode_readout_time(cfg, 'NORMAL') == 4.33
        assert smear.mode_readout_time(cfg, 'FAST') == 1.5
        assert smear.mode_readout_time(cfg, 'PED7500') is None
        assert cfg['readout_rate_time'] == {'1MHz': 4.33, '3MHz': 1.5}     # 5 and 0.05 MHz never measured
    gany = config['SPECULOOS-GANYMEDE']['andor']
    rows = smear.load_shutter_overrides()
    # 1 MHz window: READMODE decides; 3 MHz window: the override's own 1.5 s
    for date_obs, mode, expected in (('2026-02-21T03:00:00', 'NORMAL', 4.33), ('2026-02-28T03:00:00', 'FAST', 1.5),
                                     ('2026-02-28T03:00:00', 'NORMAL', 1.5), ('2025-04-08T23:00:00', 'NORMAL', 4.33)):
        h = header(None, date_obs=date_obs, telescope='SPECULOOS-GANYMEDE', readmode=mode)
        smear.correct_frame(star_field(shape=(60, 80), n=2), h, gany, overrides=rows)
        assert h['SMEARCOR'] is True and h['SMEARTRD'] == expected, date_obs


def _saturating_star(pedestal, peak):
    true = star_field(n=0, sky=500.0)
    yy, xx = np.mgrid[:true.shape[0], :true.shape[1]]
    true += peak * np.exp(-((yy - 200) ** 2 + (xx - 150) ** 2) / (2 * 2.4 ** 2))
    raw = smear.add_smear(true, TAU_20S) + pedestal
    return raw - pedestal, raw


def test_clipped_rows_use_the_frames_own_saturation():
    # Ganymede 2026-03-03 (7,546 pedestal mode): saturates ~36,400 ADU above the pedestal
    pedestal, saturate = 7546.0 + 7.0, 36400.0
    signal, raw = _saturating_star(pedestal, 45000.0)
    assert raw.max() < 62000                       # the config level would never call it clipped
    h = header(False)
    h['SATURATE'] = saturate
    smear.correct_frame(signal, h, ANDOR, raw=raw, overrides=[], pedestal=pedestal)
    assert h['SMSATLEV'] == pytest.approx(saturate + pedestal) and 'SATURATE' in h.comments['SMSATLEV']
    assert h['SMSATROW'] > 0 and 'SMTR0001' in h
    lines, first, _, _ = smear.read_trail_header(h)
    want_lines, want_first = smear.clipped_lines(raw, saturate + pedestal)
    assert list(lines) == list(want_lines) and list(first) == list(want_first)


def test_clipped_rows_fall_back_to_the_config_level():
    signal, raw = _saturating_star(300.0, 45000.0)
    h = header(False)                               # calibrated before SATURATE existed
    smear.correct_frame(signal, h, ANDOR, raw=raw, overrides=[], pedestal=300.0)
    assert h['SMSATLEV'] == 62000 and 'config' in h.comments['SMSATLEV']
    assert h['SMSATROW'] == 0
    signal, raw = _saturating_star(300.0, 70000.0)
    h = header(False)
    h['SATURATE'] = 'N'                             # the pipeline's 'absent' value: also the config level
    smear.correct_frame(signal, h, ANDOR, raw=raw, overrides=[], pedestal=300.0)
    assert h['SMSATLEV'] == 62000 and h['SMSATROW'] > 0


@pytest.mark.parametrize('value', ['open', 'UNKNOWN', '', 'shutterless', 'T?'])
def test_unreadable_camshutt_is_taken_as_shuttered(value, capsys, fresh_warnings):
    h = header(value)
    assert smear.shutter_state(h, overrides=[]) == (True, 'CAMSHUTT?', None)
    assert 'unrecognised CAMSHUTT' in capsys.readouterr().out
    img = star_field()
    assert smear.correct_frame(img, h, ANDOR, overrides=[]) is img
    assert h['SMEARCOR'] is False and 'CAMSHUTT' in h.comments['SMEARCOR']
    assert not smear.shutterless(h, ANDOR, overrides=[])


def test_unreadable_override_camshutt_is_taken_as_shuttered(tmp_path, capsys, fresh_warnings):
    path = tmp_path / 'overrides.csv'
    path.write_text('telescope,start_utc,end_utc,camshutt,readout_time,note\n'
                    'SPECULOOS-GANYMEDE,2025-04-08T12:00:00,2025-04-09T12:00:00,maybe,,typo\n')
    rows = smear.load_shutter_overrides(str(path))
    assert rows[0]['camshutt'] is True
    assert 'unrecognised camshutt' in capsys.readouterr().out


def test_shipped_override_file_has_only_readable_camshutt_values():
    import csv
    with open(smear.OVERRIDES_PATH) as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith('#')]
    for row in csv.DictReader(lines):
        assert smear._parse_camshutt(row['camshutt']) is not None, row


def test_shutterless_needs_smear_parameters():
    assert smear.shutterless(header(False), ANDOR, overrides=[])
    assert not smear.shutterless(header(True), ANDOR, overrides=[])
    assert not smear.shutterless(header(False), {'saturation_threshold': 11000}, overrides=[])


# --------------------------------------------------------------- T6 and T4 end to end

def _andor_raw(rng, pedestal, noise, signal=None):
    """A full-size SSO ANDOR raw frame: 2088 x 2048, overscan in the last 15 rows."""
    raw = np.full((2088, 2048), pedestal)
    raw[:-15] += 7.0                                 # bias structure above the overscan
    if signal is not None:
        raw[22:2066, 2:2048] += signal
    raw += rng.normal(0, noise, raw.shape)
    return np.round(raw)


def _write_raw(path, raw, date_obs, camshutt=None, exptime=20.0, imagetyp='Light Frame',
               telescope='SPECULOOS-GANYMEDE'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    hdu = fits.PrimaryHDU(raw.astype(np.float32))
    hdu.header['TELESCOP'] = telescope
    hdu.header['EXPTIME'] = exptime
    hdu.header['FILTER'] = 'I+z'
    hdu.header['DATE-OBS'] = date_obs
    hdu.header['JD'] = 2460000.5
    hdu.header['IMAGETYP'] = imagetyp
    if camshutt is not None:
        hdu.header['CAMSHUTT'] = camshutt
    hdu.writeto(path, overwrite=True)
    return path


def _reduce(tmp_path, path, night, telescope='SPECULOOS-GANYMEDE'):
    from calibration import pipered
    from calibration import readout as ro
    trimmed = (2066 - 22, 2048 - 2)
    bias = np.full(trimmed, 7.0, dtype=np.float32)
    flat = {'I+z': np.ones(trimmed, dtype=np.float32)}
    dark = {'combined': np.zeros(trimmed, dtype=np.float32)}
    with fits.open(path) as hdul:
        params = pipered.get_instrument_parameters(hdul)
    record = ro.measure_night([], [], night, telescope, params)
    outdir = str(tmp_path) + '/out/'
    os.makedirs(outdir, exist_ok=True)
    pipered.reduce_file(path, outdir, params, 'bias.fits', bias, dark, flat, None, '6.2', 0.0, 'v3',
                        readout_record=record)
    out = outdir + 'proc' + os.path.basename(path).replace('.fts', '.fits')
    with fits.open(out) as h:
        return h[0].header.copy(), np.array(h[0].data)


def test_t6_shuttered_frame_is_calibrated_as_before(tmp_path):
    rng = np.random.default_rng(3)
    signal = smear.add_smear(star_field(shape=(2044, 2046), n=30, sky=600.0), 0.0)
    raw = _andor_raw(rng, 320.0, 6.2, signal)
    path = _write_raw(os.path.join(str(tmp_path), 'images', '20260810', 'a.fts'), raw, '2026-08-11T03:00:00',
                      camshutt=True)
    h, data = _reduce(tmp_path, path, '20260810')
    # what main's T6 computes: (raw - overscan - bias - dark * t) / flat, as float32
    expected = np.float32((raw[22:2066, 2:2048] - h['OVERSCAN'] - 7.0 - 0.0) / 1.0)
    assert np.array_equal(data, expected)
    assert h['SMEARCOR'] is False
    assert [k for k in h if k.startswith('SM')] == ['SMEARCOR']
    assert h['READMODE'] == 'NORMAL'


def test_t6_desmears_by_readout_mode_and_uses_the_frames_saturation(tmp_path):
    rng = np.random.default_rng(4)
    true = star_field(shape=(2044, 2046), n=30, sky=600.0)
    tau = smear.smear_tau(4.33, 2048, 20.0)
    raw = _andor_raw(rng, 320.0, 6.2, smear.add_smear(true, tau))
    # shutterless outside every override window: CAMSHUTT and READMODE decide
    path = _write_raw(os.path.join(str(tmp_path), 'images', '20260810', 'b.fts'), raw, '2026-08-11T03:00:00',
                      camshutt=False)
    h, data = _reduce(tmp_path, path, '20260810')
    assert h['SMEARCOR'] is True and h['SMEARTRD'] == 4.33 and h['READMODE'] == 'NORMAL'
    assert h['SMSATLEV'] == pytest.approx(h['SATURATE'] + h['OVERSCAN'] + 7.0)
    resid = data - (true + 7.0 + 320.0 - h['OVERSCAN'] - 7.0)
    assert np.abs(np.median(resid)) < 1.0 and np.std(resid) < 10.0


def test_t6_fast_readout_window_and_unknown_mode(tmp_path, fresh_warnings, capsys):
    rng = np.random.default_rng(5)
    true = star_field(shape=(2044, 2046), n=30, sky=600.0)
    # 2026-02-28, 3 MHz: the override window gives 1.5 s; the overscan says FAST
    raw = _andor_raw(rng, 319.0, 9.8, smear.add_smear(true, smear.smear_tau(1.5, 2048, 20.0)))
    path = _write_raw(os.path.join(str(tmp_path), 'images', '20260227', 'c.fts'), raw, '2026-02-28T03:00:00')
    h, _ = _reduce(tmp_path, path, '20260227')
    assert h['READMODE'] == 'FAST' and h['SMEARCOR'] is True and h['SMEARTRD'] == 1.5
    assert h['SMEARSRC'] == 'override'
    assert h['SMSATLEV'] == pytest.approx(59500.0)      # FAST's level, not the config's 62,000
    # 2026-03-03, 7,546 pedestal mode, CAMSHUTT F: no readout time known, so left as it is
    raw = _andor_raw(rng, 7546.0, 18.2, true)
    path = _write_raw(os.path.join(str(tmp_path), 'images', '20260302', 'd.fts'), raw, '2026-03-03T03:00:00',
                      camshutt=False)
    h, data = _reduce(tmp_path, path, '20260302')
    assert h['READMODE'] == 'PED7500' and h['SMEARCOR'] is False and 'PED7500' in h['SMEARERR']
    assert np.array_equal(data, np.float32(raw[22:2066, 2:2048] - h['OVERSCAN'] - 7.0))
    assert 'readout smear NOT removed' in capsys.readouterr().out


def test_t4_master_flat_from_shutterless_flats(tmp_path, fresh_warnings, capsys):
    from calibration import pipeflat
    rng = np.random.default_rng(6)
    reddir = os.path.join(str(tmp_path), 'output', '20260810', 'reduction') + '/'
    os.makedirs(reddir)
    trimmed = (2044, 2046)
    fits.PrimaryHDU(np.full(trimmed, 7.0, dtype=np.float32)).writeto(reddir + 'run_MasterBias.fits')
    fits.PrimaryHDU(np.zeros(trimmed, dtype=np.float32)).writeto(reddir + 'run_MasterDark.fits')
    illum = np.ones(trimmed) * 20000.0
    names = []
    for i, (camshutt, date_obs) in enumerate([(False, '2026-08-10T23:00:%02d'), (True, '2026-08-10T23:01:%02d')]):
        for j in range(2):
            exptime = 6.0 + j
            flat = smear.add_smear(illum, smear.smear_tau(4.33, 2048, exptime)) if not camshutt else illum
            names.append(_write_raw(os.path.join(str(tmp_path), 'images', '20260810', 'f%d%d.fts' % (i, j)),
                                    _andor_raw(rng, 320.0, 6.2, flat), date_obs % j, camshutt=camshutt,
                                    exptime=exptime, imagetyp='Flat Field'))
    flist = reddir + 'run_flat_I+z.list'
    with open(flist, 'w') as f:
        f.write('\n'.join(names) + '\n')
    pipeflat.reducer([flist], 'run_MasterBias.fits', 'run_MasterDark.fits', 'run_MasterFlat.fits', '1', '1',
                     reddir, reddir)
    with fits.open(reddir + 'run_MasterFlat_I+z.fits') as h:
        master, hdr = np.array(h[0].data), h[0].header
    assert hdr['NSMEARED'] == 2 and 'NSMRLEFT' not in hdr
    # the de-smeared flats carry no readout ramp along x
    left, right = np.median(master[:, 100:300]), np.median(master[:, -300:-100])
    assert abs(right / left - 1) < 2e-3
    assert 'De-smeared shutterless flat' in capsys.readouterr().out


# ------------------------------------------------------------------ FLAGS readers

def test_flags_keep_bw_flag_for_the_promotion_check(tmp_path):
    pytest.importorskip('fitsio')
    pytest.importorskip('matplotlib')
    from diff_photometry import SPlightcurve as spl
    from utils import promote_check

    bw = np.array([0, 1, 1, 0, 0, 1])
    table = spl.flags_table(bw, np.array([0, 0, 1, 1, 0, 0], dtype=bool))
    assert list(table.dtype.names) == ['FLAG', 'BW_FLAG', 'SMEAR_FLAG']
    assert list(table['FLAG']) == list(bw) and list(table['BW_FLAG']) == list(bw)

    cat = fits.BinTableHDU.from_columns([fits.Column('TARGET', 'K', array=np.array([1, 0])),
                                         fits.Column('GAIA_DR3_ID', '19A', array=np.array(['123', '456']))],
                                        name='CATALOGUE')
    lc = fits.ImageHDU(np.vstack([1 + 1e-3 * np.arange(6), np.ones(6)]).T, name='LIGHTCURVE_5')
    path = str(tmp_path / '123_I+z_20260810_5_diff.fits')
    fits.HDUList([fits.PrimaryHDU(), cat, lc, spl.flags_hdu(table)]).writeto(path)
    assert promote_check.read_light_curve(path)['bw'] == 3


# ------------------------------------------------- T9 condense: SMEAR_AP and --update

def _condense_night(tmp_path, night, jd0, nframes, trail_line=None, readout_cards=True):
    """Fake T8 products for one night of a 3-star field, then condense_photometry.main()."""
    pytest.importorskip('fitsio')
    import argparse
    from condense import condense_photometry as cp

    night_dir = tmp_path / 'output' / night / 'TARG'
    night_dir.mkdir(parents=True, exist_ok=True)
    xs, ys = np.array([301.0, 101.0, 501.0]), np.array([301.0, 51.0, 151.0])   # imcore, 1-based
    nstars = len(xs)
    cat = tmp_path / 'stack_catalogue.fits'
    if not cat.exists():
        cols = [fits.Column(n, 'D', array=np.zeros(nstars)) for n in ('pmra', 'pmdec', 'g_rp', 'bp_rp', 'parallax',
                                                                       'teff', 'gmag')]
        cols += [fits.Column('gaia_dr2_id', '19A', array=np.array(['111', '222', '333'])),
                 fits.Column('gaia_dr3_id', '19A', array=np.array(['111', '222', '333'])),
                 fits.Column('target_role', '10A', array=np.array(['primary', '', '']))]
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols, name='Gaia_Crossmatch')]).writeto(cat)
    phots = []
    for k in range(nframes):
        jd = jd0 + k / 1440.
        frame = night_dir / ('proc%s_%03d.fits' % (night, k))
        h = fits.Header()
        h['SMEARCOR'] = trail_line is not None
        if trail_line is not None:
            smear.write_trail_header(h, np.array([trail_line]), np.array([9]), 1, 1)
        fits.PrimaryHDU(np.zeros((8, 8), dtype=np.float32), header=h).writeto(frame)
        cols = [fits.Column('dummy%02d' % i, 'E', array=np.zeros(nstars)) for i in range(1, 20)]
        for a in range(1, 14):
            cols += [fits.Column('Aper_flux_%d' % a, 'E', array=np.full(nstars, 1e4 * a)),
                     fits.Column('Aper_flux_%d_err' % a, 'E', array=np.full(nstars, 10.0))]
        cols += [fits.Column('X_coordinate', 'E', array=xs), fits.Column('Y_coordinate', 'E', array=ys),
                 fits.Column('Sky_level', 'E', array=np.full(nstars, 500.0)),
                 fits.Column('Peak_height', 'E', array=np.full(nstars, 3e3)),
                 fits.Column('ra', 'D', array=np.radians([10.0, 10.01, 10.02])),
                 fits.Column('dec', 'D', array=np.radians([-20.0, -20.01, -20.02]))]
        table = fits.BinTableHDU.from_columns(cols)
        for key, value in (('jd-obs', jd), ('bjd-obs', jd), ('hjd-obs', jd), ('exptime', 20.0), ('rcore', 4.0),
                           ('TCTYP3', 'RA---ZPN'), ('TCTYP5', 'DEC--ZPN'), ('TCRPX3', 1.0), ('TCRPX5', 1.0),
                           ('TCRVL3', 10.0), ('TCRVL5', -20.0), ('TC3_3', 1e-4), ('TC3_5', 0.0), ('TC5_3', 0.0),
                           ('TC5_5', 1e-4)):
            table.header[key] = value
        if readout_cards:
            for key, value in (('RDNOISE', 6.5), ('SATURATE', 61500.0), ('OSNOISE', 6.2), ('READMODE', 'NORMAL'),
                               ('GAINSRC', 'night')):
                table.header[key] = value
        primary = fits.PrimaryHDU()
        primary.header['jd-obs'] = jd
        phot = str(frame) + '.phot'
        fits.HDUList([primary, table]).writeto(phot)
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols[:1])]).writeto(
            str(frame).replace('.fits', '.cat'))
        phots.append(phot)
    args = argparse.Namespace(filename=phots, outputdir=str(tmp_path / 'output'), update=False, verbose=False,
                              sha=None, date=night, targname='TARG', outname='output.fits', filter='I+z',
                              ext='fits', version='v3', cat=str(cat), tlist=None)
    cp.main(args)
    with fits.open(str(night_dir / ('111_I+z_%s_output.fits' % night))) as h:
        nightly = {x.name: (x.data.copy() if x.data is not None else None) for x in h[1:]}
        last = h[-1].name
    return nightly, last, str(tmp_path / 'output' / '111_I+z_output.fits')


def test_condense_smear_ap_with_readout_columns_and_update(tmp_path):
    # night A, shuttered: SMEAR_AP all 0, the last HDU; the readout columns last in the IMAGELIST
    a, last, glob_path = _condense_night(tmp_path, '20260810', 2461263.5, 3)
    assert last == 'SMEAR_AP' and a['SMEAR_AP'].shape == (3, 3) and not a['SMEAR_AP'].any()
    assert list(a['IMAGELIST'].dtype.names[-5:]) == ['RDNOISE', 'SATURATE', 'OSNOISE', 'READMODE', 'GAINSRC']

    # make the global file what one condensed before #32 and the smear flag looks like
    with fits.open(glob_path) as h:
        old_cols = [c for c in h['IMAGELIST'].columns if c.name not in ('RDNOISE', 'SATURATE', 'OSNOISE',
                                                                          'READMODE', 'GAINSRC')]
        hdus = [h[0].copy(), fits.BinTableHDU.from_columns(old_cols, name='IMAGELIST')]
        hdus += [x.copy() for x in h[1:] if x.name not in ('IMAGELIST', 'SMEAR_AP')]
        fits.HDUList(hdus).writeto(glob_path, overwrite=True)

    # night B, shutterless with star 1 (0-based y = 50) on a trail row: updates that old global file
    b, _, _ = _condense_night(tmp_path, '20260811', 2461264.5, 2, trail_line=50)
    assert b['SMEAR_AP'][:, 0].tolist() == [0, 1, 0]
    with fits.open(glob_path) as h:
        assert h['SMEAR_AP'].data.shape == h['JD'].data.shape == (3, 5)
        assert (h['SMEAR_AP'].data[:, :3] == smear.SMEAR_NOT_ASSESSED).all()     # never assessed
        np.testing.assert_array_equal(h['SMEAR_AP'].data[:, 3:], b['SMEAR_AP'])
        il = h['IMAGELIST'].data
        assert 'READMODE' not in il.dtype.names and len(il) == 5
        np.testing.assert_allclose(il['JD-OBS'][3:], b['IMAGELIST']['JD-OBS'])   # positional zip kept the rows
        np.testing.assert_allclose(il['EXPTIME'], 20.0)

    # night C on the now current-format global file: appended, no second backfill
    c, _, _ = _condense_night(tmp_path, '20260812', 2461265.5, 2)
    with fits.open(glob_path) as h:
        assert [x.name for x in h].count('SMEAR_AP') == 1
        assert h['SMEAR_AP'].data.shape == h['JD'].data.shape == (3, 7)
        np.testing.assert_array_equal(h['SMEAR_AP'].data[:, 5:], c['SMEAR_AP'])
        np.testing.assert_array_equal(h['SMEAR_AP'].data[:, 3:5], b['SMEAR_AP'])

    # a job still pinned to an older release updates the global file: every HDU but SMEAR_AP gains its two
    # frames. The next update pads SMEAR_AP with "not assessed" for them, so it lines up with JD again.
    with fits.open(glob_path, mode='update') as h:
        for x in h:
            if isinstance(x, fits.ImageHDU) and x.name != 'SMEAR_AP':
                x.data = np.append(x.data, x.data[:, -2:] + (1.0 if x.name in ('JD', 'HJD', 'BJD') else 0.0), axis=1)
    e, _, _ = _condense_night(tmp_path, '20260814', 2461267.5, 2, trail_line=50)
    with fits.open(glob_path) as h:
        sa, jd = h['SMEAR_AP'].data, h['JD'].data
        assert sa.shape == jd.shape == (3, 11)
        assert (sa[:, 7:9] == smear.SMEAR_NOT_ASSESSED).all()
        np.testing.assert_array_equal(sa[:, 9:], e['SMEAR_AP'])
        np.testing.assert_array_equal(sa[:, 3:5], b['SMEAR_AP'])


def test_pad_smear_ap():
    pytest.importorskip('fitsio')
    from condense import condense_photometry as cp
    na = smear.SMEAR_NOT_ASSESSED
    assert (cp.pad_smear_ap(None, (2, 3)) == na).all()
    got = cp.pad_smear_ap(np.array([[1, 0], [0, 2]], dtype=np.int16), (2, 4))
    assert got.tolist() == [[1, 0, na, na], [0, 2, na, na]] and got.dtype == np.int16
    assert (cp.pad_smear_ap(np.zeros((3, 2), dtype=np.int16), (2, 4)) == na).all()      # stars differ
    same = np.array([[4, 0]], dtype=np.int16)
    assert cp.pad_smear_ap(same, (1, 2)).tolist() == [[4, 0]]


def test_stable_weights_leave_smeared_stars_out(tmp_path):
    """#28's stable weights: a star on a trail gets no weight, and SMEAR_AP all 0 changes nothing."""
    pytest.importorskip('fitsio')
    pytest.importorskip('matplotlib')
    from diff_photometry import SPlightcurve as spl

    rng = np.random.RandomState(2)
    n, nf = 30, 200
    t = np.linspace(0, 1, nf)
    transparency = 1 - 0.1 * np.exp(-((t - 0.5) / 0.05) ** 2)
    base = np.geomspace(2e5, 5e3, n)
    flux = base[:, None] * transparency[None, :] * (1 + 0.003 * rng.randn(n, nf))
    ra = 150.0 + rng.uniform(-0.05, 0.05, n)
    dec = -30.0 + rng.uniform(-0.05, 0.05, n)

    class Tel:
        rcore = np.full(nf, 4.0)
        exposure = np.full(nf, 20.0)
        gain = np.full(nf, 1.0)

    def weights(smear_ap):
        jd = 2461095.5 + np.arange(nf) / 1440.
        nan = [np.nan] * n
        lcurves = spl.init_all_lcurves(Tel(), ['SP%03d' % i for i in range(n)], nan, nan,
                                       [spl.normalise(f) for f in flux], [f.copy() for f in flux],
                                       np.zeros(n, dtype=bool), jd, jd, jd, list(ra), list(dec), nan, nan, nan,
                                       [np.zeros(nf)] * n, [np.zeros(nf)] * n, nan, nan, nan, nan, 5, 'I+z',
                                       '20260810', np.full(nf, 20.0), smear_ap)
        lcurves[0].set_as_target('SP000')
        comp, _ = spl.stable_weights(spl.set_all_weights(np.delete(lcurves, 0)), lcurves[0])
        return np.array(spl.extract(comp, 'weight'), dtype=float), comp

    w_none, _ = weights(None)
    w_zero, _ = weights(np.zeros((n, nf), dtype=np.int16))
    assert np.array_equal(w_none, w_zero)                 # a shuttered night: bit-identical weights
    smear_ap = np.zeros((n, nf), dtype=np.int16)
    smear_ap[3, 50] = 2                                   # star 3 (comparison 2) on a trail at ap 2+, one frame
    smear_ap[4, :] = 8                                    # star 4 only from ap 8: clear at ap 5
    w_trail, comp = weights(smear_ap)
    assert w_none[2] > 0 and w_trail[2] == 0 and comp[2].smeared
    assert w_trail[3] > 0 and not comp[3].smeared
    assert abs(w_trail.sum() - 1) < 1e-9
