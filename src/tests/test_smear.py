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
    'smear': {'axis': 1, 'direction': 1, 'readout_time': 4.33, 'n_shifts': 2048},
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


def header(camshutt=None, date_obs='2026-08-01T03:00:00', telescope='SPECULOOS-EUROPA', exptime=20.0):
    h = fits.Header()
    h['TELESCOP'] = telescope
    h['DATE-OBS'] = date_obs
    h['EXPTIME'] = exptime
    if camshutt is not None:
        h['CAMSHUTT'] = camshutt
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
