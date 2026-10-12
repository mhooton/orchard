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
