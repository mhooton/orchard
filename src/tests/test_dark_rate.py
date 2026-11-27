"""
Unit tests for the dark current T5 writes to darkcurrent.dat (calibration/pipebadpixel.py).

SPIRIT's master darks are exposure-matched and stored per second, (dark - bias) / t.  Each carries an offset c
from the master bias, so its median is r + c/t and the shortest master overstates the dark current: on Callisto
20251212 the 2 s master read 39.0 ADU/s against 24.2 at 120 s, and T5 reported 4.77 x 39.0 = 186 e-/s for a
true ~114.  exposure_fit_rate() takes the slope of median x t against t instead.  Nights with one master (every
Andor night: one combined master) keep the plain median.

Run:  pytest src/tests/test_dark_rate.py -q
"""
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import pipebadpixel as pbp      # noqa: E402

SHAPE = (40, 50)
EXPOSURES = (2, 5, 10, 15, 30, 60, 120)


def _masters(rate, offset, exposures=EXPOSURES, noise=0.0, seed=1):
    """{exposure: master dark in ADU/s} for a detector with dark current `rate` and offset `offset` ADU."""
    rng = np.random.default_rng(seed)
    return {t: (rate + offset / t + noise * rng.standard_normal(SHAPE) / t).astype(np.float32) for t in exposures}


def test_the_fit_removes_the_offset_the_shortest_master_carries():
    masters = _masters(24.0, 25.0)
    rate, note = pbp.exposure_fit_rate(masters)
    assert rate == pytest.approx(24.0, rel=1e-5)
    assert np.median(masters[2]) == pytest.approx(36.5, rel=1e-5)      # what T5 used to report, in ADU/s
    assert 'slope over 7 exposures 2-120 s' in note and 'offset 25.0 ADU' in note


def test_noise_and_a_negative_offset_are_handled():
    rate, _ = pbp.exposure_fit_rate(_masters(21.9, -18.0, noise=40.0))
    assert rate == pytest.approx(21.9, rel=0.01)


def test_one_bad_master_does_not_move_the_rate():
    # Callisto 20220715's medians (ADU/s): the 68 s master reads twice the others' rate.  Least squares gave
    # 14.4 ADU/s; the median of the pairwise slopes stays with the other six masters.
    night = {5: 17.1, 10: 12.2, 11: 10.636, 15: 10.933, 30: 9.8, 60: 9.467, 68: 17.353}
    rate, _ = pbp.exposure_fit_rate({t: np.full(SHAPE, m, np.float32) for t, m in night.items()})
    assert 8.5 < rate < 9.5


def test_bad_pixels_are_left_out_of_each_median():
    masters = _masters(24.0, 25.0)
    bpm = np.zeros(SHAPE, bool)
    bpm[:, :30] = True                                     # 60% of the frame, far above the dark level
    for t in masters:
        masters[t][bpm] = 5000.0
    assert pbp.exposure_fit_rate(masters)[0] > 1000
    assert pbp.exposure_fit_rate(masters, bpm)[0] == pytest.approx(24.0, rel=1e-5)


def test_a_night_with_fewer_than_two_exposure_matched_masters_is_left_to_the_caller():
    one = np.full(SHAPE, 30.0, np.float32)
    assert pbp.exposure_fit_rate({10: one}) is None
    assert pbp.exposure_fit_rate({'combined': one}) is None
    assert pbp.exposure_fit_rate({'combined': one, 120: one}) is None   # MasterDark.fits copies the longest


def test_an_unusable_slope_falls_back_to_the_longest_master():
    # totals falling with exposure: no dark current can do that, so the longest master's median is used
    masters = {10: np.full(SHAPE, 50.0, np.float32), 20: np.full(SHAPE, 10.0, np.float32),
               40: np.full(SHAPE, 2.0, np.float32)}
    rate, note = pbp.exposure_fit_rate(masters)
    assert rate == pytest.approx(2.0) and note.startswith('longest master (40 s)')


# --- bad_pixel_maps() end to end on a tiny night -------------------------------------------------------------

def _night(tmp_path, masters, combined=None):
    cal = tmp_path / 'reduction'
    cal.mkdir()
    frame = tmp_path / 'frame.fits'
    fits.writeto(frame, np.zeros(SHAPE, np.float32), fits.Header({'EXPTIME': 10.0}))
    (cal / '1_image_T.list').write_text(str(frame) + '\n')
    fits.writeto(cal / '1_MasterBias.fits', np.zeros(SHAPE, np.float32))
    fits.writeto(cal / '1_MasterFlat_zYJ.fits', np.ones(SHAPE, np.float32))
    darknames = []
    for t, dark in masters.items():
        fits.writeto(cal / '1_MasterDark_{}s.fits'.format(t), dark)
        darknames.append(str(cal / '1_MasterDark_{}s.fits'.format(t)))
    if combined is not None:
        fits.writeto(cal / '1_MasterDark.fits', combined)
        darknames.append(str(cal / '1_MasterDark.fits'))
    return cal, darknames


def _patch(monkeypatch, params, gain):
    monkeypatch.setattr(pbp, 'get_instrument_parameters', lambda hdul, **kw: params)
    monkeypatch.setattr(pbp, 'frame_telescope', lambda hdul: 'SPECULOOS-CALLISTO')
    monkeypatch.setattr(pbp.readout, 'load_or_measure', lambda *a, **kw: {})
    monkeypatch.setattr(pbp.readout, 'calib_values', lambda values: {'GAIN': gain})


def test_bad_pixel_maps_writes_the_fitted_dark_current_for_spirit(tmp_path, monkeypatch):
    params = {'gain': 4.77, 'saturation_threshold': 11000,
              'bad_pixel_correction': {'preferred_flat_filter': 'zYJ', 'hot_sigma_threshold': 5,
                                       'cold_sigma_threshold': 5, 'flat_threshold': 0.1, 'overwrite_existing': True}}
    _patch(monkeypatch, params, 4.77)
    masters = _masters(24.0, 25.0, noise=5.0)
    cal, darknames = _night(tmp_path, masters, combined=masters[120])
    pbp.bad_pixel_maps(str(cal / '1_image_T.list'), str(cal) + '/', darknames, [str(cal / '1_MasterFlat_zYJ.fits')],
                       '1_BadPixelMap.fits', '1_MasterBias.fits')
    darkc = float((cal / 'darkcurrent.dat').read_text())
    assert darkc == pytest.approx(4.77 * 24.0, rel=0.01)                # ~114.5, not 4.77 x 36.5 = 174
    assert (cal / '1_BadPixelMap.fits').exists()


def test_a_combined_master_keeps_the_plain_median(tmp_path, monkeypatch):
    # the Andor route: no bad pixel correction, one combined master; the value must not change
    _patch(monkeypatch, {'gain': 1.0029, 'saturation_threshold': 62000}, 1.0029)
    rng = np.random.default_rng(3)
    combined = (0.02 + 0.5 * rng.standard_normal(SHAPE)).astype(np.float32)
    cal, darknames = _night(tmp_path, {}, combined=combined)
    pbp.bad_pixel_maps(str(cal / '1_image_T.list'), str(cal) + '/', darknames, [str(cal / '1_MasterFlat_zYJ.fits')],
                       '1_BadPixelMap.fits', '1_MasterBias.fits')
    assert float((cal / 'darkcurrent.dat').read_text()) == pytest.approx(float(1.0029) * np.median(combined), rel=1e-6)
