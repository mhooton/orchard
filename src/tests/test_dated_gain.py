"""
Unit tests for dated gains ("gain_epochs" in instrument_config.json).

SPIRIT has no overscan, so readout.py cannot measure its readout per night and every stage
(T2 read noise, T3/T5 dark current, T6 GAIN) takes the configured gain.  That gain changed
between 2022 and 2023, so the config lists it by date and get_instrument_parameters() picks
the epoch from the frame's own DATE-OBS.

Run:  pytest src/tests/test_dated_gain.py -q
"""
import datetime
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import pipeutils as pu            # noqa: E402
from calibration import readout as ro              # noqa: E402

EPOCHS = {'gain': 4.8, 'naxis2': 1280,
          'gain_epochs': [{'from': '2022-01-01T12:00:00', 'gain': 5.0},
                          {'from': '2023-01-01T12:00:00', 'gain': 4.8}]}


def _hdr(**cards):
    h = fits.Header()
    for k, v in cards.items():
        h[k] = v
    return h


def test_epoch_follows_the_frame_date():
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2022-10-16T03:12:44.123'}))['gain'] == 5.0
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2023-02-14T01:00:00'}))['gain'] == 4.8
    # the boundary is at midday UTC, between nights: a frame taken before it keeps the old epoch
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2023-01-01T09:59:59'}))['gain'] == 5.0
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2023-01-01T12:00:00'}))['gain'] == 4.8


def test_jd_is_used_when_date_obs_is_missing():
    # 2022-10-16T00:00 UTC = JD 2459868.5
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'JD-OBS': 2459868.5}))['gain'] == 5.0
    assert pu.resolve_dated_gain(EPOCHS, _hdr(JD=2460300.5))['gain'] == 4.8


def test_undated_or_too_early_frames_keep_the_plain_gain():
    assert pu.resolve_dated_gain(EPOCHS, _hdr())['gain'] == 4.8
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': 'not a date'}))['gain'] == 4.8
    assert pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2021-06-01T00:00:00'}))['gain'] == 4.8


def test_parameters_without_epochs_are_returned_unchanged():
    andor = {'gain': 1.0029, 'naxis2': 2088}
    assert pu.resolve_dated_gain(andor, _hdr(**{'DATE-OBS': '2022-10-16T03:12:44'})) is andor


def test_the_config_is_never_modified():
    before = dict(EPOCHS)
    resolved = pu.resolve_dated_gain(EPOCHS, _hdr(**{'DATE-OBS': '2022-10-16T03:12:44'}))
    assert resolved is not EPOCHS and EPOCHS == before and EPOCHS['gain'] == 4.8


def test_configured_epochs_are_ordered_and_plausible():
    config = pu.load_instrument_config()
    found = 0
    for tel, instruments in config['telescopes'].items():
        for name, params in instruments.items():
            epochs = params.get('gain_epochs')
            if not epochs:
                continue
            found += 1
            starts = [datetime.datetime.strptime(e['from'], '%Y-%m-%dT%H:%M:%S') for e in epochs]
            assert starts == sorted(starts), (tel, name)
            assert all(s.hour == 12 and s.minute == 0 for s in starts), 'boundaries sit at midday UTC'
            assert all(0.5 < float(e['gain']) < 10 for e in epochs)
            assert float(params['gain']) == float(epochs[-1]['gain']), 'plain gain = latest epoch'
    assert found >= 1


def _spirit_frame(path, date_obs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = fits.PrimaryHDU(data=np.zeros((1280, 1024), dtype=np.float32))
    h.header['TELESCOP'] = 'SPECULOOS-CALLISTO'
    h.header['DATE-OBS'] = date_obs
    h.writeto(str(path))
    return str(path)


def test_spirit_frames_get_the_gain_of_their_epoch(tmp_path):
    config = pu.load_instrument_config()
    epochs = config['telescopes']['SPECULOOS-CALLISTO']['spirit']['gain_epochs']
    first, last = epochs[0], epochs[-1]
    when_first = datetime.datetime.strptime(first['from'], '%Y-%m-%dT%H:%M:%S') + datetime.timedelta(days=30)
    when_last = datetime.datetime.strptime(last['from'], '%Y-%m-%dT%H:%M:%S') + datetime.timedelta(days=30)
    a = _spirit_frame(tmp_path / 'Observations/Callisto/images/x/a.fits', when_first.strftime('%Y-%m-%dT%H:%M:%S'))
    b = _spirit_frame(tmp_path / 'Observations/Callisto/images/y/b.fits', when_last.strftime('%Y-%m-%dT%H:%M:%S'))
    with fits.open(a) as ha, fits.open(b) as hb:
        assert pu.detect_instrument(ha) == 'spirit'
        pa, pb = pu.get_instrument_parameters(ha), pu.get_instrument_parameters(hb)
    assert pa['gain'] == first['gain'] and pb['gain'] == last['gain']
    # and the readout path (T2 -> readout.json -> T6 GAIN) carries it through for SPIRIT
    record = ro.measure_night([a], [], when_first.strftime('%Y%m%d'), 'SPECULOOS-CALLISTO', pa)
    assert ro.frame_values(record, None, None, legacy_ron_e=79.0)['GAIN'] == first['gain']


def test_andor_parameters_are_untouched(tmp_path):
    path = tmp_path / 'Observations/Io/images/20221016/i.fits'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = fits.PrimaryHDU(data=np.zeros((2088, 8), dtype=np.uint16))
    h.header['TELESCOP'] = 'SPECULOOS-IO'
    h.header['DATE-OBS'] = '2022-10-16T03:12:44'
    h.writeto(str(path))
    with fits.open(str(path)) as hdul:
        params = pu.get_instrument_parameters(hdul)
    assert params is pu.load_instrument_config()['telescopes']['SPECULOOS-IO']['andor'] or \
        params == pu.load_instrument_config()['telescopes']['SPECULOOS-IO']['andor']
    assert params['gain'] == 1.0029
