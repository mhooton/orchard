"""
Tests for per-night gain, read noise and saturation (calibration/readout.py)
and the stages that consume them.

Synthetic frames carry Poisson signal at a known gain and Gaussian read noise
at a known level, in a small detector geometry (the measurement only needs the
overscan strip and the central 512x512 pixels).  Covered:
  - the overscan statistics and the readout-mode classification, including
    the episodes found in the archive (Europa 2019, Artemis 2020/2025-26,
    Ganymede 2026-02/03, Artemis 2026-08, Artemis's noisy 2023);
  - gain and read noise recovered from bias and flat pairs, per mode, with
    saturated flats and borrowed frames left out;
  - the fallbacks: dated table, estimate, config, legacy instruments;
  - T6 writing the keywords into a calibrated frame;
  - T10's per-frame saturation test and the PDF report's readout cell.

Run:  pytest src/tests/test_readout.py -q
"""
import importlib
import json
import os
import sys
import types

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import readout as ro            # noqa: E402

# A small detector laid out like the Andors: 15 overscan rows at the bottom of
# the array, the image trimmed to rows 22..578 and columns 2..560.
PARAMS = {
    'gain': 1.0029,
    'saturation_threshold': 62000,
    'naxis2': 600,
    'overscan': {'top_row': -15, 'bottom_row': None, 'left_col': 4, 'right_col': None},
    'trim': {'top_row': 22, 'bottom_row': 578, 'left_col': 2, 'right_col': 560},
}
SHAPE = (600, 562)


def frame(rng, level_adu=0.0, gain=1.0, ron_adu=6.0, pedestal=300.0, bias_offset=7.0,
          saturate=None, prnu=None):
    """Raw frame in ADU: pedestal + Poisson(level*gain e-)/gain + N(0, ron)."""
    ny, nx = SHAPE
    img = np.full(SHAPE, pedestal, dtype=float)
    img[:-15] += bias_offset
    if level_adu > 0:
        yy, xx = np.mgrid[0:ny - 15, 0:nx]
        # vignetting-like gradient, plus pixel response non-uniformity
        illum = level_adu * (1.0 + 0.04 * (xx / nx - 0.5) + 0.03 * (yy / ny - 0.5))
        if prnu is not None:
            illum = illum * prnu[:ny - 15]
        img[:-15] += rng.poisson(illum * gain) / gain
    img += rng.normal(0.0, ron_adu, SHAPE)
    if saturate is not None:
        img = np.minimum(img, saturate)
    return np.round(img)


def write(path, data, imagetyp, t, filt='I+z', exptime=10.0, telescope='SPECULOOS-IO'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    hdu = fits.PrimaryHDU(data.astype(np.float32))
    hdu.header['IMAGETYP'] = imagetyp
    hdu.header['DATE-OBS'] = t
    hdu.header['FILTER'] = filt
    hdu.header['EXPTIME'] = exptime
    hdu.header['TELESCOP'] = telescope
    hdu.writeto(path, overwrite=True)
    return path


def make_night(tmp_path, night, rng, gain, ron_adu, flat_levels, nbias=6, pedestal=300.0,
               saturate=None, telescope='SPECULOOS-IO', prefix='', t0=10):
    """Bias and flat frames for one night; returns (bias paths, flat paths)."""
    base = os.path.join(str(tmp_path), 'images', night, 'Calibration')
    prnu = 1.0 + 0.01 * rng.standard_normal(SHAPE)
    bias, flats = [], []
    for i in range(nbias):
        bias.append(write(os.path.join(base, '%sBias-%03d.fits' % (prefix, i)),
                          frame(rng, 0.0, gain, ron_adu, pedestal),
                          'Bias Frame', '2026-01-01T%02d:00:%02d' % (t0, i), telescope=telescope))
    for i, level in enumerate(flat_levels):
        flats.append(write(os.path.join(base, '%sFlat-%03d.fits' % (prefix, i)),
                           frame(rng, level, gain, ron_adu, pedestal, saturate=saturate, prnu=prnu),
                           'Flat Field', '2026-01-01T%02d:01:%02d' % (t0, i), telescope=telescope))
    return bias, flats


# ------------------------------------------------------------ frame statistics

def test_tile_sigma_recovers_gaussian_noise():
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 5.0, (512, 512))
    noise[100:110, 200:260] += 4000.0          # a star-like defect is clipped
    assert ro.tile_sigma(noise) == pytest.approx(5.0, rel=0.03)


def test_tile_sigma_matches_clipping_tile_by_tile():
    from astropy.stats import sigma_clipped_stats
    rng = np.random.default_rng(11)
    data = rng.normal(0, 3.0, (150, 330)) + rng.standard_exponential((150, 330)) ** 4
    expected = np.median([sigma_clipped_stats(data[y:y + 64, x:x + 64], sigma=4, maxiters=3)[2]
                          for y in range(0, 150 - 63, 64) for x in range(0, 330 - 63, 64)])
    assert ro.tile_sigma(data) == pytest.approx(expected, rel=1e-9)
    strip = data[:15]
    expected = np.median([sigma_clipped_stats(strip[:, x:x + 64], sigma=4, maxiters=3)[2]
                          for x in range(0, 330 - 63, 64)])
    assert ro.tile_sigma(strip, 15, 64) == pytest.approx(expected, rel=1e-9)


def test_overscan_signature_reads_the_configured_strip():
    rng = np.random.default_rng(2)
    level, noise = ro.overscan_signature(frame(rng, 20000.0, 1.0, 2.45, 300.0), PARAMS)
    assert level == pytest.approx(300.0, abs=1.0)
    assert noise == pytest.approx(2.45, rel=0.06)


def test_frames_without_overscan_have_no_signature():
    # Artemis wrote 2048 x 2048 frames in 2021-04/06: the strip would be sky
    data = np.zeros((SHAPE[0] - 40, SHAPE[1]))
    assert ro.overscan_signature(data, PARAMS) == (None, None)
    assert ro.overscan_signature(np.zeros(SHAPE), {'gain': 5.09, 'naxis2': 600}) == (None, None)


# ------------------------------------------------------------- classification

@pytest.mark.parametrize('telescope, night, noise, level, expected', [
    ('SPECULOOS-IO', '20240616', 5.9, 297, 'NORMAL'),
    ('SPECULOOS-GANYMEDE', '20240616', 7.4, 300, 'NORMAL'),        # noisy but normal
    # Europa switched during 2019-07-25: evening frames normal, later ones not
    ('SPECULOOS-EUROPA', '20190725', 5.8, 300, 'NORMAL'),
    ('SPECULOOS-EUROPA', '20190725', 2.43, 301, 'GAIN3.8'),
    ('SPECULOOS-EUROPA', '20190815', 2.55, 300, 'GAIN3.8'),
    # Artemis was already in the 3.8 e-/ADU mode on 2025-10-29 (no biases)
    ('Artemis', '20251029', 2.42, 301, 'GAIN3.8'),
    ('Artemis', '20200806', 2.43, 300, 'GAIN3.8'),
    ('Artemis', '20200805', 6.6, 303, 'NORMAL'),
    # 2026-08-06/07: bias level ~2,900, then back to 3.8 during 08-07
    ('Artemis', '20260806', 6.25, 2890, 'PED2900'),
    ('Artemis', '20260807', 6.40, 2810, 'PED2900'),
    ('Artemis', '20260807', 2.43, 300, 'GAIN3.8'),
    # Artemis 2023: up to 10.7 ADU at the normal gain is still the normal mode
    ('Artemis', '20230815', 10.7, 301, 'NORMAL'),
    ('Artemis', '20230815', 6.5, 301, 'NORMAL'),
    ('Artemis', '20230815', 12.0, 301, 'UNKNOWN'),
    # Ganymede 2026-02-27..03-01 and 03-02/03
    ('SPECULOOS-GANYMEDE', '20260228', 10.1, 321, 'FAST'),
    ('SPECULOOS-GANYMEDE', '20260228', 11.1, 320, 'FAST'),         # bright flat
    ('SPECULOOS-GANYMEDE', '20260305', 6.7, 298, 'NORMAL'),
    ('SPECULOOS-GANYMEDE', '20260303', 18.0, 7525, 'PED7500'),
    ('SPECULOOS-GANYMEDE', '20260303', 29.0, 7354, 'PED7500'),     # bright flat
    # a mode never seen on the telescope is still recognised if distinctive
    ('SPECULOOS-IO', '20270101', 2.45, 300, 'GAIN3.8'),
    ('SPECULOOS-IO', '20270101', 12.0, 318, 'UNKNOWN'),
    # Europa's 2017-12 readout fault
    ('SPECULOOS-EUROPA', '20171210', 100.0, 300, 'UNKNOWN'),
])
def test_classify(telescope, night, noise, level, expected):
    assert ro.classify(noise, level, telescope, night) == expected


def test_classify_without_overscan_follows_the_table():
    assert ro.classify(None, None, 'Artemis', '20260420') == 'GAIN3.8'
    assert ro.classify(None, None, 'Artemis', '20250714') == 'NORMAL'
    assert ro.classify(None, None, 'SAINT-EX', '20250714') == 'NORMAL'


def test_ganymede_normal_gain_changed_in_2022_01():
    assert ro.regime_for_mode('SPECULOOS-GANYMEDE', 'NORMAL', '20220109').gain == 0.99
    assert ro.regime_for_mode('SPECULOOS-GANYMEDE', 'NORMAL', '20220111').gain == 1.11
    old = ro.frame_values(_empty('SPECULOOS-GANYMEDE', '20210815'), 300.0, 6.0)
    assert (old['GAIN'], old['GAINSRC']) == (0.99, 'table')


def test_regime_lookup():
    assert ro.expected_regime('Artemis', '20260806').mode == 'PED2900'
    assert ro.expected_regime('SPECULOOS-IO', '20260806').mode == 'NORMAL'
    assert ro.expected_regime('SPECULOOS-IO', '20090101') is None
    # two 3.8 e-/ADU spells either side of 2026-08-06/07: the nearer one
    assert ro.regime_for_mode('Artemis', 'GAIN3.8', '20260807').first == '20260807'
    assert ro.regime_for_mode('Artemis', 'GAIN3.8', '20260801').first == '20251029'
    assert ro.regime_for_mode('SPECULOOS-IO', 'GAIN3.8', '20270101') is None


def test_acp_era_artemis_frames_use_the_artemis_regimes():
    assert ro.canonical_telescope('ACP->Artemis') == 'Artemis'
    assert ro.expected_regime('ACP->Artemis', '20200806').mode == 'GAIN3.8'
    assert ro.classify(2.45, 300, 'ACP->Artemis', '20200806') == 'GAIN3.8'
    assert ro.measure_night([], [], '20200806', 'ACP->Artemis', PARAMS)['telescope'] == 'Artemis'


def test_own_night():
    assert ro.own_night('/data/Observations/Io/images/20260423/Calibration/Bias.fits', '20260423')
    assert not ro.own_night('/data/Observations/Io/images/20260422/Calibration/Bias.fits', '20260423')


# --------------------------------------------------------------- measurement

def test_measure_night_recovers_gain_and_read_noise(tmp_path):
    rng = np.random.default_rng(3)
    levels = [26000, 23500, 21300, 19300, 17500, 15900, 14400, 13000, 11800]
    bias, flats = make_night(tmp_path, '20260423', rng, gain=1.08, ron_adu=6.0, flat_levels=levels)
    record = ro.measure_night(bias, flats, '20260423', 'SPECULOOS-IO', PARAMS)
    mode = record['modes']['NORMAL']
    assert record['calib_mode'] == 'NORMAL'
    assert mode['ron_adu'] == pytest.approx(6.0, rel=0.03)
    assert mode['gain_accepted'] and mode['gain_pairs'] >= 6
    assert mode['gain'] == pytest.approx(1.08, rel=0.03)
    values = ro.calib_values(record)
    assert values['GAINSRC'] == 'night' and values['RONSRC'] == 'bias'
    assert values['RDNOISE'] == pytest.approx(6.0 * 1.08, rel=0.05)
    json.dumps(record)                       # the record must be writable


def test_saturated_flats_are_not_used(tmp_path):
    """Flats pinned at the full well understate the variance and inflate the
    gain (5.4 instead of 3.83 on Artemis 2025-12-02); they must be left out."""
    rng = np.random.default_rng(4)
    levels = [45000, 42000, 40000, 26000, 24000, 22000, 20000, 18000]
    bias, flats = make_night(tmp_path, '20260420', rng, gain=3.83, ron_adu=2.45, flat_levels=levels,
                             saturate=39600, telescope='Artemis')
    record = ro.measure_night(bias, flats, '20260420', 'Artemis', PARAMS)
    mode = record['modes']['GAIN3.8']
    assert mode['gain_accepted']
    assert mode['gain'] == pytest.approx(3.83, rel=0.04)
    assert all(p['m1'] < 33000 and p['m2'] < 33000 for p in mode['pairs'])


def test_borrowed_frames_are_not_measured(tmp_path):
    rng = np.random.default_rng(5)
    own_bias, own_flats = make_night(tmp_path, '20260423', rng, 1.04, 5.8, [24000, 21800, 19800, 18000])
    old_bias, old_flats = make_night(tmp_path, '20260301', rng, 3.83, 2.45, [24000, 21800, 19800, 18000],
                                     prefix='old')
    record = ro.measure_night(old_bias + own_bias, old_flats + own_flats, '20260423', 'SPECULOOS-IO', PARAMS)
    assert record['n_bias_own'] == len(own_bias) and record['n_flats_own'] == len(own_flats)
    assert list(record['modes']) == ['NORMAL']
    assert record['modes']['NORMAL']['ron_adu'] == pytest.approx(5.8, rel=0.04)
    assert any('borrowed' in n for n in record['notes'])


def test_mode_switch_during_a_night_is_measured_per_mode(tmp_path):
    """Europa 2019-07-25: evening flats in the normal mode, biases and later
    flats in the 3.8 e-/ADU mode."""
    rng = np.random.default_rng(6)
    night = '20190725'
    bias, flats38 = make_night(tmp_path, night, rng, 3.80, 2.45, [24000, 21800, 19800, 18000, 16400],
                               telescope='SPECULOOS-EUROPA', t0=11)
    _, flats_n = make_night(tmp_path, night, rng, 1.08, 6.0, [24000, 21800, 19800, 18000, 16400],
                            nbias=0, telescope='SPECULOOS-EUROPA', prefix='eve', t0=1)
    record = ro.measure_night(bias, flats_n + flats38, night, 'SPECULOOS-EUROPA', PARAMS)
    assert set(record['modes']) == {'NORMAL', 'GAIN3.8'}
    assert record['modes']['GAIN3.8']['gain'] == pytest.approx(3.80, rel=0.05)
    assert record['modes']['NORMAL']['gain'] == pytest.approx(1.08, rel=0.04)
    assert record['calib_mode'] == 'GAIN3.8'
    # a science frame of each mode gets its own values
    normal = ro.frame_values(record, 300.0, 6.1)
    fast = ro.frame_values(record, 300.0, 2.5)
    assert (normal['READMODE'], fast['READMODE']) == ('NORMAL', 'GAIN3.8')
    assert normal['GAIN'] == pytest.approx(1.08, rel=0.04) and normal['GAINSRC'] == 'night'
    assert normal['RONSRC'] == 'overscan'           # no normal-mode biases that night
    assert fast['RONSRC'] == 'bias'
    assert fast['SATURATE_RAW'] == 34000.0          # Europa's 3.8 e-/ADU level


def test_too_few_flat_pairs_fall_back_to_the_table(tmp_path):
    rng = np.random.default_rng(7)
    bias, flats = make_night(tmp_path, '20260420', rng, 3.83, 2.45, [24000, 22000], telescope='Artemis')
    record = ro.measure_night(bias, flats, '20260420', 'Artemis', PARAMS)
    assert not record['modes']['GAIN3.8'].get('gain_accepted')
    values = ro.frame_values(record, 300.0, 2.45)
    assert (values['GAIN'], values['GAINSRC']) == (3.83, 'table')
    assert values['RONSRC'] == 'bias'
    assert values['SATURATE_RAW'] == 37300.0


# ---------------------------------------------------------------- resolution

def _empty(telescope, night, params=PARAMS):
    return ro.measure_night([], [], night, telescope, params)


def test_night_without_calibration_frames_uses_frame_overscan_and_table():
    record = _empty('Artemis', '20251105')
    values = ro.frame_values(record, 300.0, 2.42)
    assert values['READMODE'] == 'GAIN3.8'
    assert (values['GAIN'], values['GAINSRC']) == (3.83, 'table')
    assert values['RDNOISE'] == pytest.approx(2.42 * 3.83)
    assert values['RONSRC'] == 'overscan'


def test_ped2900_gain_is_an_estimate():
    values = ro.frame_values(_empty('Artemis', '20260806'), 2890.0, 6.3)
    assert values['READMODE'] == 'PED2900'
    assert (values['GAIN'], values['GAINSRC']) == (6.5, 'estimate')
    assert values['SATURATE_RAW'] == 22800.0


def test_telescope_without_regimes_keeps_its_config_gain():
    params = dict(PARAMS, gain=3.48)
    values = ro.frame_values(ro.measure_night([], [], '20250101', 'SAINT-EX', params), 300.0, 6.2)
    assert (values['READMODE'], values['GAIN'], values['GAINSRC']) == ('NORMAL', 3.48, 'config')
    # a telescope the table knows, in a mode it has never used: the mode's gain
    io = ro.frame_values(_empty('SPECULOOS-IO', '20270101'), 300.0, 2.45)
    assert (io['READMODE'], io['GAIN'], io['GAINSRC']) == ('GAIN3.8', 3.80, 'table')


def test_unknown_readout_keeps_the_expected_modes_gain_and_saturation():
    values = ro.frame_values(_empty('SPECULOOS-EUROPA', '20171210'), 300.0, 100.0)
    assert values['READMODE'] == 'UNKNOWN'
    assert (values['GAIN'], values['GAINSRC']) == (1.08, 'table')
    assert values['RDNOISE'] == pytest.approx(108.0)
    assert values['SATURATE_RAW'] == 62000.0


def test_frame_noise_far_from_the_bias_pairs_uses_the_frames_own(tmp_path):
    """Artemis 2023-08-15: biases at ~10 ADU, some science frames at 6.5."""
    rng = np.random.default_rng(8)
    bias, flats = make_night(tmp_path, '20230815', rng, 1.03, 10.3, [24000, 21800, 19800, 18000],
                             telescope='Artemis')
    record = ro.measure_night(bias, flats, '20230815', 'Artemis', PARAMS)
    assert record['calib_mode'] == 'NORMAL'
    noisy = ro.frame_values(record, 301.0, 10.4)
    quiet = ro.frame_values(record, 301.0, 6.5)
    assert noisy['READMODE'] == quiet['READMODE'] == 'NORMAL'
    assert noisy['GAINSRC'] == 'night' and noisy['GAIN'] == pytest.approx(1.03, rel=0.05)
    assert noisy['RONSRC'] == 'bias' and noisy['RDNOISE'] == pytest.approx(10.3 * noisy['GAIN'], rel=0.05)
    assert quiet['RONSRC'] == 'overscan' and quiet['RDNOISE'] == pytest.approx(6.5 * noisy['GAIN'])


def test_instrument_without_overscan_keeps_the_legacy_values():
    params = {'gain': 5.09, 'saturation_threshold': 11000, 'naxis2': 1280}
    record = ro.measure_night(['/x/20250815/b.fits'], [], '20250815', 'SPECULOOS-CALLISTO', params)
    assert record['andor'] is False and record['modes'] == {}
    values = ro.frame_values(record, None, None, legacy_ron_e=79.3)
    assert values == {'READMODE': 'DEFAULT', 'GAIN': 5.09, 'GAINSRC': 'config', 'RDNOISE': 79.3,
                      'RONSRC': 'legacy', 'SATURATE_RAW': 11000.0, 'OSNOISE': None}
    assert ro.calib_values(record, legacy_ron_e=79.3)['RDNOISE'] == 79.3


def test_header_cards_put_saturation_in_calibrated_units():
    values = ro.frame_values(_empty('Artemis', '20260420'), 300.0, 2.42)
    cards = {k: v for k, v, _ in ro.header_cards(values, pedestal=301.5)}
    assert cards['SATURATE'] == pytest.approx(37300.0 - 301.5)
    assert cards['RON'] == cards['RDNOISE']
    assert cards['READMODE'] == 'GAIN3.8' and cards['GAINSRC'] == 'table'
    assert cards['OSNOISE'] == pytest.approx(2.42)
    for key, value, comment in ro.header_cards(values, pedestal=301.5):
        fits.Card(key, value, comment).verify('exception')       # all are valid FITS cards


# --------------------------------------------------------------- record I/O

def test_load_or_measure_writes_reuses_and_refreshes(tmp_path):
    rng = np.random.default_rng(9)
    night = '20260423'
    bias, flats = make_night(tmp_path, night, rng, 1.04, 5.8, [24000, 21800, 19800, 18000, 16400])
    reddir = os.path.join(str(tmp_path), 'output', night, 'reduction')
    os.makedirs(reddir)
    with open(os.path.join(reddir, '1_bias.list'), 'w') as f:
        f.write('\n'.join(bias) + '\n')
    with open(os.path.join(reddir, '1_flat_I+z.list'), 'w') as f:
        f.write('\n'.join(flats) + '\n')

    record = ro.load_or_measure(reddir + '/', '1', PARAMS, 'SPECULOOS-IO')
    assert record['night'] == night and os.path.exists(os.path.join(reddir, 'readout.json'))
    assert record['modes']['NORMAL']['gain_accepted']

    # reused as written...
    stored = ro.read_record(reddir)
    stored['notes'].append('marker')
    ro.write_record(stored, reddir)
    assert 'marker' in ro.load_or_measure(reddir, '1', PARAMS, 'SPECULOOS-IO')['notes']
    # ...but not for another telescope, and T2 always measures afresh
    assert 'marker' not in ro.load_or_measure(reddir, '1', PARAMS, 'SPECULOOS-EUROPA', write=False)['notes']
    assert 'marker' not in ro.load_or_measure(reddir, '1', PARAMS, 'SPECULOOS-IO', refresh=True)['notes']


def test_a_reduction_dir_outside_the_pipeline_layout_is_reported(tmp_path):
    reddir = os.path.join(str(tmp_path), 'replay', 'reduction')
    os.makedirs(reddir)
    messages = []
    record = ro.load_or_measure(reddir, '1', PARAMS, 'Artemis', log=messages.append)
    assert record['night'] == 'replay'
    assert any('not YYYYMMDD' in n for n in record['notes']) and messages
    assert ro.night_from_reddir('/data/x/output/20251202/reduction/') == '20251202'


def test_a_failed_measurement_still_gives_a_usable_record(tmp_path, monkeypatch):
    reddir = os.path.join(str(tmp_path), 'output', '20251105', 'reduction')
    os.makedirs(reddir)
    monkeypatch.setattr(ro, 'night_lists', lambda *a: (_ for _ in ()).throw(OSError('disk gone')))
    record = ro.load_or_measure(reddir, '1', PARAMS, 'Artemis', log=lambda *a: None)
    assert any('disk gone' in n for n in record['notes'])
    assert ro.frame_values(record, 300.0, 2.42)['READMODE'] == 'GAIN3.8'


# ------------------------------------------------------- T6 calibrated frames

def test_reduce_file_writes_the_readout_keywords(tmp_path):
    """pipered.reduce_file with the real instrument config (Io, 2088 rows)."""
    from calibration import pipered
    rng = np.random.default_rng(10)
    raw = np.full((2088, 2048), 300.0)
    raw[:-15] += 7.0 + rng.poisson(1500.0, (2073, 2048))
    raw += rng.normal(0, 2.45, raw.shape)
    path = os.path.join(str(tmp_path), 'images', '20251105', 'Sp0000-S001-R001-C001-I+z.fts')
    os.makedirs(os.path.dirname(path))
    hdu = fits.PrimaryHDU(np.round(raw).astype(np.float32))
    hdu.header['TELESCOP'] = 'Artemis'
    hdu.header['EXPTIME'] = 30.0
    hdu.header['FILTER'] = 'I+z'
    hdu.writeto(path)

    trimmed_shape = (2066 - 22, 2048 - 2)
    bias = np.full(trimmed_shape, 7.0, dtype=np.float32)
    flat = {'I+z': np.ones(trimmed_shape, dtype=np.float32)}
    with fits.open(path) as hdul:
        params = pipered.get_instrument_parameters(hdul)
    record = ro.measure_night([], [], '20251105', 'Artemis', params)
    outdir = str(tmp_path) + '/out/'
    os.makedirs(outdir)
    dark = {'combined': np.zeros(trimmed_shape, dtype=np.float32)}
    pipered.reduce_file(path, outdir, params, 'bias.fits', bias, dark, flat, None, '9.2', 0.0, 'v3',
                        readout_record=record)
    with fits.open(outdir + 'procSp0000-S001-R001-C001-I+z.fits') as out:
        h = out[0].header
    assert h['READMODE'] == 'GAIN3.8' and h['GAINSRC'] == 'table'
    assert h['GAIN'] == pytest.approx(3.83)
    assert h['OSNOISE'] == pytest.approx(2.45, rel=0.08)
    assert h['RDNOISE'] == pytest.approx(h['OSNOISE'] * 3.83) and h['RON'] == h['RDNOISE']
    assert h['SATURATE'] == pytest.approx(37300.0 - h['OVERSCAN'] - 7.0, abs=1.0)


# ------------------------------------------------- T10 saturation and products

def _load_splightcurve():
    """SPlightcurve, with stand-ins for modules this machine may lack."""
    try:
        return importlib.import_module('diff_photometry.SPlightcurve')
    except ImportError:
        pass
    saved = dict(sys.modules)
    try:
        for name in ('fitsio', 'water_vapour', 'water_vapour.water_vapour', 'water_vapour.pwvGrid'):
            try:
                importlib.import_module(name)
            except ImportError:
                stub = types.ModuleType(name)
                stub.getLHATPROdata = None
                stub.FITS = object
                sys.modules[name] = stub
        sys.modules.pop('diff_photometry.SPlightcurve', None)
        return importlib.import_module('diff_photometry.SPlightcurve')
    finally:
        for name in list(sys.modules):
            if name not in saved and name != 'diff_photometry.SPlightcurve' \
                    and not name.startswith('diff_photometry'):
                sys.modules.pop(name, None)


@pytest.fixture(scope='module')
def splc():
    try:
        return _load_splightcurve()
    except Exception as e:                       # pragma: no cover - environment
        pytest.skip('SPlightcurve cannot be imported here: %s' % e)


def test_saturation_uses_each_frames_level_and_adds_the_sky(splc):
    # three stars x two frames, SATURATE in calibrated ADU of each frame
    peak = np.array([[30000., 36800.], [36000., 1000.], [61000., 61500.]])
    sky = np.full_like(peak, 800.)
    saturate = np.array([37000., 37000.])
    readmode = np.array(['GAIN3.8', 'GAIN3.8'])
    flags = splc.saturated_stars(peak, sky, saturate, readmode, 62000.)
    # star 0 reaches 37,600 with its sky; star 1 reaches 36,800 only
    assert list(flags) == [True, False, True]


def test_saturation_keeps_the_old_test_without_frame_levels(splc):
    peak = np.array([[61000., 62500.], [61900., np.nan]])
    flags_old = splc.saturated_stars(peak, None, None, None, 62000.)
    assert list(flags_old) == [True, False]
    # DEFAULT (instruments readout.py does not handle) and blank modes too
    level = np.array([1000., 1000.])
    assert list(splc.saturated_stars(peak, None, level, np.array(['DEFAULT', '']), 62000.)) == [True, False]
    # NaN SATURATE in some frames: those frames keep the old test
    mixed = splc.saturated_stars(np.array([[100., 62100.]]), np.zeros((1, 2)), np.array([50000., np.nan]),
                                 np.array(['NORMAL', 'NORMAL']), 62000.)
    assert list(mixed) == [True]


def test_init_all_lcurves_marks_the_flagged_stars(splc):
    """The per-star flags reach the light curves (a loop-local 'saturated'
    once shadowed them)."""
    nstar, nframe = 3, 20
    rng = np.random.default_rng(12)
    flux = [1000.0 + rng.normal(0, 5, nframe) for _ in range(nstar)]
    tel = types.SimpleNamespace(rcore=np.full(nframe, 4.0), exposure=np.full(nframe, 30.0),
                                gain=np.full(nframe, 3.83))
    ids = ['SP%06d' % i for i in range(nstar)]
    zeros = np.zeros(nstar)
    lcs = splc.init_all_lcurves(tel, ids, ids, ids, [f / np.median(f) for f in flux], flux,
                                np.array([False, True, False]), np.arange(nframe), np.arange(nframe),
                                np.arange(nframe), zeros, zeros, zeros, zeros, zeros,
                                [np.zeros(nframe)] * nstar, [np.zeros(nframe)] * nstar, zeros, zeros, zeros,
                                zeros, 5, 'I+z', '20251202', np.full(nframe, 30.0))
    assert [lc.saturated for lc in lcs] == [False, True, False]


def test_readout_summary_cards(splc):
    tel = types.SimpleNamespace(readmode=np.array(['GAIN3.8', 'GAIN3.8', 'PED2900']),
                                gain=np.array([3.83, 3.83, 6.5]), ron=np.array([9.3, 9.3, 41.0]))
    cards = {k: v for k, v, _ in splc.readout_summary_cards(tel)}
    assert cards == {'READMODE': 'GAIN3.8+PED2900', 'GAIN': 3.83, 'RDNOISE': 9.3}
    old = types.SimpleNamespace(readmode=np.array(['', '']), gain=np.array([1.0029, 1.0029]),
                                ron=np.array([6.3, 6.3]))
    assert {k: v for k, v, _ in splc.readout_summary_cards(old)}['READMODE'] == 'N'
