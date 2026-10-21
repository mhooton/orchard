"""
Tests for monitoring/monitoring.py, the 06:30 cron job that rewrites PipelineOutput/v2/<telescope>/monitoring.

Covered:
  - Callisto's camera on each night: the measured SPIRIT spells in jobqueue/config.py, not the hardcoded
    2022-06-01 / 2023-03-17 / 2025-03-10 the script used before, which put ~5 weeks under the wrong camera;
  - a night's readout mode, gain, gain source and bias level from readout.json, and the notes that disagree
    with readout.REGIMES; a record written for another night is ignored;
  - CCD temperature and PHOT_ZP from the night's *_5_diff.fits, PHOT_ZP only in the camera's survey filter;
  - the cache that spares unchanged _diff.fits files, re-reads rewritten ones and survives a corrupt one;
  - readout-mode changes between nights;
  - a whole run into a separate output directory, leaving the input tree alone.

Run:  pytest src/tests/test_monitoring.py -q
"""
import csv
import json
import os
import sys

import numpy as np
import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from monitoring import monitoring as mon          # noqa: E402

try:
    import fitsio  # noqa: F401
    HAVE_FITSIO = True
except ImportError:
    HAVE_FITSIO = False
needs_fitsio = pytest.mark.skipif(not HAVE_FITSIO, reason='fitsio is not installed (it is in orchard-server)')


# ------------------------------------------------------------------ helpers

def write_night(top, telescope, night, ron=6.5, overscan=330.0, dark=0.002, dat=('readoutnoise', 'overscan',
                                                                                   'darkcurrent')):
    red = os.path.join(top, telescope, 'output', night, 'reduction')
    os.makedirs(red, exist_ok=True)
    values = {'readoutnoise': ron, 'overscan': overscan, 'darkcurrent': dark}
    for name in dat:
        with open(os.path.join(red, name + '.dat'), 'w') as f:
            f.write('%r\n' % values[name])
    return red


def write_record(reddir, telescope, night, calib_mode, modes, expected=None, andor=True, notes=()):
    record = {
        'version': 1, 'telescope': telescope, 'night': night, 'andor': andor, 'config_gain': 1.0029,
        'config_saturation': 62000.0, 'expected_mode': expected or calib_mode, 'calib_mode': calib_mode,
        'modes': modes, 'notes': list(notes),
    }
    with open(os.path.join(reddir, 'readout.json'), 'w') as f:
        json.dump(record, f)


def write_diff(top, telescope, night, target, filt, temps, phot_zp=None):
    from astropy.io import fits
    d = os.path.join(top, telescope, 'output', night, target)
    os.makedirs(d, exist_ok=True)
    primary = fits.PrimaryHDU()
    if phot_zp is not None:
        primary.header['PHOT_ZP'] = phot_zp
    cols = fits.ColDefs([fits.Column(name='JD-OBS', format='D', array=2460000.5 + np.arange(len(temps)) / 1440.),
                         fits.Column(name='CCD-TEMP', format='E', array=np.asarray(temps, dtype=float))])
    imagelist = fits.BinTableHDU.from_columns(cols, name='IMAGELIST')
    path = os.path.join(d, '%s_%s_%s_5_diff.fits' % (target, filt, night))
    fits.HDUList([primary, imagelist]).writeto(path)
    return path


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


# ------------------------------------------------------------------ camera

@pytest.mark.parametrize('night,expected', [
    ('20220110', 'ANDOR'), ('20220111', 'SPIRIT'), ('20220112', 'ANDOR'), ('20220115', 'SPIRIT'),
    ('20220509', 'ANDOR'),
    ('20220510', 'SPIRIT'),     # the old dates said ANDOR until 2022-06-01
    ('20220531', 'SPIRIT'),
    ('20230315', 'SPIRIT'),
    ('20230316', 'ANDOR'),      # the old dates said SPIRIT until 2023-03-17
    ('20240601', 'ANDOR'),
    ('20250225', 'ANDOR'),
    ('20250226', 'SPIRIT'),     # the old dates said ANDOR until 2025-03-10
    ('20250309', 'SPIRIT'),
    ('20261019', 'SPIRIT'),
])
def test_callisto_camera_follows_the_measured_spells(night, expected):
    assert mon.get_instrument('Callisto', night) == expected


def test_other_telescopes_are_always_andor():
    for telescope in ('Artemis', 'Europa', 'Ganymede', 'Io'):
        for night in ('20220510', '20230101', '20250301'):
            assert mon.get_instrument(telescope, night) == 'ANDOR'


def test_split_by_instrument():
    rows = [{'instrument': 'ANDOR', 'night': '20230316'}, {'instrument': 'SPIRIT', 'night': '20250301'}]
    assert mon.split_by_instrument('Europa', rows) == {'': rows}
    assert mon.split_by_instrument('Callisto', rows) == {'_ANDOR': rows[:1], '_SPIRIT': rows[1:]}


# ------------------------------------------------------------------ readout.json

def test_readout_measured_gain_and_bias_level(tmp_path):
    red = write_night(str(tmp_path), 'Artemis', '20251030')
    write_record(red, 'Artemis', '20251030', 'GAIN3.8', {
        'GAIN3.8': {'gain': 3.81, 'gain_accepted': True, 'os_level': 300.0, 'bias_offset': 2.5,
                    'os_noise': 2.45, 'ron_adu': 2.4, 'n_bias': 10}})
    v = mon.read_readout(red, '20251030')
    assert v['readout_mode'] == 'GAIN3.8'
    assert v['gain'] == pytest.approx(3.81)
    assert v['gain_source'] == 'night'
    assert v['bias_level'] == pytest.approx(302.5)
    assert v['modes_seen'] == 'GAIN3.8'
    assert v['readout_notes'] == []


def test_readout_table_gain_and_notes(tmp_path):
    red = write_night(str(tmp_path), 'Europa', '20190801')
    write_record(red, 'SPECULOOS-EUROPA', '20190801', 'GAIN3.8', {
        'GAIN3.8': {'gain': None, 'gain_accepted': False, 'os_level': 310.0, 'os_noise': 2.5, 'n_bias': 6},
        'NORMAL': {'gain': None, 'gain_accepted': False, 'n_bias': 0}},
        expected='NORMAL',
        notes=['calibration frames borrowed from other nights are not measured',
               'calibration frames in GAIN3.8+NORMAL; table expects NORMAL'])
    v = mon.read_readout(red, '20190801')
    assert v['readout_mode'] == 'GAIN3.8'
    assert v['gain'] == pytest.approx(3.80)        # readout.REGIMES, Europa 2019-07-25..10-26
    assert v['gain_source'] == 'table'
    assert v['bias_level'] == pytest.approx(310.0)  # no bias offset recorded
    assert v['modes_seen'] == 'GAIN3.8+NORMAL'
    assert v['readout_notes'] == ['calibration frames in GAIN3.8+NORMAL; table expects NORMAL']


def test_readout_spirit_is_default(tmp_path):
    red = write_night(str(tmp_path), 'Callisto', '20250301')
    write_record(red, 'SPECULOOS-CALLISTO', '20250301', 'DEFAULT', {}, andor=False)
    v = mon.read_readout(red, '20250301')
    assert v['readout_mode'] == 'DEFAULT'
    assert v['gain_source'] == 'config'
    assert v['bias_level'] is None


def test_readout_record_of_another_night_is_ignored(tmp_path):
    red = write_night(str(tmp_path), 'Io', '20240510')
    write_record(red, 'SPECULOOS-IO', '20240509', 'NORMAL', {'NORMAL': {'gain': 1.04, 'gain_accepted': True}})
    assert mon.read_readout(red, '20240510') == {}
    os.remove(os.path.join(red, 'readout.json'))
    assert mon.read_readout(red, '20240510') == {}


def test_mode_changes():
    rows = [{'night': '20251001', 'readout_mode': 'NORMAL'},
            {'night': '20251010'},                                   # no readout.json
            {'night': '20251030', 'readout_mode': 'GAIN3.8'},
            {'night': '20251031', 'readout_mode': 'GAIN3.8'},
            {'night': '20260806', 'readout_mode': 'PED2900'},
            {'night': '20260808', 'readout_mode': 'GAIN3.8'}]
    assert mon.mode_changes(rows) == [('20251030', '20251001', 'NORMAL', 'GAIN3.8'),
                                      ('20260806', '20251031', 'GAIN3.8', 'PED2900'),
                                      ('20260808', '20260806', 'PED2900', 'GAIN3.8')]
    assert mon.mode_changes([{'night': '20251001'}]) == []


# ------------------------------------------------------------------ nights and .dat files

def test_collect_rows_and_missing_dates(tmp_path):
    top = str(tmp_path)
    write_night(top, 'Europa', '20230101', ron=6.1, overscan=320.0, dark=0.001)
    write_night(top, 'Europa', '20230102', dat=('readoutnoise', 'overscan'))       # no darkcurrent.dat
    os.makedirs(os.path.join(top, 'Europa', 'output', '20230103'))                  # no reduction directory
    os.makedirs(os.path.join(top, 'Europa', 'output', 'notanight', 'reduction'))
    nights = mon.list_nights('Europa', top)
    assert nights == ['20230101', '20230102', '20230103']
    rows, missing = mon.collect_monitoring_data('Europa', nights, top)
    assert [r['night'] for r in rows] == ['20230101']
    assert rows[0]['readout_noise'] == pytest.approx(6.1)
    assert rows[0]['zero_point'] == pytest.approx(320.0)
    assert rows[0]['dark_current'] == pytest.approx(0.001)
    assert rows[0]['instrument'] == 'ANDOR'
    assert 'readout_mode' not in rows[0] and 'ccd_temp' not in rows[0]
    assert missing == ['20230102']


def test_list_nights_of_a_telescope_with_no_output(tmp_path):
    assert mon.list_nights('Artemis', str(tmp_path)) == []


# ------------------------------------------------------------------ _diff.fits values

@needs_fitsio
def test_night_diff_values_survey_filter_only(tmp_path):
    top = str(tmp_path)
    write_diff(top, 'Europa', '20230101', '111', 'I+z', [-60.0, -60.2, -59.8], phot_zp=23.9)
    write_diff(top, 'Europa', '20230101', '222', 'r', [-61.0], phot_zp=22.0)
    date_dir = os.path.join(top, 'Europa', 'output', '20230101')
    v = mon.night_diff_values(date_dir, 'ANDOR', {}, set())
    assert v['ccd_temp'] == pytest.approx(-60.5)    # median of the two files' medians, any filter
    assert v['phot_zp'] == pytest.approx(23.9)      # I+z only
    assert v['n_diff'] == 2


@needs_fitsio
def test_night_diff_values_spirit_uses_zyj_and_tolerates_missing_zp(tmp_path):
    top = str(tmp_path)
    write_diff(top, 'Callisto', '20250301', '333', 'zYJ', [-70.0], phot_zp=21.4)
    write_diff(top, 'Callisto', '20250301', '444', 'I+z', [-70.0], phot_zp=23.0)
    write_diff(top, 'Callisto', '20250301', '555', 'zYJ', [-70.0])                  # no PHOT_ZP
    date_dir = os.path.join(top, 'Callisto', 'output', '20250301')
    v = mon.night_diff_values(date_dir, 'SPIRIT', {}, set())
    assert v['phot_zp'] == pytest.approx(21.4)
    assert v['ccd_temp'] == pytest.approx(-70.0)


def test_target_ids_with_underscores_keep_their_filter():
    m = mon.diff_name.search('TOI_4336_c_I+z_20240101_5_diff.fits')
    assert m.group('filter') == 'I+z' and m.group('night') == '20240101'


@needs_fitsio
def test_cache_spares_unchanged_files_and_rereads_rewritten_ones(tmp_path, monkeypatch):
    top = str(tmp_path)
    a = write_diff(top, 'Io', '20240101', '111', 'I+z', [-60.0], phot_zp=23.9)
    write_diff(top, 'Io', '20240102', '222', 'I+z', [-61.0], phot_zp=23.8)
    for night in ('20240101', '20240102'):
        write_night(top, 'Io', night)
    calls = []
    real = mon.read_diff_values
    monkeypatch.setattr(mon, 'read_diff_values', lambda p: calls.append(p) or real(p))
    nights = mon.list_nights('Io', top)
    cache = {}
    rows, _ = mon.collect_monitoring_data('Io', nights, top, cache)
    assert len(calls) == 2 and len(cache) == 2
    assert [r['phot_zp'] for r in rows] == pytest.approx([23.9, 23.8])

    calls.clear()
    cache = json.loads(json.dumps(cache))          # as saved and loaded between runs
    rows, _ = mon.collect_monitoring_data('Io', nights, top, cache)
    assert calls == []
    assert [r['ccd_temp'] for r in rows] == pytest.approx([-60.0, -61.0])

    os.remove(a)
    write_diff(top, 'Io', '20240101', '111', 'I+z', [-55.0], phot_zp=24.1)
    st = os.stat(a)
    os.utime(a, (st.st_atime, st.st_mtime + 10))
    rows, _ = mon.collect_monitoring_data('Io', nights, top, cache)
    assert calls == [a]
    assert rows[0]['ccd_temp'] == pytest.approx(-55.0) and rows[0]['phot_zp'] == pytest.approx(24.1)

    os.remove(a)
    mon.collect_monitoring_data('Io', nights, top, cache)
    assert a not in cache and len(cache) == 1      # a file that has gone leaves the cache


@needs_fitsio
def test_corrupt_diff_file_is_recorded_and_not_reread(tmp_path, monkeypatch):
    top = str(tmp_path)
    write_night(top, 'Ganymede', '20260228')
    d = os.path.join(top, 'Ganymede', 'output', '20260228', '999')
    os.makedirs(d)
    bad = os.path.join(d, '999_I+z_20260228_5_diff.fits')
    with open(bad, 'wb') as f:
        f.write(b'not a fits file')
    cache = {}
    rows, _ = mon.collect_monitoring_data('Ganymede', ['20260228'], top, cache)
    assert rows[0]['ccd_temp'] is None and rows[0]['phot_zp'] is None and rows[0]['n_diff'] == 1
    assert 'error' in cache[bad]
    monkeypatch.setattr(mon, 'read_diff_values', lambda p: pytest.fail('re-read an unchanged file'))
    mon.collect_monitoring_data('Ganymede', ['20260228'], top, cache)


# ------------------------------------------------------------------ whole runs

def test_run_callisto_into_a_separate_directory(tmp_path):
    top, out = str(tmp_path / 'v2'), str(tmp_path / 'out')
    write_night(top, 'Callisto', '20230314', ron=40.0, overscan=0.0, dark=120.0)    # SPIRIT
    write_night(top, 'Callisto', '20230316', ron=6.3, overscan=330.0, dark=0.003)   # ANDOR
    write_night(top, 'Callisto', '20230317', dat=('overscan',))
    mon.main(['--top-dir', top, '--out-dir', out, '--telescopes', 'Callisto', '--no-diff'])

    mdir = os.path.join(out, 'Callisto', 'monitoring')
    for suffix in ('_ANDOR', '_SPIRIT'):
        for name in ('monitoring_timeseries%s.csv', 'monitoring_timeseries%s.txt', 'readoutnoise%s_log.pdf',
                     'readoutnoise%s_linear.pdf', 'zeropoint%s.pdf', 'darkcurrent%s_log.pdf',
                     'darkcurrent%s_linear.pdf'):
            assert os.path.exists(os.path.join(mdir, name % suffix)), name % suffix
    andor = read_csv(os.path.join(mdir, 'monitoring_timeseries_ANDOR.csv'))
    spirit = read_csv(os.path.join(mdir, 'monitoring_timeseries_SPIRIT.csv'))
    assert [r['date'] for r in andor] == ['2023-03-16']
    assert [r['date'] for r in spirit] == ['2023-03-14']
    assert list(andor[0]) == mon.csv_columns
    assert float(andor[0]['readout_noise']) == pytest.approx(6.3)
    assert andor[0]['gain'] == '' and andor[0]['ccd_temp'] == ''
    with open(os.path.join(mdir, 'missing_data.txt')) as f:
        assert '20230317' in f.read()
    assert not os.path.exists(os.path.join(mdir, 'gain_ANDOR.pdf'))       # no readout.json anywhere
    assert not os.path.exists(os.path.join(mdir, mon.cache_name))         # --no-diff
    assert not os.path.exists(os.path.join(top, 'Callisto', 'monitoring'))


def test_run_with_readout_records_writes_gain_and_changes(tmp_path):
    top, out = str(tmp_path / 'v2'), str(tmp_path / 'out')
    normal = {'NORMAL': {'gain': 1.02, 'gain_accepted': True, 'os_level': 330.0, 'bias_offset': 1.0}}
    gain38 = {'GAIN3.8': {'gain': 3.82, 'gain_accepted': True, 'os_level': 300.0, 'bias_offset': 1.0}}
    write_record(write_night(top, 'Artemis', '20251020'), 'Artemis', '20251020', 'NORMAL', normal)
    write_night(top, 'Artemis', '20251025')
    write_record(write_night(top, 'Artemis', '20251030'), 'Artemis', '20251030', 'GAIN3.8', gain38,
                 notes=['GAIN3.8 gain 3.820 measured, table says 3.400'])
    mon.main(['--top-dir', top, '--out-dir', out, '--telescopes', 'Artemis'])

    mdir = os.path.join(out, 'Artemis', 'monitoring')
    assert os.path.exists(os.path.join(mdir, 'gain.pdf'))
    rows = read_csv(os.path.join(mdir, 'monitoring_timeseries.csv'))
    assert [(r['date'], r['readout_mode'], r['gain_source']) for r in rows] == [
        ('2025-10-20', 'NORMAL', 'night'), ('2025-10-25', '', ''), ('2025-10-30', 'GAIN3.8', 'night')]
    assert float(rows[2]['bias_level']) == pytest.approx(301.0)
    with open(os.path.join(mdir, 'readout_changes.txt')) as f:
        text = f.read()
    assert '2 of 3 nights have a readout.json' in text
    assert '20251030  NORMAL -> GAIN3.8  (previous record 20251020)' in text
    assert '20251030  GAIN3.8 gain 3.820 measured, table says 3.400' in text
    with open(os.path.join(mdir, mon.cache_name)) as f:
        assert json.load(f) == {}                                         # no _diff.fits in this tree


@needs_fitsio
def test_run_writes_ccd_temperature_and_zero_point_plots(tmp_path):
    top, out = str(tmp_path / 'v2'), str(tmp_path / 'out')
    for night, temp, zp in (('20240101', -60.0, 23.9), ('20240102', -60.5, 23.85)):
        write_night(top, 'Europa', night)
        write_diff(top, 'Europa', night, '111', 'I+z', [temp], phot_zp=zp)
    mon.main(['--top-dir', top, '--out-dir', out, '--telescopes', 'Europa'])
    mdir = os.path.join(out, 'Europa', 'monitoring')
    assert os.path.exists(os.path.join(mdir, 'ccdtemp.pdf'))
    assert os.path.exists(os.path.join(mdir, 'photzp.pdf'))
    rows = read_csv(os.path.join(mdir, 'monitoring_timeseries.csv'))
    assert [float(r['phot_zp']) for r in rows] == pytest.approx([23.9, 23.85])
    assert [r['n_diff'] for r in rows] == ['1', '1']
    with open(os.path.join(mdir, mon.cache_name)) as f:
        assert len(json.load(f)) == 2
