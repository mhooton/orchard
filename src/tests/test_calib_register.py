"""
Tests for the calibration register (calibration/calregister.py), the choice of a night's calibrations from it
(calibration/calselect.py), T1's use of both (main/createlists.py) and pipedark.py's exposure-matched masters.

Covered:
  - a frame's register fields: configuration (shape, binning, orientation), overscan, readout mode, level,
    saturation; SPIRIT's measured pixel ceiling;
  - the usability rules: saturated and faint flats, bias and dark outliers;
  - the builder's prefilter: ESO types, frames unpacked from datacubes (both naming forms), ACP names;
  - camera runs (ANDOR -> SPIRIT -> ANDOR) and hand-maintained breaks that apply to flats only;
  - the choice: nearest night either side, the past on a tie, never another camera, configuration or readout
    mode, SPIRIT darks per exposure and temperature with a covering night, flags and errors, stale rows,
    compressed frames;
  - createlists: Callisto 20250407 (job 1043), where the old search took ANDOR flats for SPIRIT frames;
  - pipedark: exposure keys rounded as pipered.py looks them up, MasterDark.fits the highest exposure;
  - the queue's hook only updates an existing register and never fails a job.

Run:  pytest src/tests/test_calib_register.py -q
"""
import json
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import calregister as cr      # noqa: E402
from calibration import calselect as cs        # noqa: E402

ANDOR = (2088, 2048)       # (NAXIS2, NAXIS1) of a stored Andor frame
SPIRIT = (1280, 1024)


def write_frame(path, imagetyp, level, shape=ANDOR, exptime=10.0, filt='I+z', telescop='SPECULOOS-IO',
                xpixsz=13.5, flip='Rotate 90 CCW', astra=None, ccd_temp=-60.0, pedestal=300.0, noise=6.0,
                obj=None, rng=None, binning=1):
    """A raw frame: pedestal + level in the image, pedestal alone in the 15 overscan rows (Andor)."""
    rng = rng or np.random.default_rng(abs(hash(path)) % (2 ** 32))
    ny, nx = shape
    data = rng.normal(pedestal + level, noise, shape)
    if shape == ANDOR:
        data[-15:] = rng.normal(pedestal, noise, (15, nx))
    h = fits.Header()
    h['IMAGETYP'] = imagetyp
    h['EXPTIME'] = exptime
    h['FILTER'] = filt
    h['TELESCOP'] = telescop
    h['XBINNING'] = binning
    h['YBINNING'] = binning
    h['XPIXSZ'] = xpixsz * binning
    h['CCD-TEMP'] = ccd_temp
    h['DATE-OBS'] = '2025-04-07T23:00:00'
    if obj:
        h['OBJECT'] = obj
    if astra:
        for k, v in astra.items():
            h[k] = v
    elif flip:
        h['FLIPSTAT'] = flip
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fits.writeto(path, np.clip(np.round(data), 0, 65535).astype('uint16'), h, overwrite=True)
    return path


# ------------------------------------------------------------------------------------- frames

def test_andor_frame_fields_and_statistics(tmp_path):
    p = write_frame(str(tmp_path / 'Observations/Io/images/20240101/b.fits'), 'Bias Frame', 0.0)
    f = cr.read_frame(p, 'Io', '20240101')
    assert f['error'] is None and f['kind'] == 'bias'
    assert f['config'] == '2048x2088 b1x1 maxim:Rotate 90 CCW'
    assert f['detector'] == '13.5um' and f['instrument'] == 'andor'
    assert f['readmode'] == 'NORMAL'
    assert f['os_level'] == pytest.approx(300.0, abs=1.0)
    assert f['os_noise'] == pytest.approx(6.0, rel=0.15)
    assert abs(f['level']) < 1.0
    assert f['saturation'] == pytest.approx(62000.0)


def test_orientation_labels_keep_software_eras_apart():
    assert cr.orientation({'FLIPSTAT': 'Rotate 90 CCW'}) == 'maxim:Rotate 90 CCW'
    assert cr.orientation({'ASTRAROT': 90}) == 'astra:rot90'
    assert cr.orientation({'ASTRAMIR': 'Y'}) == 'astra:mirY'
    assert cr.orientation({'FLIPSTAT': 'Mirror'}) == 'maxim:Mirror'
    assert cr.orientation({}) == 'none'


def test_spirit_frames_use_the_measured_ceiling(tmp_path):
    d = tmp_path / 'Observations/Callisto/images/20250407'
    flat = write_frame(str(d / 'f.fits'), 'FLAT', 8000.0, shape=SPIRIT, telescop='SPECULOOS-CALLISTO',
                       xpixsz=12.0, astra={'ASTRAMIR': 'Y'}, pedestal=2810.0, filt='zYJ')
    f = cr.read_frame(flat, 'Callisto', '20250407')
    assert f['instrument'] == 'spirit' and f['readmode'] == 'DEFAULT'
    assert f['config'] == '1024x1280 b1x1 astra:mirY' and f['detector'] == '12.0um'
    assert f['saturation'] == cr.PIXEL_CEILING['spirit']
    # a 10,800 ADU flat is above the configured 11,000 "saturation" only if the pedestal is ignored
    rows = [f]
    cr.flag_night(rows)
    assert f['usable'] == 1


def test_binned_frames_have_their_own_configuration(tmp_path):
    p = write_frame(str(tmp_path / 'Observations/Artemis/images/20250115/AutoFlat/f.fts'), 'FLAT', 20000.0,
                    shape=(1044, 1024), telescop='Artemis', binning=2)
    f = cr.read_frame(p, 'Artemis', '20250115')
    assert f['config'].startswith('1024x1044 b2x2') and f['detector'] == '13.5um'


# ------------------------------------------------------------------------------------- usability

def _row(kind, level, **kw):
    r = {'kind': kind, 'level': level, 'config': 'C', 'readmode': 'NORMAL', 'os_level': 300.0, 'sat_frac': 0.0,
         'saturation': 62000.0, 'exptime': kw.pop('exptime', 10.0), 'error': None}
    r.update(kw)
    return r


def test_flag_night_rules():
    rows = [_row('bias', 1.0), _row('bias', 2.0), _row('bias', 1.5), _row('bias', 300.0),
            _row('dark', 5.0), _row('dark', 5.5), _row('dark', 6.0), _row('dark', 900.0),
            _row('flat', 30000.0), _row('flat', 500.0), _row('flat', 55000.0), _row('flat', 30000.0, sat_frac=0.2),
            _row('flat', 30000.0, error='bad')]
    cr.flag_night(rows)
    reasons = [r['reason'] for r in rows]
    assert reasons[:3] == [None] * 3 and reasons[3] == 'level outlier'
    assert reasons[4:7] == [None] * 3 and reasons[7] == 'level outlier'
    assert reasons[8:] == [None, 'too faint', 'near saturation', 'saturated', 'unreadable']
    assert rows[8]['zero'] == pytest.approx(1.75) and rows[8]['signal'] == pytest.approx(30000.0 - 1.75)


def test_darks_are_judged_within_their_exposure():
    rows = [_row('dark', 20.0, exptime=10), _row('dark', 21.0, exptime=10), _row('dark', 22.0, exptime=10),
            _row('dark', 2200.0, exptime=120), _row('dark', 2210.0, exptime=120), _row('dark', 2205.0, exptime=120)]
    cr.flag_night(rows)
    assert all(r['usable'] for r in rows)


# ------------------------------------------------------------------------------------- the builder's prefilter

def test_candidates_from_eso_types_and_acp_names():
    rel = ['SPECU1.2024-04-10T22:53:25.660.fits', 'SPECU1.2024-04-10T10:30:28.340.fits',
           'Sp0746+2000-S001-R001-C001-I+z.fts', 'AutoFlat/AutoFlat-Dusk-I+z-Bin1-001.fts',
           'Calibration/Bias-001.fts', 'Bias-002.fts']
    types = {cr.dp_key('SPECU1.2024-04-10T22:53:25.660'): 'OBJECT', cr.dp_key('SPECU1.2024-04-10T10:30:28.340'): 'FLAT'}
    got = cr.candidate_paths('Io', '20240410', rel, types, {})
    assert got == ['SPECU1.2024-04-10T10:30:28.340.fits', 'AutoFlat/AutoFlat-Dusk-I+z-Bin1-001.fts',
                   'Calibration/Bias-001.fts', 'Bias-002.fts']
    assert cr.candidate_paths('Saint-Ex', '20210702', rel, {}, {}) == rel


def test_candidates_on_a_datacube_night_follow_the_cube_starts():
    start = lambda s: cr._mjd(__import__('datetime').datetime.strptime(s, '%Y%m%dT%H%M%S'))
    cubes = {('Callisto', '20260310'): [(start('20260310T112318'), 'DARK'), (start('20260310T230000'), 'OBJECT'),
                                        (start('20260311T100000'), 'FLAT')]}
    rel = ['SPECU4.2026-03-10T11:23:18.438.0.fits',     # in the dark cube, "<time>.<k>" naming
           'SPECU4.2026-03-10T11:23:19.38.0.fits',
           'SPECU4.2026-03-10T23:10:00.100.fits',       # in the OBJECT cube
           'SPECU4.2026-03-11T10:00:05.000.fits',       # in the flat cube
           'SPECU4.2026-03-10T10:00:00.000.fits']       # before any cube: read
    got = cr.candidate_paths('Callisto', '20260310', rel, {}, cubes)
    assert got == [rel[0], rel[1], rel[3], rel[4]]


def test_load_eso_index_takes_cube_starts_from_origfile(tmp_path):
    p = tmp_path / 'eso.csv'
    p.write_text('telescope,night,dp_id,dp_type,mjd_obs,exposure,filter_path,origfile,cube\n'
                 'Callisto,20260601,SPECU4.2026-06-01T10:41:22.992,FLAT,61192.44540,60,ZYJ,'
                 'SPECU4.20260601T103622_C_flat_zYJ_60s.fits,1\n')
    by_id, cubes = cr.load_eso_index(str(p))
    assert by_id == {'SPECU4.2026-06-01T10:41:22.992': 'FLAT'}
    assert cr.dp_key('SPECU4.2026-02-04T23:48:53.95.0') == 'SPECU4.2026-02-04T23:48:53.095'
    assert cr.dp_key('SPECU4.2026-06-01T22:19:24.932') == 'SPECU4.2026-06-01T22:19:24.932'
    assert cr.dp_key('SPECU4.2026-06-01T22:19:24.9') == 'SPECU4.2026-06-01T22:19:24.900'
    assert cr.candidate_paths('Callisto', '20260204', ['SPECU4.2026-02-04T23:48:53.95.0.fits'],
                              {'SPECU4.2026-02-04T23:48:53.095': 'FLAT'}, cubes) == ['SPECU4.2026-02-04T23:48:53.95.0.fits']
    (s, t), = cubes[('Callisto', '20260601')]
    assert t == 'FLAT' and s == pytest.approx(cr._name_mjd('SPECU4.2026-06-01T10:36:22'))


# ------------------------------------------------------------------------------------- the register

def test_register_replaces_a_night_and_refuses_to_be_created_by_t1(tmp_path):
    root = tmp_path
    d = root / 'Observations/Io/images/20240101'
    for i in range(3):
        write_frame(str(d / 'b{}.fits'.format(i)), 'Bias Frame', 0.0)
    path = str(root / 'reg.sqlite')
    with pytest.raises(cr.RegisterError):
        cr.open_existing(path)
    reg = cr.Register(path, create=True, check_local=False)
    row, rows = cr.register_night(reg, str(root), 'Io', '20240101')
    assert row['n_bias'] == 3 and len(reg.frames('Io')) == 3
    os.remove(str(d / 'b2.fits'))
    cr.register_night(reg, str(root), 'Io', '20240101')
    assert len(reg.frames('Io')) == 2
    reg.close()


def test_build_reads_candidates_with_a_process_pool_and_resumes(tmp_path):
    root = tmp_path
    for n in ('20180301', '20180302'):
        d = root / 'Observations/Io/images' / n
        for i in range(2):
            write_frame(str(d / 'Calibration/Bias-{}.fts'.format(i)), 'Bias Frame', 0.0)
        write_frame(str(d / 'Sp0000+0000-S001-R001-C001-I+z.fts'), 'Light Frame', 500.0, obj='Sp0000+0000')
    reg = cr.Register(str(root / 'reg.sqlite'), create=True, check_local=False)
    t = cr.build(reg, str(root), ['Io'], workers=2, log=lambda *a: None)
    assert t['nights'] == 2 and t['cal'] == 4 and len(reg.frames('Io')) == 4
    sci = json.loads(reg.night('Io', '20180301')['science'])
    assert list(sci) == ['2048x2088 b1x1 maxim:Rotate 90 CCW'] and reg.night('Io', '20180301')['detector'] == '13.5um'
    assert cr.build(reg, str(root), ['Io'], workers=2, log=lambda *a: None)['skipped'] == 2
    reg.close()


def test_readonly_register_is_not_written(tmp_path, monkeypatch):
    path = str(tmp_path / 'reg.sqlite')
    cr.Register(path, create=True, check_local=False).close()
    monkeypatch.setenv(cr.READONLY_ENV, '1')
    reg = cr.open_existing(path)
    assert reg.readonly
    with pytest.raises(cr.RegisterError):
        reg.put_night({'telescope': 'Io', 'night_dir': '20240101', 'night': '20240101'}, [])
    reg.close()


def test_hook_needs_an_existing_register_and_never_raises(tmp_path, monkeypatch):
    monkeypatch.setenv(cr.REGISTER_ENV, str(tmp_path / 'elsewhere.sqlite'))     # ignored by the queue's hook
    msgs = []
    assert cr.hook_register_night(str(tmp_path), 'Io', '20240101', log=msgs.append) is None
    assert not os.path.exists(cr.default_register_path(str(tmp_path)))
    path = cr.default_register_path(str(tmp_path))
    cr.Register(path, create=True, check_local=False).close()
    monkeypatch.setattr(cr, 'register_night', lambda *a, **k: (_ for _ in ()).throw(OSError('disk gone')))
    assert cr.hook_register_night(str(tmp_path), 'Io', '20240101', log=msgs.append) is None
    assert any('could not register' in m for m in msgs)


# ------------------------------------------------------------------------------------- epochs

def test_camera_runs_and_flat_only_breaks():
    nights = [('20220301', '13.5um'), ('20220302', '13.5um'), ('20220510', '12.0um'), ('20230315', '12.0um'),
              ('20230318', '13.5um'), ('20230908', '13.5um'), ('20250226', '12.0um'), ('20250407', '12.0um')]
    e = cs.Epochs(nights, [('20230901', 'flat')])
    assert e.epoch('20220302', 'bias') != e.epoch('20230318', 'bias')          # ANDOR either side of SPIRIT
    assert e.epoch('20220510', 'dark') != e.epoch('20250407', 'dark')          # SPIRIT spell 1 vs 2
    assert e.epoch('20250226', 'flat') == e.epoch('20250407', 'flat')
    assert e.epoch('20230318', 'bias') == e.epoch('20230908', 'bias')          # the break is for flats only
    assert e.epoch('20230318', 'flat') != e.epoch('20230908', 'flat')


def test_breaks_file_parses():
    breaks = cs.load_breaks()
    assert ('20231211', 'flat') in breaks['Io']


# ------------------------------------------------------------------------------------- choosing

def reg_row(night, kind, n=0, config='C', readmode='NORMAL', exptime=10.0, filt='I+z', ccd_temp=-60.0,
            usable=1, compressed=0, instrument='andor'):
    rel = '{}_{}_{}.fits'.format(kind, n, int(exptime))
    return {'telescope': 'T', 'night_dir': night, 'night': night, 'relpath': rel, 'stem': rel[:-5], 'kind': kind,
            'config': config, 'readmode': readmode, 'exptime': exptime, 'filter': filt, 'ccd_temp': ccd_temp,
            'usable': usable, 'compressed': compressed, 'instrument': instrument}


def many(night, kind, k, **kw):
    return [reg_row(night, kind, n=i, **kw) for i in range(k)]


PROFILE = {'config': 'C', 'detector': '13.5um', 'instrument': 'andor', 'telescop': 'X', 'n': 50,
           'filters': {'I+z': 50}, 'exposures': {10: 50}, 'ccd_temp': -60.0, 'readmode': 'NORMAL',
           'readmodes': {'NORMAL': 5}, 'other_configs': {}}


def chooser(rows, profile=PROFILE, night='20250407', dark_matching=False, nights=None, breaks=(),
            resolve=None):
    by_kind = {}
    for r in rows:
        by_kind.setdefault(r['kind'], []).append(r)
    nights = nights or [(r['night'], '13.5um') for r in rows] + [(night, '13.5um')]
    return cs.Chooser('T', night, dict(profile), by_kind, cs.Epochs(nights, breaks),
                      resolve or (lambda r: '/x/{}/{}'.format(r['night_dir'], r['relpath'])), dark_matching)


def test_nearest_night_either_side_past_wins_a_tie():
    rows = many('20250405', 'bias', 5) + many('20250409', 'bias', 5) + many('20250410', 'bias', 9)
    assert chooser(rows).bias()['night'] == '20250405'
    rows = many('20250404', 'bias', 5) + many('20250409', 'bias', 5)
    e = chooser(rows).bias()
    assert e['night'] == '20250409' and e['delta_days'] == 2           # the future when it is nearer


def test_own_night_wins_when_it_has_enough():
    rows = many('20250407', 'flat', 4) + many('20250406', 'flat', 30)
    e = chooser(rows).flats()['I+z']
    assert e['night'] == '20250407' and e['own_night'] and e['delta_days'] == 0


def test_too_few_frames_or_unusable_frames_are_passed_over():
    rows = many('20250407', 'bias', 4) + many('20250406', 'bias', 5, usable=0) + many('20250403', 'bias', 5)
    assert chooser(rows).bias()['night'] == '20250403'


def test_never_another_camera_configuration_or_readout_mode():
    rows = (many('20250406', 'bias', 9, config='ANDOR') + many('20250405', 'bias', 9, readmode='FAST')
            + many('20250301', 'bias', 9))
    assert chooser(rows).bias()['night'] == '20250301'
    # another camera run between the nights: nothing qualifies, the night fails
    nights = [('20250301', '13.5um'), ('20250320', '12.0um'), ('20250407', '13.5um')]
    c = chooser(many('20250301', 'bias', 9), nights=nights)
    assert c.bias() is None and c.errors


def test_flats_never_combine_nights_and_missing_filters_are_flagged():
    profile = dict(PROFILE, filters={'zYJ': 40, 'J': 10})
    rows = many('20250407', 'flat', 2, filt='zYJ') + many('20250401', 'flat', 2, filt='zYJ') \
        + many('20250330', 'flat', 6, filt='zYJ')
    c = chooser(rows, profile=profile)
    flats = c.flats()
    assert flats['zYJ']['night'] == '20250330' and flats['zYJ']['n'] == 6
    assert 'J' not in flats and any("in J" in f for f in c.flags) and not c.errors
    c = chooser([], profile=profile)
    assert c.flats() == {} and c.errors


def test_andor_darks_outside_the_temperature_tolerance_are_flagged():
    rows = many('20250406', 'dark', 5, ccd_temp=-40.0) + many('20250301', 'dark', 5)
    c = chooser(rows)
    assert c.darks()['night'] == '20250301'
    c = chooser(many('20250406', 'dark', 5, ccd_temp=-40.0))
    assert c.darks()['night'] == '20250406' and any('outside' in f for f in c.flags)


SPIRIT_PROFILE = dict(PROFILE, instrument='spirit', exposures={10: 100, 60: 50}, ccd_temp=-60.0, readmode='DEFAULT')


def spirit_darks(night, exposures, k=3, ccd_temp=-60.0):
    return sum((many(night, 'dark', k, exptime=e, readmode='DEFAULT', ccd_temp=ccd_temp, instrument='spirit')
                for e in exposures), [])


def test_spirit_needs_fewer_bias_frames_than_andor():
    rows = many('20250407', 'bias', 2, readmode='DEFAULT') + many('20250401', 'bias', 9, readmode='DEFAULT')
    assert chooser(rows, profile=SPIRIT_PROFILE).bias()['night'] == '20250407'
    rows = many('20250407', 'bias', 2) + many('20250401', 'bias', 9)
    assert chooser(rows).bias()['night'] == '20250401'


def test_spirit_darks_per_exposure_and_temperature():
    rows = (spirit_darks('20250407', [10]) + spirit_darks('20250405', [60]) + spirit_darks('20250406', [60],
                                                                                            ccd_temp=-50.0))
    c = chooser(rows, profile=SPIRIT_PROFILE, dark_matching=True)
    d = c.darks()
    assert d['by_exposure']['10']['night'] == '20250407'
    assert d['by_exposure']['60']['night'] == '20250405'               # -50 C darks are not used for -60 C
    assert d['covering_night'] is None and not d['missing_exposures']


def test_spirit_prefers_one_night_covering_every_exposure_within_a_few_days():
    rows = spirit_darks('20250407', [10]) + spirit_darks('20250406', [60]) + spirit_darks('20250404', [10, 60])
    d = chooser(rows, profile=SPIRIT_PROFILE, dark_matching=True).darks()
    assert d['covering_night'] == '20250404'
    assert {e['night'] for e in d['by_exposure'].values()} == {'20250404'}
    rows = spirit_darks('20250407', [10]) + spirit_darks('20250406', [60]) + spirit_darks('20250320', [10, 60])
    d = chooser(rows, profile=SPIRIT_PROFILE, dark_matching=True).darks()
    assert d['covering_night'] is None and d['by_exposure']['10']['night'] == '20250407'
    # an exposure no night has darks at does not stop a covering night being used for the others
    profile = dict(SPIRIT_PROFILE, exposures={3: 20, 10: 100, 60: 50})
    rows = spirit_darks('20250407', [10]) + spirit_darks('20250406', [60]) + spirit_darks('20250405', [10, 60])
    d = chooser(rows, profile=profile, dark_matching=True).darks()
    assert d['covering_night'] == '20250405' and d['missing_exposures'] == [3]


def test_spirit_keeps_the_primary_nights_other_exposures_but_never_mixes_nights():
    rows = spirit_darks('20250407', [10, 60, 120]) + spirit_darks('20250406', [30])
    profile = dict(SPIRIT_PROFILE, exposures={10: 100, 60: 50})
    d = chooser(rows, profile=profile, dark_matching=True).darks()
    assert d['by_exposure']['extra'] == {'night': '20250407', 'exposures': [120]}
    assert len(d['frames']) == 9 and not any('20250406' in f for f in d['frames'])
    # an exposure taken from another night is not topped up from the primary night
    rows = spirit_darks('20250407', [10, 60]) + spirit_darks('20250406', [10, 60, 120]) \
        + many('20250407', 'dark', 3, exptime=60, readmode='DEFAULT', instrument='spirit', usable=0)
    d = chooser(rows, profile=dict(profile, exposures={10: 100, 60: 50}), dark_matching=True).darks()
    assert {e['night'] for k, e in d['by_exposure'].items() if k != 'extra'} == {'20250407'}


def test_spirit_exposure_without_darks_is_flagged_not_dropped():
    rows = spirit_darks('20250407', [10, 120])
    c = chooser(rows, profile=SPIRIT_PROFILE, dark_matching=True)
    d = c.darks()
    assert d['missing_exposures'] == [60] and any('60 s' in f for f in c.flags) and not c.errors
    assert len(d['frames']) == 6 and d['by_exposure']['extra']['exposures'] == [120]   # the scaled fallback
    c = chooser(spirit_darks('20250407', [120]), profile=SPIRIT_PROFILE, dark_matching=True)
    d = c.darks()
    assert 'any' in d['by_exposure'] and d['missing_exposures'] == [10, 60] and not c.errors
    c = chooser(spirit_darks('20250407', [120], ccd_temp=-50.0), profile=SPIRIT_PROFILE, dark_matching=True)
    assert c.darks() is None and c.errors


def test_stale_rows_and_compressed_frames_are_skipped():
    rows = many('20250406', 'bias', 5) + many('20250405', 'bias', 5, compressed=1) + many('20250401', 'bias', 5)
    c = chooser(rows, resolve=lambda r: None if r['night'] == '20250406' else '/x')
    assert c.bias()['night'] == '20250401' and c.stale == {'20250406'}


def test_select_records_errors_without_science_and_describe_lines(tmp_path):
    rec = cs.select(None, 'T', '20250407', str(tmp_path), [], [], [])
    assert rec['errors'] == ['no readable science frames']
    rec = {'bias': {'night': '20250405', 'delta_days': -2, 'n': 7, 'own_night': False},
           'dark': {'night': '20250407', 'delta_days': 0, 'n': 9, 'own_night': True},
           'flat': {'zYJ': {'night': '20250409', 'delta_days': 2, 'n': 6, 'own_night': False}},
           'flags': ['x'], 'errors': []}
    text = '\n'.join(cs.describe(rec))
    assert 'bias           20250405 (-2 d) 7 frames  [borrowed, |delta| = 2 d]' in text
    assert 'flat zYJ       20250409 (+2 d)' in text and 'WARNING: x' in text
    line = cs.report_line(rec)
    assert line == 'Calibrations: bias 20250405 (-2 d); dark 20250407; flat zYJ 20250409 (+2 d)  [1 warning]'


# ------------------------------------------------------------------------------------- createlists (T1)

def spirit(path, imagetyp, level, **kw):
    kw.setdefault('telescop', 'SPECULOOS-CALLISTO')
    return write_frame(path, imagetyp, level, shape=SPIRIT, xpixsz=12.0, astra={'ASTRAMIR': 'Y'},
                       pedestal=2810.0, noise=15.0, flip=None, **kw)


def andor_callisto(path, imagetyp, level, **kw):
    return write_frame(path, imagetyp, level, telescop='', astra={'ASTRAROT': 90}, flip=None, **kw)


@pytest.fixture
def callisto_archive(tmp_path):
    """Callisto around 20250407 (SPIRIT): ANDOR Clear flats on 20230908, SPIRIT zYJ flats on 20250401 and two
    on the night itself, darks at 10 s only on the night, at 60 s on 20250405; science in zYJ, Clear and J."""
    root = tmp_path / 'root'
    img = root / 'Observations/Callisto/images'
    for i in range(6):
        andor_callisto(str(img / '20230908/SPECU4.2023-09-08T23:01:{:02d}.000.fits'.format(i)), 'FLAT', 30000.0,
                       filt='Clear')
        andor_callisto(str(img / '20230908/SPECU4.2023-09-08T10:00:{:02d}.000.fits'.format(i)), 'Bias Frame', 0.0)
    for i in range(5):
        spirit(str(img / '20250401/SPECU4.2025-04-01T10:18:{:02d}.000.fits'.format(i)), 'FLAT', 8000.0,
               filt='zYJ', exptime=11.0)
        spirit(str(img / '20250405/SPECU4.2025-04-05T11:00:{:02d}.000.fits'.format(i)), 'Dark Frame', 1100.0,
               filt='Dark', exptime=60.0)
    night = img / '20250407'
    for i in range(5):
        spirit(str(night / 'SPECU4.2025-04-08T10:48:{:02d}.000.fits'.format(i)), 'Bias Frame', 0.0, filt='Dark',
               exptime=0.0)
        spirit(str(night / 'SPECU4.2025-04-08T10:49:{:02d}.000.fits'.format(i)), 'Dark Frame', 200.0,
               filt='Dark', exptime=10.0)
    for i in range(2):
        spirit(str(night / 'SPECU4.2025-04-08T10:20:{:02d}.000.fits'.format(i)), 'FLAT', 8000.0, filt='zYJ',
               exptime=10.0)
    for i in range(6):
        spirit(str(night / 'SPECU4.2025-04-08T01:00:{:02d}.000.fits'.format(i)), 'Light Frame', 500.0,
               filt='zYJ', exptime=10.0, obj='Sp1614-0250')
        spirit(str(night / 'SPECU4.2025-04-08T03:00:{:02d}.000.fits'.format(i)), 'Light Frame', 500.0,
               filt='Clear', exptime=60.0, obj='TOI-1227')
        spirit(str(night / 'SPECU4.2025-04-08T05:00:{:02d}.000.fits'.format(i)), 'Light Frame', 500.0,
               filt='J', exptime=60.0, obj='TOI-1227B')
    reg = cr.Register(cr.default_register_path(str(root)), create=True, check_local=False)
    for d in ('20230908', '20250401', '20250405'):
        cr.register_night(reg, str(root), 'Callisto', d)
    reg.close()
    out = root / 'PipelineOutput/v3/Callisto'
    (out / 'output/20250407/reduction').mkdir(parents=True)
    (out / 'reports').mkdir(parents=True)          # QC writes the file reports/QC
    return root


def run_t1(root, night='20250407', tel='Callisto'):
    pytest.importorskip('tabulate')                     # createlists imports reporting.QC
    from main import createlists
    img = str(root / 'Observations' / tel / 'images' / night)
    out = str(root / 'PipelineOutput/v3' / tel)
    createlists.main(['createlists.py', img, out, out + '/output', out + '/reports', 'IMAGE', 'fits', '1',
                      night, tel])
    red = os.path.join(out, 'output', night, 'reduction')
    lists = {}
    for name in os.listdir(red):
        if name.endswith('.list'):
            with open(os.path.join(red, name)) as f:
                lists[name] = [line.strip() for line in f if line.strip()]
    return red, lists


def test_createlists_job_1043_never_takes_andor_flats(callisto_archive, monkeypatch):
    monkeypatch.setattr(cr, 'mount_fstype', lambda path, mounts_file=None: None)
    monkeypatch.delenv(cr.REGISTER_ENV, raising=False)
    red, lists = run_t1(callisto_archive)
    with open(os.path.join(red, '1_calib_selection.json')) as f:
        rec = json.load(f)
    # zYJ: tonight has only two flats, so the five of 20250401 (6 days back); never 20250407's two added to them
    assert [os.path.basename(os.path.dirname(p)) for p in lists['1_flat_zYJ.list']] == ['20250401'] * 5
    # Clear: the only Clear flats are ANDOR ones: flagged, not used; J has none at all
    assert '1_flat_Clear.list' not in lists and '1_flat_J.list' not in lists
    assert any('in Clear' in f for f in rec['flags']) and any('in J' in f for f in rec['flags'])
    assert not any('20230908' in p for v in lists.values() for p in v)
    # darks: 10 s from the night itself, 60 s from 20250405
    assert rec['dark']['by_exposure']['10']['night'] == '20250407'
    assert rec['dark']['by_exposure']['60']['night'] == '20250405'
    assert len(lists['1_dark.list']) == 10 and len(lists['1_bias.list']) == 5
    assert 'extra' not in rec['dark']['by_exposure']
    # the night's own frames were registered by T1
    reg = cr.open_existing(cr.default_register_path(str(callisto_archive)), readonly=True)
    assert reg.night('Callisto', '20250407')['method'] == 't1'
    reg.close()


def test_createlists_fails_clearly_without_a_register_or_a_bias(callisto_archive, monkeypatch):
    monkeypatch.setattr(cr, 'mount_fstype', lambda path, mounts_file=None: None)
    monkeypatch.setenv(cr.REGISTER_ENV, str(callisto_archive / 'nowhere.sqlite'))
    with pytest.raises(SystemExit) as e:
        run_t1(callisto_archive)
    assert 'no calibration register' in str(e.value)
    monkeypatch.delenv(cr.REGISTER_ENV)
    night = callisto_archive / 'Observations/Callisto/images/20250407'
    for p in night.glob('SPECU4.2025-04-08T10:48:*'):
        p.unlink()                                       # no bias left anywhere in the SPIRIT run
    with pytest.raises(SystemExit) as e:
        run_t1(callisto_archive)
    assert 'ERROR: T1 cannot calibrate Callisto 20250407' in str(e.value) and 'bias' in str(e.value)


def test_createlists_readonly_register_is_left_unchanged(callisto_archive, monkeypatch):
    monkeypatch.setattr(cr, 'mount_fstype', lambda path, mounts_file=None: None)
    monkeypatch.setenv(cr.READONLY_ENV, '1')
    red, lists = run_t1(callisto_archive)
    assert len(lists['1_bias.list']) == 5
    reg = cr.open_existing(cr.default_register_path(str(callisto_archive)), readonly=True)
    assert reg.night('Callisto', '20250407') is None
    reg.close()


# ------------------------------------------------------------------------------------- pipedark

def test_pipedark_groups_rounded_exposures_and_keeps_the_longest_as_fallback(tmp_path):
    pytest.importorskip('tabulate')                     # pipedark imports reporting.QC
    sys.path.insert(0, os.path.join(SRC, 'calibration'))  # and pipeutils as a top-level module
    from calibration import pipedark
    d = tmp_path / 'Observations/Callisto/images/20260601'
    out = tmp_path / 'red'
    out.mkdir()
    files = []
    for i, e in enumerate([42.407, 42.408, 42.407, 10.0, 10.0, 120.0, 120.0]):
        files.append(spirit(str(d / 'd{}.fits'.format(i)), 'Dark Frame', 18.5 * e, filt='Dark', exptime=e))
    fits.writeto(str(out / '1_MasterBias.fits'), np.full(SPIRIT, 2810.0, dtype=np.float32))
    lst = tmp_path / 'dark.list'
    lst.write_text('\n'.join(files) + '\n')
    made = pipedark.darkmaker(str(lst), '1_MasterBias.fits', '1_MasterDark.fits', str(out) + '/')
    names = sorted(os.path.basename(m) for m in made)
    assert names == ['1_MasterDark_10s.fits', '1_MasterDark_120s.fits', '1_MasterDark_42s.fits']
    with fits.open(str(out / '1_MasterDark_42s.fits')) as h:
        assert 'from 3 images' in str(h[0].header['HISTORY'])
    with fits.open(str(out / '1_MasterDark.fits')) as h:
        assert h[0].header['EXPTIME'] == 120
