"""Disk-side helpers: the cube-aware ESO diff, the FITS header reader, transformation_check, targets, products."""
import datetime as dt
import os

from jobqueue.nights import (artemis_targets, diff_eso_disk, products_state, read_header, read_transfer_log,
                             transformation_check_kind, transformation_check_would_delete)
from jobqueue.tests.fakes import write_fits
from jobqueue.util import last_night, mjd_to_night, night_bounds, to_mjd


def row(dp_id, dp_type='OBJECT', origfile=''):
    return {'dp_id': dp_id, 'dp_type': dp_type, 'origfile': origfile}


def test_andor_rows_match_dp_id_files():
    rows = [row('SPECU2.2026-10-08T07:58:06.000'), row('SPECU2.2026-10-08T07:58:26.000'),
            row('SPECU2.2026-10-08T07:58:46.000', 'BIAS')]
    names = ['SPECU2.2026-10-08T07:58:06.000.fits', 'SPECU2.2026-10-08T07:58:46.000.fits.fz']
    present, missing, _ = diff_eso_disk(rows, names)
    assert [r['dp_id'] for r in missing] == ['SPECU2.2026-10-08T07:58:26.000']
    assert len(present) == 2


def test_spirit_cubes_match_frames_by_time():
    # real shapes from Callisto 2026-10-07: the cube's dp_id time is not its last frame
    rows = [row('SPECU4.2026-10-08T00:32:25.987', origfile='SPECU4.20261007T235657_S_Sp2205-1104_zYJ_6s.fits'),
            row('SPECU4.2026-10-08T06:19:43.162', origfile='SPECU4.20261008T045509_S_Sp0218-0617_zYJ_15s.fits'),
            row('SPECU4.2026-10-08T09:45:21.859', 'FLAT', origfile='SPECU4.20261008T094421_C_flat_zYJ_60s.fits'),
            row('SPECU4.2026-10-08T10:05:01.392', 'BIAS', origfile='SPECU4.20261008T100500_C_bias.fits')]
    names = ['SPECU4.2026-10-07T23:56:57.087.fits', 'SPECU4.2026-10-08T01:10:00.000.fits',   # cube 1, after its dp_id
             'SPECU4.2026-10-08T04:55:09.500.fits',                                         # cube 2
             'SPECU4.2026-10-08T09:44:21.100.fits']                                         # the flat cube
    present, missing, counts = diff_eso_disk(rows, names)
    assert [r['dp_type'] for r in missing] == ['BIAS']
    assert counts['SPECU4.2026-10-08T00:32:25.987'] == 2


def test_header_reader_and_transformation_check(tmp_path):
    p = write_fits(str(tmp_path / 'SPECU2.2026-10-08T07:58:06.000.fits'), OBJECT="Sp0246+1625", IMAGETYP='Light Frame',
                   OBSERVER='Astra', ASTRAROT=90, EXPTIME=14.5, DATE_OBS='2026-10-08T07:58:06.000')
    h = read_header(p)
    assert h['OBJECT'] == 'Sp0246+1625' and h['ASTRAROT'] == 90 and h['EXPTIME'] == 14.5 and h['SIMPLE'] is True
    assert h['DATE-OBS'] == '2026-10-08T07:58:06.000'
    name = os.path.basename(p)
    assert not transformation_check_would_delete(name, h, 'ANDOR')
    assert transformation_check_would_delete(name, {'OBSERVER': 'speculoos'}, 'ANDOR')     # pre-Astra: deleted
    assert transformation_check_would_delete(name, h, 'SPIRIT')                            # needs ASTRAMIR
    assert transformation_check_would_delete('SPECU4.2019-04-15T01:00:00.000.fts', {}, 'ANDOR')
    assert not transformation_check_would_delete('Sp1056+0700-S001-R001-C001-i.fts', {}, 'ANDOR')


def test_which_nights_sso_download_checks():
    assert transformation_check_kind('Callisto', '20240910') == 'ANDOR'
    assert transformation_check_kind('Callisto', '20261007') == 'SPIRIT'
    assert transformation_check_kind('Callisto', '20221001') == 'SPIRIT'
    assert transformation_check_kind('Ganymede', '20250101') is None
    assert transformation_check_kind('Ganymede', '20251001') == 'ANDOR'
    assert transformation_check_kind('Europa', '20261008') is None


def test_artemis_targets_from_names():
    names = ['Sp0020+3305-S001-R001-C397-i.fts', 'Sp0449+5138-S002-R001-C001-i.fts',
             'Dark-S001-R004-C007-B1.fts', 'AutoFlat-Dusk-I+z-Bin1-007.fts', 'guider.log']
    assert artemis_targets(names) == ['Sp0020+3305', 'Sp0449+5138']


def test_transfer_log_and_products(tmp_path):
    base = str(tmp_path)
    os.makedirs(os.path.join(base, 'Observations', 'Callisto'))
    with open(os.path.join(base, 'Observations', 'Callisto', 'transfer_log.txt'), 'w') as f:
        f.write('20261007 19\n20261008 17\n')
    assert read_transfer_log(base, 'Callisto', '20261008') == 17
    assert read_transfer_log(base, 'Callisto', '20261009') is None
    assert read_transfer_log(base, 'Io', '20261008') is None
    assert products_state(base, 'Europa', '20261008') == 'none'
    v3 = os.path.join(base, 'PipelineOutput', 'v3', 'Europa', 'output', '20261008', 'Sp0001')
    os.makedirs(v3)
    open(os.path.join(v3, 'Sp0001_I+z_4_diff.fits'), 'w').close()
    assert products_state(base, 'Europa', '20261008') == 'v3_only'
    v2 = v3.replace('/v3/', '/v2/')
    os.makedirs(v2)
    open(os.path.join(v2, 'Sp0001_I+z_5_diff.fits'), 'w').close()
    assert products_state(base, 'Europa', '20261008') == 'v2'


def test_night_arithmetic():
    a, b = night_bounds('20261008')
    assert (a, b) == (dt.datetime(2026, 10, 8, 15), dt.datetime(2026, 10, 9, 15))
    assert str(mjd_to_night(to_mjd(dt.datetime(2026, 10, 9, 2, 11)))) == '2026-10-08'
    assert str(mjd_to_night(to_mjd(dt.datetime(2026, 10, 8, 16, 0)))) == '2026-10-08'
    assert str(last_night(dt.datetime(2026, 10, 9, 14, 0))) == '2026-10-08'
    assert str(last_night(dt.datetime(2026, 10, 10, 0, 30))) == '2026-10-08'
    assert str(last_night(dt.datetime(2026, 10, 9, 10, 59))) == '2026-10-07'
