"""The 30-night look-back and the add-only fetch."""
import datetime as dt
import json
import os

from jobqueue.fetch import add_only_fetch
from jobqueue.lookback import Lookback
from jobqueue.tests.fakes import PROGS, FakeEso, write_fits

NOW = dt.datetime(2026, 10, 9, 16, 30)


def europa_rows(night, n, science=True):
    base = dt.datetime.strptime(night, '%Y%m%d') + dt.timedelta(hours=23)
    return [{'dp_id': 'SPECU2.{}.000'.format((base + dt.timedelta(seconds=30 * k)).strftime('%Y-%m-%dT%H:%M:%S')),
             'dp_type': 'OBJECT' if science else 'BIAS', 'object': 'Sp0246+1625', 'origfile': 'x.fits',
             'mjd_obs': '0', 'night': night} for k in range(n)]


def put_on_disk(cfg, tel, night, rows):
    d = os.path.join(cfg['basedir'], 'Observations', tel, 'images', night)
    for r in rows:
        write_fits(os.path.join(d, r['dp_id'] + '.fits'), IMAGETYP='Light Frame', OBJECT=r['object'])


def only_europa(cfg):
    return dict(cfg, telescopes={'Europa': cfg['telescopes']['Europa']})


def make_lc(cfg, version, tel, night):
    d = os.path.join(cfg['basedir'], 'PipelineOutput', version, tel, 'output', night, 'Sp0246+1625')
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, 'Sp0246+1625_I+z_5_diff.fits'), 'w').close()


def scenario(cfg):
    eso = FakeEso()
    complete, partial, unproc, v3only = (europa_rows(n, 4) for n in ('20261005', '20261004', '20261003', '20261002'))
    eso.rows[PROGS['Europa']] = complete + partial + unproc + v3only
    put_on_disk(cfg, 'Europa', '20261005', complete)
    put_on_disk(cfg, 'Europa', '20261004', partial[:2])
    put_on_disk(cfg, 'Europa', '20261003', unproc)
    put_on_disk(cfg, 'Europa', '20261002', v3only)
    make_lc(cfg, 'v2', 'Europa', '20261005')
    make_lc(cfg, 'v2', 'Europa', '20261004')
    make_lc(cfg, 'v3', 'Europa', '20261002')
    return eso, partial


def by_night(results):
    return {r['night']: r for r in results}


def test_lookback_fetches_missing_frames_and_reruns_nights_without_products(cfg, store):
    cfg = only_europa(cfg)
    eso, partial = scenario(cfg)
    res = by_night(Lookback(cfg, store, eso, now=NOW, log=lambda m: None).run())
    assert eso.calls[0] == ('frames', PROGS['Europa'], '20260908', '20261007')     # 30 nights, not last night
    assert res['20261005']['action'] == 'ok'
    assert res['20261004']['action'].startswith('fetch 2 rows (2 science)')
    assert res['20261003']['action'].startswith('rerun (no light curves')
    assert 'T12 promotion' in res['20261002']['action']
    fetch = store.jobs(kind='fetch')[0]
    assert fetch['class'] == 'P1' and fetch['night'] == '20261004' and 'sem:eso' in fetch['locks']
    rows_file = fetch['argv'][fetch['argv'].index('--rows') + 1]
    assert [r['dp_id'] for r in json.load(open(rows_file))] == [r['dp_id'] for r in partial[2:]]
    rerun = store.jobs(kind='pipeline')[0]
    assert rerun['class'] == 'P1' and rerun['night'] == '20261003' and 'target:SP0246+1625' in rerun['locks']


def test_lookback_leaves_busy_nights_alone_and_limits_reruns(cfg, store):
    cfg = only_europa(cfg)
    eso, _ = scenario(cfg)
    Lookback(cfg, store, eso, now=NOW, log=lambda m: None).run()
    res = by_night(Lookback(cfg, store, eso, now=NOW, log=lambda m: None).run())
    assert res['20261003']['action'].startswith('skipped')            # its rerun is still queued
    rerun = store.jobs(kind='pipeline')[0]
    for _ in range(2):
        store.mark_running(rerun['id'])
        store.finish(rerun['id'], 1, signature='boom')
        res = by_night(Lookback(cfg, store, eso, now=NOW, log=lambda m: None).run())
        rerun = (store.jobs(kind='pipeline', night='20261003', states=['queued']) or [rerun])[0]
    assert 'automatic reruns already tried' in res['20261003']['action']


def test_lookback_dry_run_records_nothing(cfg, store):
    cfg = only_europa(cfg)
    eso, _ = scenario(cfg)
    res = by_night(Lookback(cfg, store, eso, now=NOW, log=lambda m: None, dry_run=True).run())
    assert res['20261004']['action'].startswith('would fetch 2 rows')
    assert res['20261003']['action'].startswith('would rerun')
    assert store.jobs() == [] and store.lookback_rows() == []


def test_lookback_survives_an_eso_failure(cfg, store):
    eso = FakeEso()
    eso.fail = True
    res = Lookback(only_europa(cfg), store, eso, now=NOW, log=lambda m: None).run()
    assert 'ESO query failed' in res[0]['action'] and store.jobs() == []


# ---------------------------------------------------------------- add-only fetch


def test_fetch_only_adds_and_never_overwrites(tmp_path):
    night = tmp_path / 'images' / '20261004'
    night.mkdir(parents=True)
    existing = night / 'SPECU2.A.fits'
    existing.write_bytes(b'original')
    staging = tmp_path / 'staging'

    def download(dp_ids, dest):
        out = []
        for d in dp_ids:
            if d == 'SPECU2.LOST':
                continue
            p = os.path.join(dest, d + '.fits')
            write_fits(p, IMAGETYP='Light Frame' if d != 'SPECU2.CUBE' else 'Light Frame', OBSERVER='Astra')
            out.append(p)
        return out

    def is_cube(path):
        return path.endswith('SPECU2.CUBE.fits')

    def unpack(path, dest):
        frames = [write_fits(os.path.join(dest, 'SPECU2.F{}.fits'.format(k)), IMAGETYP='Light Frame',
                             OBSERVER='Astra') for k in (1, 2)]
        os.remove(path)
        return frames

    transformed = []
    rows = [{'dp_id': d} for d in ('SPECU2.A', 'SPECU2.B', 'SPECU2.CUBE', 'SPECU2.LOST')]
    res = add_only_fetch(rows, str(night), str(staging), download, unpack=unpack, is_cube=is_cube,
                         transform=lambda ids, d: transformed.extend(ids), log=lambda m: None)
    assert existing.read_bytes() == b'original'                         # never overwritten
    assert sorted(os.listdir(night)) == ['SPECU2.A.fits', 'SPECU2.B.fits', 'SPECU2.F1.fits', 'SPECU2.F2.fits']
    assert res['added'] == 3 and res['added_science'] == 3 and res['already_present'] == 1
    assert res['not_downloaded'] == ['SPECU2.LOST'] and res['settled'] == ['SPECU2.A']
    assert sorted(transformed) == ['SPECU2.A', 'SPECU2.B', 'SPECU2.F1', 'SPECU2.F2']   # Astra frames only
    assert not staging.exists()


def test_fetch_creates_no_empty_night_directory(tmp_path):
    night = tmp_path / 'images' / '20261001'
    res = add_only_fetch([{'dp_id': 'SPECU2.X'}], str(night), str(tmp_path / 'st'), lambda ids, d: [],
                         log=lambda m: None)
    assert res['added'] == 0 and not night.exists()


def test_fetch_skips_the_astra_transform_for_pre_astra_frames(tmp_path):
    def download(dp_ids, dest):
        return [write_fits(os.path.join(dest, d + '.fits'), IMAGETYP='Light Frame', OBSERVER='speculoos')
                for d in dp_ids]
    transformed = []
    add_only_fetch([{'dp_id': 'SPECU4.OLD'}], str(tmp_path / 'n'), str(tmp_path / 's'), download,
                   transform=lambda ids, d: transformed.extend(ids), log=lambda m: None)
    assert transformed == [] and (tmp_path / 'n' / 'SPECU4.OLD.fits').exists()


def test_a_cube_that_fails_to_unpack_is_not_settled(tmp_path):
    def download(dp_ids, dest):
        return [write_fits(os.path.join(dest, d + '.fits'), IMAGETYP='Light Frame') for d in dp_ids]
    res = add_only_fetch([{'dp_id': 'SPECU4.BADCUBE'}], str(tmp_path / 'n'), str(tmp_path / 's'), download,
                         is_cube=lambda p: True, unpack=lambda p, d: [], log=lambda m: None)
    assert res['unpack_failed'] == ['SPECU4.BADCUBE'] and res['settled'] == [] and res['added'] == 0
    assert not (tmp_path / 'n').exists()
