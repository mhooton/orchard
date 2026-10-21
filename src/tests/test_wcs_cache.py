"""
T7's WCS cache (astrom/wcs_cache.py, docs/wcs-cache.md), and T7 no longer writing raw frames.

T7 used to copy every solution back into the raw frame as a cache, but SSO always forces a re-solve (the
"already solved" test could not tell Astra's WCS from orchard's), so ~379,000 raw-frame writes in 30 days were never
read. Now nothing in T7 writes a raw frame; a night's solutions go to
<BASEDIR>/PipelineOutput/wcs_cache/<TEL>/<night>.ecsv and a later run applies them, header only, when the raw frame,
the calibrated frame's shape, the solver, the Gaia database and the libraries are unchanged and the row came from a
release or the running commit.

The solves here are real pointer_wcs + twirl solves; only the image cleaning, source detection and the Gaia query
are faked (the pattern of test_sexagesimal_ra.py): the catalogue is the detections projected through a known WCS.

Run:  pytest src/tests/test_wcs_cache.py -q
"""
import os
import re
import shutil
import subprocess
import sys
from argparse import Namespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pw = pytest.importorskip('astrom.pointer_wcs')
from astrom import astrometry, header_write, wcs_cache   # noqa: E402
from calibration.pipeutils import proc_name              # noqa: E402
from download import raw_provenance                      # noqa: E402
from utils import orchard_version as ov                  # noqa: E402

NIGHT = '20261018'
COMMIT_A = 'a' * 40
COMMIT_B = 'b' * 40
SHAPE = (200, 200)
# where the fake detector finds stars, (x, y) in pixels
STARS = np.array([[20.0 + 13 * k, 25.0 + (37 * k) % 150] for k in range(12)])
SCALE_DEG = 0.3480718605416447 / 3600       # the plate scale pointer_wcs works out for these headers


def _code(monkeypatch, commit=COMMIT_A, release='', version='4.0.2'):
    monkeypatch.setattr(ov, 'info', lambda: {'version': version, 'commit': commit, 'release': release,
                                            'date': None, 'source': 'archive'})


def _truth(center):
    """The WCS a frame really has: pointed 20 arcsec off the header position, rotated by 1 degree."""
    ra, dec = center
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [ra + 20 / 3600 / np.cos(np.radians(dec)), dec - 20 / 3600]
    w.wcs.crpix = [100.5, 100.5]
    c, s = np.cos(np.radians(1.0)), np.sin(np.radians(1.0))
    w.wcs.cd = SCALE_DEG * np.array([[-c, s], [s, c]])
    return w


@pytest.fixture
def solver(monkeypatch):
    """Real pointer_wcs solves on synthetic frames. Returns the list of frame centres the Gaia query was asked for
    (one per solve attempted in this process)."""
    asked = []

    def fake_query(center, **kwargs):
        asked.append(center)
        sky = _truth(center).pixel_to_world(STARS[:, 0], STARS[:, 1])
        return np.column_stack([sky.ra.deg, sky.dec.deg])

    monkeypatch.setattr(pw, 'clean_image', lambda data: data)
    monkeypatch.setattr(pw, 'find_stars_multiscale', lambda *a, **k: STARS.copy())
    monkeypatch.setattr(pw, 'gaia_db_query', fake_query)
    _code(monkeypatch)
    return asked


def _raw_header(k, astra_wcs=True):
    h = fits.Header()
    h.update(TELESCOP='SPECULOOS-IO', IMAGETYP='Light Frame', RA=313.6875 + 0.001 * k, DEC=9.1747 - 0.001 * k,
             XPIXSZ=13.5, FOCALLEN=8000.0, EXPTIME=30.0, FILTER='I+z')
    h['DATE-OBS'] = '2026-10-18T23:{:02d}:00'.format(k)
    if astra_wcs:       # what an SSO frame arrives with: Astra's on-site solution, which T7 must not trust
        h.update(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CRVAL1=313.0, CRVAL2=9.0, CRPIX1=1024.0, CRPIX2=1024.0,
                 CD1_1=-9.7e-5, CD1_2=0.0, CD2_1=0.0, CD2_2=9.7e-5)
    return h


class Night:
    """A data root with one night of raw frames, T1's raw list, and T6's calibrated frames and list. With
    `checksums`, the raw frames carry CHECKSUM/DATASUM as ESO's do, and the calibrated frames are T6's: new pixels in
    the opened raw frame, written by `t6_writer`. That is header_write.writeto as T6 writes now, or by default plain
    writeto() as T6 wrote before, which keeps the raw frame's cards (stale)."""

    def __init__(self, root, n=3, names=None, checksums=False, t6_writer=None):
        self.root = root
        self.checksums = checksums
        self.t6_writer = t6_writer or (lambda hdulist, path: hdulist.writeto(path))
        self.rawdir = root / 'Observations' / 'Io' / 'images' / NIGHT
        self.outdir = root / 'PipelineOutput' / 'v3' / 'Io' / 'output' / NIGHT
        self.procdir = self.outdir / 'T' / '1'
        self.rawdir.mkdir(parents=True)
        self.procdir.mkdir(parents=True)
        (self.outdir / 'reduction').mkdir()
        names = names or ['Io_{:03d}.fits'.format(k) for k in range(n)]
        self.raws = []
        for k, name in enumerate(names):
            path = self.rawdir / name
            data = (1000 + np.arange(np.prod(SHAPE)) % 97).reshape(SHAPE).astype(np.uint16)
            fits.PrimaryHDU(data=data, header=_raw_header(k)).writeto(str(path), checksum=checksums)
            self.raws.append(str(path))
        self.raw_list = self.outdir / 'reduction' / '1_image_T.list'
        self.raw_list.write_text('\n'.join(self.raws) + '\n')
        self.db = root / 'gaia.db'
        self.db.write_bytes(b'not really a database')
        self.cache = root / 'PipelineOutput' / 'wcs_cache' / 'Io' / '{}.ecsv'.format(NIGHT)
        self.reduce()

    def reduce(self):
        """What T6 does to the WCS: a calibrated copy of each raw frame that inherits its header."""
        self.procs = []
        for raw in self.raws:
            out = str(self.procdir / proc_name(raw))
            if self.checksums:          # calibration/pipered.py: new pixels in the opened raw frame, written
                if os.path.exists(out):
                    os.remove(out)
                with fits.open(raw) as h:
                    h[0].data = h[0].data.astype(np.float32)
                    self.t6_writer(h, out)
                self.procs.append(out)
                continue
            with fits.open(raw) as h:
                header = h[0].header.copy()
                for key in ('BZERO', 'BSCALE'):
                    header.remove(key, ignore_missing=True)
                fits.PrimaryHDU(data=h[0].data.astype(np.float32), header=header).writeto(out, overwrite=True)
            self.procs.append(out)
        self.filelist = self.procdir / 'I+z_processed.dat'
        self.filelist.write_text('\n'.join(self.procs))
        return self.procs

    def args(self, **kw):
        a = dict(filelist=str(self.filelist), nproc=1, ext='fits', db_path=str(self.db), raw_list=str(self.raw_list),
                 wcs_cache=str(self.cache), force_platesolve=False, timeout=60)
        a.update(kw)
        return Namespace(**a)

    def t7(self, **kw):
        astrometry.main(self.args(**kw))

    def rows(self):
        return wcs_cache.load(str(self.cache))


def _cards(path):
    """Every card of a frame's header but HISTORY, as written: (keyword, image)."""
    h = fits.getheader(path)
    return [(c.keyword, c.image) for c in h.cards if c.keyword != 'HISTORY']


def _history(path):
    return [str(v) for v in fits.getheader(path).get('HISTORY', [])]


def _raw_state(night):
    out = {}
    for raw in night.raws:
        with open(raw, 'rb') as f:
            out[raw] = (f.read(), os.stat(raw).st_mtime_ns)
    return out


# ------------------------------------------------------------------------------------------------ name map


@pytest.mark.parametrize('raw, proc', [
    ('/x/Io_001.fits', 'procIo_001.fits'),
    ("/x/Sp1056+0700 'I+z' 001.fts", 'procSp1056+0700--I+z--001.fits'),
    ('/x/IMAGE123.fits.bz2', 'procIMAGE123.fits.bz2'),
])
def test_proc_name_is_t6s_rule(raw, proc):
    assert proc_name(raw) == proc


def test_raw_frames_are_found_by_t6s_rule(tmp_path):
    night = Night(tmp_path, names=["Sp1056+0700 'I+z' 001.fts", 'Io_002.fits'])
    by_proc = astrometry.raw_frames_by_proc(str(night.raw_list))
    assert by_proc == {os.path.basename(p): r for p, r in zip(night.procs, night.raws)}
    assert astrometry.raw_frames_by_proc(str(tmp_path / 'missing.list')) == {}


# ------------------------------------------------------------------------------------------------ T7 runs


def test_first_run_solves_every_frame_and_fills_the_cache(tmp_path, solver):
    night = Night(tmp_path)
    night.t7()
    assert len(solver) == 3
    rows = night.rows()
    assert sorted(rows) == sorted(os.path.basename(r) for r in night.raws)
    for raw, proc in zip(night.raws, night.procs):
        row = rows[os.path.basename(raw)]
        h = fits.getheader(proc)
        assert 'CD1_1' not in h and h['CTYPE1'] == 'RA---TAN'           # Astra's WCS gone, orchard's in its place
        assert 'orchard 4.0.2 T7: WCS solved' in _history(proc)
        st = os.stat(raw)
        assert (row['raw_size'], row['raw_mtime_ns']) == (st.st_size, st.st_mtime_ns)
        assert (row['naxis1'], row['naxis2']) == (200, 200)
        assert (row['trim_x'], row['trim_y']) == (2, 22)
        assert row['solver_epoch'] == astrometry.SOLVER_EPOCH
        assert row['gaia_db'] == str(night.db) and row['gaia_db_size'] == night.db.stat().st_size
        assert (row['orchver'], row['orchrel'], row['commit']) == ('4.0.2', '', COMMIT_A)
        assert row['matches'] >= 10 and row['proc'] == proc
        assert row['cards'] == wcs_cache.cards_from_header(h, astrometry.WCS_KEYWORDS)
        assert 'MJDREF' in wcs_cache.cards_header(row['cards'])


def test_second_run_applies_the_cache_with_the_cards_of_a_fresh_solve(tmp_path, solver):
    night = Night(tmp_path)
    night.t7()
    fresh = {p: _cards(p) for p in night.procs}
    night.reduce()                                  # a rerun recalibrates every frame from its raw frame
    del solver[:]
    night.t7()
    assert solver == []                             # nothing solved
    for proc in night.procs:
        assert _cards(proc) == fresh[proc]          # every card, in order, as the fresh solve wrote it
        hist = _history(proc)
        assert 'orchard 4.0.2 T7: WCS from cache' in hist and 'orchard 4.0.2 T7: WCS solved' not in hist
        assert 'orchard 4.0.2 T7: inherited WCS removed' in hist


def test_raw_frames_are_never_written(tmp_path, solver):
    night = Night(tmp_path)
    before = _raw_state(night)
    night.t7()                                      # solves
    night.reduce()
    night.t7()                                      # from the cache
    night.reduce()
    night.t7(force_platesolve=True)                 # forced
    assert _raw_state(night) == before


def test_parent_collects_the_pool_results(tmp_path, solver, monkeypatch):
    # the workers fork (astrom/frame_pool.py), so they inherit the fake catalogue
    stored = []
    real_store = wcs_cache.store
    monkeypatch.setattr(wcs_cache, 'store', lambda path, new: stored.append(sorted(new)) or real_store(path, new))
    night = Night(tmp_path, n=4)
    night.t7(nproc=2)
    assert stored == [sorted(os.path.basename(r) for r in night.raws)]   # one write, by the parent, of every solve
    rows = night.rows()
    for raw, proc in zip(night.raws, night.procs):
        assert rows[os.path.basename(raw)]['cards'] == wcs_cache.cards_from_header(fits.getheader(proc),
                                                                                    astrometry.WCS_KEYWORDS)
    night.reduce()
    night.t7(nproc=2)
    assert len(stored) == 1                         # second run: every frame from the cache, nothing to store


def test_force_ignores_the_cache_and_refreshes_its_rows(tmp_path, solver, monkeypatch):
    night = Night(tmp_path)
    night.t7()
    night.reduce()
    del solver[:]
    _code(monkeypatch, version='4.0.3')
    night.t7(force_platesolve=True)
    assert len(solver) == 3
    assert {row['orchver'] for row in night.rows().values()} == {'4.0.3'}
    for proc in night.procs:
        assert 'orchard 4.0.3 T7: WCS solved' in _history(proc)


def test_failures_and_timeouts_are_not_cached(tmp_path, solver, monkeypatch):
    night = Night(tmp_path)
    good = pw.find_stars_multiscale
    bad = os.path.basename(night.procs[1])
    monkeypatch.setattr(pw, 'find_stars_multiscale',
                        lambda data, *a, **k: STARS[:3].copy() if data[0, 0] == -1 else good(data, *a, **k))
    with fits.open(night.procs[1], mode='update') as h:
        h[0].data[0, 0] = -1                                    # this frame's detections are too few to solve
    night.t7()
    rows = night.rows()
    assert sorted(rows) == sorted(os.path.basename(r) for r in night.raws if proc_name(r) != bad)

    def timeout(path, *a, **k):                                 # what the SIGALRM timeout returns
        return astrometry._partial_result({'stage': 'gaia_query'}, 'Plate solving timed out after 60s')
    monkeypatch.setattr(astrometry, 'pointer_wcs', timeout)
    night.reduce()
    night.t7(force_platesolve=True)
    assert night.rows() == rows                                 # nothing refreshed, nothing added


def test_a_header_wcs_is_not_a_solution(tmp_path, solver):
    """Astra's WCS (or PinPoint's, or orchard's copied into old raw frames) used to make T7 skip a frame unless
    forced. Without a cache row the frame is now always solved."""
    night = Night(tmp_path, n=2)
    night.t7(wcs_cache=None)
    assert len(solver) == 2
    assert not night.cache.exists()
    for proc in night.procs:
        h = fits.getheader(proc)
        assert 'CD1_1' not in h and h['CRVAL1'] != 313.0


def test_an_empty_list_is_not_an_error(tmp_path, solver):
    night = Night(tmp_path, n=1)
    night.filelist.write_text('')
    night.t7()
    assert solver == [] and not night.cache.exists()


# ------------------------------------------------------------------------------------------------ validity


def _fields():
    return {'raw_size': 1000, 'raw_mtime_ns': 1760000000123456789, 'naxis1': 2046, 'naxis2': 2044, 'trim_x': 2,
            'trim_y': 22, 'solver_epoch': 1, 'gaia_db': '/gaia_database/gaia_tmass_16_jm_cut.db',
            'gaia_db_size': 123456789, 'gaia_db_mtime_ns': 1750000000000000000, 'astropy': '5.3.4', 'twirl': '0.4.2'}


def _changed(value):
    return value + 1 if isinstance(value, int) else value + 'x'


@pytest.mark.parametrize('field', wcs_cache.VALIDITY)
def test_each_validity_field_invalidates_a_row(field, monkeypatch):
    _code(monkeypatch)
    row = wcs_cache.make_entry('Io_001.fits', _fields(), 'CTYPE1  = ...')
    assert wcs_cache.check(row, _fields(), wcs_cache.code_identity()) is None
    expected = dict(_fields(), **{field: _changed(_fields()[field])})
    assert wcs_cache.check(row, expected, wcs_cache.code_identity()) == field


def test_validity_fields_survive_the_table(tmp_path, monkeypatch):
    _code(monkeypatch)
    path = str(tmp_path / 'c.ecsv')
    row = wcs_cache.make_entry('Io 001.fits', _fields(), 'CTYPE1  = ...', matches=None, proc='/p/procIo--001.fits')
    wcs_cache.store(path, {row['raw']: row})
    back = wcs_cache.load(path)['Io 001.fits']
    assert back == row
    assert wcs_cache.check(back, _fields(), wcs_cache.code_identity()) is None


@pytest.mark.parametrize('change', ['touch raw', 'gaia db', 'solver epoch', 'twirl', 'shape'])
def test_a_changed_input_solves_again(tmp_path, solver, monkeypatch, change):
    night = Night(tmp_path, n=2)
    night.t7()
    if change == 'touch raw':                       # a repair or re-download rewrites the raw frame
        st = os.stat(night.raws[0])
        os.utime(night.raws[0], ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
    elif change == 'gaia db':
        night.db.write_bytes(b'a rebuilt database')
    elif change == 'solver epoch':
        monkeypatch.setattr(astrometry, 'SOLVER_EPOCH', astrometry.SOLVER_EPOCH + 1)
    elif change == 'twirl':
        real = wcs_cache._package_version
        monkeypatch.setattr(wcs_cache, '_package_version', lambda name: '9.9' if name == 'twirl' else real(name))
    night.reduce()
    if change == 'shape':                           # e.g. a binned frame of the same name
        with fits.open(night.procs[0]) as h:
            fits.PrimaryHDU(data=h[0].data[:100], header=h[0].header).writeto(night.procs[0], overwrite=True)
    del solver[:]
    night.t7()
    assert len(solver) == (1 if change in ('touch raw', 'shape') else 2)
    assert 'orchard 4.0.2 T7: WCS solved' in _history(night.procs[0])


# ------------------------------------------------------------------------------------------------ provenance


@pytest.mark.parametrize('writer, reader, usable', [
    (dict(release='20261019T120000Z-abc1234', commit=COMMIT_B), dict(commit=COMMIT_A), True),   # a release
    (dict(release='', commit=COMMIT_A), dict(commit=COMMIT_A), True),                            # same commit
    (dict(release='', commit=COMMIT_B), dict(commit=COMMIT_A), False),                           # a scratch tree
    (dict(release='', commit=COMMIT_B), dict(commit=COMMIT_A, release='20261019T1-a'), False),   # ... into production
    (dict(release='', commit=COMMIT_A, version='4.0.2-dirty'), dict(commit=COMMIT_A), False),   # uncommitted code
    (dict(release='', commit=''), dict(commit=''), False),                                      # no known commit
])
def test_only_a_release_or_the_running_commit_feeds_a_run(monkeypatch, writer, reader, usable):
    _code(monkeypatch, **writer)
    row = wcs_cache.make_entry('Io_001.fits', _fields(), 'CTYPE1  = ...')
    _code(monkeypatch, **reader)
    assert wcs_cache.check(row, _fields(), wcs_cache.code_identity()) == (None if usable else 'provenance')


def test_a_scratch_trees_rows_do_not_feed_another_commit(tmp_path, solver, monkeypatch):
    night = Night(tmp_path, n=2)
    _code(monkeypatch, commit=COMMIT_B)             # a scratch tree run against this data root
    night.t7()
    _code(monkeypatch, commit=COMMIT_A, release='20261020T000000Z-aaaaaaa')
    night.reduce()
    del solver[:]
    night.t7()
    assert len(solver) == 2                         # production solved for itself ...
    assert {r['orchrel'] for r in night.rows().values()} == {'20261020T000000Z-aaaaaaa'}
    _code(monkeypatch, commit=COMMIT_B)
    night.reduce()
    del solver[:]
    night.t7()
    assert solver == []                             # ... and its rows feed anyone


# ------------------------------------------------------------------------------------------------ damage


def test_a_corrupt_cache_falls_back_to_solving(tmp_path, solver, capsys):
    night = Night(tmp_path)
    night.cache.parent.mkdir(parents=True)
    night.cache.write_text('# %ECSV 1.0\nthis is not a table\n"unbalanced\n')
    night.t7()
    assert len(solver) == 3
    assert 'WARNING: WCS cache not used' in capsys.readouterr().out
    assert len(night.rows()) == 3                                      # replaced by a good table
    assert (night.cache.parent / (night.cache.name + '.corrupt')).exists()


def test_a_bad_row_falls_back_to_solving_that_frame(tmp_path, solver, capsys):
    night = Night(tmp_path)
    night.t7()
    rows = night.rows()
    broken = os.path.basename(night.raws[0])
    rows[broken]['cards'] = rows[broken]['cards'][:160]                # CTYPEs lost
    shutil.move(str(night.cache), str(night.cache) + '.good')
    wcs_cache.store(str(night.cache), rows)
    night.reduce()
    del solver[:]
    night.t7()
    assert len(solver) == 1
    assert 'WARNING: WCS cache not used for ' + os.path.basename(night.procs[0]) in capsys.readouterr().out
    assert 'orchard 4.0.2 T7: WCS solved' in _history(night.procs[0])
    assert night.rows()[broken]['cards'] == wcs_cache.cards_from_header(fits.getheader(night.procs[0]),
                                                                         astrometry.WCS_KEYWORDS)


def test_a_failed_header_write_falls_back_to_solving_that_frame(tmp_path, solver, capsys, monkeypatch):
    night = Night(tmp_path)
    night.t7()
    night.reduce()
    real = wcs_cache.apply
    broken = night.procs[1]

    def apply(path, cards, keywords):
        if path == broken:
            raise OSError('disk went away')
        return real(path, cards, keywords)
    monkeypatch.setattr(wcs_cache, 'apply', apply)
    del solver[:]
    night.t7()
    assert len(solver) == 1
    out = capsys.readouterr().out
    assert 'WARNING: WCS cache not used for {}: disk went away'.format(os.path.basename(broken)) in out
    assert 'WCS cache: 2 frames from the cache, 1 to solve (unusable 1)' in out
    assert 'orchard 4.0.2 T7: WCS solved' in _history(broken)
    assert all('orchard 4.0.2 T7: WCS from cache' in _history(p) for p in night.procs if p != broken)


def test_an_unwritable_cache_does_not_stop_t7(tmp_path, solver, capsys):
    night = Night(tmp_path, n=2)
    night.t7(wcs_cache=str(night.raws[0]) + '/cannot/be/here.ecsv')   # a file where a directory should be
    assert len(solver) == 2
    assert 'WARNING: WCS cache not updated' in capsys.readouterr().out
    for proc in night.procs:
        assert 'orchard 4.0.2 T7: WCS solved' in _history(proc)


def test_a_missing_gaia_database_turns_the_cache_off(tmp_path, solver, capsys):
    night = Night(tmp_path, n=1)
    night.t7(db_path=str(tmp_path / 'nowhere.db'))
    assert len(solver) == 1 and not night.cache.exists()
    assert 'WARNING: WCS cache not used' in capsys.readouterr().out


# ------------------------------------------------------------------------------------------------ the header write


def _verifies(path):
    with fits.open(path, do_not_scale_image_data=True) as h:
        return h[0].verify_datasum() == 1 and h[0].verify_checksum() == 1


def _data_bytes(path):
    raw, _ = raw_provenance.read_header_bytes(path)
    with open(path, 'rb') as f:
        return f.read()[len(raw):]


def _no_data_sum(*args):
    raise AssertionError('the data were summed although the DATASUM could be trusted')


def test_checksums_verify_after_a_solve_and_after_a_cache_hit(tmp_path, solver, monkeypatch):
    """Frames T6 calibrated before it set CHECKSUM/DATASUM carry the raw frame's, stale. T7 sums the data once and
    leaves both valid, whether it solves a frame or applies its cached solution; a frame T7 has written is trusted
    after that."""
    night = Night(tmp_path, checksums=True)
    assert not any(_verifies(p) for p in night.procs)
    pixels = {p: _data_bytes(p) for p in night.procs}
    night.t7()                                      # solves
    assert len(solver) == 3
    assert all(_verifies(p) and _data_bytes(p) == pixels[p] for p in night.procs)
    night.reduce()
    del solver[:]
    night.t7()                                      # from the cache
    assert solver == []
    assert all(_verifies(p) and _data_bytes(p) == pixels[p] for p in night.procs)
    assert all('orchard 4.0.2 T7: WCS from cache' in _history(p) for p in night.procs)

    monkeypatch.setattr(header_write, 'data_sum', _no_data_sum)
    for run in ({}, {'force_platesolve': True}):    # from the cache again, then solved, without T6 in between
        headers = {p: raw_provenance.read_header_bytes(p)[0] for p in night.procs}
        night.t7(**run)
        for p in night.procs:
            assert raw_provenance.read_header_bytes(p)[0] != headers[p]         # written
            assert _verifies(p) and _data_bytes(p) == pixels[p]
    assert len(solver) == 3


def test_t7_never_sums_the_data_of_frames_t6_wrote(tmp_path, solver, monkeypatch):
    """T6 now writes CHECKSUM/DATASUM that verify, so every T7 write, a solve's or a cache hit's, trusts them."""
    night = Night(tmp_path, checksums=True, t6_writer=header_write.writeto)
    assert all(_verifies(p) for p in night.procs)
    pixels = {p: _data_bytes(p) for p in night.procs}
    monkeypatch.setattr(header_write, 'data_sum', _no_data_sum)
    night.t7()                                      # solves
    night.reduce()
    night.t7()                                      # from the cache
    assert len(solver) == 3
    assert all(_verifies(p) and _data_bytes(p) == pixels[p] for p in night.procs)
    assert all('orchard 4.0.2 T7: WCS from cache' in _history(p) for p in night.procs)


def test_applying_a_row_keeps_a_16_bit_frame_16_bit(tmp_path, solver):
    night = Night(tmp_path, n=1, checksums=True)
    night.t7()
    row = night.rows()[os.path.basename(night.raws[0])]
    f = str(tmp_path / 'int16.fits')
    shutil.copy(night.raws[0], f)                                       # uint16 data: BITPIX 16, BZERO 32768
    before, pixels = fits.getdata(f), _data_bytes(f)
    wcs_cache.apply(f, row['cards'], astrometry.WCS_KEYWORDS)
    h = fits.getheader(f)
    assert (h['BITPIX'], h['BZERO']) == (16, 32768)
    np.testing.assert_array_equal(fits.getdata(f), before)
    assert _data_bytes(f) == pixels and _verifies(f)
    assert WCS(h).has_celestial and 'CD1_1' not in h


def test_a_cache_hit_on_a_stale_checksum_sums_the_data(tmp_path, solver, monkeypatch):
    night = Night(tmp_path, n=1, checksums=True)
    night.t7()
    night.reduce()                                  # T6 again, as it wrote before: the raw frame's cards, stale
    calls = []
    monkeypatch.setattr(header_write, 'data_sum', lambda *a: calls.append(a) or raw_provenance.data_sum(*a))
    row = night.rows()[os.path.basename(night.raws[0])]
    wcs_cache.apply(night.procs[0], row['cards'], astrometry.WCS_KEYWORDS)
    assert len(calls) == 1 and _verifies(night.procs[0])


def test_a_row_that_does_not_parse_leaves_the_frame_untouched(tmp_path):
    f = str(tmp_path / 'f.fits')
    fits.PrimaryHDU(data=np.zeros((4, 4), np.float32), header=_raw_header(0)).writeto(f)
    before = open(f, 'rb').read()
    with pytest.raises(wcs_cache.CacheError):
        wcs_cache.apply(f, 'CRVAL1  = 1.0', astrometry.WCS_KEYWORDS)
    assert open(f, 'rb').read() == before


# ------------------------------------------------------------------------------------------------ writers


WRITERS = 6
ROWS_PER_WRITER = 15


def _env():
    env = dict(os.environ)
    env['PYTHONPATH'] = SRC + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    return env


def test_concurrent_writers_lose_no_rows(tmp_path):
    """Two targets of one night (or a rerun beside a nightly run) finishing T7 at once: one read-change-write
    each, under the table's lock."""
    path = str(tmp_path / 'wcs_cache' / 'Io' / (NIGHT + '.ecsv'))
    script = """
import sys
from astrom import wcs_cache
fields = {{f: 1 for f in wcs_cache.VALIDITY}}
fields.update(gaia_db='/g.db', astropy='5.3.4', twirl='0.4.2')
for k in range({n}):
    name = 'w{w}_{{:03d}}.fits'.format(k)
    wcs_cache.store(sys.argv[1], {{name: wcs_cache.make_entry(name, fields, 'CTYPE1  = x' * 20)}})
"""
    procs = [subprocess.Popen([sys.executable, '-c', script.format(w=w, n=ROWS_PER_WRITER), path], env=_env(),
                              stderr=subprocess.PIPE) for w in range(WRITERS)]
    for p in procs:
        _, err = p.communicate(timeout=300)
        assert p.returncode == 0, err.decode()
    rows = wcs_cache.load(path)
    assert sorted(rows) == sorted('w{}_{:03d}.fits'.format(w, k) for w in range(WRITERS)
                                  for k in range(ROWS_PER_WRITER))
    assert not [f for f in os.listdir(os.path.dirname(path)) if '.tmp' in f]


# ------------------------------------------------------------------------------------------------ ZLP


ZLP = os.path.join(SRC, 'main', 'ZLP_pipeline.sh')


def _zlp_cache_path(env):
    line = [ln.strip() for ln in open(ZLP) if ln.strip().startswith('WCS_CACHE=')]
    assert len(line) == 1
    script = 'set -o nounset; BASEDIR=/data/SPECULOOSPipeline; TEL=Io; DATE={}; {}; echo "$WCS_CACHE"'.format(
        NIGHT, line[0])
    full = {k: v for k, v in os.environ.items() if k != 'WCS_CACHE_DIR'}
    full.update(env)
    return subprocess.run(['bash', '-c', script], env=full, stdout=subprocess.PIPE, check=True,
                          universal_newlines=True).stdout.strip()


def test_zlp_derives_the_cache_path():
    assert _zlp_cache_path({}) == '/data/SPECULOOSPipeline/PipelineOutput/wcs_cache/Io/{}.ecsv'.format(NIGHT)
    assert _zlp_cache_path({'WCS_CACHE_DIR': '/data/x/cache'}) == '/data/x/cache/Io/{}.ecsv'.format(NIGHT)


def test_zlp_passes_the_cache_to_t7_and_its_own_command_line_is_unchanged():
    text = open(ZLP).read()
    start = text.index('check_astrometry(){')
    t7 = text[start:text.index('\n}\n', start)]
    assert '--wcs-cache ${WCS_CACHE}' in t7 and '--raw-list ' in t7
    assert not re.search(r'--wcs-cache\)', text)     # not a ZLP option: jobqueue's argv pinning is untouched
