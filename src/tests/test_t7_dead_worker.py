"""
T7 finishes when a plate-solving worker dies or hangs (astrom/astrometry.py on astrom/frame_pool.py).

On appct (2026-10-20/21) T7 workers died of OpenBLAS's SIGINT when orchard-server ran out of threads, and Pool.map
waited for their frames until the queue's job timeout, 5-15 h later (job 1217: 55 dead workers). Now the frame a
worker had is tried again in a new worker, and fails, with the reason in its log line, if it is lost again. The other
frames, their log lines and their WCS cache rows are those of a healthy run.

The solves are real; test_wcs_cache's `solver` fakes only the image cleaning, detection and Gaia query. Needs the
pipeline's own environment (astrometry.py imports twirl): run inside orchard-server.
Run:  pytest src/tests/test_t7_dead_worker.py -q
"""
import os
import re
import signal
import sys
import time

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pytest.importorskip('twirl')
from astrom import astrometry as am, wcs_cache              # noqa: E402
from astropy.io import fits                                 # noqa: E402
from utils import promote_policy                            # noqa: E402
from test_wcs_cache import Night, solver                    # noqa: E402,F401
from test_neighbour_pass import _hard_frame                 # noqa: E402

pytestmark = pytest.mark.skipif(not hasattr(os, 'fork'), reason='T7 forks its workers')


@pytest.fixture
def deadline():
    """Fail, rather than hang the suite, if T7 waits for ever, as Pool.map did for a dead worker's frames."""
    def expired(signum, frame):
        raise AssertionError('T7 still running after 120 s')
    old = signal.signal(signal.SIGALRM, expired)
    signal.alarm(120)
    yield
    signal.alarm(0)
    signal.signal(signal.SIGALRM, old)


def _sabotage(monkeypatch, tmp_path, frames):
    """Run frames[name](marker) in the worker before a frame's solve: it kills or blocks the worker."""
    real = am.pointer_wcs

    def pointer_wcs(filepath, *args, **kwargs):
        name = os.path.basename(filepath)
        if name in frames:
            frames[name](str(tmp_path / (name + '.sabotaged')), kwargs)
        return real(filepath, *args, **kwargs)
    monkeypatch.setattr(am, 'pointer_wcs', pointer_wcs)


def _first_time(marker):
    try:
        os.close(os.open(marker, os.O_CREAT | os.O_EXCL))
    except FileExistsError:
        return False
    return True


def sigint_once(marker, kwargs):
    """OpenBLAS's SIGINT when it cannot start its threads: a KeyboardInterrupt the worker does not catch."""
    if _first_time(marker):
        os.kill(os.getpid(), signal.SIGINT)
        time.sleep(10)


def sigkill(marker, kwargs):
    """The OOM killer, every time this frame is tried."""
    os.kill(os.getpid(), signal.SIGKILL)


def hang(marker, kwargs):
    """Stuck where SIGALRM cannot reach it, like the workers spinning in OpenBLAS's sched_yield."""
    signal.pthread_sigmask(signal.SIG_BLOCK, [signal.SIGALRM])
    time.sleep(3600)


def _rows_are_the_frames_solutions(night, names):
    rows = night.rows()
    assert sorted(rows) == sorted(names)
    for raw, proc in zip(night.raws, night.procs):
        if os.path.basename(raw) in names:
            assert rows[os.path.basename(raw)]['cards'] == wcs_cache.cards_from_header(fits.getheader(proc),
                                                                                        am.WCS_KEYWORDS)


def _summaries(out):
    return {m.group('name'): (m.group('state'), m.group('rest')) for line in out.splitlines()
            for m in promote_policy.SUMMARY.finditer(line)}


def test_a_dead_worker_costs_only_its_frame(tmp_path, solver, monkeypatch, capfd, deadline):
    night = Night(tmp_path, n=6)
    _sabotage(monkeypatch, tmp_path, {'procIo_001.fits': sigint_once, 'procIo_004.fits': sigkill})
    night.t7(nproc=3)
    out = capfd.readouterr().out
    summaries = _summaries(out)
    assert summaries['procIo_001.fits'][0] == 'SUCCESS'                # tried again in a new worker
    assert summaries['procIo_004.fits'] == ('FAILED', ' Error: Lost with its worker process: killed by SIGKILL '
                                                      '(attempt 2 of 2)')
    assert [summaries['procIo_00{}.fits'.format(k)][0] for k in (0, 2, 3, 5)] == ['SUCCESS'] * 4
    assert 'on procIo_001.fits; trying it again in a new worker (attempt 2 of 2)' in out
    assert 'T7 pool: 3 worker processes lost (3 died, 0 hung); 1 frames failed with them' in out
    assert 'Successfully plate solved: 5' in out and 'Failed plate solving: 1' in out
    _rows_are_the_frames_solutions(night, ['Io_00{}.fits'.format(k) for k in (0, 1, 2, 3, 5)])


def test_a_hung_worker_is_killed(tmp_path, solver, monkeypatch, capfd, deadline):
    night = Night(tmp_path, n=4)
    _sabotage(monkeypatch, tmp_path, {'procIo_002.fits': hang})
    monkeypatch.setattr(am, 'HANG_GRACE', 1)
    started = time.monotonic()
    night.t7(nproc=2, timeout=3)                    # killed 4 s after it was handed out, twice
    assert time.monotonic() - started < 60
    out = capfd.readouterr().out
    assert _summaries(out)['procIo_002.fits'] == (
        'FAILED', ' Error: Lost with its worker process: gave no result in 4 s and was killed (attempt 2 of 2)')
    assert 'Successfully plate solved: 3' in out and 'Failed plate solving: 1' in out
    _rows_are_the_frames_solutions(night, ['Io_000.fits', 'Io_001.fits', 'Io_003.fits'])


def test_the_second_pass_survives_a_dead_worker(tmp_path, solver, monkeypatch, capfd, deadline):
    night = Night(tmp_path, n=3)
    _hard_frame(monkeypatch, night.procs[1])        # fails alone, solves from its neighbour (#69)

    def sigint_in_the_second_pass(marker, kwargs):
        if kwargs.get('neighbour'):
            sigint_once(marker, kwargs)
    _sabotage(monkeypatch, tmp_path, {'procIo_001.fits': sigint_in_the_second_pass})
    night.t7(nproc=2)
    out = capfd.readouterr().out
    assert 'Second pass: 1 frames tried again' in out and 'Second pass: 1 of 1 solved' in out
    assert 'on procIo_001.fits; trying it again in a new worker' in out
    assert 'Method: neighbour' in out
    assert 'Successfully plate solved: 3' in out and 'Failed plate solving: 0' in out
    _rows_are_the_frames_solutions(night, ['Io_000.fits', 'Io_001.fits', 'Io_002.fits'])


def test_the_pool_writes_what_a_serial_run_writes(tmp_path, solver, monkeypatch, capfd):
    """A healthy run: the same log lines, headers and cache rows from the pool as from one process."""
    def run(root, nproc):
        night = Night(root, n=5)
        _hard_frame(monkeypatch, night.procs[2])    # so the second pass runs too
        night.t7(nproc=nproc)
        out = capfd.readouterr().out
        lines = sorted(re.sub(r'^\[\d+/\d+\] ', '', line.strip()) for line in out.splitlines()
                       if ' SUCCESS - ' in line or ' FAILED - ' in line)   # the serial run numbers its lines
        headers = [{k: v for k, v in fits.getheader(p).items() if k in am.WCS_KEYWORDS} for p in night.procs]
        rows = {name: row['cards'] for name, row in night.rows().items()}
        return lines, headers, rows, 'T7 pool' in out
    (tmp_path / 'serial').mkdir()
    (tmp_path / 'pool').mkdir()
    serial = run(tmp_path / 'serial', 1)
    pool = run(tmp_path / 'pool', 3)
    assert len(serial[0]) == 6 and serial[0] == pool[0]                # 5 frames and the second pass
    assert serial[1:3] == pool[1:3]
    assert not pool[3]                                                   # nothing from the pool in a healthy run
