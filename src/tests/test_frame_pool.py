"""
T7's process pool (astrom/frame_pool.py): a worker that dies or hangs costs only the frame it had.

multiprocessing.Pool, which T7 used before, loses a dead worker's chunk and Pool.map then waits for ever. On appct
(2026-10-20/21) workers killed by OpenBLAS's SIGINT when orchard-server ran out of threads left T7 waiting until the
queue's job timeout, 5-15 h later (job 1217: 55 dead workers). Here the workers are killed in every way that happened
or could: SIGINT (an uncaught KeyboardInterrupt), SIGKILL (the OOM killer), from inside and from another process, and
a worker stuck where no Python signal handler runs.

The jobs here are toys, so this runs anywhere with fork (Linux, macOS); test_t7_dead_worker.py runs T7 itself.
Run:  pytest src/tests/test_frame_pool.py -q
"""
import multiprocessing
import os
import signal
import subprocess
import sys
import time

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from astrom import frame_pool       # noqa: E402

pytestmark = pytest.mark.skipif(not hasattr(os, 'fork'), reason='T7 forks its workers')


def lost(job, why):
    return ('lost', job, why)


def run(func, jobs, nproc=3, hang_after=30, **kw):
    started = time.monotonic()
    results = frame_pool.run(func, jobs, nproc, hang_after=hang_after, lost=lost, label=repr, **kw)
    assert not multiprocessing.active_children()        # every worker ended, none left behind
    return results, time.monotonic() - started


def square(x):
    time.sleep(0.02)
    return x * x


def once(marker):
    """True the first time it is called with this marker file, in any process."""
    try:
        fd = os.open(marker, os.O_CREAT | os.O_EXCL)
    except FileExistsError:
        return False
    os.close(fd)
    return True


class Dies:
    """square(x), except that a worker given `victim` dies there: every time, or (`once`) only the first time."""

    def __init__(self, victim, sig, tmp_path, every_time=False):
        self.victim, self.sig, self.every_time = victim, sig, every_time
        self.marker = str(tmp_path / 'died')

    def __call__(self, x):
        if x == self.victim and (self.every_time or once(self.marker)):
            time.sleep(0.2)                              # mid-job
            os.kill(os.getpid(), self.sig)
            time.sleep(10)
        return square(x)


def test_results_are_in_order_and_a_healthy_run_prints_nothing(capfd):
    results, _ = run(square, list(range(40)), nproc=4)
    assert results == [x * x for x in range(40)]
    assert capfd.readouterr().out == ''


def test_fewer_jobs_than_workers():
    assert run(square, [3], nproc=20)[0] == [9]
    assert run(square, [], nproc=20)[0] == []


@pytest.mark.parametrize('sig', [signal.SIGINT, signal.SIGKILL, signal.SIGSEGV])
def test_a_worker_that_dies_once_costs_nothing(tmp_path, capfd, sig):
    # SIGINT is OpenBLAS's when it cannot start its threads: a KeyboardInterrupt the worker does not catch
    results, elapsed = run(Dies(7, sig, tmp_path), list(range(20)))
    assert results == [x * x for x in range(20)]        # the lost frame solved in a new worker, the rest untouched
    out = capfd.readouterr().out
    assert 'on 7; trying it again in a new worker (attempt 2 of 2)' in out
    assert '1 worker processes lost (1 died, 0 hung); 0 frames failed with them' in out
    assert elapsed < 20


def test_a_frame_that_kills_every_worker_fails_alone(tmp_path, capfd):
    results, elapsed = run(Dies(5, signal.SIGKILL, tmp_path, every_time=True), list(range(20)))
    assert results[5] == ('lost', 5, 'Lost with its worker process: killed by SIGKILL (attempt 2 of 2)')
    assert results[:5] + results[6:] == [x * x for x in range(20) if x != 5]
    out = capfd.readouterr().out
    assert 'killed by SIGKILL on 5; frame failed (attempt 2 of 2)' in out
    assert '2 worker processes lost (2 died, 0 hung); 1 frames failed with them' in out
    assert elapsed < 20


def stuck(x):
    """Frame 2 never returns and SIGALRM cannot reach it, like a worker spinning in OpenBLAS's sched_yield."""
    if x == 2:
        signal.pthread_sigmask(signal.SIG_BLOCK, [signal.SIGALRM, signal.SIGINT, signal.SIGTERM])
        time.sleep(3600)
    return square(x)


def test_a_hung_worker_is_killed_and_its_frame_fails(capfd):
    results, elapsed = run(stuck, list(range(12)), hang_after=1.0, poll=0.2)
    assert results[2] == ('lost', 2, 'Lost with its worker process: gave no result in 1 s and was killed '
                                     '(attempt 2 of 2)')
    assert [r for k, r in enumerate(results) if k != 2] == [x * x for x in range(12) if x != 2]
    assert '2 worker processes lost (0 died, 2 hung); 1 frames failed with them' in capfd.readouterr().out
    assert elapsed < 10


class KillsASibling:
    """A worker killed by another process while it works: job 0 records its pid and works for 2 s; job 3, in another
    worker, SIGKILLs it."""

    def __init__(self, tmp_path):
        self.pidfile = str(tmp_path / 'pid0')
        self.marker = str(tmp_path / 'killed')

    def __call__(self, x):
        if x == 0 and not os.path.exists(self.marker):
            with open(self.pidfile + '.tmp', 'w') as f:
                f.write(str(os.getpid()))
            os.rename(self.pidfile + '.tmp', self.pidfile)
            time.sleep(2)
        if x == 3:
            while not os.path.exists(self.pidfile):
                time.sleep(0.01)
            open(self.marker, 'w').close()
            os.kill(int(open(self.pidfile).read()), signal.SIGKILL)
        return square(x)


def test_a_worker_killed_from_outside(tmp_path, capfd):
    results, elapsed = run(KillsASibling(tmp_path), list(range(10)), nproc=2)
    assert results == [x * x for x in range(10)]
    assert 'killed by SIGKILL on 0; trying it again' in capfd.readouterr().out
    assert elapsed < 10


def test_frames_lost_over_and_over_stop_the_run(tmp_path, capfd):
    def dies(x):
        os.kill(os.getpid(), signal.SIGKILL)
    started = time.monotonic()
    with pytest.raises(frame_pool.WorkersLost, match='4 frames lost with their worker processes'):
        frame_pool.run(dies, list(range(100)), 4, hang_after=30, lost=lost, loss_limit=3)
    assert time.monotonic() - started < 20
    for _ in range(50):                                 # the rest are killed, not left running
        if not multiprocessing.active_children():
            break
        time.sleep(0.1)
    assert not multiprocessing.active_children()


def test_a_job_that_raises_fails_without_a_retry(capfd):
    def raises(x):
        if x == 1:
            raise ValueError('bad frame')
        return square(x)
    results, _ = run(raises, [0, 1, 2])
    assert results == [0, ('lost', 1, 'ValueError: bad frame'), 4]
    assert 'T7 pool: 1 raised ValueError: bad frame' in capfd.readouterr().out


def test_the_worker_log_lines_come_out_as_each_frame_finishes(tmp_path):
    # A worker's print() is flushed after each frame, so a worker that dies later takes none of its lines with it.
    script = tmp_path / 'lines.py'
    script.write_text(
        'import os, signal, sys, time\n'
        'sys.path.insert(0, {src!r})\n'
        'from astrom import frame_pool\n'
        'def job(x):\n'
        '    if x == 4:\n'
        '        time.sleep(0.5)\n'
        '        os.kill(os.getpid(), signal.SIGKILL)\n'
        '    print("frame", x, "done")\n'
        '    return x\n'
        'frame_pool.run(job, list(range(5)), 1, hang_after=30, lost=lambda j, w: None)\n'.format(src=SRC))
    out = subprocess.run([sys.executable, str(script)], stdout=subprocess.PIPE, universal_newlines=True,
                         timeout=60).stdout                 # stdout a pipe, so block-buffered, as in a job's log
    assert [line for line in out.splitlines() if line.startswith('frame')] == \
        ['frame {} done'.format(x) for x in range(4)]


def test_workers_exit_when_t7_is_killed(tmp_path):
    # The queue's timeout or a person can kill T7; its workers must not live on.
    script = tmp_path / 'killed.py'
    pids = tmp_path / 'pids'
    script.write_text(
        'import os, sys, time\n'
        'sys.path.insert(0, {src!r})\n'
        'from astrom import frame_pool\n'
        'def job(x):\n'
        '    with open({pids!r}, "a") as f:\n'
        '        f.write("%d\\n" % os.getpid())\n'
        '    time.sleep(1)\n'
        '    return x\n'
        'frame_pool.run(job, list(range(100)), 2, hang_after=30, lost=lambda j, w: None)\n'.format(
            src=SRC, pids=str(pids)))
    t7 = subprocess.Popen([sys.executable, str(script)])
    for _ in range(100):
        if pids.exists() and len(set(pids.read_text().split())) >= 2:
            break
        time.sleep(0.1)
    t7.kill()
    t7.wait()
    workers = {int(p) for p in pids.read_text().split()}
    deadline = time.monotonic() + 3 * frame_pool.POLL + 5
    while time.monotonic() < deadline and any(_alive(p) for p in workers):
        time.sleep(0.2)
    assert not [p for p in workers if _alive(p)]


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # an orphan that has exited may linger as a zombie until init reaps it
    try:
        with open('/proc/{}/stat'.format(pid)) as f:
            return f.read().split(')')[-1].split()[0] != 'Z'
    except OSError:
        return True
