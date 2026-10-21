"""
T7's process pool: one plate solve at a time in each worker, and a worker that dies or stops answering costs only
the frame it had.

multiprocessing.Pool loses a dead worker's chunk of tasks, and Pool.map then waits for them for ever. On appct
(2026-10-20/21) T7 workers killed by OpenBLAS's SIGINT when orchard-server ran out of threads, and workers spinning in
sched_yield out of reach of the 60 s SIGALRM, left T7 waiting until the queue's job timeout, 5-15 h later.

Here the parent hands each worker one frame at a time, so it knows which frame a worker had when it died (an exit
of any kind: a signal, an uncaught KeyboardInterrupt, an OOM kill) or gave no result within `hang_after` seconds (it is
then killed). That frame is tried again in a new worker, up to `attempts` times in all, and then fails with the reason.
The other workers' frames carry on. A run in which frames keep being lost is not a few bad frames: when more than
`loss_limit` frames have been lost on every attempt, run() kills the workers and raises WorkersLost, which fails T7.
Nothing is printed in a run that loses no worker.
"""
import collections
import multiprocessing
import multiprocessing.connection
import os
import signal
import sys
import time

ATTEMPTS = 2            # a frame whose worker died or hung is tried once more, in a new worker
LOSS_LIMIT = 10         # frames lost on every attempt before the run stops with WorkersLost
POLL = 5.0              # seconds between checks for workers past hang_after
START_TRIES = 12        # failed starts in a row, POLL apart, with no worker left, before WorkersLost


class WorkersLost(RuntimeError):
    """Frames keep being lost with their workers, or no worker process can be started."""


def _serve(conn, func, parent):
    """A worker: run func on each job the parent sends, until told to stop or the parent has gone."""
    while True:
        # A worker forked after this one holds a copy of the parent's end of this pipe, so if the parent is killed
        # recv() may never see EOF: check the parent while waiting.
        while not conn.poll(POLL):
            if os.getppid() != parent:
                return
        try:
            task = conn.recv()
        except EOFError:
            return
        if task is None:
            return
        k, job = task
        try:
            outcome = (True, func(job))
        except Exception as e:      # func's own failure: the parent fails the frame, as Pool.map would have failed T7
            outcome = (False, '{}: {}'.format(type(e).__name__, e))
        sys.stdout.flush()          # its log lines now, not at exit: a worker that dies later takes none with it
        try:
            conn.send((k, outcome))
        except OSError:             # the parent has gone
            return


class _Worker(object):
    def __init__(self, ctx, func):
        self.conn, child = ctx.Pipe()
        self.proc = ctx.Process(target=_serve, args=(child, func, os.getpid()), daemon=True)
        self.proc.start()
        child.close()
        self.task = None            # index of the job it has, or None when idle
        self.started = None

    def give(self, k, job):
        self.conn.send((k, job))
        self.task, self.started = k, time.monotonic()

    def end(self, wait=0.0):
        """Wait up to `wait` seconds for the process to exit, then kill it. Returns its exit code."""
        self.proc.join(wait)
        if self.proc.exitcode is None:
            self.proc.kill()
            self.proc.join()
        self.conn.close()
        return self.proc.exitcode


def _exit_reason(code):
    if code is not None and code < 0:
        try:
            return 'killed by {}'.format(signal.Signals(-code).name)
        except ValueError:
            return 'killed by signal {}'.format(-code)
    return 'exited with status {}'.format(code)


def run(func, jobs, nproc, hang_after, lost, label=str, attempts=ATTEMPTS, loss_limit=LOSS_LIMIT, poll=POLL):
    """
    [func(job) for job in jobs], in up to nproc worker processes, in the order of jobs.

    A job whose worker died or gave no result within hang_after seconds is tried again in a new worker, up to
    `attempts` times in all; after that, and for a job whose func raised, its result is lost(job, why). label(job)
    names the job in the log. Raises WorkersLost when more than loss_limit jobs have been lost on every attempt, or
    when no worker can be started.
    """
    ctx = multiprocessing.get_context('fork')    # as multiprocessing.Pool on Linux: the workers inherit T7's state
    results = [None] * len(jobs)
    todo = collections.deque(range(len(jobs)))
    tries = [0] * len(jobs)
    finished = 0
    workers = []
    died = hung = given_up = start_failures = 0

    def settle(k, value):
        nonlocal finished
        results[k] = value
        finished += 1

    def lose(w, why):
        """w has died or been killed: retire it, and try its job again or give the job up."""
        nonlocal given_up
        workers.remove(w)
        k = w.task
        if k is None:
            return
        name = label(jobs[k])
        if tries[k] < attempts:
            print('T7 pool: worker {} {} on {}; trying it again in a new worker (attempt {} of {})'.format(
                w.proc.pid, why, name, tries[k] + 1, attempts), flush=True)
            todo.appendleft(k)
            return
        print('T7 pool: worker {} {} on {}; frame failed (attempt {} of {})'.format(
            w.proc.pid, why, name, tries[k], attempts), flush=True)
        given_up += 1
        settle(k, lost(jobs[k], 'Lost with its worker process: {} (attempt {} of {})'.format(why, tries[k],
                                                                                            attempts)))
        if given_up > loss_limit:
            raise WorkersLost('T7: {} frames lost with their worker processes, the last one {}: {} died, {} hung. '
                              'Stopping rather than failing every frame.'.format(given_up, why, died, hung))

    try:
        while finished < len(jobs):
            # Start workers while there are frames waiting and no idle worker for them.
            while todo and len(workers) < nproc and sum(w.task is None for w in workers) < len(todo):
                try:
                    workers.append(_Worker(ctx, func))
                    start_failures = 0
                except OSError as e:            # fork's EAGAIN when the container is out of threads or processes
                    start_failures += 1
                    print('T7 pool: could not start a worker ({}); {} running'.format(e, len(workers)), flush=True)
                    if not workers and start_failures >= START_TRIES:
                        raise WorkersLost('T7: could not start a worker process {} times in a row: {}'.format(
                            start_failures, e))
                    break
            if not workers:
                time.sleep(poll)
                continue

            for w in list(workers):
                if w.task is None and todo:
                    k = todo.popleft()
                    try:
                        w.give(k, jobs[k])
                    except OSError:             # it died idle: the job was never handed over
                        todo.appendleft(k)
                        w.end()
                        workers.remove(w)
                        continue
                    tries[k] += 1

            now = time.monotonic()
            due = [w.started + hang_after - now for w in workers if w.task is not None]
            ready = set(multiprocessing.connection.wait(
                [w.conn for w in workers] + [w.proc.sentinel for w in workers], timeout=max(0.0, min([poll] + due))))

            for w in list(workers):
                if w.conn not in ready and w.proc.sentinel not in ready:
                    continue
                try:
                    message = w.conn.recv() if w.conn.poll() else None
                except (EOFError, OSError):
                    message = None              # its end of the pipe has closed: it has exited
                if message is not None:         # a result (one sent just before exiting is still read)
                    k, (ok, value) = message
                    if not ok:
                        print('T7 pool: {} raised {}'.format(label(jobs[k]), value), flush=True)
                    settle(k, value if ok else lost(jobs[k], value))
                    w.task = None
                    continue
                code = w.end(wait=10.0)
                if w.task is not None:
                    died += 1
                lose(w, _exit_reason(code))

            now = time.monotonic()
            for w in list(workers):
                if w.task is not None and now - w.started > hang_after and not w.conn.poll():
                    hung += 1
                    w.end()
                    lose(w, 'gave no result in {:.0f} s and was killed'.format(hang_after))
    except BaseException:
        for w in workers:
            w.end()
        raise

    for w in workers:
        try:
            w.conn.send(None)
        except OSError:
            pass
    for w in workers:
        w.end(wait=30.0)
    if died or hung:
        print('T7 pool: {} worker processes lost ({} died, {} hung); {} frames failed with them'.format(
            died + hung, died, hung, given_up), flush=True)
    return results
