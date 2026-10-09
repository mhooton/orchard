"""Wrapper that runs one job and records how it ended.

The dispatcher starts this in a new session, so the job and everything it spawns share one process
group (killed together only on a timeout or an explicit cancel). The wrapper writes
<status>.started as soon as it runs and <status> when the command exits, so a dispatcher that
restarts in the meantime can still learn the exit code of a job it did not see finish.

    python -m jobqueue.runner --log LOG --status STATUS [--cwd DIR] [--job ID] -- command args...
"""
import argparse
import os
import signal
import socket
import subprocess
import sys

from .resources import proc_start_time
from .util import iso, utcnow, write_json_atomic


def main(argv=None):
    p = argparse.ArgumentParser(prog='jobqueue.runner')
    p.add_argument('--log', required=True)
    p.add_argument('--status', required=True)
    p.add_argument('--cwd')
    p.add_argument('--job')
    p.add_argument('--sha', default='')
    p.add_argument('cmd', nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd
    if not cmd:
        p.error('no command')

    started = iso(utcnow())
    write_json_atomic(a.status + '.started', {'pid': os.getpid(), 'proc_start': proc_start_time(os.getpid()),
                                              'started_at': started, 'job': a.job})
    terminated = []

    def on_term(signum, frame):
        # The signal went to the whole process group, so the command has it too; wait for it to exit.
        terminated.append(signum)

    signal.signal(signal.SIGTERM, on_term)
    signal.signal(signal.SIGHUP, signal.SIG_IGN)

    os.makedirs(os.path.dirname(os.path.abspath(a.log)), exist_ok=True)
    with open(a.log, 'ab', buffering=0) as log:
        log.write('=== jobqueue job {} | started {} | host {} | pid {} | code {}\n=== cwd {}\n=== cmd {}\n'.format(
            a.job, started, socket.gethostname(), os.getpid(), a.sha or '-', a.cwd or os.getcwd(),
            ' '.join(cmd)).encode())
        try:
            child = subprocess.Popen(cmd, cwd=a.cwd or None, stdin=subprocess.DEVNULL, stdout=log,
                                     stderr=subprocess.STDOUT)
        except OSError as e:
            log.write('=== could not start: {}\n'.format(e).encode())
            rc = 127
        else:
            while True:
                try:
                    rc = child.wait()
                    break
                except InterruptedError:
                    continue
        finished = iso(utcnow())
        log.write('=== jobqueue job {} | finished {} | exit {}\n'.format(a.job, finished, rc).encode())

    write_json_atomic(a.status, {'exit_code': rc, 'started_at': started, 'finished_at': finished,
                                 'terminated_by_signal': terminated[0] if terminated else None})
    return 0


if __name__ == '__main__':
    sys.exit(main())
