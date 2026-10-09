"""Host readings for admission, from /proc (the container sees the host's load, memory and disks).

- load: 5-minute load average from /proc/loadavg
- memory: MemAvailable from /proc/meminfo
- disk: %util per device from the io_ticks column of /proc/diskstats, averaged since the previous
  reading (the dispatcher's poll interval), as iostat computes it
- processes: liveness with PID-reuse protection (start time from /proc/<pid>/stat), and a scan for
  pipeline or download processes that the queue did not start (the legacy cron, manual runs)
"""
import os
import time
from collections import namedtuple

Snapshot = namedtuple('Snapshot', 'load5 mem_available_gb disk_util cpu_count')


class ProcResources:
    def __init__(self, devices, proc='/proc', clock=time.monotonic, sleep=time.sleep, first_sample_seconds=3.0):
        self.devices = list(devices)
        self.proc = proc
        self.clock = clock
        self.sleep = sleep
        self.first_sample_seconds = first_sample_seconds
        self._last = None  # (t, {dev: io_ticks})

    def _read(self, name):
        with open(os.path.join(self.proc, name)) as f:
            return f.read()

    def load5(self):
        try:
            return float(self._read('loadavg').split()[1])
        except (OSError, ValueError, IndexError):
            return None

    def mem_available_gb(self):
        try:
            for line in self._read('meminfo').splitlines():
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) / 1024.0 / 1024.0
        except (OSError, ValueError):
            pass
        return None

    def io_ticks(self):
        ticks = {}
        try:
            for line in self._read('diskstats').splitlines():
                parts = line.split()
                if len(parts) >= 13 and parts[2] in self.devices:
                    ticks[parts[2]] = int(parts[12])
        except (OSError, ValueError):
            pass
        return ticks

    def disk_util(self):
        """{device: %util} since the previous call; the first call samples for a few seconds."""
        now, ticks = self.clock(), self.io_ticks()
        if self._last is None:
            if not ticks:
                return {d: None for d in self.devices}
            self._last = (now, ticks)
            self.sleep(self.first_sample_seconds)
            now, ticks = self.clock(), self.io_ticks()
        t0, prev = self._last
        self._last = (now, ticks)
        elapsed_ms = (now - t0) * 1000.0
        out = {}
        for d in self.devices:
            if d in ticks and d in prev and elapsed_ms > 0:
                out[d] = max(0.0, min(100.0, 100.0 * (ticks[d] - prev[d]) / elapsed_ms))
            else:
                out[d] = None
        return out

    def snapshot(self):
        return Snapshot(self.load5(), self.mem_available_gb(), self.disk_util(), os.cpu_count())


# ------------------------------------------------------------------------- processes

def proc_start_time(pid, proc='/proc'):
    """Start time of a process in clock ticks since boot, as a string, or None (no /proc or no process)."""
    try:
        with open(os.path.join(proc, str(pid), 'stat')) as f:
            stat = f.read()
    except OSError:
        return None
    fields = stat[stat.rfind(')') + 2:].split()
    return fields[19] if len(fields) > 19 else None


def _proc_state(pid, proc='/proc'):
    try:
        with open(os.path.join(proc, str(pid), 'stat')) as f:
            stat = f.read()
    except OSError:
        return None
    return stat[stat.rfind(')') + 2:].split()[0]


def reap(pid):
    """Collect a finished child so it does not linger as a zombie.

    Returns 'reaped' (it had exited), 'running' (our child, still running) or 'gone' (not our child).
    """
    try:
        done, _ = os.waitpid(pid, os.WNOHANG)
    except ChildProcessError:
        return 'gone'
    except OSError:
        return 'gone'
    return 'reaped' if done == pid else 'running'


def pid_alive(pid, proc_start=None, proc='/proc'):
    """True while the process exists, is not a zombie and (where /proc exists) has the recorded start time."""
    if not pid:
        return False
    if os.path.isdir(proc):
        state = _proc_state(pid, proc)
        if state is None or state == 'Z':
            return False
        if proc_start and proc_start_time(pid, proc) != str(proc_start):
            return False
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def group_alive(pgid):
    """True while any process of the group remains (a job's command can outlive a runner that was killed)."""
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def session_id(pid, proc='/proc'):
    try:
        with open(os.path.join(proc, str(pid), 'stat')) as f:
            stat = f.read()
        return int(stat[stat.rfind(')') + 2:].split()[3])
    except (OSError, ValueError, IndexError):
        return None


ZLP_FLAGS_WITH_VALUE = {'--cores'}


def parse_zlp_argv(argv):
    """Locks held by a ZLP_pipeline.sh invocation: positional RUNNAME BASEDIR DATES CTHRESH STHRESH TEL [TARGETS]."""
    try:
        i = next(k for k, a in enumerate(argv) if a.endswith('ZLP_pipeline.sh'))
    except StopIteration:
        return None
    pos, args = [], argv[i + 1:]
    k = 0
    while k < len(args):
        a = args[k]
        if a in ZLP_FLAGS_WITH_VALUE:
            k += 2
            continue
        if a.startswith('--'):
            k += 1
            continue
        pos.append(a)
        k += 1
    if len(pos) < 6:
        return None
    tel, dates = pos[5], pos[2].split()
    locks = {'tel:' + tel} | {'night:{}:{}'.format(tel, d) for d in dates}
    targets = ' '.join(pos[6:]).split()
    locks |= {target_lock(t) for t in targets}
    return locks


def parse_sso_download_argv(argv):
    """Locks held by an SSO_download.py run: the telescope-nights it writes into."""
    if not any(a.endswith('SSO_download.py') for a in argv):
        return None
    opts = {}
    for k, a in enumerate(argv):
        if a in ('--telescope', '--sdate', '--edate') and k + 1 < len(argv):
            opts[a] = argv[k + 1]
    tel, s, e = opts.get('--telescope'), opts.get('--sdate'), opts.get('--edate')
    if not tel or not s:
        return None
    import datetime as dt
    try:
        d0 = dt.datetime.strptime(s, '%Y%m%d')
        d1 = dt.datetime.strptime(e, '%Y%m%d') if e else d0 + dt.timedelta(days=1)
    except ValueError:
        return None
    tels = ['Io', 'Europa', 'Ganymede', 'Callisto'] if tel == 'all' else [tel]
    nights = [(d0 + dt.timedelta(days=n)).strftime('%Y%m%d') for n in range(max(1, (d1 - d0).days))]
    return {'night:{}:{}'.format(t, n) for t in tels for n in nights}


def target_lock(name):
    return 'target:' + name.replace(' ', '--').upper()


def external_locks(own_sessions, proc='/proc'):
    """Locks held by pipeline and download processes the queue did not start.

    Returns {lock: description}. Our own jobs run in their own sessions, so any matching process
    whose session is not one of own_sessions belongs to someone else (cron, a person, a replay).
    """
    held = {}
    if not os.path.isdir(proc):
        return held
    own = set(s for s in own_sessions if s)
    for name in os.listdir(proc):
        if not name.isdigit():
            continue
        try:
            with open(os.path.join(proc, name, 'cmdline'), 'rb') as f:
                argv = [a.decode('utf-8', 'replace') for a in f.read().split(b'\0') if a]
        except OSError:
            continue
        if not argv:
            continue
        locks = parse_zlp_argv(argv) or parse_sso_download_argv(argv)
        if not locks:
            continue
        sid = session_id(int(name), proc)
        if sid in own:
            continue
        desc = 'pid {} {}'.format(name, ' '.join(argv)[:160])
        for lock in locks:
            held.setdefault(lock, desc)
    return held
