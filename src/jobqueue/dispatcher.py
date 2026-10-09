"""The dispatcher daemon: one copy, held by flock, polling every poll_seconds.

Each poll it
  1. checks running jobs: collects exit codes, kills a job that has run past timeout_factor times its
     estimate (and one the operator cancelled with --kill), requeues once a job whose process vanished;
  2. reads load, memory and disk %util from /proc, and the locks held by pipeline or download
     processes the queue did not start;
  3. starts whatever scheduler.plan admits, each job in its own session with its own log.

Running jobs are never paused or killed to make room: the plate solver's 60 s SIGALRM and the
pipeline's open files make that unsafe. In shadow mode nothing is launched: an admitted job is
recorded as started and finishes after its estimated run time, so the shadow queue holds cores and
locks the way the real one would.
"""
import datetime as dt
import errno
import fcntl
import hashlib
import json
import os
import signal
import socket
import subprocess
import sys
import time
import traceback

from .config import camera, ensure_dirs, queue_paths
from .resources import ProcResources, external_locks, group_alive, pid_alive, proc_start_time, reap
from .scheduler import plan
from .signatures import extract, is_transient, read_tail
from .util import iso, last_night, night_str, parse_iso, read_json, utcnow

FITS_SUFFIXES = ('.fits', '.fts', '.fit', '.fits.fz', '.fits.gz')


class AlreadyRunning(Exception):
    pass


def acquire_single_instance(lock_path):
    """Hold an exclusive flock for the life of the process; raise AlreadyRunning if another holds it."""
    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o664)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as e:
        os.close(fd)
        if e.errno in (errno.EAGAIN, errno.EACCES, errno.EWOULDBLOCK):
            raise AlreadyRunning(lock_path)
        raise
    os.ftruncate(fd, 0)
    os.write(fd, '{}\n'.format(os.getpid()).encode())
    return fd


def lock_is_held(lock_path):
    try:
        fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o664)
    except OSError:
        return False
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return True
    else:
        fcntl.flock(fd, fcntl.LOCK_UN)
        return False
    finally:
        os.close(fd)


def deployed_sha(src_dir):
    """The deployed commit from src/VERSION, else a fingerprint of the tree's .py and .sh files.

    The production tree is a working copy, not a git checkout, so `git rev-parse` there is meaningless.
    """
    version = os.path.join(src_dir, 'VERSION')
    try:
        with open(version) as f:
            first = f.readline().split()
        if first:
            return first[0]
    except OSError:
        pass
    h = hashlib.sha1()
    for root, dirs, files in os.walk(src_dir):
        dirs[:] = sorted(d for d in dirs if not d.startswith('.') and d not in ('__pycache__', 'tests', 'catcache'))
        for name in sorted(files):
            if not name.endswith(('.py', '.sh')) or name.startswith('._'):
                continue
            path = os.path.join(root, name)
            try:
                with open(path, 'rb') as f:
                    digest = hashlib.sha1(f.read()).hexdigest()
            except OSError:
                continue
            h.update(os.path.relpath(path, src_dir).encode())
            h.update(digest.encode())
    return 'tree-' + h.hexdigest()[:12]


def count_frames(path):
    try:
        return sum(1 for n in os.listdir(path) if n.endswith(FITS_SUFFIXES))
    except OSError:
        return 0


def pipeline_estimate(cfg, telescope, night, frames):
    rate = cfg['minutes_per_1000_frames'][camera(cfg, telescope, night)]
    return max(float(cfg['min_estimate_minutes']), frames / 1000.0 * rate)


class Dispatcher:
    def __init__(self, cfg, store, resources=None, launcher=None, now_fn=None, log=None, proc='/proc'):
        self.cfg = cfg
        self.store = store
        self.paths = ensure_dirs(cfg)
        self.shadow = cfg['mode'] == 'shadow'
        self.resources = resources or ProcResources(cfg['disk_devices'])
        self.launcher = launcher or self._launch
        self.now = now_fn or utcnow
        self.log = log or (lambda msg: print('{} {}'.format(iso(utcnow()), msg), flush=True))
        self.proc = proc
        self.stop = False
        self._children = set()
        self._waiting = {}
        self._sha = (None, 0.0)
        self._drain_logged = False

    # ------------------------------------------------------------------ helpers

    def code_sha(self):
        sha, at = self._sha
        if sha is None or time.time() - at > 600:
            sha = deployed_sha(self.cfg['src_dir'])
            self._sha = (sha, time.time())
        return sha

    def _job_paths(self, job, attempt, now):
        d = os.path.join(self.paths['job_logs'], now.strftime('%Y-%m'))
        os.makedirs(d, exist_ok=True)
        stem = '{:06d}-{}-{}-{}-a{}'.format(job['id'], job['kind'], job['telescope'] or 'x', job['night'] or 'x',
                                            attempt)
        log = os.path.join(d, stem + '.log')
        return log, log + '.status.json'

    def _estimate(self, job):
        if job['kind'] != 'pipeline' or not job['telescope'] or not job['night']:
            return float(job['est_minutes'])
        frames = count_frames(os.path.join(self.cfg['basedir'], 'Observations', job['telescope'], 'images',
                                           job['night']))
        if frames == 0:
            return float(job['est_minutes'])
        return pipeline_estimate(self.cfg, job['telescope'], job['night'], frames)

    def _launch(self, job, log_path, status_path, sha):
        env = os.environ.copy()
        env.update({k: str(v) for k, v in job['env'].items()})
        env['JOBQUEUE_JOB_ID'] = str(job['id'])
        if job['kind'] == 'pipeline':
            env['N_CORES'] = str(job['cores'])
        cwd = job['cwd'] if job['cwd'] and os.path.isdir(job['cwd']) else None
        cmd = [sys.executable, '-m', 'jobqueue.runner', '--log', log_path, '--status', status_path,
               '--job', str(job['id']), '--sha', sha or '']
        if cwd:
            cmd += ['--cwd', cwd]
        cmd += ['--'] + list(job['argv'])
        with open(log_path, 'ab') as errlog:
            p = subprocess.Popen(cmd, start_new_session=True, stdin=subprocess.DEVNULL, stdout=errlog,
                                 stderr=errlog, close_fds=True, env=env, cwd=cwd)
        return p.pid

    def _kill(self, job, sig, reason):
        if os.path.isdir(self.proc) and not job['proc_start']:
            # without the start time a reused PID could be someone else's process group
            self.log('job {}: not signalling pid {}: its start time is unknown'.format(job['id'], job['pid']))
            return
        try:
            os.killpg(job['pid'], sig)
        except ProcessLookupError:
            pass
        except PermissionError as e:
            self.log('job {}: cannot signal process group {}: {}'.format(job['id'], job['pid'], e))
            return
        self.store.mark_kill_sent(job['id'], reason, signal.Signals(sig).name)
        self.log('job {}: {} sent ({})'.format(job['id'], signal.Signals(sig).name, reason))

    def timeout_minutes(self, job):
        return max(float(self.cfg['min_timeout_minutes']), float(self.cfg['timeout_factor']) * float(job['est_minutes']))

    # ------------------------------------------------------------------ running jobs

    def check_running(self, job, now):
        started = parse_iso(job['started_at']) or now
        if job['shadow']:
            if now >= started + dt.timedelta(minutes=float(job['est_minutes'])):
                self.store.finish(job['id'], 0, shadow=True, message='simulated run of {:.0f} min ended'.format(
                    float(job['est_minutes'])))
            return
        pid = job['pid']
        if job['status_path'] and (not pid or not job['proc_start']):
            info = read_json(job['status_path'] + '.started')
            if info and info.get('pid'):
                pid = info['pid']
                self.store.set_process(job['id'], pid, info.get('proc_start') or job['proc_start'])
                job = self.store.get_job(job['id'])
        if pid:
            reap(pid)
        # The runner writes its status file as its last act, so a status file means the job is over, whatever
        # now holds its PID.
        status = read_json(job['status_path']) if job['status_path'] else None
        leader = bool(pid) and pid_alive(pid, job['proc_start'], self.proc)
        if status is None and pid and not leader and group_alive(pid):
            # The runner is gone but the command it started is not: still running, exit status unknown. It is
            # not signalled any more (its group can no longer be verified); q status shows it overrunning.
            if self._waiting.get(('orphan', job['id'])) is None:
                self._waiting[('orphan', job['id'])] = True
                self.log('job {}: runner {} gone but its process group is alive; waiting for it'.format(job['id'], pid))
                self.store.event('runner-gone', job_id=job['id'], source='dispatcher',
                                 message='runner {} gone, process group still alive'.format(pid))
            return
        if status is None and leader:
            runtime = (now - started).total_seconds() / 60.0
            if job['cancel_requested'] and not job['kill_sent_at']:
                self._kill(job, signal.SIGTERM, 'cancelled')
            elif runtime > self.timeout_minutes(job) and not job['kill_sent_at']:
                self._kill(job, signal.SIGTERM, 'timeout')
            elif job['kill_sent_at'] and not (job['kill_reason'] or '').endswith(':KILL'):
                waited = (now - parse_iso(job['kill_sent_at'])).total_seconds()
                if waited > float(self.cfg['kill_grace_seconds']):
                    self._kill(job, signal.SIGKILL, '{}:KILL'.format(job['kill_reason'] or 'timeout'))
            return

        if status is not None:
            self._finalise(job, status.get('exit_code'))
        elif job['kill_reason']:
            self._finalise(job, -9)
        elif not pid and (now - started).total_seconds() < 30:
            return  # still launching
        else:
            out = self.store.vanish(job['id'], int(self.cfg['vanished_requeues']))
            self.log('job {}: process gone without an exit status -> {}'.format(job['id'], out['state']))

    def _finalise(self, job, exit_code):
        tail = read_tail(job['log_path']) if job['log_path'] else ''
        reason = (job['kill_reason'] or '').split(':')[0] or None
        if job['cancel_requested']:
            reason = 'cancelled'
        if reason == 'timeout':
            sig = 'timeout: killed after {:.0f} min ({:g}x estimate of {:.0f} min)'.format(
                self.timeout_minutes(job), float(self.cfg['timeout_factor']), float(job['est_minutes']))
            out = self.store.finish(job['id'], exit_code, failure_kind='timeout', signature=sig)
        elif reason == 'cancelled':
            out = self.store.finish(job['id'], exit_code, failure_kind='cancelled', signature='cancelled by operator')
        elif exit_code == 0:
            out = self.store.finish(job['id'], 0)
        else:
            sig = extract(tail, exit_code)
            transient = is_transient(sig, tail, exit_code, self.cfg)
            out = self.store.finish(job['id'], exit_code, signature=sig, transient=transient,
                                    retry_delay_minutes=float(self.cfg['retry_delay_minutes']))
        self.log('job {} ({} {} {} {}) exit {} -> {}{}'.format(
            job['id'], job['class'], job['kind'], job['telescope'], job['night'], exit_code, out['state'],
            ': ' + out['failure_signature'] if out.get('failure_signature') else ''))

    # ------------------------------------------------------------------ starting jobs

    def start(self, job, now):
        attempt = job['attempts'] + 1
        sha = self.code_sha()
        est = self._estimate(job)
        if self.shadow:
            if self.store.mark_running(job['id'], code_sha=sha, est_minutes=est, shadow=True,
                                       message='would start now ({} cores, est {:.0f} min)'.format(job['cores'], est)):
                self.log('[shadow] would start job {} {} {} {} {}'.format(job['id'], job['class'], job['kind'],
                                                                      job['telescope'], job['night']))
            return
        log_path, status_path = self._job_paths(job, attempt, now)
        if not self.store.mark_running(job['id'], log_path=log_path, status_path=status_path, code_sha=sha,
                                       est_minutes=est):
            return
        try:
            pid = self.launcher(job, log_path, status_path, sha)
        except Exception as e:
            self.store.finish(job['id'], 127, signature='launch failed: {}'.format(e), transient=True,
                              retry_delay_minutes=float(self.cfg['retry_delay_minutes']))
            self.log('job {}: launch failed: {}'.format(job['id'], e))
            return
        self._children.add(pid)
        self.store.set_process(job['id'], pid, proc_start_time(pid, self.proc))
        self.log('started job {} {} {} {} {} pid {} est {:.0f} min log {}'.format(
            job['id'], job['class'], job['kind'], job['telescope'], job['night'], pid, est, log_path))

    def _note_waiting(self, job, waiting):
        key = tuple(sorted(waiting))
        if self._waiting.get(job['id']) == key:
            return
        self._waiting[job['id']] = key
        self.store.event('waiting', job_id=job['id'], source='dispatcher', telescope=job['telescope'],
                         night=job['night'], message=', '.join(key), data={'waiting': list(key)})

    def p0_settled(self, now):
        night = night_str(last_night(now))
        settled = set()
        for tel in self.cfg['telescopes']:
            if self.store.jobs(classes=['P0'], kind='pipeline', telescope=tel, night=night, limit=1):
                settled.add(tel)
                continue
            w = self.store.watch_get(tel, night)
            if w and w['state'] == 'no_data' and w['noted']:
                settled.add(tel)
        return settled

    # ------------------------------------------------------------------ the poll

    def heartbeat(self, now):
        info = {'pid': os.getpid(), 'proc_start': proc_start_time(os.getpid(), self.proc),
                'host': socket.gethostname(), 'mode': self.cfg['mode'], 'heartbeat': iso(now),
                'started_at': getattr(self, '_started_at', iso(now))}
        self.store.kv_set('dispatcher', json.dumps(info, sort_keys=True))

    def poll(self):
        now = self.now()
        self.heartbeat(now)
        for pid in list(self._children):
            if reap(pid) != 'running':
                self._children.discard(pid)
        for job in self.store.jobs(states=['running']):
            try:
                self.check_running(job, now)
            except Exception:
                self.log('job {}: check failed:\n{}'.format(job['id'], traceback.format_exc()))

        if self.store.draining():
            if not self._drain_logged:
                self.log('draining: no new jobs start')
                self._drain_logged = True
            return []
        self._drain_logged = False

        queued = self.store.jobs(states=['queued'])
        dep_states = self.store.states_of([j['depends_on'] for j in queued])
        for job in queued:
            dep = job['depends_on']
            if dep is not None and job['dep_requires_success'] and dep_states.get(dep) in ('failed', 'cancelled'):
                self.store.cancel(job['id'], source='dispatcher')
        queued = self.store.jobs(states=['queued'])
        running = self.store.jobs(states=['running'])
        external = external_locks([j['pid'] for j in running if j['pid']], self.proc)
        snapshot = self.resources.snapshot()
        decisions = plan(queued, running, snapshot, self.cfg, now, paused=self.store.paused(),
                         external=external, dep_states=dep_states, p0_settled=self.p0_settled(now))
        for d in decisions:
            if d.start:
                self.start(d.job, now)
            else:
                self._note_waiting(d.job, d.waiting)
        return decisions

    def run(self, once=False):
        lock_fd = acquire_single_instance(self.paths['lock'])
        self._started_at = iso(self.now())
        self.log('dispatcher started: pid {} mode {} db {}'.format(os.getpid(), self.cfg['mode'], self.store.path))
        self.store.event('dispatcher-start', source='dispatcher', message='pid {} mode {}'.format(
            os.getpid(), self.cfg['mode']))

        def on_signal(signum, frame):
            self.stop = True

        signal.signal(signal.SIGTERM, on_signal)
        signal.signal(signal.SIGINT, on_signal)
        try:
            while not self.stop:
                try:
                    self.poll()
                except Exception:
                    self.log('poll failed:\n' + traceback.format_exc())
                if once:
                    break
                deadline = time.monotonic() + float(self.cfg['poll_seconds'])
                while not self.stop and time.monotonic() < deadline:
                    time.sleep(1)
        finally:
            self.store.event('dispatcher-stop', source='dispatcher', message='pid {}'.format(os.getpid()))
            self.log('dispatcher stopping; running jobs carry on in their own sessions')
            os.close(lock_fd)


# ---------------------------------------------------------------------- watchdog

def watchdog(cfg, store, now=None, proc='/proc', kill=os.kill, log=print):
    """Called from cron every 5 minutes. Returns 0 if a healthy dispatcher holds the lock, 3 if one must be
    started (the cron script then runs `docker exec -d orchard-server python -m jobqueue dispatcher`),
    4 if a hung dispatcher was sent SIGTERM (the next run restarts it)."""
    paths = ensure_dirs(cfg)
    now = now or utcnow()
    if not lock_is_held(paths['lock']):
        log('{} no dispatcher holds {}; start one'.format(iso(now), paths['lock']))
        return 3
    info = json.loads(store.kv_get('dispatcher') or '{}')
    beat = parse_iso(info.get('heartbeat'))
    stale = float(cfg['heartbeat_stale_seconds'])
    if beat is None or (now - beat).total_seconds() <= stale:
        return 0
    pid = info.get('pid')
    if pid and pid_alive(pid, info.get('proc_start'), proc):
        log('{} dispatcher pid {} has not polled since {}; sending SIGTERM'.format(iso(now), pid, info['heartbeat']))
        store.event('watchdog-kill', source='watchdog', message='pid {} heartbeat {}'.format(pid, info['heartbeat']))
        try:
            kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        return 4
    return 0


def spawn_detached(cfg, extra_env=None):
    """Start a dispatcher in its own session (for use outside docker exec -d)."""
    paths = queue_paths(cfg)
    env = os.environ.copy()
    env.update(extra_env or {})
    with open(os.path.join(paths['logs'], 'dispatcher.log'), 'ab') as out:
        p = subprocess.Popen([sys.executable, '-m', 'jobqueue', 'dispatcher'], start_new_session=True,
                             stdin=subprocess.DEVNULL, stdout=out, stderr=out, close_fds=True, env=env)
    return p.pid
