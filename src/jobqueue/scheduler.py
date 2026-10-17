"""Which queued jobs may start now. Pure: no I/O, so the rules are easy to test.

Jobs are considered in priority order: class (P0 first), then priority within the class, then age.
A job starts only if all of these hold:

  cores     running declared cores + its cores <= max_cores (112 of 128 threads); outside P0, any
            p0_reserve_cores (0 by default) are taken off that cap during the nightly window
  load      5-minute load average (plus cores started earlier in this pass) < load_max
  mem       MemAvailable > mem_min_gb
  disk      every array's %util < disk_util_max, for disk-heavy jobs outside P0; and at most
            max_disk_heavy_starts_per_poll such starts per pass, since %util lags a new job
  locks     none of its locks is held by a running job or by a process the queue did not start;
            'sem:NAME' locks are counted against semaphores[NAME]; 'tel:TEL' locks count only while
            telescope_lock is set (off: any number of nights of one telescope may run at once)

Strict priority: when a job cannot start, the resources it is waiting for are reserved for its class,
and no job of a lower class may take a reserved resource. Jobs of the same class may still start
alongside it if they do not need what it waits for (a download is not held up by a pipeline job
waiting for cores). Nothing running is ever stopped to make room, and nothing is held back for a
job that has not been queued yet: last night's data goes to the front when it arrives and takes the next
capacity that frees up.
"""
from collections import namedtuple

from . import CLASS_RANK
from .util import in_daily_window, parse_iso

Decision = namedtuple('Decision', 'job start waiting')

GLOBAL = ('cores', 'load', 'mem')


def order_key(job):
    return (CLASS_RANK[job['class']], -float(job.get('priority') or 0), job['created_at'], job['id'])


def needs(job):
    out = set(GLOBAL)
    if job['disk_heavy'] and job['class'] != 'P0':
        out.add('disk')
    out.update(job['locks'])
    return out


def eligible(queued, now, paused, dep_states):
    """Queued jobs that may be considered at all: due, class not paused, dependency settled."""
    out = []
    for job in queued:
        if job['state'] != 'queued' or job['class'] in paused:
            continue
        nb = parse_iso(job.get('not_before'))
        if nb and nb > now:
            continue
        dep = job.get('depends_on')
        if dep is not None and dep_states.get(dep) not in ('done', 'failed', 'cancelled'):
            continue
        out.append(job)
    return out


def plan(queued, running, snapshot, cfg, now, paused=(), draining=False, external=None, dep_states=None):
    """Return a Decision for every eligible queued job, in the order considered."""
    if draining:
        return []
    external = external or {}
    dep_states = dep_states or {}
    telescope_lock = bool(cfg.get('telescope_lock'))
    held = set(lock for lock in external if telescope_lock or not lock.startswith('tel:'))
    sem_count = {}
    used_cores = 0
    for job in running:
        used_cores += int(job['cores'])
        for lock in job['locks']:
            if lock.startswith('tel:') and not telescope_lock:
                continue
            if lock.startswith('sem:'):
                sem_count[lock] = sem_count.get(lock, 0) + 1
            else:
                held.add(lock)
    semaphores = cfg.get('semaphores') or {}

    reserve = 0
    if cfg.get('p0_reserve_cores') and cfg.get('p0_reserve_window_utc'):
        a, b = cfg['p0_reserve_window_utc']
        if in_daily_window(now, a, b):
            reserve = int(cfg['p0_reserve_cores'])

    reserved = {}            # resource -> best (lowest) class rank waiting for it
    started_cores = 0
    disk_heavy_started = 0
    decisions = []
    utils = list((snapshot.disk_util or {}).values()) if snapshot else []

    for job in sorted(eligible(queued, now, set(paused), dep_states), key=order_key):
        rank = CLASS_RANK[job['class']]
        cores = int(job['cores'])
        waiting = []

        cap = int(cfg['max_cores']) - (reserve if job['class'] != 'P0' else 0)
        if used_cores + cores > cap:
            waiting.append('cores')
        if snapshot and snapshot.load5 is not None and snapshot.load5 + started_cores >= float(cfg['load_max']):
            waiting.append('load')
        if snapshot and snapshot.mem_available_gb is not None and snapshot.mem_available_gb < float(cfg['mem_min_gb']):
            waiting.append('mem')
        if job['disk_heavy'] and job['class'] != 'P0':
            if not utils or any(u is None or u >= float(cfg['disk_util_max']) for u in utils):
                waiting.append('disk')
            elif disk_heavy_started >= int(cfg.get('max_disk_heavy_starts_per_poll') or 1):
                waiting.append('disk')
        for lock in job['locks']:
            if lock.startswith('tel:') and not telescope_lock:
                continue
            if lock.startswith('sem:'):
                cap_n = int(semaphores.get(lock[4:], 1))
                if sem_count.get(lock, 0) >= cap_n:
                    waiting.append(lock)
            elif lock in held:
                waiting.append(lock)
        for r in needs(job):
            if r in reserved and reserved[r] < rank and r not in waiting:
                waiting.append(r)

        if waiting:
            for r in waiting:
                reserved[r] = min(reserved.get(r, rank), rank)
            decisions.append(Decision(job, False, waiting))
            continue

        decisions.append(Decision(job, True, []))
        used_cores += cores
        started_cores += cores
        if job['disk_heavy'] and job['class'] != 'P0':
            disk_heavy_started += 1
        for lock in job['locks']:
            if lock.startswith('sem:'):
                sem_count[lock] = sem_count.get(lock, 0) + 1
            else:
                held.add(lock)
    return decisions
