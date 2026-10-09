"""`q`: the queue's command line. Run inside orchard-server on appct:

    docker exec orchard-server python -m jobqueue status
    docker exec orchard-server python -m jobqueue add --class P2 --telescope Europa --night 20250101
    docker exec orchard-server python -m jobqueue pause P2        # resume P2 | resume all
    docker exec orchard-server python -m jobqueue drain           # start nothing new; undrain to undo
    docker exec orchard-server python -m jobqueue cancel 42 [--kill]
    docker exec orchard-server python -m jobqueue requeue --signature "KeyError: 'gaia_dr3_id'"

src/jobqueue/cron/q is a host-side wrapper that does the docker exec for you.
"""
import argparse
import datetime as dt
import json
import os
import sys

from . import CLASSES, jobspec
from .config import ensure_dirs, load_config, queue_paths
from .util import iso, parse_iso, utcnow


def _store(cfg):
    from .store import Store
    ensure_dirs(cfg)
    return Store(queue_paths(cfg)['db'])


def _ago(ts, now):
    t = parse_iso(ts)
    if t is None:
        return '-'
    s = (now - t).total_seconds()
    if s < 120:
        return '{:.0f}s'.format(s)
    if s < 7200:
        return '{:.0f}m'.format(s / 60)
    return '{:.1f}h'.format(s / 3600)


def _waiting_for(store, job_id):
    ev = store.conn.execute("SELECT message FROM events WHERE job_id = ? AND kind = 'waiting' ORDER BY id DESC LIMIT 1",
                            (job_id,)).fetchone()
    return ev['message'] if ev else ''


def cmd_status(cfg, store, a):
    from .dispatcher import lock_is_held
    now = utcnow()
    paths = queue_paths(cfg)
    info = json.loads(store.kv_get('dispatcher') or '{}')
    alive = lock_is_held(paths['lock'])
    print('queue      {} (mode {})'.format(store.path, cfg['mode']))
    if alive:
        print('dispatcher pid {} on {}, last poll {} ago, up since {}'.format(
            info.get('pid'), info.get('host'), _ago(info.get('heartbeat'), now), info.get('started_at')))
    else:
        print('dispatcher NOT RUNNING (last heartbeat {})'.format(info.get('heartbeat', 'never')))
    paused = sorted(store.paused())
    print('drain      {}   paused: {}'.format('ON' if store.draining() else 'off', ', '.join(paused) or 'none'))
    if a.resources:
        from .resources import ProcResources
        s = ProcResources(cfg['disk_devices']).snapshot()
        print('host       load5 {} | MemAvailable {:.0f} GB | {}'.format(
            s.load5, s.mem_available_gb or 0,
            ' '.join('{} {}%'.format(d, '?' if u is None else '{:.0f}'.format(u)) for d, u in s.disk_util.items())))
    since = iso(now - dt.timedelta(hours=24))
    print('\nclass  queued running  done/24h failed/24h')
    for c in CLASSES:
        row = store.conn.execute(
            "SELECT SUM(state = 'queued') q, SUM(state = 'running') r, "
            "SUM(state = 'done' AND finished_at >= ?) d, SUM(state = 'failed' AND finished_at >= ?) f "
            'FROM jobs WHERE class = ?', (since, since, c)).fetchone()
        print('{:6} {:6} {:7} {:9} {:10}'.format(c, row['q'] or 0, row['r'] or 0, row['d'] or 0, row['f'] or 0))

    running = store.jobs(states=['running'])
    print('\nRUNNING ({})'.format(len(running)))
    for j in running:
        print('  {:>6} {} {:8} {:9} {:8} {:>3}c  {:>6} of est {:.0f}m  {}{}'.format(
            j['id'], j['class'], j['kind'], j['telescope'] or '-', j['night'] or '-', j['cores'],
            _ago(j['started_at'], now), j['est_minutes'], 'pid {}'.format(j['pid']) if j['pid'] else '',
            ' [shadow]' if j['shadow'] else ''))
    queued = store.jobs(states=['queued'])
    queued.sort(key=lambda j: (CLASSES.index(j['class']), -j['priority'], j['created_at'], j['id']))
    print('\nQUEUED ({}{})'.format(len(queued), ', first {}'.format(a.limit) if len(queued) > a.limit else ''))
    for j in queued[:a.limit]:
        extra = ''
        if j['not_before'] and parse_iso(j['not_before']) > now:
            extra = 'retry after {}'.format(j['not_before'])
        elif j['depends_on']:
            extra = 'after job {}'.format(j['depends_on'])
        else:
            extra = _waiting_for(store, j['id'])
        print('  {:>6} {} {:8} {:9} {:8} {:>3}c est {:>4.0f}m  {}'.format(
            j['id'], j['class'], j['kind'], j['telescope'] or '-', j['night'] or '-', j['cores'], j['est_minutes'],
            extra))
    failed = store.conn.execute("SELECT * FROM jobs WHERE state = 'failed' AND finished_at >= ? ORDER BY id DESC",
                                (since,)).fetchall()
    if failed:
        print('\nFAILED in the last 24 h ({})'.format(len(failed)))
        for j in failed[:a.limit]:
            print('  {:>6} {} {:8} {:9} {:8} exit {} [{}] {}'.format(
                j['id'], j['class'], j['kind'], j['telescope'] or '-', j['night'] or '-', j['exit_code'],
                j['failure_kind'], j['failure_signature'] or ''))
    return 0


def cmd_show(cfg, store, a):
    j = store.get_job(a.id)
    if not j:
        print('no job {}'.format(a.id))
        return 1
    for k in sorted(j):
        print('{:22} {}'.format(k, j[k]))
    print('\nevents:')
    for e in store.events(job_id=a.id):
        print('  {} {:18} {:10} {}'.format(e['ts'], e['kind'], e['source'], e['message'] or ''))
    return 0


def cmd_add(cfg, store, a):
    if a.kind == 'pipeline':
        if not (a.telescope and a.night):
            print('a pipeline job needs --telescope and --night')
            return 2
        run_targets = a.targets.split() if a.targets else []
        lock_targets = a.lock_targets.split() if a.lock_targets else []
        if not lock_targets and not run_targets:
            from .nights import night_targets
            lock_targets = night_targets(cfg['basedir'], a.telescope, a.night)
        from .dispatcher import count_frames
        from .nights import night_dir
        spec = jobspec.pipeline_job(
            cfg, a.cls, a.telescope, a.night, lock_targets=lock_targets, run_targets=run_targets,
            frames=count_frames(night_dir(cfg['basedir'], a.telescope, a.night)), est_minutes=a.estimate,
            source='cli', depends_on=a.depends_on, priority=a.priority, no_t12=a.no_T12, note=a.note, cores=a.cores,
            extra_flags=a.flag or ())
    else:
        cmd = a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd
        if not cmd:
            print('a command job needs a command after --')
            return 2
        locks = list(a.lock or [])
        if a.telescope and a.night:
            locks.append('night:{}:{}'.format(a.telescope, a.night))
        spec = dict(cls=a.cls, kind=a.kind, telescope=a.telescope, night=a.night, argv=cmd, cwd=a.cwd or cfg['src_dir'],
                    env=jobspec.job_env(cfg), cores=a.cores or 1, disk_heavy=not a.no_disk_heavy,
                    est_minutes=a.estimate or 60, locks=locks, depends_on=a.depends_on, priority=a.priority,
                    source='cli', note=a.note, max_retries=cfg['max_retries'])
    if a.no_disk_heavy:
        spec['disk_heavy'] = False
    jid, created = jobspec.add(store, spec)
    print('{} job {}: {} {} {} {} cores={} est={:.0f}m locks={}'.format(
        'added' if created else 'already queued as', jid, spec['cls'], spec['kind'], spec.get('telescope') or '-',
        spec.get('night') or '-', spec['cores'], spec['est_minutes'], ' '.join(spec['locks'])))
    return 0


def cmd_pause(cfg, store, a):
    store.pause(a.cls, a.reason or '')
    print('paused {}'.format(a.cls))
    return 0


def cmd_resume(cfg, store, a):
    store.resume(a.cls)
    print('resumed {}'.format(a.cls))
    return 0


def cmd_drain(cfg, store, a):
    store.set_drain(not a.off)
    print('drain {}: {}'.format('off' if a.off else 'ON', 'jobs start again' if a.off else
                                'running jobs finish, nothing new starts'))
    return 0


def cmd_cancel(cfg, store, a):
    for jid in a.ids:
        try:
            print('job {}: {}'.format(jid, store.cancel(jid, kill=a.kill)))
        except (KeyError, ValueError) as e:
            print('job {}: {}'.format(jid, e))
    return 0


def cmd_requeue(cfg, store, a):
    if not a.ids and a.signature is None:
        print('give job ids or --signature')
        return 2
    out = store.requeue(job_ids=a.ids or None, signature=a.signature, exact=a.exact, dry_run=a.dry_run)
    for jid, what in out:
        print('job {}: {}'.format(jid, what))
    print('{} job(s)'.format(len(out)))
    return 0


def cmd_signatures(cfg, store, a):
    since = iso(utcnow() - dt.timedelta(days=a.days))
    rows = store.conn.execute(
        "SELECT failure_signature s, failure_kind k, COUNT(*) n, GROUP_CONCAT(id) ids FROM jobs "
        "WHERE state = 'failed' AND finished_at >= ? GROUP BY s, k ORDER BY n DESC", (since,)).fetchall()
    for r in rows:
        print('{:4} [{}] {}\n     jobs {}'.format(r['n'], r['k'], r['s'], r['ids']))
    return 0


def cmd_dispatcher(cfg, store, a):
    from .dispatcher import AlreadyRunning, Dispatcher
    if a.log:
        # started detached (docker exec -d), so keep its output in the queue's log directory
        path = os.path.join(queue_paths(cfg)['logs'], 'dispatcher.log')
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o664)
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(fd, 1)
        os.dup2(fd, 2)
        os.close(fd)
    try:
        Dispatcher(cfg, store).run(once=a.once)
    except AlreadyRunning:
        print('a dispatcher already holds {}'.format(queue_paths(cfg)['lock']))
        return 1
    return 0


def cmd_watchdog(cfg, store, a):
    from .dispatcher import spawn_detached, watchdog
    rc = watchdog(cfg, store)
    if rc == 3 and a.spawn:
        pid = spawn_detached(cfg)
        print('{} started dispatcher pid {}'.format(iso(utcnow()), pid))
        return 0
    return rc


def cmd_stop(cfg, store, a):
    """SIGTERM the dispatcher (running jobs carry on). The watchdog starts a new one within 5 minutes unless
    its cron line is gone, so use this for a restart, or after removing the cron lines to stop the queue."""
    import signal
    from .resources import pid_alive
    info = json.loads(store.kv_get('dispatcher') or '{}')
    pid = info.get('pid')
    if not pid or not pid_alive(pid, info.get('proc_start')):
        print('no dispatcher is running')
        return 0
    os.kill(pid, signal.SIGTERM)
    print('sent SIGTERM to dispatcher pid {}; running jobs carry on in their own sessions'.format(pid))
    return 0


def _eso(cfg):
    from .eso import EsoArchive
    return EsoArchive(cfg['eso_env_file'])


def cmd_watch(cfg, store, a):
    from .watcher import Watcher
    now = utcnow()
    print('{} watcher ({} mode)'.format(iso(now), cfg['mode']))
    for line in Watcher(cfg, store, _eso(cfg), now=now).run():
        print('  ' + line)
    return 0


def cmd_lookback(cfg, store, a):
    from .lookback import Lookback, format_report, results_json
    now = utcnow()
    print('{} look-back over {} nights ({}{})'.format(iso(now), cfg['lookback']['nights'], cfg['mode'],
                                                     ', dry run' if a.dry_run else ''))
    res = Lookback(cfg, store, _eso(cfg), now=now, dry_run=a.dry_run).run()
    print(format_report(res))
    if a.json:
        with open(a.json, 'w') as f:
            f.write(results_json(res))
    return 0


def cmd_shadow_report(cfg, store, a):
    from .shadow_report import build_report
    text = build_report(cfg, store, days=a.days)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(text)
    print(text)
    return 0


def cmd_probe(cfg, store, a):
    from .dispatcher import deployed_sha
    from .resources import ProcResources, external_locks
    s = ProcResources(cfg['disk_devices']).snapshot()
    print('load5 {}  MemAvailable {:.0f} GB  cpus {}'.format(s.load5, s.mem_available_gb or 0, s.cpu_count))
    print('disk %util over 3 s: {}'.format(s.disk_util))
    print('deployed code: {}'.format(deployed_sha(cfg['src_dir'])))
    held = external_locks([])
    print('locks held by processes the queue did not start: {}'.format(len(held)))
    for lock, desc in sorted(held.items()):
        print('  {:34} {}'.format(lock, desc))
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog='python -m jobqueue', description='orchard job queue on appct')
    p.add_argument('--config', help='JSON config (default $ORCHARD_QUEUE_CONFIG, then <queue_root>/config.json)')
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('status', help='dispatcher, classes, running and queued jobs, recent failures')
    s.add_argument('--limit', type=int, default=25)
    s.add_argument('--resources', action='store_true', help='also sample load, memory and disk (3 s)')
    s.set_defaults(fn=cmd_status)

    s = sub.add_parser('show', help='one job and its events')
    s.add_argument('id', type=int)
    s.set_defaults(fn=cmd_show)

    s = sub.add_parser('add', help='queue a job (a pipeline night by default)')
    s.add_argument('--class', dest='cls', required=True, choices=CLASSES)
    s.add_argument('--kind', default='pipeline', help='pipeline (default) or a free-form label for a command job')
    s.add_argument('--telescope')
    s.add_argument('--night')
    s.add_argument('--targets', help='targets to run, space separated (default: every target of the night)')
    s.add_argument('--lock-targets', help='targets to lock (default: read from the night\'s frames)')
    s.add_argument('--no-T12', dest='no_T12', action='store_true', help='keep products in v3')
    s.add_argument('--flag', action='append', help='extra ZLP_pipeline.sh flag, e.g. --flag=--no_T7')
    s.add_argument('--cores', type=int)
    s.add_argument('--estimate', type=float, help='minutes (default: from the frame count)')
    s.add_argument('--priority', type=float, default=0.0, help='higher runs first within the class')
    s.add_argument('--depends-on', type=int)
    s.add_argument('--no-disk-heavy', action='store_true')
    s.add_argument('--lock', action='append', help='extra lock for a command job')
    s.add_argument('--cwd')
    s.add_argument('--note')
    s.add_argument('cmd', nargs=argparse.REMAINDER, help='-- command for a non-pipeline job')
    s.set_defaults(fn=cmd_add)

    s = sub.add_parser('pause', help='stop starting jobs of a class (or all)')
    s.add_argument('cls', choices=CLASSES + ('all',))
    s.add_argument('--reason')
    s.set_defaults(fn=cmd_pause)
    s = sub.add_parser('resume', help='undo pause')
    s.add_argument('cls', choices=CLASSES + ('all',))
    s.set_defaults(fn=cmd_resume)

    s = sub.add_parser('drain', help='start nothing new; running jobs finish (drain --off to undo)')
    s.add_argument('--off', action='store_true')
    s.set_defaults(fn=cmd_drain)
    s = sub.add_parser('undrain', help='same as drain --off')
    s.set_defaults(fn=cmd_drain, off=True)

    s = sub.add_parser('cancel', help='cancel queued jobs; a running job only with --kill')
    s.add_argument('ids', type=int, nargs='+')
    s.add_argument('--kill', action='store_true', help='also stop a running job (SIGTERM, then SIGKILL)')
    s.set_defaults(fn=cmd_cancel)

    s = sub.add_parser('requeue', help='requeue failed jobs, by id or by failure signature')
    s.add_argument('ids', type=int, nargs='*')
    s.add_argument('--signature', help='substring of the failure signature (see `signatures`)')
    s.add_argument('--exact', action='store_true')
    s.add_argument('--dry-run', action='store_true')
    s.set_defaults(fn=cmd_requeue)

    s = sub.add_parser('signatures', help='failed jobs grouped by failure signature')
    s.add_argument('--days', type=float, default=30)
    s.set_defaults(fn=cmd_signatures)

    s = sub.add_parser('dispatcher', help='run the dispatcher in the foreground (the watchdog starts it)')
    s.add_argument('--once', action='store_true')
    s.add_argument('--log', action='store_true', help='append output to <queue_root>/logs/dispatcher.log')
    s.set_defaults(fn=cmd_dispatcher)

    s = sub.add_parser('stop', help='stop the dispatcher (not the jobs); the watchdog restarts it if still in cron')
    s.set_defaults(fn=cmd_stop)

    s = sub.add_parser('watchdog', help='cron: exit 0 if the dispatcher is healthy, 3 if it must be started')
    s.add_argument('--spawn', action='store_true', help='start it here instead of leaving that to docker exec -d')
    s.set_defaults(fn=cmd_watchdog)

    s = sub.add_parser('watch', help='cron: one ESO watcher poll')
    s.set_defaults(fn=cmd_watch)

    s = sub.add_parser('lookback', help='cron: the daily 30-night look-back')
    s.add_argument('--dry-run', action='store_true', help='report only; queue nothing, record nothing')
    s.add_argument('--json', help='also write the per-night results here')
    s.set_defaults(fn=cmd_lookback)

    s = sub.add_parser('shadow-report', help='compare shadow decisions with what the cron actually ran')
    s.add_argument('--days', type=int, default=3)
    s.add_argument('--out')
    s.set_defaults(fn=cmd_shadow_report)

    s = sub.add_parser('probe', help='read-only: load, memory, disk %%util, code version, external locks')
    s.set_defaults(fn=cmd_probe)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'jobs':
        from .jobs import main as jobs_main
        return jobs_main(argv[1:])
    a = build_parser().parse_args(argv)
    cfg = load_config(a.config)
    store = _store(cfg)
    try:
        return a.fn(cfg, store, a) or 0
    finally:
        store.close()
