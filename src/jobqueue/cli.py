"""`q`: the queue's command line. Run inside orchard-server on appct:

    docker exec orchard-server python -m jobqueue status
    docker exec orchard-server python -m jobqueue add --class P2 --telescope Europa --night 20250101
    docker exec orchard-server python -m jobqueue add --class P2 --telescope Europa --night 20250101 --download
    docker exec orchard-server python -m jobqueue add --class P2 --from-file backlog.csv --dry-run
    docker exec orchard-server python -m jobqueue pause P2        # resume P2 | resume all
    docker exec orchard-server python -m jobqueue drain           # start nothing new; undrain to undo
    docker exec orchard-server python -m jobqueue cancel 42 [--kill]
    docker exec orchard-server python -m jobqueue requeue --signature "KeyError: 'gaia_dr3_id'"
    docker exec orchard-server python -m jobqueue add --class P2 --telescope Europa --night 20250101 --promote
    docker exec orchard-server python -m jobqueue promote list      # = python utils/promote.py list

src/jobqueue/cron/q is a host-side wrapper that does the docker exec for you.
"""
import argparse
import datetime as dt
import json
import os
import sys

from . import CLASSES, jobspec
from .config import ensure_dirs, load_config, queue_paths
from .util import iso, parse_iso, parse_night, utcnow


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
    from .releases import current, current_id
    release, rel = current_id(cfg), current(cfg)
    print('code       {}'.format('release {} (orchard {})'.format(release, (rel or {}).get('version') or '-')
                                 if release else 'live tree (no release)'))
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
        print('  {:>6} {} {:8} {:9} {:8} {:>3}c  {:>6} of est {:.0f}m  {}{}{}'.format(
            j['id'], j['class'], j['kind'], j['telescope'] or '-', j['night'] or '-', j['cores'],
            _ago(j['started_at'], now), j['est_minutes'], 'pid {}'.format(j['pid']) if j['pid'] else '',
            ' [shadow]' if j['shadow'] else '',
            '  code {}'.format(j['code_sha']) if release and j['code_sha'] != release else ''))
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


# --------------------------------------------------------------------------- nights: download and/or process

BACKLOG_FORMAT = 'TEL,DATE,DOWNLOAD,DELETE,PROCESS[,TARGETS]'


class EsoNights:
    """ESO rows of a night, fetched one telescope-month at a time (a backlog file lists many nights)."""

    def __init__(self, cfg, eso):
        self.cfg = cfg
        self.eso = eso
        self.cache = {}

    def rows(self, telescope, night):
        from .eso import group_by_night
        prog = self.cfg['telescopes'][telescope]['prog_id']
        month = night[:6]
        if (prog, month) not in self.cache:
            y, m = int(month[:4]), int(month[4:])
            last = dt.date(y + (m == 12), m % 12 + 1, 1) - dt.timedelta(days=1)
            self.cache[(prog, month)] = group_by_night(self.eso.frames(prog, month + '01', last.strftime('%Y%m%d')))
        return self.cache[(prog, month)].get(night, [])


def parse_backlog(path, cfg):
    """Requests from a file of TEL,DATE,DOWNLOAD,DELETE,PROCESS[,TARGETS] lines (manual_runs/backlog/v3_backlog.txt).
    DELETE=1 is taken as a refetch: the old files are kept in quarantine, never deleted. Returns (requests, errors)."""
    tels = {t.lower(): t for t in cfg['telescopes']}
    reqs, errors = [], []
    with open(path) as f:
        lines = f.read().splitlines()
    for no, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = [x.strip() for x in line.split(',')]
        if len(fields) < 5:
            errors.append('line {}: {!r}: expected {}'.format(no, line, BACKLOG_FORMAT))
            continue
        tel = tels.get(fields[0].lower())
        night = fields[1]
        flags = fields[2:5]
        bad = []
        if tel is None:
            bad.append('unknown telescope {!r}'.format(fields[0]))
        try:
            parse_night(night)
        except ValueError:
            bad.append('DATE {!r} is not YYYYMMDD'.format(night))
        if any(x not in ('0', '1') for x in flags):
            bad.append('DOWNLOAD, DELETE and PROCESS must be 0 or 1')
        if bad:
            errors.append('line {}: {!r}: {}'.format(no, line, '; '.join(bad)))
            continue
        download, delete, process = (x == '1' for x in flags)
        if delete and not download:
            errors.append('line {}: {!r}: DELETE needs DOWNLOAD'.format(no, line))
        elif download and cfg['telescopes'][tel].get('source') != 'eso':
            errors.append('line {}: {!r}: {} frames do not come from ESO'.format(no, line, tel))
        elif not (download or process):
            errors.append('line {}: {!r}: nothing to do'.format(no, line))
        else:
            targets = ' '.join(fields[5:]).replace('"', ' ').replace("'", ' ').split()
            reqs.append({'telescope': tel, 'night': night, 'download': download, 'process': process,
                         'replace': 'all' if delete else 'never', 'run_targets': targets, 'line': no})
    return reqs, errors


def queue_night(cfg, store, eso_nights, req, cls, priority=0.0, note=None, no_t12=False, extra_flags=(), cores=None,
                dry_run=False, log=print, campaign=None, promote=False, from_stage=None):
    """Queue one night: an add-only fetch (if req['download']) and the pipeline run after it (if req['process']).
    The pipeline job waits for the fetch and is cancelled if the fetch fails. promote: the pipeline runs with
    --no_T12 and a promote job follows it. from_stage: a rerun from that stage instead (queue_rerun)."""
    from .dispatcher import count_frames
    from .jobs import plan_fetch
    from .nights import night_dir, night_targets, normalise_target
    tel, night = req['telescope'], req['night']
    if from_stage:
        if req['download']:
            log('{} {} {}: not queued: a rerun from T{} reprocesses what is on disk; it does not download'.format(
                cls, tel, night, from_stage))
            return None
        queue_rerun(cfg, store, cls, tel, night, from_stage, req.get('run_targets'), priority=priority, note=note,
                    extra_flags=extra_flags, cores=cores, dry_run=dry_run, log=log, campaign=campaign, promote=promote)
        return None
    path = night_dir(cfg['basedir'], tel, night)
    fetch_id, lock_targets, what = None, [], []
    if req['download']:
        rows = eso_nights.rows(tel, night)
        todo, summary, _ = plan_fetch(rows, path, req['replace'])
        lock_targets = sorted({normalise_target(r['object']) for r in rows
                               if (r.get('dp_type') or '').upper() == 'OBJECT' and (r.get('object') or '').strip()})
        desc = 'ESO {} rows, on disk {} frames (+{} under other names), {} to fetch ({} science, {:.1f} GB{})'.format(
            summary['eso_rows'], summary['disk_files'], summary['foreign_files'], summary['todo'],
            summary['todo_science'], summary['todo_bytes'] / 1e9,
            '' if req['replace'] == 'never' else ', replace {}'.format(req['replace']))
        if not todo:
            what.append('nothing to fetch ({})'.format(desc))
        elif dry_run:
            what.append('would fetch: {}'.format(desc))
        else:
            fetch_id, created = jobspec.add(store, jobspec.fetch_job(
                cfg, cls, tel, night, n_rows=summary['todo'], n_science=summary['todo_science'], source='cli',
                rerun=False, replace=req['replace'], nbytes=summary['todo_bytes'], priority=priority,
                note=note or desc), campaign=campaign)
            what.append('{} fetch job {}: {}'.format('queued' if created else 'already queued as', fetch_id, desc))
    if req['process']:
        run_targets = list(req.get('run_targets') or [])
        scan = not lock_targets and not run_targets
        if dry_run:
            # no target scan here: for SSO it reads the header of every raw frame of the night, and a dry run of a
            # long backlog file once read 20 GB that way and held the queue's disk gate for an hour
            what.append('would run the pipeline{}{}{}{}'.format(
                ' after the fetch' if req['download'] else '',
                ' for ' + ' '.join(run_targets) if run_targets else '',
                ' (targets read from its frames when queued)' if scan else '',
                ' with --no_T12, then a promote job' if promote else ''))
        else:
            if scan:
                lock_targets = night_targets(cfg['basedir'], tel, night)
            jid, created, promote_id = jobspec.add_pipeline(store, cfg, jobspec.pipeline_job(
                cfg, cls, tel, night, lock_targets=lock_targets, run_targets=run_targets, frames=count_frames(path),
                source='cli', depends_on=fetch_id, dep_requires_success=fetch_id is not None, priority=priority,
                no_t12=no_t12 or promote, extra_flags=extra_flags, note=note, cores=cores), promote=promote,
                campaign=campaign)
            what.append('{} pipeline job {}{}{}'.format(
                'queued' if created else 'already queued as', jid, ' after job {}'.format(fetch_id) if fetch_id else '',
                ', then promote job {}'.format(promote_id) if promote_id else ''))
    log('{} {} {}: {}'.format(cls, tel, night, '; '.join(what)))
    return fetch_id


def queue_rerun(cfg, store, cls, tel, night, first, targets=None, priority=0.0, note=None, extra_flags=(), cores=None,
                dry_run=False, log=print, campaign=None, promote=False):
    """Queue a rerun of a night from stage `first` (docs/rerun-from-stage.md): a restage job (utils/restage.py), then
    ZLP_pipeline.sh --from_T<first> --no_T12 for the restaged targets, then, with promote, a promote job, which
    overlays v2 with what the rerun wrote where the night was already promoted. The targets default to every target
    of the night with the rerun's inputs in v3 or v2, and the pipeline job is always given them: its PDF then says
    which targets it covers. Returns the pipeline job id, or None."""
    from utils import restage
    from .dispatcher import count_frames
    from .nights import night_dir
    targets = list(targets or restage.night_targets(cfg['basedir'], tel, night, first))
    if not targets:
        log('{} {} {}: not queued: no target of the night has the inputs of a rerun from T{} in v3 or v2'.format(
            cls, tel, night, first))
        return None
    if dry_run:
        r = restage.Restager(cfg['basedir'], apply=False)
        _, lines = r.night(tel, night, first, targets)
        log('{} {} {}: would restage, then run ZLP_pipeline.sh --from_T{} --no_T12 for {}{}'.format(
            cls, tel, night, first, ' '.join(targets), ', then a promote job' if promote else ''))
        for line in lines:
            log('  ' + line)
        return None
    frames = count_frames(night_dir(cfg['basedir'], tel, night))
    spec = jobspec.pipeline_job(
        cfg, cls, tel, night, lock_targets=targets, run_targets=targets, frames=frames,
        est_minutes=jobspec.rerun_estimate(cfg, tel, night, frames, first), source='cli', priority=priority,
        no_t12=True, extra_flags=list(extra_flags or ()) + ['--from_T{}'.format(first)], note=note,
        cores=cores or min(jobspec.RERUN_CORES.get(first) or int(cfg['pipeline_cores']), int(cfg['pipeline_cores'])))
    rid, jid, created, pid = jobspec.add_rerun(store, cfg, spec, first, promote=promote, campaign=campaign)
    log('{} {} {}: restage job {}, then {} pipeline job {} (--from_T{}, {}){}'.format(
        cls, tel, night, rid, 'queued' if created else 'already queued as', jid, first, ' '.join(targets),
        ', then promote job {}'.format(pid) if pid else ''))
    return jid


def parse_from_stage(value):
    """--from: a stage (T10, 10) or a cause the stage map knows (pwv, teff, pdf: utils/stages.py CAUSES)."""
    from utils import stages
    v = value.strip()
    if v.upper().lstrip('T').isdigit():
        return int(v.upper().lstrip('T'))
    try:
        return stages.cause_stage(v)
    except KeyError as e:
        raise argparse.ArgumentTypeError(e.args[0])


def cmd_add_nights(cfg, store, a):
    if a.from_file:
        reqs, errors = parse_backlog(a.from_file, cfg)
        if errors:
            print('\n'.join(errors))
            print('{} bad line(s): nothing queued'.format(len(errors)))
            return 2
    else:
        if not (a.telescope and a.night):
            print('--download needs --telescope and --night')
            return 2
        if cfg['telescopes'].get(a.telescope, {}).get('source') != 'eso':
            print('{} frames do not come from ESO'.format(a.telescope))
            return 2
        reqs = [{'telescope': a.telescope, 'night': a.night, 'download': True, 'process': not a.download_only,
                 'replace': 'all' if a.refetch else 'invalid' if a.repair else 'never',
                 'run_targets': a.targets.split() if a.targets else []}]
    if a.repair or a.refetch:
        for r in reqs:
            if r['download'] and r['replace'] == 'never':
                r['replace'] = 'all' if a.refetch else 'invalid'
    if a.from_stage and any(r['download'] for r in reqs):
        print('--from reruns what is on disk: a line with DOWNLOAD=1 needs a full run; nothing queued')
        return 2
    eso_nights = EsoNights(cfg, _eso(cfg)) if any(r['download'] for r in reqs) else None
    for r in reqs:
        try:
            queue_night(cfg, store, eso_nights, r, a.cls, priority=a.priority, note=a.note, no_t12=a.no_T12,
                        extra_flags=a.flag or (), cores=a.cores, dry_run=a.dry_run, campaign=a.campaign,
                        promote=a.promote, from_stage=a.from_stage)
        except Exception as e:  # one bad night (ESO unreachable, say) must not stop the rest
            print('{} {} {}: not queued: {}: {}'.format(a.cls, r['telescope'], r['night'], type(e).__name__, e))
    return 0


def cmd_add_promote(cfg, store, a):
    """Promote jobs: one after each pipeline job of --after (finished ones included), or one for a night."""
    if a.after:
        for jid in sorted(_parse_ids(a.after)):
            job = store.get_job(jid)
            if job is None or job['kind'] != 'pipeline':
                print('job {}: not a pipeline job; nothing queued for it'.format(jid))
                continue
            pid, created = jobspec.add(store, jobspec.promote_job(
                cfg, a.cls, job['telescope'], job['night'], after=jid, targets=(job['targets'] or '').split(),
                source='cli', note=a.note or 'after pipeline job {}'.format(jid)),
                campaign=a.campaign or (job['meta'] or {}).get('campaign'))
            print('{} promote job {} after pipeline job {} ({} {}, {}){}'.format(
                'added' if created else 'already queued as', pid, jid, job['telescope'], job['night'], job['state'],
                '' if '--no_T12' in job['argv'] else '; note: that job runs T12, which sweeps v3 into v2 first'))
        return 0
    if not (a.telescope and a.night):
        print('a promote job needs --after JOB_IDS, or --telescope and --night')
        return 2
    pid, created = jobspec.add(store, jobspec.promote_job(cfg, a.cls, a.telescope, a.night, source='cli', note=a.note),
                               campaign=a.campaign)
    print('{} promote job {}: {} {} {} (every v3 target directory of the night, against its latest run)'.format(
        'added' if created else 'already queued as', pid, a.cls, a.telescope, a.night))
    return 0


def check_rerun_args(a):
    """What is wrong with a --from request, or ''."""
    from utils import stages
    if a.kind != 'pipeline':
        return '--from goes with pipeline jobs'
    if a.from_stage not in stages.RESTAGE:
        return ('a rerun from T{} is not supported: from {} only (utils/stages.py says why); for an earlier stage, '
                'rerun the whole night'.format(a.from_stage, ' or '.join('T{}'.format(n) for n in stages.RESTAGE)))
    if a.download or a.download_only:
        return '--from reruns what is on disk: it does not go with --download'
    if not (a.promote or a.no_T12):
        return ('--from needs --promote (the promotion overlays v2 with what the rerun wrote) or --no-T12 (keep it in '
                'v3): T12 would sweep the restaged inputs into v2')
    return ''


def cmd_add(cfg, store, a):
    if getattr(a, 'from_stage', None):
        why = check_rerun_args(a)
        if why:
            print(why)
            return 2
        if not a.from_file:
            if not (a.telescope and a.night):
                print('a rerun needs --telescope and --night, or --from-file')
                return 2
            queue_rerun(cfg, store, a.cls, a.telescope, a.night, a.from_stage, a.targets.split() if a.targets else None,
                        priority=a.priority, note=a.note, extra_flags=a.flag or (), cores=a.cores, dry_run=a.dry_run,
                        campaign=a.campaign, promote=a.promote)
            return 0
    if a.kind == 'promote':
        return cmd_add_promote(cfg, store, a)
    if a.after:
        print('--after goes with --kind promote')
        return 2
    if a.promote and a.kind != 'pipeline':
        print('--promote goes with pipeline jobs')
        return 2
    if a.from_file or a.download or a.download_only:
        return cmd_add_nights(cfg, store, a)
    if a.repair or a.refetch or a.dry_run:
        print('--repair, --refetch and --dry-run go with --download, --download-only or --from-file')
        return 2
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
            source='cli', depends_on=a.depends_on, priority=a.priority, no_t12=a.no_T12 or a.promote, note=a.note,
            cores=a.cores, extra_flags=a.flag or ())
        if a.promote:
            if a.no_disk_heavy:
                spec['disk_heavy'] = False
            jid, created, pid = jobspec.add_pipeline(store, cfg, spec, promote=True, campaign=a.campaign)
            print('{} job {}: {} pipeline {} {} with --no_T12, then promote job {}'.format(
                'added' if created else 'already queued as', jid, spec['cls'], spec['telescope'], spec['night'], pid))
            return 0
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
    jid, created = jobspec.add(store, spec, campaign=a.campaign)
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


# --------------------------------------------------------------------------- code releases

def _release_or_current(cfg, rid):
    from .releases import ReleaseError, current_id
    rid = rid or current_id(cfg)
    if not rid:
        raise ReleaseError('no release given and none is current')
    return rid


def _print_blocked(blocked):
    for pid, cmd in sorted(blocked.items()):
        print('  pid {} {}'.format(pid, cmd))


def _after_sync(cfg, store, a, res):
    from .releases import queue_code_changed
    for path in res['written']:
        print('  wrote ' + path)
    if res.get('backup'):
        print('replaced files are backed up in {}'.format(res['backup']))
    if queue_code_changed(res):
        print('the queue\'s own code changed: restarting the dispatcher (the watchdog starts the new one within 5 min)')
        cmd_stop(cfg, store, a)


def _deploy_sync(cfg, store, a):
    """After activating: bring the live tree up now if nothing runs from it, else leave it to the dispatcher."""
    from .releases import sync_pending
    res = sync_pending(cfg, store)
    if res is None:
        return 0
    if 'blocked' in res:
        print('the live tree is still used by {} process(es) started before releases or outside the queue; the '
              'dispatcher brings it up to date once they end (q release shows when):'.format(len(res['blocked'])))
        _print_blocked(res['blocked'])
        return 0
    _after_sync(cfg, store, a, res)
    print('live tree matches {}'.format(res['release']))
    return 0


def cmd_release(cfg, store, a):
    from .releases import (SYNC_KEY, ReleaseError, activate, build, current_id, list_releases, live_diff, run_tests,
                           sync_live)
    try:
        if a.action == 'list':
            cur = current_id(cfg)
            running = {}
            for j in store.jobs(states=['running']):
                running.setdefault(j['code_sha'], []).append(j['id'])
            print('current: {}'.format(cur or 'none (jobs run from the live tree)'))
            releases = list_releases(cfg)
            for r in releases[-a.limit:]:
                print('{} {:28} {:>8}  orchard {:18} built {}  {:>3} changed  {:>3} differed from live{}{}'.format(
                    '*' if r['id'] == cur else ' ', r['id'], (r.get('sha') or '-')[:8], r.get('version') or '-',
                    r['built_at'],
                    len(r.get('changed_since_previous') or []), len(r.get('differs_from_live') or []),
                    '  running jobs {}'.format(running[r['id']]) if r['id'] in running else '',
                    '  tests {}'.format(r['tests']) if r.get('tests') else ''))
            known = {r['id'] for r in releases}
            for sha, ids in sorted(((k, v) for k, v in running.items() if k not in known), key=lambda kv: str(kv[0])):
                print('  live tree {:20} running jobs {}'.format(sha or '-', ids))
            if store.kv_get(SYNC_KEY):
                print('live tree sync pending: the dispatcher does it once nothing runs from the live tree')
            if cur:
                diff = live_diff(cfg, cur)
                print('live tree {}: {}'.format(cfg['src_dir'], 'matches {}'.format(cur) if not diff else
                                                '{} file(s) differ from {}{}'.format(
                                                    len(diff), cur, ': ' + ', '.join(diff[:8]) +
                                                    (' ...' if len(diff) > 8 else ''))))
            return 0
        if a.action == 'build':
            rel = build(cfg, stage=a.stage, sha=a.sha)
            for path in rel['differs_from_live']:
                print('  differs from live: ' + path)
            print('next: q release test {0}; q release activate {0}'.format(rel['id']))
            return 0
        if a.action == 'test':
            return 1 if run_tests(cfg, _release_or_current(cfg, a.id)) else 0
        if a.action == 'activate':
            activate(cfg, a.id, store=store)
            return _deploy_sync(cfg, store, a)
        if a.action == 'sync-live':
            rid = _release_or_current(cfg, a.id)
            res = sync_live(cfg, rid, dry_run=a.dry_run, force=a.force)
            if 'blocked' in res:
                print('live tree not updated: {} process(es) still run from it:'.format(len(res['blocked'])))
                _print_blocked(res['blocked'])
                return 2
            if a.dry_run:
                print('{} file(s) of {} differ in the live tree'.format(len(res['written']), rid))
                for path in res['written']:
                    print('  ' + path)
                return 0
            if rid == current_id(cfg):
                store.kv_delete(SYNC_KEY)
            _after_sync(cfg, store, a, res)
            print('live tree matches {}'.format(rid))
            return 0
        if a.action == 'deploy':
            rel = build(cfg, stage=a.stage, sha=a.sha)
            for path in rel['differs_from_live']:
                print('  differs from live: ' + path)
            if not rel['differs_from_live'] and a.stage:
                print('note: nothing in {} differs from the live tree'.format(a.stage))
            if a.skip_tests:
                print('tests skipped')
            elif run_tests(cfg, rel['id']):
                print('tests failed: {} is built but NOT active; read {}/tests.log, then '
                      'q release activate {} if the failures are known'.format(rel['id'], rel['path'], rel['id']))
                return 1
            activate(cfg, rel['id'], store=store)
            return _deploy_sync(cfg, store, a)
    except ReleaseError as e:
        print('error: {}'.format(e))
        return 1
    raise ValueError(a.action)



def _eso(cfg):
    from .eso import EsoArchive
    return EsoArchive(cfg['eso_env_file'])


def _eso_or_record(cfg, store, now, source, record=True):
    """The ESO client, or None after recording a failed poll: a client that cannot even be set up (credentials
    unreadable, the downloader failing to import) is ESO unreachable as far as the nightly work is concerned."""
    from . import esoreach
    try:
        return _eso(cfg)
    except Exception as e:                       # noqa: BLE001 - reported and recorded, like a failed query
        print('{} cannot set up the ESO client: {}'.format(iso(now), e))
        if record:
            esoreach.record(store, now, ok=False, error='cannot set up the ESO client: {}'.format(e), source=source)
        return None


def cmd_watch(cfg, store, a):
    from .watcher import Watcher
    now = utcnow()
    print('{} watcher ({} mode)'.format(iso(now), cfg['mode']))
    eso = _eso_or_record(cfg, store, now, 'watcher')
    if eso is None:
        return 1
    for line in Watcher(cfg, store, eso, now=now).run():
        print('  ' + line)
    return 0


def cmd_lookback(cfg, store, a):
    from .lookback import Lookback, format_report, results_json
    now = utcnow()
    print('{} look-back over {} nights ({}{})'.format(iso(now), cfg['lookback']['nights'], cfg['mode'],
                                                     ', dry run' if a.dry_run else ''))
    eso = _eso_or_record(cfg, store, now, 'look-back', record=not a.dry_run)
    if eso is None:
        return 1
    res = Lookback(cfg, store, eso, now=now, dry_run=a.dry_run).run()
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


def _parse_ids(spec):
    ids = set()
    for part in (spec or '').split(','):
        part = part.strip()
        if '-' in part:
            first, last = part.split('-', 1)
            ids.update(range(int(first), int(last) + 1))
        elif part:
            ids.add(int(part))
    return ids


def cmd_campaign(cfg, store, a):
    """Label jobs already queued (or finished) with a campaign, by id and/or by a note substring."""
    if not (a.ids or a.note):
        print('give --ids and/or --note to choose the jobs')
        return 2
    if not (a.label or a.clear):
        print('give a label, or --clear')
        return 2
    ids = _parse_ids(a.ids)
    text = (a.note or '').lower()
    rows = [r for r in store.conn.execute('SELECT id, meta, note FROM jobs ORDER BY id')
            if (not ids or r['id'] in ids) and (not text or text in (r['note'] or '').lower())]
    for r in rows:
        meta = json.loads(r['meta'] or '{}')
        if a.clear:
            meta.pop('campaign', None)
        else:
            meta['campaign'] = a.label
        if not a.dry_run:
            store.update_job(r['id'], meta=meta)
    print('{} {} job(s){}{}'.format('would update' if a.dry_run else 'updated', len(rows),
                                    ': campaign cleared' if a.clear else ': campaign "{}"'.format(a.label),
                                    ' ({}-{})'.format(rows[0]['id'], rows[-1]['id']) if rows else ''))
    return 0


def cmd_notify(cfg, store, a):
    from .notify import Notifier, NotifyError, format_incidents
    try:
        n = Notifier(cfg, store, outbox_only=getattr(a, 'outbox', False))
        if a.action == 'digest':
            print(n.digest(night=a.night, dry_run=a.print))
        elif a.action == 'check':
            print('\n'.join(n.check(dry_run=a.print)))
        elif a.action == 'incidents':
            rows = n.incidents()
            if a.all:
                rows += list(reversed(n.incidents(state='resolved', limit=30)))
            print(format_incidents(rows, n.now))
        elif a.action == 'ack':
            r = n.ack(a.id, a.until, a.note)
            print('incident {} ({}) acknowledged until {}'.format(r['id'], r['title'], a.until))
        elif a.action == 'unack':
            r = n.unack(a.id)
            print('incident {} ({}) no longer acknowledged'.format(r['id'], r['title']))
        elif a.action == 'messages':
            for m in n.messages(a.last):
                print('{}  {:<8} {:<7} {}  {}{}'.format(m['ts'], m['kind'], m['status'], m['subject'], m['path'] or '',
                                                       '  ({})'.format(m['error']) if m['error'] else ''))
        elif a.action == 'recipients':
            from .notify import config_path, edit_recipients
            if a.add or a.remove:
                to, path = edit_recipients(cfg, add=a.add or (), remove=a.remove or ())
            else:
                to, path = list(n.delivery.to or cfg['notify'].get('to') or []), config_path(cfg)
            print('recipients ({}): {}'.format(path, ', '.join(to) or 'none'))
            print('from: {}; transport: {}{}'.format(
                n.delivery.sender, cfg['notify']['transport'],
                '' if cfg['notify']['transport'] == 'smtp' else ' (nothing is sent until it is smtp)'))
    except (NotifyError, ValueError) as e:
        print('notify: {}'.format(e), file=sys.stderr)
        return 2
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
    s.add_argument('--kind', default='pipeline', help='pipeline (default), promote (utils/promote.py for a night, or '
                                                       'after --after jobs), or a free-form label for a command job')
    s.add_argument('--telescope')
    s.add_argument('--night')
    s.add_argument('--targets', help='targets to run, space separated (default: every target of the night)')
    s.add_argument('--lock-targets', help='targets to lock (default: read from the night\'s frames)')
    s.add_argument('--no-T12', dest='no_T12', action='store_true', help='keep products in v3')
    s.add_argument('--promote', action='store_true',
                   help='run the pipeline with --no_T12 and queue a promote job after it, which publishes into v2 '
                        'what passes promote_check and keeps the rest for review (utils/promote.py)')
    s.add_argument('--after', help='with --kind promote: pipeline job ids (1544-1670,1700) to promote after')
    s.add_argument('--flag', action='append', help='extra ZLP_pipeline.sh flag, e.g. --flag=--no_T7')
    s.add_argument('--cores', type=int)
    s.add_argument('--estimate', type=float, help='minutes (default: from the frame count)')
    s.add_argument('--priority', type=float, default=0.0, help='higher runs first within the class')
    s.add_argument('--depends-on', type=int)
    s.add_argument('--no-disk-heavy', action='store_true')
    s.add_argument('--lock', action='append', help='extra lock for a command job')
    s.add_argument('--cwd')
    s.add_argument('--note')
    s.add_argument('--campaign', help='short label the morning digest groups this work under, e.g. "WCS audit reruns"')
    s.add_argument('--download', action='store_true',
                   help='fetch what ESO holds for the night and the disk lacks (add-only), then run the pipeline')
    s.add_argument('--download-only', action='store_true', help='as --download, without the pipeline run')
    s.add_argument('--repair', action='store_true',
                   help='also re-fetch frames on disk that fail verification (truncated); old files go to quarantine')
    s.add_argument('--refetch', action='store_true',
                   help='fetch the whole night again and swap in every frame that differs; old files go to quarantine')
    s.add_argument('--from-file', help='queue each line of a {} file (the v3_backlog.txt format); DELETE=1 means '
                                       '--refetch'.format(BACKLOG_FORMAT))
    s.add_argument('--from', dest='from_stage', type=parse_from_stage,
                   help='rerun from this stage, or for this cause (T10 or pwv/teff/light-curves: T10 and T11; T11 or '
                        'pdf: the PDF): a restage job puts back into v3 what the rerun reads (from v2 for a promoted '
                        'night), then ZLP_pipeline.sh --from_T<N> runs; with --promote the promotion overlays v2 with '
                        'what the rerun wrote (docs/rerun-from-stage.md). Needs --promote or --no-T12')
    s.add_argument('--dry-run', action='store_true',
                   help='with --download, --from-file or --from: only say what would be queued')
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

    for name, text in (('report', 'daily report: per night, when it was ready, its jobs and results (live), or '
                                  'shadow decisions next to what the cron ran (shadow)'),
                       ('shadow-report', 'same as report')):
        s = sub.add_parser(name, help=text)
        s.add_argument('--days', type=int, default=3)
        s.add_argument('--out')
        s.set_defaults(fn=cmd_shadow_report)

    s = sub.add_parser('release', help='code releases: every job runs from a frozen copy of the deployed tree, '
                                       'so deploying never waits for the queue (see docs/job-queue.md)')
    rs = s.add_subparsers(dest='action')
    s.set_defaults(fn=cmd_release, action='list', limit=10)
    r = rs.add_parser('list', help='releases, the current one, what runs on each, and whether the live tree matches')
    r.add_argument('--limit', type=int, default=10)
    for name, text in (('build', 'build a release from the live tree with a staged commit laid over it'),
                       ('deploy', 'build, test, activate, then bring the live tree up to it')):
        r = rs.add_parser(name, help=text)
        r.add_argument('--from', dest='stage', help='src/ of `git archive <sha> src`, as the container sees it '
                                                   '(none: snapshot the live tree as it is)')
        r.add_argument('--sha', help='the staged commit')
        if name == 'deploy':
            r.add_argument('--skip-tests', action='store_true')
    r = rs.add_parser('test', help='run the test suite from a release (default the current one)')
    r.add_argument('id', nargs='?')
    r = rs.add_parser('activate', help='make a release current: jobs started from now on run it, and the live '
                                       'tree follows as soon as nothing runs from it (also a rollback)')
    r.add_argument('id')
    r = rs.add_parser('sync-live', help='bring the live tree up to a release now (default the current one; the '
                                        'dispatcher does this itself after activate): writes what differs, backs '
                                        'up what it replaces, deletes nothing')
    r.add_argument('id', nargs='?')
    r.add_argument('--dry-run', action='store_true', help='list the files that differ')
    r.add_argument('--force', action='store_true',
                   help='write even while pipeline processes run from the live tree (they may mix old and new code)')

    s = sub.add_parser('campaign', help='label jobs already in the queue with a campaign, for the morning digest')
    s.add_argument('label', nargs='?', help='e.g. "WCS audit reruns"')
    s.add_argument('--ids', help='job ids: 1200-1350,1400')
    s.add_argument('--note', help='jobs whose note contains this text (case-insensitive)')
    s.add_argument('--clear', action='store_true', help='remove the label instead')
    s.add_argument('--dry-run', action='store_true', help='only count the jobs it would change')
    s.set_defaults(fn=cmd_campaign)

    s = sub.add_parser('notify', help='the notifier: morning digest and incident alerts, the only email the '
                                      'pipeline sends (see docs/job-queue.md)')
    ns = s.add_subparsers(dest='action')
    s.set_defaults(fn=cmd_notify, action='incidents', all=False)
    for name, text in (('digest', 'cron, each morning: bring incidents up to date, deliver their alerts, then the '
                                  "digest of last night"),
                       ('check', 'cron, every 30 minutes: bring incidents up to date and deliver their alerts')):
        r = ns.add_parser(name, help=text)
        if name == 'digest':
            r.add_argument('--night', help='YYYYMMDD instead of last night')
        r.add_argument('--print', action='store_true', help='show what it would say; change and deliver nothing')
        r.add_argument('--outbox', action='store_true', help='only write to the outbox, even if notify.transport '
                                                             'is smtp')
    r = ns.add_parser('incidents', help='open incidents (--all: the last 30 resolved too)')
    r.add_argument('--all', action='store_true')
    r = ns.add_parser('ack', help="silence an incident's alert until a date; it stays in the digest")
    r.add_argument('id', help='incident id, or the key of an open incident')
    r.add_argument('--until', required=True, help='YYYY-MM-DD')
    r.add_argument('--note', default='')
    r = ns.add_parser('unack', help='undo ack')
    r.add_argument('id')
    r = ns.add_parser('messages', help='messages delivered, newest first, with their outbox files')
    r.add_argument('--last', type=int, default=20)
    r = ns.add_parser('recipients', help='the recipient list (notify.to in the queue config): show, --add, --remove')
    r.add_argument('--add', action='append', metavar='ADDRESS')
    r.add_argument('--remove', action='append', metavar='ADDRESS')

    s = sub.add_parser('probe', help='read-only: load, memory, disk %%util, code version, external locks')
    s.set_defaults(fn=cmd_probe)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'jobs':
        from .jobs import main as jobs_main
        return jobs_main(argv[1:])
    if argv and argv[0] == 'promote':
        from utils.promote import main as promote_main
        return promote_main(argv[1:])
    a = build_parser().parse_args(argv)
    cfg = load_config(a.config)
    store = _store(cfg)
    try:
        return a.fn(cfg, store, a) or 0
    finally:
        store.close()
