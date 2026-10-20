"""How the watcher, the look-back and `q add` describe pipeline, download, fetch and promote jobs.

Locks:
  night:<TEL>:<NIGHT>   every job that reads or writes that night's raw or output directories
  tel:<TEL>             pipeline jobs; enforced only while config telescope_lock is set (one job per telescope)
  target:<NAME>         pipeline jobs, one per target: v2/StackImages is shared across telescopes and nights
  sem:eso               download and fetch jobs, counted against semaphores.eso
  sem:eso_bulk          P2/P3 fetch jobs as well, counted against semaphores.eso_bulk (1), so bulk downloads
                        never take both ESO slots from the nightly ones
"""
import os

from .resources import target_lock
from .dispatcher import pipeline_estimate


def job_env(cfg):
    env = {'ORCHARD_QUEUE_ROOT': cfg['queue_root']}
    if cfg.get('_config_path'):
        env['ORCHARD_QUEUE_CONFIG'] = cfg['_config_path']
    return env


def pipeline_locks(cfg, telescope, night, targets):
    locks = {'night:{}:{}'.format(telescope, night), 'tel:' + telescope}
    locks.update(target_lock(t) for t in targets)
    return sorted(locks)


def pipeline_argv(cfg, telescope, night, run_targets=(), no_t12=False, extra_flags=()):
    tel_cfg = cfg['telescopes'].get(telescope, {})
    argv = ['./main/ZLP_pipeline.sh']
    if tel_cfg.get('force_platesolve'):
        argv.append('--force-platesolve')
    if no_t12:
        argv.append('--no_T12')
    argv += list(extra_flags)
    argv += [cfg['pipeline_runname'], cfg['basedir'], night, cfg['pipeline_cthresh'], cfg['pipeline_sthresh'],
             telescope]
    if run_targets:
        argv.append(' '.join(run_targets))  # ZLP_pipeline.sh takes the targets as one space-separated argument
    return argv


def pipeline_job(cfg, cls, telescope, night, lock_targets=(), run_targets=(), frames=None, est_minutes=None,
                 source='manual', depends_on=None, priority=0.0, no_t12=False, extra_flags=(), note=None, meta=None,
                 cores=None, max_retries=None, dep_requires_success=False):
    if est_minutes is None:
        est_minutes = pipeline_estimate(cfg, telescope, night, frames or 0)
    lock_targets = sorted(set(lock_targets) | set(run_targets))
    meta = dict(meta or {})
    meta.update({'lock_targets': lock_targets, 'frames_at_enqueue': frames})
    return dict(
        cls=cls, kind='pipeline', telescope=telescope, night=night, targets=' '.join(run_targets),
        argv=pipeline_argv(cfg, telescope, night, run_targets, no_t12, extra_flags), cwd=cfg['src_dir'],
        env=job_env(cfg), cores=cores or cfg['pipeline_cores'], disk_heavy=True, est_minutes=est_minutes,
        locks=pipeline_locks(cfg, telescope, night, lock_targets), depends_on=depends_on,
        dep_requires_success=dep_requires_success, priority=priority,
        dedupe_key='pipeline:{}:{}:{}'.format(telescope, night, ' '.join(run_targets) or '*'),
        max_retries=cfg['max_retries'] if max_retries is None else max_retries, source=source, meta=meta, note=note)


def download_job(cfg, telescope, night, eso_rows, source='watcher', note=None):
    """SSO_download.py for one recent night, the way the nightly cron does it (see jobs.download)."""
    argv = [cfg['python'], '-m', 'jobqueue', 'jobs', 'download', '--telescope', telescope, '--night', night]
    return dict(
        cls='P0', kind='download', telescope=telescope, night=night, argv=argv, cwd=cfg['src_dir'],
        env=job_env(cfg), cores=cfg['download_cores'], disk_heavy=False,
        est_minutes=max(10.0, eso_rows / 100.0), locks=['night:{}:{}'.format(telescope, night), 'sem:eso'],
        dedupe_key='download:{}:{}'.format(telescope, night), max_retries=cfg['max_retries'], source=source,
        meta={'eso_rows': eso_rows}, note=note)


def fetch_estimate(cfg, n_rows, nbytes=None):
    """Minutes for a fetch: ESO's compressed size at fetch.mb_per_s, plus unpacking and placing."""
    if nbytes:
        return max(10.0, nbytes / (float(cfg['fetch']['mb_per_s']) * 1e6) / 60.0 + n_rows * 0.02)
    return max(10.0, n_rows / 60.0)


def fetch_job(cfg, cls, telescope, night, rows_file=None, n_rows=0, n_science=0, source='lookback', rerun=True,
              replace='never', nbytes=None, est_minutes=None, priority=0.0, note=None):
    """Add-only fetch of one night (see jobs.fetch and fetch.py): the listed rows, or, without rows_file, whatever
    ESO holds that is not on disk. rerun=False: download only, never queue a pipeline run."""
    argv = [cfg['python'], '-m', 'jobqueue', 'jobs', 'fetch', '--telescope', telescope, '--night', night]
    if rows_file:
        argv += ['--rows', rows_file]
    argv += ['--rerun-class', cls if rerun else 'none']
    if replace != 'never':
        argv += ['--replace', replace]
    locks = ['night:{}:{}'.format(telescope, night), 'sem:eso']
    if cls in ('P2', 'P3'):
        locks.append('sem:eso_bulk')
    return dict(
        cls=cls, kind='fetch', telescope=telescope, night=night, argv=argv, cwd=cfg['src_dir'], env=job_env(cfg),
        cores=cfg['download_cores'], disk_heavy=True,
        est_minutes=est_minutes or fetch_estimate(cfg, n_rows, nbytes), locks=locks,
        dedupe_key='fetch:{}:{}'.format(telescope, night), max_retries=cfg['max_retries'], source=source,
        priority=priority, note=note,
        meta={'rows': n_rows, 'science': n_science, 'rows_file': os.path.basename(rows_file) if rows_file else None,
              'bytes': nbytes, 'replace': replace, 'rerun': rerun})


def promote_job(cfg, cls, telescope, night, after=None, targets=(), source='cli', note=None, meta=None):
    """utils/promote.py for one night (see its docstring): it judges what the run left in v3 with promote_check and
    publishes what passes into v2, keeping holds and rejects for a person. after: the pipeline job it follows,
    whether that job succeeds or fails (a cancelled one cancels it); it then judges that job's targets. Without
    `after` it judges every v3 target directory of the night against the night's latest run."""
    pcfg = cfg['promote']
    argv = [cfg['python'], 'utils/promote.py', 'run', '--apply', '--telescope', telescope, '--night', night]
    if after:
        argv += ['--job', str(after)]
    # one BLAS/OpenMP thread, as ZLP_pipeline.sh does for its report: numpy otherwise starts a pool sized to the host
    # (64 threads) in every promote job, and many processes at once can exhaust RLIMIT_NPROC
    env = dict(job_env(cfg), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    return dict(
        cls=cls, kind='promote', telescope=telescope, night=night, targets=' '.join(targets), argv=argv,
        cwd=cfg['src_dir'], env=env, cores=int(pcfg['cores']), disk_heavy=False,
        est_minutes=float(pcfg['est_minutes']), locks=['night:{}:{}'.format(telescope, night)], depends_on=after,
        dep_requires_success=False, priority=float(pcfg['priority']),
        dedupe_key='promote:{}:{}:{}'.format(telescope, night, after or '*'), max_retries=cfg['max_retries'],
        source=source, meta=dict(meta or {}, after=after), note=note)


# A rerun from a late stage costs this share of a full run's estimate (minutes per 1000 frames): T10 is six
# SPlightcurve.py runs per target, T11 the PDF (which aligns every frame for its stack panel). The scratch trial of
# 2026-10-20 (docs/rerun-from-stage.md) took 10-32% of the original full runs for T10 + T11, 2-7% for T11.
RERUN_SHARE = {10: 0.3, 11: 0.1}
# A rerun from T10 reserves fewer cores than a full run: SPlightcurve.py is one process, and only T11's alignment of the
# frames for its stack panel runs in parallel (the trial ran with 8)
RERUN_CORES = 8


def rerun_estimate(cfg, telescope, night, frames, first):
    return max(float(cfg['min_estimate_minutes']), pipeline_estimate(cfg, telescope, night, frames or 0) *
               RERUN_SHARE.get(first, 1.0))


def restage_job(cfg, cls, telescope, night, first, targets=(), source='cli', note=None, priority=0.0):
    """utils/restage.py for one night (docs/rerun-from-stage.md): puts back into v3 what a rerun from stage `first`
    reads (from v2, as hard links, where promotion already moved the night) and records it, so that the promotion
    after the rerun overlays v2 with only what the rerun wrote. Holds the night lock; links and a few small copies."""
    argv = [cfg['python'], 'utils/restage.py', 'stage', '--apply', '--telescope', telescope, '--night', night,
            '--from', str(first)]
    for t in targets:
        argv += ['--target', t]
    return dict(
        cls=cls, kind='restage', telescope=telescope, night=night, targets=' '.join(targets), argv=argv,
        cwd=cfg['src_dir'], env=job_env(cfg), cores=1, disk_heavy=False, est_minutes=5.0,
        locks=['night:{}:{}'.format(telescope, night)], priority=priority,
        dedupe_key='restage:{}:{}:{}:T{}'.format(telescope, night, ' '.join(targets) or '*', first),
        max_retries=cfg['max_retries'], source=source, note=note, meta={'from_stage': first})


def add_rerun(store, cfg, spec, first, promote=False, campaign=None):
    """Queue a rerun from stage `first` of a night: a restage job, then the pipeline job described by pipeline_job()
    with --from_T<first>, which runs only if the restage succeeded, and, with promote, a promote job after it.
    Returns (restage job id, pipeline job id, created, promote job id or None)."""
    flag = '--from_T{}'.format(first)
    if flag not in spec['argv']:
        raise ValueError('a rerun from T{} must run ZLP_pipeline.sh with {}'.format(first, flag))
    if '--no_T12' not in spec['argv']:
        raise ValueError('a rerun from T{} must run with --no_T12: T12 would sweep the restaged inputs into v2'.format(
            first))
    rid, _ = add(store, restage_job(cfg, spec['cls'], spec['telescope'], spec['night'], first,
                                    (spec.get('targets') or '').split(), source=spec.get('source') or 'cli',
                                    note=spec.get('note'), priority=spec.get('priority') or 0.0), campaign=campaign)
    spec = dict(spec, depends_on=rid, dep_requires_success=True,
                dedupe_key='{}:from_T{}'.format(spec['dedupe_key'], first))
    jid, created, pid = add_pipeline(store, cfg, spec, promote=promote, campaign=campaign)
    return rid, jid, created, pid


def nightly_promote(cfg):
    """Whether the watcher's, the look-back's and fetch reruns' pipeline jobs run with --no_T12 and a promote job."""
    return bool((cfg.get('promote') or {}).get('nightly'))


def add_pipeline(store, cfg, spec, promote=False, campaign=None):
    """Queue a pipeline job described by pipeline_job() and, with promote, a promote job after it. Returns
    (pipeline job id, created, promote job id or None). A promoted run must not run T12, which would sweep v3 into
    v2 before anything is judged."""
    if promote and '--no_T12' not in spec['argv']:
        raise ValueError('a pipeline job followed by a promote job must run with --no_T12')
    jid, created = add(store, spec, campaign=campaign)
    pid = None
    if promote:
        pid, _ = add(store, promote_job(cfg, spec['cls'], spec['telescope'], spec['night'], after=jid,
                                        targets=(spec.get('targets') or '').split(), source=spec.get('source') or 'cli',
                                        note='after pipeline job {}'.format(jid)), campaign=campaign)
    return jid, created, pid


def add(store, spec, campaign=None):
    """Queue a described job. campaign: a short label (`q add --campaign`) the morning digest groups work by; it is
    kept in the job's meta."""
    spec = dict(spec)
    if campaign:
        spec['meta'] = dict(spec.get('meta') or {}, campaign=campaign)
    return store.add_job(spec.pop('cls'), spec.pop('kind'), spec.pop('argv'), **spec)
