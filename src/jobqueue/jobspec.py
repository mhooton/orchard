"""How the watcher, the look-back and `q add` describe pipeline, download and fetch jobs.

Locks:
  night:<TEL>:<NIGHT>   every job that reads or writes that night's raw or output directories
  tel:<TEL>             pipeline jobs, while ZLP_pipeline.sh T8/T9 `rm -rf ${DATDIR}/catcache` share one
                        directory per telescope (config telescope_lock; off once a per-run catcache lands)
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
    locks = {'night:{}:{}'.format(telescope, night)}
    if cfg.get('telescope_lock', True):
        locks.add('tel:' + telescope)
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


def add(store, spec):
    spec = dict(spec)
    return store.add_job(spec.pop('cls'), spec.pop('kind'), spec.pop('argv'), **spec)
