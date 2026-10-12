"""Job bodies that are not the pipeline itself: the nightly download and the add-only fetch.

download  SSO_download.py for one recent night, exactly as the nightly cron runs it (log to
          ESO_logs/<NIGHT><TEL>.log, its usual summary email), with --max-retries 1 because the watcher
          only queues a night once ESO has it. Refuses nights older than watcher.recent_days, and refuses
          (falling back to the add-only fetch) if transformation_check() would delete frames already on
          disk. Exits 0 whenever frames are on disk, so the pipeline job behind it runs even if a few are
          missing (the look-back fetches those); exits 75 (retry later) if nothing landed.
fetch     add-only fetch (fetch.py) of any night: the rows the look-back listed (--rows), or, without a list,
          whatever ESO holds for the night that is not on disk. With --rerun-class P1 (the look-back) it queues
          a rerun if science frames were added; with --rerun-class none it only downloads (`q add --download`
          queues the pipeline run itself, as a job that waits for this one). --replace invalid also swaps
          frames on disk that fail verification (truncated); --replace all fetches the whole night again and
          swaps every frame that differs. Swapped files are kept in <queue_root>/quarantine/<TEL>/<NIGHT>.
          Every file added or swapped goes into <queue_root>/fetch/manifest.csv.
"""
import argparse
import datetime as dt
import os
import subprocess
import sys
import threading

from . import jobspec
from .config import ensure_dirs, load_config, queue_paths
from .nights import (FitsError, diff_eso_disk, foreign_times, is_cube_layout, list_frames, night_dir, read_header,
                     row_files, science_count, transformation_check_kind, transformation_check_would_delete,
                     verify_frame)
from .util import parse_night, read_json, utcnow, write_json_atomic

EX_TEMPFAIL = 75


def _archive(cfg):
    from .eso import EsoArchive
    return EsoArchive(cfg['eso_env_file'])


def eso_token_fn(archive):
    """token(force=False) for EsoFileFetcher: the archive's bearer token, renewed when due or when forced."""
    lock = threading.Lock()

    def token(force=False):
        with lock:
            if force:
                archive._d.token = None
            archive._ensure_token()
            return archive._d.token
    return token


def _verify_placed(path):
    """What a frame must pass before it goes into a night directory: complete FITS, and transformed if Astra."""
    from .fetch import _needs_astra
    verify_frame(path)
    if _needs_astra(path):
        raise FitsError('an Astra frame without ASTRAROT/ASTRAMIR: the transform did not run')


def _is_cube(path):
    try:
        return is_cube_layout(verify_frame(path))
    except (OSError, FitsError):
        return False


def fetch_rows(cfg, archive, telescope, night, rows, job_id=None, dest=None, replace='never', foreign=None,
               record=None, log=print):
    """Fetch ESO rows add-only into the night directory (or dest, for a trial). Records every file placed in
    the manifest unless dest is given."""
    from download.unpack_datacubes import unpack_datacube
    from download.astra_transform import astra_transform
    from .fetch import EsoFileFetcher, Manifest, add_only_fetch

    def unpack(path, dest):
        ids, _ = unpack_datacube(path, dest)
        return [os.path.join(dest, i + '.fits') for i in ids]

    fcfg = cfg['fetch']
    paths = queue_paths(cfg)
    job_id = job_id or os.getenv('JOBQUEUE_JOB_ID')
    if record is None and not dest:
        record = Manifest(paths['manifest'], telescope=telescope, night=night, job=job_id or '')
    fetcher = EsoFileFetcher(eso_token_fn(archive), workers=fcfg['workers'], attempts=fcfg['attempts'],
                             backoff=fcfg['backoff_seconds'], breaker=fcfg['breaker'], log=log)
    staging = os.path.join(paths['staging'], 'job-{}-{}-{}'.format(job_id or os.getpid(), telescope, night))
    res = add_only_fetch(rows, dest or night_dir(cfg['basedir'], telescope, night), staging, fetcher,
                         unpack=unpack, is_cube=_is_cube, transform=lambda ids, d: astra_transform(ids, d, 1),
                         verify=_verify_placed, replace=replace, foreign=foreign, record=record, log=log,
                         quarantine=os.path.join(paths['quarantine'], telescope, night), batch=int(fcfg['batch']))
    res['failed'] = dict(fetcher.failed)
    res['retry_later'] = sorted(fetcher.retry_later)
    res['bytes_downloaded'] = fetcher.bytes
    return res


def invalid_rows(path, rows, names):
    """Rows whose frames on disk fail verification (truncated or corrupt): what --replace invalid re-fetches."""
    files = row_files(rows, names)
    out = []
    for r in rows:
        for name in files.get(r['dp_id'], []):
            if not name.endswith('.fits'):
                continue  # .fits.fz and others are not swapped
            try:
                verify_frame(os.path.join(path, name))
            except (OSError, FitsError):
                out.append(r)
                break
    return out


def plan_fetch(rows, path, replace='never'):
    """(todo, summary) for a night: the rows to fetch given what is on disk (cube- and ACP-name-aware)."""
    names = list_frames(path)
    foreign = foreign_times(path)
    present, missing, _ = diff_eso_disk(rows, names, foreign)
    if replace == 'all':
        todo = list(rows)
    elif replace == 'invalid':
        todo = missing + invalid_rows(path, present, names)
    else:
        todo = missing
    summary = {'eso_rows': len(rows), 'eso_science': science_count(rows), 'disk_files': len(names),
               'foreign_files': len(foreign), 'missing': len(missing), 'missing_science': science_count(missing),
               'todo': len(todo), 'todo_science': science_count(todo),
               'todo_bytes': sum(int(float(r.get('access_estsize') or 0)) for r in todo) * 1024}
    return todo, summary, foreign


def download(cfg, telescope, night):
    age = (utcnow().date() - parse_night(night)).days
    recent = int(cfg['watcher']['recent_days'])
    if age > recent:
        print('refusing: {} is {} days old; SSO_download.py is only run on nights up to {} days old '
              '(use the look-back\'s add-only fetch for older nights)'.format(night, age, recent))
        return 2
    prog = cfg['telescopes'][telescope]['prog_id']
    path = night_dir(cfg['basedir'], telescope, night)
    archive = _archive(cfg)

    kind = transformation_check_kind(telescope, night)
    names = list_frames(path)
    doomed = []
    if kind:
        for n in names:
            try:
                header = read_header(os.path.join(path, n), max_blocks=24)
            except OSError:
                header = {}
            if transformation_check_would_delete(n, header, kind):
                doomed.append(n)
    if doomed:
        print('NOT running SSO_download.py: its transformation_check() would delete {} frames already on disk '
              '(e.g. {}); fetching missing frames add-only instead'.format(len(doomed), doomed[0]))
        rows = archive.frames(prog, night, night)
        _, missing, _ = diff_eso_disk(rows, names)
        fetch_rows(cfg, archive, telescope, night, missing, os.getenv('JOBQUEUE_JOB_ID'))
    else:
        log_path = os.path.join(cfg['basedir'], 'ESO_logs', '{}{}.log'.format(night, telescope))
        edate = (parse_night(night) + dt.timedelta(days=1)).strftime('%Y%m%d')
        cmd = [cfg['python'], 'download/SSO_download.py', '--dir', cfg['basedir'], '--telescope', telescope,
               '--sdate', night, '--edate', edate, '--max-retries', '1']
        print('running: {} > {}'.format(' '.join(cmd), log_path), flush=True)
        with open(log_path, 'w') as out:
            rc = subprocess.call(cmd, cwd=cfg['src_dir'], stdout=out, stderr=subprocess.STDOUT)
        print('SSO_download.py exit {}'.format(rc), flush=True)

    names = list_frames(path)
    try:
        rows = archive.frames(prog, night, night)
    except Exception as e:  # the check is informational; the frames on disk decide
        print('could not list the night at ESO to check the download ({}); {} files on disk'.format(e, len(names)))
        return 0 if names else EX_TEMPFAIL
    present, missing, _ = diff_eso_disk(rows, names)
    print('{} {}: ESO {} rows ({} science); on disk {} files; missing {} rows ({} science)'.format(
        telescope, night, len(rows), science_count(rows), len(names), len(missing), science_count(missing)))
    if rows and not names:
        print('nothing landed on disk; exit {} so the queue retries later'.format(EX_TEMPFAIL))
        return EX_TEMPFAIL
    if missing:
        print('{} rows still missing; processing what is here, the look-back will fetch the rest'.format(
            len(missing)))
    return 0


def fetch(cfg, store, telescope, night, rows_file=None, rerun_class='P1', dest=None, replace='never'):
    archive = _archive(cfg)
    if rows_file:
        rows = read_json(rows_file)
        if rows is None:
            print('cannot read {}'.format(rows_file))
            return 2
    else:
        rows = archive.frames(cfg['telescopes'][telescope]['prog_id'], night, night)
    path = dest or night_dir(cfg['basedir'], telescope, night)
    todo, summary, foreign = plan_fetch(rows, path, replace)
    print('{} {}: ESO {eso_rows} rows ({eso_science} science); on disk {disk_files} frames (+{foreign_files} under '
          'other names); missing {missing} rows ({missing_science} science); fetching {todo} rows, about '
          '{gb:.1f} GB{r}'.format(telescope, night, gb=summary['todo_bytes'] / 1e9,
                                  r='' if replace == 'never' else ' (replace: {})'.format(replace), **summary),
          flush=True)
    res = fetch_rows(cfg, archive, telescope, night, todo, dest=dest, replace=replace, foreign=foreign)
    res['plan'] = summary
    if dest:  # a trial into a scratch directory: record nothing, queue nothing
        print(res)
        return 0
    result_path = rows_file + '.result.json' if rows_file else os.path.join(
        queue_paths(cfg)['fetch'], '{}_{}_job{}.result.json'.format(telescope, night, os.getenv('JOBQUEUE_JOB_ID', 'x')))
    os.makedirs(os.path.dirname(result_path), exist_ok=True)
    write_json_atomic(result_path, res)
    if res.get('failed'):
        reasons = {}
        for why in res['failed'].values():
            reasons[why] = reasons.get(why, 0) + 1
        for why, n in sorted(reasons.items(), key=lambda x: -x[1])[:5]:
            print('{} rows not downloaded: {}'.format(n, why))
    if rerun_class == 'none':
        if res.get('retry_later'):
            print('{} rows may download later; exit {} so the queue retries'.format(len(res['retry_later']),
                                                                                   EX_TEMPFAIL))
            return EX_TEMPFAIL
        return 0

    lb = store.lookback_get(telescope, night) or {'telescope': telescope, 'night': night, 'fetches': 0,
                                                   'reruns': 0, 'data': {}}
    data = dict(lb.get('data') or {})
    data['settled'] = sorted(set(data.get('settled', [])) | set(res['settled']))
    data['last_fetch'] = {k: res[k] for k in ('requested', 'downloaded', 'added', 'added_science')}
    lb['data'] = data
    if res['added_science'] > 0:
        if int(lb.get('reruns') or 0) < int(cfg['lookback']['max_reruns_per_night']):
            objects = sorted({r.get('object', '').strip() for r in rows
                              if (r.get('dp_type') or '').upper() == 'OBJECT' and (r.get('object') or '').strip()})
            jid, created = jobspec.add(store, jobspec.pipeline_job(
                cfg, rerun_class, telescope, night, lock_targets=objects,
                frames=len(list_frames(night_dir(cfg['basedir'], telescope, night))), source='fetch',
                priority=float(night), note='{} science frames added by the look-back'.format(res['added_science'])))
            if created:
                lb['reruns'] = int(lb.get('reruns') or 0) + 1
            print('queued rerun job {} ({} science frames added)'.format(jid, res['added_science']))
        else:
            print('{} science frames added, but the automatic rerun limit is reached'.format(res['added_science']))
    store.lookback_put(lb)
    if todo and res['downloaded'] == 0:
        return EX_TEMPFAIL
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog='python -m jobqueue jobs')
    p.add_argument('--config')
    sub = p.add_subparsers(dest='job', required=True)
    d = sub.add_parser('download')
    d.add_argument('--telescope', required=True)
    d.add_argument('--night', required=True)
    f = sub.add_parser('fetch')
    f.add_argument('--telescope', required=True)
    f.add_argument('--night', required=True)
    f.add_argument('--rows', help='JSON list of ESO rows to fetch (default: everything ESO holds that is not on disk)')
    f.add_argument('--rerun-class', default='P1', choices=('P0', 'P1', 'P2', 'P3', 'none'),
                   help='queue a rerun of this class if science frames were added; none: download only')
    f.add_argument('--replace', default='never', choices=('never', 'invalid', 'all'),
                   help='invalid: also swap frames that fail verification; all: swap every frame that differs')
    f.add_argument('--dest', help='trial run: put the frames here instead of the night directory, record nothing')
    a = p.parse_args(argv)
    cfg = load_config(a.config)
    ensure_dirs(cfg)
    if a.job == 'download':
        return download(cfg, a.telescope, a.night)
    from .store import Store
    store = Store(queue_paths(cfg)['db'])
    return fetch(cfg, store, a.telescope, a.night, a.rows, a.rerun_class, dest=a.dest, replace=a.replace)


if __name__ == '__main__':
    sys.exit(main())
