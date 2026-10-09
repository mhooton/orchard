"""Job bodies that are not the pipeline itself: the nightly download and the look-back's add-only fetch.

download  SSO_download.py for one recent night, exactly as the nightly cron runs it (log to
          ESO_logs/<NIGHT><TEL>.log, its usual summary email), with --max-retries 1 because the watcher
          only queues a night once ESO has it. Refuses nights older than watcher.recent_days, and refuses
          (falling back to the add-only fetch) if transformation_check() would delete frames already on
          disk. Exits 0 whenever frames are on disk, so the pipeline job behind it runs even if a few are
          missing (the look-back fetches those); exits 75 (retry later) if nothing landed.
fetch     add-only fetch of the rows listed by the look-back; queues a rerun if science frames were added.
"""
import argparse
import datetime as dt
import io
import os
import subprocess
import sys
import contextlib

from . import jobspec
from .config import ensure_dirs, load_config, queue_paths
from .nights import (diff_eso_disk, list_frames, night_dir, read_header, science_count,
                     transformation_check_kind, transformation_check_would_delete)
from .util import parse_night, read_json, utcnow, write_json_atomic

EX_TEMPFAIL = 75


def _archive(cfg):
    from .eso import EsoArchive
    return EsoArchive(cfg['eso_env_file'])


def eso_download_fn(archive):
    """download(dp_ids, dest) through ESODownloader.download_files, its output scrubbed of credentials."""
    from .eso import scrub

    def download(dp_ids, dest):
        archive._ensure_token()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            got = archive._d.download_files([[d] for d in dp_ids], dest)
        print(scrub(buf.getvalue(), archive._secrets), end='', flush=True)
        return got
    return download


def fetch_rows(cfg, archive, telescope, night, rows, job_id=None, dest=None):
    from download.unpack_datacubes import is_datacube, unpack_datacube
    from download.astra_transform import astra_transform
    from .fetch import add_only_fetch

    def unpack(path, dest):
        ids, _ = unpack_datacube(path, dest)
        return [os.path.join(dest, i + '.fits') for i in ids]

    staging = os.path.join(queue_paths(cfg)['staging'], 'job-{}-{}'.format(job_id or os.getpid(), night))
    return add_only_fetch(rows, dest or night_dir(cfg['basedir'], telescope, night), staging, eso_download_fn(archive),
                          unpack=unpack, is_cube=is_datacube, transform=lambda ids, d: astra_transform(ids, d, 1))


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


def fetch(cfg, store, telescope, night, rows_file, rerun_class='P1', dest=None):
    rows = read_json(rows_file)
    if rows is None:
        print('cannot read {}'.format(rows_file))
        return 2
    archive = _archive(cfg)
    res = fetch_rows(cfg, archive, telescope, night, rows, os.getenv('JOBQUEUE_JOB_ID'), dest=dest)
    if dest:  # a trial into a scratch directory: record nothing, queue nothing
        print(res)
        return 0
    write_json_atomic(rows_file + '.result.json', res)

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
    if rows and res['downloaded'] == 0:
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
    f.add_argument('--rows', required=True)
    f.add_argument('--rerun-class', default='P1')
    f.add_argument('--dest', help='trial run: put the frames here instead of the night directory, record nothing')
    a = p.parse_args(argv)
    cfg = load_config(a.config)
    ensure_dirs(cfg)
    if a.job == 'download':
        return download(cfg, a.telescope, a.night)
    from .store import Store
    store = Store(queue_paths(cfg)['db'])
    return fetch(cfg, store, a.telescope, a.night, a.rows, a.rerun_class, dest=a.dest)


if __name__ == '__main__':
    sys.exit(main())
