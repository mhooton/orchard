#!/usr/bin/env python3
"""
promote.py
==========

The promotion step: publish pipeline output from PipelineOutput/v3 into PipelineOutput/v2, the tree the
SPECULOOS portal reads, only when the promotion check (promote_check.py) says PROMOTE.  Phase 1 of the
replacement for T12's sweep (perform_version_migration in main/ZLP_pipeline.sh), which moves any v3 target
directory holding a *4_diff.fits over the v2 one without looking at either.  The 2026-10-09 audit found that
the sweep nests on a night's second supersession and fails on the third, swaps the whole night's reduction
on a single-target rerun, strands good results when a later target crashes, and made 13 products worse.

Pipeline jobs run with --no_T12 and a 'promote' job after them (`q add ... --promote`), which runs this.

For each target-night, the candidate is a v3 target directory holding a *_diff.fits.  The step

  judges it  with promote_check.assess() (A1-A5 and R1-R3, over every star in the directory, since the
             directory is replaced as a whole), plus three checks of its own:
               S1  the run being promoted wrote every file in the directory: no file's inode last changed
                   before the run started.  Reruns reuse v3 directories, so an earlier run's files can sit
                   beside the new ones, or make up the whole directory.  It uses ctime, not mtime, because
                   T8 restores stacks from v2/StackImages with copy2, which keeps their old mtime
               S2  the run is known: the latest pipeline run of the night that covers the target, from the
                   queue's jobs or the logs in v3/ and v2/<TEL>/logs
               A3  with --identity table: the plan file and the K1 target table agree on the star
  PROMOTE    1. renames the v2 incumbent, if any, to PipelineOutput/superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>/,
                a directory made for it, so a night can be superseded any number of times without nesting
             2. renames the candidate into v2
             3. touches its *_5_diff.fits, so that the portal's 90-day mtime window ingests them
             4. installs a copy of the run's reduction when v2 has none for the night, or one without
                masters.  It replaces v2's only when the run covered every target of the night, and the
                reduction has masters and 1_image_*.list files and was made by this run.  v2 is never left
                without a reduction
             5. copies the run's log and PDF into v2/<TEL>/logs and v2/<TEL>/reports under names that carry the
                run (<NIGHT>_1_v3.run<start>[-j<job>].log, <TEL>_<NIGHT>.run<...>.pdf), and links the plain
                name to it when v2 has none, since completeness_check.py looks for <TEL>_<NIGHT>.pdf by name
             6. drops the superseded copy's frames and per-frame files: everything directly inside its run
                directories (1/ and any other all-digit one), apart from the frame lists, and a superseded
                reduction's flats/.  *_diff.fits, *_output.fits, lightcurves/, plots and stacks stay, with
                thinned.tsv listing what went.  The frames go last, so that v2 lacks the target only for the
                moment between two renames, not for the time the deletions take
  HOLD and   With --apply the verdict and its reasons go into the journal (promote_journal.py), and the approval
  REJECT     policy (promote_policy.py) decides it: promote (the PROMOTE steps, with the files the run did not
             write set aside first), reject (v2 keeps what it has), reprocess (the pipeline, not the data, is at
             fault: listed for a rerun), or a person.  `accept <id> --reason` performs the PROMOTE steps by
             hand, and `reject <id> --reason` records a person's decision.  --no-policy leaves every one for a
             person.

A candidate restaged for a rerun from a late stage (utils/restage.py, docs/rerun-from-stage.md) carries a record of
its inputs. S1 then judges only what the rerun wrote (its manifest), and, when the inputs came from v2, the PROMOTE
steps are an overlay instead: the manifest's files replace or join v2's, the v2 files they replace are set aside
under superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-overlay/, and v2's directory, reduction and every other file stay as
they are.  O1-O3 hold a restaged candidate whose inputs, or v2, changed after the restage.

Each promotion is journalled before and after each step.  One cut short by a crash is listed by `list`
and in the queue's daily report, and finished by `resume`; `run --apply` finishes those of a night before
judging the night again.  Everything is a dry run unless --apply is given.  A night is left alone while a
queue job, or a pipeline or promotion run outside the queue, works on it.  Every process that changes the
tree carries --apply --telescope T --night N on its command line, so that the queue's dispatcher sees the
night as locked (jobqueue.resources.parse_promote_argv); commands that span nights run one child per night.

Usage (inside orchard-server; `q promote ...` is the same)
    python utils/promote.py run [--telescope T] [--night N] [--target X] [--job ID] [--apply]
    python utils/promote.py list [--all] [--grep TEXT] | show ID
    python utils/promote.py accept ID[-ID][,ID...] --reason TEXT [--apply]
    python utils/promote.py reject ID[-ID][,ID...] --reason TEXT
    python utils/promote.py resume [--promotion ID] [--apply]
    python utils/promote.py abandon PROMOTION_ID --reason TEXT
"""

import argparse
import csv
import datetime as dt
import getpass
import glob
import hashlib
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import stat
import subprocess
import sys
import time
from collections import Counter, OrderedDict

import numpy as np

from utils import promote_check
from utils import promote_policy
from utils import restage as restage_mod
from utils import stages as stage_map
from utils import target_management
from utils.promote_check import Identity, RawFrames, assess
from utils.promote_journal import Journal, iso, parse_iso, utcnow
from utils.target_list import load_target_list

RUNNAME = '1'                       # the queue's pipeline_runname; names the logs and reduction files
SUPERSEDED = 'superseded'           # PipelineOutput/superseded/<TEL>/<NIGHT>/<TARGET or reduction>/<stamp>/
KEEP_IN_RUN_DIRS = ('.dat', '.dat_phot')   # frame lists: tiny, and they say which frames a product used
TOLERANCE = 1.0                     # seconds: run starts are logged to the second
EX_TEMPFAIL = 75                    # the queue retries a job with this exit code an hour later
PROC = '/proc'                      # where night_busy looks for pipeline runs outside the queue (tests point it away)
STEPS = ('check', 'aside', 'supersede', 'publish', 'touch', 'reduction', 'records', 'thin', 'unstage')
# a night restaged from v2 for a rerun from a late stage (utils/restage.py): only what the rerun wrote goes into v2
OVERLAY_STEPS = ('check', 'install', 'touch', 'records', 'leftover', 'thin')
OVERLAY = '-overlay'                # superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-overlay/: v2 files an overlay replaced
RESTAGED = '-restage'               # superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-restage/: the restaged inputs, after
STALE = '-stale'                    # superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-stale/: files set aside before publishing
RANK = {'PROMOTE': 0, 'HOLD': 1, 'REJECT': 2}
EPOCH = dt.datetime(1970, 1, 1)

DATE_LINE = re.compile(r'^(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) (\w{3}) +(\d{1,2}) (\d\d:\d\d:\d\d) (\w+) (\d{4})$')
LOG_NAME = re.compile(r'^\d{8}_[^_]+_v3(?:_(.+))?\.log$')
COPY_STAMP = re.compile(r'\.run\d{8}T\d{6}Z(?:-j\d+)?(?:\.\d+)?$')


class PromoteError(Exception):
    """A promotion cannot go on safely: something is not as it was when the promotion was planned."""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def norm_target(name):
    return name.strip().replace('--', ' ').upper()


def epoch(t):
    return (t - EPOCH).total_seconds()


def lino(path):
    """(st_dev, st_ino) of path itself, or None: a renamed directory keeps them, which is how a resumed step
    tells whether its rename already happened."""
    try:
        st = os.lstat(path)
    except OSError:
        return None
    return [st.st_dev, st.st_ino]


def nearest_existing(path):
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:
            break
        path = parent
    return path


def same_device(*paths):
    """Refuse unless every path (or its nearest existing parent) is on one filesystem: os.rename never copies,
    and nothing here may fall back to copy-then-delete."""
    devs = {os.stat(nearest_existing(p)).st_dev for p in paths if p}
    if len(devs) > 1:
        raise PromoteError('{} are not on one filesystem'.format(', '.join(p for p in paths if p)))


def file_sha1(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def same_content(a, b):
    try:
        return os.path.getsize(a) == os.path.getsize(b) and file_sha1(a) == file_sha1(b)
    except OSError:
        return False


def copy_no_overwrite(src, dst):
    """Copy src to dst without ever replacing a file. An identical file already at dst counts as the copy; a
    different one sends the copy to dst with .2, .3, ... before the extension. Returns (path, copied)."""
    stem, ext = os.path.splitext(dst)
    for k in range(1, 100):
        path = dst if k == 1 else '{}.{}{}'.format(stem, k, ext)
        if os.path.lexists(path):
            if same_content(src, path):
                return path, False
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = os.path.join(os.path.dirname(path), '.{}.promote-{}'.format(os.path.basename(path), os.getpid()))
        shutil.copy2(src, tmp)
        try:
            os.link(tmp, path)          # fails if the name was taken meanwhile: never overwrites
        except FileExistsError:
            os.unlink(tmp)
            continue
        os.unlink(tmp)
        return path, True
    raise PromoteError('no free name for a copy of {} beside {}'.format(src, dst))


def jsonable(row):
    """An assess() row as plain JSON values (NaN -> None, numpy scalars -> Python)."""
    out = {}
    for k, v in row.items():
        if isinstance(v, (np.generic,)):
            v = v.item()
        if isinstance(v, float) and not math.isfinite(v):
            v = None
        out[k] = v
    return out


def run_json(run):
    if run is None:
        return None
    return {k: (iso(v) if isinstance(v, dt.datetime) else v) for k, v in run.items()}


def run_from_json(run):
    if not run:
        return None
    out = dict(run)
    for k in ('start', 'end'):
        out[k] = parse_iso(run.get(k)) if run.get(k) else None
    return out


def describe_run(run):
    if run is None:
        return 'no run'
    who = 'queue job {}'.format(run['job']) if run.get('job') else os.path.basename(run.get('log') or '?')
    return '{}, started {}'.format(who, iso(run['start']))


def run_stamp(run):
    """'run<start>[-j<job>]': names the copies of a run's log and PDF in v2."""
    job = '-j{}'.format(run['job']) if run.get('job') else ''
    return 'run{}{}'.format(run['start'].strftime('%Y%m%dT%H%M%SZ'), job)


# ---------------------------------------------------------------------------
# The tree
# ---------------------------------------------------------------------------

class Tree(object):
    def __init__(self, base):
        self.base = base
        self.po = os.path.join(base, 'PipelineOutput')

    def output(self, version, tel, night, name=None):
        path = os.path.join(self.po, version, tel, 'output', night)
        return os.path.join(path, name) if name else path

    def logs(self, version, tel):
        return os.path.join(self.po, version, tel, 'logs')

    def reports(self, version, tel):
        return os.path.join(self.po, version, tel, 'reports')

    def superseded(self, tel, night, name, stamp):
        return os.path.join(self.po, SUPERSEDED, tel, night, name, stamp)

    def under_superseded(self, path):
        root = os.path.join(os.path.realpath(self.po), SUPERSEDED) + os.sep
        return os.path.realpath(path).startswith(root)

    def candidates(self, tels=None, nights=None):
        """{(tel, night): [target, ...]} for every v3 target directory holding a *_diff.fits."""
        out = OrderedDict()
        for tel, night, target, _, _, _ in promote_check.v3_target_nights(self.base, tels, nights):
            out.setdefault((tel, night), []).append(target)
        return out


# ---------------------------------------------------------------------------
# Runs: which pipeline run a candidate should come from
# ---------------------------------------------------------------------------

def _date_line(line):
    m = DATE_LINE.match(line.strip())
    if not m or m.group(4) not in ('UTC', 'GMT'):
        return None
    try:
        return dt.datetime.strptime(' '.join(m.group(1, 2, 3, 5)), '%b %d %H:%M:%S %Y')
    except ValueError:
        return None


def parse_log(path):
    """What a ZLP_pipeline.sh log says about its run: start and end (its first and last `date` lines, UTC in the
    container), the targets it ran and whether they were given, and whether it completed. None without a start."""
    run = {'source': 'log', 'log': path, 'job': None, 'state': None, 'start': None, 'end': None, 'targets': None,
           'specified': None, 'complete': False, 'science': True, 'stages': None}
    try:
        with open(path, 'rb') as f:
            for k, raw in enumerate(f):
                if k > 5000:
                    break
                line = raw.decode('utf-8', 'replace').rstrip('\r\n')
                if run['start'] is None:
                    run['start'] = _date_line(line)
                if run['stages'] is None:
                    run['stages'] = stage_map.parse_stages_line(line)       # 'Stages run: T10 T11' (since 2026-10)
                for label, given in (('Using specified targets:', True), ('Using discovered targets:', False)):
                    if line.startswith(label):
                        run['targets'], run['specified'] = line[len(label):].split(), given
                if run['start'] is not None and run['targets'] is not None:
                    break
            f.seek(0, os.SEEK_END)
            f.seek(max(0, f.tell() - 65536))
            tail = f.read().decode('utf-8', 'replace').splitlines()
    except OSError:
        return None
    if run['start'] is None:
        return None
    ends = [t for t in (_date_line(line) for line in tail) if t]
    run['end'] = ends[-1] if ends and ends[-1] > run['start'] else None
    run['complete'] = any(line.strip() == 'PIPELINE COMPLETE' for line in tail)
    if run['stages'] is not None:
        run['science'] = 10 in run['stages']
    if run['targets'] is None:
        # a log with neither line: its name says whether targets were given (<NIGHT>_<RUN>_v3_<T1_T2>.log)
        m = LOG_NAME.match(os.path.basename(path))
        run['specified'] = bool(m and m.group(1))
        if run['specified']:
            run['targets'] = m.group(1).split('_')
    return run


def runs_from_logs(tree, tel, night):
    """Runs recorded by logs in v3/ and v2/<TEL>/logs (T12 moves a night's log to v2 when it sweeps any target
    of it), the v3 copy first; copies of one log count once."""
    out, seen = [], set()
    for version in ('v3', 'v2'):
        for path in sorted(glob.glob(os.path.join(tree.logs(version, tel), '{}_*_v3*.log'.format(night)))):
            run = parse_log(path)
            if run is None:
                continue
            key = (run['start'], tuple(sorted(norm_target(t) for t in run['targets'] or ())), run['specified'])
            if key not in seen:
                seen.add(key)
                out.append(run)
    return out


def open_queue(path):
    """The queue database, read-only, or None."""
    if not path or not os.path.exists(path):
        return None
    conn = sqlite3.connect('file:{}?mode=ro'.format(path), uri=True, timeout=60)
    conn.row_factory = sqlite3.Row
    return conn


def queue_job(conn, job_id):
    if conn is None:
        return None
    r = conn.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone()
    return dict(r) if r else None


def science_argv(argv):
    """Whether a ZLP_pipeline.sh invocation runs T10, the stage that writes the light curves. An --only_T11 run (the
    PDF and frame packing) or an --only_T12 sweep makes no products, so S1 must not judge files against it."""
    return 10 in stage_map.stages_from_argv(argv)


def runs_from_queue(conn, tel, night):
    """Pipeline jobs of the night that have started, as runs (the latest attempt of each)."""
    if conn is None:
        return []
    out = []
    rows = conn.execute("SELECT * FROM jobs WHERE kind = 'pipeline' AND telescope = ? AND night = ? "
                        'AND started_at IS NOT NULL ORDER BY id', (tel, night)).fetchall()
    for r in rows:
        targets = (r['targets'] or '').split()
        try:
            argv = json.loads(r['argv'] or '[]')
        except ValueError:
            argv = []
        out.append({'source': 'job', 'job': r['id'], 'state': r['state'], 'start': parse_iso(r['started_at']),
                    'end': None if r['state'] == 'running' else parse_iso(r['finished_at']),
                    'targets': targets or None, 'specified': bool(targets), 'complete': None, 'log': None,
                    'job_log': r['log_path'], 'exit_code': r['exit_code'], 'code': r['code_sha'],
                    'science': science_argv(argv), 'stages': stage_map.stages_from_argv(argv)})
    return out


def discover_runs(tree, conn, tel, night):
    """Every run of the night that left a trace: queue jobs, each with its own log while that log is still there
    (a later run of the night overwrites it), and the runs outside the queue that only a log records."""
    logs = runs_from_logs(tree, tel, night)
    jobs = runs_from_queue(conn, tel, night)
    used = set()
    for job in jobs:
        for k, run in enumerate(logs):
            if k in used or bool(run['specified']) != job['specified']:
                continue
            if job['specified'] and {norm_target(t) for t in run['targets'] or ()} != \
                    {norm_target(t) for t in job['targets']}:
                continue
            if job['start'] - dt.timedelta(seconds=120) <= run['start'] <= (job['end'] or utcnow()):
                used.add(k)
                job['log'] = run['log']
                job['complete'] = run['complete']
                if not job['specified']:
                    job['targets'] = run['targets']
                break
    return jobs + [run for k, run in enumerate(logs) if k not in used]


def covers(run, target):
    """A run given targets covers those; one without covers the targets it discovered, or the whole night if its
    log is gone (a queue job whose log a later run overwrote)."""
    return run['targets'] is None or norm_target(target) in {norm_target(t) for t in run['targets']}


def latest_run(runs, target):
    cover = [r for r in runs if covers(r, target)]
    return max(cover, key=lambda r: (r['start'], r['job'] or 0)) if cover else None


def science_run(runs, target):
    """The run a candidate's light curves come from: the latest run covering the target that ran T10. The runs after
    it that made no light curves (an --only_T11 PDF run, an --only_T12 sweep) are listed in its 'later', oldest first,
    because their PDF and log belong to the night's records. None if no such run is known."""
    order = lambda r: (r['start'], r['job'] or 0)
    cover = sorted((r for r in runs if covers(r, target)), key=order)
    made = [r for r in cover if r.get('science', True)]
    if not made:
        return None
    run = dict(made[-1])
    run['later'] = [run_json(r) for r in cover if order(r) > order(run)]
    return run


def night_busy(conn, tel, night, own_job=None, proc=None):
    """Why the night must be left alone now ('' if it need not): a queue job that holds its lock is running, or a
    pipeline, download or promotion process outside the queue works on it."""
    lock = 'night:{}:{}'.format(tel, night)
    why, sessions = [], [os.getsid(0)]
    if conn is not None:
        for r in conn.execute("SELECT id, kind, locks, pid FROM jobs WHERE state = 'running'"):
            sessions.append(r['pid'])           # each queue job runs in its own session, led by its runner
            if r['id'] != own_job and lock in json.loads(r['locks'] or '[]'):
                why.append('queue job {} ({}) is running on it'.format(r['id'], r['kind']))
    from jobqueue.resources import external_locks
    held = external_locks(sessions, proc or PROC)
    if lock in held:
        why.append('a process outside the queue holds it: {}'.format(held[lock]))
    return '; '.join(why)


# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------

def scan_dir(path, since=None):
    """A directory tree's fingerprint (relative path, size and mtime of every file), its size and identity, and
    the files whose inode last changed before `since` (UTC datetime): files the run did not write."""
    out = {'fp': None, 'files': 0, 'bytes': 0, 'stale': [], 'ino': lino(path)}
    if not path or not os.path.isdir(path) or os.path.islink(path):
        return out
    limit = None if since is None else epoch(since) - TOLERANCE
    h = hashlib.sha1()
    for root, dirs, names in os.walk(path):
        dirs.sort()
        for name in sorted(names):
            p = os.path.join(root, name)
            try:
                st = os.lstat(p)
            except OSError:
                continue
            rel = os.path.relpath(p, path)
            h.update('{}\0{}\0{}\n'.format(rel, st.st_size, st.st_mtime_ns).encode('utf-8', 'surrogateescape'))
            out['files'] += 1
            out['bytes'] += st.st_size
            if limit is not None and st.st_ctime < limit:
                out['stale'].append(rel)
    out['fp'] = h.hexdigest()
    return out


def pwv_state(path):
    """Whether a *_diff.fits carries T10's PWV correction: 'applied' (a PWV table with data and a finite
    PWV_LIGHTCURVE), 'missing' (T10 wrote an empty table or an all-NaN curve, as when the LHATPRO fetch fails), or
    None (no file, or unreadable). Artemis, SAINT-EX and targets without a usable Teff never get one."""
    if not path or not os.path.isfile(path):
        return None
    try:
        from astropy.io import fits
        with fits.open(path, memmap=False) as h:
            names = [x.name for x in h]
            rows = len(h['PWV'].data) if 'PWV' in names and h['PWV'].data is not None else 0
            curves = [x for x in names if x.startswith('PWV_LIGHTCURVE_')]
            finite = bool(curves) and h[curves[0]].data is not None and \
                bool(np.isfinite(np.asarray(h[curves[0]].data, dtype=float)).any())
    except Exception:
        return None
    return 'applied' if rows and finite else 'missing'


def stale_summary(stale, examples=3):
    """'12 in 1/, 2 in lightcurves/, 1 *_diff.fits (e.g. ...)'."""
    def kind(rel):
        parts = rel.split(os.sep, 1)
        if len(parts) == 2:
            return 'in {}/'.format(parts[0])
        for suffix in ('_diff.fits', '_output.fits'):
            if rel.endswith(suffix):
                return '*' + suffix
        return 'other top-level'
    counts = Counter(kind(r) for r in stale)
    text = ', '.join('{} {}'.format(n, k) for k, n in counts.most_common())
    return '{} (e.g. {})'.format(text, ', '.join(sorted(stale)[:examples]))


def reduction_info(path):
    """A night's reduction directory: does it exist, its masters and image lists, its top-level files."""
    info = {'exists': os.path.isdir(path) and not os.path.islink(path), 'masters': [], 'lists': [], 'files': {}}
    if not info['exists']:
        return info
    for name in sorted(os.listdir(path)):
        p = os.path.join(path, name)
        try:
            st = os.lstat(p)
        except OSError:
            continue
        if not stat.S_ISREG(st.st_mode):
            continue
        info['files'][name] = [st.st_size, st.st_mtime_ns, st.st_ctime]
        if name.startswith(RUNNAME + '_Master') and name.endswith('.fits'):
            info['masters'].append(name)
        elif name.startswith(RUNNAME + '_image_') and name.endswith('.list'):
            info['lists'].append(name)
    return info


def list_targets(info):
    """The night's targets according to its reduction: one 1_image_<target>.list each (T1 writes them)."""
    n = len(RUNNAME + '_image_')
    return sorted({norm_target(name[n:-len('.list')]) for name in info['lists']})


def thin(path, reduction=False):
    """Delete a superseded copy's frames and per-frame files: every file directly inside its run directories (1/
    and any other all-digit directory) except the frame lists, or, for a reduction, inside flats/. Never
    descends further. Each file is listed in <path>/thinned.tsv before it goes. Returns the counts."""
    if reduction:
        dirs = [os.path.join(path, 'flats')]
    else:
        dirs = [os.path.join(path, d) for d in sorted(os.listdir(path)) if d.isdigit()]
    manifest = os.path.join(path, 'thinned.tsv')
    out = {'deleted': 0, 'bytes': 0, 'kept_lists': 0, 'manifest': manifest}
    with open(manifest, 'a') as m:
        for d in dirs:
            if not os.path.isdir(d) or os.path.islink(d):
                continue
            for name in sorted(os.listdir(d)):
                p = os.path.join(d, name)
                try:
                    st = os.lstat(p)
                except FileNotFoundError:
                    continue
                if stat.S_ISDIR(st.st_mode):
                    continue
                if not reduction and name.endswith(KEEP_IN_RUN_DIRS):
                    out['kept_lists'] += 1
                    continue
                m.write('{}\t{}\t{}\n'.format(os.path.relpath(p, path), st.st_size, int(st.st_mtime)))
                m.flush()
                os.unlink(p)
                out['deleted'] += 1
                out['bytes'] += st.st_size
    return out


# ---------------------------------------------------------------------------
# promote_check's helpers with what the step needs on top
# ---------------------------------------------------------------------------

class CachedRawFrames(RawFrames):
    """promote_check.RawFrames with its per-night counts kept in SQLite, which concurrent promote jobs can share,
    and invalidated when the night's raw directory, or one of its subdirectories, changed after the count was
    made (the ESO backfill added frames to 845 nights after the 2026-10-09 cache). seed: a cache file in
    RawFrames' own JSON format, whose entries count as made when the file was last written."""

    VERSION = 1

    def __init__(self, base, path=None, seed=None, readonly=False, threads=8):
        super(CachedRawFrames, self).__init__(base, cache_path=None, threads=threads)
        self.readonly = readonly
        self.stats = Counter()
        self.db = None
        if path and (not readonly or os.path.exists(path)):
            if readonly:
                self.db = sqlite3.connect('file:{}?mode=ro'.format(path), uri=True, timeout=60)
            else:
                os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
                self.db = sqlite3.connect(path, timeout=60)
                self.db.execute('PRAGMA journal_mode=WAL')
                self.db.execute('CREATE TABLE IF NOT EXISTS raw_counts (telescope TEXT, night TEXT, counted_at REAL, '
                                'version INTEGER, counts TEXT, PRIMARY KEY (telescope, night))')
                self.db.commit()
        self.seed, self.seed_at = {}, None
        if seed and os.path.exists(seed):
            with open(seed) as f:
                self.seed = json.load(f)
            self.seed_at = os.path.getmtime(seed)

    def changed_at(self, tel, date):
        """Latest mtime of the night's raw directory and its subdirectories (RawFrames reads both levels), or
        None if there is no raw directory."""
        imgdir = os.path.join(self.base, 'Observations', tel, 'images', date)
        try:
            latest = os.stat(imgdir).st_mtime
            for entry in os.scandir(imgdir):
                if entry.is_dir(follow_symlinks=False):
                    latest = max(latest, entry.stat(follow_symlinks=False).st_mtime)
        except OSError:
            return None
        return latest

    @staticmethod
    def _valid(counts, counted_at, changed):
        if changed is None:
            return counts is None
        return counts is not None and counted_at is not None and changed < counted_at

    def night(self, tel, date):
        key = '%s/%s' % (tel, date)
        if key in self.cache:
            return self.cache[key]
        changed = self.changed_at(tel, date)
        if self.db is not None:
            try:
                row = self.db.execute('SELECT counted_at, version, counts FROM raw_counts WHERE telescope = ? AND '
                                      'night = ?', (tel, date)).fetchone()
            except sqlite3.OperationalError:     # a read-only cache file without the table yet
                row = None
            if row and row[1] == self.VERSION and self._valid(json.loads(row[2]), row[0], changed):
                self.cache[key] = json.loads(row[2])
                self.stats['cached'] += 1
                return self.cache[key]
        if key in self.seed and self._valid(self.seed[key], self.seed_at, changed):
            self.cache[key] = self.seed[key]
            self.stats['seed'] += 1
            return self.cache[key]
        if key in self.seed:
            self.stats['seed_invalidated'] += 1
        counted_at = time.time()             # before listing: a change while counting invalidates the count
        counts = super(CachedRawFrames, self).night(tel, date)
        self.stats['counted'] += 1
        if self.db is not None and not self.readonly:
            self.db.execute('INSERT OR REPLACE INTO raw_counts VALUES (?, ?, ?, ?, ?)',
                            (tel, date, counted_at, self.VERSION, json.dumps(counts)))
            self.db.commit()
        return counts


class TableIdentity(Identity):
    """A3 from the K1 target table (ml_40pc_v2.csv) for Sp targets, with both of its Gaia IDs, and from the plan
    file for every other target, as before. A name can have several rows (the components of a resolved binary
    share one): the plan's row when the plan names one of them, else the one the light curves carry. Where the
    plan names a star that is in none of the name's rows, the disagreement is kept in .conflicts, and promote.py
    holds the candidate."""

    def __init__(self, base, table_path):
        super(TableIdentity, self).__init__(base, target_list=table_path)
        self.table = load_target_list(table_path)
        self.rows = {}
        if self.table is not None and self.table.has('Sp_ID'):
            for i in range(len(self.table)):
                if self.table.ids(i):
                    self.rows.setdefault(norm_target(self.table.sp_id(i)), []).append(i)
        self.conflicts = {}

    def expected(self, tel, date, target, catalogue_ids=()):
        gid, src, dec = super(TableIdentity, self).expected(tel, date, target, catalogue_ids)
        if self.table is None or not norm_target(target).startswith('SP'):
            return gid, src, dec
        rows = self.rows.get(norm_target(target), [])
        if not rows:
            return gid, src + ', not in the target table', dec
        plan_id = gid if src.startswith('plan') else ''
        key = (tel, date, target)
        if plan_id:
            same = [i for i in rows if plan_id in self.table.ids(i)]
            if not same:
                other = super(TableIdentity, self).counterpart(plan_id, dec)     # the plan may use the other release
                same = [i for i in rows if other & set(self.table.ids(i))]
            if not same:
                self.conflicts[key] = 'the plan file names {} ({}), the target table {} for {}'.format(
                    plan_id, src, ' or '.join(' / '.join(self.table.ids(i)) for i in rows), target)
            row = (same or rows)[0]
        else:
            hit = [i for i in rows if set(self.table.ids(i)) & set(catalogue_ids or ())]
            row = (hit or rows)[0]
        ids = self.table.ids(row)
        _, table_dec = self.table.coords(row)
        how = 'target table' + ('' if not plan_id else ', plan disagrees' if key in self.conflicts else ', plan agrees')
        return (plan_id if plan_id in ids else ids[0]), how, table_dec

    def counterpart(self, gid, dec):
        row = self.table.row_for_id(gid) if self.table is not None else None
        if row is not None:
            return set(self.table.ids(row)) - {gid}
        return super(TableIdentity, self).counterpart(gid, dec)


# ---------------------------------------------------------------------------
# The promotion step
# ---------------------------------------------------------------------------

class Promoter(object):
    def __init__(self, base, journal=None, queue=None, identity=None, raw=None, th=None, apply=False, own_job=None,
                 invocation=None, proc=None, now_fn=utcnow, out=print, policy=False, policy_th=None):
        self.tree = Tree(base)
        self.journal = journal
        self.queue = queue
        self.identity = identity or Identity(base)
        self.raw = raw
        self.th = dict(th or promote_check.THRESHOLDS)
        self.policy = policy                 # decide HOLDs and REJECTs with promote_policy.py
        self.policy_th = dict(policy_th or {})
        self.apply = apply
        self.own_job = own_job
        self.invocation = invocation
        self.proc = proc
        self.now = now_fn
        self.out = out
        self._dry_reductions = {}

    # ------------------------------------------------------------------ judging

    def restage_record(self, tel, night, target):
        """The restage record of a candidate (utils/restage.py), if it was restaged for a rerun from a late stage."""
        rec = restage_mod.load_record(self.tree.base, tel, night, target)
        if rec is None or os.path.abspath(rec['dest']) != os.path.abspath(self.tree.output('v3', tel, night, target)):
            return None
        return rec

    @staticmethod
    def run_stale(stale, run, rec, night, target):
        """The files of a candidate that predate its run and should not: S1 limited to what the run had to write.
        A restaged candidate's unchanged inputs (its record says which) are exempt; so, after a run that started at a
        later stage (--from_T10, or --no_T1 ... --no_T9), are the products of the stages it did not run (the calibrated
        frames, the stacks, the night's *_output.fits), which it read and must predate it. Returns (stale, the
        record's check or None)."""
        if rec is not None:
            st = restage_mod.check(rec)
            inputs = set(st['unchanged'])
            return [r for r in stale if r not in inputs], st
        first = stage_map.first_stage(run.get('stages')) if run else None
        if first and first > 1:
            return [r for r in stale if not stage_map.made_before(r, first, night=night, target=target)], None
        return list(stale), None

    def judge(self, tel, night, target, runs):
        cand = self.tree.output('v3', tel, night, target)
        inc = self.tree.output('v2', tel, night, target)
        run = science_run(runs, target)
        rec = self.restage_record(tel, night, target)
        c = scan_dir(cand, run['start'] if run else None)
        c['stale'], st = self.run_stale(c['stale'], run, rec, night, target)
        i = scan_dir(inc)
        row = assess(tel, night, target, cand, inc if i['fp'] is not None else '', self.identity, self.raw, self.th)
        hold = []
        overlay = rec is not None and rec['source'] == 'v2'
        if rec is not None:
            bad = st['changed'] + st['missing']
            if rec.get('state') != 'staged':
                hold.append('O1 the restage of {} did not finish ({} is still {})'.format(
                    target, os.path.basename(rec['_path']), rec.get('state')))
            if bad:
                linked = {e['rel'] for e in rec['files'] if e['how'] == 'link'} & set(bad)
                hold.append('O1 {} of the {} restaged inputs changed or went after the restage ({}){}'.format(
                    len(bad), len(rec['files']), ', '.join(sorted(bad)[:3]),
                    "; {} were hard links to v2's files, so v2's changed too".format(len(linked)) if linked else ''))
            if st['reduction_changed']:
                hold.append("O1 the night's v3 reduction changed after the restage ({})".format(
                    ', '.join(st['reduction_changed'][:3])))
            if overlay and st['src_changed']:
                hold.append('O2 v2 changed after the restage: {} of the files the rerun read are no longer what v2 '
                            'holds ({})'.format(len(st['src_changed']), ', '.join(sorted(st['src_changed'])[:3])))
            if overlay and run and run.get('code'):
                # files the rerun made (a restaged frame it edited keeps the ORCHREL of the run that made it)
                m = restage_mod.manifest(rec, sha=False)
                restaged = {e['rel'] for e in rec['files']}
                foreign = sorted(e['rel'] for e in m['written'] if e['rel'] not in restaged and e.get('orchrel') and
                                 e['orchrel'] != run['code'])
                if foreign:
                    hold.append('O3 {} file(s) of the rerun were made by another code release than its run ({}): '
                                '{}'.format(len(foreign), run['code'], ', '.join(foreign[:3])))
        if run is None and latest_run(runs, target) is not None:
            hold.append('S2 no run of the night covering {} made light curves (the latest, {}, did not run T10), so '
                        'its files cannot be checked against the run that wrote them'.format(
                            target, describe_run(latest_run(runs, target))))
        elif run is None:
            hold.append('S2 no pipeline run of the night covering {} was found (no queue job, no log in v3/ or v2/{}/'
                        'logs), so its files cannot be checked against the run'.format(target, tel))
        elif c['stale'] and len(c['stale']) == c['files']:
            hold.append('S1 the run ({}) wrote none of its {} files: an earlier run made this directory'.format(
                describe_run(run), c['files']))
        elif c['stale']:
            hold.append('S1 {} of its {} files predate the run ({}): {}'.format(
                len(c['stale']), c['files'], describe_run(run), stale_summary(c['stale'])))
        conflict = getattr(self.identity, 'conflicts', {}).get((tel, night, target))
        if conflict:
            hold.append('A3 ' + conflict)
        # P1: a replacement must not lose the PWV correction v2 has (T10's fetch fails when appct's containers
        # cannot resolve names, as from 2026-10-18); a gap fill without one is still better than nothing
        notes = [row['notes']] if row['notes'] else []
        c_pwv = pwv_state(os.path.join(cand, row['c_file'])) if row.get('c_file') else None
        i_pwv = pwv_state(os.path.join(inc, row['i_file'])) if row.get('i_file') and i['fp'] is not None else None
        if c_pwv == 'missing' and i_pwv == 'applied':
            hold.append('P1 the candidate has no PWV correction ({} has an empty PWV table or curve) and v2\'s '
                        'curve has one'.format(row['c_file']))
        elif c_pwv == 'missing':
            notes.append('no PWV correction in the candidate (nor in v2)' if i_pwv else
                         'no PWV correction in the candidate')
        if rec is not None:
            notes.append('restaged from {} for a rerun from T{} ({} inputs){}'.format(
                rec['source'], rec['from_stage'], len(rec['files']),
                ': promoted as an overlay of what the rerun wrote' if overlay else ''))
        verdict = row['verdict']
        if hold and verdict == 'PROMOTE':
            verdict = 'HOLD'
        return {'telescope': tel, 'night': night, 'target': target, 'candidate': cand,
                'incumbent': inc if i['fp'] is not None else '', 'verdict': verdict,
                'checker_verdict': row['verdict'], 'reasons': '; '.join(x for x in [row['reasons']] + hold if x),
                'notes': '; '.join(notes), 'metrics': jsonable(row), 'run': run, 'cand_fp': c['fp'],
                'inc_fp': i['fp'], 'cand_files': c['files'], 'stale_files': len(c['stale']), 'stale': c['stale'],
                'mode': 'overlay' if overlay else 'whole'}

    # ------------------------------------------------------------------ the journal

    def record(self, j):
        """Journal a judged candidate unless the journal holds the same verdict on the same files. Returns
        (verdict row, None) to act on it, or (row, why) to leave it."""
        prev = self.journal.latest_verdict(j['telescope'], j['night'], j['target'])
        same = bool(prev) and prev['candidate'] == j['candidate'] and prev['cand_fp'] == j['cand_fp'] and \
            prev['inc_fp'] == j['inc_fp']
        if prev and prev['state'] == 'promoting' and prev['promotion_id'] is None:
            # a crash between journalling a PROMOTE verdict and starting its promotion: nothing was moved
            self.journal.set_verdict_state(prev['id'], 'outdated', expect=('promoting',))
            prev = dict(prev, state='outdated')
        if prev and prev['state'] == 'promoting':
            return prev, 'promotion {} of verdict {} has not finished'.format(prev['promotion_id'], prev['id'])
        if same and prev['state'] == 'rejected':
            return prev, 'rejected by {} on {}: {}'.format(prev['decided_by'], prev['decided_at'],
                                                           prev['decision_reason'])
        if same and prev['state'] == 'pending' and prev['verdict'] == j['verdict'] and \
                prev['reasons'] == j['reasons']:
            self.journal.touch_verdict(prev['id'])
            return prev, 'unchanged since verdict {} of {}'.format(prev['id'], prev['created_at'])
        if prev and prev['state'] == 'pending':
            self.journal.set_verdict_state(prev['id'], 'outdated', expect=('pending',))
        vid = self.journal.add_verdict(
            self.invocation, j['telescope'], j['night'], j['target'], j['candidate'], j['incumbent'], j['verdict'],
            j['reasons'], j['notes'], j['metrics'], run_json(j['run']), j['cand_fp'], j['inc_fp'],
            'promoting' if j['verdict'] == 'PROMOTE' else 'pending')
        return self.journal.verdict(vid), None

    def close_gone(self, tel, night, present):
        """Pending verdicts of the night whose candidate is no longer in v3 (promoted by T12, rerun, removed)."""
        for v in self.journal.verdicts(states=['pending'], telescope=tel, night=night):
            if v['target'] not in present or not os.path.isdir(v['candidate']):
                self.journal.set_verdict_state(v['id'], 'outdated', expect=('pending',))

    # ------------------------------------------------------------------ planning

    def plan(self, v, stamp):
        tel, night, target = v['telescope'], v['night'], v['target']
        run = run_from_json(v['run'])
        rec = self.restage_record(tel, night, target)
        if rec is not None and os.path.abspath(v['candidate']) == os.path.abspath(rec['dest']) and \
                rec['source'] == 'v2':
            return self.plan_overlay(v, stamp, rec, run)
        dest = self.tree.output('v2', tel, night, target)
        c = scan_dir(v['candidate'], run['start'] if run and run.get('start') else None)
        c['stale'], _ = self.run_stale(c['stale'], run, rec, night, target)
        i = scan_dir(dest)
        # files the run did not write go aside, so v2 holds only the promoted run's own (an accept of an S1 hold,
        # or the policy's)
        aside = {'dir': self.tree.superseded(tel, night, target, stamp + STALE), 'files': c['stale']} \
            if c['stale'] else None
        return {'stamp': stamp, 'candidate': v['candidate'], 'dest': dest, 'cand_fp': v['cand_fp'],
                'cand_ino': c['ino'], 'inc_fp': v['inc_fp'], 'inc_ino': i['ino'],
                'superseded': self.tree.superseded(tel, night, target, stamp) if i['ino'] else None, 'aside': aside,
                'reduction': self.plan_reduction(tel, night, run, stamp), 'records': self.plan_records(tel, night, run),
                'run': v['run'],
                'record': {'path': rec['_path'], 'dir': self.tree.superseded(tel, night, target, stamp + RESTAGED)}
                if rec is not None else None}

    def plan_overlay(self, v, stamp, rec, run):
        """A candidate restaged from v2 (utils/restage.py): v2 keeps its directory, and only the files the rerun wrote
        (its manifest: everything in the candidate but the unchanged inputs, less any that predate the run) replace
        or join v2's. v2's files they replace, and v2's outputs of the rerun's stages that the rerun did not write
        again, are set aside under superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-overlay/; the inputs go to
        <stamp>-restage/ with the record."""
        tel, night, target = v['telescope'], v['night'], v['target']
        dest = self.tree.output('v2', tel, night, target)
        c = scan_dir(v['candidate'])
        i = scan_dir(dest)
        m = restage_mod.manifest(rec, run['start'] if run and run.get('start') else None)
        early = set(m['early'])
        install = [{k: e[k] for k in ('rel', 'size', 'mtime_ns', 'sha1')} for e in m['written']
                   if e['rel'] not in early]
        names = {e['rel'] for e in install}
        have = set(restage_mod.walk_files(dest)) if os.path.isdir(dest) else set()
        orphans = sorted(rel for rel in have - names
                         if stage_map.is_output(rel, rec['from_stage'], night=night, target=target))
        records = self.plan_records(tel, night, run)
        covered = self.covers_night(tel, night, run)
        for r in records:
            if r['kind'] == 'pdf' and r.get('src'):
                # the night's report shows the rerun's curves: it takes the plain name when it covers every target
                r['replace_plain'] = covered
                r['plain_aside'] = os.path.join(self.tree.po, SUPERSEDED, tel, night, 'reports', stamp)
        return {'mode': 'overlay', 'stamp': stamp, 'candidate': v['candidate'], 'dest': dest,
                'cand_fp': v['cand_fp'], 'cand_ino': c['ino'], 'inc_fp': v['inc_fp'], 'inc_ino': i['ino'],
                'install': install, 'replace': sorted(names & have), 'orphans': orphans, 'early': sorted(early),
                'inputs': len(m['inputs']), 'from_stage': rec['from_stage'],
                'overlay_dir': self.tree.superseded(tel, night, target, stamp + OVERLAY),
                'leftover': self.tree.superseded(tel, night, target, stamp + RESTAGED),
                'record': {'path': rec['_path'], 'dir': self.tree.superseded(tel, night, target, stamp + RESTAGED)},
                'links': sorted(e['rel'] for e in rec['files'] if e['how'] == 'link'),
                'superseded': None, 'aside': None,
                'reduction': {'action': 'keep', 'why': 'an overlay leaves the reduction alone: the rerun only read it'},
                'records': records, 'run': v['run']}

    def covers_night(self, tel, night, run):
        """Whether a rerun given its targets covered every target of the night with products (a v2 or v3 directory
        with the night's *_output.fits), so that its PDF can stand for the night. A run that discovered its targets
        could have rerun some of them on nothing, so it never does."""
        if run is None or not run.get('specified'):
            return False
        made = set()
        for version in ('v2', 'v3'):
            d = self.tree.output(version, tel, night)
            for name in os.listdir(d) if os.path.isdir(d) else []:
                if name != 'reduction' and glob.glob(os.path.join(d, name, '*_{}_output.fits'.format(night))):
                    made.add(norm_target(name))
        return made <= {norm_target(t) for t in run.get('targets') or ()}

    def plan_reduction(self, tel, night, run, stamp):
        src = self.tree.output('v3', tel, night, 'reduction')
        dst = self.tree.output('v2', tel, night, 'reduction')
        r = {'action': 'keep', 'src': src, 'dst': dst, 'superseded': None, 'old_ino': lino(dst),
             'tmp': os.path.join(os.path.dirname(dst), '.reduction.promote-' + stamp)}
        new, old = reduction_info(src), reduction_info(dst)
        v2_state = 'v2 has none for the night' if not old['exists'] else \
            "v2's has no masters" if not old['masters'] else 'v2 keeps its own'
        if not new['exists']:
            r['why'] = 'there is no v3 reduction; ' + v2_state
            return r
        if not new['masters'] or not new['lists']:
            r['why'] = 'the v3 reduction has no {}; {}'.format('masters' if not new['masters'] else '1_image_*.list',
                                                               v2_state)
            return r
        if run is None:
            r['why'] = 'the run is unknown; ' + v2_state
            return r
        lo, hi = epoch(run['start']) - TOLERANCE, epoch(run['end']) + 60 if run.get('end') else None
        made = [new['files'][m][2] for m in new['masters']]
        if min(made) < lo or (hi is not None and max(made) > hi):
            r['why'] = 'the v3 reduction was not made by this run ({}): its masters changed {}; {}'.format(
                describe_run(run), iso(dt.datetime.fromtimestamp(max(made), dt.timezone.utc).replace(
                    tzinfo=None, microsecond=0)), v2_state)
            return r
        if old['exists'] and all(old['files'].get(m, [None])[:2] == new['files'][m][:2] for m in new['masters']) \
                and set(old['lists']) == set(new['lists']):
            r['why'] = "v2 already has this run's reduction"
            return r
        if not old['exists']:
            r['action'], r['why'] = 'install', 'v2 has no reduction for the night'
        elif not old['masters']:
            r['action'], r['why'] = 'replace', "v2's reduction has no masters"
        else:
            night_targets = set(list_targets(new))
            ran = None if not run['specified'] else {norm_target(t) for t in run['targets'] or ()}
            missed = sorted(night_targets - ran) if ran is not None else []
            if missed:
                r['why'] = 'the run covered {} of the {} targets of the night (not {}); v2 keeps its reduction'.format(
                    len(night_targets) - len(missed), len(night_targets), ' '.join(missed))
                return r
            r['action'], r['why'] = 'replace', 'the run covered every target of the night ({})'.format(
                ' '.join(sorted(night_targets)))
        if old['exists']:
            r['superseded'] = self.tree.superseded(tel, night, 'reduction', stamp)
        return r

    def plan_records(self, tel, night, run):
        """The logs of the run and of the later runs that made no light curves (an --only_T11 run writes the PDF),
        and the night's PDF if one of them wrote it."""
        out = []
        if run is None:
            return out
        runs = [run] + [run_from_json(r) for r in run.get('later') or ()]
        logs = self.tree.logs('v2', tel)
        for r in runs:
            stamp = run_stamp(r)
            log = r.get('log') if r.get('log') and os.path.isfile(r['log']) else None
            job_log = r.get('job_log') if r.get('job_log') and os.path.isfile(r['job_log']) else None
            if log and os.path.dirname(os.path.abspath(log)) == os.path.abspath(logs):
                out.append({'kind': 'log', 'missing': 'the log is in v2 already ({})'.format(os.path.basename(log))})
            elif log or job_log:
                if log:
                    stem = COPY_STAMP.sub('', os.path.basename(log)[:-len('.log')])
                else:
                    stem = '{}_{}_v3{}'.format(night, RUNNAME, '_' + '_'.join(r['targets']) if r['specified'] else '')
                out.append({'kind': 'log', 'src': log or job_log,
                            'dst': os.path.join(logs, '{}.{}.log'.format(stem, stamp)),
                            'plain': os.path.join(logs, stem + '.log')})
            else:
                out.append({'kind': 'log', 'missing': 'no log of {} was found'.format(describe_run(r))})
        pdf = os.path.join(self.tree.reports('v3', tel), '{}_{}.pdf'.format(tel, night))
        reports = self.tree.reports('v2', tel)
        if not os.path.isfile(pdf):
            out.append({'kind': 'pdf', 'missing': 'there is no {}'.format(pdf)})
            return out
        t = os.path.getmtime(pdf)
        for r in reversed(runs):
            end = epoch(r['end']) + 60 if r.get('end') else None
            if epoch(r['start']) - TOLERANCE <= t and (end is None or t <= end):
                out.append({'kind': 'pdf', 'src': pdf,
                            'dst': os.path.join(reports, '{}_{}.{}.pdf'.format(tel, night, run_stamp(r))),
                            'plain': os.path.join(reports, '{}_{}.pdf'.format(tel, night))})
                return out
        out.append({'kind': 'pdf', 'missing': '{} was written by another run'.format(os.path.basename(pdf))})
        return out

    # ------------------------------------------------------------------ promoting

    def promote(self, v, how, by=None, reason=None):
        """The PROMOTE steps for a verdict: an automatic one (state promoting), one a person accepted (pending), or
        one the approval policy accepted (pending, how 'policy'). Returns (state, error): state is done, aborted
        (nothing moved) or failed."""
        stamp = [None]

        def make_plan(pid):
            stamp[0] = '{}-p{}'.format(self.now().strftime('%Y%m%dT%H%M%SZ'), pid)
            return self.plan(v, stamp[0])

        pid = self.journal.start_promotion(v['id'], self.invocation, v['telescope'], v['night'], v['target'], how,
                                           make_plan)
        if pid is None:
            return 'aborted', 'verdict {} is no longer open'.format(v['id'])
        if how in ('accept', 'policy'):
            self.journal.set_verdict_state(v['id'], 'promoting', decided_by=by, decided_at=iso(self.now()),
                                           decision_reason=reason)
        return self.run_steps(self.journal.promotion(pid))

    # ------------------------------------------------------------------ the approval policy

    def apply_policy(self, j, v):
        """Ask the approval policy (promote_policy.py) about a pending HOLD or REJECT and act on its answer: promote
        it, reject it, or leave it pending for reprocessing or for a person. Every answer is journalled. Returns
        the action text for the night's listing."""
        run = j['run']
        log = None
        if run:
            log = run.get('log') if run.get('log') and os.path.isfile(run['log']) else run.get('job_log')
        losses = promote_policy.frame_losses(log, j['target'], dict(promote_policy.THRESHOLDS, **self.policy_th))
        d = promote_policy.decide(j['verdict'], j['metrics'], j['reasons'], bool(j['incumbent']),
                                  j['stale'] if run else None, j['cand_files'], losses,
                                  dict(self.th, **self.policy_th))
        evidence = {'losses': losses, 'stale': promote_policy.stale_summary(j['stale']) if run else None,
                    'files': j['cand_files']}
        if not self.apply:
            return 'policy would {} ({}{}): {}'.format(d.outcome, d.rule, ', rerun' if d.rerun else '', d.why)
        self.journal.add_policy_decision(v['id'], self.invocation, promote_policy.VERSION, d.outcome, d.rule, d.why,
                                         evidence, rerun=d.rerun)
        if d.outcome == promote_policy.PROMOTE:
            state, err = self.promote(v, 'policy', by=d.by, reason=d.why)
            p = self.journal.promotion(self.journal.verdict(v['id'])['promotion_id'] or 0)
            j['plan'] = p['plan'] if p else None
            j['error'] = state == 'failed'
            return ('promoted by policy ({}): {}'.format(d.rule, d.why) if state == 'done'
                    else 'policy promotion {}: {}'.format(state, err))
        if d.outcome == promote_policy.REJECT:
            self.journal.set_verdict_state(v['id'], 'rejected', expect=('pending',), decided_by=d.by,
                                           decided_at=iso(self.now()), decision_reason=d.why)
            return 'rejected by policy ({}): {}'.format(d.rule, d.why)
        return 'left for {} ({}): {}'.format('reprocessing' if d.outcome == promote_policy.REPROCESS else 'a person',
                                            d.rule, d.why)

    def run_steps(self, prom):
        """Run a promotion's remaining steps. Returns (state, why): done, aborted (nothing moved), failed, or
        waiting (a resumed promotion whose night is busy now; nothing changed)."""
        pid, plan = prom['id'], prom['plan']
        overlay = plan.get('mode') == 'overlay'
        done = {s['step'] for s in self.journal.steps(pid) if s['state'] in ('done', 'skipped')}
        if 'check' in done:
            # resuming: the pre-flight check ran before the interruption, and the night may be busy again since
            busy = night_busy(self.queue, prom['telescope'], prom['night'], self.own_job, self.proc)
            if busy:
                return 'waiting', 'the night is busy: ' + busy
        for name in steps_of(plan):
            if name in done:
                continue
            sid = self.journal.step_intent(pid, name, {})
            try:
                fn = getattr(self, 'step_overlay_' + name, None) if overlay else None
                detail = (fn or getattr(self, 'step_' + name))(plan, prom)
            except Exception as e:
                err = '{} step: {}: {}'.format(name, type(e).__name__, e)
                note = {'error': err}
                if name == 'check':
                    # nothing moved: an automatic or policy verdict is left to the next judgement, an accepted one to
                    # the person
                    self.journal.step_end(sid, 'failed', note)
                    self.journal.finish_promotion(pid, 'aborted', error=err,
                                                  verdict_state='pending' if prom['how'] == 'accept' else 'outdated')
                    return 'aborted', err
                if name == 'publish' and 'supersede' in done:
                    note['restored'] = self.undo_supersede(plan)
                self.journal.step_end(sid, 'failed', note)
                self.journal.finish_promotion(pid, 'failed', error=err + (
                    '; v2 restored from {}'.format(plan['superseded']) if note.get('restored') else ''))
                return 'failed', err
            self.journal.step_end(sid, 'skipped' if detail.pop('skipped', False) else 'done', detail)
            done.add(name)
        self.journal.finish_promotion(pid, 'done', verdict_state='accepted' if prom['how'] in ('accept', 'policy')
                                      else 'promoted')
        return 'done', ''

    def undo_supersede(self, plan):
        """After a failed publish: put the incumbent back, if v2 is still free and the candidate still in v3."""
        try:
            if plan['superseded'] and lino(plan['superseded']) == plan['inc_ino'] and \
                    not os.path.lexists(plan['dest']) and lino(plan['candidate']) == plan['cand_ino']:
                os.rename(plan['superseded'], plan['dest'])
                return True
        except OSError:
            pass
        return False

    def step_check(self, plan, prom):
        busy = night_busy(self.queue, prom['telescope'], prom['night'], self.own_job, self.proc)
        if busy:
            raise PromoteError('the night is busy: ' + busy)
        c = scan_dir(plan['candidate'])
        if c['fp'] != plan['cand_fp'] or c['ino'] != plan['cand_ino']:
            raise PromoteError('the candidate {} changed after it was judged'.format(plan['candidate']))
        i = scan_dir(plan['dest'])
        if (plan['superseded'] is None) != (i['ino'] is None) or i['fp'] != plan['inc_fp']:
            raise PromoteError('v2 {} changed after the candidate was judged'.format(plan['dest']))
        if os.path.islink(plan['dest']) or os.path.islink(plan['candidate']):
            raise PromoteError('a symbolic link where a directory should be')
        red = plan['reduction']
        aside = (plan.get('aside') or {}).get('dir')
        same_device(plan['candidate'], plan['dest'], plan['superseded'], red.get('dst'), red.get('superseded'), aside)
        problems = access_problems(plan)
        if problems:
            raise PromoteError('; '.join(problems))
        return {'candidate_files': c['files'], 'incumbent_files': i['files']}

    def step_aside(self, plan, prom):
        """Move the files the run did not write out of the candidate, into superseded/.../<stamp>-stale/ with their
        relative paths, so that v2 gets only the promoted run's own files."""
        aside = plan.get('aside')
        if not aside or not aside.get('files'):
            return {'skipped': True, 'why': "every file in the candidate is the run's own"}
        root = aside['dir']
        if not self.tree.under_superseded(root):
            raise PromoteError('refusing to set files aside outside {}: {}'.format(SUPERSEDED, root))
        moved = already = gone = nbytes = 0
        for rel in aside['files']:
            src, dst = os.path.join(plan['candidate'], rel), os.path.join(root, rel)
            if not os.path.lexists(src):
                if os.path.lexists(dst):
                    already += 1        # moved before an interruption
                else:
                    gone += 1
                continue
            if os.path.lexists(dst):
                raise PromoteError('{} is set aside already, and {} is still in the candidate'.format(dst, rel))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            nbytes += os.lstat(src).st_size
            os.rename(src, dst)
            moved += 1
        write_note(root, {'set_aside_at': iso(self.now()), 'from': plan['candidate'], 'promotion': prom['id'],
                          'verdict': prom['verdict_id'], 'why': 'these files predate the run that was promoted',
                          'run': plan.get('run')})
        return {'dir': root, 'moved': moved, 'bytes': nbytes, 'resumed': already, 'gone': gone}

    def step_supersede(self, plan, prom):
        src, dst = plan['dest'], plan['superseded']
        if dst is None:
            return {'skipped': True, 'why': 'v2 had no directory for the target'}
        if lino(dst) == plan['inc_ino']:
            return {'moved': [src, dst], 'resumed': True}
        if lino(src) != plan['inc_ino'] or not os.path.isdir(src):
            raise PromoteError('the v2 incumbent {} is not the one that was judged'.format(src))
        claim_and_rename(src, dst)
        write_note(dst, {'superseded_at': iso(self.now()), 'was': src, 'promotion': prom['id'],
                         'verdict': prom['verdict_id'], 'replaced_by_run': plan.get('run')})
        return {'moved': [src, dst]}

    def step_publish(self, plan, prom):
        src, dst = plan['candidate'], plan['dest']
        if lino(dst) == plan['cand_ino']:
            return {'moved': [src, dst], 'resumed': True}
        if os.path.lexists(dst):
            raise PromoteError('{} exists: v2 should be free for the target at this point'.format(dst))
        if lino(src) != plan['cand_ino']:
            raise PromoteError('the candidate {} is gone or replaced'.format(src))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        same_device(src, os.path.dirname(dst))
        os.rename(src, dst)
        return {'moved': [src, dst]}

    def step_touch(self, plan, prom):
        names = sorted(n for n in os.listdir(plan['dest']) if n.endswith('_5_diff.fits'))
        now = time.time()
        for name in names:
            os.utime(os.path.join(plan['dest'], name), (now, now))
        return {'touched': names}

    def step_reduction(self, plan, prom):
        r = plan['reduction']
        if r['action'] == 'keep':
            return {'skipped': True, 'why': r['why']}
        tmp, dst, sup = r['tmp'], r['dst'], r['superseded']
        if r.get('tmp_ino') and lino(dst) == r['tmp_ino']:
            return {'installed': dst, 'superseded': sup, 'why': r['why'], 'resumed': True}
        if not (r.get('tmp_ino') and lino(tmp) == r['tmp_ino']):
            if os.path.lexists(tmp):        # an incomplete copy from before a crash: ours, by its name
                if not os.path.basename(tmp).startswith('.reduction.promote-'):
                    raise PromoteError('refusing to remove {}'.format(tmp))
                shutil.rmtree(tmp)
            shutil.copytree(r['src'], tmp, symlinks=True)     # copy2: masters keep their mtimes
            r['tmp_ino'] = lino(tmp)
            self.journal.update_plan(prom['id'], plan)
        if sup:
            if lino(sup) != r['old_ino']:
                if lino(dst) != r['old_ino']:
                    raise PromoteError('the v2 reduction {} changed after the promotion was planned'.format(dst))
                claim_and_rename(dst, sup)
                write_note(sup, {'superseded_at': iso(self.now()), 'was': dst, 'promotion': prom['id'],
                                 'replaced_by_run': plan.get('run')})
        elif os.path.lexists(dst):
            raise PromoteError('{} appeared after the promotion was planned'.format(dst))
        os.rename(tmp, dst)
        return {'installed': dst, 'superseded': sup, 'why': r['why'], 'files': len(os.listdir(dst))}

    def step_records(self, plan, prom):
        done = []
        for rec in plan['records']:
            if rec.get('missing') or not os.path.isfile(rec['src']):
                done.append({'kind': rec['kind'], 'missing': rec.get('missing') or 'gone: ' + rec['src']})
                continue
            path, copied = copy_no_overwrite(rec['src'], rec['dst'])
            plain = None
            if not os.path.lexists(rec['plain']):
                try:
                    os.link(path, rec['plain'])
                    plain = rec['plain']
                except FileExistsError:
                    pass
            elif rec.get('replace_plain') and not os.path.samefile(rec['plain'], path):
                plain = self.replace_plain(rec['plain'], path, rec['plain_aside'])
            done.append({'kind': rec['kind'], 'copy': path, 'copied': copied, 'plain': plain})
        return {'records': done}

    def replace_plain(self, plain, copy, aside_dir):
        """Point a plain record name (<TEL>_<NIGHT>.pdf) at a new copy. The file it named is kept: under its own .run
        name if it has one (another hard link), else set aside into aside_dir first."""
        if os.lstat(plain).st_nlink < 2:
            if not self.tree.under_superseded(aside_dir):
                raise PromoteError('refusing to set {} aside outside {}'.format(plain, SUPERSEDED))
            os.makedirs(aside_dir, exist_ok=True)
            kept = os.path.join(aside_dir, os.path.basename(plain))
            if not os.path.lexists(kept):
                os.link(plain, kept)
        tmp = os.path.join(os.path.dirname(plain), '.{}.promote-{}'.format(os.path.basename(plain), os.getpid()))
        if os.path.lexists(tmp):
            os.unlink(tmp)
        os.link(copy, tmp)
        os.replace(tmp, plain)
        return plain

    def step_thin(self, plan, prom):
        out = {}
        red = plan['reduction'].get('superseded') if plan['reduction'].get('action') != 'keep' else None
        aside = (plan.get('aside') or {}).get('dir')
        # an overlay of a rerun from T9 replaced v2's frames and their per-frame files: the old ones, under
        # <stamp>-overlay/1/, go as a superseded copy's do
        for key, path in (('target', plan['superseded']), ('reduction', red), ('aside', aside),
                          ('restaged', plan.get('leftover')), ('overlay', plan.get('overlay_dir'))):
            if not path or not os.path.isdir(path):
                continue
            if key == 'overlay' and not any(n.isdigit() for n in os.listdir(path)):
                continue                    # an overlay from T10 on replaced no frames
            if not self.tree.under_superseded(path):
                raise PromoteError('refusing to delete anything outside {}: {}'.format(SUPERSEDED, path))
            out[key] = thin(path, reduction=(key == 'reduction'))
        if not out:
            return {'skipped': True, 'why': 'nothing was superseded'}
        return out

    def step_unstage(self, plan, prom):
        """After a restaged candidate is published whole (restaged in place in v3): its record goes to the superseded
        tree, so that the night is no longer taken for one restaged."""
        r = plan.get('record')
        if not r:
            return {'skipped': True, 'why': 'the candidate was not restaged'}
        dst = os.path.join(r['dir'], 'restage.json')
        if not os.path.lexists(r['path']):
            return {'record': dst, 'resumed': True}
        os.makedirs(r['dir'], exist_ok=True)
        os.rename(r['path'], dst)
        _rmdir_empty(os.path.dirname(r['path']))
        return {'record': dst}

    # ------------------------------------------------------------------ the overlay of a restaged rerun

    def step_overlay_check(self, plan, prom):
        busy = night_busy(self.queue, prom['telescope'], prom['night'], self.own_job, self.proc)
        if busy:
            raise PromoteError('the night is busy: ' + busy)
        c = scan_dir(plan['candidate'])
        if c['fp'] != plan['cand_fp'] or c['ino'] != plan['cand_ino']:
            raise PromoteError('the candidate {} changed after it was judged'.format(plan['candidate']))
        i = scan_dir(plan['dest'])
        if i['fp'] != plan['inc_fp'] or i['ino'] != plan['inc_ino'] or i['ino'] is None:
            raise PromoteError('v2 {} changed after the candidate was judged'.format(plan['dest']))
        rec = restage_mod.load_record(self.tree.base, prom['telescope'], prom['night'], prom['target'])
        if rec is None or rec['_path'] != plan['record']['path']:
            raise PromoteError('the restage record {} is gone'.format(plan['record']['path']))
        st = restage_mod.check(rec)
        bad = st['changed'] + st['missing'] + st['src_changed']
        if bad:
            raise PromoteError('{} restaged inputs (or the v2 files they came from) changed: {}'.format(
                len(bad), ', '.join(bad[:3])))
        for e in plan['install']:
            p = os.path.join(plan['candidate'], e['rel'])
            if os.path.getsize(p) != e['size'] or os.lstat(p).st_mtime_ns != e['mtime_ns']:
                raise PromoteError('{} is not the file the manifest lists'.format(p))
        if os.path.islink(plan['dest']) or os.path.islink(plan['candidate']):
            raise PromoteError('a symbolic link where a directory should be')
        same_device(plan['candidate'], plan['dest'], plan['overlay_dir'], plan['leftover'])
        problems = overlay_access_problems(plan)
        if problems:
            raise PromoteError('; '.join(problems))
        return {'install': len(plan['install']), 'replace': len(plan['replace']), 'orphans': len(plan['orphans'])}

    def step_install(self, plan, prom):
        """Install the rerun's files into v2. Every v2 file it replaces, and every v2 output of the rerun's stages it
        did not write again, gets a hard link under <stamp>-overlay/ first, so v2 never lacks a file and nothing is
        lost; then each new file is renamed over v2's (one filesystem: atomic) and read back against its sha1."""
        root, dest, cand = plan['overlay_dir'], plan['dest'], plan['candidate']
        if not self.tree.under_superseded(root):
            raise PromoteError('refusing to set files aside outside {}: {}'.format(SUPERSEDED, root))
        kept = installed = resumed = 0
        for rel in plan['replace'] + plan['orphans']:
            v2f, a = os.path.join(dest, rel), os.path.join(root, rel)
            if os.path.lexists(a):
                continue                                    # kept before an interruption
            if not os.path.lexists(v2f):
                raise PromoteError('{} went from v2 during the overlay'.format(v2f))
            os.makedirs(os.path.dirname(a), exist_ok=True)
            os.link(v2f, a)
            kept += 1
        for e in plan['install']:
            src, dst = os.path.join(cand, e['rel']), os.path.join(dest, e['rel'])
            if not os.path.lexists(src):
                if os.path.lexists(dst) and os.path.getsize(dst) == e['size'] and file_sha1(dst) == e['sha1']:
                    resumed += 1                            # installed before an interruption
                    continue
                raise PromoteError('{} went from the candidate during the overlay'.format(src))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            os.replace(src, dst)
            if os.path.getsize(dst) != e['size'] or file_sha1(dst) != e['sha1']:
                raise PromoteError('{} does not read back as the file the manifest lists'.format(dst))
            installed += 1
        gone = 0
        for rel in plan['orphans']:
            v2f = os.path.join(dest, rel)
            if os.path.lexists(v2f) and os.path.samefile(v2f, os.path.join(root, rel)):
                os.unlink(v2f)                              # its data stays under <stamp>-overlay/
                gone += 1
        if plan['replace'] or plan['orphans']:
            write_note(root, {'overlaid_at': iso(self.now()), 'v2': dest, 'promotion': prom['id'],
                              'verdict': prom['verdict_id'], 'from_stage': plan['from_stage'],
                              'why': "v2's files that a rerun from T{} replaced, or wrote no more".format(
                                  plan['from_stage']),
                              'replaced': plan['replace'], 'no_longer_written': plan['orphans'],
                              'manifest': plan['install'], 'run': plan.get('run')})
        return {'installed': installed, 'resumed': resumed, 'kept': kept, 'removed_orphans': gone,
                'dir': root if (plan['replace'] or plan['orphans']) else None}

    def step_leftover(self, plan, prom):
        """The restaged directory, now holding only the inputs (and anything that predated the run), goes to
        <stamp>-restage/ with its record; its hard links to v2's own files are dropped there (v2 keeps them)."""
        src, dst = plan['candidate'], plan['leftover']
        if lino(dst) != plan['cand_ino']:
            if lino(src) != plan['cand_ino']:
                raise PromoteError('the restaged directory {} is gone or replaced'.format(src))
            claim_and_rename(src, dst)
        rec_path, kept = plan['record']['path'], os.path.join(dst, 'restage.json')
        rec = None
        if os.path.lexists(rec_path):
            with open(rec_path) as f:
                rec = json.load(f)
            os.rename(rec_path, kept)
            _rmdir_empty(os.path.dirname(rec_path))
        elif os.path.lexists(kept):
            with open(kept) as f:
                rec = json.load(f)
        dropped = restage_mod.drop_links(dst, rec) if rec else 0
        write_note(dst, {'restaged_inputs_of': plan['dest'], 'promotion': prom['id'], 'verdict': prom['verdict_id'],
                         'why': 'the inputs a rerun from T{} read, restaged from v2; its own files went into v2'.format(
                             plan['from_stage']), 'links_to_v2_dropped': dropped})
        return {'moved': [src, dst], 'links_dropped': dropped}

    # ------------------------------------------------------------------ one night

    def night(self, tel, night, targets, only=None):
        """Judge (and with apply, act on) the candidates of one night. Returns one result dict per candidate."""
        present = set(targets)
        if only:
            want = {norm_target(t) for t in only}
            targets = [t for t in targets if norm_target(t) in want]
        busy = night_busy(self.queue, tel, night, self.own_job, self.proc)
        if busy:
            return [self._result(tel, night, t, action='skipped: night busy ({})'.format(busy), busy=True)
                    for t in targets]
        if self.apply:
            self.close_gone(tel, night, present)
        if not targets:
            return []
        if self.apply:
            for prom in self.journal.promotions(states=['in_progress'], telescope=tel, night=night):
                state, err = self.run_steps(prom)
                self.out('{} {} {}: resumed promotion {}: {}{}'.format(tel, night, prom['target'], prom['id'], state,
                                                                    ' ({})'.format(err) if err else ''))
            stuck = self.journal.promotions(states=['in_progress', 'failed'], telescope=tel, night=night)
            if stuck:
                why = 'promotion {} ({}) did not finish; see promote.py show'.format(stuck[0]['id'], stuck[0]['target'])
                return [self._result(tel, night, t, action='skipped: ' + why) for t in targets]
            # a resumed promotion has taken its candidate out of v3
            targets = [t for t in targets if glob.glob(os.path.join(self.tree.output('v3', tel, night, t),
                                                                    '*_diff.fits'))]
        runs = discover_runs(self.tree, self.queue, tel, night)
        results = []
        for target in targets:
            waiting = self.rerun_pending(tel, night, target, runs)
            if waiting:
                results.append(self._result(tel, night, target, action='skipped: ' + waiting))
                continue
            try:
                j = self.judge(tel, night, target, runs)
            except Exception as e:      # one unreadable directory must not stop the night
                results.append(self._result(tel, night, target, action='error: {}: {}'.format(type(e).__name__, e),
                                            error=True))
                continue
            j['action'], j['verdict_id'], j['plan'] = '', None, None
            if self.apply:
                v, why = self.record(j)
                j['verdict_id'] = v['id']
                if why and not (self.policy and v['state'] == 'pending' and self._policy_due(v)):
                    j['action'] = 'left: ' + why
                elif j['verdict'] == 'PROMOTE':
                    state, err = self.promote(v, 'auto')
                    p = self.journal.promotion(self.journal.verdict(v['id'])['promotion_id'] or 0)
                    j['plan'] = p['plan'] if p else None
                    j['action'] = 'promoted' if state == 'done' else '{}: {}'.format(state, err)
                    j['error'] = state == 'failed'
                elif self.policy:
                    j['action'] = 'verdict {}: {}'.format(v['id'], self.apply_policy(j, v))
                else:
                    j['action'] = 'recorded {} as verdict {}'.format(j['verdict'], v['id'])
            else:
                j['action'] = self._dry_action(j)
            results.append(j)
            self.out('{:8} {} {:16} {:7} {}'.format(tel, night, target, j['verdict'], j['action']))
        return results

    def rerun_pending(self, tel, night, target, runs):
        """Why a restaged candidate is not judged yet ('' if it is): its restage did not finish, or no run covering
        it has started since (a candidate restaged for a rerun from T11 holds T10's curves from the start)."""
        rec = self.restage_record(tel, night, target)
        if rec is None:
            return ''
        if rec.get('state') != 'staged':
            return 'the restage of {} for a rerun from T{} did not finish (utils/restage.py stage again)'.format(
                target, rec['from_stage'])
        last = latest_run(runs, target)
        staged = parse_iso(rec.get('staged_at') or rec['made_at'])
        if last is None or epoch(last['start']) < epoch(staged) - TOLERANCE:
            return 'restaged at {} for a rerun from T{}, which has not run yet'.format(
                rec.get('staged_at') or rec['made_at'], rec['from_stage'])
        return ''

    def _policy_due(self, v):
        """Whether the policy should decide a pending verdict the journal already held unchanged: not when this
        version of it has already left the verdict for reprocessing or for a person."""
        prev = self.journal.latest_policy_decision(v['id'])
        return not (prev and prev['version'] == promote_policy.VERSION)

    def _dry_action(self, j):
        if j['verdict'] != 'PROMOTE':
            prev = self.journal.latest_verdict(j['telescope'], j['night'], j['target']) if self.journal else None
            if prev and prev['state'] == 'rejected' and prev['cand_fp'] == j['cand_fp']:
                return 'would leave: rejected by {} on {}'.format(prev['decided_by'], prev['decided_at'])
            if self.policy:
                return 'would record; ' + self.apply_policy(j, None)
            return 'would record for review'
        plan = self.plan({'telescope': j['telescope'], 'night': j['night'], 'target': j['target'],
                          'candidate': j['candidate'], 'cand_fp': j['cand_fp'], 'inc_fp': j['inc_fp'],
                          'run': run_json(j['run'])}, 'DRYRUN')
        key = (j['telescope'], j['night'])
        red = plan['reduction']
        if red['action'] != 'keep' and key in self._dry_reductions:
            red = dict(red, action='keep', why="v2 gets this run's reduction from the promotion of {}".format(
                self._dry_reductions[key]))
            plan['reduction'] = red
        elif red['action'] != 'keep':
            self._dry_reductions[key] = j['target']
        j['plan'] = plan
        if plan.get('mode') == 'overlay':
            problems = overlay_access_problems(plan)
            return 'would overlay v2 with the {} files the rerun from T{} wrote ({} replacing v2\'s, {} new; {} of ' \
                'v2\'s no longer written set aside); the {} restaged inputs go to the superseded tree{}'.format(
                    len(plan['install']), plan['from_stage'], len(plan['replace']),
                    len(plan['install']) - len(plan['replace']), len(plan['orphans']), plan['inputs'],
                    '; BUT {}'.format('; '.join(problems)) if problems else '')
        problems = access_problems(plan)
        return 'would promote{}; reduction: {} ({}){}'.format(
            ', superseding the v2 directory' if plan['superseded'] else ' (new in v2)', red['action'], red['why'],
            '; BUT {}'.format('; '.join(problems)) if problems else '')

    def _result(self, tel, night, target, action, busy=False, error=False):
        self.out('{:8} {} {:16} {:7} {}'.format(tel, night, target, '-', action))
        return {'telescope': tel, 'night': night, 'target': target, 'verdict': '', 'action': action, 'busy': busy,
                'error': error, 'reasons': '', 'notes': '', 'metrics': {}, 'run': None, 'plan': None}


def steps_of(plan):
    """The steps of a promotion: the whole-directory ones, or an overlay's."""
    return OVERLAY_STEPS if (plan or {}).get('mode') == 'overlay' else STEPS


def _rmdir_empty(path):
    try:
        os.rmdir(path)
    except OSError:
        pass


def overlay_access_problems(plan):
    """What would stop an overlay: writing in v2's directory and its subdirectories, renaming the restaged
    directory out of v3, and making the directories under superseded/."""
    want = [plan['dest'], os.path.dirname(plan['candidate']), plan['candidate'], nearest_existing(plan['overlay_dir']),
            nearest_existing(plan['leftover'])]
    want += sorted({os.path.join(plan['dest'], os.path.dirname(e['rel'])) for e in plan['install']})
    want += sorted({os.path.join(plan['candidate'], os.path.dirname(e['rel'])) for e in plan['install']})
    want += sorted({os.path.join(plan['dest'], os.path.dirname(rel)) for rel in plan['orphans']})
    return ['cannot write in {}'.format(d) for d in want if os.path.isdir(d) and not os.access(d, os.W_OK | os.X_OK)]


def access_problems(plan):
    """What would stop the renames: each parent directory written to, and each directory moved to a new parent
    (Linux rewrites its '..' entry, so the directory itself must be writable too)."""
    red = plan['reduction']
    want = [os.path.dirname(plan['candidate']), plan['candidate'], nearest_existing(os.path.dirname(plan['dest']))]
    if plan['superseded']:
        want += [plan['dest'], nearest_existing(plan['superseded'])]
    if (plan.get('aside') or {}).get('files'):
        want += [nearest_existing(plan['aside']['dir'])]
        want += sorted({os.path.join(plan['candidate'], os.path.dirname(rel)) for rel in plan['aside']['files']})
    if red.get('action') != 'keep':
        want.append(nearest_existing(os.path.dirname(red['dst'])))
        if red.get('superseded'):
            want += [red['dst'], nearest_existing(red['superseded'])]
    return ['cannot write in {}'.format(d) for d in want if os.path.isdir(d) and not os.access(d, os.W_OK | os.X_OK)]


def claim_and_rename(src, dst):
    """Rename directory src to dst, a name made fresh for it: dst is created empty first (failing if it exists and
    is not empty), then src renamed onto it, so a supersession can never land inside an older one."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        os.mkdir(dst)
    except FileExistsError:
        if not os.path.isdir(dst) or os.listdir(dst):
            raise PromoteError('{} already exists'.format(dst))
    same_device(src, dst)
    os.rename(src, dst)     # replaces the empty directory just made


def write_note(path, info):
    """superseded.json beside a superseded copy: what it was and what replaced it."""
    with open(os.path.join(path, 'superseded.json'), 'w') as f:
        json.dump(info, f, indent=1, sort_keys=True, default=str)


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

EXTRA = ['action', 'verdict_id', 'checker_verdict', 'run_source', 'run_job', 'run_log', 'run_start', 'stale_files',
         'cand_files', 'supersedes', 'reduction', 'reduction_why']


def write_csv(results, path):
    cols = ['tel', 'date', 'target', 'verdict', 'reasons', 'notes'] + EXTRA + \
        [c for c in promote_check.COLUMNS if c not in ('tel', 'date', 'target', 'verdict', 'reasons', 'notes')]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in results:
            row = dict(r.get('metrics') or {})
            run, plan = r.get('run') or {}, r.get('plan') or {}
            row.update(tel=r['telescope'], date=r['night'], target=r['target'], verdict=r['verdict'],
                       reasons=r['reasons'], notes=r['notes'], action=r['action'], verdict_id=r.get('verdict_id'),
                       checker_verdict=r.get('checker_verdict', ''), run_source=run.get('source', ''),
                       run_job=run.get('job') or '', run_log=run.get('log') or '',
                       run_start=iso(run['start']) if isinstance(run.get('start'), dt.datetime) else
                       run.get('start', ''),
                       stale_files=r.get('stale_files', ''), cand_files=r.get('cand_files', ''),
                       supersedes='overlay' if plan.get('mode') == 'overlay' else
                       'yes' if plan.get('superseded') else ('no' if plan else ''),
                       reduction=(plan.get('reduction') or {}).get('action', ''),
                       reduction_why=(plan.get('reduction') or {}).get('why', ''))
            w.writerow({k: promote_check._fmt(v) if v is not None else '' for k, v in row.items()})


def write_markdown(results, path, header):
    by = Counter(r['verdict'] or 'skipped' for r in results)
    lines = ['# Promotion {}'.format(header['title']), '']
    for k, v in header['facts']:
        lines.append('- {}: {}'.format(k, v))
    lines += ['', '{} target-nights: {}.'.format(len(results), ', '.join(
        '{} {}'.format(by[v], v) for v in ('PROMOTE', 'HOLD', 'REJECT', 'skipped') if by[v])), '']
    tels = sorted({r['telescope'] for r in results})
    lines += ['| Telescope | PROMOTE | HOLD | REJECT | skipped |', '|---|---:|---:|---:|---:|']
    for tel in tels:
        c = Counter(r['verdict'] or 'skipped' for r in results if r['telescope'] == tel)
        lines.append('| {} | {} | {} | {} | {} |'.format(tel, c['PROMOTE'], c['HOLD'], c['REJECT'], c['skipped']))
    lines.append('')
    for verdict, title in (('PROMOTE', 'Promote'), ('HOLD', 'Hold'), ('REJECT', 'Reject'), ('', 'Skipped')):
        rows = [r for r in results if r['verdict'] == verdict]
        if not rows:
            continue
        lines += ['## {} ({})'.format(title, len(rows)), '']
        if verdict == 'PROMOTE':
            lines += ['| Telescope | Night | Target | Points | p2p | v2 now | Reduction | Run |',
                      '|---|---|---|---:|---:|---|---|---|']
            for r in rows:
                m, plan, run = r['metrics'], r.get('plan') or {}, r.get('run') or {}
                inc = 'none' if not r.get('incumbent') else '{} points, p2p {}'.format(
                    m.get('i_valid'), promote_check._fmt(m.get('i_p2p')) if m.get('i_p2p') is not None else '-')
                red = plan.get('reduction') or {}
                lines.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
                    r['telescope'], r['night'], r['target'], m.get('c_valid'),
                    promote_check._fmt(m.get('c_p2p')) if m.get('c_p2p') is not None else '-', inc,
                    '{}: {}'.format(red.get('action', '-'), red.get('why', '')) if red else '-',
                    describe_run(run) if run else '-'))
        else:
            lines += ['| Telescope | Night | Target | Why | Action |', '|---|---|---|---|---|']
            for r in rows:
                lines.append('| {} | {} | {} | {} | {} |'.format(r['telescope'], r['night'], r['target'],
                                                                 (r['reasons'] or '').replace('|', '/'),
                                                                 r['action'].replace('|', '/')))
        lines.append('')
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _queue_defaults(config=None):
    from jobqueue.config import load_config, queue_paths
    cfg = load_config(config)
    return cfg, queue_paths(cfg)


def _by(a):
    if a.by:
        return a.by
    try:
        return getpass.getuser()
    except Exception:           # no passwd entry for the container's uid
        return 'uid{}'.format(os.getuid())


def _code():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _child(a, tel, night, extra=()):
    """Run this command again for one night, so that the night is on the command line of the process that
    changes the tree, where the dispatcher sees it as a lock."""
    argv = [sys.executable, os.path.abspath(__file__)] + list(a.argv) + ['--telescope', tel, '--night', night] + \
        list(extra)
    return subprocess.call(argv)


def _ids(specs):
    """Verdict ids from '12', '12-20' and '12,15' (any mix)."""
    out = set()
    for spec in specs:
        for part in str(spec).split(','):
            part = part.strip()
            if '-' in part:
                first, last = part.split('-', 1)
                out.update(range(int(first), int(last) + 1))
            elif part:
                out.add(int(part))
    return sorted(out)


def _spawn(a, argv, tel, night):
    """This script for one verdict or night, with the same paths, the night on its command line."""
    paths = ['--base', a.base, '--db', a.db, '--queue-db', a.queue_db] + (['--config', a.config] if a.config else []) + \
        (['--by', a.by] if a.by else [])
    return subprocess.call([sys.executable, os.path.abspath(__file__)] + list(argv) + paths +
                           ['--telescope', tel, '--night', night])


def make_promoter(a, cfg, paths, apply, journal, queue=None):
    th = dict(promote_check.THRESHOLDS)
    th.update(cfg.get('promote', {}).get('thresholds') or {})
    for k in promote_check.THRESHOLDS:
        if getattr(a, k, None) is not None:
            th[k] = getattr(a, k)
    identity = Identity(a.base)
    if a.identity == 'table':
        table = a.target_table or os.path.join(a.base, 'ml_40pc_v2.csv')
        if os.path.exists(table):
            identity = TableIdentity(a.base, table)
        else:
            print('note: --identity table, but {} does not exist: using plan files'.format(table))
    raw_path = a.raw_cache or (paths['promote_raw'] if apply or os.path.exists(paths['promote_raw']) else None)
    raw = None if a.no_raw else CachedRawFrames(a.base, raw_path, seed=a.raw_seed,
                                                readonly=not apply and not a.raw_cache)
    own = os.environ.get('JOBQUEUE_JOB_ID')
    pcfg = cfg.get('promote') or {}
    policy = bool(pcfg.get('policy', True)) and not getattr(a, 'no_policy', False)
    return Promoter(a.base, journal=journal, queue=queue or open_queue(a.queue_db), identity=identity, raw=raw,
                    th=th, apply=apply, own_job=int(own) if own and own.isdigit() else None, policy=policy,
                    policy_th=pcfg.get('policy_thresholds') or {})


def cmd_run(a, cfg, paths):
    conn = open_queue(a.queue_db)
    only = a.target
    tels, nights = ([a.telescope] if a.telescope else None), ([a.night] if a.night else None)
    if a.job:
        job = queue_job(conn, a.job)
        if job is None or job['kind'] != 'pipeline':
            print('no pipeline job {} in {}'.format(a.job, a.queue_db))
            return 2
        if (a.telescope and a.telescope != job['telescope']) or (a.night and a.night != job['night']):
            print('job {} is {} {}, not {} {}'.format(a.job, job['telescope'], job['night'], a.telescope, a.night))
            return 2
        tels, nights = [job['telescope']], [job['night']]
        only = (job['targets'] or '').split() or only
        if job['state'] == 'running':
            print('job {} is still running'.format(a.job))
            return EX_TEMPFAIL
    tree = Tree(a.base)
    found = tree.candidates(tels, nights)
    if a.apply and not (a.telescope and a.night):
        # one child per night, each with the night on its command line (see the module docstring), including nights
        # whose pending verdicts may have lost their candidate
        nights = list(found)
        if os.path.exists(a.db) and not a.job:
            j = Journal(a.db, readonly=True)
            nights += sorted({(v['telescope'], v['night']) for v in j.verdicts(states=['pending'])
                              if (not tels or v['telescope'] in tels) and (not a.night or v['night'] == a.night)}
                             - set(found))
            j.close()
        if a.csv or a.report:
            print('note: --csv and --report go with one night (or a dry run); the journal holds every verdict '
                  '(promote.py list --all)')
        rc = 0
        for tel, night in nights:
            rc = max(rc, _child(a, tel, night))
        return rc
    journal = Journal(a.db) if a.apply else (Journal(a.db, readonly=True) if os.path.exists(a.db) else None)
    inv = journal.invocation_start('run', sys.argv[1:], _code(), os.environ.get('JOBQUEUE_JOB_ID'), a.job) \
        if a.apply else None
    p = make_promoter(a, cfg, paths, a.apply, journal, queue=conn)
    p.invocation = inv
    started = time.time()
    results = []
    if a.apply and a.telescope and a.night and (a.telescope, a.night) not in found:
        p.close_gone(a.telescope, a.night, set())
    for (tel, night), targets in found.items():
        results += p.night(tel, night, targets, only)
    want = {norm_target(x) for x in only or ()}
    for tel in tels or ():
        for night in nights or ():
            for t in restage_mod.restaged_targets(a.base, tel, night):
                if t not in found.get((tel, night), ()) and (not want or norm_target(t) in want):
                    print('{:8} {} {:16} -       restaged for a rerun, but v3 holds no light curves of it to promote '
                          '(the rerun has not run, or failed): requeue it, or retire the restage with '
                          'utils/restage.py clear'.format(tel, night, t))
    counts = Counter(r['verdict'] or 'skipped' for r in results)
    acted = Counter(r['action'].split(':')[0] for r in results)
    summary = {'target_nights': len(results), 'verdicts': dict(counts), 'actions': dict(acted),
               'raw': dict(p.raw.stats) if p.raw else {}, 'seconds': round(time.time() - started, 1)}
    told = ', '.join('{} {}'.format(counts[v], v) for v in ('PROMOTE', 'HOLD', 'REJECT', 'skipped') if counts[v])
    print('{} target-nights: {}; {} mode'.format(len(results), told, 'apply' if a.apply else 'dry-run'))
    if a.csv:
        write_csv(results, a.csv)
    if a.report:
        facts = [('when', iso(utcnow())), ('mode', 'apply' if a.apply else 'dry run: nothing moved, nothing recorded'),
                 ('base', a.base), ('code', _code()), ('identity', type(p.identity).__name__),
                 ('thresholds', ', '.join('{} {}'.format(k, v) for k, v in sorted(p.th.items()))),
                 ('raw-frame counts',
                  ', '.join('{} {}'.format(k, v) for k, v in sorted(summary['raw'].items())) or '-'),
                 ('seconds', summary['seconds'])]
        write_markdown(results, a.report, {'title': 'dry run' if not a.apply else 'run', 'facts': facts})
    if journal:
        if inv:
            journal.invocation_end(inv, summary)
        journal.close()
    if any(r.get('error') for r in results):
        return 1
    if a.job and any(r.get('busy') for r in results):
        return EX_TEMPFAIL
    return 0


def _journal_or_none(a, write):
    if write:
        return Journal(a.db)
    if not os.path.exists(a.db):
        print('no journal at {} yet'.format(a.db))
        return None
    return Journal(a.db, readonly=True)


def cmd_list(a, cfg, paths):
    j = _journal_or_none(a, False)
    if j is None:
        return 0
    rows = j.verdicts(states=None if a.all or a.needs == 'rerun' else ['pending', 'promoting'],
                      telescope=a.telescope, night=a.night)
    if a.grep:
        rows = [v for v in rows if a.grep.lower() in (v['reasons'] or '').lower()]
    policy = j.latest_policy_decisions([v['id'] for v in rows])
    if a.needs == 'rerun':
        # published or not: the nights whose processing, not data, the policy found at fault
        rows = [v for v in rows if v['state'] not in ('outdated',) and (policy.get(v['id']) or {}).get('rerun')]
    elif a.needs:
        rows = [v for v in rows if v['state'] == 'pending' and
                (policy.get(v['id']) or {}).get('outcome', promote_policy.PERSON) == a.needs]
    for v in rows:
        d = policy.get(v['id'])
        print('{:>6} {:9} {:8} {} {:16} {:7} {:9} {:10} {}'.format(
            v['id'], v['state'], v['telescope'], v['night'], v['target'], v['verdict'], v['created_at'][:10],
            (d['outcome'] if d else '-') if v['state'] == 'pending' else '', (d['why'] if d else v['reasons'])[:150]))
    for p in j.promotions(states=['in_progress', 'failed']):
        print('promotion {} ({} {} {}, verdict {}) is {}: {}'.format(p['id'], p['telescope'], p['night'], p['target'],
                                                                  p['verdict_id'], p['state'], p['error'] or ''))
    print('{} verdict(s)'.format(len(rows)))
    j.close()
    return 0


def cmd_show(a, cfg, paths):
    j = _journal_or_none(a, False)
    if j is None:
        return 1
    v = j.verdict(a.id)
    if v is None:
        print('no verdict {}'.format(a.id))
        return 1
    for k in ('id', 'state', 'verdict', 'telescope', 'night', 'target', 'candidate', 'incumbent', 'created_at',
              'updated_at', 'decided_by', 'decided_at', 'decision_reason', 'promotion_id'):
        print('{:16} {}'.format(k, v[k]))
    print('{:16} {}'.format('run', describe_run(run_from_json(v['run']))))
    print('reasons:\n  ' + '\n  '.join((v['reasons'] or '-').split('; ')))
    print('notes:\n  ' + '\n  '.join((v['notes'] or '-').split('; ')))
    m = v['metrics']
    f = promote_check._fmt
    print('candidate: {} of {} points, p2p {}, {} raw frames; v2: {} points, p2p {}'.format(
        m.get('c_valid'), m.get('c_frames'), f(m.get('c_p2p')) or '-', f(m.get('raw_frames')) or '-',
        m.get('i_valid'), f(m.get('i_p2p')) or '-'))
    for d in j.policy_decisions(v['id']):
        print('policy v{} {}: {} ({}): {}'.format(d['version'], d['decided_at'], d['outcome'], d['rule'], d['why']))
        losses = (d['evidence'] or {}).get('losses')
        if losses:
            print('  plate solves: {} solved, {} lost to the sky, {} to the pipeline; {}'.format(
                losses.get('solved'), losses.get('data'), losses.get('processing'), losses.get('causes')))
    for p in j.promotions():
        if p['verdict_id'] != v['id']:
            continue
        print('promotion {} ({}, {}) {} -> {}{}'.format(p['id'], p['how'], p['state'], p['started_at'],
                                                       p['finished_at'], ': ' + p['error'] if p['error'] else ''))
        for s in j.steps(p['id']):
            print('  {:9} {:7} {}'.format(s['step'], s['state'], json.dumps(s['detail'], sort_keys=True)[:300]))
    j.close()
    return 0


def cmd_accept(a, cfg, paths):
    ids = _ids(a.ids)
    if len(ids) > 1 or (a.apply and not (a.telescope and a.night)):
        # one process per verdict, its night on the command line (where the dispatcher sees the lock)
        j = _journal_or_none(a, False)
        if j is None:
            return 1
        rc = 0
        for vid in ids:
            v = j.verdict(vid)
            if v is None:
                print('no verdict {}'.format(vid))
                rc = max(rc, 2)
                continue
            argv = ['accept', str(vid), '--reason', a.reason] + (['--apply'] if a.apply else [])
            rc = max(rc, _spawn(a, argv, v['telescope'], v['night']))
        j.close()
        return rc
    a.id = ids[0]
    j = _journal_or_none(a, a.apply)
    if j is None:
        return 1
    v = j.verdict(a.id)
    if v is None or v['state'] != 'pending':
        print('verdict {} is {}: only a pending HOLD or REJECT can be accepted'.format(
            a.id, v['state'] if v else 'absent'))
        return 2
    if a.apply and (a.telescope, a.night) != (v['telescope'], v['night']):
        j.close()
        return _spawn(a, ['accept', str(v['id']), '--reason', a.reason, '--apply'], v['telescope'], v['night'])
    p = make_promoter(a, cfg, paths, a.apply, j)
    busy = night_busy(p.queue, v['telescope'], v['night'], p.own_job, p.proc)
    if busy:
        print('not now: {}'.format(busy))
        return EX_TEMPFAIL
    c, i = scan_dir(v['candidate']), scan_dir(p.tree.output('v2', v['telescope'], v['night'], v['target']))
    if c['fp'] != v['cand_fp']:
        print('the candidate changed after verdict {} (or is gone); judge it again with promote.py run --telescope {} '
              '--night {}'.format(v['id'], v['telescope'], v['night']))
        return 2
    if i['fp'] != v['inc_fp']:
        print('v2 changed after verdict {}; judge the candidate again'.format(v['id']))
        return 2
    if not a.apply:
        plan = p.plan(v, 'DRYRUN')
        print('dry run: accepting verdict {} ({} {} {}, {}) would promote {}'.format(
            v['id'], v['telescope'], v['night'], v['target'], v['verdict'], v['candidate']))
        if plan.get('mode') == 'overlay':
            print('  as an overlay: {} files of the rerun from T{} into v2 ({} replacing, {} no longer written set '
                  'aside under {})'.format(len(plan['install']), plan['from_stage'], len(plan['replace']),
                                           len(plan['orphans']), os.path.dirname(plan['overlay_dir'])))
        print('  v2 directory: {}'.format('superseded into ' + os.path.dirname(plan['superseded']) + '/<stamp>'
                                          if plan['superseded'] else 'kept' if plan.get('mode') == 'overlay'
                                          else 'none yet'))
        print('  reduction: {} ({})'.format(plan['reduction']['action'], plan['reduction']['why']))
        for rec in plan['records']:
            print('  {}: {}'.format(rec['kind'], rec.get('missing') or 'copy to ' + rec['dst']))
        print('add --apply to do it')
        return 0
    p.invocation = j.invocation_start('accept', sys.argv[1:], _code())
    state, err = p.promote(v, 'accept', by=_by(a), reason=a.reason)
    j.invocation_end(p.invocation, {'verdict': v['id'], 'state': state, 'error': err})
    print('verdict {}: {}{}'.format(v['id'], state, ' ({})'.format(err) if err else ''))
    return 0 if state == 'done' else 1


def cmd_reject(a, cfg, paths):
    j = _journal_or_none(a, True)
    inv = j.invocation_start('reject', sys.argv[1:], _code())
    rc, done = 0, []
    for vid in _ids(a.ids):
        v = j.verdict(vid)
        if v is None or v['state'] != 'pending':
            print('verdict {} is {}: only a pending HOLD or REJECT can be rejected'.format(
                vid, v['state'] if v else 'absent'))
            rc = 2
            continue
        j.set_verdict_state(v['id'], 'rejected', expect=('pending',), decided_by=_by(a), decided_at=iso(utcnow()),
                            decision_reason=a.reason)
        done.append(v['id'])
        print('verdict {} ({} {} {}) rejected; nothing moved'.format(v['id'], v['telescope'], v['night'], v['target']))
    j.invocation_end(inv, {'verdicts': done})
    j.close()
    return rc


def cmd_resume(a, cfg, paths):
    j = _journal_or_none(a, a.apply)
    if j is None:
        return 0
    states = ['in_progress', 'failed'] if a.promotion else ['in_progress']
    proms = [p for p in j.promotions(states=states) if not a.promotion or p['id'] == a.promotion]
    if a.telescope and a.night:
        proms = [p for p in proms if (p['telescope'], p['night']) == (a.telescope, a.night)]
    if not proms:
        print('nothing to resume')
        return 0
    if not a.apply:
        for p in proms:
            done = [s['step'] for s in j.steps(p['id']) if s['state'] in ('done', 'skipped')]
            print('promotion {} ({} {} {}, {}): done {}; would run {}'.format(
                p['id'], p['telescope'], p['night'], p['target'], p['state'], ' '.join(done) or 'nothing',
                ' '.join(s for s in steps_of(p['plan']) if s not in done)))
        print('add --apply to finish them')
        return 0
    nights = sorted({(p['telescope'], p['night']) for p in proms})
    if not (a.telescope and a.night):
        j.close()
        return max(_child(a, tel, night) for tel, night in nights)
    pr = make_promoter(a, cfg, paths, True, j)
    pr.invocation = j.invocation_start('resume', sys.argv[1:], _code())
    rc = 0
    for p in proms:
        state, err = pr.run_steps(p)
        print('promotion {} ({} {} {}): {}{}'.format(p['id'], p['telescope'], p['night'], p['target'], state,
                                                   ' ({})'.format(err) if err else ''))
        rc = rc or (0 if state == 'done' else EX_TEMPFAIL if state == 'waiting' else 1)
    j.invocation_end(pr.invocation, {'promotions': [p['id'] for p in proms]})
    return rc


def cmd_abandon(a, cfg, paths):
    j = _journal_or_none(a, True)
    p = j.promotion(a.id)
    if p is None or p['state'] not in ('in_progress', 'failed'):
        print('promotion {} is {}: only an unfinished one can be abandoned'.format(a.id, p['state'] if p else 'absent'))
        return 2
    published = any(s['step'] == 'publish' and s['state'] == 'done' for s in j.steps(p['id']))
    inv = j.invocation_start('abandon', sys.argv[1:], _code())
    j.finish_promotion(p['id'], 'abandoned', error='abandoned by {}: {}'.format(_by(a), a.reason),
                       verdict_state=('accepted' if p['how'] in ('accept', 'policy') else 'promoted') if published
                       else 'outdated')
    j.invocation_end(inv, {'promotion': p['id'], 'published': published})
    print('promotion {} abandoned ({}); nothing on disk was changed by this'.format(
        p['id'], 'the candidate is in v2' if published else 'the candidate was not published'))
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog='promote.py', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    def common(s):
        s.add_argument('--config', help='queue config (default $ORCHARD_QUEUE_CONFIG, then <queue_root>/config.json)')
        s.add_argument('--base', help='data root (default: the queue config basedir)')
        s.add_argument('--db', help='the promotion journal (default <queue_root>/promote.sqlite)')
        s.add_argument('--queue-db', help='the queue database, read-only (default <queue_root>/queue.sqlite)')
        s.add_argument('--telescope')
        s.add_argument('--night')
        s.add_argument('--by', help='who decides (accept, reject, abandon; default the user)')

    def judging(s):
        s.add_argument('--identity', choices=('plan', 'table'), help='A3: schedule plan files (as before), or the K1 '
                       'target table ml_40pc_v2.csv in the data root for Sp targets, holding where the plan disagrees')
        s.add_argument('--target-table', help='with --identity table: the table, if not <base>/ml_40pc_v2.csv')
        s.add_argument('--raw-cache', help='SQLite cache of raw-frame counts (default <queue_root>/promote_raw.sqlite; '
                                           'a dry run only reads it)')
        s.add_argument('--raw-seed', help="promote_check's JSON raw cache to start from (invalidated per night)")
        s.add_argument('--no-raw', action='store_true', help='skip A4')
        for k, v in promote_check.THRESHOLDS.items():
            s.add_argument('--' + k.replace('_', '-'), type=type(v), help='default {}'.format(v))

    s = sub.add_parser('run', help='judge v3 target-nights and, with --apply, promote what passes')
    common(s)
    judging(s)
    s.add_argument('--target', action='append', help='only this target (repeatable)')
    s.add_argument('--job', type=int, help="judge the targets of this queue pipeline job's night")
    s.add_argument('--apply', action='store_true', help='act: promote, and journal holds and rejects')
    s.add_argument('--no-policy', action='store_true',
                   help='leave every HOLD and REJECT for a person (default: the approval policy, promote_policy.py, '
                        'decides them unless promote.policy is false in the queue config)')
    s.add_argument('--csv', help='write one row per target-night here')
    s.add_argument('--report', help='write a markdown summary here')
    s.set_defaults(fn=cmd_run)

    s = sub.add_parser('list', help='verdicts waiting for a person, and promotions that did not finish')
    common(s)
    s.add_argument('--all', action='store_true', help='every verdict, decided ones included')
    s.add_argument('--grep', help='only verdicts whose reasons contain this text, e.g. S1')
    s.add_argument('--needs', choices=(promote_policy.PERSON, promote_policy.REPROCESS, 'rerun'),
                   help='only pending verdicts the policy left for a person, or for reprocessing; rerun: every '
                        'verdict, published or not, whose night the policy found the processing at fault for')
    s.set_defaults(fn=cmd_list)

    s = sub.add_parser('show', help='one verdict, its reasons, and its promotion step by step')
    common(s)
    s.add_argument('id', type=int)
    s.set_defaults(fn=cmd_show)

    for name, fn, text in (('accept', cmd_accept, 'promote a held or rejected candidate (dry run without --apply)'),
                           ('reject', cmd_reject, 'close a held or rejected candidate; nothing is moved')):
        s = sub.add_parser(name, help=text)
        common(s)
        judging(s)
        s.add_argument('ids', nargs='+', help='verdict ids (promote.py list): 12, 12-20 or 12,15')
        s.add_argument('--reason', required=True)
        if name == 'accept':
            s.add_argument('--apply', action='store_true')
        s.set_defaults(fn=fn)

    s = sub.add_parser('resume', help='finish promotions a crash cut short (dry run without --apply)')
    common(s)
    judging(s)
    s.add_argument('--promotion', type=int, help='only this one (also retries a failed one)')
    s.add_argument('--apply', action='store_true')
    s.set_defaults(fn=cmd_resume)

    s = sub.add_parser('abandon', help='close an unfinished promotion after sorting it out by hand')
    common(s)
    s.add_argument('id', type=int, help='promotion id')
    s.add_argument('--reason', required=True)
    s.set_defaults(fn=cmd_abandon)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        sys.stdout.reconfigure(line_buffering=True)     # a job log shows each target-night as it is judged
    except (AttributeError, ValueError):
        pass
    a = build_parser().parse_args(argv)
    # what a child gets: the command as given, without --telescope/--night (it gets its own) or --csv/--report
    # (each child would overwrite the last one's)
    keep, skip = [], False
    for k, x in enumerate(argv):
        if skip:
            skip = False
            continue
        if x in ('--telescope', '--night', '--csv', '--report'):
            skip = True
            continue
        if x.startswith(('--telescope=', '--night=', '--csv=', '--report=')):
            continue
        keep.append(x)
    a.argv = keep
    cfg, paths = _queue_defaults(a.config)
    paths['promote_db'] = (cfg.get('promote') or {}).get('db') or os.path.join(cfg['queue_root'], 'promote.sqlite')
    paths['promote_raw'] = (cfg.get('promote') or {}).get('raw_cache') or os.path.join(cfg['queue_root'],
                                                                                     'promote_raw.sqlite')
    a.base = a.base or cfg['basedir']
    a.db = a.db or paths['promote_db']
    a.queue_db = a.queue_db or paths['db']
    if hasattr(a, 'identity') and not a.identity:
        a.identity = (cfg.get('promote') or {}).get('identity') or 'plan'
    logging.getLogger(target_management.__name__).setLevel(logging.ERROR)
    return a.fn(a, cfg, paths)


if __name__ == '__main__':
    sys.exit(main())
