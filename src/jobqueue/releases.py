"""Code releases: every job runs from a frozen copy of the deployed tree, so a deploy never waits for the queue.

A release is a directory <releases_dir>/<id>/src holding a complete copy of the production tree as it was
when the release was built, with a commit's files laid over it. Its files are read-only and never change.
`current` in <releases_dir> is a relative symlink to the release new jobs use.

When the dispatcher starts a job it pins the job to the current release: the job's working directory, its
ORCHARD_PATH and the live-tree entries of its PYTHONPATH are pointed into the release, so ZLP_pipeline.sh
(SCRIPTDIR=$ORCHARD_PATH), every Python script it runs and `python -m jobqueue jobs ...` all read the
release, from the first stage to the last. An attempt keeps its release until it ends; a retry starts on
whatever release is current then. The job's code_sha records the release.

A deploy is then:

    build     copy the live tree (minus caches, logs, core dumps) and lay the staged commit over it. Files
              identical to the previous release are hard links to it, so a release costs only what changed.
    test      run the suite from the release (optional, needs a pytest pylib), its temporary files on the data
              disk next to the release, removed after
    activate  repoint `current`: every job started from now on runs the new code; running jobs carry on,
              each on the release it started with
    sync      bring the live tree (/opt/orchard/src) up to the release, for what does not go through the
              queue (the dispatcher, the watcher, cron, people). Only files that differ are written, each
              through a temporary file and a rename; replaced files are backed up; nothing is deleted. It
              waits while a pipeline or download process still runs from the live tree (one started outside
              the queue, or by the queue before releases existed): activating a release marks the sync
              pending, and the dispatcher does it on the first poll the live tree is free.

Without a `current` release the dispatcher launches from the live tree as it always did.

Credentials are the exception to the copy (`release_link`, by default every `.env`, such as reporting/.env). A
release holds a symlink to the live tree's file instead, so the file keeps its own mode, exists once, and a
rotated password reaches every release at once; neither a sync nor a rollback to an older release ever writes one
into the live tree. Copied files keep their own read bits, so anything only its owner may read stays that way.
"""
import fnmatch
import importlib.util
import json
import os
import shutil
import stat
import subprocess
import time

from .util import iso, one_blas_thread, utcnow

CURRENT = 'current'
RELEASE_JSON = 'RELEASE.json'
TESTS_JSON = 'tests.json'
TEST_TMP = 'pytest-tmp'            # <release>/pytest-tmp: the suite's temporary files while it runs
SYNC_KEY = 'release_sync'          # kv: the live tree still has to be brought up to the current release
HISTORY = 'history.jsonl'
# the tree's own record of which orchard it is, filled in by `git archive` (utils/orchard_version.py)
VERSION_MODULE = os.path.join('utils', 'orchard_version.py')

# Processes that do pipeline work and so read code for as long as they run.
WORK_MARKERS = ('ZLP_pipeline.sh', 'SSO_download.py', 'jobqueue.runner')


class ReleaseError(Exception):
    pass


# --------------------------------------------------------------------------- locations

def releases_dir(cfg):
    return cfg.get('releases_dir') or os.path.join(cfg['queue_root'], 'releases')


def release_path(cfg, rid):
    if not rid or os.sep in rid or rid.startswith('.') or rid == CURRENT:
        raise ReleaseError('not a release id: {!r}'.format(rid))
    return os.path.join(releases_dir(cfg), rid)


def read_release(cfg, rid):
    path = release_path(cfg, rid)
    try:
        with open(os.path.join(path, RELEASE_JSON)) as f:
            info = json.load(f)
    except (OSError, ValueError):
        return None
    if not os.path.isdir(os.path.join(path, 'src')):
        return None
    info['path'] = path
    info['src'] = os.path.join(path, 'src')
    try:
        with open(os.path.join(path, TESTS_JSON)) as f:
            info['tests'] = json.load(f).get('summary')
    except (OSError, ValueError):
        pass
    return info


def current_id(cfg):
    try:
        target = os.readlink(os.path.join(releases_dir(cfg), CURRENT))
    except OSError:
        return None
    return os.path.basename(target.rstrip(os.sep))


def current(cfg):
    """The current release's info (with 'src'), or None if there is none or it is unreadable."""
    rid = current_id(cfg)
    return read_release(cfg, rid) if rid else None


def list_releases(cfg):
    d = releases_dir(cfg)
    try:
        names = sorted(n for n in os.listdir(d) if not n.startswith('.') and n != CURRENT and
                       os.path.isfile(os.path.join(d, n, RELEASE_JSON)))
    except OSError:
        return []
    return [r for r in (read_release(cfg, n) for n in names) if r]


# --------------------------------------------------------------------------- pinning a job

def _under(path, root):
    return path == root or path.startswith(root + os.sep)


def remap(path, cfg, rel_src):
    """A path into the live tree, or into any release, as the same path in rel_src; anything else unchanged."""
    if not path:
        return path
    p = os.path.normpath(path)
    live = os.path.normpath(cfg['src_dir'])
    if _under(p, live):
        return rel_src + p[len(live):]
    rdir = os.path.normpath(releases_dir(cfg))
    if _under(p, rdir):
        parts = p[len(rdir) + 1:].split(os.sep)
        if len(parts) >= 2 and parts[1] == 'src':
            return os.path.join(rel_src, *parts[2:])
    return path


def pin(cfg, cwd, env, release=None):
    """Point a job's cwd and environment at a release (default: the current one).

    Returns (cwd, env, release_id); release_id is None, and nothing changes, if there is no usable release.
    Paths outside the live tree and the releases (a scratch checkout, say) are left alone.
    """
    rel = release if release is not None else current(cfg)
    if not rel:
        return cwd, env, None
    src = rel['src']
    env = dict(env)
    new_cwd = remap(cwd, cfg, src) if cwd else src
    if not os.path.isdir(new_cwd):
        new_cwd = cwd
    env['ORCHARD_PATH'] = remap(env.get('ORCHARD_PATH') or cfg['src_dir'], cfg, src)
    pp = env.get('PYTHONPATH')
    env['PYTHONPATH'] = os.pathsep.join(remap(x, cfg, src) for x in pp.split(os.pathsep)) if pp else src
    env['JOBQUEUE_RELEASE'] = rel['id']
    env['JOBQUEUE_SRC'] = src
    return new_cwd, env, rel['id']


def code_dir(cfg):
    """The tree a job's own code should run things from: its release if pinned, else the live tree."""
    return os.environ.get('JOBQUEUE_SRC') or cfg['src_dir']


# --------------------------------------------------------------------------- building

def _excluded(rel, name, patterns):
    for pat in patterns:
        if pat.startswith('/'):
            if fnmatch.fnmatch(rel, pat[1:]):
                return True
        elif fnmatch.fnmatch(name, pat):
            return True
    return False


def tree_files(root, patterns=()):
    """{relative path: absolute path} of the regular files under root, skipping excluded names."""
    out = {}
    root = os.path.normpath(root)
    for dirpath, dirs, files in os.walk(root):
        reldir = os.path.relpath(dirpath, root)
        reldir = '' if reldir == '.' else reldir
        dirs[:] = sorted(d for d in dirs if not _excluded(os.path.join(reldir, d), d, patterns) and
                         not os.path.islink(os.path.join(dirpath, d)))
        for name in files:
            rel = os.path.join(reldir, name)
            path = os.path.join(dirpath, name)
            if _excluded(rel, name, patterns) or os.path.islink(path) or not os.path.isfile(path):
                continue
            out[rel] = path
    return out


def same_bytes(a, b, chunk=1 << 20):
    try:
        sa, sb = os.stat(a), os.stat(b)
    except OSError:
        return False
    if sa.st_size != sb.st_size:
        return False
    if (sa.st_dev, sa.st_ino) == (sb.st_dev, sb.st_ino):
        return True
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        while True:
            x, y = fa.read(chunk), fb.read(chunk)
            if x != y:
                return False
            if not x:
                return True


def _release_mode(path):
    """A release copy's mode: the source's read and execute bits, readable by its owner, writable by nobody."""
    return (stat.S_IMODE(os.stat(path).st_mode) | 0o400) & ~0o222


def _same_file(source, copy):
    """True if copy (in an earlier release) is what copying source makes now: same bytes and the same mode, so a
    file whose mode has narrowed is copied afresh, not hard-linked to the wider inode."""
    return same_bytes(source, copy) and _release_mode(source) == stat.S_IMODE(os.stat(copy).st_mode)


def _linked(rel, patterns):
    return _excluded(rel, os.path.basename(rel), patterns)


def _readonly(path):
    os.chmod(path, stat.S_IMODE(os.stat(path).st_mode) & ~0o222)


def release_id(sha, now):
    return '{}-{}'.format(now.strftime('%Y%m%dT%H%M%SZ'), (sha or 'live')[:7])


def build(cfg, stage=None, sha=None, now=None, log=print):
    """Build a release from the live tree with `stage` (a `git archive <sha> src` tree) laid over it.

    With no stage the release is a snapshot of the live tree as it is. Returns the release info. The release
    is complete or absent: it is assembled under a hidden name and renamed into place at the end.
    """
    now = now or utcnow()
    live = os.path.normpath(cfg['src_dir'])
    patterns = list(cfg.get('release_exclude') or ())
    links = list(cfg.get('release_link') or ())
    if stage:
        stage = os.path.normpath(stage)
        for probe in ('main/ZLP_pipeline.sh', 'jobqueue/__init__.py'):
            if not os.path.isfile(os.path.join(stage, probe)):
                raise ReleaseError('{} has no {}: give the src directory of the staged commit'.format(stage, probe))
        if not sha:
            raise ReleaseError('a staged tree needs --sha')
    root = releases_dir(cfg)
    os.makedirs(root, exist_ok=True)
    rid = base = release_id(sha, now)
    n = 1
    while os.path.lexists(os.path.join(root, rid)) or os.path.lexists(os.path.join(root, '.building-' + rid)):
        n += 1
        rid = '{}.{}'.format(base, n)
    final = os.path.join(root, rid)

    prev = current(cfg)
    prev_src = prev['src'] if prev else None
    live_files = tree_files(live, patterns)
    stage_files = tree_files(stage, patterns) if stage else {}
    sources = dict(live_files)
    sources.update(stage_files)
    sources.pop('VERSION', None)
    if stage and VERSION_MODULE not in stage_files:
        sources.pop(VERSION_MODULE, None)       # a commit from before versioning must not claim the live tree's
    # credentials: a symlink to the live file, not a copy (a staged one, if any, gives way to the live one)
    live_links = sorted(rel for rel in live_files if _linked(rel, links))
    for rel in live_links:
        del sources[rel]

    tmp = os.path.join(root, '.building-' + rid)
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    src = os.path.join(tmp, 'src')
    os.makedirs(src)
    linked = copied = new_bytes = 0
    changed_since_prev = []
    try:
        for rel in sorted(sources):
            source = sources[rel]
            dst = os.path.join(src, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            before = os.path.join(prev_src, rel) if prev_src else None
            if before and os.path.isfile(before) and not os.path.islink(before) and _same_file(source, before):
                os.link(before, dst)
                linked += 1
                continue
            shutil.copyfile(source, dst)
            os.chmod(dst, _release_mode(source))
            copied += 1
            new_bytes += os.path.getsize(dst)
            changed_since_prev.append(rel)
        for rel in live_links:
            dst = os.path.join(src, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            os.symlink(os.path.join(live, rel), dst)
        with open(os.path.join(src, 'VERSION'), 'w') as f:
            f.write('{} {}\n'.format(rid, sha or 'live-snapshot'))
        _readonly(os.path.join(src, 'VERSION'))
        version = tree_version(src)
        if stage and version and version['source'] == 'archive' and not version['commit'].startswith(sha.lower()):
            raise ReleaseError('{} is commit {} by its own {}, not --sha {}: stage `git archive {} src`'
                               .format(stage, version['commit'][:12], VERSION_MODULE, sha, sha))
        if stage and version and version['source'] != 'archive':
            log('warning: {} in {} was not filled in by git archive, so this release cannot say which orchard it '
                'is and its products are stamped {!r}'.format(VERSION_MODULE, stage, version['version']))

        # content only, as sync_live compares: the live tree's modes are its own (deploys never copied them)
        differs_from_live = sorted(rel for rel in stage_files if rel not in live_links and
                                   (rel not in live_files or not same_bytes(stage_files[rel], live_files[rel])))
        removed_since_prev = []
        if prev_src:
            removed_since_prev = sorted(set(tree_files(prev_src, patterns)) - set(sources) - set(live_links) -
                                        {'VERSION'})
        info = {'id': rid, 'sha': sha, 'built_at': iso(now), 'stage': stage, 'live_src': live,
                'version': version['version'] if version else None,
                'commit': version['commit'] if version else None,
                'previous': prev['id'] if prev else None, 'files': len(sources),
                'linked_from_previous': linked, 'copied': copied, 'new_bytes': new_bytes,
                'live_links': live_links,
                'differs_from_live': differs_from_live,
                'changed_since_previous': changed_since_prev if prev_src else [],
                'removed_since_previous': removed_since_prev}
        with open(os.path.join(tmp, RELEASE_JSON), 'w') as f:
            json.dump(info, f, indent=1, sort_keys=True)
            f.write('\n')
        os.rename(tmp, final)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    log('built release {} (orchard {}): {} files ({} new or changed, {:.1f} MB; {} linked to {}); {} differ from '
        'the live tree{}'
        .format(rid, info['version'] or '-', len(sources), copied, new_bytes / 1e6, linked, prev['id'] if prev else '-',
                len(differs_from_live),
                '; {} point{} at the live tree'.format(', '.join(live_links), 's' if len(live_links) == 1 else '')
                if live_links else ''))
    return read_release(cfg, rid)


def tree_version(src):
    """What a tree's own utils/orchard_version.py says it is ({'version', 'commit', 'source', ...}), or None for a
    tree without one (a commit from before versioning). The tree's module answers, so a release reports itself
    the way its jobs stamp their products."""
    path = os.path.join(src, VERSION_MODULE)
    if not os.path.isfile(path):
        return None
    try:
        spec = importlib.util.spec_from_file_location('_orchard_version_of_tree', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.compute(src=src)
    except Exception:
        return None


# --------------------------------------------------------------------------- testing

def suite_tmp_dir(cfg, rel):
    """Where a release's suite keeps its temporary files: on the data disk, never the container's /tmp. That is
    appct's 20 GB /var; a full suite left ~1.1 GB in /tmp/pytest-of-unknown, pytest keeps three runs, and on
    2026-10-19 /var filled and production jobs died with Errno 28."""
    base = cfg.get('test_tmp')
    return os.path.join(base, rel['id']) if base else os.path.join(rel['path'], TEST_TMP)


def remove_tree(path):
    """rm -rf a directory, also through read-only directories a suite leaves behind. True if it is gone."""
    def writable(p):
        try:
            os.chmod(p, stat.S_IMODE(os.lstat(p).st_mode) | stat.S_IRWXU)
        except OSError:
            pass
    if os.path.isdir(path) and not os.path.islink(path):
        writable(path)
        for dirpath, dirs, _ in os.walk(path):      # top-down: a subdirectory is made readable before it is listed
            for d in dirs:
                if not os.path.islink(os.path.join(dirpath, d)):
                    writable(os.path.join(dirpath, d))
        shutil.rmtree(path, ignore_errors=True)
    return not os.path.lexists(path)


def run_tests(cfg, rid, log=print):
    """Run the test suite from the release, as the deploy scripts did from a staged tree. Returns the exit code.

    pytest's --basetemp and TMPDIR point into suite_tmp_dir (outside the release's src, so no later build copies
    it), which is removed when the run ends, however it ends. tests.log and tests.json stay in the release.
    """
    rel = read_release(cfg, rid)
    if not rel:
        raise ReleaseError('no release {}'.format(rid))
    pylib = cfg.get('test_pylib')
    if not pylib or not os.path.isdir(pylib):
        raise ReleaseError('no pytest pylib at {!r} (config test_pylib); the container has no pytest of its own'
                           .format(pylib))
    paths = [p for p in cfg.get('test_paths') or ('tests',) if os.path.isdir(os.path.join(rel['src'], p))]
    tmp = suite_tmp_dir(cfg, rel)
    remove_tree(tmp)                    # what a killed run left
    basetemp, tmpdir = os.path.join(tmp, 'basetemp'), os.path.join(tmp, 'tmp')
    os.makedirs(tmpdir, exist_ok=True)  # pytest empties basetemp when it starts, so TMPDIR is not inside it
    env = one_blas_thread(os.environ)   # under load, OpenBLAS's default pools hung the suite (2026-10-20)
    env.update({'PYTHONPATH': os.pathsep.join([pylib, rel['src']]), 'ORCHARD_PATH': rel['src'],
                'PYTHONDONTWRITEBYTECODE': '1', 'TMPDIR': tmpdir})
    out = os.path.join(rel['path'], 'tests.log')
    cmd = ['nice', '-n', '15', cfg['python'], '-m', 'pytest', '-q', '-p', 'no:cacheprovider', '-rs',
           '--basetemp', basetemp] + paths
    log('running the suite from {} ({}); log {}; temporary files in {}'.format(rel['src'], ' '.join(paths), out, tmp))
    try:
        with open(out, 'w') as f:
            rc = subprocess.call(cmd, cwd=rel['src'], env=env, stdout=f, stderr=subprocess.STDOUT)
    finally:
        if not remove_tree(tmp):
            log('warning: could not remove {}: delete it by hand'.format(tmp))
    with open(out) as f:
        lines = [line.strip() for line in f if line.strip()]
    # pytest's last line is the count ("308 passed in 128.57s"); the lines before it are progress dots
    summary = 'exit {}: {}'.format(rc, lines[-1].strip('= ') if lines else 'no output')
    with open(os.path.join(rel['path'], TESTS_JSON), 'w') as f:
        json.dump({'exit_code': rc, 'at': iso(utcnow()), 'summary': summary, 'paths': paths}, f)
    log('suite ' + summary)
    return rc


# --------------------------------------------------------------------------- activating

def activate(cfg, rid, store=None, source='cli', log=print):
    rel = read_release(cfg, rid)
    if not rel:
        raise ReleaseError('no release {} (or it has no {} or src/)'.format(rid, RELEASE_JSON))
    root = releases_dir(cfg)
    before = current_id(cfg)
    tmp = os.path.join(root, '.{}.{}'.format(CURRENT, os.getpid()))
    if os.path.lexists(tmp):
        os.unlink(tmp)
    os.symlink(rid, tmp)
    os.replace(tmp, os.path.join(root, CURRENT))
    with open(os.path.join(root, HISTORY), 'a') as f:
        f.write(json.dumps({'ts': iso(utcnow()), 'from': before, 'to': rid, 'source': source}) + '\n')
    if store is not None:
        store.event('release', source=source, message='{} -> {}'.format(before or 'live tree', rid),
                    data={'from': before, 'to': rid})
        store.kv_set(SYNC_KEY, rid)
    log('current release: {} (was {}); jobs started from now on run it'.format(rid, before or 'the live tree'))
    return before


# --------------------------------------------------------------------------- the live tree

def _environ(pid, proc):
    try:
        with open(os.path.join(proc, str(pid), 'environ'), 'rb') as f:
            raw = f.read()
    except OSError:
        return None
    env = {}
    for item in raw.split(b'\0'):
        k, sep, v = item.partition(b'=')
        if sep:
            env[k.decode('utf-8', 'replace')] = v.decode('utf-8', 'replace')
    return env


def _is_work(argv):
    if any(a.endswith(m) for a in argv for m in WORK_MARKERS):
        return True
    return any(argv[k:k + 2] == ['jobqueue', 'jobs'] for k in range(len(argv)))


def live_tree_users(cfg, proc='/proc', exclude=()):
    """Pipeline, download and job processes that read code from the live tree: {pid: command}.

    A process whose ORCHARD_PATH is not the live tree is pinned (a release, or someone's scratch copy). One whose
    environment cannot be read is counted, to be safe.
    """
    users = {}
    if not os.path.isdir(proc):
        return users
    live = os.path.normpath(cfg['src_dir'])
    skip = set(exclude) | {os.getpid()}
    for name in os.listdir(proc):
        if not name.isdigit() or int(name) in skip:
            continue
        try:
            with open(os.path.join(proc, name, 'cmdline'), 'rb') as f:
                argv = [a.decode('utf-8', 'replace') for a in f.read().split(b'\0') if a]
        except OSError:
            continue
        if not argv or not _is_work(argv):
            continue
        env = _environ(name, proc)
        orchard = (env or {}).get('ORCHARD_PATH')
        if env is not None and orchard and os.path.normpath(orchard) != live:
            continue
        users[int(name)] = ' '.join(argv)[:160]
    return users


def live_diff(cfg, rid):
    """Files of the release missing from, or different in, the live tree (the live tree's extra files don't count).

    Credentials never count: the live tree's are the ones in use, and a release built before they were linked holds
    a copy that may be a revoked password by now.
    """
    rel = read_release(cfg, rid)
    if not rel:
        raise ReleaseError('no release {}'.format(rid))
    live = os.path.normpath(cfg['src_dir'])
    links = cfg.get('release_link') or ()
    out = []
    for relpath, path in sorted(tree_files(rel['src'], cfg.get('release_exclude') or ()).items()):
        if _linked(relpath, links):
            continue
        if not same_bytes(path, os.path.join(live, relpath)):
            out.append(relpath)
    return out


def sync_live(cfg, rid, dry_run=False, force=False, proc='/proc', now=None, log=print):
    """Bring the live tree up to a release: write the files that differ, back up what they replace, delete nothing.

    Returns {'blocked': {pid: cmd}} if a pipeline or download process still runs from the live tree (unless
    force), else {'written': [...], 'backup': dir or None}.
    """
    rel = read_release(cfg, rid)
    if not rel:
        raise ReleaseError('no release {}'.format(rid))
    users = {} if force else live_tree_users(cfg, proc)
    if users:
        return {'blocked': users}
    todo = live_diff(cfg, rid)
    if dry_run or not todo:
        return {'written': todo, 'backup': None, 'dry_run': dry_run}
    now = now or utcnow()
    live = os.path.normpath(cfg['src_dir'])
    backup = os.path.join(releases_dir(cfg), 'live_backups', '{}-before-{}'.format(now.strftime('%Y%m%dT%H%M%SZ'), rid))
    for relpath in todo:
        src = os.path.join(rel['src'], relpath)
        dst = os.path.join(live, relpath)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.lexists(dst):
            b = os.path.join(backup, relpath)
            os.makedirs(os.path.dirname(b), exist_ok=True)
            shutil.copy2(dst, b)
            mode = stat.S_IMODE(os.stat(dst).st_mode)
        else:
            mode = stat.S_IMODE(os.stat(src).st_mode) | 0o644
        tmp = os.path.join(os.path.dirname(dst), '.{}.{}.release-tmp'.format(os.path.basename(dst), os.getpid()))
        shutil.copyfile(src, tmp)
        os.chmod(tmp, mode)
        os.replace(tmp, dst)
    log('live tree {}: {} files written{}'.format(live, len(todo),
                                                   '; replaced files backed up in ' + backup if os.path.isdir(backup)
                                                   else ''))
    return {'written': todo, 'backup': backup if os.path.isdir(backup) else None}


def sync_pending(cfg, store, proc='/proc', log=print):
    """Do the pending live-tree sync if nothing runs from the live tree any more. The dispatcher calls this
    every poll; `q release deploy` calls it once straight away.

    Returns None if no sync is pending, {'blocked': {pid: cmd}} while it has to wait, else sync_live's result
    with 'release' added. The live tree is brought up to the current release, whichever that is by then.
    """
    if not store.kv_get(SYNC_KEY):
        return None
    rid = current_id(cfg)
    if not rid or not read_release(cfg, rid):
        store.kv_delete(SYNC_KEY)
        return None
    res = sync_live(cfg, rid, proc=proc, log=log)
    if 'blocked' in res:
        return res
    store.kv_delete(SYNC_KEY)
    res['release'] = rid
    if res['written']:
        store.event('release-sync', source='release', message='{}: {} files into the live tree'.format(
            rid, len(res['written'])), data={'release': rid, 'written': res['written'], 'backup': res['backup']})
    return res


def queue_code_changed(res):
    """True if a sync rewrote the queue's own code, so the running dispatcher is out of date."""
    return any(path.startswith('jobqueue' + os.sep) for path in (res or {}).get('written') or ())
