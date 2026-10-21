#!/usr/bin/env python3
"""
rawpack.py
==========

Let the pipeline run on a raw night whose frames are stored compressed (``.fz``), and store them compressed again
afterwards.  The pipeline only ever reads plain frames (decision D2 of the frame compression inventory, 16 October
2026); ZLP_pipeline.sh calls this before T1 and after the night's last stage.

    python utils/rawpack.py unpack <Observations/TEL/images/NIGHT> [--workers N]
    python utils/rawpack.py repack <Observations/TEL/images/NIGHT> [--workers N]
    python utils/rawpack.py status <Observations/TEL/images/NIGHT>

``unpack``
    Every ``X.fz`` in the night, at the top level and one directory down (where T1 lists frames), is decompressed to
    ``X``: the downloaded file byte for byte for a byte-exact ``.fz`` (``framepack rawpack``), otherwise the same
    pixels and header cards.  Each working copy is verified, given the original's mtime (so T7's WCS cache, keyed on
    raw name, size and mtime, still matches), and the ``.fz`` is moved into a hidden ``.rawpack/`` beside it, where
    no pipeline glob looks.  ``.rawpack/unpacked.json`` records each working copy's inode, size, mtime and ctime.

``repack``
    Each ``.rawpack/X.fz`` goes back beside ``X`` and ``X`` is deleted, if ``X`` is unchanged since ``unpack``: the
    same inode, size, mtime and ctime, or else the same bytes as the ``.fz`` rebuilds.  Nothing in the pipeline
    writes into raw frames since PR #64; a working copy that changed anyway is compressed again (byte-exact, as the
    new original) and the earlier ``.fz`` is kept as ``.rawpack/X.fz.<UTC stamp>``, which nothing here deletes.

Both are safe to repeat after a run that was killed half-way: ``unpack`` takes over the working copies a killed
run left, and ``repack`` finishes a repack that was cut short.  A night left unpacked is consistent for every
reader (plain frames only), it just takes twice the space until the next run on it.

The caller holds the night (the queue's night lock): nothing else may unpack or repack the same night meanwhile.
"""
import argparse
import datetime
import json
import os
import stat
import sys
from multiprocessing import Pool

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import framepack as fp  # noqa: E402

STATE = '.rawpack'              # hidden: bash and Python globs skip it, so T1 never lists the .fz files
RECORD = 'unpacked.json'
PLAIN = fp.FITS_SUFFIXES


def frame_dirs(night):
    """The night directory and the directories one level down, as T1's get_liste reads them."""
    out = [night]
    for name in sorted(os.listdir(night)):
        p = os.path.join(night, name)
        if not name.startswith('.') and os.path.isdir(p) and not os.path.islink(p):
            out.append(p)
    return out


def _is_frame_fz(name):
    low = name.lower()
    return low.endswith(fp.COMPRESSED) and low[:-len(fp.COMPRESSED)].endswith(PLAIN) and not name.startswith('.')


def _signature(path):
    s = os.lstat(path)
    return {'ino': s.st_ino, 'size': s.st_size, 'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns}


def _fsync_dir(d):
    fd = os.open(d, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _load_record(d):
    try:
        with open(os.path.join(d, STATE, RECORD)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _save_record(d, record):
    path = os.path.join(d, STATE, RECORD)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(record, f, indent=1, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def _same_bytes(plain, fz):
    """Whether the plain file holds exactly what the .fz rebuilds."""
    with open(plain, 'rb') as f:
        return f.read() == fp.plain_bytes(fz)


# --------------------------------------------------------------------------- unpack

def _working_copy(fz, plain):
    """Decompress ``fz`` to ``plain`` through a verified hidden temporary file."""
    d = os.path.dirname(plain)
    tmp = os.path.join(d, '.' + os.path.basename(plain) + '.rawpack-{}'.format(os.getpid()))
    try:
        fp.decompress(fz, tmp, record=False)
        info = fp.original_info(fz)
        if info is not None:
            fp.check_rebuilt(fz, tmp)
            mtime = info['ORIGMTIM']
        else:
            fp.check_lossless(fz, tmp)
            mtime = os.stat(fz).st_mtime_ns
        st = os.stat(fz)
        os.chmod(tmp, stat.S_IMODE(st.st_mode))
        os.utime(tmp, ns=(st.st_atime_ns, mtime))
        fd = os.open(tmp, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        os.replace(tmp, plain)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _unpack_one(args):
    """One frame: returns (dir, plain name, ready, signature to record or None, message or '')."""
    d, name, where = args
    plain_name = name[:-len(fp.COMPRESSED)]
    plain = os.path.join(d, plain_name)
    kept = os.path.join(d, STATE, name)
    try:
        if where == 'kept':
            # an earlier run unpacked this frame and never repacked it (killed): keep its working copy and the
            # signature recorded then, so repack still sees anything that run changed; remake only a lost copy
            if os.path.exists(plain):
                return d, plain_name, True, None, ''
            _working_copy(kept, plain)
            return d, plain_name, True, _signature(plain), ''
        fz = os.path.join(d, name)
        if os.path.exists(plain):
            # both visible: a repack cut short, or a frame stored twice; only the first is safe to take over
            if not _same_bytes(plain, fz):
                return d, plain_name, False, None, 'both {} and its .fz exist and differ: left as they are'.format(
                    plain_name)
        else:
            _working_copy(fz, plain)
        os.replace(fz, kept)
        return d, plain_name, True, _signature(plain), ''
    except Exception as e:                      # noqa: BLE001 - one bad frame must not stop the night
        return d, plain_name, False, None, '{}: {}: {}'.format(name, type(e).__name__, e)


def unpack(night, workers=1):
    """Decompress every .fz frame of a raw night for a pipeline run; returns (frames unpacked, failure messages)."""
    jobs, records = [], {}
    for d in frame_dirs(night):
        state = os.path.join(d, STATE)
        names = sorted(n for n in os.listdir(d) if _is_frame_fz(n))
        kept = sorted(n for n in os.listdir(state) if _is_frame_fz(n)) if os.path.isdir(state) else []
        if not names and not kept:
            continue
        os.makedirs(state, exist_ok=True)
        records[d] = _load_record(d)
        jobs += [(d, n, 'visible') for n in names] + [(d, n, 'kept') for n in kept if n not in names]
    results = _map(_unpack_one, jobs, workers)
    failures = []
    for d, plain_name, ready, sig, msg in results:
        if sig is not None:
            records[d][plain_name] = sig
        if msg:
            failures.append(os.path.join(os.path.relpath(d, night), msg) if d != night else msg)
    for d, record in records.items():
        _save_record(d, record)
        _fsync_dir(d)
        _fsync_dir(os.path.join(d, STATE))
    return sum(1 for r in results if r[2]), failures


# --------------------------------------------------------------------------- repack

def _stamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')


def _repack_one(args):
    """One frame: returns (dir, plain name, outcome, message)."""
    d, name, sig = args
    plain_name = name[:-len(fp.COMPRESSED)]
    plain, fz, kept = os.path.join(d, plain_name), os.path.join(d, name), os.path.join(d, STATE, name)
    try:
        if not os.path.exists(plain):
            os.replace(kept, fz)
            return d, plain_name, 'restored', ''
        if sig == _signature(plain) or _same_bytes(plain, kept):
            os.replace(kept, fz)                # the .fz first: the frame is never without a stored copy
            os.remove(plain)
            return d, plain_name, 'repacked', ''
        # the run changed the raw frame: store what it wrote, keep what was there
        tmp = os.path.join(d, '.' + name + '.rawpack-{}'.format(os.getpid()))
        try:
            fp.compress_exact(plain, tmp)
            fp.check_exact(plain, tmp)
            st = os.stat(plain)
            os.chmod(tmp, stat.S_IMODE(st.st_mode))
            os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
            history = '{}.{}'.format(kept, _stamp())
            os.replace(kept, history)
            os.replace(tmp, fz)
            os.remove(plain)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
        return d, plain_name, 'recompressed', '{} was changed by the run: compressed again, the earlier .fz kept ' \
                                              'as {}'.format(plain_name, os.path.relpath(history, d))
    except Exception as e:                      # noqa: BLE001
        return d, plain_name, 'failed', '{}: {}: {}'.format(plain_name, type(e).__name__, e)


def repack(night, workers=1):
    """Store an unpacked night compressed again; returns ({outcome: count}, messages)."""
    jobs, dirs = [], []
    for d in frame_dirs(night):
        state = os.path.join(d, STATE)
        if not os.path.isdir(state):
            continue
        dirs.append(d)
        record = _load_record(d)
        jobs += [(d, n, record.get(n[:-len(fp.COMPRESSED)])) for n in sorted(os.listdir(state)) if _is_frame_fz(n)]
    results = _map(_repack_one, jobs, workers)
    counts, messages = {}, []
    for d, plain_name, outcome, msg in results:
        counts[outcome] = counts.get(outcome, 0) + 1
        if msg:
            messages.append(msg)
    for d in dirs:
        state = os.path.join(d, STATE)
        if not any(_is_frame_fz(n) for n in os.listdir(state)):
            try:
                os.remove(os.path.join(state, RECORD))
            except FileNotFoundError:
                pass
            if not os.listdir(state):
                os.rmdir(state)
        _fsync_dir(d)
    return counts, messages


def status(night):
    out = {'compressed': 0, 'unpacked': 0, 'plain': 0}
    for d in frame_dirs(night):
        for n in os.listdir(d):
            if _is_frame_fz(n):
                out['compressed'] += 1
            elif n.lower().endswith(PLAIN) and not n.startswith('.'):
                out['plain'] += 1
        state = os.path.join(d, STATE)
        if os.path.isdir(state):
            out['unpacked'] += sum(1 for n in os.listdir(state) if _is_frame_fz(n))
    return out


def _map(fn, jobs, workers):
    if workers > 1 and len(jobs) > 1:
        with Pool(min(workers, len(jobs))) as pool:
            return list(pool.imap_unordered(fn, jobs, chunksize=8))
    return [fn(j) for j in jobs]


def main(argv=None):
    ap = argparse.ArgumentParser(description='Decompress a raw night for a pipeline run, or store it compressed again.')
    ap.add_argument('op', choices=('unpack', 'repack', 'status'))
    ap.add_argument('night', help='Observations/<TEL>/images/<NIGHT>')
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args(argv)
    if not os.path.isdir(a.night):
        print('rawpack: no such night directory {}'.format(a.night))
        return 2
    if a.op == 'status':
        s = status(a.night)
        print('rawpack status {}: {} compressed, {} unpacked for a run, {} plain'.format(
            a.night, s['compressed'], s['unpacked'], s['plain']))
        return 0
    if a.op == 'unpack':
        n, failures = unpack(a.night, a.workers)
        for m in failures:
            print('rawpack unpack FAILED: ' + m)
        print('rawpack unpack {}: {} frames decompressed for the run, {} failed'.format(a.night, n, len(failures)))
        return 1 if failures else 0
    counts, messages = repack(a.night, a.workers)
    for m in messages:
        print('rawpack repack: ' + m)
    print('rawpack repack {}: '.format(a.night) + (', '.join('{} {}'.format(v, k) for k, v in sorted(counts.items()))
                                                   or 'nothing unpacked'))
    return 1 if counts.get('failed') else 0


if __name__ == '__main__':
    sys.exit(main())
