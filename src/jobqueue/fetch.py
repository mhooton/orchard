"""Add-only download of missing ESO frames into a night directory.

SSO_download.py must not be looped over old nights: for Callisto (all dates) and Ganymede (after
2025-03-03) its transformation_check() deletes every existing frame without ASTRAROT/ASTRAMIR before
re-downloading it, which destroys pre-Astra Callisto frames; it also creates empty night directories and
rewrites download_log.csv. This path only ever adds:

  - frames are downloaded into a private staging directory (same filesystem as the archive);
  - datacubes are unpacked there, and the Astra transform is applied there, only to frames whose
    OBSERVER is Astra and that do not already carry ASTRAROT/ASTRAMIR;
  - each frame is hard-linked into the night directory only if no file of that name exists, so nothing
    on disk is ever overwritten; the night directory is created only when there is a frame to put in it;
  - nothing outside the staging directory is deleted, and download_log.csv is not touched.
"""
import errno
import os
import shutil

from .nights import read_header, strip_suffix


def _link_new(src, dst):
    """Hard-link src to dst unless dst exists. Returns True if added."""
    if os.path.lexists(dst):
        return False
    try:
        os.link(src, dst)
        return True
    except FileExistsError:
        return False
    except OSError as e:
        if e.errno != errno.EXDEV:
            raise
    tmp = os.path.join(os.path.dirname(dst), '.{}.jobqueue-tmp{}'.format(os.path.basename(dst), os.getpid()))
    shutil.copyfile(src, tmp)
    try:
        os.link(tmp, dst)
        return True
    except FileExistsError:
        return False
    finally:
        os.unlink(tmp)


def _needs_astra(path):
    try:
        h = read_header(path, max_blocks=24)
    except OSError:
        return False
    return str(h.get('OBSERVER', '')).strip() == 'Astra' and 'ASTRAROT' not in h and 'ASTRAMIR' not in h


def _is_science(path):
    try:
        return str(read_header(path, max_blocks=24).get('IMAGETYP', '')).lower().startswith('light')
    except OSError:
        return False


def add_only_fetch(rows, night_path, staging, download, unpack=None, is_cube=None, transform=None, log=print,
                   batch=25):
    """Fetch ESO rows (dicts with dp_id) into night_path without touching anything already there.

    download(dp_ids, dest) -> paths written in dest; is_cube(path) -> bool; unpack(path, dest) -> frame paths;
    transform(stems, dirname) applies the Astra transform in place (download.astra_transform.astra_transform).
    """
    os.makedirs(staging, exist_ok=True)
    unpacked_dir = os.path.join(staging, 'frames')
    os.makedirs(unpacked_dir, exist_ok=True)
    res = {'requested': len(rows), 'downloaded': 0, 'frames': 0, 'added': 0, 'added_science': 0,
           'already_present': 0, 'not_downloaded': [], 'unpack_failed': [], 'settled': []}
    try:
        frames_by_row = {}
        for k in range(0, len(rows), batch):
            chunk = [r['dp_id'] for r in rows[k:k + batch]]
            got = download(chunk, staging) or []
            res['downloaded'] += len(got)
            landed = {strip_suffix(os.path.basename(p)): p for p in got}
            for dp_id in chunk:
                path = landed.get(dp_id)
                if path is None or not os.path.exists(path):
                    res['not_downloaded'].append(dp_id)
                    continue
                if is_cube and unpack and is_cube(path):
                    frames = list(unpack(path, unpacked_dir) or [])
                    if not frames:
                        res['unpack_failed'].append(dp_id)  # not settled: ask again another day
                        continue
                    frames_by_row[dp_id] = frames
                else:
                    dst = os.path.join(unpacked_dir, os.path.basename(path))
                    os.replace(path, dst)
                    frames_by_row[dp_id] = [dst]
        frames = [p for ps in frames_by_row.values() for p in ps if os.path.exists(p)]
        res['frames'] = len(frames)
        if transform:
            astra = [p for p in frames if _needs_astra(p)]
            if astra:
                transform([strip_suffix(os.path.basename(p)) for p in astra], unpacked_dir)
        for dp_id, paths in frames_by_row.items():
            added_here = 0
            for p in paths:
                if not os.path.exists(p):
                    continue
                dst = os.path.join(night_path, os.path.basename(p))
                if os.path.lexists(dst):
                    res['already_present'] += 1
                    continue
                os.makedirs(night_path, exist_ok=True)
                if _link_new(p, dst):
                    res['added'] += 1
                    added_here += 1
                    if _is_science(dst):
                        res['added_science'] += 1
                else:
                    res['already_present'] += 1
            if added_here == 0:
                res['settled'].append(dp_id)  # downloaded, nothing new: do not ask for it again
        log('fetch: {requested} rows requested, {downloaded} downloaded, {frames} frames, {added} added '
            '({added_science} science), {already_present} already present, {n} not downloaded, {u} cubes did not '
            'unpack'.format(n=len(res['not_downloaded']), u=len(res['unpack_failed']), **res))
        return res
    finally:
        shutil.rmtree(staging, ignore_errors=True)
