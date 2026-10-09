"""Add-only download of ESO frames into a night directory. The look-back's fetch jobs, `q add --download` and
the ESO backfill all come through here.

SSO_download.py must not be looped over old nights: for Callisto (all dates) and Ganymede (after
2025-03-03) its transformation_check() deletes every existing frame without ASTRAROT/ASTRAMIR before
re-downloading it, which destroys pre-Astra Callisto frames; it counts only .fits/.fts files as present, so it
fetches the .fits.fz nights and the unpacked SPIRIT nights all over again; it also creates empty night
directories and rewrites download_log.csv. This path only ever adds:

  - each file is downloaded under a temporary name into a private staging directory on the archive's
    filesystem, at most `workers` files at a time, with exponential backoff on errors; gzip and compress
    data are decompressed as request_eso.py does;
  - it is verified before anything else happens (nights.verify_frame: SIMPLE = T, a size that is a multiple of
    2,880, sane BITPIX and NAXISn, every HDU accounted for by the size); a file that fails is fetched again;
  - datacubes are unpacked there (download.unpack_datacubes), and the Astra transform is applied there, only to
    frames whose OBSERVER is Astra and that lack ASTRAROT/ASTRAMIR; the results are verified again;
  - each frame is hard-linked into the night directory only if no file with its stem exists there (any FITS
    suffix) and no frame stored under another name (the 2018 ACP names) has its DATE-OBS; the night
    directory is created only when there is a frame to put in it;
  - nothing outside staging is deleted, and download_log.csv is not touched;
  - every file added is appended to a manifest (time, telescope, night, file, ESO row, bytes, md5, job), so
    the additions can be audited or undone.

The one exception to adding is a swap, and only on request: replace='invalid' (repair) swaps an existing file
of the same name that fails verification, e.g. a truncated frame (decision G3); replace='all' (refetch) swaps one
that differs from the fresh copy. The old file is first hard-linked into a quarantine directory, then the
verified fresh copy is renamed over it, so the swap is atomic and the old file is kept.
"""
import concurrent.futures
import csv
import errno
import fcntl
import gzip
import hashlib
import os
import shutil
import subprocess
import threading
import time

from .nights import FITS_SUFFIXES, FitsError, ForeignMatcher, read_header, strip_suffix, verify_frame
from .util import iso, utcnow

DATAPORTAL = 'https://dataportal.eso.org/dataPortal/file/{}'
REPLACE_MODES = ('never', 'invalid', 'all')
MANIFEST_FIELDS = ('time', 'telescope', 'night', 'action', 'file', 'dp_id', 'dp_type', 'object', 'bytes', 'md5',
                   'job', 'old_bytes', 'old_md5', 'quarantine')


# --------------------------------------------------------------------------- downloading


class FetchError(RuntimeError):
    def __init__(self, message, retry=True):
        super().__init__(message)
        self.retry = retry          # worth another attempt in this job


def _requests_get(url, headers, timeout):
    import requests
    return requests.get(url, headers=headers, stream=True, timeout=timeout)


class EsoFileFetcher:
    """download(dp_ids, dest) -> paths of verified files written in dest, at most `workers` at once.

    token(force=False) returns a bearer token (force: get a new one); get(url, headers, timeout) returns a
    streaming requests response. Each file gets `attempts` tries with exponential backoff; after
    `breaker` consecutive network failures nothing more is tried in this job. Afterwards self.failed maps each
    dp_id that did not land to its reason, and self.retry_later holds those that may work in a later run
    (network trouble, as opposed to a file ESO does not have or keeps serving broken).
    """

    def __init__(self, token, get=None, workers=3, attempts=4, backoff=30.0, timeout=(30, 600), breaker=8,
                 log=print, sleep=time.sleep):
        self.token = token
        self.get = get or _requests_get
        self.workers = max(1, int(workers))
        self.attempts = max(1, int(attempts))
        self.backoff = float(backoff)
        self.timeout = timeout
        self.breaker = int(breaker)
        self.log = log
        self.sleep = sleep
        self.failed = {}
        self.retry_later = set()
        self.bytes = 0
        self.done = 0
        self._lock = threading.Lock()
        self._consecutive = 0
        self._tripped = False

    def __call__(self, dp_ids, dest):
        os.makedirs(dest, exist_ok=True)
        out = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as pool:
            futures = {pool.submit(self._one, d, dest): d for d in dp_ids}
            for fut in concurrent.futures.as_completed(futures):
                dp_id = futures[fut]
                try:
                    path = fut.result()
                except FetchError as e:
                    with self._lock:
                        self.failed[dp_id] = str(e)
                        if e.retry:
                            self.retry_later.add(dp_id)
                    self.log('fetch: {} not downloaded: {}'.format(dp_id, e))
                    continue
                out.append(path)
        return out

    def _one(self, dp_id, dest):
        last = None
        for attempt in range(self.attempts):
            if self._tripped:
                raise FetchError('not tried: {} downloads in a row failed'.format(self.breaker))
            if attempt:
                self.sleep(self.backoff * 2 ** (attempt - 1))
            try:
                path = self._fetch(dp_id, dest)
            except FetchError as e:
                last = e
                with self._lock:
                    if e.retry:
                        self._consecutive += 1
                        if self._consecutive >= self.breaker:
                            self._tripped = True
                if not e.retry:
                    break
                continue
            with self._lock:
                self._consecutive = 0
                self.done += 1
                self.bytes += os.path.getsize(path)
            return path
        if last.retry and str(last).startswith('not valid FITS'):
            last.retry = False          # ESO keeps serving a broken file: another run will not help
        raise last

    def _fetch(self, dp_id, dest):
        part = os.path.join(dest, '.{}.part'.format(dp_id))
        tmp = os.path.join(dest, '.{}.fits.tmp'.format(dp_id))
        try:
            try:
                headers = {'Authorization': 'Bearer {}'.format(self._token())}
                resp = self.get(DATAPORTAL.format(dp_id), headers, self.timeout)
            except FetchError:
                raise
            except Exception as e:  # requests' connection and timeout errors
                raise FetchError('{}: {}'.format(type(e).__name__, e))
            try:
                self._save(resp, part)
            finally:
                close = getattr(resp, 'close', None)
                if close:
                    close()
            with open(part, 'rb') as f:
                magic = f.read(2)
            if magic == b'\x1f\x8b':                       # gzip, as request_eso.py handles it
                with gzip.open(part, 'rb') as src, open(tmp, 'wb') as out:
                    shutil.copyfileobj(src, out, 1 << 20)
            elif magic == b'\x1f\x9d':                     # compress (.Z): request_eso.py uses uncompress too
                with open(tmp, 'wb') as out:
                    if subprocess.call(['uncompress', '-c', part], stdout=out, stderr=subprocess.DEVNULL):
                        raise FetchError('not valid FITS: uncompress failed')
            else:
                os.replace(part, tmp)
            try:
                verify_frame(tmp)
            except FitsError as e:
                raise FetchError('not valid FITS: {}'.format(e))
            path = os.path.join(dest, dp_id + '.fits')
            os.replace(tmp, path)
            return path
        except OSError as e:
            raise FetchError('{}: {}'.format(type(e).__name__, e))
        finally:
            for p in (part, tmp):
                if os.path.exists(p):
                    os.unlink(p)

    def _token(self, force=False):
        try:
            return self.token(force=True) if force else self.token()
        except Exception as e:  # ESO's login service: try again later
            raise FetchError('no ESO token: {}'.format(e))

    def _save(self, resp, part):
        status = resp.status_code
        if status in (401, 403):
            self._token(force=True)
            raise FetchError('HTTP {} (token renewed)'.format(status))
        if status == 404:
            raise FetchError('HTTP 404: not in the archive', retry=False)
        if status != 200:
            raise FetchError('HTTP {}'.format(status), retry=status == 429 or status >= 500)
        if 'text/html' in (resp.headers.get('Content-Type') or ''):
            self._token(force=True)
            raise FetchError('HTML instead of a file (token renewed)')
        expected = resp.headers.get('Content-Length')
        n = 0
        with open(part, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
                    n += len(chunk)
        if expected is not None and str(expected).isdigit() and int(expected) != n:
            raise FetchError('connection closed after {} of {} bytes'.format(n, expected))


# --------------------------------------------------------------------------- placing frames


def md5sum(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


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


def _link_or_copy(src, dst):
    try:
        os.link(src, dst)
    except OSError as e:
        if e.errno != errno.EXDEV:
            raise
        shutil.copy2(src, dst)


def _swap(src, dst, quarantine_dir):
    """Keep dst in quarantine_dir, then atomically rename a link to src over it. Returns the kept copy's path."""
    os.makedirs(quarantine_dir, exist_ok=True)
    kept = os.path.join(quarantine_dir, os.path.basename(dst))
    if os.path.lexists(kept):
        kept = '{}.{}'.format(kept, time.strftime('%Y%m%dT%H%M%S'))
    _link_or_copy(dst, kept)
    tmp = os.path.join(os.path.dirname(dst), '.{}.jobqueue-new{}'.format(os.path.basename(dst), os.getpid()))
    _link_or_copy(src, tmp)
    try:
        os.rename(tmp, dst)
    except OSError:
        os.unlink(tmp)
        raise
    return kept


def existing_frame(night_path, stem):
    """The file in night_path with this stem and any FITS suffix, or None."""
    for suffix in FITS_SUFFIXES:
        p = os.path.join(night_path, stem + suffix)
        if os.path.lexists(p):
            return p
    return None


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


def _date_obs(path):
    try:
        return str(read_header(path, max_blocks=24).get('DATE-OBS', '')).strip()
    except OSError:
        return ''


class Manifest:
    """Append-only CSV of every file the fetch added or swapped, shared by all jobs (flock per write)."""

    def __init__(self, path, **fixed):
        self.path = path
        self.fixed = fixed

    def __call__(self, entry):
        row = dict(self.fixed, time=iso(utcnow()), **entry)
        os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
        with open(self.path, 'a', newline='') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                f.seek(0, 2)
                w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS, extrasaction='ignore')
                if f.tell() == 0:
                    w.writeheader()
                w.writerow({k: row.get(k, '') for k in MANIFEST_FIELDS})
                f.flush()
                os.fsync(f.fileno())
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)


def add_only_fetch(rows, night_path, staging, download, unpack=None, is_cube=None, transform=None, log=print,
                   batch=25, verify=None, replace='never', quarantine=None, foreign=None, record=None):
    """Fetch ESO rows (dicts with dp_id) into night_path without touching anything already there.

    download(dp_ids, dest) -> paths written in dest (EsoFileFetcher, or anything with that signature);
    is_cube(path) -> bool; unpack(path, dest) -> frame paths; transform(stems, dirname) applies the Astra
    transform in place (download.astra_transform.astra_transform); verify(path) raises FitsError for a broken
    file; foreign: {DATE-OBS: name} of frames stored under other names, which are never duplicated;
    record(entry) is called for every file added or swapped; replace and quarantine: see the module docstring.
    Works through the rows `batch` at a time, so staging never holds more than one batch.
    """
    if replace not in REPLACE_MODES:
        raise ValueError('replace must be one of {}'.format(', '.join(REPLACE_MODES)))
    if replace != 'never' and not quarantine:
        raise ValueError('replace={} needs a quarantine directory'.format(replace))
    os.makedirs(staging, exist_ok=True)
    unpacked_dir = os.path.join(staging, 'frames')
    match = ForeignMatcher(foreign) if foreign else None
    res = {'requested': len(rows), 'downloaded': 0, 'frames': 0, 'added': 0, 'added_science': 0, 'replaced': 0,
           'already_present': 0, 'foreign_duplicates': 0, 'bytes_added': 0, 'not_downloaded': [],
           'unpack_failed': [], 'invalid': [], 'settled': []}
    try:
        for k in range(0, len(rows), batch):
            chunk = rows[k:k + batch]
            os.makedirs(unpacked_dir, exist_ok=True)
            got = download([r['dp_id'] for r in chunk], staging) or []
            res['downloaded'] += len(got)
            landed = {strip_suffix(os.path.basename(p)): p for p in got}
            frames_by_row = {}
            for r in chunk:
                dp_id = r['dp_id']
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
            res['frames'] += len(frames)
            if transform:
                astra = [p for p in frames if _needs_astra(p)]
                if astra:
                    transform([strip_suffix(os.path.basename(p)) for p in astra], unpacked_dir)
            by_id = {r['dp_id']: r for r in chunk}
            for dp_id, paths in frames_by_row.items():
                placed, invalid = _place_row(by_id[dp_id], paths, night_path, res, verify, replace, quarantine, match,
                                             record, log)
                if placed == 0 and invalid == 0:
                    res['settled'].append(dp_id)  # downloaded, nothing new: do not ask for it again
            shutil.rmtree(unpacked_dir, ignore_errors=True)
            for name in os.listdir(staging):  # anything a failed unpack left behind
                p = os.path.join(staging, name)
                if os.path.isfile(p):
                    os.unlink(p)
            log('fetch: {} of {} rows done, {} frames added, {} replaced'.format(
                min(k + batch, len(rows)), len(rows), res['added'], res['replaced']))
        log('fetch: {requested} rows requested, {downloaded} downloaded, {frames} frames, {added} added '
            '({added_science} science, {gb:.2f} GB), {replaced} replaced, {already_present} already present, '
            '{foreign_duplicates} already present under another name, {n} not downloaded, {u} cubes did not '
            'unpack, {i} frames failed verification'.format(
                n=len(res['not_downloaded']), u=len(res['unpack_failed']), i=len(res['invalid']),
                gb=res['bytes_added'] / 1e9, **res))
        return res
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _place_row(row, paths, night_path, res, verify, replace, quarantine, match, record, log):
    """Put one row's frames into the night directory. Returns (frames added or swapped, frames that failed
    verification)."""
    placed = invalid = 0
    for p in paths:
        if not os.path.exists(p):
            continue
        name = os.path.basename(p)
        try:
            if verify:
                verify(p)
        except FitsError as e:
            res['invalid'].append(name)
            invalid += 1
            log('fetch: {} not placed: {}'.format(name, e))
            continue
        stem = strip_suffix(name)
        dst = os.path.join(night_path, name)
        old = existing_frame(night_path, stem)
        if old is None:
            if match is not None and match.find(_date_obs(p)) is not None:
                res['foreign_duplicates'] += 1
                continue
            os.makedirs(night_path, exist_ok=True)
            if not _link_new(p, dst):
                res['already_present'] += 1
                continue
            action, kept, old_size, old_md5 = 'added', '', '', ''
            res['added'] += 1
        elif replace == 'never' or old != dst:
            res['already_present'] += 1  # never across formats (.fits.fz stays as it is)
            continue
        else:
            old_size = os.path.getsize(old)
            old_md5 = ''
            if replace == 'invalid':
                try:
                    verify_frame(old)
                    res['already_present'] += 1
                    continue
                except FitsError as e:
                    log('fetch: replacing {}: {}'.format(name, e))
            else:
                old_md5 = md5sum(old)
                if old_md5 == md5sum(p):
                    res['already_present'] += 1
                    continue
            kept = _swap(p, dst, quarantine)
            action = 'replaced'
            res['replaced'] += 1
        placed += 1
        size = os.path.getsize(dst)
        res['bytes_added'] += size
        science = _is_science(dst)
        if science and action == 'added':
            res['added_science'] += 1
        if record:
            record({'action': action, 'file': name, 'dp_id': row.get('dp_id', ''), 'dp_type': row.get('dp_type', ''),
                    'object': (row.get('object') or '').strip(), 'bytes': size, 'md5': md5sum(dst),
                    'old_bytes': old_size, 'old_md5': old_md5, 'quarantine': kept})
    return placed, invalid

