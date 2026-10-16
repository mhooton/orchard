#!/usr/bin/env python3
"""
framepack.py
============

Store SPECULOOS frames compactly, and restore them, without losing data.

Three conversions, decided on 16 October 2026 (the frame compression
inventory):

* ``int16``       a calibrated ``proc*.fits`` frame, float32 -> 16-bit
                  integers.  Pixels are rounded to whole ADU, half away from
                  zero as CASU's NINT does, so no pixel moves by more than
                  0.5 ADU apart from the few clipped ones counted below.
* ``compress``    any 16-bit frame, calibrated or raw -> lossless RICE_1 tile
                  compression, written as ``<name>.fz`` with an empty primary
                  HDU and the image in HDU 1 (the same layout fpack writes).
* ``decompress``  ``<name>.fz`` -> ``<name>``: the same integers, the same
                  scaling and the same header cards.

``pack`` is ``int16`` followed by ``compress``, for processed frames in one
pass.  The pipeline itself only ever reads decompressed frames: a job
decompresses a night before it runs and compresses it again afterwards.

How the 16-bit window is chosen
-------------------------------
A 16-bit integer holds 65,536 values.  FITS stores physical = stored + BZERO
(BSCALE 1), so BZERO places each frame's window on the ADU scale.  Stored
-32,768 is reserved for BLANK, which leaves a window 65,534 ADU wide.

The window starts at the frame's faintest finite pixel, so sky is never
clipped; only saturated cores above floor + 65,534 can be.  A frame spanning
more than that has its floor raised to min(max - 65,534, -1,001): CASU
already treats pixels below -1,000 as bad, so clipping those changes nothing.

BZERO is capped at 32,767 and kept away from three values:

* 32,768, which makes astropy read the frame as uint16, so ``data - sky``
  wraps round;
* 0, which makes astropy and fitsio return raw int16 arrays, so arithmetic on
  them can overflow;
* above 32,767 or below -32,667, where CASU imstack overflows: it copies the
  frame's BZERO onto its 16-bit confidence map (cir_imstack_cat.c:579).

Non-finite pixels (NaN, +-inf) are stored as BLANK = -32,768, which astropy
reads back as NaN.

Safety
------
Nothing on the servers is backed up, so every conversion writes a hidden
temporary file beside the original (a dot-file, which the pipeline's globs
ignore), reads it back and checks it, gives it the original's mode, group and
mtime, flushes it to disk, and only then renames it into place.  When the name
changes, the original is removed only after that.  Compression and
decompression must reproduce the physical pixel values exactly, with the same
scaling and the same header cards; the 16-bit conversion must keep every
unclipped pixel within 0.5 ADU and every blank pixel blank.  A file that fails
any check is left untouched.

Nothing here takes the job queue's locks: callers must not convert a night
while the pipeline is working on it.

Usage
-----
::

    python -m utils.framepack pack <target>/1/ --manifest framepack.csv
    python -m utils.framepack decompress <night dir> --workers 8

Directories are expanded to ``proc*.fits`` (int16, pack), ``*.fits``,
``*.fts`` and ``*.fit`` (compress) or ``*.fz`` (decompress), not recursively.
"""
import argparse
import csv
import datetime
import glob
import math
import os
import stat
import sys
import warnings
from dataclasses import asdict, dataclass, fields
from multiprocessing import Pool

import numpy as np
from astropy.io import fits

BLANK = -32768            # stored value of a NaN or infinite pixel
STORED_MIN = -32767       # the usable stored range, BLANK excluded
STORED_MAX = 32767
SPAN = STORED_MAX - STORED_MIN    # 65,534 ADU: the widest frame stored unclipped
MAX_BZERO = 32767         # 32,768 reads as uint16; above 32,767 CASU imstack's confidence map overflows
MIN_BZERO = -32667        # below this the confidence map overflows the other way
CASU_BAD_BELOW = -1000    # CASU imcore_list flags pixels below -1,000 as bad

SCALING = ('BSCALE', 'BZERO', 'BLANK')
STRUCTURAL = ('SIMPLE', 'XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND', 'PCOUNT', 'GCOUNT',
              'CHECKSUM', 'DATASUM')
NOT_COMPARED = set(STRUCTURAL) | set(SCALING)
FITS_SUFFIXES = ('.fits', '.fts', '.fit')
HISTORY_TAG = 'framepack:'
OPS = ('int16', 'compress', 'decompress', 'pack')
DIRECTORY_GLOBS = {'int16': ('proc*.fits',), 'pack': ('proc*.fits',),
                   'compress': ('*.fits', '*.fts', '*.fit'), 'decompress': ('*.fz',)}


class FramePackError(Exception):
    """A frame that can't be converted safely.  The original is left untouched."""


# --------------------------------------------------------------------------- 16-bit

def choose_bzero(lo, hi):
    """BZERO for a frame whose finite pixels run from ``lo`` to ``hi`` ADU."""
    floor = math.floor(lo)
    top = math.ceil(hi)
    if top - floor > SPAN:
        floor = max(floor, min(top - SPAN, CASU_BAD_BELOW - 1))
    bzero = min(MAX_BZERO, floor - STORED_MIN)
    if bzero < MIN_BZERO:
        raise FramePackError(f'pixels {lo:.1f} to {hi:.1f} ADU would need BZERO {bzero}, below {MIN_BZERO}')
    if bzero == 0:
        bzero = -1        # never 0: readers would return raw int16
    return bzero


def quantise(data):
    """Round a float frame to whole ADU and store it as 16-bit integers.

    Returns ``(stored, bzero, n_clipped, n_blank)``, ``stored`` being int16
    with physical = stored + bzero and BLANK where the input was not finite.
    """
    x = np.asarray(data, dtype=np.float64)
    finite = np.isfinite(x)
    n_blank = int(x.size - np.count_nonzero(finite))
    if n_blank == x.size:
        raise FramePackError('the frame has no finite pixels')
    values = x[finite]
    bzero = choose_bzero(values.min(), values.max())
    with np.errstate(invalid='ignore'):
        whole = np.sign(x) * np.floor(np.abs(x) + 0.5)
        s = whole - bzero
        clipped = finite & ((s < STORED_MIN) | (s > STORED_MAX))
    s = np.clip(s, STORED_MIN, STORED_MAX)
    s[~finite] = BLANK
    return s.astype(np.int16), bzero, int(np.count_nonzero(clipped)), n_blank


def _set_scaling(header, bscale, bzero, blank):
    """Put the scaling cards where FITS expects them, after the last NAXISn."""
    last = f"NAXIS{header['NAXIS']}"
    header.set('BSCALE', bscale, 'physical = stored * BSCALE + BZERO', after=last)
    header.set('BZERO', bzero, after='BSCALE')
    if blank is not None:
        header.set('BLANK', blank, 'stored value of a blank (NaN) pixel', after='BZERO')


def _strip(header, keywords):
    header = header.copy()
    for k in keywords:
        header.remove(k, ignore_missing=True, remove_all=True)
    return header


def to_int16(src, dst):
    """Write ``src`` (one float image) to ``dst`` as 16-bit integers; returns (bzero, n_clipped, n_blank)."""
    with fits.open(src) as h:
        if len(h) != 1 or h[0].header.get('NAXIS', 0) == 0:
            raise FramePackError('not a single-image FITS file')
        if h[0].header['BITPIX'] > 0:
            raise FramePackError(f"already an integer image (BITPIX {h[0].header['BITPIX']})")
        header = h[0].header.copy()
        stored, bzero, n_clipped, n_blank = quantise(h[0].data)
    out = fits.PrimaryHDU(data=stored, header=_strip(header, STRUCTURAL + SCALING), do_not_scale_image_data=True)
    _set_scaling(out.header, 1.0, bzero, BLANK if n_blank else None)
    # one card (<= 72 characters), so the header comparison can recognise and skip it
    out.header.add_history(f'{HISTORY_TAG} float32 -> int16, {n_clipped} clipped, {n_blank} blank, '
                           f'{datetime.date.today().isoformat()}')
    out.writeto(dst, checksum=True)
    return bzero, n_clipped, n_blank


# --------------------------------------------------------------------------- Rice

def compress(src, dst):
    """Write the 16-bit image in ``src`` to ``dst`` with lossless RICE_1 tile compression."""
    with fits.open(src, do_not_scale_image_data=True) as h:
        if len(h) != 1 or h[0].header.get('NAXIS', 0) == 0:
            raise FramePackError('not a single-image FITS file (already compressed?)')
        bitpix = h[0].header['BITPIX']
        if bitpix != 16:
            raise FramePackError(f'BITPIX {bitpix}: only 16-bit images are compressed losslessly')
        header = _strip(h[0].header, ('CHECKSUM', 'DATASUM'))
        c = fits.CompImageHDU(data=h[0].data, header=header, compression_type='RICE_1',
                              do_not_scale_image_data=True)
        fits.HDUList([fits.PrimaryHDU(), c]).writeto(dst, checksum=True)


def decompress(src, dst):
    """Write the Rice-compressed image in ``src`` to ``dst`` as a plain 16-bit FITS file."""
    with fits.open(src, do_not_scale_image_data=True) as h:
        if len(h) != 2 or not isinstance(h[1], fits.CompImageHDU):
            raise FramePackError('not a tile-compressed frame (primary + one compressed image)')
        header = h[1].header
        scaling = (header.get('BSCALE', 1.0), header.get('BZERO', 0.0), header.get('BLANK'))
        stored = h[1].data
        out = fits.PrimaryHDU(data=stored, header=_strip(header, STRUCTURAL + SCALING), do_not_scale_image_data=True)
    _set_scaling(out.header, *scaling)
    out.writeto(dst, checksum=True)


# --------------------------------------------------------------------------- checks

def _image(path, verify_checksum=False):
    """(header, physical pixels as float64, stored integers) of the one image in a plain or compressed file.

    Only files written here are checked against their CHECKSUM/DATASUM: current frames carry stale ones
    copied from the raw frames.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with fits.open(path, checksum=verify_checksum) as h:
            hdu = h[1] if len(h) == 2 and isinstance(h[1], fits.CompImageHDU) else h[0]
            header = hdu.header.copy()
            physical = np.array(hdu.data, dtype=np.float64)
    if verify_checksum:
        bad = [str(w.message) for w in caught if 'checksum' in str(w.message).lower()
               or 'datasum' in str(w.message).lower()]
        if bad:
            raise FramePackError(f'checksum check failed: {bad[0]}')
    with fits.open(path, do_not_scale_image_data=True) as h:
        hdu = h[1] if len(h) == 2 and isinstance(h[1], fits.CompImageHDU) else h[0]
        stored = np.array(hdu.data)
    return header, physical, stored


def _cards(header):
    return [str(c).rstrip() for c in header.cards
            if c.keyword not in NOT_COMPARED and not str(c.value).startswith(HISTORY_TAG)]


def _scaling(header):
    return header.get('BSCALE', 1.0), header.get('BZERO', 0.0), header.get('BLANK')


def check_lossless(src, dst):
    """Raise unless ``dst`` holds exactly the same frame as ``src``: pixels, scaling and header cards."""
    hs, ps, ss = _image(src)
    hd, pd, sd = _image(dst, verify_checksum=True)
    if ps.shape != pd.shape or not np.array_equal(ps, pd, equal_nan=True):
        raise FramePackError('physical pixel values differ')
    if not np.array_equal(ss.astype(np.int64), sd.astype(np.int64)):
        raise FramePackError('stored integers differ')
    if _scaling(hs) != _scaling(hd):
        raise FramePackError(f'scaling differs: {_scaling(hs)} vs {_scaling(hd)}')
    if _cards(hs) != _cards(hd):
        raise FramePackError('header cards differ')


def check_int16(src, dst, bzero, n_clipped, n_blank):
    """Raise unless ``dst`` is ``src`` rounded to whole ADU, within 0.5 ADU except the counted clipped pixels."""
    hs, ps, _ = _image(src)
    hd, pd, _ = _image(dst, verify_checksum=True)
    if hd['BITPIX'] != 16 or hd.get('BZERO') != bzero or hd.get('BSCALE', 1.0) != 1.0:
        raise FramePackError('16-bit scaling cards are wrong')
    finite = np.isfinite(ps)
    if not np.array_equal(np.isnan(pd), ~finite):
        raise FramePackError('blank pixels do not match the non-finite input pixels')
    if int(np.count_nonzero(~finite)) != n_blank:
        raise FramePackError('blank pixel count is wrong')
    far = np.abs(pd[finite] - ps[finite]) > 0.5 + 1e-9
    if int(np.count_nonzero(far)) > n_clipped:
        raise FramePackError(f'{int(np.count_nonzero(far))} pixels moved more than 0.5 ADU, '
                             f'but only {n_clipped} were clipped')
    edge = (pd[finite][far] == bzero + STORED_MIN) | (pd[finite][far] == bzero + STORED_MAX)
    if not edge.all():
        raise FramePackError('a pixel moved more than 0.5 ADU without being clipped')
    if _cards(hs) != _cards(hd):
        raise FramePackError('header cards differ')


# --------------------------------------------------------------------------- in place

@dataclass
class Report:
    path: str
    op: str
    final: str = ''
    status: str = ''          # done, skipped, failed, dry-run
    reason: str = ''
    bytes_in: int = 0
    bytes_out: int = 0
    bzero: str = ''
    n_clipped: str = ''
    n_blank: str = ''
    nlink: int = 0
    when: str = ''


def final_name(path, op):
    if op in ('compress', 'pack'):
        return path + '.fz'
    if op == 'decompress':
        return path[:-3]
    return path


def _fsync(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _check_name(path, op):
    name = os.path.basename(path).lower()
    if op == 'decompress':
        if not name.endswith(tuple(s + '.fz' for s in FITS_SUFFIXES)):
            raise FramePackError('not a .fits.fz, .fts.fz or .fit.fz file')
    elif not name.endswith(FITS_SUFFIXES):
        raise FramePackError('not a .fits, .fts or .fit file')


def _convert(path, op, out, tmp_dir):
    """Run ``op`` from ``path`` into ``out``, verifying each step; returns (bzero, n_clipped, n_blank) or Nones."""
    if op == 'int16':
        bzero, n_clipped, n_blank = to_int16(path, out)
        check_int16(path, out, bzero, n_clipped, n_blank)
        return bzero, n_clipped, n_blank
    if op == 'compress':
        compress(path, out)
        check_lossless(path, out)
        return None, None, None
    if op == 'decompress':
        decompress(path, out)
        check_lossless(path, out)
        return None, None, None
    # pack: 16-bit into a second hidden file, then compress that; a frame already 16-bit is only compressed
    if fits.getheader(path)['BITPIX'] == 16:
        compress(path, out)
        check_lossless(path, out)
        return None, None, None
    mid = os.path.join(tmp_dir, '.' + os.path.basename(path) + f'.framepack-int16-{os.getpid()}')
    try:
        bzero, n_clipped, n_blank = to_int16(path, mid)
        check_int16(path, mid, bzero, n_clipped, n_blank)
        compress(mid, out)
        check_lossless(mid, out)
    finally:
        if os.path.exists(mid):
            os.remove(mid)
    return bzero, n_clipped, n_blank


def convert_in_place(path, op, dry_run=False):
    """Convert one file in place (renaming it if ``op`` changes its suffix) and return a Report."""
    rep = Report(path=path, op=op, final=final_name(path, op),
                 when=datetime.datetime.now().isoformat(timespec='seconds'))
    tmp = os.path.join(os.path.dirname(path) or '.',
                       '.' + os.path.basename(rep.final) + f'.framepack-{os.getpid()}')
    try:
        if op not in OPS:
            raise FramePackError(f'unknown operation {op}')
        _check_name(path, op)
        st = os.lstat(path)
        if not stat.S_ISREG(st.st_mode):
            raise FramePackError('not a regular file')
        rep.bytes_in, rep.nlink = st.st_size, st.st_nlink
        if st.st_uid != os.getuid():
            raise FramePackError('owned by another account')
        if not os.access(os.path.dirname(path) or '.', os.W_OK):
            raise FramePackError('directory not writable')
        if rep.final != path and os.path.exists(rep.final):
            raise FramePackError(f'{os.path.basename(rep.final)} already exists')
        if dry_run:
            rep.status = 'dry-run'
            return rep
        bzero, n_clipped, n_blank = _convert(path, op, tmp, os.path.dirname(path) or '.')
        if bzero is not None:
            rep.bzero, rep.n_clipped, rep.n_blank = str(bzero), str(n_clipped), str(n_blank)
        os.chmod(tmp, stat.S_IMODE(st.st_mode))
        try:
            os.chown(tmp, -1, st.st_gid)
        except OSError:
            pass
        os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
        _fsync(tmp)
        os.replace(tmp, rep.final)
        if rep.final != path:
            os.remove(path)
        _fsync(os.path.dirname(path) or '.')
        rep.bytes_out = os.path.getsize(rep.final)
        rep.status = 'done'
    except FramePackError as e:
        rep.status, rep.reason = ('skipped' if 'already' in str(e) or 'owned' in str(e) else 'failed'), str(e)
    except Exception as e:        # noqa: BLE001 - one bad frame must not stop a sweep
        rep.status, rep.reason = 'failed', f'{type(e).__name__}: {e}'
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return rep


# --------------------------------------------------------------------------- command line

def expand(paths, op):
    files = []
    for p in paths:
        if os.path.isdir(p):
            found = set()
            for pattern in DIRECTORY_GLOBS[op]:
                found.update(glob.glob(os.path.join(p, pattern)))
            files.extend(sorted(found))
        else:
            files.append(p)
    return files


def _run(args):
    return convert_in_place(*args)


def write_manifest(path, reports):
    new = not os.path.exists(path)
    with open(path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[x.name for x in fields(Report)])
        if new:
            w.writeheader()
        for r in reports:
            w.writerow(asdict(r))


def main(argv=None):
    ap = argparse.ArgumentParser(description='Store frames as 16-bit and/or Rice-compressed, or restore them.')
    ap.add_argument('op', choices=OPS)
    ap.add_argument('paths', nargs='+', help='files, or directories (not searched recursively)')
    ap.add_argument('--manifest', help='append one CSV row per file here')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--dry-run', action='store_true', help='check names and permissions only')
    ap.add_argument('--verbose', action='store_true', help='print a line per file')
    args = ap.parse_args(argv)

    files = expand(args.paths, args.op)
    jobs = [(f, args.op, args.dry_run) for f in files]
    if args.workers > 1:
        with Pool(args.workers) as pool:
            reports = list(pool.imap(_run, jobs, chunksize=4))
    else:
        reports = [_run(j) for j in jobs]
    if args.manifest:
        write_manifest(args.manifest, reports)

    counts = {}
    for r in reports:
        counts[r.status] = counts.get(r.status, 0) + 1
        if args.verbose or r.status == 'failed':
            print(f'{r.status:8s} {r.path}' + (f'  ({r.reason})' if r.reason else ''))
    before = sum(r.bytes_in for r in reports if r.status == 'done')
    after = sum(r.bytes_out for r in reports if r.status == 'done')
    print(f'{args.op}: {len(reports)} files, ' + ', '.join(f'{v} {k}' for k, v in sorted(counts.items()))
          + (f'; {before / 1e9:.2f} GB -> {after / 1e9:.2f} GB' if before else ''))
    return 1 if counts.get('failed') else 0


if __name__ == '__main__':
    sys.exit(main())
