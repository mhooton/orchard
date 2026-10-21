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
pass.

* ``rawpack``     a raw frame -> lossless RICE_1, byte-exact: the ``.fz``
                  also keeps the original header block verbatim (extension
                  ``ORIGHDR``, with the original file's sha256, size and
                  mtime), so ``decompress`` rebuilds the file as it was
                  downloaded, byte for byte, and ESO's own CHECKSUM verifies
                  on it.  It is a storage encoding, not an edit: the image
                  header is left exactly as it was (Matthew, 2026-10-21: all
                  raw, ESO downloads included, on these terms).

Readers
-------
The pipeline works on plain frames.  For raw frames, ZLP_pipeline.sh
decompresses the night before T1 and restores the compressed files afterwards
(``utils/rawpack.py``); a rerun from T7-T9 decompresses its calibrated frames
first.  Code that reads a frame that may be stored either way, for example a
calibration frame from an earlier night, uses ``open_frame``: it presents a
``.fz`` frame as the plain file it stands for (HDU 0 is the image, with the
dtype astropy gives the plain file), and follows a list entry that still names
``X.fits`` to ``X.fits.fz`` once that frame has been compressed.  It never
opens anything for writing.

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
    python -m utils.framepack rawpack <raw night dir> --workers 8
    python -m utils.framepack decompress <night dir> --workers 8

Directories are expanded to ``proc*.fits`` (int16, pack), ``*.fits``,
``*.fts`` and ``*.fit`` (compress, rawpack) or ``*.fz`` (decompress), not
recursively.
"""
import argparse
import csv
import datetime
import glob
import hashlib
import io
import math
import os
import stat
import sys
import warnings
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields
from multiprocessing import Pool

import numpy as np
from astropy.io import fits

try:
    from utils import orchard_version
except ImportError:         # run as a bare script
    import orchard_version

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
# ORCHMOD: each conversion records itself as the file's latest edit, so it differs from the source by design
NOT_COMPARED = set(STRUCTURAL) | set(SCALING) | {orchard_version.EDIT_KEYWORD}
FITS_SUFFIXES = ('.fits', '.fts', '.fit')
COMPRESSED = '.fz'
ORIGINAL_HEADER = 'ORIGHDR'       # extension of a byte-exact .fz: the original header block, verbatim
HISTORY_TAG = 'framepack:'
OPS = ('int16', 'compress', 'decompress', 'pack', 'rawpack')
DIRECTORY_GLOBS = {'int16': ('proc*.fits',), 'pack': ('proc*.fits',),
                   'compress': ('*.fits', '*.fts', '*.fit'), 'rawpack': ('*.fits', '*.fts', '*.fit'),
                   'decompress': ('*.fz',)}
BLOCK = 2880


class FramePackError(Exception):
    """A frame that can't be converted safely.  The original is left untouched."""


# --------------------------------------------------------------------------- names and reading

def is_compressed(path):
    return str(path).lower().endswith(COMPRESSED)


def plain_name(path):
    """The name a frame has when it is not compressed: ``X.fits`` for ``X.fits.fz``."""
    return path[:-len(COMPRESSED)] if is_compressed(path) else path


def stored_path(path):
    """Where the frame ``path`` names is stored now: ``path`` itself, or ``path + '.fz'`` once it has been compressed.

    Frame lists written by T1 and T6 keep the plain names, so a later rerun follows them here."""
    if not is_compressed(path) and not os.path.exists(path) and os.path.exists(path + COMPRESSED):
        return path + COMPRESSED
    return path


def image_hdu(hdul):
    """The HDU that holds a frame's image: HDU 1 of a tile-compressed file, whose primary HDU is empty, else HDU 0."""
    if len(hdul) > 1 and isinstance(hdul[1], fits.CompImageHDU) and not hdul[0].header.get('NAXIS', 0):
        return hdul[1]
    return hdul[0]


def plain_bytes(path):
    """The plain single-image FITS file a tile-compressed frame stands for, as bytes.

    A byte-exact .fz (``rawpack``) gives back the original file byte for byte, checked against the sha256 recorded
    when it was compressed.  Any other .fz (the processed frames, the 2019-20 raw nights compressed with fpack in
    2022) gives a file with the same stored integers, scaling and header cards."""
    with fits.open(path, do_not_scale_image_data=True) as h:
        hdu = image_hdu(h)
        if hdu is h[0]:
            raise FramePackError('not a tile-compressed frame')
        stored = np.array(hdu.data)
        if ORIGINAL_HEADER in h:
            return _exact_bytes(h[ORIGINAL_HEADER], stored)
        header = hdu.header
        scaling = (header.get('BSCALE', 1.0), header.get('BZERO', 0.0), header.get('BLANK'))
        out = fits.PrimaryHDU(data=stored, header=_strip(header, STRUCTURAL + SCALING), do_not_scale_image_data=True)
    _set_scaling(out.header, *scaling)
    buf = io.BytesIO()
    out.writeto(buf, checksum=True)
    return buf.getvalue()


@contextmanager
def open_frame(path, **kwargs):
    """Open a frame read-only, however it is stored, as an HDUList whose HDU 0 is the image.

    A plain file is opened as it is.  A .fz frame is decompressed in memory and opened as the plain file it stands
    for, so ``hdul[0].data`` has the dtype astropy gives the plain file (uint16 for raw frames) and every
    ``hdul[0]`` reader works on it unchanged; ``hdul.filename()`` is the .fz file.  A name that no longer exists is
    followed to its compressed copy (``stored_path``)."""
    path = stored_path(path)
    if not is_compressed(path):
        with fits.open(path, **kwargs) as h:
            yield h
        return
    buf = io.BytesIO(plain_bytes(path))
    buf.name = path
    with fits.open(buf, **kwargs) as h:
        yield h


def getheader(path):
    """The image header of a frame stored plain or compressed (``fits.getheader`` gives a .fz's empty primary)."""
    path = stored_path(path)
    if is_compressed(path):
        with fits.open(path) as h:
            return image_hdu(h).header.copy()
    return fits.getheader(path)


# --------------------------------------------------------------------------- byte-exact raw frames

def _sha256(data):
    return hashlib.sha256(data).hexdigest()


def _padded(n):
    return -(-n // BLOCK) * BLOCK


def _split_plain(raw):
    """(header block, header, shape) of a plain single-image FITS file held in ``raw``.

    Refuses anything that would not come back byte for byte from its header and pixels alone: extensions or other
    bytes after the image, non-zero padding, or an image that is not 16-bit integers."""
    end = next((k for k in range(0, min(len(raw), 2000 * 80), 80) if raw[k:k + 8] == b'END     '), None)
    if end is None:
        raise FramePackError('no END card: not a FITS file')
    nhead = _padded(end + 80)
    header = fits.Header.fromstring(raw[:nhead].decode('ascii', errors='replace'))
    naxis = header.get('NAXIS', 0)
    if not naxis or header.get('XTENSION'):
        raise FramePackError('not a single-image FITS file (already compressed?)')
    if header['BITPIX'] != 16:
        raise FramePackError(f"BITPIX {header['BITPIX']}: only 16-bit raw frames are compressed byte-exact")
    shape = tuple(header[f'NAXIS{k}'] for k in range(naxis, 0, -1))
    nbytes = 2 * int(np.prod(shape))
    if len(raw) != nhead + _padded(nbytes):
        raise FramePackError(f'{len(raw) - nhead - _padded(nbytes):+d} bytes after the image: extensions or a '
                             'truncated file, not stored byte-exact')
    if raw[nhead + nbytes:].strip(b'\0'):
        raise FramePackError('the data padding is not zeros, not stored byte-exact')
    return raw[:nhead], header, shape


def _exact_bytes(orig, stored):
    """The original file from its ``ORIGHDR`` extension and the stored integers; raises unless its sha256 matches."""
    out = orig.data.tobytes() + np.ascontiguousarray(stored, dtype='>i2').tobytes()
    out += b'\0' * (_padded(len(out)) - len(out))
    if len(out) != orig.header['ORIGSIZE'] or _sha256(out) != orig.header['ORIGSHA']:
        raise FramePackError('the rebuilt file does not match the sha256 of the original')
    return out


def compress_exact(src, dst):
    """Write the raw frame ``src`` to ``dst`` as a byte-exact .fz: RICE_1 image + the original header block."""
    with open(src, 'rb') as f:
        raw = f.read()
    head, header, shape = _split_plain(raw)
    stored = np.frombuffer(raw, dtype='>i2', count=int(np.prod(shape)), offset=len(head)).reshape(shape)
    image = fits.CompImageHDU(data=stored.astype(np.int16), header=header, compression_type='RICE_1',
                              do_not_scale_image_data=True)
    orig = fits.ImageHDU(data=np.frombuffer(head, dtype=np.uint8), name=ORIGINAL_HEADER)
    orig.header['ORIGSHA'] = (_sha256(raw), 'sha256 of the original file')
    orig.header['ORIGSIZE'] = (len(raw), 'bytes of the original file')
    orig.header['ORIGMTIM'] = (os.stat(src).st_mtime_ns, 'mtime of the original file, ns since 1970')
    orig.header['ORIGNAME'] = (os.path.basename(src)[:68], 'name of the original file')
    orig.header.add_comment('framepack rawpack: this is the original header, byte for byte; the image is in HDU 1')
    orchard_version.stamp(orig.header)
    fits.HDUList([fits.PrimaryHDU(), image, orig]).writeto(dst, checksum=True)


def check_exact(src, dst):
    """Raise unless the .fz ``dst`` rebuilds ``src`` byte for byte and its own checksums verify."""
    with open(src, 'rb') as f:
        raw = f.read()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with fits.open(dst, checksum=True) as h:
            for hdu in h:
                hdu.data        # noqa: B018 - reading the data checks DATASUM
    bad = [str(w.message) for w in caught if 'checksum' in str(w.message).lower() or 'datasum' in str(w.message).lower()]
    if bad:
        raise FramePackError(f'checksum check failed: {bad[0]}')
    if plain_bytes(dst) != raw:
        raise FramePackError('the compressed file does not rebuild the original byte for byte')


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
    orchard_version.record_edit(out.header, 'framepack', history=False)    # the card above is the record
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
        header = orchard_version.record_edit(_strip(h[0].header, ('CHECKSUM', 'DATASUM')), 'framepack',
                                             history=False)
        c = fits.CompImageHDU(data=h[0].data, header=header, compression_type='RICE_1',
                              do_not_scale_image_data=True)
        fits.HDUList([fits.PrimaryHDU(), c]).writeto(dst, checksum=True)


def original_info(path):
    """The ORIGHDR cards (ORIGSHA, ORIGSIZE, ORIGMTIM, ORIGNAME) of a byte-exact .fz, or None for any other file."""
    with fits.open(path) as h:
        return h[ORIGINAL_HEADER].header.copy() if ORIGINAL_HEADER in h else None


def decompress(src, dst, record=True):
    """Write the Rice-compressed image in ``src`` to ``dst`` as a plain FITS file.

    A byte-exact .fz gives back the original file.  Any other gives the same stored integers, scaling and cards,
    with framepack recorded as the file's last edit unless ``record`` is False (a working copy the pipeline reads
    and that is deleted afterwards, utils/rawpack.py)."""
    with fits.open(src, do_not_scale_image_data=True) as h:
        if len(h) < 2 or image_hdu(h) is h[0]:
            raise FramePackError('not a tile-compressed frame (empty primary + a compressed image)')
        exact = ORIGINAL_HEADER in h
        if not exact:
            if len(h) != 2:
                raise FramePackError('not a tile-compressed frame (primary + one compressed image)')
            header = h[1].header
            scaling = (header.get('BSCALE', 1.0), header.get('BZERO', 0.0), header.get('BLANK'))
            stored = h[1].data
            out = fits.PrimaryHDU(data=stored, header=_strip(header, STRUCTURAL + SCALING),
                                  do_not_scale_image_data=True)
    if exact:
        with open(dst, 'wb') as f:
            f.write(plain_bytes(src))
        return
    _set_scaling(out.header, *scaling)
    if record:
        orchard_version.record_edit(out.header, 'framepack', history=False)
    out.writeto(dst, checksum=True)


def check_rebuilt(src, dst):
    """Raise unless the plain file ``dst`` is, byte for byte, the original that the byte-exact .fz ``src`` holds."""
    info = original_info(src)
    with open(dst, 'rb') as f:
        out = f.read()
    if len(out) != info['ORIGSIZE'] or _sha256(out) != info['ORIGSHA']:
        raise FramePackError('the decompressed file does not match the sha256 of the original')


# --------------------------------------------------------------------------- checks

def _image(path, verify_checksum=False):
    """(header, physical pixels as float64, stored integers) of the one image in a plain or compressed file.

    Only files written here are checked against their CHECKSUM/DATASUM: current frames carry stale ones
    copied from the raw frames.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with fits.open(path, checksum=verify_checksum) as h:
            hdu = image_hdu(h)
            header = hdu.header.copy()
            physical = np.array(hdu.data, dtype=np.float64)
    if verify_checksum:
        bad = [str(w.message) for w in caught if 'checksum' in str(w.message).lower()
               or 'datasum' in str(w.message).lower()]
        if bad:
            raise FramePackError(f'checksum check failed: {bad[0]}')
    with fits.open(path, do_not_scale_image_data=True) as h:
        hdu = image_hdu(h)
        stored = np.array(hdu.data)
    return header, physical, stored


def _cards(header):
    # Blank cards are padding: astropy leaves them where it drops BSCALE/BZERO/BLANK as T6 turns an integer raw
    # frame into float32, and reuses them for the next keyword it appends, so the output never has the same number
    return [str(c).rstrip() for c in header.cards
            if c.keyword not in NOT_COMPARED and not str(c.value).startswith(HISTORY_TAG)
            and not (c.keyword == '' and not str(c.value).strip())]


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
    if op in ('compress', 'pack', 'rawpack'):
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
    if op == 'rawpack':
        compress_exact(path, out)
        check_exact(path, out)
        return None, None, None
    if op == 'decompress':
        decompress(path, out)
        if original_info(path) is not None:
            check_rebuilt(path, out)
        else:
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
        mtime = st.st_mtime_ns
        if op == 'decompress':
            info = original_info(path)
            mtime = info['ORIGMTIM'] if info is not None else mtime
        os.utime(tmp, ns=(st.st_atime_ns, mtime))
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
