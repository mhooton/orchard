"""
Writes to a calibrated frame that leave CHECKSUM and DATASUM valid without astropy re-summing the pixels: T6's
calibrated frame (writeto) and T7's edits to its header (updating: pointer_wcs after a solve, wcs_cache.apply on a
cache hit). T7's header is read and written as bytes and the pixels are never decoded, so a 16-bit frame stays 16-bit.

astropy recomputes CHECKSUM and DATASUM when it writes a header that has them (update mode on close, writeto with
checksum=True). For DATASUM it sums the data unit in 2880-byte steps, in Python and holding the GIL: on appct, 136 ms
for a 16.7 MB Andor frame and 48 ms for a 5.3 MB SPIRIT frame. Instead:

  - T6 sums the float32 array it has just written, in memory (raw_provenance.array_datasum), and sets both cards in
    the file it wrote. Before, it wrote with checksum=False, which kept the raw frame's CHECKSUM/DATASUM: ESO's, for
    the raw header and pixels, so stale. Frames without the cards (Artemis) now get them too.
  - T7 takes DATASUM from the header when the old CHECKSUM agrees with it (raw_provenance.claimed_datasum, as
    repair_astra_wcs does for raw frames), which reads the header only. That is every frame T6 or T7 has written
    since. Otherwise (frames calibrated by older code) it sums the data bytes in 4 MB chunks (raw_provenance.data_sum):
    ~10 ms for 16.7 MB.
  - CHECKSUM is computed from the header and that DATASUM (raw_provenance.set_checksum).
"""
import os
import shutil
from contextlib import contextmanager

import numpy as np

from download.raw_provenance import (CHUNK, array_datasum, claimed_datasum, data_sum, read_header_bytes,
                                     set_checksum)


def writeto(hdulist, path):
    """
    T6: hdulist.writeto(path), with the primary HDU's CHECKSUM and DATASUM valid. Room for the two cards is made in
    the header first (placeholders where the header lacks them), so setting them afterwards is a header write in place.
    DATASUM comes from the data in memory when it is a plain float32 array, as T6's is, and otherwise from the file.
    """
    header = hdulist[0].header
    if 'DATASUM' not in header:
        header['DATASUM'] = ('0', 'data unit checksum')
    if 'CHECKSUM' not in header:
        header.set('CHECKSUM', '0' * 16, 'HDU checksum', before='DATASUM')
    hdulist.writeto(path)
    data = hdulist[0].data
    set_checksums(path, array_datasum(data) if type(data) is np.ndarray and data.dtype == np.float32 else None)


def set_checksums(path, datasum=None):
    """Set the primary header's CHECKSUM and DATASUM for a data unit whose sum is `datasum` (known to the caller) or,
    when it is None, the sum of the data bytes. Adds the cards if the header lacks them. Header only."""
    old_raw, header = read_header_bytes(path)
    if datasum is None:
        datasum = data_sum(path, len(old_raw), header)
    set_checksum(header, datasum)
    _replace(path, header.tostring(padding=True).encode('ascii'), old_raw)


@contextmanager
def updating(path):
    """Yield the primary header of the FITS file `path` for editing, in place of astropy's update mode. On leaving the
    block, also through a return or an exception, the header is written back if it changed (write())."""
    old_raw, header = read_header_bytes(path)
    trusted = claimed_datasum(old_raw, header)          # before the block edits the header
    try:
        yield header
    finally:
        write(path, header, old_raw, trusted)


def write(path, header, old_raw, datasum=None):
    """
    Write `header` over the primary header of `path`, whose bytes on disk are `old_raw`. Returns False, having written
    nothing, if the header is unchanged.

    If the header has CHECKSUM or DATASUM, both are set for `datasum` (the data unit's sum, which the caller trusts)
    or, when it is None, for the sum of the data bytes.
    """
    new_raw = header.tostring(padding=True).encode('ascii')
    if new_raw == old_raw:
        return False
    if 'CHECKSUM' in header or 'DATASUM' in header:
        if datasum is None:     # a header-only edit leaves BITPIX and NAXISn, so the data unit, as they were
            datasum = data_sum(path, len(old_raw), header)
        set_checksum(header, datasum)
        new_raw = header.tostring(padding=True).encode('ascii')
    _replace(path, new_raw, old_raw)
    return True


def _replace(path, new_raw, old_raw):
    """Replace the header bytes `old_raw` at the start of `path` with `new_raw`. In place when both fill the same
    number of 2880-byte blocks. Otherwise the file is rewritten through a temporary file, with the bytes after the
    header copied unchanged, and swapped in with os.replace, keeping its mode."""
    if len(new_raw) == len(old_raw):
        with open(path, 'r+b') as f:
            f.write(new_raw)
        return
    size = os.path.getsize(path)
    tmp = '{}.tmp{}'.format(path, os.getpid())
    try:
        with open(path, 'rb') as src, open(tmp, 'wb') as dst:
            dst.write(new_raw)
            src.seek(len(old_raw))
            shutil.copyfileobj(src, dst, CHUNK)
        if os.path.getsize(tmp) != size - len(old_raw) + len(new_raw):
            raise IOError('{}: the rewritten file has the wrong size; left untouched'.format(path))
        os.chmod(tmp, os.stat(path).st_mode & 0o7777)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
