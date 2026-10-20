"""
T7's writes to a calibrated frame's header (pointer_wcs after a solve, wcs_cache.apply on a cache hit). The header is
read and written as bytes. CHECKSUM and DATASUM are left valid, and the pixels are never decoded, so a 16-bit frame
stays 16-bit.

astropy's update mode recomputes CHECKSUM and DATASUM on close when the header has them. For DATASUM it sums the data
unit in 2880-byte steps, in Python: 47 ms for a 16 MB frame on appct, most of the 20-60 ms a cache hit cost (PR #64),
and paid again after every solve. A header edit cannot change the data, so here:

  - DATASUM is taken from the header when the old CHECKSUM agrees with it (download/raw_provenance.claimed_datasum,
    as repair_astra_wcs does for raw frames). This reads the header only.
  - Otherwise DATASUM is summed from the data bytes, in large chunks (raw_provenance.data_sum), ~10 ms for 16 MB.
  - CHECKSUM is computed from the new header and that DATASUM (raw_provenance.set_checksum).

The second case is the usual one on a fresh night. T6 (calibration/pipered.py) writes with checksum=False, which
keeps the raw frame's CHECKSUM/DATASUM: ESO's, for the raw header and pixels, so they do not agree with the calibrated
header. A frame T7 has already written verifies, so a second write trusts its DATASUM. A header with neither card
(Artemis) gets none, as with astropy.
"""
import os
import shutil
from contextlib import contextmanager

from download.raw_provenance import CHUNK, claimed_datasum, data_sum, read_header_bytes, set_checksum


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
    or, when it is None, for the sum of the data bytes. A header that still fills the same number of 2880-byte blocks
    is written in place. Otherwise the file is rewritten through a temporary file, with the bytes after the header
    copied unchanged, and swapped in with os.replace, keeping its mode.
    """
    new_raw = header.tostring(padding=True).encode('ascii')
    if new_raw == old_raw:
        return False
    if 'CHECKSUM' in header or 'DATASUM' in header:
        if datasum is None:     # a header-only edit leaves BITPIX and NAXISn, so the data unit, as they were
            datasum = data_sum(path, len(old_raw), header)
        set_checksum(header, datasum)
        new_raw = header.tostring(padding=True).encode('ascii')
    if len(new_raw) == len(old_raw):
        with open(path, 'r+b') as f:
            f.write(new_raw)
        return True
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
    return True
