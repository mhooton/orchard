"""
Provenance of the raw frames orchard writes into Observations/ (docs/versioning.md, "Raw frames").

  - A frame orchard writes at ingest is stamped as new: ORCHVER (and ORCHREL from a release) and one HISTORY card
    'orchard <version> <script>'. These are frames unpacked from a datacube and Astra frames rotated or mirrored
    on arrival (stamp_ingest).
  - A later in-place edit of a raw frame records itself with orchard_version.record_edit: ORCHMOD and one card
    'orchard <version> <script>: <what>'. These edits are the astra_transform CLI, repair_astra_wcs and
    repair_orientation.
  - A download is left byte for byte as ESO delivered it. A raw frame with neither ORCHVER nor ORCHMOD is ESO's
    content as delivered, or was written before versioning (ASTRAROT/ASTRAMIR and HISTORY tell that story).

Every orchard write leaves CHECKSUM and DATASUM valid. A writer that writes the whole file passes checksum=True to
astropy. A writer that replaces header bytes only (repair_astra_wcs) uses set_checksum() below, which recomputes
CHECKSUM from the new header and a DATASUM it can trust. T7's header writes to calibrated frames use the same
helpers (astrom/header_write.py).

Until this module, every Astra frame kept the CHECKSUM/DATASUM ESO computed for the untransformed image
(astropy's writeto(checksum=False) leaves cards it was given), so they did not verify. The exceptions were light
frames that T7 later wrote a WCS into: astropy's update mode recomputes cards that are present.
"""
import datetime
import os
import sys

import numpy as np
from astropy.io import fits

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import orchard_version  # noqa: E402

NEGATIVE_ZERO = 0xFFFFFFFF          # the ones' complement sum of an HDU whose CHECKSUM is valid
CHUNK = 1 << 22                     # bytes read at a time when summing a data unit (a multiple of 4)
BLOCK = 2880                        # a FITS block
# ASCII punctuation the FITS checksum encoding steps around (Seaman, Pence & Rots, the FITS checksum convention)
_EXCLUDE = tuple(range(0x3a, 0x41)) + tuple(range(0x5b, 0x61))


def stamp_ingest(header, script):
    """Stamp an astropy Header of a raw frame written at ingest, and return it: orchard_version.stamp() (ORCHVER,
    ORCHREL from a release, no ORCHMOD) and the card 'orchard <version> <script>' unless the header holds it."""
    orchard_version.stamp(header)
    card = orchard_version.edit_card(script)
    if card not in [str(h) for h in header.get('HISTORY', [])]:
        header.add_history(card)
    return header


# --- FITS checksums by hand, for a header written as bytes ---------------------------------------------------------

def header_bytes(path):
    """The primary header exactly as stored: to the END card, padded to whole blocks."""
    chunks = []
    with open(path, 'rb') as f:
        while True:
            block = f.read(BLOCK)
            if len(block) < BLOCK:
                raise ValueError('{}: no END card before the end of the file'.format(path))
            chunks.append(block)
            if any(block[i:i + 8] == b'END     ' for i in range(0, BLOCK, 80)):
                return b''.join(chunks)


def read_header_bytes(path):
    """The primary header exactly as stored (header_bytes), and its parse."""
    raw = header_bytes(path)
    return raw, fits.Header.fromstring(raw.decode('ascii'))


def ones_add(a, b):
    """32-bit ones' complement addition."""
    s = a + b
    while s >> 32:
        s = (s & 0xFFFFFFFF) + (s >> 32)
    return s


def ones_sum(buf, start=0):
    """The 32-bit ones' complement sum of buf (big-endian words; a length that is not a multiple of 4 is padded
    with zeros, as the FITS block padding would be), added to start."""
    if len(buf) % 4:
        buf = bytes(buf) + b'\0' * (4 - len(buf) % 4)
    return ones_add(start, int(np.frombuffer(buf, dtype='>u4').sum(dtype=np.uint64)))


def encode_checksum(total):
    """The 16-character CHECKSUM value that makes an HDU whose bytes sum to `total` (with CHECKSUM set to sixteen
    '0') sum to negative zero."""
    value = ~total & 0xFFFFFFFF
    asc = [0] * 16
    for i in range(4):
        byte = (value >> (24 - 8 * i)) & 0xFF
        quotient, remainder = byte // 4 + 0x30, byte % 4
        ch = [quotient + remainder, quotient, quotient, quotient]
        moved = True
        while moved:
            moved = False
            for k in _EXCLUDE:
                for j in (0, 2):
                    if ch[j] == k or ch[j + 1] == k:
                        ch[j] += 1
                        ch[j + 1] -= 1
                        moved = True
        for j in range(4):
            asc[4 * j + i] = ch[j]
    text = ''.join(map(chr, asc))
    return text[-1] + text[:-1]


def data_bytes(header):
    """Size in bytes of a primary HDU's data unit, without the block padding."""
    naxis = int(header.get('NAXIS', 0))
    if naxis == 0:
        return 0
    n = abs(int(header['BITPIX'])) // 8
    for k in range(1, naxis + 1):
        n *= int(header['NAXIS{}'.format(k)])
    return n


def data_sum(path, offset, header):
    """The DATASUM of the primary data unit of `path`, which starts at byte `offset`. Reads the data as bytes;
    nothing is decoded."""
    remaining, total = data_bytes(header), 0
    with open(path, 'rb') as f:
        f.seek(offset)
        while remaining:
            chunk = f.read(min(remaining, CHUNK))
            if not chunk:
                raise IOError('{}: data unit ends early'.format(path))
            remaining -= len(chunk)
            total = ones_sum(chunk, total)
    return total


def array_datasum(data):
    """The DATASUM of the data unit astropy writes for `data`, a numpy array of floats or of signed integers wider
    than a byte (which astropy writes unscaled), from memory: nothing is read from disk. A 32-bit word has the same
    value whichever byte order holds it, so 4-byte items are summed as they are; other sizes go through their
    big-endian bytes (a copy)."""
    a = np.ascontiguousarray(data)
    if a.dtype.kind not in 'fi' or a.dtype.itemsize < 2:
        raise ValueError('astropy writes {} data scaled; sum the file instead'.format(a.dtype))
    if a.dtype.itemsize == 4:
        words = a.view(np.dtype('u4').newbyteorder(a.dtype.byteorder))
        return ones_add(0, int(words.sum(dtype=np.uint64)))
    return ones_sum(a.astype(a.dtype.newbyteorder('>')).tobytes())


def claimed_datasum(raw_header, header):
    """The DATASUM a header claims, when its CHECKSUM agrees with it, else None. Reads header bytes only, so it
    cannot tell a valid pair from one that went stale together, which no writer produces: astropy computes both
    from the data, and a writer that changed the data without them also changed the header."""
    if 'CHECKSUM' not in header:
        return None
    try:
        datasum = int(str(header['DATASUM']).strip())
    except (KeyError, ValueError):
        return None
    return datasum if ones_sum(raw_header, datasum) == NEGATIVE_ZERO else None


def _now():
    return datetime.datetime.now().isoformat()[:19]       # as astropy's comments have it


def set_checksum(header, datasum):
    """Set DATASUM and CHECKSUM in an astropy Header for a data unit whose sum is `datasum`, so that the header as
    header.tostring(padding=True) serialises it, followed by that data unit, verifies. A DATASUM card that already
    holds this value is kept as it is. Returns the header."""
    if str(header.get('DATASUM', '')).strip() != str(datasum):
        header['DATASUM'] = (str(datasum), 'data unit checksum updated ' + _now())
    header['CHECKSUM'] = ('0' * 16, 'HDU checksum updated ' + _now())
    header['CHECKSUM'] = encode_checksum(ones_sum(header.tostring(padding=True).encode('ascii'), datasum))
    if ones_sum(header.tostring(padding=True).encode('ascii'), datasum) != NEGATIVE_ZERO:
        raise ValueError('CHECKSUM did not verify after it was set')
    return header
