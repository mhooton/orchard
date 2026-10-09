"""What is on disk for a telescope-night, compared with what ESO holds.

- frames: top-level FITS files in Observations/<TEL>/images/<NIGHT>
- the ESO-to-disk diff is cube-aware: an ANDOR row is present when <dp_id>.fits exists (or .fits.fz, ...); a
  SPIRIT datacube row (origfile SPECU<n>.<YYYYMMDDTHHMMSS>_...) is unpacked on download into frames named by
  their own timestamps, so a cube counts as present when any frame on disk falls between its start
  and the next cube's start. Checked on Callisto 2026-10-05..07: every unpacked frame lands in a cube.
- frames stored under other names (the 2017-2018 ACP names, Sp1609-3431-S001-R001-C001-I+z.fts, also in the
  AutoFlat/ and Calibration/ subdirectories) are matched by DATE-OBS: ESO's dp_id is SPECU<n>.<DATE-OBS>
  (checked on Io 2018-09-05: DATE-OBS '2018-09-05T23:11:08.710' is SPECU1.2018-09-05T23:11:08.710).
- a stdlib FITS primary-header reader, for OBJECT, IMAGETYP, OBSERVER and the Astra keywords, and a check
  that a file is complete FITS (every HDU's header and data accounted for by its size)
- products: whether a night has light curves in PipelineOutput v2, only in v3, or nowhere
"""
import bisect
import datetime as dt
import glob
import os
import re

FITS_SUFFIXES = ('.fits', '.fts', '.fit', '.fits.fz', '.fits.gz', '.fits.Z')
DP_ID_TIME = re.compile(r'^SPECU\d+\.(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)')
DP_ID_STAMP = re.compile(r'^SPECU\d+\.(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?)$')
CUBE_ORIGFILE = re.compile(r'^SPECU\d+\.(\d{8}T\d{6})_')
FOREIGN_MATCH_SECONDS = 0.5
ARTEMIS_NAME = re.compile(r'^(?P<target>.+?)-S\d{3}-R\d{3}-C\d{3,}-')
CALIBRATION_NAMES = re.compile(r'^(bias|dark|flat|autoflat|dusk|dawn|calib)', re.I)


def night_dir(basedir, telescope, night):
    return os.path.join(basedir, 'Observations', telescope, 'images', night)


def list_frames(path):
    try:
        return sorted(n for n in os.listdir(path) if n.endswith(FITS_SUFFIXES) and not n.startswith('.'))
    except OSError:
        return []


def strip_suffix(name):
    for s in sorted(FITS_SUFFIXES, key=len, reverse=True):
        if name.endswith(s):
            return name[:-len(s)]
    return name


def frame_time(name):
    m = DP_ID_TIME.match(name)
    return dt.datetime.strptime(m.group(1), '%Y-%m-%dT%H:%M:%S') if m else None


def cube_start(row):
    m = CUBE_ORIGFILE.match(row.get('origfile') or '')
    return dt.datetime.strptime(m.group(1), '%Y%m%dT%H%M%S') if m else None


def parse_stamp(s):
    """'2018-09-05T23:11:08.710' (or without fraction) as a datetime, else None."""
    s = str(s or '').strip()
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return dt.datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


def dp_id_stamp(dp_id):
    m = DP_ID_STAMP.match(dp_id or '')
    return m.group(1) if m else None


def is_dp_id_name(name):
    """True for SPECU<n>.<timestamp> names: ESO dp_ids and frames unpacked from SPIRIT datacubes."""
    return DP_ID_TIME.match(os.path.basename(name)) is not None


def diff_eso_disk(rows, names, foreign=None):
    """Split ESO rows into (present, missing) given the frame names on disk.

    foreign maps the DATE-OBS of frames stored under other names (see foreign_frames) to their names; a row
    whose dp_id time matches one, exactly or within FOREIGN_MATCH_SECONDS, is present.
    """
    stems = {strip_suffix(n) for n in names}
    cubes = sorted((cube_start(r), r['dp_id']) for r in rows if cube_start(r) is not None)
    counts = {}
    if cubes:
        starts = [c[0] for c in cubes]
        for n in names:
            t = frame_time(n)
            if t is None:
                continue
            i = bisect.bisect_right(starts, t + dt.timedelta(seconds=1)) - 1
            if i >= 0:
                counts[cubes[i][1]] = counts.get(cubes[i][1], 0) + 1
    match = ForeignMatcher(foreign) if foreign else None
    present, missing = [], []
    for r in rows:
        ok = r['dp_id'] in stems or (cube_start(r) is not None and counts.get(r['dp_id'], 0) > 0)
        if not ok and match is not None:
            ok = match.find(dp_id_stamp(r['dp_id'])) is not None
        (present if ok else missing).append(r)
    return present, missing, counts


def row_files(rows, names):
    """{dp_id: [names on disk]}: the file named after the row, or, for a datacube row, the frames unpacked
    from it (those between its start and the next cube's start, as in diff_eso_disk)."""
    by_stem = {}
    for n in names:
        by_stem.setdefault(strip_suffix(n), []).append(n)
    out = {r['dp_id']: list(by_stem.get(r['dp_id'], [])) for r in rows}
    cubes = sorted((cube_start(r), r['dp_id']) for r in rows if cube_start(r) is not None)
    if cubes:
        starts = [c[0] for c in cubes]
        for n in names:
            t = frame_time(n)
            if t is None:
                continue
            i = bisect.bisect_right(starts, t + dt.timedelta(seconds=1)) - 1
            if i >= 0 and n not in out[cubes[i][1]]:
                out[cubes[i][1]].append(n)
    return out


class ForeignMatcher:
    """Find the frame stored under another name that has a given DATE-OBS. Each frame matches one row at most."""

    def __init__(self, foreign):
        self.exact = dict(foreign)
        pairs = sorted((t, n) for t, n in ((parse_stamp(s), n) for s, n in foreign.items()) if t is not None)
        self.times = [p[0] for p in pairs]
        self.names = [p[1] for p in pairs]
        self.used = set()

    def find(self, stamp):
        if not stamp:
            return None
        name = self.exact.get(stamp)
        if name is None or name in self.used:
            name = self._nearest(stamp)
        if name is not None:
            self.used.add(name)
        return name

    def _nearest(self, stamp):
        t = parse_stamp(stamp)
        if t is None or not self.times:
            return None
        i = bisect.bisect_left(self.times, t)
        best = None
        for j in (i - 2, i - 1, i, i + 1):
            if 0 <= j < len(self.times) and self.names[j] not in self.used:
                d = abs((self.times[j] - t).total_seconds())
                if d <= FOREIGN_MATCH_SECONDS and (best is None or d < best[0]):
                    best = (d, self.names[j])
        return best[1] if best else None


def foreign_frames(path):
    """FITS files in a night directory, and one level of subdirectories (AutoFlat/, Calibration/), whose names
    are not SPECU<n>.<timestamp>: the 2017-2018 ACP names. Paths relative to the night directory."""
    out = []
    try:
        entries = sorted(os.listdir(path))
    except OSError:
        return out
    for name in entries:
        if name.startswith('.'):
            continue
        full = os.path.join(path, name)
        if os.path.isdir(full):
            try:
                sub = sorted(os.listdir(full))
            except OSError:
                continue
            out += [os.path.join(name, s) for s in sub
                    if s.endswith(FITS_SUFFIXES) and not s.startswith('.') and not is_dp_id_name(s)]
        elif name.endswith(FITS_SUFFIXES) and not is_dp_id_name(name):
            out.append(name)
    return out


def foreign_times(path, names=None):
    """{DATE-OBS: name} for the frames stored under other names (read from their headers)."""
    out = {}
    for rel in foreign_frames(path) if names is None else names:
        if not rel.endswith(('.fits', '.fts', '.fit')):
            continue  # compressed: the header is not readable without decompressing
        try:
            stamp = read_header(os.path.join(path, rel), max_blocks=12).get('DATE-OBS')
        except OSError:
            continue
        if stamp:
            out.setdefault(str(stamp).strip(), rel)
    return out


# --------------------------------------------------------------------------- FITS headers

def _card_value(raw):
    raw = raw.strip()
    if raw.startswith("'"):
        out, i = [], 1
        while i < len(raw):
            if raw[i] == "'":
                if i + 1 < len(raw) and raw[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                break
            out.append(raw[i])
            i += 1
        return ''.join(out).rstrip()
    value = raw.split('/', 1)[0].strip()
    if value in ('T', 'F'):
        return value == 'T'
    for cast in (int, float):
        try:
            return cast(value)
        except ValueError:
            pass
    return value


def read_header(path, max_blocks=72):
    """Primary-HDU header keywords of a FITS file (uncompressed), with the standard library only."""
    out = {}
    with open(path, 'rb') as f:
        for _ in range(max_blocks):
            block = f.read(2880)
            if len(block) < 2880:
                break
            text = block.decode('ascii', 'replace')
            for k in range(36):
                card = text[80 * k:80 * (k + 1)]
                key = card[:8].strip()
                if key == 'END':
                    return out
                if card[8:10] == '= ' and key not in out:
                    out[key] = _card_value(card[10:])
    return out


FITS_BLOCK = 2880
BITPIX_VALUES = (8, 16, 32, 64, -32, -64)
MAX_AXIS = 100000


class FitsError(ValueError):
    pass


def _header_at(f, offset):
    """Cards of the header starting at offset, and the header's length in bytes."""
    f.seek(offset)
    cards, nblocks = {}, 0
    while True:
        block = f.read(FITS_BLOCK)
        if len(block) < FITS_BLOCK:
            raise FitsError('header at byte {} has no END card'.format(offset))
        nblocks += 1
        text = block.decode('ascii', 'replace')
        for k in range(36):
            card = text[80 * k:80 * (k + 1)]
            key = card[:8].strip()
            if key == 'END':
                return cards, nblocks * FITS_BLOCK
            if card[8:10] == '= ' and key not in cards:
                cards[key] = _card_value(card[10:])
        if nblocks > 1000:
            raise FitsError('header at byte {} runs past 1000 blocks'.format(offset))


def _int_card(cards, key, hdu):
    v = cards.get(key)
    if not isinstance(v, int) or isinstance(v, bool):
        raise FitsError('HDU {}: {} is {!r}'.format(hdu, key, v))
    return v


def fits_layout(path):
    """Check that a file is complete, well-formed FITS and describe its HDUs.

    Raises FitsError unless the size is a multiple of 2,880, the first card is SIMPLE = T, every HDU has sane
    BITPIX and NAXISn, and the HDUs' headers and data account for every byte of the file. Returns
    [{'extname', 'bitpix', 'naxis': [...], 'header_bytes', 'data_bytes'}, ...].
    """
    size = os.path.getsize(path)
    if size == 0 or size % FITS_BLOCK:
        raise FitsError('size {} is not a positive multiple of 2880'.format(size))
    hdus = []
    with open(path, 'rb') as f:
        card = f.read(80).decode('ascii', 'replace')
        if card[:8].strip() != 'SIMPLE' or card[8:10] != '= ' or _card_value(card[10:]) is not True:
            raise FitsError('does not start with SIMPLE = T')
        pos = 0
        while pos < size:
            f.seek(pos)
            first = f.read(FITS_BLOCK)
            if hdus and not first.startswith(b'XTENSION'):
                if not first.strip(b'\x00 '):
                    pos += FITS_BLOCK  # trailing padding block
                    continue
                raise FitsError('no XTENSION card at byte {}'.format(pos))
            n = len(hdus)
            cards, hbytes = _header_at(f, pos)
            bitpix = _int_card(cards, 'BITPIX', n)
            if bitpix not in BITPIX_VALUES:
                raise FitsError('HDU {}: BITPIX {}'.format(n, bitpix))
            naxis = _int_card(cards, 'NAXIS', n)
            if not 0 <= naxis <= 999:
                raise FitsError('HDU {}: NAXIS {}'.format(n, naxis))
            dims = [_int_card(cards, 'NAXIS{}'.format(i), n) for i in range(1, naxis + 1)]
            if any(d < 0 or d > MAX_AXIS * 100 for d in dims):
                raise FitsError('HDU {}: NAXISn {}'.format(n, dims))
            count = 1
            for d in dims:
                count *= d
            if not naxis:
                count = 0
            if hdus:
                count = int(cards.get('GCOUNT', 1)) * (int(cards.get('PCOUNT', 0)) + count)
            dbytes = abs(bitpix) // 8 * count
            hdus.append({'extname': str(cards.get('EXTNAME', '')).strip(), 'bitpix': bitpix, 'naxis': dims,
                         'header_bytes': hbytes, 'data_bytes': dbytes})
            pos += hbytes + -(-dbytes // FITS_BLOCK) * FITS_BLOCK
    if pos != size:
        raise FitsError('HDUs need {} bytes, the file has {}'.format(pos, size))
    return hdus


def verify_frame(path):
    """fits_layout, plus what a raw frame or datacube needs: a 2- or 3-D primary image of sane size.
    Returns the layout; raises FitsError."""
    hdus = fits_layout(path)
    dims = hdus[0]['naxis']
    if len(dims) not in (2, 3) or any(d < 1 or d > MAX_AXIS for d in dims):
        raise FitsError('primary image is {}'.format('x'.join(str(d) for d in dims) or 'empty'))
    return hdus


def is_cube_layout(hdus):
    """A datacube from create_datacubes.py: it carries a METADATA table (as unpack_datacubes.is_datacube)."""
    return any(h['extname'].upper() == 'METADATA' for h in hdus[1:])


def transformation_check_would_delete(name, header, kind):
    """Mirror of download.astra_transform.transformation_check: True if SSO_download.py would delete this
    frame before re-downloading it. kind is 'ANDOR' (needs ASTRAROT) or 'SPIRIT' (needs ASTRAMIR)."""
    stem = strip_suffix(name)
    if 'SPECU' not in stem or '.Z' in stem:
        return False
    if not name.endswith('.fits'):
        return True  # transformation_check opens <stem>.fits, which does not exist: it would crash
    key = 'ASTRAROT' if kind == 'ANDOR' else 'ASTRAMIR'
    return not (key in header and len(stem.split('.')[-1]) == 3)


def transformation_check_kind(telescope, night):
    """Which check SSO_download.py applies to frames already on disk, or None (it applies none)."""
    n = int(night)
    if telescope == 'Callisto':
        return 'SPIRIT' if (20220509 < n < 20230317) or n > 20250225 else 'ANDOR'
    if telescope == 'Ganymede' and n > 20250303:
        return 'ANDOR'
    return None


# --------------------------------------------------------------------------- targets

def normalise_target(name):
    return name.strip().replace(' ', '--')


def artemis_targets(names):
    out = set()
    for n in names:
        m = ARTEMIS_NAME.match(n)
        if m and not CALIBRATION_NAMES.match(m.group('target')):
            out.add(normalise_target(m.group('target')))
    return sorted(out)


def targets_from_headers(path, names, limit=None):
    """Distinct OBJECT values of light frames (slow for big nights: one header read per frame)."""
    out = set()
    for n in names[:limit] if limit else names:
        try:
            h = read_header(os.path.join(path, n), max_blocks=12)
        except OSError:
            continue
        if str(h.get('IMAGETYP', '')).lower().startswith('light') and h.get('OBJECT'):
            out.add(normalise_target(str(h['OBJECT'])))
    return sorted(out)


def night_targets(basedir, telescope, night, eso_objects=None):
    if eso_objects:
        return sorted({normalise_target(o) for o in eso_objects})
    path = night_dir(basedir, telescope, night)
    names = list_frames(path)
    if telescope == 'Artemis':
        return artemis_targets(names)
    return targets_from_headers(path, names)


# --------------------------------------------------------------------------- products and logs

def products_state(basedir, telescope, night):
    """'v2' if any light curve is in v2 (what the portal reads), 'v3_only' if only v3 has one, else 'none'."""
    for version in ('v2', 'v3'):
        pattern = os.path.join(basedir, 'PipelineOutput', version, telescope, 'output', night, '*', '*_diff.fits')
        if glob.glob(pattern):
            return version if version == 'v2' else 'v3_only'
    return 'none'


def pipeline_logs(basedir, telescope, night):
    """Logs of earlier pipeline runs of this night (PipelineOutput/v3 or, after T12, v2)."""
    out = []
    for version in ('v3', 'v2'):
        out += glob.glob(os.path.join(basedir, 'PipelineOutput', version, telescope, 'logs', '{}_*.log'.format(night)))
    return sorted(out)


def read_transfer_log(basedir, telescope, night):
    """Frame count the observatory reports for a night (cubes, for Callisto SPIRIT), or None."""
    path = os.path.join(basedir, 'Observations', telescope, 'transfer_log.txt')
    count = None
    try:
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == night:
                    try:
                        count = int(parts[1])
                    except ValueError:
                        pass
    except OSError:
        return None
    return count


def science_count(rows):
    return sum(1 for r in rows if (r.get('dp_type') or '').upper() == 'OBJECT')
