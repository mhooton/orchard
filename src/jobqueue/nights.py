"""What is on disk for a telescope-night, compared with what ESO holds.

- frames: top-level FITS files in Observations/<TEL>/images/<NIGHT>
- the ESO-to-disk diff is cube-aware: an ANDOR row is present when <dp_id>.fits exists; a SPIRIT
  datacube row (origfile SPECU<n>.<YYYYMMDDTHHMMSS>_...) is unpacked on download into frames named by
  their own timestamps, so a cube counts as present when any frame on disk falls between its start
  and the next cube's start. Checked on Callisto 2026-10-05..07: every unpacked frame lands in a cube.
- a stdlib FITS primary-header reader, for OBJECT, IMAGETYP, OBSERVER and the Astra keywords
- products: whether a night has light curves in PipelineOutput v2, only in v3, or nowhere
"""
import bisect
import datetime as dt
import glob
import os
import re

FITS_SUFFIXES = ('.fits', '.fts', '.fit', '.fits.fz', '.fits.gz', '.fits.Z')
DP_ID_TIME = re.compile(r'^SPECU\d+\.(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)')
CUBE_ORIGFILE = re.compile(r'^SPECU\d+\.(\d{8}T\d{6})_')
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


def diff_eso_disk(rows, names):
    """Split ESO rows into (present, missing) given the frame names on disk."""
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
    present, missing = [], []
    for r in rows:
        ok = r['dp_id'] in stems or (cube_start(r) is not None and counts.get(r['dp_id'], 0) > 0)
        (present if ok else missing).append(r)
    return present, missing, counts


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
