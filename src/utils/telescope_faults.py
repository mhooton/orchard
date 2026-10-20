"""
telescope_faults.py
===================

Frames the telescope lost, not the pipeline: for each target of a night, the light frames exposed while the mount
had stopped tracking, and those taken while the pointing had walked away from where it sat for the rest of the block.
Read-only and stdlib only (the job queue's digest imports it), from what the pipeline already has: the raw frames'
headers, T7's solutions (the calibrated frames, its WCS cache) and T7's per-frame lines in the night's logs.
docs/telescope-faults.md has the measurements behind it and the archive census.

Tracking stopped
    Astra (SSO: Callisto from 2023-12-15, Ganymede from 2025-02-27, Io and Europa from 2026-01-28/29) writes into
    every frame the share of the exposure the mount spent tracking (TRACKING) and slewing (SLEWING), and whether it
    was parked (TELPARK). A light frame for which the mount did neither for IDLE_MIN of the exposure was exposed on
    a stopped mount, and its stars trail 15 arcsec per second. Its header RA, the mount's position, then runs at the
    sidereal rate while ALTITUDE, AZIMUTH and AIRMASS stay fixed, and Astra's own WCS is empty (CRVAL 0). Measured:
    Ganymede 2025-07-16 TOI-4336 c 00:30-03:11 UT (570 frames) and 2025-08-25 Sp1954+1801 00:32-04:42 UT (333).
    MaxIm DL frames (Artemis, and SSO before Astra) have no TRACKING and compute ALTITUDE from the target
    coordinates, so there only T7 can tell, by failing frames as "Sources trailed" (astrom/pointer_wcs.py, v4.1.0).

Pointing drift
    Each solved frame's WCS puts the block's pointing (its header target position) somewhere on the detector.
    Guided, it stays within a few pixels for the whole block. A drift is the pointing more than DRIFT_ARCMIN from
    home, where it sat longest, on DRIFT_MIN_SOLVED solved frames in a row (a running median of three, so one wrong
    solution is not a drift); the unsolved frames among them count with them. Artemis 2024-02-11 Sp1356+4342 jumped
    4.3' at 01:34 UT and nothing solved after 01:49. Solutions turned or scaled unlike the rest of the block are
    wrong ones and left out; each side of a meridian flip has its own home; an episode over within SETTLE_MINUTES of
    the block's start is the telescope settling on the field.

night() reads one night; describe(), flags() and digest_text() put its episodes into words for the night report
(T11) and the morning digest, and write_summary() leaves them where the digest reads them.
"""
import csv
import datetime as dt
import glob
import json
import math
import os
import re
from collections import Counter

from utils import promote_policy

BLOCK, CARD = 2880, 80
RAW_EXTENSIONS = ('.fits', '.fts', '.fit')

IDLE_MIN = 0.5            # share of an exposure with the mount neither tracking nor slewing: a stopped mount
TRAIL_MIN_FRAMES = 3      # T7 'Sources trailed' frames in a row that make an episode (one can be a bright star)
DRIFT_ARCMIN = 1.0        # pointing this far from home: a drift (guided, it stays within seconds of arc)
DRIFT_MIN_SOLVED = 3      # solved frames in a row beyond DRIFT_ARCMIN that make a drift
REF_RADIUS_ARCMIN = 0.5   # the block's home: where the most solved frames sat within this of each other
SETTLE_MINUTES = 10.0     # an episode over within this of a block's first frame is the pointing settling
BLOCK_GAP_MINUTES = 20.0  # without the night's other frames, a gap this long between frames starts a block
ROT_TOL = 3.0             # deg: a solution turned further from the block's orientation (or its flip) is wrong
SCALE_TOL = 0.03          # a solution whose plate scale is this far from the block's is wrong
JUMP_MINUTES = 5.0        # a drift that arrives between solved frames this close together is a jump
SIDEREAL_DEG_PER_H = 15.041
LOW_ALTITUDE = 25.0       # deg: a stop this low may be the mount's altitude limit; the report says the altitude

TRAILED = 'Sources trailed'
CRPIX = re.compile(r'CRPIX=\(([-\d.eE+]+),([-\d.eE+]+)\)')
CRVAL = re.compile(r'CRVAL=\(([-\d.eE+]+),([-\d.eE+]+)\)')
GAIA = re.compile(r' Gaia: (\d+)')


# ---------------------------------------------------------------------------
# FITS headers, without astropy
# ---------------------------------------------------------------------------

def _value(text):
    """A card's value field as str, bool, int or float; None if empty."""
    s = text.strip()
    if not s:
        return None
    if s[0] == "'":
        out, i = [], 1
        while i < len(s):
            if s[i] == "'":
                if s[i + 1:i + 2] == "'":
                    out.append("'")
                    i += 2
                    continue
                break
            out.append(s[i])
            i += 1
        return ''.join(out).rstrip()
    s = s.split('/', 1)[0].strip()
    if s in ('T', 'F'):
        return s == 'T'
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s.replace('D', 'E'))
    except ValueError:
        return s or None


def _read_cards(f):
    """The header at the file's position as {KEY: value}, leaving the file at the start of the HDU's data."""
    cards = {}
    while True:
        block = f.read(BLOCK)
        if len(block) < BLOCK:
            raise ValueError('header ends before END')
        for i in range(0, BLOCK, CARD):
            card = block[i:i + CARD].decode('ascii', 'replace')
            key = card[:8].rstrip()
            if key == 'END':
                return cards
            if card.startswith('HIERARCH ') and '=' in card:
                k, v = card[9:].split('=', 1)
                cards.setdefault(k.strip(), _value(v))
            elif card[8:10] == '= ' and key not in cards:
                cards[key] = _value(card[10:])


def read_header(path):
    """The image header of a FITS file as {KEY: value}: the primary HDU's or, for a tile-compressed frame
    (.fits.fz), the compressed image's, with NAXISn the image's own. Raises OSError or ValueError."""
    with open(path, 'rb') as f:
        cards = _read_cards(f)
        if cards.get('NAXIS') == 0 and cards.get('EXTEND'):
            try:
                ext = _read_cards(f)            # no primary data, so the first extension follows at once
            except ValueError:
                return cards
            if ext.get('ZIMAGE'):
                for n in (1, 2):
                    if 'ZNAXIS{}'.format(n) in ext:
                        ext['NAXIS{}'.format(n)] = ext['ZNAXIS{}'.format(n)]
                cards.update(ext)
    return cards


# ---------------------------------------------------------------------------
# Header values
# ---------------------------------------------------------------------------

def _num(x):
    if isinstance(x, bool) or x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _sexagesimal(text):
    parts = re.split(r'[\s:]+', text.strip())
    if not parts or not parts[0]:
        return None
    sign = -1.0 if parts[0].startswith('-') else 1.0
    try:
        values = [abs(float(parts[0]))] + [float(p) for p in parts[1:3]]
    except ValueError:
        return None
    return sign * sum(v / 60.0 ** i for i, v in enumerate(values))


def angle(value, hours=False):
    """Degrees from a header angle: a number of degrees, or 'HH MM SS.ss' / 'DD:MM:SS' text (RA text is hours)."""
    if isinstance(value, str):
        s = value.strip()
        if re.search(r'\d[\s:]+\d', s):
            x = _sexagesimal(s)
            return None if x is None else x * (15.0 if hours else 1.0)
        return _num(s)
    return _num(value)


def frame_time(header):
    """Exposure start (UTC, naive) from DATE-OBS (+ TIME-OBS when DATE-OBS holds the date only), else JD/MJD-OBS."""
    d = header.get('DATE-OBS')
    if isinstance(d, str) and d.strip():
        s = d.strip()
        if 'T' not in s and isinstance(header.get('TIME-OBS'), str):
            s = s + 'T' + header['TIME-OBS'].strip()
        m = re.match(r'(\d{4})-(\d\d)-(\d\d)(?:[T ](\d\d):(\d\d):(\d\d)(?:\.(\d+))?)?$', s)
        if m:
            y, mo, da, h, mi, se, frac = m.groups()
            t = dt.datetime(int(y), int(mo), int(da), int(h or 0), int(mi or 0), int(se or 0))
            if frac:
                t += dt.timedelta(microseconds=int(round(float('0.' + frac) * 1e6)))
            return t
    for key, zero in (('JD-OBS', 2440587.5), ('JD', 2440587.5), ('MJD-OBS', 40587.0)):
        x = _num(header.get(key))
        if x:
            return dt.datetime(1970, 1, 1) + dt.timedelta(days=x - zero)
    return None


def is_light(header):
    kind = str(header.get('IMAGETYP') or '').strip().lower()
    return 'light' in kind or kind in ('object', 'science') or \
        str(header.get('ESO DPR CATG') or '').strip().upper() == 'SCIENCE'


def frame_key(name):
    """What the raw frame and its processed frame share: the raw name without its extension
    ('procSPECU3.2025-07-17T02:00:50.000.fits' and 'SPECU3.2025-07-17T02:00:50.000.fits' -> 'SPECU3.2025-...')."""
    base = os.path.basename(name)
    if base.startswith('proc'):
        base = base[4:]
    for ext in ('.fits.fz', '.fits', '.fts', '.fit'):
        if base.endswith(ext):
            return base[:-len(ext)]
    return base


# ---------------------------------------------------------------------------
# WCS: the TAN projection by hand
# ---------------------------------------------------------------------------

def header_wcs(header):
    """{'crval', 'crpix', 'cd', 'shape'} from a header with a celestial TAN WCS, else None. Astra's empty solutions
    (CRVAL 0, CTYPE 'nan') and placeholder ones (CTYPE PIXEL) are None. SIP terms are ignored: they move a star by a
    few pixels at most, nothing next to the arcminutes measured here."""
    if not str(header.get('CTYPE1') or '').startswith('RA---') or \
            not str(header.get('CTYPE2') or '').startswith('DEC--'):
        return None
    crval = (_num(header.get('CRVAL1')), _num(header.get('CRVAL2')))
    crpix = (_num(header.get('CRPIX1')), _num(header.get('CRPIX2')))
    if None in crval or None in crpix or crval == (0.0, 0.0):
        return None
    if any(k in header for k in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2')):
        cd = [[_num(header.get('CD1_1')) or 0.0, _num(header.get('CD1_2')) or 0.0],
              [_num(header.get('CD2_1')) or 0.0, _num(header.get('CD2_2')) or 0.0]]
    else:
        c1, c2 = _num(header.get('CDELT1')), _num(header.get('CDELT2'))
        if c1 is None or c2 is None:
            return None
        pc = [[_num(header.get('PC1_1')), _num(header.get('PC1_2'))],
              [_num(header.get('PC2_1')), _num(header.get('PC2_2'))]]
        pc = [[1.0 if v is None and i == j else (v or 0.0) for j, v in enumerate(row)] for i, row in enumerate(pc)]
        cd = [[c1 * pc[0][0], c1 * pc[0][1]], [c2 * pc[1][0], c2 * pc[1][1]]]
    det = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0]
    if not det or not math.isfinite(det):
        return None
    shape = (_num(header.get('NAXIS2')), _num(header.get('NAXIS1')))
    return {'crval': crval, 'crpix': crpix, 'cd': cd, 'shape': shape}


def with_reference(cd, crpix, crval, shape=(None, None)):
    """A WCS from a T7 log line's CRPIX and CRVAL and the matrix of another frame of the block."""
    return {'crval': tuple(crval), 'crpix': tuple(crpix), 'cd': cd, 'shape': shape}


def scale_arcsec(wcs):
    cd = wcs['cd']
    return math.sqrt(abs(cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0])) * 3600.0


def rotation_deg(wcs):
    """Position angle of the detector's +y axis on the sky, east of north: a meridian flip turns it by 180 deg."""
    cd = wcs['cd']
    return math.degrees(math.atan2(cd[0][1], cd[1][1])) % 360.0


def world_to_pixel(wcs, ra, dec):
    """(x, y), 1-based FITS pixels, of a sky position; None behind the projection."""
    a0, d0 = math.radians(wcs['crval'][0]), math.radians(wcs['crval'][1])
    a, d = math.radians(ra), math.radians(dec)
    cosc = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    if cosc <= 0:
        return None
    xi = math.degrees(math.cos(d) * math.sin(a - a0) / cosc)
    eta = math.degrees((math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / cosc)
    (c11, c12), (c21, c22) = wcs['cd']
    det = c11 * c22 - c12 * c21
    return (wcs['crpix'][0] + (c22 * xi - c12 * eta) / det,
            wcs['crpix'][1] + (-c21 * xi + c11 * eta) / det)


def pixel_to_world(wcs, x, y):
    """(ra, dec) in degrees of a 1-based FITS pixel position."""
    (c11, c12), (c21, c22) = wcs['cd']
    dx, dy = x - wcs['crpix'][0], y - wcs['crpix'][1]
    xi, eta = math.radians(c11 * dx + c12 * dy), math.radians(c21 * dx + c22 * dy)
    a0, d0 = math.radians(wcs['crval'][0]), math.radians(wcs['crval'][1])
    rho = math.hypot(xi, eta)
    if rho == 0:
        return wcs['crval']
    c = math.atan(rho)
    dec = math.asin(math.cos(c) * math.sin(d0) + eta * math.sin(c) * math.cos(d0) / rho)
    ra = a0 + math.atan2(xi * math.sin(c), rho * math.cos(d0) * math.cos(c) - eta * math.sin(d0) * math.sin(c))
    return math.degrees(ra) % 360.0, math.degrees(dec)


def separation_arcmin(ra1, dec1, ra2, dec2):
    a1, d1, a2, d2 = (math.radians(v) for v in (ra1, dec1, ra2, dec2))
    h = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((a2 - a1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(h)))) * 60.0


# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------

class Frame(object):
    """One raw frame's header, as far as the faults need it."""

    def __init__(self, path, header):
        self.path = path
        self.name = os.path.basename(path)
        self.key = frame_key(self.name)
        self.time = frame_time(header)
        self.object = str(header.get('OBJECT') or '').strip()
        self.light = is_light(header)
        self.exptime = _num(header.get('EXPTIME'))
        ra = header.get('RA') if header.get('RA') is not None else header.get('OBJCTRA')
        dec = header.get('DEC') if header.get('DEC') is not None else header.get('OBJCTDEC')
        self.ra, self.dec = angle(ra, hours=isinstance(ra, str)), angle(dec)
        self.airmass = _num(header.get('AIRMASS'))
        self.alt, self.az = _num(header.get('ALTITUDE')), _num(header.get('AZIMUTH'))
        self.tracking, self.slewing = _num(header.get('TRACKING')), _num(header.get('SLEWING'))
        self.parked = bool(_num(header.get('TELPARK')))
        self.wcs = header_wcs(header)
        self.solved = None          # the WCS used: processed frame, WCS cache, raw header or log line
        self.t7 = None              # T7's last attempt at the frame, if the log has it
        self.block = None           # which run of one OBJECT in the night's light frames (night() sets it)

    @property
    def idle(self):
        """Share of the exposure the mount neither tracked nor slewed (Astra), None where the header cannot say."""
        if self.tracking is None:
            return None
        if self.parked:
            return 1.0
        return max(0.0, 1.0 - self.tracking - (self.slewing or 0.0))

    @property
    def stopped(self):
        idle = self.idle
        return idle is not None and idle >= IDLE_MIN

    @property
    def trailed(self):
        return bool(self.t7 and self.t7['state'] == 'FAILED' and self.t7['error'].startswith(TRAILED))

    def __repr__(self):
        return 'Frame({!r})'.format(self.name)


def raw_frames(raw_dir):
    """Every raw frame of a night directory, light or not, sorted by exposure start. Unreadable files are skipped."""
    try:
        names = sorted(n for n in os.listdir(raw_dir) if n.endswith(RAW_EXTENSIONS))
    except OSError:
        return []
    out = []
    for n in names:
        path = os.path.join(raw_dir, n)
        try:
            out.append(Frame(path, read_header(path)))
        except (OSError, ValueError):
            continue
    return sorted(out, key=lambda f: (f.time or dt.datetime.max, f.name))


def processed_frames(output_night_dirs):
    """{frame key: path} of the processed frames in <output>/<night>/<target>/<run>/ under each directory given,
    later directories (and later runs) taking precedence."""
    out = {}
    for d in output_night_dirs:
        for path in sorted(glob.glob(os.path.join(d, '*', '*', 'proc*.fits')) +
                           glob.glob(os.path.join(d, '*', '*', 'proc*.fits.fz'))):
            if os.path.basename(os.path.dirname(path)).isdigit():
                out[frame_key(path)] = path
    return out


def cache_path(basedir, telescope, night_):
    """T7's WCS cache table of a night (astrom/wcs_cache.py; docs/wcs-cache.md)."""
    return os.path.join(basedir, 'PipelineOutput', 'wcs_cache', telescope, '{}.ecsv'.format(night_))


def wcs_cache(path):
    """{frame key: WCS} from a night's WCS cache table, read without astropy: an ECSV table (space-separated, quoted
    with "), one row per raw frame T7 solved, holding the solution's header cards as 80-character card images in the
    calibrated frame's pixels. {} if the table is missing or unreadable."""
    out = {}
    try:
        with open(path, newline='') as f:
            rows = csv.reader((line for line in f if not line.startswith('#')), delimiter=' ', quotechar='"')
            head = next(rows, None)
            if not head or 'raw' not in head or 'cards' not in head:
                return out
            col = {name: i for i, name in enumerate(head)}
            for row in rows:
                if len(row) != len(head):
                    continue
                cards, hdr = row[col['cards']], {}
                for i in range(0, len(cards), CARD):
                    card = cards[i:i + CARD]
                    if card[8:10] == '= ':
                        hdr.setdefault(card[:8].rstrip(), _value(card[10:]))
                for n in ('naxis1', 'naxis2'):
                    if n in col:
                        hdr[n.upper()] = _num(row[col[n]])
                wcs = header_wcs(hdr)
                if wcs is not None:
                    out[frame_key(row[col['raw']])] = wcs
    except (OSError, csv.Error, UnicodeDecodeError):
        return {}
    return out


def t7_outcomes(log_paths):
    """{frame key: {'state', 'sources', 'gaia', 'error', 'crpix', 'crval'}}: each frame's last attempt in T7's
    summary lines (astrometry.py) of the logs given, oldest log first. Unreadable logs are skipped."""
    out = {}
    for path in log_paths:
        try:
            with open(path, 'rb') as f:
                for raw in f:
                    if b' - Sources: ' not in raw:
                        continue
                    line = raw.decode('utf-8', 'replace').rstrip('\r\n')
                    for s in promote_policy.SUMMARY.finditer(line):
                        text, rest = s.group(0), s.group('rest')
                        gaia, crpix, crval = GAIA.search(text), CRPIX.search(rest), CRVAL.search(rest)
                        out[frame_key(s.group('name'))] = {
                            'state': s.group('state'),
                            'sources': int(s.group('src')) if s.group('src').isdigit() else None,
                            'gaia': int(gaia.group(1)) if gaia else None,
                            'error': rest.split(' Error: ', 1)[1].strip() if ' Error: ' in rest else '',
                            'crpix': (float(crpix.group(1)), float(crpix.group(2))) if crpix else None,
                            'crval': (float(crval.group(1)), float(crval.group(2))) if crval else None}
        except OSError:
            continue
    return out


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def _median(values):
    v = sorted(values)
    n = len(v)
    if not n:
        return None
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def _unwrap(ras):
    """RAs made continuous across 0/360, about the first."""
    if not ras:
        return []
    out, ref = [], ras[0]
    for r in ras:
        out.append(ref + ((r - ref + 180.0) % 360.0 - 180.0))
    return out


def blocks(frames):
    """Runs of consecutive light frames (in time) with the same OBJECT: [(object, [Frame])]."""
    out = []
    for f in frames:
        if not f.light:
            continue
        if out and out[-1][0] == f.object:
            out[-1][1].append(f)
        else:
            out.append((f.object, [f]))
    return out


def target_blocks(frames):
    """One target's frames (in time) split into its blocks: where another target's frames came between (Frame.block,
    which night() sets), or, without the night's other frames, at a gap of BLOCK_GAP_MINUTES."""
    out = []
    for f in frames:
        if out:
            prev = out[-1][-1]
            if f.block is not None and prev.block is not None:
                same = f.block == prev.block
            else:
                same = not (f.time and prev.time and f.time - prev.time > dt.timedelta(minutes=BLOCK_GAP_MINUTES))
            if same:
                out[-1].append(f)
                continue
        out.append([f])
    return out


def _runs(flags):
    """[(i0, i1)] inclusive index ranges where flags is true."""
    out, start = [], None
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(flags) - 1))
    return out


def _episode(kind, target, frames, i0, i1, index, detail):
    a, b = frames[i0], frames[i1]
    return {'kind': kind, 'target': target, 'first': a.name, 'last': b.name, 'i0': index[a.key] + 1,
            'i1': index[b.key] + 1, 'n': i1 - i0 + 1, 'start': a.time.isoformat() if a.time else None,
            'end': b.time.isoformat() if b.time else None, 'detail': detail}


def _tracking_detail(run):
    detail = {'parked': sum(1 for f in run if f.parked)}
    timed = [f for f in run if f.time and f.ra is not None]
    if len(timed) >= 2:
        hours = (timed[-1].time - timed[0].time).total_seconds() / 3600.0
        ras = _unwrap([f.ra for f in timed])
        if hours >= 1.0 / 60:
            detail['ra_rate'] = round((ras[-1] - ras[0]) / hours, 2)
    airmass = [f.airmass for f in run if f.airmass]
    if len(airmass) >= 2:
        detail['airmass'] = round(_median(airmass), 4)
        detail['airmass_span'] = round(max(airmass) - min(airmass), 6)
    alt = [f.alt for f in run if f.alt is not None]
    if alt:
        detail['altitude'] = round(_median(alt), 1)
    failed = [f for f in run if f.t7 and f.t7['state'] == 'FAILED']
    if any(f.t7 for f in run):
        detail['t7_failed'] = len(failed)
        detail['t7_trailed'] = sum(1 for f in run if f.trailed)
    return detail


def tracking_episodes(target, frames, index):
    """Frames exposed with the mount stopped (Astra's TRACKING/SLEWING/TELPARK) and, where the headers cannot say,
    runs of TRAIL_MIN_FRAMES or more frames T7 failed as trailed."""
    out = []
    for i0, i1 in _runs([f.stopped for f in frames]):
        out.append(_episode('tracking', target, frames, i0, i1, index, _tracking_detail(frames[i0:i1 + 1])))
    for i0, i1 in _runs([f.trailed and not f.stopped for f in frames]):
        if i1 - i0 + 1 >= TRAIL_MIN_FRAMES:
            run = frames[i0:i1 + 1]
            detail = _tracking_detail(run)
            detail['header_tracking'] = all(f.tracking is not None and not f.stopped for f in run)
            out.append(_episode('trailed', target, frames, i0, i1, index, detail))
    return out


def _pointing(frames):
    """The block's pointing: the median header RA/DEC of the frames the mount was not stopped for."""
    use = [f for f in frames if f.ra is not None and f.dec is not None and not f.stopped]
    if not use:
        return None
    ras = _unwrap([f.ra for f in use])
    return _median(ras) % 360.0, _median([f.dec for f in use])


def solve_frames(frames, processed=None, cached=None):
    """Set each frame's .solved, the first WCS of: its processed frame (`processed`, {key: path}: T7's solution
    for the run), T7's WCS cache (`cached`, {key: WCS}), the raw header (Astra's on-site solution, or orchard's
    from before v4.1.1), T7's log line (CRPIX, CRVAL) with the matrix of the nearest frame that has one. A
    frame the mount was stopped for gets none: its pointing means nothing."""
    processed, cached = processed or {}, cached or {}
    matrix = None
    pending = []
    for f in frames:
        f.solved = None
        if f.stopped:
            continue
        path = processed.get(f.key)
        if path:
            try:
                f.solved = header_wcs(read_header(path))
            except (OSError, ValueError):
                f.solved = None
        if f.solved is None:
            f.solved = cached.get(f.key) or f.wcs
        if f.solved is not None:
            matrix = f.solved
            for g in pending:
                g.solved = with_reference(matrix['cd'], g.t7['crpix'], g.t7['crval'], matrix['shape'])
            pending = []
        elif f.t7 and f.t7['state'] == 'SUCCESS' and f.t7['crpix'] and f.t7['crval']:
            if matrix is not None:
                f.solved = with_reference(matrix['cd'], f.t7['crpix'], f.t7['crval'], matrix['shape'])
            else:
                pending.append(f)


def _turn(a, b, period=360.0):
    """The smaller angle between two directions (period 360) or two axes (period 180), degrees."""
    d = (a - b) % period
    return min(d, period - d)


def _orientation(rotations):
    """The block's orientation, an axis in degrees mod 180 (so both sides of a meridian flip agree): the 1-degree
    bin with most solutions within ROT_TOL, refined as the mean of those solutions' doubled angles."""
    counts = [0] * 180
    for r in rotations:
        counts[int(r % 180.0) % 180] += 1
    width = int(math.ceil(ROT_TOL))
    peak = max(range(180), key=lambda b: sum(counts[(b + d) % 180] for d in range(-width, width + 1)))
    near = [r for r in rotations if _turn(r, peak + 0.5, 180.0) <= ROT_TOL + 0.5]
    s = sum(math.sin(math.radians(2 * r)) for r in near)
    c = sum(math.cos(math.radians(2 * r)) for r in near)
    return (math.degrees(math.atan2(s, c)) / 2.0) % 180.0


def _home(points, scale):
    """Where the pointing sat longest: the median of the points in the densest REF_RADIUS_ARCMIN cell and its
    neighbours (the earliest such cell on a tie). Pixels, from (x, y) points and the plate scale in arcsec/px."""
    r = REF_RADIUS_ARCMIN * 60.0 / scale
    cell = [(int(math.floor(x / r)), int(math.floor(y / r))) for x, y in points]
    counts = Counter(cell)
    best = max(dict.fromkeys(cell), key=lambda c: sum(counts.get((c[0] + a, c[1] + b), 0)
                                                      for a in (-1, 0, 1) for b in (-1, 0, 1)))
    near = [p for p, c in zip(points, cell) if abs(c[0] - best[0]) <= 1 and abs(c[1] - best[1]) <= 1]
    return _median([p[0] for p in near]), _median([p[1] for p in near])


def drift_track(frames, point):
    """[(frame index, arcmin from home, x, y, shape)]: where `point` falls on each of the block's solved frames,
    measured from where it sat longest (_home), so frames taken while the telescope still settled on the field do not
    make the rest of the block look displaced. Solutions turned more than ROT_TOL from the block's orientation, or
    with a plate scale SCALE_TOL off the block's, are wrong (old twirl solutions, accepted before the plate-scale
    check, land at random angles) and left out. Each side of a meridian flip has its own home."""
    pts = []
    for i, f in enumerate(frames):
        if f.solved is None:
            continue
        xy = world_to_pixel(f.solved, point[0], point[1])
        if xy is not None:
            pts.append((i, xy[0], xy[1], scale_arcsec(f.solved), rotation_deg(f.solved), f.solved['shape']))
    if not pts:
        return []
    axis = _orientation([p[4] for p in pts])
    scale = _median([p[3] for p in pts])
    good = [p for p in pts if _turn(p[4], axis, 180.0) <= ROT_TOL and abs(p[3] / scale - 1.0) <= SCALE_TOL]
    sides = {}
    for p in good:
        sides.setdefault(_turn(p[4], axis) > 90.0, []).append((p[1], p[2]))
    home = {side: _home(xy, scale) for side, xy in sides.items()}
    out = []
    for p in good:
        hx, hy = home[_turn(p[4], axis) > 90.0]
        out.append((p[0], math.hypot(p[1] - hx, p[2] - hy) * p[3] / 60.0, p[1], p[2], p[5]))
    return out


def drift_episodes(target, frames, index):
    """Where the block's pointing sat more than DRIFT_ARCMIN from home (drift_track), on DRIFT_MIN_SOLVED solved
    frames in a row (in a running median of three); the unsolved frames among them count with them, and an episode
    with no solved frame after it runs to the end of the block. One over within SETTLE_MINUTES of the block's first
    frame is the telescope settling on the field, not a fault. Returns (episodes, solved frames measured)."""
    point = _pointing(frames)
    if point is None:
        return [], 0
    track = drift_track(frames, point)
    if len(track) < DRIFT_MIN_SOLVED:
        return [], len(track)
    drift = [t[1] for t in track]
    smooth = [_median(drift[max(0, k - 1):k + 2]) for k in range(len(drift))]
    out = []
    for k0, k1 in _runs([d >= DRIFT_ARCMIN for d in smooth]):
        while k0 < k1 and drift[k0] < DRIFT_ARCMIN:   # the median can start a run on its last frame below
            k0 += 1
        if k1 - k0 + 1 < DRIFT_MIN_SOLVED:
            continue
        i0 = track[k0][0]
        last = track[k1][0]
        after = track[k1 + 1][0] if k1 + 1 < len(track) else None
        i1 = len(frames) - 1 if after is None else after - 1
        if frames[0].time and frames[i1].time and \
                frames[i1].time - frames[0].time <= dt.timedelta(minutes=SETTLE_MINUTES):
            continue
        peak = max(range(k0, k1 + 1), key=lambda k: drift[k])
        detail = {'max_arcmin': round(drift[peak], 1), 'max_at': _iso(frames[track[peak][0]]),
                  'solved': k1 - k0 + 1, 'last_solved_at': _iso(frames[last]),
                  'unsolved_after': after is None and last < len(frames) - 1}
        if k0 > 0:
            gap = frames[i0].time - frames[track[k0 - 1][0]].time if frames[i0].time and \
                frames[track[k0 - 1][0]].time else None
            step = drift[k0] - drift[k0 - 1]
            if step >= DRIFT_ARCMIN and gap is not None and gap <= dt.timedelta(minutes=JUMP_MINUTES):
                detail['jump_arcmin'] = round(step, 1)
        for k in range(k0, k1 + 1):
            ny, nx = track[k][4]
            x, y = track[k][2], track[k][3]
            if nx and ny and not (0.5 <= x <= nx + 0.5 and 0.5 <= y <= ny + 0.5):
                detail['off_field_at'] = _iso(frames[track[k][0]])
                break
        out.append(_episode('drift', target, frames, i0, i1, index, detail))
    return out, len(track)


def _iso(frame):
    return frame.time.isoformat() if frame.time else None


def target_faults(target, frames, processed=None, t7=None, cached=None):
    """{'frames', 'tracking_known', 'solved', 'episodes'} for one target's light frames of a night, in time order:
    how many frames, whether their headers say if the mount tracked, on how many frames the pointing was measured,
    and the episodes (tracking, trailed, drift)."""
    t7 = t7 or {}
    for f in frames:
        f.t7 = t7.get(f.key)
    index = {f.key: i for i, f in enumerate(frames)}
    episodes, solved = [], 0
    for block in target_blocks(frames):
        episodes += tracking_episodes(target, block, index)
        solve_frames(block, processed, cached)
        found, n = drift_episodes(target, block, index)
        episodes += found
        solved += n
    episodes.sort(key=lambda e: (e['i0'], e['kind']))
    return {'frames': len(frames), 'tracking_known': any(f.tracking is not None for f in frames),
            'solved': solved, 'episodes': episodes}


def night(raw_dir, log_paths=(), processed=None, cache=None, frames=None):
    """The faults of one telescope-night: {target (OBJECT): target_faults()}. raw_dir is the night's raw frame
    directory; log_paths the night's pipeline logs, oldest first; processed {key: path} (processed_frames());
    cache the night's WCS cache table (cache_path()); frames the night's raw_frames(), if already read."""
    frames = raw_frames(raw_dir) if frames is None else frames
    t7 = t7_outcomes(log_paths) if log_paths else {}
    cached = wcs_cache(cache) if cache else {}
    by_target = {}
    for b, (target, run) in enumerate(blocks(frames)):
        for f in run:
            f.block = b
        by_target.setdefault(target, []).extend(run)
    return {target: target_faults(target, fs, processed, t7, cached) for target, fs in by_target.items()}


def episodes(result):
    """Every episode of a night() result, in time order."""
    return sorted((e for r in result.values() for e in r['episodes']), key=lambda e: (e['start'] or '', e['kind']))


# ---------------------------------------------------------------------------
# The night's record, for the morning digest
# ---------------------------------------------------------------------------

SUMMARY_VERSION = 1


def summary_path(datdir, night_):
    """Where T11 leaves a night's faults: <PipelineOutput>/v3/<TEL>/reports/telescope_faults/<night>.json. The name
    must not hold '<TEL>_<night>': T12 moves the first file under reports/ matching that, which must be the PDF."""
    return os.path.join(datdir, 'reports', 'telescope_faults', '{}.json'.format(night_))


def write_summary(path, result, telescope, night_, orchard=''):
    """Write a night() result as JSON (through a temporary file, so a reader never sees half of one)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tracking, drift = flags(result)
    doc = {'version': SUMMARY_VERSION, 'telescope': telescope, 'night': night_, 'orchard': orchard,
           'written': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
           'tracking_lost': tracking, 'pointing_drift': drift, 'targets': result}
    tmp = '{}.tmp{}'.format(path, os.getpid())
    try:
        with open(tmp, 'w') as f:
            json.dump(doc, f, indent=1, sort_keys=True)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_summary(path):
    """The record write_summary() left, or None if there is none (or it cannot be read)."""
    try:
        with open(path) as f:
            doc = json.load(f)
    except (OSError, ValueError):
        return None
    return doc if isinstance(doc, dict) and isinstance(doc.get('targets'), dict) else None


# ---------------------------------------------------------------------------
# Words
# ---------------------------------------------------------------------------

def _hm(iso):
    return iso[11:16] if iso and len(iso) >= 16 else '?'


def _span(e):
    return '{}–{} UT'.format(_hm(e['start']), _hm(e['end']))


def describe(e, frames_total=None):
    """One episode in a sentence for the night report."""
    d = e['detail']
    of = ' of {}'.format(frames_total) if frames_total else ''
    frames = 'frames {}–{}{} ({} frame{})'.format(e['i0'], e['i1'], of, e['n'], '' if e['n'] == 1 else 's')
    if e['kind'] == 'tracking':
        why = []
        if d.get('parked'):
            why.append('{} with the telescope parked'.format(d['parked']))
        if d.get('ra_rate') is not None and abs(d['ra_rate'] - SIDEREAL_DEG_PER_H) < 1.5:
            why.append('header RA ran {:.1f} deg/h'.format(d['ra_rate']) +
                       (' at constant airmass {:.3f}'.format(d['airmass']) if d.get('airmass_span', 1) < 1e-3
                        else ''))
        if d.get('altitude') is not None and d['altitude'] < LOW_ALTITUDE:
            why.append('at altitude {:.0f} deg'.format(d['altitude']))
        return '{}: tracking stopped {}, {}{}'.format(e['target'] or '?', _span(e), frames,
                                                     '; ' + ', '.join(why) if why else '')
    if e['kind'] == 'trailed':
        return '{}: stars trailed (T7) {}, {}'.format(e['target'] or '?', _span(e), frames)
    if d.get('jump_arcmin'):
        how = "jumped {:.1f}' at {} UT".format(d['jump_arcmin'], _hm(e['start']))
        if d['max_arcmin'] >= d['jump_arcmin'] + DRIFT_ARCMIN:
            how += ", up to {:.1f}' away by {} UT".format(d['max_arcmin'], _hm(d.get('max_at')))
    else:
        how = "drifted, up to {:.1f}' away by {} UT".format(d['max_arcmin'], _hm(d.get('max_at')))
    text = '{}: pointing {}; {}, {}'.format(e['target'] or '?', how, _span(e), frames)
    if d.get('off_field_at'):
        text += '; target off the field from {} UT'.format(_hm(d['off_field_at']))
    if d.get('unsolved_after'):
        text += '; no frame solved after {} UT'.format(_hm(d.get('last_solved_at')))
    return text


def flags(result):
    """('Y'|'N'|'unknown', 'Y'|'N'|'unknown'): tracking lost, pointing drift, for the night report's events table.
    Tracking is unknown without Astra's TRACKING (MaxIm DL frames) unless T7 found trailed frames; drift is unknown
    when no frame had a solution."""
    eps = episodes(result)
    tracking = 'Y' if any(e['kind'] in ('tracking', 'trailed') for e in eps) else \
        'N' if any(r['tracking_known'] for r in result.values()) else 'unknown'
    drift = 'Y' if any(e['kind'] == 'drift' for e in eps) else \
        'N' if any(r['solved'] >= DRIFT_MIN_SOLVED for r in result.values()) else 'unknown'
    return tracking, drift


def digest_text(result, limit=2):
    """The night's telescope faults for the digest's telescope line, '' if none."""
    eps = episodes(result)
    if not eps:
        return ''
    def words(e):
        if e['kind'] != 'drift':
            return {'tracking': 'tracking stopped', 'trailed': 'stars trailed'}[e['kind']]
        d = e['detail']
        return "pointing jumped {:.1f}'".format(d['jump_arcmin']) if d.get('jump_arcmin') else \
            "pointing drifted {:.1f}'".format(d['max_arcmin'])
    parts = ['{} on {} {} ({} frame{})'.format(words(e), e['target'] or '?', _span(e), e['n'],
                                               '' if e['n'] == 1 else 's') for e in eps[:limit]]
    if len(eps) > limit:
        parts.append('{} more'.format(len(eps) - limit))
    return 'telescope: ' + '; '.join(parts)
