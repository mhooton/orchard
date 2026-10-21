"""
Register of calibration frames.

One row for every bias, dark and flat frame under Observations/<Tel>/images/<night>/, with what T1
(main/createlists.py) needs to choose a night's calibrations without opening any other night's frames:
type, filter, exposure time, detector configuration (shape, binning, window, pixel size, orientation),
detector temperature, readout mode, and the level of the frame (median of the central rows above the
overscan, the fraction of saturated pixels).  calselect.py queries it for the best single night.

Why.  T1 used to look for missing calibrations one night at a time, backwards, re-reading every header of every
night it passed, for up to 3 years (flats) or 100 nights (bias, darks), and took what it found first - from any
camera, at any exposure time, and for flats adding nights together until it had four.  Callisto 20250407 (SPIRIT)
spent 2h20m in T1 and got its Clear master flat from ANDOR flats of 20230908; v2 SPIRIT nights 20220111 and
20220607 took ANDOR bias and darks.

Storage: one SQLite file (stdlib sqlite3), by default <data root>/calibration_register/register.sqlite.  The
register has ~600,000 rows, is queried by every T1 (several nights per telescope can run at once) and gains a
night's rows from T1, from the queue's download and fetch jobs and from the builder, so it needs indexed queries
and atomic updates from concurrent writers; a CSV would have to be parsed whole by every T1 and rewritten under a
lock for every night.  Like the queue's database it must sit on appct's local disk and be opened only from
orchard-server on appct (appcs sees that disk over NFS, where SQLite locking is unsafe); `export` writes a CSV
copy to read anywhere.

Frames are keyed by (telescope, night directory, stem), the path in the night directory without its FITS
suffixes, so a frame that later gains a .fz suffix keeps its row; calselect.py resolves the file on disk.

The register only describes frames.  Whether a frame is usable is decided by flag_night() from its statistics
and the other frames of its night (RULES_VERSION); `reflag` re-applies the rules without reading any frame.

Commands (python -m calibration.calregister ...; read-only on the frames, writes only the register):
    build      one-off scan of the archive (prefiltered, see candidate_paths), resumable
    eso-index  the ESO archive's calibration rows, which tell the builder which ESO-named frames to read
    night      (re-)register one night from its own headers: T1, download and fetch jobs do this
    update     register nights that are missing or whose directory changed (e.g. the last N days)
    reflag     re-apply the usability rules
    validate   compare the builder's prefilter with a full header scan of some nights
    summary / export / epochs
"""
import argparse
import contextlib
import csv
import datetime as dt
import functools
import json
import os
import re
import sqlite3
import sys
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np

REGISTER_VERSION = 1
RULES_VERSION = 1
REGISTER_ENV = 'CALIB_REGISTER'
READONLY_ENV = 'CALIB_REGISTER_READONLY'
DEFAULT_RELPATH = os.path.join('calibration_register', 'register.sqlite')

KINDS = ('bias', 'dark', 'flat')
# the suffixes T1 can be handed (find_ext: fits, fts, fz) plus the ones the archive also holds
FITS_SUFFIXES = ('.fits.fz', '.fts.fz', '.fit.fz', '.fits.bz2', '.fts.bz2', '.fits', '.fts', '.fit')
COMPRESSED_SUFFIXES = ('.fz',)
NIGHT_DIR = re.compile(r'^\d{8}$')
# an ESO-named frame: SPECUn.<time>, or SPECUn.<time>.<k> for frames unpacked from datacubes (2026-02/03)
DP_ID_NAME = re.compile(r'^SPECU\d+\.(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)(\.\d+)?(\.\d+)?$')
CALIBRATION_NAME = re.compile(r'^(bias|dark|flat|autoflat|skyflat|domeflat|twilight|dusk|dawn|zero|calib)', re.I)

# rows of the central band of a frame whose level is measured (each side of the middle row)
BAND_HALF = 32
# Raw level (ADU) at which pixels stop responding, for instruments whose instrument_config.json
# saturation_threshold is something else.  SPIRIT's pixels fill at 14,300-14,700 ADU raw (bias pedestal ~2,810;
# a zYJ flat with a median of 14,410 is flat-topped, the 14-bit ADC tops out at ~16,200; Callisto 2025-04-01,
# mh_scratch/calib_register/ceiling.py), so its configured 11,000 would call most of its twilight flats and
# every dark longer than 30 s saturated.  The Andor modes take readout.mode_saturation().
PIXEL_CEILING = {'spirit': 14400.0}
# usability rules (RULES_VERSION 1); levels in ADU above the overscan, or above the night's bias level
SATURATED_FRACTION = 0.05      # of the band's pixels at or above the ceiling (hot pixels and stars are fewer)
OUTLIER_MIN_ADU = 30.0         # a bias or dark further than this, and OUTLIER_SIGMA, from its night's median
OUTLIER_SIGMA = 6.0
FLAT_MIN_SIGNAL = 1000.0       # a flat below this is too faint to add anything
FLAT_MAX_FRACTION = 0.85       # of the raw ceiling: a flat whose median is above it is near saturation

SCHEMA = """
CREATE TABLE IF NOT EXISTS frames (
    telescope TEXT NOT NULL,     -- directory under Observations (Io, Europa, Ganymede, Callisto, Artemis, ...)
    night_dir TEXT NOT NULL,     -- the night's directory name
    night TEXT NOT NULL,         -- YYYYMMDD
    stem TEXT NOT NULL,          -- path in the night directory without FITS suffixes
    relpath TEXT NOT NULL,       -- path in the night directory as last seen
    kind TEXT NOT NULL,          -- bias | dark | flat
    imagetyp TEXT,
    filter TEXT,
    exptime REAL,
    date_obs TEXT,
    naxis1 INTEGER,
    naxis2 INTEGER,
    xbin INTEGER,
    ybin INTEGER,
    xorg INTEGER,
    yorg INTEGER,
    pixsize REAL,                -- unbinned pixel size in microns (XPIXSZ / XBINNING)
    orient TEXT,                 -- how the stored frame was turned from the camera's readout
    config TEXT,                 -- naxis1 x naxis2, binning, window, orientation: must match the science frames
    detector TEXT,               -- the camera: unbinned pixel size, else unbinned shape
    ccd_temp REAL,
    set_temp REAL,
    telescop TEXT,               -- configured telescope name (instrument_config.json)
    instrument TEXT,             -- andor | spirit | ...
    readmode TEXT,
    os_level REAL,
    os_noise REAL,
    level REAL,                  -- median of the central band minus the overscan level
    p999 REAL,                   -- 99.9th percentile of the central band, raw ADU
    sat_frac REAL,
    saturation REAL,             -- raw ADU
    zero REAL,                   -- the night's bias level for this configuration (flag_night)
    signal REAL,                 -- level - zero
    compressed INTEGER NOT NULL DEFAULT 0,
    usable INTEGER NOT NULL DEFAULT 1,
    reason TEXT,
    error TEXT,
    size INTEGER,
    mtime REAL,
    scanned_at TEXT,
    PRIMARY KEY (telescope, night_dir, stem)
);
CREATE INDEX IF NOT EXISTS frames_by_kind ON frames (telescope, kind, night);
CREATE TABLE IF NOT EXISTS nights (
    telescope TEXT NOT NULL,
    night_dir TEXT NOT NULL,
    night TEXT NOT NULL,
    n_files INTEGER,
    n_bias INTEGER,
    n_dark INTEGER,
    n_flat INTEGER,
    detector TEXT,               -- majority detector of the night's frames read (calibration + sampled science)
    science TEXT,                -- JSON: the science frames read (configs, filters, exposure times, temperatures)
    method TEXT,                 -- build | night | t1 | hook
    dir_mtime REAL,
    rules_version INTEGER,
    scanned_at TEXT,
    PRIMARY KEY (telescope, night_dir)
);
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
"""

FRAME_COLUMNS = ('telescope', 'night_dir', 'night', 'stem', 'relpath', 'kind', 'imagetyp', 'filter', 'exptime',
                 'date_obs', 'naxis1', 'naxis2', 'xbin', 'ybin', 'xorg', 'yorg', 'pixsize', 'orient', 'config',
                 'detector', 'ccd_temp', 'set_temp', 'telescop', 'instrument', 'readmode', 'os_level', 'os_noise',
                 'level', 'p999', 'sat_frac', 'saturation', 'zero', 'signal', 'compressed', 'usable', 'reason',
                 'error', 'size', 'mtime', 'scanned_at')
NIGHT_COLUMNS = ('telescope', 'night_dir', 'night', 'n_files', 'n_bias', 'n_dark', 'n_flat', 'detector', 'science',
                 'method', 'dir_mtime', 'rules_version', 'scanned_at')
NETWORK_FS = ('nfs', 'nfs4', 'cifs', 'smbfs', 'smb3', 'fuse.sshfs', 'afpfs')


def utcnow_iso():
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0, tzinfo=None).isoformat() + 'Z'


def night_of(night_dir):
    """YYYYMMDD of a night directory ('20251124_rotated' -> '20251124'), or None."""
    m = re.match(r'^(\d{8})', str(night_dir))
    return m.group(1) if m else None


def night_days(night):
    return dt.datetime.strptime(night_of(night), '%Y%m%d').toordinal()


def strip_suffix(name):
    for s in FITS_SUFFIXES:
        if name.endswith(s):
            return name[:-len(s)]
    return name


def is_fits_name(name):
    return name.endswith(FITS_SUFFIXES) and not os.path.basename(name).startswith('.')


def is_compressed(name):
    return name.endswith(COMPRESSED_SUFFIXES)


def default_register_path(basedir):
    return os.path.join(basedir, DEFAULT_RELPATH)


def basedir_from_imageroot(imageroot):
    """<root> of <root>/Observations/<Tel>/images/<night>."""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(imageroot)))))


def resolve_register_path(imageroot=None, basedir=None):
    if os.getenv(REGISTER_ENV):
        return os.environ[REGISTER_ENV]
    if basedir is None:
        basedir = basedir_from_imageroot(imageroot)
    return default_register_path(basedir)


# ------------------------------------------------------------------------------------- reading frames

@functools.lru_cache(maxsize=1)
def _instrument_config():
    from calibration.pipeutils import load_instrument_config
    return load_instrument_config()


@functools.lru_cache(maxsize=256)
def _resolve_telescope(telescop, telescope_dir):
    from calibration.pipeutils import resolve_telescope
    config = _instrument_config()
    name = resolve_telescope(telescop, config, None)
    if name not in config.get('telescopes', {}):
        name = config.get('telescope_directories', {}).get(telescope_dir, name)
    return name


@functools.lru_cache(maxsize=256)
def _instrument(telescop, naxis2):
    """(instrument name, its parameters) as pipeutils.detect_instrument picks them: by NAXIS2, else the first."""
    tel_cfg = _instrument_config().get('telescopes', {}).get(telescop)
    if not tel_cfg:
        return None, {}
    for name, params in tel_cfg.items():
        if params.get('naxis2') == naxis2:
            return name, params
    name = next(iter(tel_cfg))
    return name, tel_cfg[name]


def _image_hdu(hdul):
    """The HDU holding the image: the primary one, or the first extension of a compressed or empty-primary file."""
    hdu = hdul[0]
    if int(hdu.header.get('NAXIS', 0) or 0) == 0 and len(hdul) > 1:
        return hdul[1]
    return hdu


def _section(hdu, r0, r1, c0, c1):
    """Rows r0:r1, columns c0:c1 of the image, reading whole rows (one contiguous read: astropy reads a
    column-trimmed section one row at a time)."""
    try:
        rows = hdu.section[r0:r1, :]
    except Exception:
        rows = hdu.data[r0:r1, :]
    return np.asarray(rows[:, c0:c1], dtype=float)


def frame_kind(imagetyp):
    """T1's classification (createlists.classify_file): LIGHT, DARK, BIAS, FLAT anywhere in IMAGETYP."""
    t = str(imagetyp or '').upper()
    if 'LIGHT' in t:
        return 'light'
    if 'DARK' in t:
        return 'dark'
    if 'BIAS' in t:
        return 'bias'
    if 'FLAT' in t:
        return 'flat'
    return None


def _int(v, default=None):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def _float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def orientation(hdr):
    """How the stored frame was turned from the camera's readout: Astra's cards, else MaxIm's FLIPSTAT.

    MaxIm-era Andor frames say 'Rotate 90 CCW', Astra-era ones ASTRAROT = 90, SPIRIT 'Mirror' (MaxIm) or
    ASTRAMIR = 'Y' (Astra); Callisto's 'Flip/Rotate 90 CW' nights (2022-01-08, 01-17) are mirrored the other
    way.  The labels are kept apart, so the change of control software is an epoch break (it came with camera
    rotations on Io, Europa and Ganymede) and a frame left untransformed never matches transformed ones.
    """
    parts = []
    if 'ASTRAROT' in hdr:
        parts.append('rot{}'.format(_int(hdr.get('ASTRAROT'), hdr.get('ASTRAROT'))))
    if 'ASTRAMIR' in hdr:
        parts.append('mir{}'.format(str(hdr.get('ASTRAMIR')).strip()))
    if parts:
        return 'astra:' + ':'.join(parts)
    flip = str(hdr.get('FLIPSTAT', '') or '').strip()
    if flip:
        return 'maxim:' + flip
    return 'none'


def theta_y(hdr):
    """Position angle (deg, N through E) of the frame's +Y axis from its WCS, or None without one.

    Only used to date camera rotations (`rotations`), never to match frames."""
    try:
        if 'CD1_2' in hdr or 'CD2_2' in hdr:
            c12, c22 = float(hdr.get('CD1_2', 0.0)), float(hdr.get('CD2_2', 0.0))
        elif 'CDELT1' in hdr and 'CDELT2' in hdr:
            c12 = float(hdr['CDELT1']) * float(hdr.get('PC1_2', 0.0))
            c22 = float(hdr['CDELT2']) * float(hdr.get('PC2_2', 1.0))
        else:
            return None
    except (TypeError, ValueError):
        return None
    if c12 == 0.0 and c22 == 0.0:
        return None
    return round(float(np.degrees(np.arctan2(c12, c22))) % 360.0, 2)


def header_fields(hdr, telescope_dir):
    """Everything the register keeps from a header (no pixel data)."""
    naxis1, naxis2 = _int(hdr.get('NAXIS1')), _int(hdr.get('NAXIS2'))
    xbin, ybin = _int(hdr.get('XBINNING'), 1) or 1, _int(hdr.get('YBINNING'), 1) or 1
    xorg, yorg = _int(hdr.get('XORGSUBF'), 0) or 0, _int(hdr.get('YORGSUBF'), 0) or 0
    xpix = _float(hdr.get('XPIXSZ'))
    pixsize = round(xpix / xbin, 2) if xpix else None
    orient = orientation(hdr)
    window = '' if (xorg, yorg) == (0, 0) else ' org{},{}'.format(xorg, yorg)
    config = '{}x{} b{}x{}{} {}'.format(naxis1, naxis2, xbin, ybin, window, orient)
    if pixsize:
        detector = '{:.1f}um'.format(pixsize)
    else:
        detector = '{}x{}'.format((naxis1 or 0) * xbin, (naxis2 or 0) * ybin)
    telescop = _resolve_telescope(str(hdr.get('TELESCOP', '') or '').strip(), telescope_dir)
    instrument, _ = _instrument(telescop, naxis2)
    return {
        'imagetyp': str(hdr.get('IMAGETYP', '') or '').strip(),
        'kind': frame_kind(hdr.get('IMAGETYP')),
        'filter': str(hdr.get('FILTER', '') or '').replace("'", '').strip(),
        'exptime': _float(hdr.get('EXPTIME', hdr.get('EXPOSURE'))),
        'date_obs': str(hdr.get('DATE-OBS', '') or '').strip(),
        'naxis1': naxis1, 'naxis2': naxis2, 'xbin': xbin, 'ybin': ybin, 'xorg': xorg, 'yorg': yorg,
        'pixsize': pixsize, 'orient': orient, 'config': config, 'detector': detector,
        'ccd_temp': _float(hdr.get('CCD-TEMP')), 'set_temp': _float(hdr.get('SET-TEMP')),
        'telescop': telescop, 'instrument': instrument, 'theta': theta_y(hdr),
    }


def pixel_stats(hdu, fields, night):
    """Overscan signature, readout mode, level and saturation of a frame from two small reads.

    Only the overscan strip and BAND_HALF rows either side of the middle of the trimmed image are read, so a
    frame costs ~0.3 MB of I/O rather than the whole 8.6 MB.
    """
    from calibration import readout
    _, params = _instrument(fields['telescop'], fields['naxis2'])
    n1, n2 = fields['naxis1'], fields['naxis2']
    out = {'readmode': 'DEFAULT', 'os_level': None, 'os_noise': None}
    full_height = bool(params) and params.get('naxis2') == n2
    if readout.has_overscan(params) and full_height:
        cfg = params['overscan']
        r0 = cfg.get('top_row')
        r0 = n2 + r0 if r0 is not None and r0 < 0 else (r0 or 0)
        r1 = cfg.get('bottom_row') or n2
        strip = _section(hdu, r0, r1, cfg.get('left_col') or 0, cfg.get('right_col') or n1)
        level, noise = readout.strip_signature(strip)
        out['os_level'], out['os_noise'] = level, noise
        out['readmode'] = readout.classify(noise, level, fields['telescop'], night)
    trim = params.get('trim') if full_height else None
    if trim:
        t0, t1 = trim.get('top_row') or 0, trim.get('bottom_row') or n2
        c0, c1 = trim.get('left_col') or 0, trim.get('right_col') or n1
    else:
        t0, t1, c0, c1 = 0, n2, 0, n1
    mid = (t0 + t1) // 2
    band = _section(hdu, max(t0, mid - BAND_HALF), min(t1, mid + BAND_HALF), c0, c1)
    band = band[np.isfinite(band)]
    if band.size == 0:
        out['error'] = 'no finite pixels in the central band'
        return out
    if out['readmode'] != 'DEFAULT':
        saturation = readout.mode_saturation(fields['telescop'], out['readmode'], night)
    else:
        saturation = PIXEL_CEILING.get(fields.get('instrument')) or _float(params.get('saturation_threshold')) \
            or 62000.0
    out['saturation'] = saturation
    out['level'] = float(np.median(band)) - (out['os_level'] or 0.0)
    out['p999'] = float(np.percentile(band, 99.9))
    out['sat_frac'] = float(np.mean(band >= saturation))
    return out


def read_frame(path, telescope_dir, night, stats=True):
    """Register fields of one frame (header, and pixel statistics if it is a calibration frame, or whatever it
    is with stats='always').

    Never raises: a frame that cannot be read comes back with 'error' set (and kind None unless its header
    said otherwise), so a bad file is recorded rather than stopping a night.
    """
    from astropy.io import fits
    out = {'error': None}
    try:
        st = os.stat(path)
        out['size'], out['mtime'] = st.st_size, st.st_mtime
    except OSError as e:
        out['error'] = 'stat: {}'.format(e.strerror)
        out['kind'] = None
        return out
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if path.endswith('.bz2'):
                from calibration.pipeutils import open_fits_file
                ctx = open_fits_file(path)
            else:
                ctx = fits.open(path, memmap=False, lazy_load_hdus=True)
            with ctx as hdul:
                hdu = _image_hdu(hdul)
                hdr = hdu.header
                out.update(header_fields(hdr, telescope_dir))
                if int(hdr.get('NAXIS', 0) or 0) != 2:
                    out['error'] = 'not a 2-D image (NAXIS={})'.format(hdr.get('NAXIS'))
                elif stats == 'always' or (stats and out['kind'] in KINDS):
                    out.update(pixel_stats(hdu, out, night))
    except Exception as e:
        out['error'] = '{}: {}'.format(type(e).__name__, str(e)[:200])
        out.setdefault('kind', None)
    return out


def night_files(night_path):
    """Relative paths of the FITS files T1 looks at: the night directory and its subdirectories, one level."""
    out = []
    try:
        entries = sorted(os.listdir(night_path))
    except OSError:
        return out
    for name in entries:
        if name.startswith('.'):
            continue
        full = os.path.join(night_path, name)
        if is_fits_name(name):
            out.append(name)
        elif os.path.isdir(full):
            try:
                for sub in sorted(os.listdir(full)):
                    if is_fits_name(sub):
                        out.append(name + '/' + sub)
            except OSError:
                continue
    return out


# ------------------------------------------------------------------------------------- usability

def _robust_sigma(values):
    v = np.asarray(values, dtype=float)
    if v.size < 3:
        return None
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def flag_night(rows):
    """Set zero, signal, usable and reason on one night's rows (RULES_VERSION 1), in place.

    zero: the median level of the night's bias frames of the same configuration and readout mode (0 for a
    camera whose overscan already removed it and that took no bias).  A bias or dark is an outlier when it
    lies further than max(OUTLIER_MIN_ADU, OUTLIER_SIGMA robust sigma) from the median of its group (dark:
    same exposure time); a flat must have a signal of at least FLAT_MIN_SIGNAL
    in the central band and a raw median below FLAT_MAX_FRACTION x the ceiling.  Any frame with
    SATURATED_FRACTION of its band at the ceiling is unusable.
    """
    def group_key(r, *extra):
        return (r.get('config'), r.get('readmode')) + extra

    zeros = {}
    for r in rows:
        if r.get('kind') == 'bias' and not r.get('error') and r.get('level') is not None:
            zeros.setdefault(group_key(r), []).append(r['level'])
    zero_of = {k: float(np.median(v)) for k, v in zeros.items()}
    groups = {}
    for r in rows:
        r['usable'], r['reason'] = 1, None
        if r.get('error'):
            r['usable'], r['reason'] = 0, 'unreadable'
            continue
        if r.get('level') is None:
            r['usable'], r['reason'] = 0, 'no level'
            continue
        zero = zero_of.get(group_key(r))
        if zero is None and r.get('os_level') is not None:
            zero = 0.0
        r['zero'] = zero
        r['signal'] = None if zero is None else r['level'] - zero
        if (r.get('sat_frac') or 0.0) >= SATURATED_FRACTION:
            r['usable'], r['reason'] = 0, 'saturated'
            continue
        if r['kind'] == 'flat':
            sig = r['signal'] if r['signal'] is not None else r['level']
            raw = r['level'] + (r.get('os_level') or 0.0)
            if sig < FLAT_MIN_SIGNAL:
                r['usable'], r['reason'] = 0, 'too faint'
            elif raw >= FLAT_MAX_FRACTION * (r.get('saturation') or 62000.0):
                r['usable'], r['reason'] = 0, 'near saturation'
        elif r['kind'] == 'bias':
            groups.setdefault(('bias',) + group_key(r), []).append(r)
        elif r['kind'] == 'dark':
            groups.setdefault(('dark',) + group_key(r, exposure_key(r.get('exptime'))), []).append(r)
    for members in groups.values():
        levels = [m['level'] for m in members]
        if len(levels) < 3:
            continue
        med, sig = float(np.median(levels)), _robust_sigma(levels) or 0.0
        limit = max(OUTLIER_MIN_ADU, OUTLIER_SIGMA * sig)
        for m in members:
            if abs(m['level'] - med) > limit:
                m['usable'], m['reason'] = 0, 'level outlier'
    return rows


def exposure_key(exptime):
    """The exposure time a dark is matched on: pipered.py rounds the science frame's to whole seconds."""
    if exptime is None:
        return None
    return int(round(float(exptime)))


# ------------------------------------------------------------------------------------- storage

def mount_fstype(path, mounts_file='/proc/self/mounts'):
    path = os.path.realpath(path)
    best, fstype = '', None
    try:
        with open(mounts_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                mnt = parts[1].replace('\\040', ' ')
                if (path == mnt or path.startswith(mnt.rstrip('/') + '/')) and len(mnt) > len(best):
                    best, fstype = mnt, parts[2]
    except OSError:
        return None
    return fstype


class RegisterError(RuntimeError):
    pass


class Register:
    """The SQLite register.  create=False refuses to make a new file (T1 and the queue's hooks never do)."""

    def __init__(self, path, create=False, readonly=False, check_local=True):
        self.path = path
        self.readonly = readonly
        exists = os.path.exists(path)
        if not exists and not create:
            raise RegisterError('no calibration register at {}'.format(path))
        directory = os.path.dirname(os.path.abspath(path))
        if check_local:
            fstype = mount_fstype(directory)
            if fstype in NETWORK_FS:
                raise RegisterError('refusing to open the calibration register {} on a network filesystem ({}); '
                                    'open it from orchard-server on appct'.format(path, fstype))
        if readonly:
            self.conn = sqlite3.connect('file:{}?mode=ro'.format(os.path.abspath(path)), uri=True, timeout=60,
                                        check_same_thread=False)
        else:
            os.makedirs(directory, exist_ok=True)
            self.conn = sqlite3.connect(path, timeout=60, isolation_level=None, check_same_thread=False)
            self.conn.execute('PRAGMA busy_timeout=60000')
            self.conn.execute('PRAGMA journal_mode=WAL')
            self.conn.execute('PRAGMA synchronous=NORMAL')
            self.conn.executescript(SCHEMA)
            if self.meta_get('register_version') is None:
                self.meta_set('register_version', str(REGISTER_VERSION))
                self.meta_set('created_at', utcnow_iso())
        self.conn.row_factory = sqlite3.Row
        self._lock = threading.Lock()

    def close(self):
        self.conn.close()

    @contextlib.contextmanager
    def tx(self):
        with self._lock:
            self.conn.execute('BEGIN IMMEDIATE')
            try:
                yield self.conn
            except BaseException:
                self.conn.execute('ROLLBACK')
                raise
            self.conn.execute('COMMIT')

    def meta_get(self, key):
        row = self.conn.execute('SELECT value FROM meta WHERE key = ?', (key,)).fetchone()
        return row[0] if row else None

    def meta_set(self, key, value):
        self.conn.execute('INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)', (key, value))

    def put_night(self, night_row, frame_rows):
        """Replace everything registered for (telescope, night_dir) with these rows, atomically."""
        if self.readonly:
            raise RegisterError('register opened read-only')
        key = (night_row['telescope'], night_row['night_dir'])
        with self.tx() as c:
            c.execute('DELETE FROM frames WHERE telescope = ? AND night_dir = ?', key)
            c.executemany('INSERT INTO frames ({}) VALUES ({})'.format(
                ', '.join(FRAME_COLUMNS), ', '.join('?' * len(FRAME_COLUMNS))),
                [tuple(r.get(k) for k in FRAME_COLUMNS) for r in frame_rows])
            c.execute('INSERT OR REPLACE INTO nights ({}) VALUES ({})'.format(
                ', '.join(NIGHT_COLUMNS), ', '.join('?' * len(NIGHT_COLUMNS))),
                tuple(night_row.get(k) for k in NIGHT_COLUMNS))

    def night(self, telescope, night_dir):
        row = self.conn.execute('SELECT * FROM nights WHERE telescope = ? AND night_dir = ?',
                                (telescope, night_dir)).fetchone()
        return dict(row) if row else None

    def nights(self, telescope=None):
        q, args = 'SELECT * FROM nights', ()
        if telescope:
            q, args = q + ' WHERE telescope = ?', (telescope,)
        return [dict(r) for r in self.conn.execute(q + ' ORDER BY telescope, night, night_dir', args)]

    def frames(self, telescope, kind=None, night_dir=None, usable_only=False):
        q, args = 'SELECT * FROM frames WHERE telescope = ?', [telescope]
        if kind:
            q, args = q + ' AND kind = ?', args + [kind]
        if night_dir:
            q, args = q + ' AND night_dir = ?', args + [night_dir]
        if usable_only:
            q += ' AND usable = 1'
        return [dict(r) for r in self.conn.execute(q + ' ORDER BY night, night_dir, stem', args)]

    def telescopes(self):
        return [r[0] for r in self.conn.execute('SELECT DISTINCT telescope FROM nights ORDER BY telescope')]


def open_existing(path, readonly=False):
    """The register for T1 and the queue's hooks: an existing file, read-write unless asked otherwise."""
    readonly = readonly or os.getenv(READONLY_ENV, '') not in ('', '0')
    return Register(path, create=False, readonly=readonly)


# ------------------------------------------------------------------------------------- registering a night

def summarise_science(science):
    """JSON-able summary of the science frames read for a night: per config, filters, exposures, temperatures."""
    out = {}
    for s in science:
        if s.get('error') or not s.get('config'):
            continue
        c = out.setdefault(s['config'], {'n': 0, 'detector': s.get('detector'), 'filters': {}, 'exptimes': {},
                                         'ccd_temp': [], 'theta': []})
        if s.get('theta') is not None:
            c['theta'].append(s['theta'])
        c['n'] += 1
        f = s.get('filter') or ''
        c['filters'][f] = c['filters'].get(f, 0) + 1
        e = exposure_key(s.get('exptime'))
        if e is not None:
            c['exptimes'][str(e)] = c['exptimes'].get(str(e), 0) + 1
        if s.get('ccd_temp') is not None:
            c['ccd_temp'].append(s['ccd_temp'])
    for c in out.values():
        temps = c.pop('ccd_temp')
        c['ccd_temp'] = round(float(np.median(temps)), 2) if temps else None
        th = c.pop('theta')
        if th:
            ref = th[0]
            th = [ref + ((t - ref + 180.0) % 360.0 - 180.0) for t in th]
            c['theta'] = round(float(np.median(th)) % 360.0, 2)
        else:
            c['theta'] = None
    return out


def majority_detector(rows):
    counts = {}
    for r in rows:
        if r.get('detector') and not r.get('error'):
            counts[r['detector']] = counts.get(r['detector'], 0) + 1
    if not counts:
        return None
    return max(sorted(counts), key=lambda k: counts[k])


_WORKER_THROTTLE = None


def _init_worker(rate):
    """Each builder process paces its own file opens (its share of the build's rate)."""
    global _WORKER_THROTTLE
    _WORKER_THROTTLE = Throttle(rate)


def _scan(args):
    path, telescope, night, stats = args
    if _WORKER_THROTTLE is not None:
        _WORKER_THROTTLE()
    return read_frame(path, telescope, night, stats=stats)


def build_night_rows(telescope, night_dir, night_path, relpaths, science_relpaths=(), method='night', workers=4,
                     throttle=None, n_files=None, executor=None, preread=None, log=print):
    """(night row, frame rows) for one night: reads the given candidate files (and science samples).

    Frames are read by `workers` threads, or by the builder's process pool (`executor`): the scan is mostly
    Python (headers, sigma-clipped overscan statistics), so threads cannot run more than ~15 frames a second.
    preread: (infos of relpaths, infos of science_relpaths) the builder already read in a batch of nights.
    """
    night = night_of(night_dir)
    scanned_at = utcnow_iso()

    def finish(rel, info):
        info.update({'telescope': telescope, 'night_dir': night_dir, 'night': night, 'relpath': rel,
                     'stem': strip_suffix(rel), 'compressed': int(is_compressed(rel)), 'scanned_at': scanned_at})
        return info

    def read_all(pool, rels, stats=True):
        rels = list(rels)
        if executor is not None:
            infos = executor.map(_scan, [(os.path.join(night_path, r), telescope, night, stats) for r in rels],
                                 chunksize=4)
        else:
            def one(rel):
                if throttle is not None:
                    throttle()
                return read_frame(os.path.join(night_path, rel), telescope, night, stats=stats)
            infos = pool.map(one, rels)
        return [finish(r, i) for r, i in zip(rels, infos)]

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        if preread is not None:
            infos = [finish(r, i) for r, i in zip(relpaths, preread[0])]
            sci = [finish(r, i) for r, i in zip(science_relpaths, preread[1])]
        else:
            infos = read_all(pool, relpaths)
            sci = read_all(pool, science_relpaths, stats=False)
        missed = [s['relpath'] for s in sci if s.get('kind') in KINDS]
        if missed:
            log('{} {}: {} sampled frame(s) the prefilter skipped are calibration frames, e.g. {}'.format(
                telescope, night_dir, len(missed), missed[0]))
            infos += read_all(pool, missed)
    cal = [i for i in infos if i.get('kind') in KINDS]
    science = [i for i in infos if i.get('kind') == 'light'] + [s for s in sci if s.get('kind') == 'light']
    # one row per stem: a frame stored twice (x.fits and x.fits.fz) keeps the uncompressed copy
    by_stem = {}
    for r in sorted(cal, key=lambda r: (r['stem'], r['compressed'])):
        by_stem.setdefault(r['stem'], r)
    cal = list(by_stem.values())
    flag_night(cal)
    try:
        dir_mtime = os.stat(night_path).st_mtime
    except OSError:
        dir_mtime = None
    night_row = {
        'telescope': telescope, 'night_dir': night_dir, 'night': night,
        'n_files': n_files if n_files is not None else len(relpaths),
        'n_bias': sum(1 for r in cal if r['kind'] == 'bias'),
        'n_dark': sum(1 for r in cal if r['kind'] == 'dark'),
        'n_flat': sum(1 for r in cal if r['kind'] == 'flat'),
        'detector': majority_detector(cal + science),
        'science': json.dumps(summarise_science(science), sort_keys=True),
        'method': method, 'dir_mtime': dir_mtime, 'rules_version': RULES_VERSION, 'scanned_at': scanned_at,
    }
    return night_row, cal


def science_samples(relpaths, cal_relpaths, n=3):
    """First, middle and last of the files that are not calibration candidates (a night's camera)."""
    cal = set(cal_relpaths)
    rest = [r for r in relpaths if r not in cal]
    if len(rest) <= n:
        return rest
    idx = sorted({0, len(rest) // 2, len(rest) - 1})
    return [rest[i] for i in idx]


def register_night(register, basedir, telescope, night_dir, method='night', workers=4, log=print):
    """Read every header of a night and register its calibration frames (T1, download/fetch hooks, `night`)."""
    night_path = os.path.join(basedir, 'Observations', telescope, 'images', night_dir)
    relpaths = night_files(night_path)
    night_row, rows = build_night_rows(telescope, night_dir, night_path, relpaths, method=method,
                                       workers=workers, log=log)
    register.put_night(night_row, rows)
    return night_row, rows


def register_rows_from_t1(register, telescope, night_dir, night_path, cal_relpaths, science_infos, n_files,
                          workers=4):
    """T1's own night: the calibration frames it classified (pixel statistics read here) and its science frames
    (already read by T1).  Returns (night row, rows) and writes them unless the register is read-only."""
    night_row, rows = build_night_rows(telescope, night_dir, night_path, cal_relpaths, method='t1',
                                       workers=workers, n_files=n_files)
    night_row['science'] = json.dumps(summarise_science(science_infos), sort_keys=True)
    night_row['detector'] = majority_detector(rows + list(science_infos))
    if register is not None and not register.readonly:
        register.put_night(night_row, rows)
    return night_row, rows


def hook_register_night(basedir, telescope, night_dir, register_path=None, log=print):
    """Best-effort registration of a night that just gained frames (the queue's download and fetch jobs).

    Only updates an existing register and never raises: the register is a convenience for T1 (which registers
    its own night again anyway), never a reason for a download to fail.  The path comes from the queue's
    configuration (or the data root), never from $CALIB_REGISTER, which is for T1 replays.
    """
    path = register_path or default_register_path(basedir)
    if not os.path.exists(path):
        log('calibration register: none at {}; not registering {} {}'.format(path, telescope, night_dir))
        return None
    try:
        reg = open_existing(path)
        if reg.readonly:
            reg.close()
            return None
        try:
            row, rows = register_night(reg, basedir, telescope, night_dir, method='hook', workers=2, log=log)
        finally:
            reg.close()
        log('calibration register: {} {}: {} bias, {} dark, {} flat frames registered'.format(
            telescope, night_dir, row['n_bias'], row['n_dark'], row['n_flat']))
        return row
    except Exception as e:  # never fail the job over the register
        log('calibration register: could not register {} {}: {}'.format(telescope, night_dir, e))
        return None


# ------------------------------------------------------------------------------------- the one-off builder

class Throttle:
    """At most `rate` file opens per second across threads (0 = unlimited)."""

    def __init__(self, rate):
        self.interval = 1.0 / rate if rate and rate > 0 else 0.0
        self.lock = threading.Lock()
        self.next = time.monotonic()

    def __call__(self):
        if not self.interval:
            return
        with self.lock:
            now = time.monotonic()
            wait = self.next - now
            self.next = max(self.next, now) + self.interval
        if wait > 0:
            time.sleep(wait)


CUBE_ORIGFILE = re.compile(r'^SPECU\d+\.(\d{8}T\d{6})_')


def _mjd(t):
    return (t - dt.datetime(1858, 11, 17)).total_seconds() / 86400.0


def dp_key(name):
    """An ESO-named frame's identifier with three-digit milliseconds, or None.  Frames unpacked from the
    2026-02/03 datacubes are named <time>.<milliseconds as an unpadded integer>.<k>: SPECU4.2026-02-04T23:48:53.95.0
    is ESO's SPECU4.2026-02-04T23:48:53.095; other names carry a decimal fraction (.932 = 932 ms)."""
    m = DP_ID_NAME.match(name)
    if not m:
        return None
    if not m.group(2):
        ms = 0
    elif m.group(3):
        ms = int(m.group(2)[1:])
    else:
        ms = int(round(float('0' + m.group(2)) * 1000))
    return '{}.{}.{:03d}'.format(name.split('.', 1)[0], m.group(1), ms)


def eso_nights(path):
    """{(telescope, night)} with any row in the ESO index (a night ESO knows nothing about is read whole)."""
    out = set()
    if path and os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                out.add((r['telescope'], r['night']))
    return out


def load_eso_index(path):
    """({dp_id: dp_type}, {(telescope, night): [(cube start mjd, dp_type)]}) from eso-index.

    A row whose origfile is SPECUn.<YYYYMMDDTHHMMSS>_... is a datacube that started then (the dp_id and
    mjd_obs are later); its frames are unpacked under their own times (jobqueue.nights.row_files).
    """
    by_id, cubes = {}, {}
    if not path or not os.path.exists(path):
        return by_id, cubes
    with open(path) as f:
        for r in csv.DictReader(f):
            by_id[dp_key(r['dp_id']) or r['dp_id']] = r['dp_type']
            m = CUBE_ORIGFILE.match(r.get('origfile') or '')
            if m:
                start = _mjd(dt.datetime.strptime(m.group(1), '%Y%m%dT%H%M%S'))
                cubes.setdefault((r['telescope'], r['night']), []).append((start, r['dp_type']))
    for v in cubes.values():
        v.sort()
    return by_id, cubes


def _name_mjd(base):
    """MJD of an ESO-named frame from its name, or None."""
    key = dp_key(base)
    if key is None:
        return None
    t = dt.datetime.strptime(key.split('.', 1)[1], '%Y-%m-%dT%H:%M:%S.%f')
    return _mjd(t)


def candidate_paths(telescope, night_dir, relpaths, eso_types, cubes, tolerance_s=1.0, known_nights=None):
    """Which files of a night the builder reads, without reading any header:

    - an ESO-named file (SPECUn.<time>[.<k>]) that ESO lists as BIAS, DARK or FLAT; on a night delivered as
      datacubes (SPIRIT from 2026-02) also any frame not inside an OBJECT cube's time span, which is how the
      frames unpacked from calibration cubes (named by their own time) are found;
    - every ESO-named file of a night the ESO index has no row for at all (`known_nights`): some nights on disk
      are not in ESO's archive (Ganymede 2024-12-19/20, ...), so their frames cannot be typed from it;
    - any other file in a subdirectory (AutoFlat/, Calibration/, ...) or named like a calibration frame
      (Bias-, Dark-, Flat-, AutoFlat-, Dusk-, ...).  Top-level ACP frames <target>-S001-R001-C001-<f>.fts are
      lights: in the 2026-10-16 header census no top-level non-ESO Io/Europa/Ganymede/Callisto/Artemis frame
      with such a name was a calibration frame (Saint-Ex and Moana-ES name theirs differently: read whole);
    - every file of a telescope without an ESO index (Saint-Ex, Moana-ES).
    `validate` measures what this misses against full header scans.
    """
    import bisect
    if telescope not in ('Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis'):
        return list(relpaths)
    out = []
    cube_list = cubes.get((telescope, night_of(night_dir)))
    starts = [c[0] for c in cube_list] if cube_list else None
    unknown_night = known_nights is not None and (telescope, night_of(night_dir)) not in known_nights
    for rel in relpaths:
        stem = strip_suffix(rel)
        base = stem.split('/')[-1]
        m = DP_ID_NAME.match(base)
        if m and unknown_night:
            out.append(rel)
            continue
        if m:
            types = [eso_types[i] for i in (dp_key(base),) if i in eso_types]
            if types:
                if types[0] in ('BIAS', 'DARK', 'FLAT'):
                    out.append(rel)
                continue
            if starts:
                # a frame unpacked from a cube: the latest cube that started by its time (+1 s)
                mjd = _name_mjd(base)
                i = bisect.bisect_right(starts, mjd + tolerance_s / 86400.0) - 1 if mjd is not None else -1
                if i < 0 or cube_list[i][1] != 'OBJECT':
                    out.append(rel)
            continue
        if '/' in rel or CALIBRATION_NAME.match(base):
            out.append(rel)
    return out


def build(register, basedir, telescopes, eso_index=None, since=None, until=None, workers=4, rate=0.0,
          rescan=False, nights=None, log=print):
    """Register every night of the telescopes (resumable: a night already registered by the builder with the
    same number of files is skipped unless rescan)."""
    from concurrent.futures import ProcessPoolExecutor
    eso_types, cubes = load_eso_index(eso_index)
    known = eso_nights(eso_index) if eso_index else None
    workers = max(1, workers)
    executor = ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                   initargs=((rate / workers) if rate else 0.0,))
    totals = {'nights': 0, 'skipped': 0, 'files_read': 0, 'cal': 0}
    t_start = time.time()
    try:
        _build(register, basedir, telescopes, eso_types, cubes, executor, since, until, rescan, nights, totals,
               t_start, log, known=known)
    finally:
        executor.shutdown()
    log('done: {}'.format(totals))
    return totals


def _build(register, basedir, telescopes, eso_types, cubes, executor, since, until, rescan, nights, totals,
           t_start, log, batch_files=400, known=None):
    batch = []

    def flush():
        # read the frames of several nights in one go, so no process waits for one night's slowest frame
        tasks = []
        for tel, nd, path, relpaths, cands, samples in batch:
            night = night_of(nd)
            tasks += [(os.path.join(path, r), tel, night, True) for r in cands]
            tasks += [(os.path.join(path, r), tel, night, False) for r in samples]
        infos = list(executor.map(_scan, tasks, chunksize=2)) if tasks else []
        i = 0
        for tel, nd, path, relpaths, cands, samples in batch:
            got_cal = infos[i:i + len(cands)]
            got_sci = infos[i + len(cands):i + len(cands) + len(samples)]
            i += len(cands) + len(samples)
            night_row, rows = build_night_rows(tel, nd, path, cands, samples, method='build',
                                               n_files=len(relpaths), executor=executor,
                                               preread=(got_cal, got_sci), log=log)
            register.put_night(night_row, rows)
            totals['nights'] += 1
            totals['files_read'] += len(cands) + len(samples)
            totals['cal'] += len(rows)
            if totals['nights'] % 50 == 0:
                el = time.time() - t_start
                log('{} {}: {} nights, {} files read ({:.0f}/s), {} calibration frames'.format(
                    tel, nd, totals['nights'], totals['files_read'],
                    totals['files_read'] / max(el, 1e-6), totals['cal']))
        del batch[:]

    for telescope in telescopes:
        images = os.path.join(basedir, 'Observations', telescope, 'images')
        try:
            dirs = sorted(d for d in os.listdir(images) if NIGHT_DIR.match(d))
        except OSError as e:
            log('{}: cannot list {}: {}'.format(telescope, images, e))
            continue
        if nights:
            dirs = [d for d in dirs if d in nights]
        if since:
            dirs = [d for d in dirs if d >= since]
        if until:
            dirs = [d for d in dirs if d <= until]
        done = {} if rescan else {n['night_dir']: n for n in register.nights(telescope)}
        log('{}: {} night directories'.format(telescope, len(dirs)))
        for night_dir in dirs:
            path = os.path.join(images, night_dir)
            relpaths = night_files(path)
            prev = done.get(night_dir)
            if prev and prev.get('n_files') == len(relpaths):
                totals['skipped'] += 1
                continue
            cands = candidate_paths(telescope, night_dir, relpaths, eso_types, cubes, known_nights=known)
            samples = science_samples(relpaths, cands)
            batch.append((telescope, night_dir, path, relpaths, cands, samples))
            if sum(len(b[4]) + len(b[5]) for b in batch) >= batch_files:
                flush()
    flush()


def reflag(register, telescope=None, log=print):
    n = 0
    for nrow in register.nights(telescope):
        rows = register.frames(nrow['telescope'], night_dir=nrow['night_dir'])
        if not rows:
            continue
        flag_night(rows)
        nrow['rules_version'] = RULES_VERSION
        register.put_night(nrow, rows)
        n += 1
    log('re-flagged {} nights'.format(n))
    return n


def update(register, basedir, telescopes, days=None, since=None, workers=4, log=print):
    """Register nights that are missing, or whose file count changed, in the last `days` days (or since)."""
    if days is not None:
        since = (dt.date.today() - dt.timedelta(days=int(days))).strftime('%Y%m%d')
    n = 0
    for telescope in telescopes:
        images = os.path.join(basedir, 'Observations', telescope, 'images')
        try:
            dirs = sorted(d for d in os.listdir(images) if NIGHT_DIR.match(d) and (not since or d >= since))
        except OSError:
            continue
        known = {r['night_dir']: r for r in register.nights(telescope)}
        for d in dirs:
            relpaths = night_files(os.path.join(images, d))
            k = known.get(d)
            if k and k.get('n_files') == len(relpaths):
                continue
            row, _ = register_night(register, basedir, telescope, d, method='night', workers=workers, log=log)
            log('{} {}: {} bias, {} dark, {} flat'.format(telescope, d, row['n_bias'], row['n_dark'], row['n_flat']))
            n += 1
    return n


# ------------------------------------------------------------------------------------- ESO index

ESO_PROGS = {'Io': '60.A-9009(A)', 'Europa': '60.A-9009(B)', 'Ganymede': '60.A-9009(C)', 'Callisto': '60.A-9009(D)'}
CUBE_ERA = {'Callisto': '20260201'}


def eso_index(out_path, env_path, telescopes=None, first_year=2017, log=print):
    """Write the ESO archive's calibration rows (and every row of cube-era nights) to a CSV for the builder."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from jobqueue.eso import EsoArchive
    from jobqueue.util import mjd_to_night, night_str, to_mjd
    archive = EsoArchive(env_path, log=log)
    fields = ['telescope', 'night', 'dp_id', 'dp_type', 'mjd_obs', 'exposure', 'filter_path', 'origfile', 'cube']
    tmp = out_path + '.tmp'
    n = 0
    with open(tmp, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for tel, prog in ESO_PROGS.items():
            if telescopes and tel not in telescopes:
                continue
            for year in range(first_year, dt.date.today().year + 1):
                for month in (1, 7):
                    a = dt.datetime(year, month, 1)
                    b = dt.datetime(year + (month == 7), 7 if month == 1 else 1, 1)
                    q = ("SELECT dp_id, dp_type, mjd_obs, exposure, filter_path, origfile FROM dbo.raw "
                         "WHERE prog_id = '{}' AND dp_id LIKE 'SPECU%' AND dp_type <> 'OBJECT' "
                         "AND mjd_obs >= {:.6f} AND mjd_obs < {:.6f}").format(prog, to_mjd(a), to_mjd(b))
                    rows = archive.tap(q)
                    if len(rows) >= 200000:
                        raise RuntimeError('ESO query hit MAXREC for {} {:%Y-%m}'.format(tel, a))
                    for r in rows:
                        r['telescope'] = tel
                        r['night'] = night_str(mjd_to_night(float(r['mjd_obs'])))
                        r['cube'] = ''
                        w.writerow({k: r.get(k, '') for k in fields})
                    n += len(rows)
                    log('{} {:%Y-%m}: {} calibration rows'.format(tel, a, len(rows)))
            if tel in CUBE_ERA:
                start = CUBE_ERA[tel]
                q = ("SELECT dp_id, dp_type, mjd_obs, exposure, filter_path, origfile FROM dbo.raw "
                     "WHERE prog_id = '{}' AND dp_id LIKE 'SPECU%' AND mjd_obs >= {:.6f}").format(
                    prog, to_mjd(dt.datetime.strptime(start, '%Y%m%d')))
                rows = archive.tap(q)
                if len(rows) >= 200000:
                    raise RuntimeError('ESO query hit MAXREC for the {} cube era'.format(tel))
                n_cubes = 0
                for r in rows:
                    # every row of the cube era, OBJECT included, so the builder can tell which frames came
                    # from which cube
                    r['telescope'] = tel
                    r['night'] = night_str(mjd_to_night(float(r['mjd_obs'])))
                    r['cube'] = '1' if CUBE_ORIGFILE.match(r.get('origfile') or '') else ''
                    n_cubes += r['cube'] == '1'
                    w.writerow({k: r.get(k, '') for k in fields})
                    n += 1
                log('{} cube era from {}: {} rows, {} of them cubes'.format(tel, start, len(rows), n_cubes))
    os.replace(tmp, out_path)
    log('wrote {} rows to {}'.format(n, out_path))
    return n


def validate(basedir, telescope, night_dirs, eso_index_path, workers=4, rate=0.0, log=print):
    """Full header scan of the nights against the prefilter: calibration frames the builder would not read."""
    eso_types, cubes = load_eso_index(eso_index_path)
    known = eso_nights(eso_index_path) if eso_index_path else None
    throttle = Throttle(rate)
    report = []
    for night_dir in night_dirs:
        path = os.path.join(basedir, 'Observations', telescope, 'images', night_dir)
        relpaths = night_files(path)
        cands = set(candidate_paths(telescope, night_dir, relpaths, eso_types, cubes, known_nights=known))

        def kind_of(rel):
            throttle()
            return rel, read_frame(os.path.join(path, rel), telescope, night_of(night_dir), stats=False).get('kind')

        with ThreadPoolExecutor(max_workers=workers) as pool:
            kinds = dict(pool.map(kind_of, relpaths))
        cal = {r for r, k in kinds.items() if k in KINDS}
        missed = sorted(cal - cands)
        extra = len(cands - cal)
        report.append({'telescope': telescope, 'night': night_dir, 'files': len(relpaths), 'cal': len(cal),
                       'candidates': len(cands), 'missed': len(missed), 'not_cal_read': extra,
                       'missed_examples': missed[:5]})
        log(json.dumps(report[-1]))
    return report


# ------------------------------------------------------------------------------------- CLI

def summary(register, log=print):
    for row in register.conn.execute(
            'SELECT telescope, kind, COUNT(*) AS n, SUM(usable) AS ok, COUNT(DISTINCT night_dir) AS nights '
            'FROM frames GROUP BY telescope, kind ORDER BY telescope, kind'):
        log('{:9s} {:5s} {:8d} frames, {:8d} usable, {:5d} nights'.format(*tuple(row)))
    for row in register.conn.execute(
            'SELECT telescope, reason, COUNT(*) FROM frames WHERE usable = 0 GROUP BY telescope, reason '
            'ORDER BY telescope, 3 DESC'):
        log('   unusable {:9s} {:20s} {}'.format(row[0], str(row[1]), row[2]))


def rotations(register, telescope, jump=1.5, log=print):
    """Nights where the camera's sky orientation (theta_Y of the sampled science frames' WCS) moves by more than
    `jump` degrees from the previous night with one: candidates for calib_epoch_breaks.csv."""
    prev = None
    out = []
    for n in register.nights(telescope):
        try:
            sci = json.loads(n.get('science') or '{}')
        except ValueError:
            continue
        thetas = [c.get('theta') for c in sci.values() if c.get('theta') is not None]
        if not thetas:
            continue
        theta = thetas[0]
        if prev is not None:
            d = (theta - prev[1] + 180.0) % 360.0 - 180.0
            if abs(d) > jump:
                out.append((prev[0], n['night_dir'], prev[1], theta))
                log('{} {} -> {}: theta_Y {:.1f} -> {:.1f} deg'.format(telescope, prev[0], n['night_dir'],
                                                                      prev[1], theta))
        prev = (n['night_dir'], theta)
    return out


def export(register, out_path):
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(FRAME_COLUMNS)
        for row in register.conn.execute('SELECT {} FROM frames ORDER BY telescope, night, stem'.format(
                ', '.join(FRAME_COLUMNS))):
            w.writerow(row)


def main(argv=None):
    p = argparse.ArgumentParser(prog='python -m calibration.calregister', description=__doc__.split('\n')[1])
    p.add_argument('--register', help='register file (default: $CALIB_REGISTER or <basedir>/{})'.format(
        DEFAULT_RELPATH))
    p.add_argument('--basedir', default='/data/SPECULOOSPipeline')
    sub = p.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build')
    b.add_argument('--telescopes', default='Io,Europa,Ganymede,Callisto,Artemis')
    b.add_argument('--eso-index')
    b.add_argument('--since')
    b.add_argument('--until')
    b.add_argument('--nights', help='comma-separated night directories')
    b.add_argument('--workers', type=int, default=4)
    b.add_argument('--rate', type=float, default=0.0, help='file opens per second (0: unlimited)')
    b.add_argument('--rescan', action='store_true')
    e = sub.add_parser('eso-index')
    e.add_argument('--out', required=True)
    e.add_argument('--env', default='/opt/orchard/src/reporting/.env')
    e.add_argument('--telescopes', default='Io,Europa,Ganymede,Callisto')
    n = sub.add_parser('night')
    n.add_argument('telescope')
    n.add_argument('night_dir')
    u = sub.add_parser('update')
    u.add_argument('--telescopes', default='Io,Europa,Ganymede,Callisto,Artemis')
    u.add_argument('--days', type=int)
    u.add_argument('--since')
    sub.add_parser('reflag')
    v = sub.add_parser('validate')
    v.add_argument('telescope')
    v.add_argument('nights', help='comma-separated night directories')
    v.add_argument('--eso-index')
    v.add_argument('--rate', type=float, default=0.0)
    v.add_argument('--out')
    sub.add_parser('summary')
    ro = sub.add_parser('rotations')
    ro.add_argument('telescope')
    ro.add_argument('--jump', type=float, default=1.5)
    x = sub.add_parser('export')
    x.add_argument('out')
    a = p.parse_args(argv)
    path = a.register or os.getenv(REGISTER_ENV) or default_register_path(a.basedir)
    log = functools.partial(print, flush=True)
    if a.cmd == 'eso-index':
        eso_index(a.out, a.env, telescopes=a.telescopes.split(','), log=log)
        return 0
    if a.cmd == 'validate':
        rep = validate(a.basedir, a.telescope, a.nights.split(','), a.eso_index, rate=a.rate, log=log)
        if a.out:
            with open(a.out, 'a') as f:
                for r in rep:
                    f.write(json.dumps(r) + '\n')
        return 0
    reg = Register(path, create=a.cmd == 'build', readonly=a.cmd in ('summary', 'export', 'rotations'))
    try:
        if a.cmd == 'build':
            build(reg, a.basedir, a.telescopes.split(','), eso_index=a.eso_index, since=a.since, until=a.until,
                  workers=a.workers, rate=a.rate, rescan=a.rescan,
                  nights=set(a.nights.split(',')) if a.nights else None, log=log)
        elif a.cmd == 'night':
            row, _ = register_night(reg, a.basedir, a.telescope, a.night_dir, log=log)
            log(json.dumps(row))
        elif a.cmd == 'update':
            update(reg, a.basedir, a.telescopes.split(','), days=a.days, since=a.since, log=log)
        elif a.cmd == 'reflag':
            reflag(reg, log=log)
        elif a.cmd == 'summary':
            summary(reg, log=log)
        elif a.cmd == 'export':
            export(reg, a.out)
        elif a.cmd == 'rotations':
            rotations(reg, a.telescope, a.jump, log=log)
    finally:
        reg.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
