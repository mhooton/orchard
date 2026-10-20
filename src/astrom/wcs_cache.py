"""
T7's plate-solve cache: one ECSV table per night, outside the raw and product trees (docs/wcs-cache.md).

    <BASEDIR>/PipelineOutput/wcs_cache/<TEL>/<YYYYMMDD>.ecsv     (ZLP passes it as astrometry.py --wcs-cache)

One row per raw frame, keyed by the raw file's name. A row holds the WCS cards a successful solve wrote into the
frame's calibrated copy, exactly as they stand in that file, so applying them to a later run's calibrated copy gives
the header a fresh solve would have. A row is used only when every validity field still matches (the raw file, the
calibrated frame's shape and trim, the solver, the Gaia database, astropy and twirl) and it was written by a code
release or by the commit that is running now: a scratch tree pointed at the production data root can then never
feed a production run.

The table sits outside v2/v3, so reruns, T12 and promotion never move it, and a scratch data root has its own.
Writers hold utils/shared_files.locked on the table for one read-change-write and replace it through a temporary
file, so concurrent T7 runs on one night lose no rows and a reader sees the old table or the new one.

Nothing here may end a night: astrometry.py turns every failure into a warning and solves the frames instead.
"""
import importlib
import os
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.table import Table

from utils import orchard_version, shared_files

SCHEMA = 1                      # the table layout; a table of another layout is unreadable (and replaced)

KEY = 'raw'
# Everything a row was solved under. A row whose fields differ from this run's is ignored.
VALIDITY = ('raw_size', 'raw_mtime_ns',           # the raw frame (a re-download, rotation or repair changes them)
            'naxis1', 'naxis2', 'trim_x', 'trim_y',  # the calibrated frame's shape and where it sits in the raw one
            'solver_epoch',                       # astrometry.SOLVER_EPOCH: bumped by any change to the solutions
            'gaia_db', 'gaia_db_size', 'gaia_db_mtime_ns',
            'astropy', 'twirl')
PROVENANCE = ('orchver', 'orchrel', 'commit', 'solved_at', 'matches', 'proc')
CARDS = 'cards'
COLUMNS = (KEY,) + VALIDITY + PROVENANCE + (CARDS,)
INT_COLUMNS = {'raw_size', 'raw_mtime_ns', 'naxis1', 'naxis2', 'trim_x', 'trim_y', 'solver_epoch',
               'gaia_db_size', 'gaia_db_mtime_ns', 'matches'}

# Cards pointer_wcs writes besides the WCS keywords it clears first: astropy's wcs.to_header() adds the time
# reference (MJDREF in the astropy versions orchard runs; the others for completeness).
TIME_KEYWORDS = ('MJDREF', 'MJDREFI', 'MJDREFF', 'DATEREF')
# A row without these cannot be a solve.
REQUIRED_CARDS = ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2')
CARD_LEN = 80


class CacheError(Exception):
    """The cache could not be read, or a row cannot be applied."""


def _package_version(name):
    try:
        from importlib import metadata
        return metadata.version(name)
    except Exception:
        pass
    try:
        return str(getattr(importlib.import_module(name), '__version__'))
    except Exception:
        return 'unknown'


def run_fields(db_path, trim_offsets, solver_epoch):
    """The validity fields shared by every frame of one T7 run. Raises OSError if the Gaia database is missing."""
    st = os.stat(db_path)
    return {'trim_x': int(trim_offsets[0]), 'trim_y': int(trim_offsets[1]), 'solver_epoch': int(solver_epoch),
            'gaia_db': os.path.abspath(db_path), 'gaia_db_size': int(st.st_size),
            'gaia_db_mtime_ns': int(st.st_mtime_ns),
            'astropy': _package_version('astropy'), 'twirl': _package_version('twirl')}


def frame_fields(raw_path, proc_header):
    """The validity fields of one frame: its raw file now, and the calibrated frame's shape."""
    st = os.stat(raw_path)
    return {'raw_size': int(st.st_size), 'raw_mtime_ns': int(st.st_mtime_ns),
            'naxis1': int(proc_header['NAXIS1']), 'naxis2': int(proc_header['NAXIS2'])}


def code_identity():
    """{'version', 'commit', 'release'} of the running code."""
    i = orchard_version.info()
    return {'version': i.get('version') or '', 'commit': i.get('commit') or '', 'release': i.get('release') or ''}


def reusable_by(entry, code):
    """True if a row's writer may feed the running code: a code release wrote it, or the very commit that is
    running now (never a checkout with uncommitted changes)."""
    if entry.get('orchrel'):
        return True
    commit = entry.get('commit')
    return bool(commit) and commit == code.get('commit') and not str(entry.get('orchver', '')).endswith('-dirty')


def check(entry, expected, code):
    """None if `entry` may be applied to a frame whose validity fields are `expected`, else why not: the first
    field that differs, 'no entry', 'provenance' or 'no cards'."""
    if entry is None:
        return 'no entry'
    for field in VALIDITY:
        if entry.get(field) != expected.get(field):
            return field
    if not reusable_by(entry, code):
        return 'provenance'
    if not entry.get(CARDS):
        return 'no cards'
    return None


def cards_from_header(header, wcs_keywords):
    """The solve's cards in a solved frame's header, as one string of 80-character card images (trailing blanks
    dropped; the table would keep them, but a reader need not)."""
    keys = set(wcs_keywords) | set(TIME_KEYWORDS)
    return ''.join(card.image for card in header.cards if card.keyword in keys).rstrip()


def cards_header(cards):
    """The Header a row's cards string stands for. Raises CacheError if it is not a complete celestial WCS."""
    try:
        text = str(cards)
        n = -(-len(text) // CARD_LEN)
        header = fits.Header.fromstring(text.ljust(n * CARD_LEN))
    except Exception as err:
        raise CacheError('cards do not parse: {}'.format(err))
    missing = [k for k in REQUIRED_CARDS if k not in header]
    if missing:
        raise CacheError('cards lack {}'.format(', '.join(missing)))
    try:
        from astropy.wcs import WCS
        if not WCS(header).has_celestial:
            raise CacheError('cards are not a celestial WCS')
    except CacheError:
        raise
    except Exception as err:
        raise CacheError('cards are not a WCS: {}'.format(err))
    return header


def apply(proc_path, cards, wcs_keywords):
    """Write a row's cards into a calibrated frame as a fresh solve would (pointer_wcs: clear the WCS keywords,
    header.update() with the solution), and record the edit. Header only: the pixels are never read, so a 16-bit
    frame stays 16-bit. The cards are checked before the file is opened, so a bad row leaves the frame untouched."""
    from astrom.pointer_wcs import clear_wcs_headers
    cached = cards_header(cards)
    with fits.open(proc_path, mode='update', do_not_scale_image_data=True) as hdul:
        header = hdul[0].header
        if clear_wcs_headers(header, wcs_keywords):
            orchard_version.record_edit(header, 'T7', 'inherited WCS removed')
        header.update(cached)
        orchard_version.record_edit(header, 'T7', 'WCS from cache')


def make_entry(raw_name, fields, cards, matches=None, proc=''):
    """A row for a frame this run solved: its validity fields (run_fields + frame_fields), who solved it and when,
    and its cards."""
    code = code_identity()
    entry = {field: fields[field] for field in VALIDITY}
    entry.update({KEY: raw_name, CARDS: cards, 'orchver': code['version'], 'orchrel': code['release'],
                  'commit': code['commit'], 'solved_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                  'matches': -1 if matches is None else int(matches), 'proc': proc})
    return entry


def _plain(value, column):
    """A table cell as a plain Python value: '' (or -1) for an empty cell, int or str otherwise."""
    if np.ma.is_masked(value):
        return -1 if column in INT_COLUMNS else ''
    if column in INT_COLUMNS:
        return int(value)
    return str(value)


def load(path):
    """{raw name: row} from a night's table; {} if there is none yet. Raises CacheError if it cannot be read."""
    if not os.path.exists(path):
        return {}
    try:
        table = Table.read(path, format='ascii.ecsv')
    except Exception as err:
        raise CacheError('cannot read {}: {}'.format(path, err))
    if table.meta.get('schema') != SCHEMA:
        raise CacheError('{} has layout {!r}, not {}'.format(path, table.meta.get('schema'), SCHEMA))
    missing = [c for c in COLUMNS if c not in table.colnames]
    if missing:
        raise CacheError('{} lacks columns {}'.format(path, ', '.join(missing)))
    entries = {}
    for row in table:
        entry = {c: _plain(row[c], c) for c in COLUMNS}
        entries[entry[KEY]] = entry
    return entries


def _table(entries):
    rows = [[entries[k][c] for c in COLUMNS] for k in sorted(entries)]
    dtype = [np.int64 if c in INT_COLUMNS else str for c in COLUMNS]
    table = Table(rows=rows, names=COLUMNS, dtype=dtype)
    table.meta['schema'] = SCHEMA
    table.meta['note'] = 'orchard T7 plate-solve cache; see docs/wcs-cache.md'
    return table


def store(path, new_entries):
    """Add or replace rows in a night's table: one read-change-write under its lock, written through a temporary
    file. An unreadable table is moved aside to <path>.corrupt and replaced. Returns the number of rows written."""
    if not new_entries:
        return 0
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with shared_files.locked(path):
        try:
            entries = load(path)
        except CacheError as err:
            print('WARNING: WCS cache {}; moving it to {}.corrupt and starting a new table'.format(err, path),
                  flush=True)
            os.replace(path, path + '.corrupt')
            entries = {}
        entries.update(new_entries)
        tmp = '{}.tmp{}'.format(path, os.getpid())
        try:
            _table(entries).write(tmp, format='ascii.ecsv', overwrite=True)
            os.replace(tmp, path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
    return len(new_entries)
