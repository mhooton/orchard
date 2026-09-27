"""
target_list.py
==============

Reading the SPECULOOS 40 pc target table.

Two formats are in circulation and everything here reads either, so the
pipeline can be pointed at one or the other with TARGET_LIST and behave
identically:

* the legacy ``ml_40pc.txt``.  Its header is comma-separated but its data
  rows are whitespace-separated, so the whitespace reader every caller uses
  turns the commas into part of the column names: ``Gaia_ID,``, ``T_eff,``
  and so on.  Only the last column, ``Program``, escapes, because nothing
  follows it.  A single ``Gaia_ID`` column holds a Gaia DR2 source ID.
* ``ml_40pc_v2.csv``.  Comma-separated throughout, so the names come back
  clean, and the one ambiguous ``Gaia_ID`` column is split into
  ``Gaia_DR2_ID`` and ``Gaia_DR3_ID``.

Column names are matched case-insensitively with any trailing comma
stripped, which is what makes both formats work through one code path.

A target is identified by Gaia ID and by nothing else, so a row
contributes every ID it carries: `id_rows()` yields one entry per
(row, ID) pair, letting a caller match a catalogue against DR2 and DR3
IDs at once while still recovering the row a match came from.
"""

import logging

from astropy.io import ascii

logger = logging.getLogger(__name__)

# The spellings an absent ID arrives in.  '--' is what astropy prints for a
# masked cell, which is what an empty field in a CSV integer column becomes;
# the legacy list instead writes a row of zeros.
_ABSENT = ('', 'nan', 'none', 'null', 'n/a', '--')

# Candidate ID columns, most specific first.  The legacy 'Gaia_ID' is a DR2
# ID despite the unqualified name.
_ID_COLUMNS = ('GAIA_DR2_ID', 'GAIA_DR3_ID', 'GAIA_ID')


def clean_id(value):
    """
    Normalise one Gaia ID cell to a string, or '' when the ID is absent.

    Absent has several spellings: an empty CSV field (astropy yields a
    masked cell, whose str() is '--'), the all-zeros sentinel the legacy
    list uses, and the usual nan/none.  Returning '' for all of them means
    callers need only one falsiness test.
    """
    v = '' if value is None else str(value).strip()
    if v.lower() in _ABSENT:
        return ''
    if set(v) <= set('0'):        # '0', '000...0' — the legacy sentinel
        return ''
    return v


def _norm(name):
    """Column name without the trailing comma the legacy header leaves."""
    return str(name).strip().rstrip(',').strip().upper()


class TargetList(object):
    """
    A parsed 40 pc target table.

    Row order is the file's own, so a row index is stable across formats
    provided the generator preserves source order (it does).
    """

    def __init__(self, table, path=None):
        self.table = table
        self.path = path
        self.columns = {_norm(c): c for c in table.colnames}
        self._id_cols = [self.columns[n] for n in _ID_COLUMNS
                         if n in self.columns]
        if not self._id_cols:
            raise KeyError("no Gaia ID column in target list %r "
                           "(looked for %s among %s)"
                           % (path, ', '.join(_ID_COLUMNS), table.colnames))
        self._by_id = None

    def __len__(self):
        return len(self.table)

    def has(self, name):
        return _norm(name) in self.columns

    def column(self, name):
        """The column named `name`, matched ignoring case and trailing comma."""
        return self.table[self.columns[_norm(name)]]

    def value(self, name, row):
        try:
            return self.column(name)[row]
        except (KeyError, IndexError):
            return None

    # -- per-row accessors -------------------------------------------------

    def sp_id(self, row):
        v = self.value('Sp_ID', row)
        return '' if v is None else str(v).strip()

    def teff(self, row):
        """Effective temperature as an int, or None when unusable."""
        v = self.value('T_eff', row)
        try:
            t = int(float(v))
        except (TypeError, ValueError):
            return None
        return t

    def coords(self, row):
        """(ra, dec) in degrees, or (None, None)."""
        try:
            return float(self.value('RA', row)), float(self.value('DEC', row))
        except (TypeError, ValueError):
            return None, None

    def ids(self, row):
        """
        Every distinct Gaia ID on this row, DR2 first, absent ones dropped.
        """
        out = []
        for col in self._id_cols:
            gid = clean_id(self.table[col][row])
            if gid and gid not in out:
                out.append(gid)
        return out

    # -- whole-table views -------------------------------------------------

    def id_rows(self, placeholder=None):
        """
        [(gaia_id, row_index), ...] — one entry per (row, ID) pair, so a
        row carrying both a DR2 and a DR3 ID appears twice and matches
        either.

        Rows with no ID at all are represented by `placeholder` when one is
        given (callers that walk the list positionally need every row to be
        present); when it is None they are omitted.
        """
        out = []
        for row in range(len(self.table)):
            gids = self.ids(row)
            if gids:
                out.extend((g, row) for g in gids)
            elif placeholder is not None:
                out.append((placeholder, row))
        return out

    def row_for_id(self, gaia_id):
        """The first row carrying `gaia_id`, or None."""
        if self._by_id is None:
            self._by_id = {}
            for gid, row in self.id_rows():
                self._by_id.setdefault(gid, row)
        return self._by_id.get(clean_id(gaia_id))

    def teff_for_id(self, gaia_id):
        row = self.row_for_id(gaia_id)
        return None if row is None else self.teff(row)


def _looks_like_csv(path):
    """
    Decide the format from the content rather than the file name, since
    TARGET_LIST may point at either under any name.

    The legacy list is the awkward case precisely because its *header* is
    comma-separated while its data rows are not, so the first data line is
    what separates the two formats.
    """
    lines = []
    with open(path) as f:
        for line in f:
            if line.strip():
                lines.append(line)
                if len(lines) == 2:
                    break
    if len(lines) < 2:                     # header only; fall back to the name
        return str(path).lower().endswith('.csv')
    return ',' in lines[1]


def read_target_list(path):
    """
    Read a 40 pc target table in either format.

    Raises whatever astropy raises; callers that must not fail should use
    `load_target_list`.
    """
    if _looks_like_csv(path):
        table = ascii.read(path, format='csv')
    else:
        table = ascii.read(path, delimiter=' ', header_start=0, data_start=1)
    return TargetList(table, path)


def load_target_list(path, logger_=None):
    """
    `read_target_list` that logs and returns None instead of raising.
    """
    log = logger_ or logger
    try:
        return read_target_list(path)
    except Exception as e:                                   # noqa: BLE001
        log.warning("Could not read target list '%s': %s", path, e)
        return None
