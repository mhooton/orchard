#!/usr/bin/env python3
"""
build_target_table.py
=====================

Regenerate the SPECULOOS 40 pc target table, ``ml_40pc_v2.csv``.

The table is a *generated artefact*, never hand-edited: a new 40 pc sample
release is a regeneration, not a merge.  It is a pure function of three
committed inputs, so running this twice on the same inputs gives the same
bytes:

* ``ml_40pc.txt``          the sample as circulated — the source of truth
                           for every astrophysical parameter
* ``gaia_dr3_resolved_ids.csv``  Gaia DR3 identifications for the rows whose
                           ``Gaia_ID`` is zero, with a match strength
* ``gaia_dr2_dr3_map.csv`` DR2 -> DR3 identifier map and DR3 parallaxes for
                           every other row

The two CSVs are rebuilt from the Gaia archive by
``resolve_gaia_dr3_map.py``; this script never touches the network.

What changes, and what does not
-------------------------------

*Format.*  ``ml_40pc.txt`` has a comma-separated header over
whitespace-separated data rows, so a whitespace reader turns the commas
into part of the column names (``Gaia_ID,``, ``T_eff,``) while ``Program``,
with nothing after it, escapes.  The output is comma-separated throughout.

*Identifiers.*  The legacy ``Gaia_ID`` column is a DR2 source ID under an
unqualified name.  It becomes two explicit columns, ``Gaia_DR2_ID`` and
``Gaia_DR3_ID``, either of which may be empty: 4 of the stars we identified
in DR3 have no DR2 entry at all (Wolf 359 among them — DR2's transit
matching could not follow a star moving 4.7"/yr), and a star absent from
Gaia has neither.

*Parameters.*  Every astrophysical column is copied across as the source
spelled it, token for token — no parse, no reformat — so values are
identical to ``ml_40pc.txt`` by construction.

*Distances.*  ``parallax``/``DR3_dist_pc`` and the three flags are added
from Gaia DR3.  ``Dis``/``e_Dis`` are deliberately left alone: ``M``, ``R``
and ``T_eff`` were derived from them and have not been recomputed, so
overwriting ``Dis`` would leave each row internally inconsistent.  Read
``Dis`` as "the distance the derived parameters assume" and
``DR3_dist_pc`` as "the best distance we have".

*Membership.*  No row is ever dropped.  ``dist_gt_40pc`` and ``flagged``
are advisory columns, not a filter: the table's purpose is to find as many
M and L dwarfs as possible, so a star that has drifted beyond the 40 pc
selection boundary stays in with the flag set.

Usage
-----
    python3 -m utils.build_target_table ml_40pc.txt -o ml_40pc_v2.csv
"""

import argparse
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESOLVED = os.path.join(HERE, '..', 'calibration',
                                'gaia_dr3_resolved_ids.csv')
DEFAULT_MAP = os.path.join(HERE, '..', 'calibration', 'gaia_dr2_dr3_map.csv')

# The sample as circulated.  A change here is a new release, not a bug, but
# it should be noticed rather than absorbed silently.
EXPECT_ROWS = 14168

# Columns copied verbatim from the source, in output order.  'Sp_ID' and
# '2MASS_ID' lead; the Gaia ID columns are inserted after them.
_PASSTHROUGH_HEAD = ('Sp_ID', '2MASS_ID')
_PASSTHROUGH_TAIL = (
    'RA', 'DEC', 'G', 'I', 'J', 'H', 'K', 'Dis', 'e_Dis', 'M', 'e_M',
    'R', 'e_R', 'T_eff', 'e_Teff', 'SpT', 'e_Spt', 'SNR_TESS_temp',
    'SNR_Spec_temp', 'SNR_TESS_HZ', 'SNR_Spec_HZ', 'SNR_JWST_HZ_tr',
    'SNR_JWST_HZ_occ', 'SNR_JWST_temp_occ', 'Program',
)
OUTPUT_COLUMNS = (
    _PASSTHROUGH_HEAD + ('Gaia_DR2_ID', 'Gaia_DR3_ID') + _PASSTHROUGH_TAIL +
    ('parallax', 'parallax_error', 'DR3_dist_pc', 'DR3_dist_err_pc',
     'low_conf_dist', 'dist_gt_40pc', 'flagged', 'dr3_match')
)


def read_legacy(path):
    """
    Read ``ml_40pc.txt`` as raw tokens.

    The header is split on commas *and* whitespace so the names come back
    clean; data rows are split on whitespace alone.  Values stay strings so
    that a float never round-trips through this script.
    """
    with open(path) as f:
        raw = f.read().replace('\r\n', '\n').replace('\r', '\n')
    lines = raw.split('\n')
    header_line = next(l for l in lines if l.strip())
    names = [h.strip() for h in re.split(r'[,\s]+', header_line.strip())
             if h.strip()]
    rows = []
    for line in lines[lines.index(header_line) + 1:]:
        if line.strip():
            rows.append(re.split(r'\s+', line.strip()))
    return names, rows


def is_absent(value):
    """True when an ID cell means 'no identifier'."""
    v = (value or '').strip().lower()
    return v in ('', 'nan', 'none', 'null', '--') or set(v) <= set('0')


def _key(sp_id, ra, dec):
    """
    Row key.  Sp_ID alone will not do: 418 names are duplicated, because a
    resolved binary gives both components the same coordinate-derived name.
    Rounding to 5 dp is ~0.04" at the equator, far below the separation of
    any two rows but tolerant of reformatting.
    """
    try:
        return (sp_id.strip(), round(float(ra), 5), round(float(dec), 5))
    except (TypeError, ValueError):
        return (sp_id.strip(), None, None)


def load_resolved(path):
    """Curated DR3 identifications for the zero-ID rows, keyed by position."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            out[_key(r['Sp_ID'], r.get('RA'), r.get('DEC'))] = r
    return out


def load_map(path):
    """DR2 -> DR3 map and DR3 parallaxes, keyed by DR2 source ID."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            out[r['dr2'].strip()] = r
    return out


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v          # NaN -> None


def distance_fields(parallax, parallax_error):
    """
    Distance and flags from a DR3 parallax, using the criteria Ben Rackham
    applied to the version of this table circulated in November 2024, so
    that our flags stay comparable with those:

      low_conf_dist  distance detected at < 3 sigma
      dist_gt_40pc   beyond 40 pc at >= 3 sigma
      flagged        either of the above

    A non-positive parallax yields no distance, and the flags are then left
    empty rather than False: unknown is not the same as not flagged.
    """
    plx, err = _f(parallax), _f(parallax_error)
    blank = ('', '', '', '', '')
    if plx is None or err is None or plx <= 0:
        return blank
    dist = 1000.0 / plx
    dist_err = dist * (err / plx)
    if dist_err <= 0:
        return blank
    snr = dist / dist_err
    low = snr < 3.0
    gt40 = (dist > 40.0) and (snr >= 3.0)
    return ('%.10g' % dist, '%.10g' % dist_err,
            str(bool(low)), str(bool(gt40)), str(bool(low or gt40)))


def build(source, resolved, id_map, allow_weak=False, expect_rows=EXPECT_ROWS,
          warn=None):
    """
    Return (header, rows) for the new table.

    `warn` is called with a message for every condition worth a human's
    attention; the caller decides whether that is a log line or stderr.
    """
    warn = warn or (lambda msg: None)
    names, src_rows = read_legacy(source)

    if len(src_rows) != expect_rows:
        warn("source row count changed: %d rows, expected %d. The 40 pc "
             "sample may have been re-released; check that "
             "gaia_dr3_resolved_ids.csv and gaia_dr2_dr3_map.csv still "
             "cover it before trusting the output."
             % (len(src_rows), expect_rows))

    idx = {n.upper(): i for i, n in enumerate(names)}
    missing = [c for c in _PASSTHROUGH_HEAD + _PASSTHROUGH_TAIL
               if c.upper() not in idx]
    if missing:
        raise KeyError("source is missing column(s): %s" % ', '.join(missing))
    legacy_id = idx['GAIA_ID']

    def cell(row, name):
        i = idx[name.upper()]
        return row[i] if i < len(row) else ''

    out_rows = []
    stats = {}
    unmapped, weak_skipped = [], []

    for row in src_rows:
        sp = cell(row, 'Sp_ID')
        ra, dec = cell(row, 'RA'), cell(row, 'DEC')
        legacy = row[legacy_id] if legacy_id < len(row) else ''

        dr2 = dr3 = ''
        plx = plx_err = ''
        how = ''

        if is_absent(legacy):
            # A zero-ID row: identified, if at all, by the curated file.
            rec = resolved.get(_key(sp, ra, dec))
            if rec is None:
                how = 'unresolved'
            elif (rec.get('match_strength') or '').strip() == 'weak' \
                    and not allow_weak:
                how = 'weak_excluded'
                weak_skipped.append(sp)
            else:
                dr3 = (rec.get('gaia_dr3_id') or '').strip()
                dr2 = (rec.get('gaia_dr2_id') or '').strip()
                plx = (rec.get('parallax') or '').strip()
                plx_err = (rec.get('parallax_error') or '').strip()
                how = (rec.get('dr3_match') or '').strip() or 'resolved'
        else:
            dr2 = legacy.strip()
            rec = id_map.get(dr2)
            if rec is None:
                how = 'unmapped'
                unmapped.append(dr2)
            else:
                how = (rec.get('kind') or '').strip()
                if how in ('identical', 'differs'):
                    dr3 = (rec.get('dr3') or '').strip()
                    plx = (rec.get('parallax') or '').strip()
                    plx_err = (rec.get('parallax_error') or '').strip()
                # 'ambiguous' and 'no_dr3' leave Gaia_DR3_ID empty on
                # purpose: an identifier we cannot pin down is worse than
                # none, because the pipeline would silently adopt it.

        stats[how] = stats.get(how, 0) + 1
        dist = distance_fields(plx, plx_err)

        out = [sp, cell(row, '2MASS_ID'), dr2, dr3]
        out += [cell(row, c) for c in _PASSTHROUGH_TAIL]
        out += [plx, plx_err] + list(dist) + [how]
        out_rows.append(out)

    if unmapped:
        warn("%d row(s) carry a DR2 ID absent from gaia_dr2_dr3_map.csv "
             "(e.g. %s); they keep their DR2 ID but get no DR3 ID or "
             "distance. Re-run resolve_gaia_dr3_map.py."
             % (len(unmapped), ', '.join(unmapped[:3])))
    if weak_skipped:
        warn("%d weak match(es) left unidentified (%s). These have neither "
             "a proper motion nor a 2MASS counterpart, so a neighbour is as "
             "likely as the target; pass --allow-weak to write them."
             % (len(weak_skipped), ', '.join(sorted(weak_skipped)[:5])))

    return list(OUTPUT_COLUMNS), out_rows, stats


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('source', help='path to ml_40pc.txt')
    p.add_argument('-o', '--out', required=True, help='output CSV')
    p.add_argument('--resolved', default=DEFAULT_RESOLVED)
    p.add_argument('--map', dest='id_map', default=DEFAULT_MAP)
    p.add_argument('--allow-weak', action='store_true',
                   help='also write the weak matches (see the report)')
    p.add_argument('--expect-rows', type=int, default=EXPECT_ROWS)
    a = p.parse_args(argv)

    warnings = []

    def warn(msg):
        warnings.append(msg)
        sys.stderr.write('WARNING: %s\n' % msg)

    header, rows, stats = build(
        a.source, load_resolved(a.resolved), load_map(a.id_map),
        allow_weak=a.allow_weak, expect_rows=a.expect_rows, warn=warn)

    with open(a.out, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)

    n_dr2 = sum(1 for r in rows if r[2])
    n_dr3 = sum(1 for r in rows if r[3])
    n_plx = sum(1 for r in rows if r[header.index('parallax')])
    n_flag = sum(1 for r in rows if r[header.index('flagged')] == 'True')
    print("wrote %s: %d rows, %d columns" % (a.out, len(rows), len(header)))
    print("  Gaia_DR2_ID present : %d" % n_dr2)
    print("  Gaia_DR3_ID present : %d" % n_dr3)
    print("  DR3 parallax        : %d" % n_plx)
    print("  flagged (advisory)  : %d  (no row is dropped)" % n_flag)
    print("  dr3_match: %s" % ', '.join(
        '%s=%d' % kv for kv in sorted(stats.items(), key=lambda kv: -kv[1])))
    return 0


if __name__ == '__main__':
    sys.exit(main())
