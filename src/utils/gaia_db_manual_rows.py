#!/usr/bin/env python3
"""
gaia_db_manual_rows.py
======================

Rows added by hand to the band tables of the local Gaia database, for real Gaia stars that the
build (gaia-tmass-sqlite db_maker.py) left out. The build keeps a DR3 source only if its DR2
counterpart has a 2MASS J < 16 through Gaia DR2's 2MASS cross-match, looked up in the same 0.01 deg
declination strip, so it drops stars DR2 never matched to 2MASS (most of them with poor astrometry:
binaries) and some fast movers whose DR2 and DR3 positions fall in different strips. Without a row
a target can never be identified, as primary or as secondary, on any night.

The source of truth is calibration/gaia_db_manual_rows.csv in this repo. Every row there says how
each value was chosen, as codes defined in that file's header (DR3 verbatim, proper motion from DR2,
proper motion from 2MASS positions, ...), with where each value came from.

`apply` inserts the rows a database lacks and records each one in the database's `manual_rows`
table, so the database itself says which of its rows were made by hand and how. It is idempotent:
a row already present with the same values is only recorded (once); a row whose ID is present with
other values, or a row with another star within 0.5", is reported and left alone. Re-run it after
any rebuild of the database.

Usage
-----
    python -m utils.gaia_db_manual_rows check [--db PATH ...] [--csv PATH]
    python -m utils.gaia_db_manual_rows apply [--db PATH ...] [--csv PATH]
    python -m utils.gaia_db_manual_rows list  [--db PATH ...]

`check` only reads. `apply` writes in place, one short transaction per band (0.5 s lock timeout,
retried), so it can run while the pipeline reads the database: a reader waits at most about 0.5 s.
Run it on appct inside orchard-server (local disk), never over NFS. In place matters:
goibniu-server-autosave bind-mounts the database file itself, so a copy and rename would leave it on
the old file. The database defaults to $GAIADATABASEPATH, then the pipeline's standard location.

Undo a row: DELETE FROM '<band>' WHERE rowid = <row_id>, and its entry from manual_rows, using
`list` for the row_id.
"""

import argparse
import csv
import math
import os
import sqlite3
import sys
import time

COLUMNS = ["ra", "dec", "pmra", "pmdec", "phot_g_mean_mag", "g_rp", "bp_rp", "parallax",
           "teff_gspphot", "source_id", "dr2_source_id", "j_m"]
CSV_FIELDS = ["sp_id", "band"] + COLUMNS + ["methods", "sources", "notes", "added"]
METHODS = {
    "DR2_TMASS_MISS", "STRIP_BOUNDARY", "SIX_PARAM", "DR2_ID_FROM_PLANS", "J_FROM_ML40_2MASS",
    "NEAREST_OF_SPLIT", "PM_FROM_DR2", "PM_FROM_2MASS", "COMPANION_IN_DB", "UNRESOLVED_PAIR",
    "NOT_UCD", "NO_J", "ML40_ZERO_ID", "SP_ID_SIBLING", "J_FROM_2MASS_POSITION", "DR2_ID_FROM_ML40",
}
SYNTHETIC_ID_MIN = 9_000_000_000_000_000_000   # supplementary_sources' range; never a real Gaia ID
SKIP_ARCSEC = 0.5      # another row this close is very likely the same star under another ID
REPORT_ARCSEC = 3.0
LOCK_TIMEOUT = 0.5
RETRIES = 40

PROVENANCE_SQL = (
    "CREATE TABLE IF NOT EXISTS manual_rows ("
    "band TEXT, row_id INTEGER, source_id TEXT, dr2_source_id TEXT, sp_id TEXT, "
    "methods TEXT, sources TEXT, notes TEXT, added TEXT, recorded TEXT)"
)

DEFAULT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "calibration", "gaia_db_manual_rows.csv")


def default_db_path():
    return os.environ.get("GAIADATABASEPATH", "/gaia_database/gaia_dr3_unified_16jcut.db")


def expected_band(dec):
    lo = max(-90, min(89, int(math.floor(dec))))
    return "%d_%d" % (lo, lo + 1)


def _num(x):
    x = (x or "").strip()
    return float(x) if x else None


def read_csv(path):
    """Read the curated CSV, skipping comment lines; validate every row. Returns dicts with a
    `values` tuple in COLUMNS order (floats or None, IDs as str)."""
    with open(path, newline="") as f:
        lines = [l for l in f if not l.lstrip().startswith("#")]
    rows = []
    for r in csv.DictReader(lines):
        if not (r.get("source_id") or "").strip():
            continue
        name = (r.get("sp_id") or "").strip()
        sid = r["source_id"].strip()
        dr2 = (r.get("dr2_source_id") or "").strip() or None
        for v in (sid, dr2):
            if v is not None and (not v.isdigit() or int(v) >= SYNTHETIC_ID_MIN):
                raise ValueError("%s: %r is not a real Gaia source_id (synthetic objects go in "
                                 "supplementary_sources.csv)" % (name, v))
        if _num(r.get("ra")) is None or _num(r.get("dec")) is None:
            raise ValueError("%s: ra and dec are required" % name)
        band = (r.get("band") or "").strip()
        if band != expected_band(_num(r["dec"])):
            raise ValueError("%s: band %r, but dec %s belongs in %s" % (name, band, r["dec"],
                                                                        expected_band(_num(r["dec"]))))
        methods = [m.strip() for m in (r.get("methods") or "").split(";") if m.strip()]
        unknown = set(methods) - METHODS
        if not methods or unknown:
            raise ValueError("%s: methods %r (unknown: %s)" % (name, r.get("methods"), sorted(unknown)))
        values = tuple(sid if c == "source_id" else dr2 if c == "dr2_source_id" else _num(r.get(c))
                       for c in COLUMNS)
        rows.append(dict(sp_id=name, band=band, values=values, methods=";".join(methods),
                         sources=(r.get("sources") or "").strip(), notes=(r.get("notes") or "").strip(),
                         added=(r.get("added") or "").strip()))
    ids = [r["values"][9] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate source_id in %s" % path)
    return rows


def sep_arcsec(ra1, de1, ra2, de2):
    d = math.radians
    s = (math.sin(d(de2 - de1) / 2) ** 2
         + math.cos(d(de1)) * math.cos(d(de2)) * math.sin(d(ra2 - ra1) / 2) ** 2)
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(s)))) * 3600


def _survey(conn, band, rows):
    """One pass over a band: rows carrying any candidate's ID or lying within REPORT_ARCSEC of one."""
    ids = sorted({x for r in rows for x in r["values"][9:11] if x})
    where, args = [], []
    if ids:
        marks = ",".join("?" * len(ids))
        where.append("source_id IN (%s) OR dr2_source_id IN (%s)" % (marks, marks))
        args += ids + ids
    for r in rows:
        ra, dec = r["values"][0], r["values"][1]
        h = REPORT_ARCSEC / 3600.0
        hr = h / max(math.cos(math.radians(dec)), 1e-6)
        where.append("(dec BETWEEN ? AND ? AND ra BETWEEN ? AND ?)")
        args += [dec - h, dec + h, ra - hr, ra + hr]
    q = "SELECT rowid, %s FROM '%s' WHERE %s" % (", ".join(COLUMNS), band, " OR ".join(where))
    return [(h[0], tuple(h[1:])) for h in conn.execute(q, args).fetchall()]


def _classify(row, hits):
    """('present', rowid) | ('conflict', rowid) | ('near', rowid) | ('missing', None), plus notes."""
    v = row["values"]
    ids = {x for x in v[9:11] if x}
    notes = []
    same = [rid for rid, h in hits if {x for x in h[9:11] if x} & ids]
    exact = [rid for rid, h in hits if h == v]
    for rid, h in hits:
        if rid not in same:
            s = sep_arcsec(v[0], v[1], h[0], h[1])
            if s <= REPORT_ARCSEC:
                notes.append('other row %.2f" away: rowid %s, source_id %s, G %s' % (s, rid, h[9], h[4]))
    if exact:
        return "present", exact[0], notes
    if same:
        return "conflict", same[0], notes
    near = sorted((sep_arcsec(v[0], v[1], h[0], h[1]), rid) for rid, h in hits if rid not in same)
    if near and near[0][0] < SKIP_ARCSEC:
        return "near", near[0][1], notes
    return "missing", None, notes


def _provenance(conn):
    try:
        return {(b, r) for b, r in conn.execute("SELECT band, row_id FROM manual_rows")}
    except sqlite3.OperationalError:
        return set()


def _retry(conn, sql, what):
    for attempt in range(RETRIES):
        try:
            conn.execute(sql)
            return
        except sqlite3.OperationalError as e:
            if "locked" not in str(e) and "busy" not in str(e):
                raise
            print("   %s busy (%d): %s" % (what, attempt, e), flush=True)
            time.sleep(2)
    raise RuntimeError("%s kept failing: database locked" % what)


def run(db_path, rows, write=False, log=print):
    """Check (write=False) or apply (write=True) the rows on one database. Returns a list of
    (sp_id, band, status, row_id) with status one of present, recorded, inserted, conflict, near,
    missing."""
    if write:
        conn = sqlite3.connect(db_path, timeout=LOCK_TIMEOUT, isolation_level=None)
    else:
        conn = sqlite3.connect("file:%s?mode=ro" % db_path, uri=True)
    out = []
    try:
        recorded = _provenance(conn)
        for band in sorted({r["band"] for r in rows}, key=lambda b: int(b.split("_")[0])):
            cands = [r for r in rows if r["band"] == band]
            if write:
                _retry(conn, "BEGIN IMMEDIATE", "BEGIN")
            try:
                t0 = time.time()
                hits = _survey(conn, band, cands)
                todo_insert, todo_record = [], []
                for r in cands:
                    status, rid, notes = _classify(r, hits)
                    for n in notes:
                        log("   %s: %s" % (r["sp_id"], n))
                    if status == "present" and (band, rid) not in recorded:
                        todo_record.append((r, rid))
                    elif status == "missing":
                        todo_insert.append(r)
                    else:
                        out.append((r["sp_id"], band, status, rid))
                if write and (todo_insert or todo_record):
                    conn.execute(PROVENANCE_SQL)
                    now = time.strftime("%Y-%m-%dT%H:%M:%S")
                    for r in todo_insert:
                        rid = conn.execute("INSERT INTO '%s' (%s) VALUES (%s)" % (
                            band, ", ".join(COLUMNS), ",".join("?" * len(COLUMNS))), r["values"]).lastrowid
                        todo_record.append((r, rid))
                    for r, rid in todo_record:
                        conn.execute("INSERT INTO manual_rows VALUES (?,?,?,?,?,?,?,?,?,?)",
                                     (band, rid, r["values"][9], r["values"][10], r["sp_id"], r["methods"],
                                      r["sources"], r["notes"], r["added"], now))
                    _retry(conn, "COMMIT", "COMMIT")
                    inserted = {id(r) for r in todo_insert}
                    for r, rid in todo_record:
                        got = conn.execute("SELECT %s FROM '%s' WHERE rowid = ?" % (", ".join(COLUMNS), band),
                                           (rid,)).fetchone()
                        if tuple(got) != r["values"]:
                            raise RuntimeError("%s: rowid %s in %s does not read back as written" % (
                                r["sp_id"], rid, band))
                        out.append((r["sp_id"], band, "inserted" if id(r) in inserted else "recorded", rid))
                        recorded.add((band, rid))
                else:
                    if write:
                        conn.execute("ROLLBACK")
                    out += [(r["sp_id"], band, "missing", None) for r in todo_insert]
                    out += [(r["sp_id"], band, "present, not yet recorded", rid) for r, rid in todo_record]
                log(" %s: %d candidate(s), %.1f s" % (band, len(cands), time.time() - t0))
            except BaseException:
                if write and conn.in_transaction:
                    conn.execute("ROLLBACK")
                raise
        already = {(r["sp_id"], r["band"]) for r in rows} - {(o[0], o[1]) for o in out}
        out += [(sp, band, "present", None) for sp, band in sorted(already)]
    finally:
        conn.close()
    return out


def list_rows(db_path):
    conn = sqlite3.connect("file:%s?mode=ro" % db_path, uri=True)
    try:
        return conn.execute("SELECT band, row_id, sp_id, source_id, methods, added, recorded "
                            "FROM manual_rows ORDER BY sp_id").fetchall()
    except sqlite3.OperationalError:
        return []
    finally:
        conn.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("check", "apply", "list"):
        s = sub.add_parser(name)
        s.add_argument("--db", action="append", help="database (repeatable; default $GAIADATABASEPATH)")
        if name != "list":
            s.add_argument("--csv", default=DEFAULT_CSV)
    a = p.parse_args(argv)
    dbs = a.db or [default_db_path()]
    if a.cmd == "list":
        for db in dbs:
            print("==", db)
            for r in list_rows(db):
                print("  ", r)
        return 0
    rows = read_csv(a.csv)
    print("%d row(s) in %s" % (len(rows), a.csv))
    bad = 0
    for db in dbs:
        print("==", db, "(apply)" if a.cmd == "apply" else "(check, read-only)", flush=True)
        for sp, band, status, rid in run(db, rows, write=a.cmd == "apply"):
            if status != "present":
                print("   %-14s %-8s %s%s" % (sp, band, status, "" if rid is None else " (rowid %s)" % rid))
            bad += status in ("conflict", "near")
    if bad:
        print("%d row(s) left alone: an ID present with other values, or another star within %.1f\""
              % (bad, SKIP_ARCSEC))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
