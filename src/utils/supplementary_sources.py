#!/usr/bin/env python3
"""
supplementary_sources.py
========================

Maintain the `custom_sources` table of the local Gaia database: hand-added
rows for objects that exist in neither Gaia DR2 nor DR3 (very faint T/Y
dwarfs, mostly).  Once loaded, the crossmatch in photometry.gaia_dr2_test
returns them alongside real Gaia rows and every later stage treats them
identically.

The source of truth is calibration/supplementary_sources.csv in this repo.
Loading is idempotent: the table is emptied and refilled from the CSV.

Usage
-----
    python -m utils.supplementary_sources load  [--db PATH] [--csv PATH]
    python -m utils.supplementary_sources list  [--db PATH]
    python -m utils.supplementary_sources lookup SOURCE_ID [--db PATH]

The database path defaults to $GAIADATABASEPATH, then the pipeline's
standard location.
"""

import argparse
import csv
import os
import sqlite3
import sys

SYNTHETIC_ID_MIN = 9_000_000_000_000_000_000  # above any real Gaia source_id
INT64_MAX = 9_223_372_036_854_775_807

COLUMNS = ["ra", "dec", "pmra", "pmdec", "phot_g_mean_mag", "g_rp", "bp_rp",
           "parallax", "teff_gspphot", "source_id", "dr2_source_id", "j_m",
           "name", "provenance"]

CREATE_SQL = (
    "CREATE TABLE IF NOT EXISTS custom_sources ("
    "ra REAL, dec REAL, pmra REAL, pmdec REAL, phot_g_mean_mag REAL, "
    "g_rp REAL, bp_rp REAL, parallax REAL, teff_gspphot REAL, "
    "source_id TEXT, dr2_source_id TEXT, j_m REAL, "
    "name TEXT, provenance TEXT)"
)

DEFAULT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "calibration", "supplementary_sources.csv")


def default_db_path():
    return os.environ.get("GAIADATABASEPATH", "/gaia_database/gaia_dr3_unified_16jcut.db")


def is_synthetic_id(value):
    """True if value is a synthetic ID from the reserved range."""
    try:
        v = int(str(value).strip())
    except (TypeError, ValueError):
        return False
    return SYNTHETIC_ID_MIN <= v <= INT64_MAX


def read_csv(path):
    """Read the curated CSV, skipping comment lines; validate every row."""
    rows = []
    with open(path, newline="") as f:
        lines = [l for l in f if not l.lstrip().startswith("#")]
    for r in csv.DictReader(lines):
        if not r.get("source_id", "").strip():
            continue
        sid = r["source_id"].strip()
        if not is_synthetic_id(sid):
            raise ValueError(
                "%s: source_id %s is outside the synthetic range [%d, %d]"
                % (r.get("name"), sid, SYNTHETIC_ID_MIN, INT64_MAX))
        for k in ("ra", "dec"):
            if not r.get(k, "").strip():
                raise ValueError("%s: %s is required" % (r.get("name"), k))
        rows.append(r)
    ids = [r["source_id"].strip() for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate synthetic source_id in %s" % path)
    return rows


def _num(x):
    x = (x or "").strip()
    return float(x) if x else None


def load(db_path, csv_path):
    rows = read_csv(csv_path)
    conn = sqlite3.connect(db_path)
    conn.execute(CREATE_SQL)
    conn.execute("DELETE FROM custom_sources")
    conn.executemany(
        "INSERT INTO custom_sources (%s) VALUES (%s)" % (
            ", ".join(COLUMNS), ", ".join("?" * len(COLUMNS))),
        [(
            _num(r["ra"]), _num(r["dec"]), _num(r.get("pmra")), _num(r.get("pmdec")),
            _num(r.get("phot_g_mean_mag")), _num(r.get("g_rp")), _num(r.get("bp_rp")),
            _num(r.get("parallax")), _num(r.get("teff_gspphot")),
            r["source_id"].strip(), None, _num(r.get("j_m")),
            r.get("name", "").strip(), r.get("provenance", "").strip(),
        ) for r in rows])
    conn.commit()
    conn.close()
    return len(rows)


def list_rows(db_path):
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.execute("SELECT name, source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, "
                           "teff_gspphot, j_m, provenance FROM custom_sources ORDER BY name")
        return cur.fetchall()
    except sqlite3.OperationalError:
        return []
    finally:
        conn.close()


def lookup(db_path, source_id):
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.execute("SELECT * FROM custom_sources WHERE source_id = ?", (str(source_id),))
        return cur.fetchone()
    except sqlite3.OperationalError:
        return None
    finally:
        conn.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("load"); s.add_argument("--db", default=None); s.add_argument("--csv", default=DEFAULT_CSV)
    s = sub.add_parser("list"); s.add_argument("--db", default=None)
    s = sub.add_parser("lookup"); s.add_argument("source_id"); s.add_argument("--db", default=None)
    a = p.parse_args(argv)
    db = a.db or default_db_path()
    if a.cmd == "load":
        n = load(db, a.csv)
        print("loaded %d supplementary source(s) into %s" % (n, db))
    elif a.cmd == "list":
        for r in list_rows(db):
            print(r)
    elif a.cmd == "lookup":
        print(lookup(db, a.source_id))
    return 0


if __name__ == "__main__":
    sys.exit(main())
