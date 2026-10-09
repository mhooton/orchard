#!/usr/bin/env python3
"""
resolve_gaia_dr3_map.py
=======================

Rebuild the two committed inputs that ``build_target_table.py`` consumes:

* ``calibration/gaia_dr2_dr3_map.csv``     DR2 -> DR3 identifier map, with
                                           the DR3 parallax, for every row
                                           of ``ml_40pc.txt`` that carries a
                                           DR2 identifier
* ``calibration/gaia_dr3_resolved_ids.csv`` DR3 identifications for the rows
                                           whose ``Gaia_ID`` is zero, built
                                           from a report produced by
                                           ``resolve_master_list_gaia_ids.py``

This is the only part of the target-table build that touches the network.
It is read-only with respect to every file it did not create.

Why the map exists
------------------

A Gaia DR2 ``source_id`` is not a DR3 ``source_id``.  Most are numerically
unchanged, but a meaningful minority are not, and a few stars have no DR2
entry at all.  Looking a DR2 identifier up directly in
``gaiadr3.gaia_source`` therefore fails in two ways: it silently returns
nothing for a changed identifier, and — worse — where DR2 reused that
number for a different source it returns *the wrong star*, with nothing to
indicate it.  ``gaiadr3.dr2_neighbourhood`` is the table that relates the
two releases properly, and is what this script uses.

Match kinds recorded in the map
-------------------------------

``identical``  the DR3 identifier equals the DR2 one
``differs``    exactly one DR3 counterpart, with a different identifier
``ambiguous``  several DR3 counterparts — the DR2 source was split.  The
               best candidate is recorded for inspection but
               ``build_target_table.py`` does not adopt it: an identifier
               we cannot pin down is worse than none, because the pipeline
               would take it at face value.
``no_dr3``     no counterpart in DR3 at all

Usage
-----
    python3 -m utils.resolve_gaia_dr3_map ml_40pc.txt --outdir calibration \\
        --report resolved.csv
"""

import argparse
import collections
import csv
import io
import os
import re
import sys
import time
import urllib.parse
import urllib.request

TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
BATCH = 500


def tap_query(adql, timeout=300, retries=4):
    """POST an ADQL query; GET cannot carry a batch of 500 identifiers."""
    body = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "QUERY": adql}).encode()
    last = None
    for attempt in range(retries):
        try:
            resp = urllib.request.urlopen(
                urllib.request.Request(TAP, data=body), timeout=timeout)
            return list(csv.DictReader(io.StringIO(resp.read().decode())))
        except Exception as e:                               # noqa: BLE001
            last = e
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("Gaia TAP query failed after %d attempts: %s"
                       % (retries, last))


def batched(ids, size=BATCH):
    for i in range(0, len(ids), size):
        yield i // size + 1, (len(ids) + size - 1) // size, ids[i:i + size]


def read_legacy(path):
    with open(path) as f:
        raw = f.read().replace("\r\n", "\n").replace("\r", "\n")
    lines = raw.split("\n")
    header = next(l for l in lines if l.strip())
    names = [h.strip() for h in re.split(r"[,\s]+", header.strip()) if h.strip()]
    rows = [re.split(r"\s+", l.strip())
            for l in lines[lines.index(header) + 1:] if l.strip()]
    return names, rows


def is_absent(value):
    v = (value or "").strip().lower()
    return v in ("", "nan", "none", "null", "--") or set(v) <= set("0")


def _f(x):
    try:
        return abs(float(x))
    except (TypeError, ValueError):
        return 9e9


def neighbourhood(dr2_ids, log=print):
    """DR2 -> [candidate DR3 rows] from gaiadr3.dr2_neighbourhood."""
    found = collections.defaultdict(list)
    for n, total, batch in batched(dr2_ids):
        rows = tap_query(
            "SELECT dr2_source_id, dr3_source_id, angular_distance, "
            "magnitude_difference FROM gaiadr3.dr2_neighbourhood "
            "WHERE dr2_source_id IN (%s)" % ",".join(batch))
        for r in rows:
            found[r["dr2_source_id"]].append(r)
        log("  neighbourhood batch %d/%d: %d ids -> %d rows (%d mapped)"
            % (n, total, len(batch), len(rows), len(found)))
        time.sleep(1)
    return found


def classify(dr2, candidates):
    """(dr3_id, kind, n_candidates, angular_distance, magnitude_difference)."""
    exact = [c for c in candidates if c["dr3_source_id"] == dr2]
    if exact:
        c = exact[0]
        return c["dr3_source_id"], "identical", len(candidates), \
            c["angular_distance"], c["magnitude_difference"]
    if not candidates:
        return "", "no_dr3", 0, "", ""
    best = sorted(candidates,
                  key=lambda c: (_f(c["angular_distance"]),
                                 _f(c["magnitude_difference"])))[0]
    kind = "differs" if len(candidates) == 1 else "ambiguous"
    return best["dr3_source_id"], kind, len(candidates), \
        best["angular_distance"], best["magnitude_difference"]


def parallaxes(dr3_ids, log=print):
    """DR3 source_id -> (parallax, parallax_error) from gaiadr3.gaia_source."""
    out = {}
    for n, total, batch in batched(dr3_ids):
        rows = tap_query(
            "SELECT source_id, parallax, parallax_error "
            "FROM gaiadr3.gaia_source WHERE source_id IN (%s)"
            % ",".join(batch))
        for r in rows:
            out[str(r["source_id"])] = (r["parallax"], r["parallax_error"])
        log("  parallax batch %d/%d: %d ids -> %d rows (%d total)"
            % (n, total, len(batch), len(rows), len(out)))
        time.sleep(1)
    return out


def build_map(source, log=print):
    names, rows = read_legacy(source)
    gi = names.index("Gaia_ID")
    dr2_ids = sorted({r[gi] for r in rows
                      if gi < len(r) and not is_absent(r[gi])})
    log("%d distinct DR2 identifiers to resolve" % len(dr2_ids))

    found = neighbourhood(dr2_ids, log=log)
    entries = []
    for dr2 in dr2_ids:
        dr3, kind, n, ang, mag = classify(dr2, found.get(dr2, []))
        entries.append({"dr2": dr2, "dr3": dr3, "kind": kind, "n_cands": n,
                        "ang_dist_arcsec": ang, "mag_diff": mag,
                        "parallax": "", "parallax_error": ""})

    # Parallaxes only for identifiers we are confident in; an ambiguous
    # match must not contribute a distance either.
    want = sorted({e["dr3"] for e in entries
                   if e["dr3"] and e["kind"] in ("identical", "differs")})
    log("%d DR3 identifiers to query for parallax" % len(want))
    plx = parallaxes(want, log=log)
    for e in entries:
        if e["dr3"] in plx:
            e["parallax"], e["parallax_error"] = plx[e["dr3"]]

    log("  " + ", ".join("%s=%d" % kv for kv in sorted(
        collections.Counter(e["kind"] for e in entries).items(),
        key=lambda kv: -kv[1])))
    return entries


def build_resolved(source, report, log=print):
    """
    Turn a ``resolve_master_list_gaia_ids.py`` report into the curated file.

    The report gives a DR3 identifier and a match strength for each zero-ID
    row; this adds the true DR2 identifier (empty where the star has no DR2
    entry) and the DR3 parallax.
    """
    names, rows = read_legacy(source)
    ni, ri, di, gi = (names.index(k) for k in ("Sp_ID", "RA", "DEC", "Gaia_ID"))

    recs = []
    with open(report, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("status") != "RESOLVED" or not r.get("gaia_dr3_id"):
                recs.append({"Sp_ID": r["Sp_ID"], "RA": r["RA"], "DEC": r["DEC"],
                             "gaia_dr2_id": "", "gaia_dr3_id": "",
                             "match_strength": "", "dr3_match": "unresolved",
                             "parallax": "", "parallax_error": "",
                             "note": (r.get("note") or "")[:200]})
                continue
            note = r.get("note", "")
            weak = ("; weak" in note) or (
                not (r.get("J_2MASS") or "").strip() and "pm=(," in note)
            recs.append({"Sp_ID": r["Sp_ID"], "RA": r["RA"], "DEC": r["DEC"],
                         "gaia_dr2_id": "", "gaia_dr3_id": r["gaia_dr3_id"],
                         "match_strength": "weak" if weak else "strong",
                         "dr3_match": "", "parallax": "", "parallax_error": "",
                         "note": note[:200]})

    dr3_ids = sorted({r["gaia_dr3_id"] for r in recs if r["gaia_dr3_id"]})
    log("%d DR3 identifications from the report" % len(dr3_ids))

    # Which of these have a DR2 counterpart?  Query the neighbourhood the
    # other way round, on dr3_source_id.
    back = collections.defaultdict(list)
    for n, total, batch in batched(dr3_ids):
        for r in tap_query(
                "SELECT dr3_source_id, dr2_source_id, angular_distance, "
                "magnitude_difference FROM gaiadr3.dr2_neighbourhood "
                "WHERE dr3_source_id IN (%s)" % ",".join(batch)):
            back[r["dr3_source_id"]].append(r)
        log("  reverse neighbourhood batch %d/%d" % (n, total))
        time.sleep(1)

    plx = parallaxes(dr3_ids, log=log)
    for r in recs:
        d3 = r["gaia_dr3_id"]
        if not d3:
            continue
        cands = back.get(d3, [])
        exact = [c for c in cands if c["dr2_source_id"] == d3]
        if exact:
            r["gaia_dr2_id"], r["dr3_match"] = d3, "identical"
        elif cands:
            best = sorted(cands, key=lambda c: (_f(c["angular_distance"]),
                                                _f(c["magnitude_difference"])))[0]
            r["gaia_dr2_id"] = best["dr2_source_id"]
            r["dr3_match"] = "differs" if len(cands) == 1 else "ambiguous"
        else:
            # No DR2 entry at all — the Wolf 359 case.
            r["gaia_dr2_id"], r["dr3_match"] = "", "no_dr2"
        if d3 in plx:
            r["parallax"], r["parallax_error"] = plx[d3]

    log("  " + ", ".join("%s=%d" % kv for kv in sorted(
        collections.Counter(r["dr3_match"] for r in recs).items(),
        key=lambda kv: -kv[1])))
    return recs


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("wrote %s (%d rows)" % (path, len(rows)))


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source", help="path to ml_40pc.txt")
    p.add_argument("--outdir", required=True,
                   help="directory to write the CSVs into (calibration/)")
    p.add_argument("--report", default=None,
                   help="report from resolve_master_list_gaia_ids.py; "
                        "without it only the DR2->DR3 map is rebuilt")
    p.add_argument("--skip-map", action="store_true",
                   help="only rebuild the curated file from --report")
    a = p.parse_args(argv)

    if not a.skip_map:
        entries = build_map(a.source)
        write_csv(os.path.join(a.outdir, "gaia_dr2_dr3_map.csv"), entries,
                  ["dr2", "dr3", "kind", "n_cands", "ang_dist_arcsec",
                   "mag_diff", "parallax", "parallax_error"])

    if a.report:
        recs = build_resolved(a.source, a.report)
        write_csv(os.path.join(a.outdir, "gaia_dr3_resolved_ids.csv"), recs,
                  ["Sp_ID", "RA", "DEC", "gaia_dr2_id", "gaia_dr3_id",
                   "match_strength", "dr3_match", "parallax",
                   "parallax_error", "note"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
