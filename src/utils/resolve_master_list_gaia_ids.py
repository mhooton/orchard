#!/usr/bin/env python3
"""
resolve_master_list_gaia_ids.py
===============================

Resolve Gaia DR3 source IDs for rows of the SPECULOOS 40 pc master list
(ml_40pc.txt) whose Gaia_ID column is zero or blank.

Why
---
The pipeline identifies targets by Gaia ID and by nothing else.  A master
list row with Gaia_ID = 0 can never become a primary or secondary target,
even if the star is the brightest thing in the field.  59 of the 14,168
rows are in that state.  Most are ultracool dwarfs; a third are bright
enough that they are certainly in Gaia DR3 and were simply never resolved,
usually because their proper motion moved them away from the listed
position.

How
---
For each zero-ID row, query the Gaia archive for DR3 sources within
--radius arcmin of the listed position, joined to 2MASS for J.  Every
candidate's DR3 position (epoch J2016.0) is propagated back to the
assumed epoch of the master-list position (--list-epoch, default 2000.0,
which is what the Wolf 359 row matches) and the nearest candidate is
accepted if it lies within --match arcsec AND its 2MASS J agrees with the
master-list J to within --jtol mag.  A J disagreement is treated as a
different star, which protects against picking a neighbour.

Rows that resolve are written back with the DR3 source_id in the Gaia_ID
column.  Rows that do not resolve are listed in the report as candidates
for the supplementary source table (they are probably not in Gaia at all).

Nothing is modified unless --write is given; the original file is copied
to <file>.bak-<date> first.

Usage
-----
    python3 resolve_master_list_gaia_ids.py /path/to/ml_40pc.txt --report resolved.csv
    python3 resolve_master_list_gaia_ids.py /path/to/ml_40pc.txt --report resolved.csv --write
"""

import argparse
import csv
import datetime as dt
import io
import math
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request

TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def tap_query(adql, timeout=120, retries=3):
    url = TAP + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql})
    last = None
    for attempt in range(retries):
        try:
            body = urllib.request.urlopen(url, timeout=timeout).read().decode()
            return list(csv.DictReader(io.StringIO(body)))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("Gaia TAP query failed: %s" % last)


def cone_candidates(ra, dec, radius_arcmin):
    """DR3 sources in a cone, with 2MASS J via Gaia's own best-neighbour table."""
    adql = (
        "SELECT g.source_id, g.ra, g.dec, g.pmra, g.pmdec, g.parallax, "
        "g.phot_g_mean_mag, g.bp_rp, g.g_rp, g.teff_gspphot, t.j_m "
        "FROM gaiadr3.gaia_source AS g "
        "LEFT JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xm ON xm.source_id = g.source_id "
        "LEFT JOIN gaiadr1.tmass_original_valid AS t ON t.designation = xm.original_ext_source_id "
        "WHERE 1=CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.7f, %.7f, %.6f))"
        % (ra, dec, radius_arcmin / 60.0)
    )
    return tap_query(adql)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def separation_arcsec(ra1, dec1, ra2, dec2):
    return 3600.0 * math.hypot((ra1 - ra2) * math.cos(math.radians(dec1)), dec1 - dec2)


def propagate(ra, dec, pmra, pmdec, from_epoch, to_epoch):
    """Linear proper-motion propagation; pmra is mu_alpha* in mas/yr."""
    if math.isnan(pmra) or math.isnan(pmdec):
        return ra, dec
    dt_yr = to_epoch - from_epoch
    return (ra + pmra * 1e-3 * dt_yr / 3600.0 / math.cos(math.radians(dec)),
            dec + pmdec * 1e-3 * dt_yr / 3600.0)


def read_master_list(path):
    """Header is comma-separated, rows are whitespace-separated."""
    with open(path) as f:
        raw = f.read().replace("\r\n", "\n").replace("\r", "\n")
    lines = [l for l in raw.split("\n")]
    header_line = next(l for l in lines if l.strip())
    header = [h.strip() for h in re.split(r"[,\s]+", header_line.strip()) if h.strip()]
    rows = []
    for l in lines[lines.index(header_line) + 1:]:
        if l.strip():
            rows.append(tokenise(l))
    return header_line, header, rows, lines


def tokenise(line):
    """Rows are whitespace-separated on the server copy and comma-separated
    in some local exports; accept either."""
    return [t for t in re.split(r"[,\s]+", line.strip()) if t != ""]


def is_zero_id(value):
    v = value.strip().lower()
    return v in ("", "nan", "--", "none") or set(v) <= set("0")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("master_list")
    p.add_argument("--report", required=True, help="CSV report of every zero-ID row")
    p.add_argument("--write", action="store_true", help="write resolved IDs back to the file")
    p.add_argument("--radius", type=float, default=3.0, help="search radius, arcmin")
    p.add_argument("--match", type=float, default=3.0, help="acceptance radius, arcsec")
    p.add_argument("--jtol", type=float, default=0.5, help="allowed |J_list - J_2MASS|, mag")
    p.add_argument("--list-epoch", type=float, default=2000.0,
                   help="assumed epoch of the master-list positions")
    a = p.parse_args()

    header_line, header, rows, lines = read_master_list(a.master_list)
    gi, ri, di, ji = (header.index(k) for k in ("Gaia_ID", "RA", "DEC", "J"))
    ni, ti = header.index("Sp_ID"), header.index("T_eff")

    targets = [(k, r) for k, r in enumerate(rows) if len(r) > gi and is_zero_id(r[gi])]
    print("%d rows, %d with no Gaia ID" % (len(rows), len(targets)), flush=True)

    report = []
    resolved = {}
    for n, (k, r) in enumerate(targets, 1):
        name, ra, dec, jlist = r[ni], _f(r[ri]), _f(r[di]), _f(r[ji])
        try:
            cands = cone_candidates(ra, dec, a.radius)
        except Exception as e:  # noqa: BLE001
            report.append([name, ra, dec, jlist, r[ti], "QUERY_FAILED", "", "", "", "", str(e)])
            print("%3d/%d %-14s query failed: %s" % (n, len(targets), name, e), flush=True)
            continue
        best = None
        for c in cands:
            cra, cdec = _f(c["ra"]), _f(c["dec"])
            pra, pdec = propagate(cra, cdec, _f(c["pmra"]), _f(c["pmdec"]), 2016.0, a.list_epoch)
            sep_list_epoch = separation_arcsec(pra, pdec, ra, dec)
            sep_2016 = separation_arcsec(cra, cdec, ra, dec)
            sep = min(sep_list_epoch, sep_2016)
            jm = _f(c["j_m"])
            jdiff = abs(jm - jlist) if not (math.isnan(jm) or math.isnan(jlist)) else float("nan")
            ok = sep <= a.match and (math.isnan(jdiff) or jdiff <= a.jtol)
            if ok and (best is None or sep < best[0]):
                best = (sep, c, jdiff, "epoch%.0f" % a.list_epoch if sep_list_epoch <= sep_2016 else "epoch2016")
        if best is None:
            near = sorted(((min(separation_arcsec(
                *propagate(_f(c["ra"]), _f(c["dec"]), _f(c["pmra"]), _f(c["pmdec"]), 2016.0, a.list_epoch),
                ra, dec), separation_arcsec(_f(c["ra"]), _f(c["dec"]), ra, dec)), c) for c in cands),
                key=lambda x: x[0])[:1]
            note = ""
            if near:
                s0, c0 = near[0]
                note = "nearest DR3 %s at %.1f arcsec G=%s J=%s" % (
                    c0["source_id"], s0, c0["phot_g_mean_mag"], c0["j_m"])
            report.append([name, ra, dec, jlist, r[ti], "UNRESOLVED", "", "", "", "", note])
            print("%3d/%d %-14s UNRESOLVED (%d candidates) %s" % (n, len(targets), name, len(cands), note), flush=True)
            continue
        sep, c, jdiff, how = best
        resolved[k] = c["source_id"]
        report.append([name, ra, dec, jlist, r[ti], "RESOLVED", c["source_id"], "%.2f" % sep,
                       c["phot_g_mean_mag"], c["j_m"], "%s; pm=(%s,%s) mas/yr; Jdiff=%.2f" % (
                           how, c["pmra"], c["pmdec"], jdiff if not math.isnan(jdiff) else -1)])
        print("%3d/%d %-14s -> %s  sep %.2f\"  G=%s J=%s (list J=%.2f) pm=(%s,%s)" % (
            n, len(targets), name, c["source_id"], sep, c["phot_g_mean_mag"], c["j_m"], jlist,
            c["pmra"], c["pmdec"]), flush=True)
        time.sleep(0.5)

    with open(a.report, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Sp_ID", "RA", "DEC", "J_list", "T_eff", "status", "gaia_dr3_id",
                    "sep_arcsec", "G", "J_2MASS", "note"])
        w.writerows(report)
    n_res = sum(1 for r in report if r[5] == "RESOLVED")
    print("resolved %d, unresolved %d, failed %d -> %s" % (
        n_res, sum(1 for r in report if r[5] == "UNRESOLVED"),
        sum(1 for r in report if r[5] == "QUERY_FAILED"), a.report), flush=True)

    if a.write and resolved:
        backup = "%s.bak-%s" % (a.master_list, dt.date.today().isoformat())
        shutil.copyfile(a.master_list, backup)
        # Rewrite only the affected lines, preserving everything else byte for byte.
        data_start = lines.index(header_line) + 1
        row_line_numbers = [i for i in range(data_start, len(lines)) if lines[i].strip()]
        out = list(lines)
        for k, sid in resolved.items():
            ln = row_line_numbers[k]
            old = rows[k][gi]
            # Replace the first whitespace-delimited occurrence of the old ID token.
            out[ln] = re.sub(r"(?<=[\s,])%s(?=[\s,])" % re.escape(old), sid, out[ln], count=1)
        with open(a.master_list, "w") as f:
            f.write("\n".join(out))
        print("wrote %d IDs into %s (backup: %s)" % (len(resolved), a.master_list, backup))


if __name__ == "__main__":
    sys.exit(main())
