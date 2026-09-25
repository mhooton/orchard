#!/usr/bin/env python3
"""
verify_wolf359_night.py
=======================

Check that one night of Sp1056+0700 (Wolf 359) reprocessed on the
missing-gaia-targets branch identified the right star and produced light
curves.  Read-only.  Run inside the pipeline container:

    python tests/verify_wolf359_night.py Europa 20220105 [--version v3] [--base /data/SPECULOOSPipeline]

Checks, in order:
  1. stack catalogue exists and is named by Wolf 359's DR3 ID
  2. Gaia_Crossmatch header: PRIM_ST / PRIM_HOW / PRIM_ID / SCHED_ID / PRIMINJ
  3. the TARGET_ROLE == 'primary' row is Wolf 359: within 2" of the
     proper-motion-propagated DR3 position, and the brightest detection
  4. photometry output ({ID}_r_{date}_output.fits) exists, and its
     CATALOGUE carries the ID and non-sky PEAK counts on that row
  5. differential products: *_diff.fits, lightcurves/*_MCMC, *_bestap.txt
Exit code 0 when every check passes, 1 otherwise.
"""
import argparse
import glob
import math
import os
import sys

import numpy as np
from astropy.io import fits

WOLF = "3864972938605115520"
RA0, DEC0, PMRA, PMDEC = 164.10319030755974, 7.002726940984864, -3866.338, -2699.215


def expected(date):
    y, m, d = int(date[:4]), int(date[4:6]), int(date[6:])
    epoch = y + (m - 1) / 12.0 + (d - 1) / 365.0
    dt = epoch - 2016.0
    return (RA0 + PMRA * 1e-3 * dt / 3600.0 / math.cos(math.radians(DEC0)),
            DEC0 + PMDEC * 1e-3 * dt / 3600.0)


def sep(ra1, dec1, ra2, dec2):
    return 3600.0 * math.hypot((ra1 - ra2) * math.cos(math.radians(dec1)), dec1 - dec2)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("telescope")
    p.add_argument("date")
    p.add_argument("--version", default="v3")
    p.add_argument("--base", default="/data/SPECULOOSPipeline")
    p.add_argument("--target", default="Sp1056+0700")
    a = p.parse_args()
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + msg)
        ok = ok and bool(cond)

    tdir = os.path.join(a.base, "PipelineOutput", a.version, a.telescope, "output", a.date, a.target)
    print("target dir:", tdir)
    check(os.path.isdir(tdir), "target directory exists")
    if not os.path.isdir(tdir):
        return 1

    # 1. stack catalogue
    cats = sorted(glob.glob(os.path.join(tdir, "*_stack_catalogue_*.fits")))
    print("stack catalogues:", [os.path.basename(c) for c in cats])
    check(cats, "stack catalogue present")
    wolf_cats = [c for c in cats if os.path.basename(c).startswith(WOLF + "_")]
    check(wolf_cats, "stack catalogue named by Wolf 359 DR3 ID")
    cat = wolf_cats[0] if wolf_cats else (cats[0] if cats else None)
    era, edec = expected(a.date)
    if cat:
        with fits.open(cat) as h:
            g = h["GAIA_CROSSMATCH"]
            hdr = g.header
            for k in ("PRIM_ST", "PRIM_HOW", "PRIM_ID", "SCHED_ID", "PRIMINJ", "N_MATCH"):
                print("   %-8s = %s" % (k, hdr.get(k)))
            check(str(hdr.get("PRIM_ID", "")).strip() == WOLF, "PRIM_ID is Wolf 359")
            check(str(hdr.get("PRIM_ST", "")).startswith("OK"), "primary status OK (%s)" % hdr.get("PRIM_ST"))
            roles = np.array([str(x).strip() for x in g.data["TARGET_ROLE"]]) if "TARGET_ROLE" in g.columns.names else np.array([])
            prim = np.where(roles == "primary")[0]
            check(len(prim) == 1, "exactly one primary row (%d)" % len(prim))
            ids3 = [str(x).strip() for x in g.data["GAIA_DR3_ID"]] if "GAIA_DR3_ID" in g.columns.names else []
            ids2 = [str(x).strip() for x in g.data["GAIA_DR2_ID"]]
            apm = h[1].data
            names = [n.upper() for n in apm.columns.names]
            ra = np.asarray(apm[apm.columns.names[names.index("RA")]], float)
            dec = np.asarray(apm[apm.columns.names[names.index("DEC")]], float)
            if np.nanmax(np.abs(ra)) <= 2 * np.pi + 1e-6:
                ra, dec = np.degrees(ra), np.degrees(dec)
            flux = np.asarray(apm[apm.columns.names[names.index("ISOPHOTAL_FLUX")]], float) if "ISOPHOTAL_FLUX" in names else None
            if len(prim) == 1:
                i = int(prim[0])
                s = sep(ra[i], dec[i], era, edec)
                print("   primary row %d: DR3=%s DR2=%s  sep from expected Wolf 359 position = %.2f arcsec" % (
                    i, ids3[i] if ids3 else "-", ids2[i], s))
                check(s < 2.0, "primary row within 2 arcsec of propagated DR3 position")
                check((ids3 and ids3[i] == WOLF) or ids2[i] == WOLF, "primary row carries Wolf 359 ID")
                if flux is not None:
                    rank = int(np.sum(np.nan_to_num(flux, nan=-1) > flux[i])) + 1
                    print("   primary brightness rank %d of %d (flux %.0f)" % (rank, len(flux), flux[i]))
                    check(rank == 1, "primary is the brightest detection")
                for col in ("PMRA", "PMDEC", "GMAG", "TEFF"):
                    if col in g.columns.names:
                        print("   %-5s = %s" % (col, g.data[col][i]))
                if "PMRA" in g.columns.names:
                    check(np.isfinite(float(g.data["PMRA"][i])), "primary row has proper motion")

    # 4. photometry output
    outs = sorted(glob.glob(os.path.join(tdir, "*_output.fits")))
    print("output files:", [os.path.basename(o) for o in outs])
    wolf_outs = [o for o in outs if os.path.basename(o).startswith(WOLF + "_")]
    check(wolf_outs, "photometry output named by Wolf 359 DR3 ID")
    if wolf_outs:
        with fits.open(wolf_outs[0]) as h:
            c = h["CATALOGUE"].data
            cols = [n.upper() for n in c.columns.names]
            ids = []
            for col in ("GAIA_DR3_ID", "GAIA_DR2_ID"):
                if col in cols:
                    ids.append([str(x).strip() for x in c[col]])
            rows = [i for i in range(len(c)) if any(l[i] == WOLF for l in ids)]
            check(len(rows) >= 1, "Wolf 359 ID present in output CATALOGUE (%d row(s))" % len(rows))
            if rows:
                i = rows[0]
                cra = float(c["RA"][i]); cdec = float(c["DEC"][i])
                if abs(cra) <= 2 * math.pi + 1e-6:
                    cra, cdec = math.degrees(cra), math.degrees(cdec)
                s = sep(cra, cdec, era, edec)
                pk = np.asarray(h["PEAK"].data[i, :], float)
                pk = pk[np.isfinite(pk) & (pk > 0)]
                print("   output row %d: sep %.2f arcsec, PEAK median %.0f max %.0f over %d frames" % (
                    i, s, np.median(pk) if len(pk) else -1, pk.max() if len(pk) else -1, pk.size))
                check(s < 2.0, "output row within 2 arcsec of Wolf 359")
                check(len(pk) and np.median(pk) > 1000, "aperture sits on the star (median peak > 1000 ADU)")
                check(len(pk) and pk.max() < 62000, "no saturation")

    # 5. differential products
    diffs = sorted(glob.glob(os.path.join(tdir, WOLF + "_*_diff.fits")))
    mcmc = sorted(glob.glob(os.path.join(tdir, "lightcurves", WOLF + "_*_MCMC")))
    bestap = sorted(glob.glob(os.path.join(tdir, "lightcurves", WOLF + "_*_bestap.txt")))
    print("diff.fits: %d, MCMC: %d, bestap: %s" % (len(diffs), len(mcmc), [os.path.basename(b) for b in bestap]))
    check(len(diffs) >= 1, "differential photometry files present")
    check(len(mcmc) >= 1, "MCMC light-curve text files present")
    if mcmc:
        n = sum(1 for _ in open(mcmc[-1])) - 1
        with open(mcmc[-1]) as f:
            hdr = f.readline().split()
            first = f.readline().split()
        print("   %s: %d points, columns %s, first row %s" % (os.path.basename(mcmc[-1]), n, hdr[:6], first[:5]))
        check(n > 50, "light curve has more than 50 points")
    if diffs:
        with fits.open(diffs[-1]) as h:
            c = h["CATALOGUE"].data
            t = np.asarray(c["TARGET"]) if "TARGET" in c.columns.names else None
            w = h[[x.name for x in h if x.name.startswith("WEIGHTS")][0]].data if any(x.name.startswith("WEIGHTS") for x in h) else None
            if t is not None:
                ti = np.where(t > 0)[0]
                print("   diff.fits TARGET rows:", ti.tolist())
                if len(ti) and w is not None:
                    print("   target weight/sat/lowflux:", w["W_TOTAL"][ti[0]], w["SATURATED"][ti[0]], w["LOWFLUX"][ti[0]])
                    check(w["SATURATED"][ti[0]] == 0, "target not flagged saturated")

    print("\nRESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
