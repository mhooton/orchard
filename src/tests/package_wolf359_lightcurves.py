#!/usr/bin/env python3
"""
package_wolf359_lightcurves.py
==============================

Collect the pipeline's per-night differential light curves of one target
into a self-describing directory for sharing: one CSV per night (best
aperture), an index table, and a README explaining the columns and the
provenance.  Read-only with respect to the pipeline tree.

    python tests/package_wolf359_lightcurves.py NIGHTS_FILE OUTDIR \
        [--target Sp1056+0700] [--gaia-id 3864972938605115520] \
        [--version v3] [--base /data/SPECULOOSPipeline] [--all-apertures]

NIGHTS_FILE: one "TELESCOPE DATE" per line.  For each night the script
reads lightcurves/{gaia_id}_{filter}_{date}_bestap.txt to find the best
aperture, then converts lightcurves/{gaia_id}_{filter}_{date}_{ap}_MCMC
(whitespace table written by SPlightcurve.save_mcmc_txt) into CSV.

Columns written (from the MCMC file):
  BJD_TDB          BJDTDBMID-2450000 + 2450000, mid-exposure, TDB
  BJD_UTC          BJDMID-2450000 + 2450000, mid-exposure, UTC scale
  JD_UTC           TMID-2450000 + 2450000
  DIFF_FLUX        normalised differential flux
  ERROR            flux uncertainty
  DIFF_FLUX_PWV    PWV-corrected differential flux (NaN if not applied)
  AIRMASS, FWHM (px), SKYLEVEL (ADU), EXPOSURE (s), RA_MOVE, DEC_MOVE (px)
"""
import argparse
import csv
import datetime as dt
import glob
import os
import sys

TELESCOPE_INFO = {
    "Io": "SPECULOOS-South (Paranal, Chile), 1.0 m, Andor iKon-L 936, 0.35\"/px",
    "Europa": "SPECULOOS-South (Paranal, Chile), 1.0 m, Andor iKon-L 936, 0.35\"/px",
    "Ganymede": "SPECULOOS-South (Paranal, Chile), 1.0 m, Andor iKon-L 936, 0.35\"/px",
    "Callisto": "SPECULOOS-South (Paranal, Chile), 1.0 m, Andor iKon-L 936, 0.35\"/px",
    "Artemis": "SPECULOOS-North (Teide, Tenerife), 1.0 m, Andor iKon-L 936, 0.35\"/px",
}


def read_nights(path):
    out = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        tel, date = s.split()[:2]
        out.append((tel, date))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("nights")
    p.add_argument("outdir")
    p.add_argument("--target", default="Sp1056+0700")
    p.add_argument("--gaia-id", default="3864972938605115520")
    p.add_argument("--version", default="v3")
    p.add_argument("--base", default="/data/SPECULOOSPipeline")
    p.add_argument("--all-apertures", action="store_true", help="also copy every aperture's file")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    index = []
    missing = []
    for tel, date in read_nights(a.nights):
        lcdir = os.path.join(a.base, "PipelineOutput", a.version, tel, "output", date, a.target, "lightcurves")
        bestaps = sorted(glob.glob(os.path.join(lcdir, a.gaia_id + "_*_" + date + "_bestap.txt")))
        if not bestaps:
            missing.append((tel, date, "no bestap file in %s" % lcdir))
            continue
        for bestap in bestaps:
            base = os.path.basename(bestap)[:-len("_bestap.txt")]          # {id}_{filter}_{date}
            filt = base[len(a.gaia_id) + 1:-(len(date) + 1)]
            try:
                ap, metric = open(bestap).read().split()[:2]
            except Exception as e:
                missing.append((tel, date, "unreadable bestap: %s" % e))
                continue
            mcmc = os.path.join(lcdir, "%s_%s_MCMC" % (base, ap))
            if not os.path.exists(mcmc):
                missing.append((tel, date, "missing %s" % mcmc))
                continue
            files = [mcmc]
            if a.all_apertures:
                files = sorted(glob.glob(os.path.join(lcdir, base + "_*_MCMC")))
            for f in files:
                this_ap = f.rsplit("_", 2)[-2]
                out = os.path.join(a.outdir, "%s_%s_%s_%s_ap%s.csv" % (a.target, tel, date, filt.replace("'", ""), this_ap))
                n = convert(f, out)
                if f == mcmc:
                    index.append((tel, date, filt, this_ap, n, os.path.basename(out), metric))
    with open(os.path.join(a.outdir, "index.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["telescope", "date", "filter", "best_aperture", "n_points", "file", "bestap_metric"])
        w.writerows(index)
    with open(os.path.join(a.outdir, "README.txt"), "w") as f:
        f.write(readme(a, index, missing))
    print("%d night(s) packaged, %d missing -> %s" % (len(index), len(missing), a.outdir))
    for m in missing:
        print("  MISSING", *m)
    return 0 if index else 1


def convert(mcmc_path, out_path):
    with open(mcmc_path) as f:
        hdr = f.readline().split()
        rows = [l.split() for l in f if l.strip()]
    col = {h: i for i, h in enumerate(hdr)}

    def g(r, name, off=0.0):
        try:
            return "%.8f" % (float(r[col[name]]) + off) if off else r[col[name]]
        except (KeyError, ValueError, IndexError):
            return ""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["BJD_TDB", "BJD_UTC", "JD_UTC", "DIFF_FLUX", "ERROR", "DIFF_FLUX_PWV",
                    "AIRMASS", "FWHM_px", "SKYLEVEL_ADU", "EXPOSURE_s", "RA_MOVE_px", "DEC_MOVE_px"])
        for r in rows:
            w.writerow([g(r, "BJDTDBMID-2450000", 2450000.0), g(r, "BJDMID-2450000", 2450000.0),
                        g(r, "TMID-2450000", 2450000.0), g(r, "DIFF_FLUX"), g(r, "ERROR"),
                        g(r, "DIFF_FLUX_PWV"), g(r, "AIRMASS"), g(r, "FWHM"), g(r, "SKYLEVEL"),
                        g(r, "EXPOSURE"), g(r, "RA_MOVE"), g(r, "DEC_MOVE")])
    return len(rows)


def readme(a, index, missing):
    tels = sorted({t for t, *_ in index})
    lines = [
        "SPECULOOS differential photometry of %s (Gaia DR3 %s)" % (a.target, a.gaia_id),
        "Packaged %s from pipeline version %s" % (dt.date.today().isoformat(), a.version),
        "",
        "One CSV per night, best photometric aperture (see index.csv). Columns:",
        "  BJD_TDB        mid-exposure barycentric JD, TDB scale",
        "  BJD_UTC        mid-exposure barycentric JD, UTC scale",
        "  JD_UTC         mid-exposure JD, UTC",
        "  DIFF_FLUX      differential flux relative to a weighted artificial comparison star,",
        "                 normalised to the night's median",
        "  ERROR          photometric uncertainty on DIFF_FLUX",
        "  DIFF_FLUX_PWV  the same after the precipitable-water-vapour correction (empty when not applied)",
        "  AIRMASS, FWHM_px, SKYLEVEL_ADU, EXPOSURE_s, RA_MOVE_px, DEC_MOVE_px  per-frame diagnostics",
        "",
        "Telescopes:",
    ]
    for t in tels:
        lines.append("  %-9s %s" % (t, TELESCOPE_INFO.get(t, "")))
    lines += [
        "",
        "Filter: Sloan r' on every night. Exposure time is in the EXPOSURE_s column.",
        "Fluxes are not absolutely calibrated. Flares appear as positive excursions in DIFF_FLUX;",
        "the comparison-star weighting is optimised for the target's low-frequency systematics, so",
        "large flares are preserved but very slow trends within a night may be partly absorbed.",
        "",
        "Nights included: %d" % len(index),
    ]
    for t, d, filt, ap, n, fn, metric in index:
        lines.append("  %-9s %s %-3s ap%-2s %5d points  %s" % (t, d, filt, ap, n, fn))
    if missing:
        lines += ["", "Nights requested but not available:"]
        for m in missing:
            lines.append("  %s %s: %s" % m)
    lines += ["", "Contact: Matthew Hooton, University of Cambridge (mh2143@cam.ac.uk)", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
