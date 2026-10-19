#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
lightcurve_noise_report.py
==========================

Noise budget for the differential light curves of one target, per night.

For each night it compares the measured scatter against an analytic
prediction built from the same terms as SPlightcurve.noise_model, plus
the term that a single-star model leaves out: the photon noise of the
weighted artificial comparison star (ALC).  That term matters whenever
the target is among the brightest stars in its field, because then the
comparison ensemble is fainter than the target and sets the noise floor.

Measured scatter is reported two ways:
  * point-to-point (p2p): 1.4826 * MAD of successive differences / sqrt(2).
    Robust to flares and to slow trends, so this is the number to compare
    against the photon budget for an active star.
  * RMS: plain standard deviation, which for a flare star legitimately
    includes astrophysical signal and should exceed p2p.

Predicted relative noise per frame:
    sigma^2 = 1/N_t + SUM_i w_i^2 / N_i          (target + ALC photon)
            + (npix*(ron^2 + sky*gain + dark*t)) * (1/N_t^2 + SUM w_i^2/N_i^2)
            + sigma_scint^2 * (1 + SUM w_i^2)     (Osborn et al. 2015)
where N are electrons in the aperture and w_i are the normalised weights
actually used by the pipeline (WEIGHTS_<ap>/W_TOTAL).

Scintillation is treated as uncorrelated between target and comparisons,
which is the conservative limit; over a ~12 arcmin field it is partly
common-mode, so a measured scatter slightly below prediction is expected
rather than alarming.

Usage (inside the pipeline container):
    python tests/lightcurve_noise_report.py NIGHTS_FILE [--csv out.csv]
        [--target Sp1056+0700] [--gaia-id 3864972938605115520]
        [--version v3] [--base /data/SPECULOOSPipeline] [--bin 5]
NIGHTS_FILE has one "TELESCOPE DATE" per line; nights without light
curves are skipped with a note.
"""
from __future__ import print_function

import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits

# aperture index -> multiple of rcore (SPlightcurve.ap_size_pixels)
AP_PIX = [0.5, 1. / np.sqrt(2), 1, np.sqrt(2), 2, 2 * np.sqrt(2), 4, 5, 6, 7, 8, 10, 12]
DARK_E_PER_S_PER_PIX = 0.02
SCINT_C = 1.27          # Osborn et al. 2015 median coefficient
DEFAULT_GAIN = {'Artemis': 1.1, 'SAINT-EX': 3.48}   # else 1.0029 (SSO andor)


def robust_sigma(x):
    """1.4826 * median absolute deviation."""
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def p2p_sigma(y):
    """Point-to-point scatter: robust sigma of successive differences / sqrt(2)."""
    y = np.asarray(y, float)
    d = np.diff(y[np.isfinite(y)])
    return robust_sigma(d) / np.sqrt(2.0)


def binned_sigma(t, y, bin_minutes):
    """Robust scatter of the light curve binned to bin_minutes."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]
    if t.size < 10:
        return np.nan, 0
    edges = np.arange(t.min(), t.max() + bin_minutes / 1440.0, bin_minutes / 1440.0)
    idx = np.digitize(t, edges)
    means = np.array([np.mean(y[idx == k]) for k in np.unique(idx)
                      if np.sum(idx == k) >= 3])
    if means.size < 3:
        return np.nan, means.size
    return robust_sigma(means), means.size


def col(rec, name):
    names = {c.upper(): c for c in rec.columns.names}
    return rec[names[name.upper()]] if name.upper() in names else None


def imagelist_value(il, name, default=np.nan):
    v = col(il, name)
    if v is None:
        return default
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return np.median(v) if v.size else default


def analyse_night(base, version, tel, date, target, gaia_id, bin_minutes):
    tdir = os.path.join(base, 'PipelineOutput', version, tel, 'output', date, target)
    lcdir = os.path.join(tdir, 'lightcurves')
    bestaps = sorted(glob.glob(os.path.join(lcdir, gaia_id + '_*_' + date + '_bestap.txt')))
    if not bestaps:
        return {'telescope': tel, 'night': date, 'note': 'no bestap file (night not finished)'}
    try:
        ap = int(open(bestaps[0]).read().split()[0])
    except Exception as e:
        return {'telescope': tel, 'night': date, 'note': 'unreadable bestap: %s' % e}

    diffs = glob.glob(os.path.join(tdir, '%s_*_%s_%d_diff.fits' % (gaia_id, date, ap)))
    outs = glob.glob(os.path.join(tdir, '%s_*_%s_output.fits' % (gaia_id, date)))
    if not diffs or not outs:
        return {'telescope': tel, 'night': date, 'note': 'missing diff/output for aperture %d' % ap}

    r = {'telescope': tel, 'night': date, 'aperture': ap, 'note': ''}

    with fits.open(diffs[0]) as h:
        cat = h['CATALOGUE'].data
        diff_objids = [str(v).strip() for v in col(cat, 'OBJ_ID')] \
            if col(cat, 'OBJ_ID') is not None else None
        il = h['IMAGELIST'].data
        lc = np.asarray(h['LIGHTCURVE_%d' % ap].data, float)       # (time, star)
        wrec = h['WEIGHTS_%d' % ap].data
        w = np.asarray(col(wrec, 'W_TOTAL'), float)
        sat = np.asarray(col(wrec, 'SATURATED'), float)
        low = np.asarray(col(wrec, 'LOWFLUX'), float)
        tflag = col(cat, 'TARGET')
        ids3 = col(cat, 'GAIA_DR3_ID')
        ids2 = col(cat, 'GAIA_DR2_ID')
        flags = h['FLAGS'].data if 'FLAGS' in [x.name for x in h] else None
        bw = np.asarray(col(flags, 'BW_FLAG'), float) if flags is not None else None
        bjd = np.asarray(col(il, 'BJD-TDB'), float)
        if bjd is None or not np.isfinite(bjd).any():
            bjd = np.asarray(col(il, 'BJD-OBS'), float)
        exptime = imagelist_value(il, 'EXPOSURE', 10.0)
        airmass = imagelist_value(il, 'AIRMASS', 1.3)
        skylevel = imagelist_value(il, 'SKYLEVEL', np.nan)
        ron = imagelist_value(il, 'RON', 6.3)
        rcore = imagelist_value(il, 'RCORE', 4.0)
        gain = imagelist_value(il, 'GAIN', DEFAULT_GAIN.get(tel, 1.0029))
        altitude = imagelist_value(il, 'ALTITUDE', np.nan)
        diameter = imagelist_value(il, 'DIAMETER', 100.0)

    # target row
    ti = None
    if tflag is not None and np.nansum(np.asarray(tflag, float) > 0):
        ti = int(np.where(np.asarray(tflag, float) > 0)[0][0])
    else:
        for ids in (ids3, ids2):
            if ids is not None:
                hit = [k for k, v in enumerate(ids) if str(v).strip() == gaia_id]
                if hit:
                    ti = hit[0]
                    break
    if ti is None:
        return {'telescope': tel, 'night': date, 'note': 'target row not found in diff catalogue'}

    # The diff and output catalogues are NOT in the same row order (the
    # differential stage puts the target first), so map rows by OBJ_ID
    # rather than reusing an index across files.
    with fits.open(outs[0]) as h:
        ocat = h['CATALOGUE'].data
        out_objids = [str(v).strip() for v in col(ocat, 'OBJ_ID')] \
            if col(ocat, 'OBJ_ID') is not None else None
        flux = np.asarray(h['FLUX_%d' % ap].data, float)           # (star, time) ADU
        peak = np.asarray(h['PEAK'].data, float)
    n_out = flux.shape[0]
    if out_objids is not None and diff_objids is not None:
        pos = {oid: k for k, oid in enumerate(out_objids)}
        d2o = np.array([pos.get(oid, -1) for oid in diff_objids])
    else:
        r['note'] = 'OBJ_ID missing; assuming identical row order'
        d2o = np.arange(len(w))
    if d2o[ti] < 0:
        return {'telescope': tel, 'night': date,
                'note': 'target row not matched between diff and output catalogues'}
    # reorder the photometry into diff-catalogue order
    def to_diff_order(arr2d):
        out = np.full((len(d2o),) + arr2d.shape[1:], np.nan, float)
        ok = d2o >= 0
        out[ok] = arr2d[d2o[ok]]
        return out
    if flux.shape[0] != len(d2o) or not np.array_equal(d2o, np.arange(len(d2o))):
        flux = to_diff_order(flux)
        peak = to_diff_order(peak)
    r['n_unmatched_rows'] = int(np.sum(d2o < 0))

    # ---- comparison ensemble -------------------------------------------------
    med_flux_chk = np.nanmedian(flux, axis=1)
    wn = np.where(np.isfinite(w) & (w > 0) & np.isfinite(med_flux_chk)
                  & (med_flux_chk > 0), w, 0.0)
    wn[ti] = 0.0
    if wn.sum() <= 0:
        return {'telescope': tel, 'night': date, 'note': 'no comparison stars with positive weight'}
    wn = wn / wn.sum()
    used = np.where(wn > 0)[0]
    r['n_comp'] = int(used.size)
    r['n_comp_eff'] = float(1.0 / np.sum(wn ** 2))
    r['n_saturated'] = int(np.nansum(sat > 0)) if sat is not None else -1
    r['n_lowflux'] = int(np.nansum(low > 0)) if low is not None else -1

    med_flux = np.nanmedian(flux, axis=1)                          # ADU per frame
    e_targ = med_flux[ti] * gain
    e_comp = med_flux * gain
    r['target_e_per_frame'] = float(e_targ)
    r['target_peak_median_adu'] = float(np.nanmedian(peak[ti, :]))
    with np.errstate(divide='ignore', invalid='ignore'):
        alc_e_eff = 1.0 / np.sum(wn[used] ** 2 / e_comp[used])     # effective ALC electrons
    r['alc_e_per_frame_eff'] = float(alc_e_eff)
    r['alc_over_target_flux'] = float(alc_e_eff / e_targ)
    top = used[np.argsort(-wn[used])][:3]
    r['top_comp_weights'] = ';'.join('%.3f' % wn[k] for k in top)
    r['top_comp_flux_ratio'] = ';'.join('%.2f' % (med_flux[k] / med_flux[ti]) for k in top)
    r['brightest_comp_flux_ratio'] = float(np.nanmax(med_flux[used]) / med_flux[ti])

    # ---- measured scatter ----------------------------------------------------
    y = lc[:, ti]
    good = np.isfinite(y)
    if bw is not None and bw.size == good.size:
        good &= ~(bw > 0)
    y = y / np.nanmedian(y[good])
    r['n_points'] = int(np.sum(good))
    r['n_bad_weather'] = int(np.nansum(bw > 0)) if bw is not None else 0
    r['p2p_ppt'] = 1e3 * p2p_sigma(y[good])
    r['rms_ppt'] = 1e3 * robust_sigma(y[good])
    sb, nb = binned_sigma(bjd[good], y[good], bin_minutes)
    r['binned_ppt'] = 1e3 * sb
    r['n_bins'] = int(nb)

    # ---- predicted noise -----------------------------------------------------
    npix = np.pi * (AP_PIX[ap - 1] * rcore) ** 2
    shot = 1.0 / e_targ + np.sum(wn[used] ** 2 / e_comp[used])
    detector_e2 = npix * (ron ** 2 + (skylevel * gain if np.isfinite(skylevel) else 0.0)
                          + DARK_E_PER_S_PER_PIX * exptime)
    det = detector_e2 * (1.0 / e_targ ** 2 + np.sum(wn[used] ** 2 / e_comp[used] ** 2))
    scint = (SCINT_C * np.sqrt(10e-6) * (diameter / 100.0) ** (-2. / 3.) * airmass ** 1.75
             * (np.exp(-altitude / 8000.0) if np.isfinite(altitude) else 1.0)
             / np.sqrt(exptime))
    scint2 = scint ** 2 * (1.0 + np.sum(wn[used] ** 2))
    pred = np.sqrt(shot + det + scint2)
    r['pred_ppt'] = 1e3 * pred
    r['pred_target_shot_ppt'] = 1e3 * np.sqrt(1.0 / e_targ)
    r['pred_alc_shot_ppt'] = 1e3 * np.sqrt(np.sum(wn[used] ** 2 / e_comp[used]))
    r['pred_detector_ppt'] = 1e3 * np.sqrt(det)
    r['pred_scint_ppt'] = 1e3 * np.sqrt(scint2)
    r['ratio_p2p_pred'] = r['p2p_ppt'] / r['pred_ppt'] if r['pred_ppt'] > 0 else np.nan
    # white-noise expectation for the binned point
    npb = max(1.0, float(r['n_points']) / max(1, r['n_bins']))
    r['binned_white_ppt'] = r['p2p_ppt'] / np.sqrt(npb)
    r['red_noise_factor'] = (r['binned_ppt'] / r['binned_white_ppt']
                             if np.isfinite(r['binned_ppt']) and r['binned_white_ppt'] > 0 else np.nan)
    r['airmass_median'] = float(airmass)
    r['exposure_s'] = float(exptime)
    return r


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('nights')
    p.add_argument('--csv', default=None)
    p.add_argument('--target', default='Sp1056+0700')
    p.add_argument('--gaia-id', default='3864972938605115520')
    p.add_argument('--version', default='v3')
    p.add_argument('--base', default='/data/SPECULOOSPipeline')
    p.add_argument('--bin', type=float, default=5.0, help='bin size in minutes')
    a = p.parse_args()

    rows = []
    for line in open(a.nights):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        tel, date = s.split()[:2]
        try:
            rows.append(analyse_night(a.base, a.version, tel, date, a.target,
                                      a.gaia_id, a.bin))
        except Exception as e:  # noqa: BLE001
            rows.append({'telescope': tel, 'night': date, 'note': 'ERROR %s' % e})

    done = [r for r in rows if r.get('p2p_ppt') is not None and np.isfinite(r.get('p2p_ppt', np.nan))]
    print('%-9s %-9s %2s %5s %6s %6s %7s %7s %6s  %5s %5s %6s %6s  %s' % (
        'tel', 'night', 'ap', 'npts', 'p2p', 'rms', 'pred', 'p2p/prd', 'bin%g' % a.bin,
        'ncmp', 'neff', 'alc/t', 'shot', 'note'))
    for r in rows:
        if not np.isfinite(r.get('p2p_ppt', np.nan)):
            print('%-9s %-9s %s' % (r['telescope'], r['night'], r.get('note', '')))
            continue
        print('%-9s %-9s %2d %5d %6.2f %6.2f %7.2f %7.2f %6.2f  %5d %5.1f %6.2f %6.2f  %s' % (
            r['telescope'], r['night'], r['aperture'], r['n_points'],
            r['p2p_ppt'], r['rms_ppt'], r['pred_ppt'], r['ratio_p2p_pred'],
            r['binned_ppt'], r['n_comp'], r['n_comp_eff'], r['alc_over_target_flux'],
            r['pred_alc_shot_ppt'], r.get('note', '')))

    if done:
        rat = np.array([r['ratio_p2p_pred'] for r in done])
        print('\n%d night(s) measured. p2p/predicted: median %.2f, range %.2f-%.2f' % (
            len(done), np.median(rat), rat.min(), rat.max()))
        print('scatter units are parts per thousand (ppt); 1 ppt = 1.086 mmag')
        alc = np.array([r['pred_alc_shot_ppt'] for r in done])
        tgt = np.array([r['pred_target_shot_ppt'] for r in done])
        print('ALC shot noise exceeds target shot noise on %d/%d night(s); '
              'median ALC/target shot ratio %.1f' % (np.sum(alc > tgt), len(done),
                                                     np.median(alc / tgt)))

    if a.csv and rows:
        keys = sorted({k for r in rows for k in r})
        order = ([k for k in ('telescope', 'night', 'aperture', 'n_points', 'p2p_ppt',
                              'rms_ppt', 'pred_ppt', 'ratio_p2p_pred', 'binned_ppt',
                              'red_noise_factor', 'n_comp', 'n_comp_eff',
                              'alc_over_target_flux', 'pred_target_shot_ppt',
                              'pred_alc_shot_ppt', 'pred_detector_ppt',
                              'pred_scint_ppt') if k in keys]
                 + [k for k in keys if k not in
                    ('telescope', 'night', 'aperture', 'n_points', 'p2p_ppt', 'rms_ppt',
                     'pred_ppt', 'ratio_p2p_pred', 'binned_ppt', 'red_noise_factor',
                     'n_comp', 'n_comp_eff', 'alc_over_target_flux',
                     'pred_target_shot_ppt', 'pred_alc_shot_ppt', 'pred_detector_ppt',
                     'pred_scint_ppt')])
        with open(a.csv, 'w') as f:
            w = csv.DictWriter(f, fieldnames=order)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, '') for k in order})
        print('\nwritten: %s' % a.csv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
