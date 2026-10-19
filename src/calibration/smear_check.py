"""
Measure readout smear on raw frames: was the shutter used, and what t_read?

For every bright, unsaturated, isolated star, the readout trail adds tau times
the star's flux in each row to every pixel downstream of it, and nothing
upstream. So the downstream-minus-upstream level in the star's rows, less the
same difference in nearby star-free rows (which removes the sky-smear ramp and
column structure), is tau * F_star. Summed over stars:

    tau_hat = sum(E_i) / sum(F_i),     t_read_hat = tau_hat * N * t_exp

A shuttered frame gives tau_hat consistent with zero. With --desmear the frame
is also corrected with the configured t_read and re-measured; the removed
fraction should be consistent with 100%, and "100% at" is the t_read that
would remove it exactly. That is the better t_read estimate: tau_hat itself
runs high when the flux box misses the PSF wings.

An asymmetric PSF along the readout axis fakes a small trail. On 2026-10-12
shuttered controls gave +0.3 s (Ganymede, 20 s), +0.4 s (Ganymede, 10 s) and
-1.5 s (Io, 39 s) of apparent t_read with the default windows, against
4.2-4.9 s on the shutterless Ganymede nights of 2025-04-08/09 and 2026-02-21.
That is ample to tell shutterless from shuttered, but a t_read measured to
better than ~0.3 s needs a shuttered control night of the same telescope,
target and exposure time.

Flats (IMAGETYP containing 'flat') have no stars; for them the readout ramp
across the frame is measured instead (ramp_t_read). Check the flats either side
of a shutter switch individually: de-smearing a shuttered 7 s flat by mistake
puts a ~60% ramp into the master flat.

Use it before adding a window to shutter_overrides.csv, and to check a night
whose CAMSHUTT header might be wrong.

Usage:
    python calibration/smear_check.py FRAME [FRAME ...] [--desmear]
    python calibration/smear_check.py @framelist.txt --desmear
"""

import argparse
import os
import sys

import numpy as np
from scipy.ndimage import maximum_filter
from astropy.stats import sigma_clip

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration.pipeutils import apply_image_slice, get_instrument_parameters, open_fits_file
from calibration.smear import _along_axis, desmear, mode_readout_time, readout_time, shutter_state, smear_tau

BOX = 10            # star flux box half-width (pixels, both axes)
GAP = 40            # windows start this far from the star along the readout axis,
REACH = 250         # beyond the PSF wings, and end this far
REF_ROWS = (16, 24)  # reference rows sit this far from the star across the axis
ISOLATION = 0.1     # reject a star with a neighbour this bright (peak ratio) near its windows


def load_signal(path, reduced=False):
    """
    Trimmed, overscan-subtracted frame plus what is needed to measure it.

    reduced=True takes a T6 product (proc*.fits) as it is: already trimmed,
    calibrated and, if it was shutterless, de-smeared, so what is left of the
    trail there is what the pipeline failed to remove.
    """
    with open_fits_file(path) as hdul:
        params = get_instrument_parameters(hdul)
        header = hdul[0].header.copy()
        if reduced:
            data = np.array(hdul[0].data, dtype=np.float64)
            return data, data, header, params
        overscan = apply_image_slice(hdul, 'overscan', 0.)
        if isinstance(overscan, np.ndarray) and overscan.size > 1:
            overscan = float(np.ma.median(sigma_clip(overscan, maxiters=None)))
        raw = np.array(apply_image_slice(hdul, 'trim', hdul[0].data), dtype=np.float64)
    return raw - overscan, raw, header, params


def find_stars(s, raw, saturation):
    """Bright unsaturated peaks, and every peak that could pollute a window."""
    bg = np.median(s)
    smoothed = maximum_filter(s, size=7)
    peaks = (s == smoothed)
    faint = np.argwhere(peaks & (s - bg > 300))
    bright = np.argwhere(peaks & (s - bg > 3000))
    stars = []
    nrow, ncol = s.shape
    margin_r = BOX + REF_ROWS[1] + 1
    margin_c = REACH + BOX + 1
    for y, k in bright:
        if not (margin_r <= y < nrow - margin_r and margin_c <= k < ncol - margin_c):
            continue
        if raw[y - BOX:y + BOX + 1, k - BOX:k + BOX + 1].max() >= 0.95 * saturation:
            continue
        near = faint[(np.abs(faint[:, 0] - y) <= REF_ROWS[1] + 2)
                     & (np.abs(faint[:, 1] - k) <= REACH + BOX)]
        others = [(yy, kk) for yy, kk in near if (yy, kk) != (y, k)
                  and s[yy, kk] - bg > ISOLATION * (s[y, k] - bg)]
        if others:
            continue
        stars.append((int(y), int(k)))
    return stars


def window_difference(s, rows, k):
    """Per-row median downstream minus median upstream, along the readout axis."""
    up = np.median(s[rows, k - REACH:k - GAP], axis=1)
    down = np.median(s[rows, k + GAP + 1:k + REACH + 1], axis=1)
    return down - up


def measure(s, stars):
    """tau_hat, its uncertainty, and the per-star (E, F) pairs."""
    pairs = []
    for y, k in stars:
        star_rows = np.arange(y - BOX, y + BOX + 1)
        ref_rows = np.concatenate([np.arange(y - REF_ROWS[1], y - REF_ROWS[0]),
                                   np.arange(y + REF_ROWS[0] + 1, y + REF_ROWS[1] + 1)])
        d_ref = window_difference(s, ref_rows, k).mean()
        excess = (window_difference(s, star_rows, k) - d_ref).sum()
        bg = np.median(s[ref_rows, k - BOX:k + BOX + 1])
        flux = (s[y - BOX:y + BOX + 1, k - BOX:k + BOX + 1] - bg).sum()
        if flux > 0:
            pairs.append((excess, flux))
    tau, err = combine(pairs)
    return tau, err, pairs


def combine(pairs):
    """Flux-weighted tau from (E, F) pairs, with a scatter-based uncertainty."""
    if len(pairs) < 3:
        return np.nan, np.nan
    e, f = np.array(pairs).T
    tau = e.sum() / f.sum()
    err = np.sqrt(((e - tau * f) ** 2).sum() * len(e) / (len(e) - 1)) / f.sum()
    return tau, err


def verdict(t_read, t_err):
    """Shuttered only if the measurement rules out any plausible readout time."""
    if not np.isfinite(t_read):
        return 'too few stars'
    if t_read > 5 * t_err and t_read > 2.0:
        return 'SHUTTERLESS'
    if t_read + 3 * t_err < 2.0:
        return 'shuttered'
    return 'inconclusive'


def describe(label, pairs, pairs_after, n, exptime, readout_time):
    """One result line from (E, F) pairs measured before and after de-smearing."""
    tau, err = combine(pairs)
    t_read, t_err = tau * n * exptime, err * n * exptime
    line = f'{label}  stars={len(pairs):4d}  t_read={t_read:6.2f} +/- {t_err:4.2f} s  {verdict(t_read, t_err)}'
    if pairs_after and verdict(t_read, t_err) == 'SHUTTERLESS':
        tau_after, _ = combine(pairs_after)
        removed = 1 - tau_after / tau
        line += (f'  removed {100 * removed:6.1f}% at t_read={readout_time} s'
                 f' (100% at {readout_time / removed:.2f} s)')
    return line


def ramp_t_read(oriented, n, exptime, edge=200):
    """
    Readout time implied by the ramp across a flat.

    Uniform light smeared on readout gives raw(x) = S * (1 + tau * x) along the
    readout axis, so two bands near either end fix tau. A twilight gradient adds
    an apparent t_read proportional to exposure time, with opposite effect in
    the evening and the morning: on Ganymede in February 2026 6 s evening flats
    gave 4.27-4.31 s and 8 s morning flats 4.57 s, either side of the 4.33 s
    from star trails (the 4.58 s from a flat in August 2026 was a morning one).
    """
    k = oriented.shape[-1]
    left = float(np.median(oriented[:, edge:2 * edge]))
    right = float(np.median(oriented[:, k - 2 * edge:k - edge]))
    x_left, x_right = 1.5 * edge, k - 1.5 * edge
    tau = (right - left) / (left * x_right - right * x_left)
    return tau * n * exptime


def flat_verdict(t_read, exptime, n=2048, span=1446):
    """
    Judge on the fractional ramp the readout would cause, not on t_read alone:
    a twilight gradient mimics a t_read that grows with exposure time, so a long
    flat can neither show a 1.5 s (3 MHz) readout nor rule one out.
    span is the distance between the two bands ramp_t_read() compares.
    """
    ramp = t_read * span / (n * exptime)
    if t_read > 0.75 and ramp > 0.05:
        return 'SHUTTERLESS'
    could_see_fast_readout = 1.4 * span / (n * exptime) > 0.03
    if abs(t_read) < 0.5 and could_see_fast_readout:
        return 'shuttered'
    return 'inconclusive'


def check(path, apply_desmear, reduced=False):
    """Measure one frame. Returns (pairs, pairs_after, n_shifts, exptime, readout_time, header)."""
    s, raw, header, params = load_signal(path, reduced)
    smear = params.get('smear')
    if not smear:
        raise ValueError('no smear parameters for this instrument')
    axis, direction = smear.get('axis', 1), smear.get('direction', 1)
    n, exptime = smear['n_shifts'], float(header['EXPTIME'])
    # the t_read the pipeline would use: an override window's, else the
    # frame's READMODE (T6 products), else the normal 1 MHz readout
    t_read = readout_time(header, smear, shutter_state(header)[2])[0] or mode_readout_time(smear, 'NORMAL')
    oriented = _along_axis(s, axis, direction)
    if 'flat' in str(header.get('IMAGETYP', '')).lower():
        return ramp_t_read(oriented, n, exptime), None, n, exptime, t_read, header
    stars = find_stars(oriented, _along_axis(raw, axis, direction), params['saturation_threshold'])
    _, _, pairs = measure(oriented, stars)
    pairs_after = []
    if apply_desmear and pairs:
        corrected = _along_axis(desmear(s, smear_tau(t_read, n, exptime),
                                        axis=axis, direction=direction), axis, direction)
        _, _, pairs_after = measure(corrected, stars)
    return pairs, pairs_after, n, exptime, t_read, header


def main():
    global GAP, REACH
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('frames', nargs='+', help='raw FITS frames, or @file listing them')
    parser.add_argument('--desmear', action='store_true',
                        help='also de-smear with the configured t_read and report the fraction removed')
    parser.add_argument('--reduced', action='store_true',
                        help='frames are T6 products (proc*.fits): measure the trail left after reduction')
    parser.add_argument('--gap', type=int, default=GAP,
                        help=f'trail windows start this many pixels from the star (default {GAP}); '
                             'smaller keeps more stars but lets PSF asymmetry bias t_read')
    parser.add_argument('--reach', type=int, default=REACH,
                        help=f'trail windows end this many pixels from the star (default {REACH})')
    args = parser.parse_args()
    GAP, REACH = args.gap, args.reach
    frames = []
    for item in args.frames:
        if item.startswith('@'):
            with open(item[1:]) as f:
                frames += [line.strip() for line in f if line.strip()]
        else:
            frames.append(item)
    groups = {}
    for path in frames:
        try:
            pairs, after, n, exptime, readout_time, header = check(path, args.desmear, args.reduced)
        except Exception as e:
            print(f'{os.path.basename(path)}: ERROR {e}', flush=True)
            continue
        label = (f'{os.path.basename(path)}  CAMSHUTT={header.get("CAMSHUTT", "-")!s:5} '
                 f'exp={exptime:6.1f}s')
        if after is None:
            # a flat: pairs is the t_read implied by its readout ramp
            print(f'{label}  FLAT  ramp t_read={pairs:6.2f} s  {flat_verdict(pairs, exptime)}', flush=True)
            continue
        print(describe(label, pairs, after, n, exptime, readout_time), flush=True)
        key = (str(header.get('TELESCOP', '')).strip(), exptime)
        g = groups.setdefault(key, {'pairs': [], 'after': [], 'frames': 0, 'n': n, 'rt': readout_time})
        g['pairs'] += pairs
        g['after'] += after
        g['frames'] += 1
    if len(frames) > 1:
        print('\nStacked over frames:')
        for (telescope, exptime), g in groups.items():
            label = f'{telescope}  exp={exptime:6.1f}s  frames={g["frames"]:3d}'
            print(describe(label, g['pairs'], g['after'], g['n'], exptime, g['rt']))


if __name__ == '__main__':
    main()
