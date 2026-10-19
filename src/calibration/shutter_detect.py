"""
Per-night detection of shutterless Andor frames.

CAMSHUTT is a static value in Astra's config, changed by hand when a camera is
switched, so it would not change if a shutter failed open, and
shutter_overrides.csv covers only the windows somebody has verified. This
module measures each night's shutter state from the night's own frames, so a
camera run shutterless without anyone saying so is noticed, and (once the
thresholds are agreed) corrected.

The order of evidence, first that gives an answer wins (agreed 2026-10-19):

1. shutter_overrides.csv, the human-verified windows. Applied frame by frame
   in smear.shutter_state(); the record here notes where the measurement
   disagrees with a window.
2. The night's OWN flats: the frames physically in
   Observations/<TEL>/images/<NIGHT>/, normally evening flats before the
   science and morning flats after it. Never the reduction's flat lists, which
   borrow from earlier nights: Ganymede 20250408-09 was reduced with
   20250407's shuttered morning flats while it was shutterless. Every flat
   gets the readout-ramp test (smear_check.ramp_t_read / flat_verdict). When
   the evening and the morning agree, that state holds for every frame of the
   night. The readout rate is the measured rate nearest the ramp's t_read,
   cross-checked against READOUTM (Astra, Andor frames, from ~2026-10-18).
3. Science frames stacked over the night, for a night whose own flats are
   missing, one-sided or disagree: the star-trail excess (smear_check.measure)
   for exposures up to 40 s, needing about 12 frames. On a night whose evening
   and morning flats disagree, bins of consecutive science frames locate the
   switch. For longer exposures the sky readout ramp against a shuttered
   reference flat is measured too; that test is EXPERIMENTAL and never counts
   as unambiguous.
4. CAMSHUTT, in smear.shutter_state(), when the measurement is inconclusive.
5. Shutter used.

The measurement goes into a per-night record, shutter.json in the reduction
directory (written at T2, read back at T4 and T6; a night borrowed for its
flats gets shutter_<night>.json, flats only). The record splits the night
into segments, each with a state, a readout time, the evidence and a
confidence:
    unambiguous   both twilights agree; or one twilight confirmed by the
                  science; or a split located between conclusive bins
    probable      one twilight alone, the science alone, or a READOUTM that
                  disagrees with the ramp
    inconclusive  nothing usable, or the evidence contradicts itself
Only an unambiguous segment can change what is de-smeared, and only where
the instrument's "smear" config sets "detect_auto_apply": true. Until then
the measurement is reported (SHUTD* cards, warnings, the digest) and the
decision stays override, CAMSHUTT, default. The risks are asymmetric:
de-smearing a shuttered 6-8 s flat puts a ramp of tens of per cent into it,
while leaving a shutterless night smeared costs about 1% in noise.

Measured reference numbers (Ganymede, 2025-26):
    shuttered flats, ramp t_read             0.05-0.2 s
    shutterless 1 MHz flats, evening         4.27-4.31 s
    shutterless 1 MHz flats, morning         4.57 s (twilight gradient)
    shutterless 3 MHz flats                  1.50 s
    1,653 SSO flats separated cleanly; ten long-exposure flats inconclusive.

    python -m calibration.shutter_detect <night dir> [--json] [--all-levels]
measures a night and prints the record without writing anything.
"""
import argparse
import datetime
import glob
import json
import os
import sys

import numpy as np
from astropy.stats import sigma_clip

if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration import readout, smear, smear_check
from calibration.pipeutils import apply_image_slice, get_instrument_parameters, open_fits_file, frame_telescope

RECORD_NAME = 'shutter.json'
RECORD_VERSION = 1

SHUTTERED, SHUTTERLESS, MIXED = 'shuttered', 'shutterless', 'mixed'
UNAMBIGUOUS, PROBABLE, INCONCLUSIVE = 'unambiguous', 'probable', 'inconclusive'

# ------------------------------------------------------------- thresholds
# Flats.  A flat is measured only between these levels: a faint one is noise,
# and one near saturation has its ramp compressed by the clip.
FLAT_MIN_SIGNAL = 2000.0          # ADU above the overscan
FLAT_MAX_FRACTION = 0.85          # of the readout mode's raw saturation level
MIN_SIDE_FLATS = 2                # conclusive flats a twilight needs; all must agree
# A measured t_read names a readout rate when it lies within this fraction of
# that rate's time in "readout_rate_time": 1 MHz flats gave 4.27-4.57 s
# against 4.33 s, 3 MHz flats 1.50 s against 1.5 s.
RATE_TOLERANCE = 0.15
# A twilight whose flats stray further than this from their t_read-exptime
# line mixes readout rates (1 and 3 MHz flats would stray 1.4-2.8 s). Noisy
# twilights stray up to 0.34 s (Ganymede 2026-02-18 morning, 02-24 evening).
RAMP_MIX_FRACTION = 0.25
RAMP_MIX_FLOOR = 0.1              # s
# Readout plus twilight make t_read climb at most 0.023 s per s of exposure
# (Ganymede 2026-02 mornings). Early morning twilights on shuttered nights
# climb 0.3 s per s and more (Ganymede 20260308: 28 s at 60 s, ~0 at 20 s):
# a sky gradient, not a readout, so such a twilight names no rate.
RAMP_MAX_SLOPE = 0.06             # s per s
# Flats are split into evening and morning at the first and last science
# frame, or without science at the largest gap if it is at least this long.
TWILIGHT_GAP_HOURS = 4.0

# Science star trails.
TRAIL_MAX_EXPTIME = 40.0          # s; longer exposures dilute the trail below the PSF bias
TRAIL_MIN_FRAMES = 12             # frames in an exposure-time group
TRAIL_SAMPLE = 24                 # frames measured for a whole-night verdict
TRAIL_BIN = 12                    # consecutive frames per bin when locating a switch
TRAIL_MIN_STARS = 20              # star measurements a stack needs
# PSF asymmetry along the readout axis fakes a trail, more so in longer
# exposures. Shuttered nights (validation, 2026-10-19) read -2.0 to +2.8 s:
# Ganymede +0.35..+0.61 s at 10 s, -0.03..+0.27 s at 20-25 s, +1.17..+1.38 s
# at 37 s (one field, both twilights shuttered); Io -1.97 +/- 0.17 s at 28 s,
# +2.83 +/- 1.02 s at 23 s. So the trail test can tell 1 MHz (4.33 s) from
# "not 1 MHz", but 3 MHz (1.5 s) from shuttered only when READOUTM names the
# rate. A shuttered stack would need a bias of about +3.6 s at 20 s to read
# as a conclusive 1 MHz.
TRAIL_SYS_FLOOR = 0.5             # s
TRAIL_SYS_PER_S = 0.02            # s of apparent t_read per s of exposure
# On shutterless frames the trail reads high when the PSF is broad: wings
# outside the +/-10 px flux box still feed the trail (Ganymede 2026-02-16,
# 39 s: 5.42 s against 4.33 s). So a rate is consistent with anything from
# TRAIL_RATE_RANGE[0] to TRAIL_RATE_RANGE[1] times its time.
TRAIL_RATE_RANGE = (0.9, 1.35)
CONSISTENT_SIGMA = 3.0            # a hypothesis is consistent within this
EXCLUDED_SIGMA = 3.0              # and excluded beyond it
# A night measured on its science alone is never unambiguous by default: the
# PSF bias has a long tail (Ganymede 20230130, shuttered, read -5.11 +/- 0.51 s
# at 31 s), and a positive one that size would pass as a conclusive 1 MHz
# trail and de-smear a whole shuttered night. A finite value makes a
# science-only verdict with this margin (sigma) to every other hypothesis
# unambiguous; the verified science-only nights had 4.0-4.9.
SCIENCE_ALONE_SIGMA = float('inf')

# Sky readout ramp (EXPERIMENTAL): exposures longer than TRAIL_MAX_EXPTIME.
SKY_SAMPLE = 16
SKY_MIN_FRAMES = 6
SKY_SYS_FLOOR = 0.3               # s; a 0.2% sky gradient across a 120 s frame fakes ~0.3 s
REF_MAX_DAYS = 10                 # how far back to look for a shuttered reference night
REF_MIN_FLATS = 4
REF_MAX_FLATS = 10
REF_MAX_EXPTIME = 20.0


# ------------------------------------------------------------- frames

FRAME_PATTERNS = ('*.fits', '*.fts', '*.fit', '*.fits.fz', '*.fts.fz', '*.fits.bz2')


def night_files(night_dir):
    """Every frame in a night directory and its subdirectories one level down,
    as createlists.get_liste finds them."""
    files = set()
    for depth in ('', '*'):
        for pattern in FRAME_PATTERNS:
            files.update(glob.glob(os.path.join(night_dir, depth, pattern)))
    # a compressed copy beside its uncompressed frame is the same frame
    return sorted(f for f in files if not (f.endswith(('.fz', '.bz2')) and f.rsplit('.', 1)[0] in files))


def frame_kind(imagetyp):
    """createlists.classify_file's classes, in its order."""
    text = str(imagetyp or '').upper()
    for key, kind in (('LIGHT', 'science'), ('DARK', 'dark'), ('BIAS', 'bias'), ('FLAT', 'flat')):
        if key in text:
            return kind
    return 'other'


def _header(hdul):
    header = hdul[0].header
    if header.get('NAXIS', 0) == 0 and len(hdul) > 1:
        return hdul[1].header
    return header


def _iso(value):
    try:
        return smear._parse_utc(str(value)).isoformat()
    except (TypeError, ValueError):
        return None


def scan_night(night_dir, log=print):
    """Header summary of every frame in the night directory, in time order."""
    frames = []
    for path in night_files(night_dir):
        try:
            with open_fits_file(path) as hdul:
                h = _header(hdul)
                frames.append({
                    'path': path,
                    'kind': frame_kind(h.get('IMAGETYP')),
                    't': _iso(h.get('DATE-OBS')),
                    'exptime': float(h.get('EXPTIME', 0) or 0),
                    'filter': str(h.get('FILTER', '')).replace("'", '').strip(),
                    'telescope': str(h.get('TELESCOP', '')).strip(),
                    'camshutt': smear._parse_camshutt(h['CAMSHUTT']) if 'CAMSHUTT' in h else 'missing',
                    'readoutm': h.get('READOUTM'),
                })
        except Exception as e:
            log('shutter_detect: cannot read %s: %s' % (path, e))
    return sorted((f for f in frames if f['t']), key=lambda f: f['t'])


def night_dir_of(path):
    """(directory, night) of a raw frame filed as .../<YYYYMMDD>/[sub/]frame, else None."""
    d = os.path.dirname(os.path.normpath(path))
    for _ in range(2):
        if readout.valid_night(os.path.basename(d)):
            return d, os.path.basename(d)
        d = os.path.dirname(d)
    return None


def own_night_dir(reddir, runname, night):
    """The night's raw directory, from any T1 list naming one of its own frames."""
    for pattern in ('%s_image_*.list', '%s_flat_*.list', '%s_bias.list', '%s_dark.list'):
        for name in sorted(glob.glob(os.path.join(reddir, pattern % runname))):
            for path in readout._read_list(name):
                found = night_dir_of(path)
                if found and found[1] == night:
                    return found[0]
    return None


def _load(path, telescope, night):
    """Overscan-subtracted trimmed frame, the raw trimmed frame, the frame's
    smear parameters and its raw saturation level."""
    with open_fits_file(path) as hdul:
        params = get_instrument_parameters(hdul)
        full = np.asarray(hdul[0].data, dtype=np.float64)
        overscan = apply_image_slice(hdul, 'overscan', 0.)
        if isinstance(overscan, np.ndarray) and overscan.size > 1:
            overscan = float(np.ma.median(sigma_clip(overscan, maxiters=None)))
        raw = np.asarray(apply_image_slice(hdul, 'trim', hdul[0].data), dtype=np.float64)
        exptime = float(hdul[0].header['EXPTIME'])
    level, noise = readout.overscan_signature(full, params)
    mode = readout.classify(noise, level, telescope, night)
    return {'signal': raw - float(overscan), 'raw': raw, 'params': params, 'smear': params.get('smear') or {},
            'exptime': exptime, 'mode': mode, 'saturation': readout.mode_saturation(telescope, mode, night)}


def _oriented(frame, key):
    s = frame['smear']
    return smear._along_axis(frame[key], s.get('axis', 1), s.get('direction', 1))


def rate_times(smear_params):
    return {k: float(v) for k, v in (smear_params.get('readout_rate_time') or {}).items()}


def identify_rate(t_read, rates):
    """The readout rate whose measured time t_read is within RATE_TOLERANCE of, or None."""
    best = None
    for rate, t in rates.items():
        if t_read is not None and abs(t_read - t) <= RATE_TOLERANCE * t:
            if best is None or abs(t_read - t) < abs(t_read - rates[best]):
                best = rate
    return best


# ------------------------------------------------------------- flats

def measure_flat(path, telescope, night):
    """The readout-ramp test on one flat: t_read, verdict and readout rate."""
    return judge_flat(_load(path, telescope, night))


def judge_flat(f):
    """measure_flat() on a frame from _load()."""
    out = {'exptime': f['exptime'], 'mode': f['mode'], 't_read': None, 'verdict': INCONCLUSIVE, 'rate': None}
    s = f['signal']
    oriented = _oriented(f, 'signal')
    k, edge = oriented.shape[-1], 200
    bands = (float(np.median(oriented[:, edge:2 * edge])), float(np.median(oriented[:, k - 2 * edge:k - edge])))
    out['signal'] = float(np.median(s))
    pedestal = float(np.median(f['raw'] - s))
    if f['exptime'] <= 0:
        out['reason'] = 'no exposure time'
    elif min(bands) < FLAT_MIN_SIGNAL:
        out['reason'] = 'too faint'
    elif max(bands) + pedestal >= FLAT_MAX_FRACTION * f['saturation']:
        out['reason'] = 'near saturation'
    else:
        n = f['smear'].get('n_shifts', oriented.shape[-1])
        t_read = smear_check.ramp_t_read(oriented, n, f['exptime'], edge=edge)
        out['t_read'] = round(float(t_read), 4)
        verdict = smear_check.flat_verdict(t_read, f['exptime'], n=n, span=k - 3 * edge)
        out['verdict'] = {'SHUTTERLESS': SHUTTERLESS, 'shuttered': SHUTTERED}.get(verdict, INCONCLUSIVE)
        if out['verdict'] == SHUTTERLESS:
            out['rate'] = identify_rate(t_read, rate_times(f['smear']))
    return out


def assign_sides(flats, science):
    """Label each flat 'evening' or 'morning'.

    With science frames: before the first is evening, after the last morning,
    and a flat taken during the science goes with the nearer end. Without:
    the largest gap of at least TWILIGHT_GAP_HOURS splits them; one cluster is
    evening if it is earlier than 03:00 UTC after the night's date, else
    morning (SSO twilight flats run ~22-01 and ~09-11 UTC).
    """
    if not flats:
        return
    t = [smear._parse_utc(f['t']) for f in flats]
    if science:
        first, last = smear._parse_utc(science[0]['t']), smear._parse_utc(science[-1]['t'])
        for f, tf in zip(flats, t):
            f['side'] = 'evening' if tf < first or (tf <= last and tf - first < last - tf) else 'morning'
        return
    order = sorted(range(len(t)), key=lambda i: t[i])
    gaps = [(t[order[i + 1]] - t[order[i]], i) for i in range(len(order) - 1)]
    big = max(gaps) if gaps else None
    if big and big[0] >= datetime.timedelta(hours=TWILIGHT_GAP_HOURS):
        cut = t[order[big[1]]]
        for f, tf in zip(flats, t):
            f['side'] = 'evening' if tf <= cut else 'morning'
        return
    night = flats[0].get('night')
    divide = None
    if night:
        divide = datetime.datetime.strptime(night, '%Y%m%d') + datetime.timedelta(days=1, hours=3)
    for f, tf in zip(flats, t):
        f['side'] = 'evening' if divide is None or tf < divide else 'morning'


def ramp_intercept(flats):
    """The readout time behind a twilight's flats: (t0, slope, worst residual).

    A flat's ramp t_read is the readout time plus an apparent term that grows
    with exposure time: the flat field's own large-scale structure (shuttered
    Ganymede flats, both twilights: +0.009 s per s) and the twilight gradient
    (opposite in the evening and the morning). So t_read = t0 + slope * exptime
    over a twilight, and t0 is the readout time: Ganymede 2026-02 gave 4.28 s
    (evening) and 4.44 s (morning) at 1 MHz, 1.48-1.50 s at 3 MHz, ~0.0 s
    shuttered. With one exposure time only, t0 is the median t_read. The line
    is a Theil-Sen fit (median of pairwise slopes), so a flat spoilt by a
    passing cloud moves it little.
    """
    pts = [(f['exptime'], f['t_read']) for f in flats if f.get('t_read') is not None]
    if not pts:
        return None, None, None
    e, t = np.array(pts, dtype=float).T
    i, j = np.triu_indices(len(e), 1)
    distinct = np.abs(e[j] - e[i]) > 0.1
    if not distinct.any():
        t0 = float(np.median(t))
        return round(t0, 3), None, round(float(np.max(np.abs(t - t0))), 3)
    slope = float(np.median((t[j] - t[i])[distinct] / (e[j] - e[i])[distinct]))
    t0 = float(np.median(t - slope * e))
    return round(t0, 3), round(slope, 5), round(float(np.max(np.abs(t - (t0 + slope * e)))), 3)


def side_summary(flats, rates=None):
    """A twilight's state: conclusive only when at least MIN_SIDE_FLATS flats
    are conclusive and every one of them agrees. Its readout rate is the one
    ramp_intercept() names, provided every shutterless flat sits on the line."""
    conclusive = [f for f in flats if f['verdict'] in (SHUTTERED, SHUTTERLESS)]
    states = sorted({f['verdict'] for f in conclusive})
    t_reads = [f['t_read'] for f in conclusive if f['t_read'] is not None]
    t0, slope, resid = ramp_intercept(conclusive)
    out = {'n': len(flats), 'n_conclusive': len(conclusive), 'state': None, 'rate': None,
           't_read_median': round(float(np.median(t_reads)), 3) if t_reads else None,
           't_read_intercept': t0, 'slope': slope, 'max_residual': resid,
           'first': flats[0]['t'] if flats else None, 'last': flats[-1]['t'] if flats else None}
    if len(states) > 1:
        out['state'] = MIXED
    elif states and len(conclusive) >= MIN_SIDE_FLATS:
        out['state'] = states[0]
        if states[0] == SHUTTERLESS:
            if rates is None:
                rates = {r: t for r, t in (('1MHz', 4.33), ('3MHz', 1.5))}
            if resid > RAMP_MIX_FRACTION * abs(t0) + RAMP_MIX_FLOOR:
                out['rate_note'] = 'flats scatter %.2f s about the line: not one readout rate' % resid
            elif slope is not None and abs(slope) > RAMP_MAX_SLOPE:
                out['rate_note'] = 't_read climbs %.3f s per s of exposure: a sky gradient, not a readout' % slope
            else:
                out['rate'] = identify_rate(t0, rates)
                if out['rate'] is None:
                    out['rate_note'] = 'intercept %.2f s names no measured rate' % t0
    return out


# ------------------------------------------------------------- science

def trail_pairs(path, telescope, night):
    """Star-trail (excess, flux) pairs of one raw science frame."""
    f = _load(path, telescope, night)
    oriented = _oriented(f, 'signal')
    stars = smear_check.find_stars(oriented, _oriented(f, 'raw'), f['saturation'])
    _, _, pairs = smear_check.measure(oriented, stars)
    return pairs, f


def hypotheses(rates):
    """Shuttered (t_read 0) and each measured readout rate."""
    out = {SHUTTERED: 0.0}
    out.update(rates)
    return out


def z_scores(t_read, err, sys_err, rates, rate_range=(1.0, 1.0)):
    """Distance in sigma (statistical and systematic together) from each
    hypothesis: 0 for shuttered, rate_range times its time for each rate."""
    sigma = float(np.hypot(err, sys_err))
    out = {}
    for name, h in hypotheses(rates).items():
        lo, hi = (h * rate_range[0], h * rate_range[1]) if name != SHUTTERED else (0.0, 0.0)
        out[name] = max(lo - t_read, t_read - hi, 0.0) / sigma
    return out


def judge(t_read, err, sys_err, rates, rate_range=(1.0, 1.0)):
    """Which hypothesis a measured t_read +/- err (+ systematic) picks out.

    Returns (verdict, rate, margin): verdict SHUTTERED, SHUTTERLESS or
    INCONCLUSIVE; margin is the distance in sigma to the nearest other
    hypothesis, when the best one is consistent.
    """
    if t_read is None or not np.isfinite(t_read) or not np.isfinite(err):
        return INCONCLUSIVE, None, None
    z = z_scores(t_read, err, sys_err, rates, rate_range)
    best = min(z, key=z.get)
    others = [v for k, v in z.items() if k != best]
    margin = round(min(others), 2) if others else None
    if z[best] > CONSISTENT_SIGMA or (margin is not None and margin < EXCLUDED_SIGMA):
        return INCONCLUSIVE, None, margin if z[best] <= CONSISTENT_SIGMA else None
    return (SHUTTERED, None, margin) if best == SHUTTERED else (SHUTTERLESS, best, margin)


def trail_stack(frames, telescope, night, log=print):
    """Stack the trail measurement over frames of one exposure time."""
    pairs, n_frames, smear_params = [], 0, None
    for fr in frames:
        try:
            p, f = trail_pairs(fr['path'], telescope, night)
        except Exception as e:
            log('shutter_detect: trail measurement failed on %s: %s' % (os.path.basename(fr['path']), e))
            continue
        pairs += p
        n_frames += 1
        smear_params = f['smear']
    exptime = frames[0]['exptime']
    out = {'method': 'trail', 'exptime': exptime, 'n_frames': n_frames, 'n_stars': len(pairs),
           'first': frames[0]['t'], 'last': frames[-1]['t'], 't_read': None, 'err': None,
           'verdict': INCONCLUSIVE, 'rate': None, 'margin': None}
    if smear_params is None or len(pairs) < TRAIL_MIN_STARS:
        out['reason'] = 'too few stars'
        return out
    n = smear_params.get('n_shifts', 2048)
    tau, err = smear_check.combine(pairs)
    sys_err = max(TRAIL_SYS_FLOOR, TRAIL_SYS_PER_S * exptime)
    rates, named = candidate_rates(frames, rate_times(smear_params))
    out.update(t_read=round(float(tau * n * exptime), 3), err=round(float(err * n * exptime), 3),
               sys=round(sys_err, 3), hypotheses=sorted(rates))
    if named:
        out['readoutm'] = named
    out['verdict'], out['rate'], out['margin'] = judge(out['t_read'], out['err'], sys_err, rates, TRAIL_RATE_RANGE)
    out['z'] = {k: round(v, 2) for k, v in
                z_scores(out['t_read'], out['err'], sys_err, rates, TRAIL_RATE_RANGE).items()}
    return out


def candidate_rates(frames, rates):
    """The readout rates a measurement must be told apart from: the one READOUTM
    names when every frame carries the card and they agree (Astra states the
    camera's rate), else every rate with a measured time. A trail of 1.5 s
    (3 MHz) sits within the PSF-asymmetry allowance of 0, so without
    READOUTM the trail test can rule out 1 MHz but rarely confirm a shutter."""
    named = {smear.readoutm_rate(f.get('readoutm')) for f in frames}
    if frames and all(f.get('readoutm') is not None for f in frames) and len(named) == 1 and None not in named:
        rate = named.pop()
        return ({rate: rates[rate]} if rate in rates else {}), rate
    return rates, None


def _sample(items, k):
    if len(items) <= k:
        return list(items)
    return [items[int(round(i))] for i in np.linspace(0, len(items) - 1, k)]


def trail_groups(science):
    """Science frames by exposure time, those short and numerous enough for the trail test."""
    groups = {}
    for f in science:
        groups.setdefault(round(f['exptime'], 1), []).append(f)
    return {e: g for e, g in sorted(groups.items()) if 0 < e <= TRAIL_MAX_EXPTIME and len(g) >= TRAIL_MIN_FRAMES}


def _bin_start(n, centre):
    return max(0, min(centre - TRAIL_BIN // 2, n - TRAIL_BIN))


def _conclusive(b):
    return (b['verdict'], b['rate']) if b['verdict'] != INCONCLUSIVE else None


def _side_of(b, rate):
    """'in' when a trail stack reads as shutterless at this readout rate, 'out'
    when it rules the rate out (shuttered, another rate, or simply too little
    trail), None when it can say neither."""
    if b['verdict'] == SHUTTERLESS and b['rate'] == rate:
        return 'in'
    if (b.get('z') or {}).get(rate, 0.0) > EXCLUDED_SIGMA:
        return 'out'
    return None


def _incompatible(a, b):
    """Do two trail stacks rule each other out (one reads as a rate the other excludes)?"""
    return any(x['verdict'] == SHUTTERLESS and _side_of(y, x['rate']) == 'out' for x, y in ((a, b), (b, a)))


def locate_switch(frames, telescope, night, log=print):
    """Bisect time-ordered frames of one exposure time for a change of state.

    Measures a bin of TRAIL_BIN consecutive frames at each end. A switch needs
    one end to read as shutterless at some rate and the other to rule that rate
    out; the bisection then sorts bins by that rate alone, because a shuttered
    bin can rarely be told from 3 MHz but always from the rate seen at the
    other end. A bin that reads like the first end shows that the switch comes
    after the bin's first frame; one like the last, that it comes no later
    than the bin's last frame; so a bin straddling the switch only leaves the
    bracket wider. A straddling bin that reads like neither (a mix of the two
    trails) is resolved by the bins either side of it.

    Returns (bins, switch). switch is None, or {'lo', 'hi', 'rate', 'early',
    'late'}: frames[:lo + 1] are in the early state, frames[hi:] in the late
    one, each a (verdict, rate) or None for an end known only to rule out the
    other end's rate.
    """
    n = len(frames)
    bins = [trail_stack(frames[:TRAIL_BIN], telescope, night, log),
            trail_stack(frames[n - TRAIL_BIN:], telescope, night, log)]
    first, last = bins
    rate = None
    for a, b in ((first, last), (last, first)):
        if a['verdict'] == SHUTTERLESS and _side_of(b, a['rate']) == 'out':
            rate = a['rate']
            break
    if rate is None:
        return bins, None
    early_side = _side_of(first, rate)

    def measure_bin(start):
        b = trail_stack(frames[start:start + TRAIL_BIN], telescope, night, log)
        bins.append(b)
        return _side_of(b, rate)

    lo, hi = 0, n - 1
    while hi - lo > 2 * TRAIL_BIN:
        start = _bin_start(n, (lo + hi) // 2)
        side = measure_bin(start)
        if side is None:
            before, after = max(0, start - TRAIL_BIN), min(n - TRAIL_BIN, start + TRAIL_BIN)
            if measure_bin(before) == early_side and measure_bin(after) not in (None, early_side):
                lo, hi = before, after + TRAIL_BIN - 1
                break
            return bins, None
        if side == early_side:
            lo = start
        else:
            hi = start + TRAIL_BIN - 1
    in_state = (SHUTTERLESS, rate)
    return bins, {'lo': lo, 'hi': hi, 'rate': rate,
                  'early': in_state if early_side == 'in' else _conclusive(first),
                  'late': in_state if early_side == 'out' else _conclusive(last)}


def _compatible(end, state, rate):
    """Is a twilight's state what a located switch says about its end of the night?"""
    if end is None:
        return state != (SHUTTERLESS, rate)
    return end == state


# ------------------------------------------------------------- sky ramp

def reference_flat(night_dir, night, telescope, filt, log=print):
    """Median of normalised shuttered flats from the nearest earlier night (within
    REF_MAX_DAYS) with at least REF_MIN_FLATS of them in this filter: (image,
    night, sides) or (None, None, None). Each flat must itself measure
    shuttered, which is all a reference needs; both twilights are used where
    the night has both, so their opposite twilight gradients partly cancel."""
    root = os.path.dirname(night_dir)
    day = datetime.datetime.strptime(night, '%Y%m%d')
    for back in range(1, REF_MAX_DAYS + 1):
        other = (day - datetime.timedelta(days=back)).strftime('%Y%m%d')
        d = os.path.join(root, other)
        if not os.path.isdir(d):
            continue
        frames = scan_night(d, log=lambda *a: None)
        candidates = [dict(f, night=other) for f in frames if f['kind'] == 'flat' and f['filter'] == filt
                      and 0 < f['exptime'] <= REF_MAX_EXPTIME]
        if len(candidates) < REF_MIN_FLATS:
            continue
        assign_sides(candidates, [f for f in frames if f['kind'] == 'science'])
        stack, sides = [], set()
        for f in _sample(candidates, 2 * REF_MAX_FLATS):
            try:
                loaded = _load(f['path'], telescope, other)
            except Exception:
                continue
            if judge_flat(loaded)['verdict'] == SHUTTERED:
                stack.append(loaded['signal'] / np.median(loaded['signal']))
                sides.add(f['side'])
            if len(stack) >= REF_MAX_FLATS:
                break
        if len(stack) >= REF_MIN_FLATS:
            return np.median(stack, axis=0), other, sorted(sides)
    return None, None, None


def sky_ramp(frames, ref, telescope, night, log=print):
    """EXPERIMENTAL. t_read from the sky ramp along the readout axis of science
    frames divided by a shuttered reference flat; the same ramp across the
    other axis is recorded as a control for real sky gradients."""
    along, across, smear_params = [], [], None
    for fr in _sample(frames, SKY_SAMPLE):
        try:
            f = _load(fr['path'], telescope, night)
        except Exception as e:
            log('shutter_detect: sky ramp failed on %s: %s' % (os.path.basename(fr['path']), e))
            continue
        if np.shape(ref) != np.shape(f['signal']):
            continue
        smear_params = f['smear']
        n = smear_params.get('n_shifts', 2048)
        ratio = f['signal'] / ref
        oriented = _oriented(dict(f, ratio=ratio), 'ratio')
        along.append(smear_check.ramp_t_read(oriented, n, f['exptime'], edge=min(200, oriented.shape[-1] // 5)))
        other = np.moveaxis(ratio, 1 - smear_params.get('axis', 1), -1)
        across.append(smear_check.ramp_t_read(other, n, f['exptime'], edge=min(200, other.shape[-1] // 5)))
    out = {'method': 'sky_ramp', 'exptime': frames[0]['exptime'], 'n_frames': len(along),
           'first': frames[0]['t'], 'last': frames[-1]['t'], 'verdict': INCONCLUSIVE, 'rate': None,
           'experimental': True}
    if len(along) < SKY_MIN_FRAMES:
        out['reason'] = 'too few frames'
        return out
    along, across = np.array(along), np.array(across)
    t_read = float(np.median(along))
    err = 1.253 * float(np.std(along)) / np.sqrt(len(along))
    # the cross-axis ramp says how large a sky gradient can be
    sys_err = max(SKY_SYS_FLOOR, float(np.median(np.abs(across))))
    out.update(t_read=round(t_read, 3), err=round(err, 3), sys=round(sys_err, 3),
               across_median=round(float(np.median(across)), 3))
    out['verdict'], out['rate'], out['margin'] = judge(t_read, err, sys_err,
                                                       candidate_rates(frames, rate_times(smear_params))[0])
    return out


# ------------------------------------------------------------- the night

def _segment(start, end, state, rate, t_read, source, confidence, note=''):
    return {'start': start, 'end': end, 'state': state, 'rate': rate, 't_read': t_read,
            'source': source, 'confidence': confidence, 'note': note}


def _t_read(rate, rates, measured):
    """A shutterless segment's readout time: its rate's measured time, else the
    ramp's own value (and then it is at best probable)."""
    return rates.get(rate) if rate else measured


def _readoutm_rates(frames):
    rates = {}
    for f in frames:
        if f['readoutm'] is not None:
            rate = smear.readoutm_rate(f['readoutm'])
            rates[str(rate)] = rates.get(str(rate), 0) + 1
    return rates


def _camshutt_counts(frames):
    counts = {}
    for f in frames:
        key = {True: 'T', False: 'F', None: 'unreadable'}.get(f['camshutt'], 'missing')
        counts[key] = counts.get(key, 0) + 1
    return counts


def detect_night(night_dir, night, telescope, science=True, all_levels=False, frames=None, overrides=None,
                 log=print):
    """Measure a night's shutter state. Reads frames only; writes nothing.

    science=False stops at the flats (a night borrowed for its flats).
    all_levels=True runs every test whatever the flats say (validation).
    """
    telescope = readout.canonical_telescope(telescope)
    if frames is None:
        frames = scan_night(night_dir, log=log)
    rec = {'version': RECORD_VERSION, 'telescope': telescope, 'night': night, 'night_dir': night_dir,
           'measured_at': datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None, microsecond=0)
           .isoformat(),
           'counts': {}, 'flats': [], 'sides': {}, 'science': [], 'segments': [], 'disagreements': [],
           'notes': []}
    for f in frames:
        rec['counts'][f['kind']] = rec['counts'].get(f['kind'], 0) + 1
    if not frames:
        rec['notes'].append('no frames')
        rec['result'] = {'state': None, 'source': None, 'confidence': INCONCLUSIVE, 'rate': None, 't_read': None}
        return rec
    rec['span'] = [frames[0]['t'], frames[-1]['t']]
    sci = [f for f in frames if f['kind'] == 'science']
    flats = [dict(f, night=night) for f in frames if f['kind'] == 'flat']

    rates = {}
    for f in flats:
        try:
            m = measure_flat(f['path'], telescope, night)
        except Exception as e:
            log('shutter_detect: flat measurement failed on %s: %s' % (os.path.basename(f['path']), e))
            m = {'t_read': None, 'verdict': INCONCLUSIVE, 'rate': None, 'reason': str(e)}
        f.update(m)
        if not rates:
            try:
                with open_fits_file(f['path']) as hdul:
                    rates = rate_times(get_instrument_parameters(hdul).get('smear') or {})
            except Exception:
                pass
    assign_sides(flats, sci)
    for f in flats:
        f.pop('night', None)
        f['file'] = os.path.basename(f['path'])
    rec['flats'] = flats
    if not rates and sci:
        try:
            with open_fits_file(sci[0]['path']) as hdul:
                rates = rate_times(get_instrument_parameters(hdul).get('smear') or {})
        except Exception:
            pass
    for side in ('evening', 'morning'):
        rec['sides'][side] = side_summary([f for f in flats if f.get('side') == side], rates)

    ev, mo = rec['sides']['evening'], rec['sides']['morning']
    good = (SHUTTERED, SHUTTERLESS)
    ev_ok, mo_ok = ev['state'] in good, mo['state'] in good
    first, last = frames[0]['t'], frames[-1]['t']
    segments = []
    result = None
    sides_agree = ev_ok and mo_ok and ev['state'] == mo['state'] and ev['rate'] == mo['rate']

    # the science, when the flats do not settle the night (or always, for validation)
    groups = trail_groups(sci) if science else {}
    need_science = science and (all_levels or not sides_agree)
    night_sci = []
    if need_science:
        for exptime, g in groups.items():
            night_sci.append(trail_stack(_sample(g, TRAIL_SAMPLE), telescope, night, log))
    # The sky ramp is experimental and never decides anything on its own, and
    # finding a reference flat reads up to REF_MAX_DAYS earlier nights, so
    # outside validation it runs only when nothing else gave an answer.
    trail_answered = any(st['verdict'] != INCONCLUSIVE for st in night_sci)
    if need_science and (all_levels or not (ev_ok or mo_ok or trail_answered)):
        long_groups = {}
        for f in sci:
            if f['exptime'] > TRAIL_MAX_EXPTIME:
                long_groups.setdefault((round(f['exptime'], 1), f['filter']), []).append(f)
        refs = {}
        for (exptime, filt), g in sorted(long_groups.items()):
            if len(g) < SKY_MIN_FRAMES:
                continue
            if filt not in refs:
                refs[filt] = reference_flat(night_dir, night, telescope, filt, log)
            ref, ref_night, ref_sides = refs[filt]
            if ref is None:
                night_sci.append({'method': 'sky_ramp', 'exptime': exptime, 'n_frames': 0, 'verdict': INCONCLUSIVE,
                                  'reason': 'no shuttered reference night within %d days' % REF_MAX_DAYS,
                                  'experimental': True})
                continue
            s = sky_ramp(g, ref, telescope, night, log)
            s['reference_night'], s['reference_sides'] = ref_night, ref_sides
            night_sci.append(s)
    rec['science'] = night_sci
    trail = [st for st in night_sci if st['method'] == 'trail' and st['verdict'] != INCONCLUSIVE]
    trail_states = {(st['verdict'], st['rate']) for st in trail}
    if len(trail_states) > 1:
        rec['notes'].append('exposure groups disagree: %s' % sorted(map(str, trail_states)))
    sci_state = next(iter(trail_states)) if len(trail_states) == 1 else None
    sci_margin = min((st['margin'] for st in trail if st['margin'] is not None), default=None)
    sky = [st for st in night_sci if st['method'] == 'sky_ramp' and st['verdict'] != INCONCLUSIVE]
    sky_states = {(st['verdict'], st['rate']) for st in sky}
    sky_state = next(iter(sky_states)) if len(sky_states) == 1 else None

    # A whole-night stack would average a switch away (a third of the night
    # at 1 MHz reads as 1.5 s, i.e. 3 MHz), so bins at both ends of the
    # largest exposure group are measured too; ends that rule each other out
    # are a switch, which locate_switch() then brackets.
    ends, switch, g = [], None, []
    if need_science and groups:
        g = max(groups.values(), key=len)
        if len(g) >= 2 * TRAIL_BIN:
            ends, switch = locate_switch(g, telescope, night, log)
            rec['science'] += [dict(b, bin=True) for b in ends]
    ends_differ = len(ends) >= 2 and _incompatible(ends[0], ends[1])
    end_states = {_conclusive(b) for b in ends[:2]} - {None}
    if sci_state is not None and (ends_differ or end_states - {sci_state}):
        rec['notes'].append('the night\'s end bins (%s) do not match the whole-night stack %s'
                            % (sorted(map(str, end_states)), sci_state))
        sci_state, sci_margin = None, None

    def whole(state, rate, measured, source, confidence, note=''):
        t_read = _t_read(rate, rates, measured) if state == SHUTTERLESS else None
        segments.append(_segment(first, last, state, rate, t_read, source, confidence, note))
        return {'state': state, 'source': source, 'confidence': confidence, 'rate': rate, 't_read': t_read}

    def measured_flats(side):
        return side['t_read_intercept'] if side['state'] == SHUTTERLESS else None

    def located_split(before, after, conf_before, conf_after, source, notes):
        """Segments either side of the located switch: before up to frames[lo],
        after from frames[hi], the open interval last (segment_for takes the
        first segment holding a time). A state of None is an end known only to
        rule out the other end's rate."""
        lo, hi = g[switch['lo']]['t'], g[switch['hi']]['t']
        for st, conf, start_t, end_t, note in ((before, conf_before, first, lo, notes[0]),
                                               (after, conf_after, hi, last, notes[1])):
            state, rate, measured = st if st is not None else (None, None, None)
            if state is None:
                conf, note = INCONCLUSIVE, note + '; rules out %s, nothing more' % switch['rate']
            segments.append(_segment(start_t, end_t, state, rate,
                                     _t_read(rate, rates, measured) if state == SHUTTERLESS else None,
                                     source, conf, note))
        segments.append(_segment(lo, hi, None, None, None, 'science', INCONCLUSIVE, 'the switch lies in here'))
        conf = UNAMBIGUOUS if conf_before == conf_after == UNAMBIGUOUS else PROBABLE
        return {'state': 'split', 'source': source, 'confidence': conf, 'rate': None, 't_read': None,
                'switch': [lo, hi]}

    def flat_segment(side, name, note):
        segments.append(_segment(side['first'], side['last'], side['state'], side['rate'],
                                 _t_read(side['rate'], rates, measured_flats(side)), 'flats', UNAMBIGUOUS,
                                 '%s flats themselves; %s' % (name, note)))

    def end_state(end):
        return None if end is None else (end[0], end[1], None)

    if sides_agree:
        result = whole(ev['state'], ev['rate'], measured_flats(ev), 'flats', UNAMBIGUOUS,
                       'evening and morning flats agree')
        if ev['state'] == SHUTTERLESS and ev['rate'] is None:
            result['confidence'] = segments[-1]['confidence'] = PROBABLE
            segments[-1]['note'] = 'flats agree, but their t_read names no measured readout rate'
    elif ev_ok and mo_ok:
        # a switch during the night: the flats fix each end, the science places it
        if switch is not None and _compatible(switch['early'], (ev['state'], ev['rate']), switch['rate']) \
                and _compatible(switch['late'], (mo['state'], mo['rate']), switch['rate']):
            result = located_split((ev['state'], ev['rate'], measured_flats(ev)),
                                   (mo['state'], mo['rate'], measured_flats(mo)), UNAMBIGUOUS, UNAMBIGUOUS,
                                   'flats+science', ('evening flats; science bins up to the switch',
                                                     'science bins from the switch; morning flats'))
        else:
            result = {'state': 'split', 'source': 'flats', 'confidence': PROBABLE, 'rate': None, 't_read': None}
            flat_segment(ev, 'evening', 'the switch is not located')
            flat_segment(mo, 'morning', 'the switch is not located')
            if switch is not None:
                rec['notes'].append('the switch the science locates does not match the flats')
            rec['notes'].append('evening and morning flats disagree and the science does not locate the switch')
    elif ev_ok or mo_ok:
        side = ev if ev_ok else mo
        name = 'evening' if ev_ok else 'morning'
        other = mo if ev_ok else ev
        if other['state'] == MIXED:
            rec['notes'].append('%s flats disagree among themselves' % ('morning' if ev_ok else 'evening'))
        if ends_differ:
            near = (switch['early'] if ev_ok else switch['late']) if switch else None
            if switch is not None and _compatible(near, (side['state'], side['rate']), switch['rate']):
                # the flats and the bin at their end agree; the other end is the science's alone
                flats_end = (side['state'], side['rate'], measured_flats(side))
                far = end_state(switch['late'] if ev_ok else switch['early'])
                before, after = (flats_end, far) if ev_ok else (far, flats_end)
                confs = (UNAMBIGUOUS, PROBABLE) if ev_ok else (PROBABLE, UNAMBIGUOUS)
                result = located_split(before, after, confs[0], confs[1], 'flats+science',
                                       ('%s flats and science bins' % name if ev_ok else 'science bins only',
                                        'science bins only' if ev_ok else 'science bins and %s flats' % name))
            else:
                result = {'state': None, 'source': 'flats+science', 'confidence': INCONCLUSIVE, 'rate': None,
                          't_read': None}
                flat_segment(side, name, 'the science changes during the night')
                rec['notes'].append('%s flats only, and the science changes during the night' % name)
        else:
            # A shutterless twilight needs a conclusive trail at the same rate. A
            # shuttered one needs every stack, the end bins included, to be
            # shuttered or at least to rule out the slowest readout (1 MHz, the
            # most visible trail): the PSF allowance rarely lets the trail test
            # exclude 3 MHz, but a short shuttered flat already has. Merely
            # "consistent with no trail" is not enough: a night half at 1 MHz
            # can be that.
            stacks = [st for st in night_sci if st['method'] == 'trail' and st.get('z')] + ends
            slowest = max(rates, key=rates.get) if rates else None
            if side['state'] == SHUTTERLESS:
                supports = sci_state == (side['state'], side['rate'])
                contradicts = sci_state is not None and not supports
            else:
                contradicts = any(st['verdict'] == SHUTTERLESS for st in stacks)
                supports = (not contradicts and bool(stacks) and slowest is not None
                            and all(st['verdict'] == SHUTTERED or _side_of(st, slowest) == 'out' for st in stacks))
            if supports:
                result = whole(side['state'], side['rate'], measured_flats(side), 'flats+science', UNAMBIGUOUS,
                               '%s flats, confirmed by the science' % name)
            elif contradicts:
                result = {'state': None, 'source': 'flats+science', 'confidence': INCONCLUSIVE, 'rate': None,
                          't_read': None}
                flat_segment(side, name, 'the science disagrees')
                rec['notes'].append('%s flats say %s, science says %s' % (
                    name, (side['state'], side['rate']), sorted(str((st['verdict'], st['rate'])) for st in trail)))
            else:
                result = whole(side['state'], side['rate'], measured_flats(side), 'flats', PROBABLE,
                               '%s flats only; science %s' % (name, 'inconclusive' if need_science and night_sci
                                                             else 'not usable'))
    elif ends_differ and switch is not None:
        result = located_split(end_state(switch['early']), end_state(switch['late']), PROBABLE, PROBABLE, 'science',
                               ('no usable flats; science bins', 'no usable flats; science bins'))
    elif sci_state is not None:
        confident = sci_margin is not None and sci_margin >= SCIENCE_ALONE_SIGMA
        result = whole(sci_state[0], sci_state[1], None, 'science', UNAMBIGUOUS if confident else PROBABLE,
                       'no usable flats; star trails alone, margin %s sigma' % sci_margin)
    elif sky_state is not None:
        result = whole(sky_state[0], sky_state[1], None, 'sky_ramp', PROBABLE,
                       'no usable flats; sky ramp (experimental)')
    else:
        result = {'state': None, 'source': None, 'confidence': INCONCLUSIVE, 'rate': None, 't_read': None}
        if ev['state'] == MIXED or mo['state'] == MIXED:
            rec['notes'].append('flats disagree among themselves')

    # Astra's READOUTM is the camera's own statement of the readout rate
    readoutm = _readoutm_rates(frames)
    rec['readoutm'] = readoutm
    for seg in segments:
        if seg['state'] == SHUTTERLESS and readoutm:
            in_seg = [f for f in frames if seg['start'] <= f['t'] <= seg['end'] and f['readoutm'] is not None]
            named = {str(smear.readoutm_rate(f['readoutm'])) for f in in_seg}
            if named and named != {str(seg['rate'])}:
                rec['disagreements'].append({'kind': 'READOUTM', 'start': seg['start'], 'end': seg['end'],
                                             'measured': seg['rate'], 'header': sorted(named)})
                if seg['confidence'] == UNAMBIGUOUS:
                    seg['confidence'] = PROBABLE
                    seg['note'] += '; READOUTM names %s' % sorted(named)
    _compare_with_headers(rec, frames, segments, overrides)
    rec['segments'] = segments
    if any(s['confidence'] != UNAMBIGUOUS for s in segments if s['state']) and result['confidence'] == UNAMBIGUOUS:
        result['confidence'] = PROBABLE
    rec['result'] = result
    rec['camshutt'] = _camshutt_counts(frames)
    return rec


def _compare_with_headers(rec, frames, segments, overrides):
    """Disagreements between the measured segments and CAMSHUTT, and the override table."""
    if overrides is None:
        overrides = smear._shipped_overrides()
    for seg in segments:
        if seg['state'] not in (SHUTTERED, SHUTTERLESS) or seg['confidence'] == INCONCLUSIVE:
            continue
        measured_used = seg['state'] == SHUTTERED
        in_seg = [f for f in frames if seg['start'] <= f['t'] <= seg['end']]
        wrong = [f for f in in_seg if f['camshutt'] in (True, False) and f['camshutt'] != measured_used]
        if wrong:
            rec['disagreements'].append({
                'kind': 'CAMSHUTT', 'start': seg['start'], 'end': seg['end'], 'measured': seg['state'],
                'confidence': seg['confidence'], 'header': 'T' if wrong[0]['camshutt'] else 'F',
                'n_frames': len(wrong), 'of': len(in_seg)})
        hits = {}
        for f in in_seg:
            when = smear._parse_utc(f['t'])
            for row in overrides:
                if row['telescope'] == rec['telescope'] and row['start'] <= when < row['end']:
                    hits.setdefault(id(row), [row, 0])[1] += 1
        if seg['state'] == SHUTTERLESS and not hits and not any(f['camshutt'] in (True, False) for f in in_seg):
            rec['disagreements'].append({'kind': 'no-header', 'start': seg['start'], 'end': seg['end'],
                                         'measured': seg['state'], 'measured_t_read': seg['t_read'],
                                         'confidence': seg['confidence'], 'n_frames': len(in_seg)})
        for row, n in hits.values():
            t_row = row.get('readout_time')
            state_differs = row['camshutt'] != measured_used
            rate_differs = (not state_differs and seg['state'] == SHUTTERLESS and t_row is not None
                            and seg['t_read'] is not None and abs(t_row - seg['t_read']) > RATE_TOLERANCE * t_row)
            if state_differs or rate_differs:
                rec['disagreements'].append({
                    'kind': 'override', 'start': seg['start'], 'end': seg['end'], 'measured': seg['state'],
                    'measured_t_read': seg['t_read'], 'confidence': seg['confidence'],
                    'override': 'T' if row['camshutt'] else 'F', 'override_t_read': t_row,
                    'window': [row['start'].isoformat(), row['end'].isoformat()], 'n_frames': n})


# ------------------------------------------------------------- frames and records

def segment_for(record, date_obs):
    """The record's segment holding DATE-OBS, or None (a frame from another night)."""
    if not record or not date_obs:
        return None
    t = _iso(date_obs)
    if t is None:
        return None
    for seg in record.get('segments', []):
        if seg['start'] <= t <= seg['end']:
            return seg
    return None


def header_cards(record, header):
    """SHUTD* cards: what the measurement says about this frame, and whether CAMSHUTT agrees."""
    seg = segment_for(record, header.get('DATE-OBS'))
    if record is None:
        return []
    if seg is None or seg['state'] is None:
        return [('SHUTDET', 'UNKNOWN', 'Measured shutter state (shutter_detect)'),
                ('SHUTDCNF', INCONCLUSIVE, 'Confidence of the measured shutter state')]
    cards = [('SHUTDET', seg['state'].upper(), 'Measured shutter state (shutter_detect)'),
             ('SHUTDSRC', seg['source'], 'Evidence for SHUTDET'),
             ('SHUTDCNF', seg['confidence'], 'Confidence of the measured shutter state')]
    if seg['state'] == SHUTTERLESS and seg.get('t_read') is not None:
        cards.append(('SHUTDTRD', float(seg['t_read']),
                      '[s] Measured readout time (%s)' % (seg['rate'] or 'no known rate')))
    if 'CAMSHUTT' in header:
        used = smear._parse_camshutt(header['CAMSHUTT'])
        if used is not None:
            cards.append(('SHUTDIS', used != (seg['state'] == SHUTTERED), 'Measured state disagrees with CAMSHUTT'))
    return cards


def record_path(reddir, night=None):
    """shutter.json for the reduction's own night, shutter_<night>.json for a borrowed one."""
    return os.path.join(reddir, RECORD_NAME if night is None else 'shutter_%s.json' % night)


def write_record(record, path):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(record, f, indent=1, sort_keys=True, default=readout._json_default)
    os.replace(tmp, path)
    return path


def read_record(path):
    try:
        with open(path) as f:
            record = json.load(f)
    except (OSError, ValueError):
        return None
    return record if record.get('version') == RECORD_VERSION else None


def summary(record):
    r = record.get('result') or {}
    lines = ['%s %s: %s (%s, %s)%s' % (record.get('telescope'), record.get('night'), r.get('state') or 'unknown',
                                        r.get('source') or 'no evidence', r.get('confidence'),
                                        ', t_read %.2f s (%s)' % (r['t_read'], r.get('rate'))
                                        if r.get('t_read') else '')]
    for side, s in sorted((record.get('sides') or {}).items()):
        if s.get('n'):
            lines.append('  %-8s flats %2d, conclusive %2d: %s%s, ramp t_read %s s at 0 s exposure (median %s)'
                         % (side, s['n'], s['n_conclusive'], s['state'] or 'inconclusive',
                            ' ' + s['rate'] if s.get('rate') else '', s.get('t_read_intercept'),
                            s['t_read_median']))
    for s in record.get('science') or []:
        if s.get('reason') and s.get('t_read') is None:
            lines.append('  %-8s %5.1f s: %s' % (s['method'], s.get('exptime') or 0, s['reason']))
            continue
        lines.append('  %-8s %s%5.1f s x %3d frames: t_read %s +/- %s s -> %s%s'
                     % (s['method'], 'bin ' if s.get('bin') else '', s.get('exptime') or 0, s.get('n_frames', 0),
                        s.get('t_read'), s.get('err'), s['verdict'], ' ' + s['rate'] if s.get('rate') else ''))
    segments = record.get('segments') or []
    if len(segments) > 1:
        for seg in sorted(segments, key=lambda x: x['start']):
            lines.append('  segment  %s to %s: %s%s (%s, %s)' % (
                seg['start'], seg['end'], seg['state'] or 'unknown', ' ' + seg['rate'] if seg.get('rate') else '',
                seg['source'], seg['confidence']))
    for d in record.get('disagreements') or []:
        lines.append('  DISAGREES with %s: %s' % (d['kind'], {k: v for k, v in d.items() if k != 'kind'}))
    for n in record.get('notes') or []:
        lines.append('  note: ' + n)
    return '\n'.join(lines)


def warn_disagreements(record, log=print):
    for d in record.get('disagreements') or []:
        if d['kind'] == 'CAMSHUTT':
            log('WARNING: %s %s: frames measure %s (%s) but CAMSHUTT = %s on %d of %d frames'
                % (record['telescope'], record['night'], d['measured'], d['confidence'], d['header'],
                   d['n_frames'], d['of']))
        elif d['kind'] == 'override':
            log('WARNING: %s %s: frames measure %s (t_read %s, %s) but the override window %s says CAMSHUTT %s'
                '%s' % (record['telescope'], record['night'], d['measured'], d['measured_t_read'],
                        d['confidence'], ' to '.join(d['window']), d['override'],
                        ', t_read %s' % d['override_t_read'] if d['override_t_read'] else ''))
        elif d['kind'] == 'no-header':
            log('WARNING: %s %s: frames measure shutterless (t_read %s, %s) with no CAMSHUTT and no override '
                'window, so they are not de-smeared unless the measurement is applied'
                % (record['telescope'], record['night'], d['measured_t_read'], d['confidence']))
        else:
            log('WARNING: %s %s: shutterless frames measure %s but READOUTM names %s'
                % (record['telescope'], record['night'], d['measured'], d['header']))


def load_or_detect(reddir, night_dir, night, telescope, own=True, refresh=False, write=True, log=print):
    """The night's record: read back if this code wrote one for the night, else
    measured now and saved. own=False is a night borrowed for its flats: flats
    only, saved as shutter_<night>.json. T2 passes refresh=True.

    Never raises: a failed measurement is a record with no segments.
    """
    path = record_path(reddir, None if own else night)
    telescope = readout.canonical_telescope(telescope)
    if not refresh:
        record = read_record(path)
        if record is not None and record.get('night') == night and record.get('telescope') == telescope:
            return record
    try:
        if not night_dir or not os.path.isdir(night_dir):
            raise OSError('no raw directory for night %s (%s)' % (night, night_dir))
        record = detect_night(night_dir, night, telescope, science=own, log=log)
    except Exception as e:
        log('shutter_detect: measurement failed for %s %s: %s' % (telescope, night, e))
        record = {'version': RECORD_VERSION, 'telescope': telescope, 'night': night, 'night_dir': night_dir,
                  'segments': [], 'disagreements': [], 'notes': ['measurement failed: %s' % e],
                  'result': {'state': None, 'source': None, 'confidence': INCONCLUSIVE, 'rate': None,
                             't_read': None}}
    log(summary(record))
    warn_disagreements(record, log)
    if write:
        try:
            write_record(record, path)
        except OSError as e:
            log('shutter_detect: could not write %s: %s' % (path, e))
    return record


class NightRecords:
    """The shutter record of whichever night a frame belongs to: the
    reduction's own night (from T2's shutter.json), or a night borrowed for
    its flats (flats only; measured afresh with refresh_borrowed, as T4 does)."""

    def __init__(self, reddir, runname, night=None, refresh_borrowed=False, log=print):
        self.reddir, self.runname, self.log = reddir, runname, log
        self.night = str(night or readout.night_from_reddir(reddir))
        self.refresh_borrowed = refresh_borrowed
        self._records = {}

    def own(self, telescope):
        return self.get(own_night_dir(self.reddir, self.runname, self.night), self.night, telescope, own=True)

    def get(self, night_dir, night, telescope, own):
        if night not in self._records:
            self._records[night] = load_or_detect(self.reddir, night_dir, night, telescope, own=own,
                                                  refresh=self.refresh_borrowed and not own, log=self.log)
        return self._records[night]

    def for_frame(self, path, telescope):
        found = night_dir_of(path)
        if found is None:
            return None
        night_dir, night = found
        return self.get(night_dir, night, telescope, own=(night == self.night))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('night_dir', help='Observations/<TEL>/images/<NIGHT>')
    parser.add_argument('--telescope', help='TELESCOP (default: from the first frame)')
    parser.add_argument('--flats-only', action='store_true', help='stop at the flats')
    parser.add_argument('--all-levels', action='store_true', help='run every test whatever the flats say')
    parser.add_argument('--json', action='store_true', help='print the whole record')
    parser.add_argument('--out', help='also write the record to this file')
    args = parser.parse_args(argv)
    night_dir = os.path.normpath(args.night_dir)
    night = os.path.basename(night_dir)
    frames = scan_night(night_dir)
    if not frames:
        print('no frames in %s' % night_dir)
        return 1
    telescope = args.telescope
    if not telescope:
        with open_fits_file(frames[0]['path']) as hdul:
            telescope = frame_telescope(hdul)
    record = detect_night(night_dir, night, telescope, science=not args.flats_only, all_levels=args.all_levels,
                          frames=frames)
    if args.out:
        write_record(record, args.out)
    print(json.dumps(record, indent=1, default=readout._json_default) if args.json else summary(record))
    return 0


if __name__ == '__main__':
    sys.exit(main())
