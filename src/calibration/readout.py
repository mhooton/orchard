"""
Per-night gain, read noise and saturation level for the Andor cameras.

No SPECULOOS or Artemis raw header records the readout mode, the gain or the
read noise, and the cameras have not always been read out in the same mode.
instrument_config.json carries one gain per telescope and a 62,000 ADU
saturation level, which are wrong for every night in a non-standard mode
(REGIMES below).

Each frame's values are chosen in three steps:

1. Its readout mode is recognised from its own overscan strip.  The noise
   (ADU) and the level of that strip differ between modes (MODES below), and
   the mode can change during a night: Europa switched to the 3.8 e-/ADU mode
   halfway through 2019-07-25, and Artemis left it during 2026-08-07.
2. The read noise of each mode seen that night is measured from the night's
   own bias pairs, and its gain from the night's own twilight-flat pairs, as
   in mh_scratch/header_audit/src/ptc.py:
       RON  = sigma(B1 - B2) / sqrt(2)
       gain = (m1 + k^2 m2) / (var(F1 - k F2) - (1 + k^2) RON^2),  k = m1 / m2
   with the variances taken as the median of sigma-clipped 64x64 tiles of the
   central 512x512 pixels, so a residual illumination gradient between the
   two flats does not count as noise.  Frames borrowed from earlier nights by
   createlists.py are not used: they may be in another mode.
3. A mode the night's own frames could not measure takes the REGIMES value for
   the telescope and date, and failing that the instrument config.  GAINSRC
   says which: night, table, estimate (a table gain that was never measured)
   or config.

T2 (pipebias.py) writes the night's measurements to readout.json in the
reduction directory; T6 (pipered.py) writes each frame's GAIN, RDNOISE, RON,
SATURATE, READMODE, GAINSRC and OSNOISE into the calibrated frame.

Instruments without a configured overscan (SPIRIT, MOANA) keep the config
gain and saturation and the pipeline's original read-noise estimate; their
READMODE is DEFAULT.

    python -m calibration.readout <reduction dir> [--night YYYYMMDD]
measures a night and prints the record without writing anything.
"""
import argparse
import functools
import glob
import json
import math
import os
import sys
import warnings
from collections import namedtuple

import numpy as np
from astropy.stats import sigma_clipped_stats

RECORD_NAME = 'readout.json'
RECORD_VERSION = 1

# Readout modes of the Andor iKon-L cameras, measured 2026-10-16/17 from bias
# pairs, flat pairs and the overscan strip of every kind of frame
# (mh_scratch/header_audit and mh_scratch/readout_modes on appcs).
#   noise      typical overscan noise in ADU (bias-pair noise is within 5%)
#   pedestal   range of the overscan level in ADU; bright flats pull it down
#              by up to ~60 ADU (~230 ADU in PED7500)
#   tolerance  accepted ratio between a frame's overscan noise and `noise`
#              when the dated table expects another mode
#   anywhere   False when the signature is not unique to the mode, so the mode
#              is only recognised on a telescope REGIMES says has used it
#   gain, saturation  for a telescope with no REGIMES entry for the mode
Mode = namedtuple('Mode', 'noise pedestal tolerance anywhere gain saturation')
MODES = {
    # 5.7-6.9 ADU bias noise, 0.97-1.13 e-/ADU by telescope and date.  The ADC
    # tops out at 62,600-63,950 ADU depending on the camera, hence the config's
    # 62,000.
    'NORMAL': Mode(6.2, (150.0, 1000.0), 1.30, True, 1.05, 62000.0),
    # 2.4-2.6 ADU, 3.6-3.9 e-/ADU.  The pixels fill (~140,000 e-) well below
    # the ADC limit: Europa's response bends above ~30,000 ADU, Artemis's
    # flats stop at 39,300-40,800.
    'GAIN3.8': Mode(2.45, (150.0, 1000.0), 1.30, True, 3.80, 34000.0),
    # Ganymede 2026-02-27..03-01: 9.0 ADU bias (overscan 9.4-11), 1.19 e-/ADU,
    # bias level 319, 1.5 s readout (the shutterless trails); stars pile up at
    # 62,885 ADU.  The noise alone is not unique: Artemis read out 6.4-10.7
    # ADU at 1.03 e-/ADU in 2023-07..2024-01.
    'FAST': Mode(9.8, (150.0, 1000.0), 1.30, False, 1.19, 59500.0),
    # Ganymede 2026-03-02..03: 17.6 ADU bias (overscan 17.6-18.6, up to 30 in
    # bright flats), 1.34 e-/ADU, bias level 7,546.  The response bends above
    # ~39,500 ADU and the centre of a flat stops rising near 43,700.
    'PED7500': Mode(18.2, (6000.0, 9000.0), 2.0, True, 1.34, 39500.0),
    # Artemis 2026-08-06..07: 6.3 ADU bias, bias level ~2,900 (flats ~2,810);
    # twilight flats stop rising at a median of ~24,100 ADU, so no flat pair
    # was below saturation.  The gain is estimated, not measured: the same
    # stars give 5.9 and 7.1 e-/ADU against 2026-08-05 and 08-08 (3.83), and
    # a ~150,000 e- full well gives 7.0.
    'PED2900': Mode(6.3, (2000.0, 4000.0), 2.0, True, None, 22800.0),
}
# The mode the dated table expects is kept up to this ratio; wider for NORMAL,
# whose noise wanders (Artemis 2023-07..2024-01 read 6.4-10.7 ADU within a
# night at 1.03 e-/ADU), so such a night keeps one label.
EXPECTED_TOLERANCE = {'NORMAL': 1.8}
EXPECTED_TOLERANCE_DEFAULT = 1.45

# Dated readout regimes per telescope (TELESCOP keyword).  Nights are the
# pipeline's YYYYMMDD night labels, both ends inclusive; an entry that is not
# NORMAL overrides the NORMAL entry for the same night.  `gain` and
# `saturation` (raw ADU) serve a mode the night's own frames could not
# measure; `measured` is False for an estimated gain.
Regime = namedtuple('Regime', 'first last mode gain saturation measured note')
REGIMES = {
    'SPECULOOS-IO': [
        Regime('20100101', '20991231', 'NORMAL', 1.04, 62000.0, True,
               'flat pairs 2018-2026: 1.02-1.05'),
    ],
    'SPECULOOS-EUROPA': [
        Regime('20100101', '20991231', 'NORMAL', 1.08, 62000.0, True,
               'flat pairs 2017-2026: 1.06-1.10'),
        Regime('20190725', '20191026', 'GAIN3.8', 3.80, 34000.0, True,
               'switched during 2019-07-25; 3.77-3.84; response bends above ~30,000'),
    ],
    'SPECULOOS-GANYMEDE': [
        Regime('20100101', '20220109', 'NORMAL', 0.99, 62000.0, True,
               'flat pairs 2019-03 to 2022-01-09: 0.97-1.00'),
        Regime('20220110', '20991231', 'NORMAL', 1.11, 62000.0, True,
               'flat pairs from 2022-01-11: 1.09-1.12'),
        Regime('20260227', '20260301', 'FAST', 1.19, 59500.0, True,
               '1.19 on all three nights'),
        Regime('20260302', '20260303', 'PED7500', 1.34, 39500.0, True,
               '1.32-1.37 below 30,000 ADU above bias'),
    ],
    'SPECULOOS-CALLISTO': [
        Regime('20100101', '20991231', 'NORMAL', 1.01, 62000.0, True,
               'flat pairs 2018, 2021: 1.01'),
    ],
    'Artemis': [
        Regime('20100101', '20991231', 'NORMAL', 1.02, 62000.0, True,
               'flat pairs 2019-2025: 0.97-1.04, mostly 1.01-1.03'),
        Regime('20200806', '20200807', 'GAIN3.8', 3.85, 37300.0, True,
               'flat pairs 3.82, 3.89'),
        Regime('20251029', '20260805', 'GAIN3.8', 3.83, 37300.0, True,
               'flat pairs 3.81-3.83 (2026-02, -04, -06); flats fill at 39,300-40,800'),
        Regime('20260806', '20260807', 'PED2900', 6.5, 22800.0, False,
               'bias level ~2,900; gain 5.9-7.1 from the same stars, not measured'),
        Regime('20260807', '20991231', 'GAIN3.8', 3.83, 37300.0, True,
               'back in the 3.8 e-/ADU mode during 2026-08-07; 3.89 on 2026-08-08'),
    ],
}

# Measurement settings (ptc.py)
REGION_HALF = 256          # central 512x512 pixels of the trimmed frame
TILE = 64
# Flat-pair levels in ADU above the zero level.  Pair gains rise 2-4% from
# 20,000 to 39,000 ADU in the normal mode (charge spreading at high signal), so
# one window keeps nights comparable; Artemis's 3.8 e-/ADU flats sit at
# 29,000-31,000.  Europa's 3.8 e-/ADU gains climb from 3.8 to 4.1 above
# ~30,000 ADU as the pixels fill, hence the cap below each mode's saturation.
FLAT_MIN_LEVEL = 3000.0
FLAT_MAX_LEVEL = 33000.0
FLAT_SAT_FRACTION = 0.85   # of (saturation - zero level)
FLAT_MAX_RATIO = 0.15      # |m1/m2 - 1| for a usable pair
MIN_GAIN_PAIRS = 3
MAX_GAIN_SCATTER = 0.10    # robust relative scatter of the pair gains
MAX_BIAS = 20
MAX_FLATS = 80
# a frame's own overscan noise replaces the night's bias-pair read noise when
# the two differ by more than this factor (elevated-noise frames)
FRAME_RON_RATIO = 1.3
# a measured gain this far from the table is used, but reported
TABLE_WARN_RATIO = 1.25


# --------------------------------------------------------------- frame stats

def has_overscan(params):
    return bool(params) and 'overscan' in params and params.get('naxis2') is not None


def overscan_strip(data, params):
    """The configured overscan strip of a raw frame, or None.

    Only a frame of the configured full height has one: Artemis wrote
    2048 x 2048 frames without overscan in 2021-04/06, and the same slice of
    those is sky.
    """
    if data is None or not has_overscan(params) or np.ndim(data) != 2:
        return None
    if data.shape[0] != params['naxis2']:
        return None
    cfg = params['overscan']
    strip = data[slice(cfg.get('top_row'), cfg.get('bottom_row')),
                 slice(cfg.get('left_col'), cfg.get('right_col'))]
    if strip.size < TILE * 2:
        return None
    return np.asarray(strip, dtype=float)


def tile_sigma(diff, tile_y=TILE, tile_x=TILE):
    """Median over tiles of the sigma-clipped standard deviation.

    All tiles are clipped at once (one row each), which gives the same result
    as clipping them one by one in a fraction of the time.
    """
    ny, nx = diff.shape
    tile_y, tile_x = min(tile_y, ny), min(tile_x, nx)
    my, mx = ny // tile_y, nx // tile_x
    if my == 0 or mx == 0:
        return float('nan')
    tiles = np.asarray(diff[:my * tile_y, :mx * tile_x], dtype=float)
    tiles = tiles.reshape(my, tile_y, mx, tile_x).transpose(0, 2, 1, 3).reshape(my * mx, tile_y * tile_x)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, _, sig = sigma_clipped_stats(tiles, sigma=4, maxiters=3, axis=1)
    sig = np.asarray(sig, dtype=float)
    sig = sig[np.isfinite(sig)]
    return float(np.median(sig)) if sig.size else float('nan')


def overscan_signature(data, params):
    """(level, noise) of a raw frame's overscan strip in ADU, or (None, None).

    The noise is taken in tiles spanning the strip's height, so a slow drift
    of the bias level along the strip is not counted.
    """
    return strip_signature(overscan_strip(data, params))


def strip_signature(strip):
    """(level, noise) of an overscan strip already cut out of a frame (the
    calibration register reads only those rows), or (None, None)."""
    if strip is None:
        return None, None
    strip = np.asarray(strip, dtype=float)
    if strip.ndim != 2 or strip.size < TILE * 2:
        return None, None
    level = float(np.median(strip))
    noise = tile_sigma(strip, strip.shape[0], TILE)
    if not np.isfinite(noise):
        return level, None
    return level, noise


def trimmed(data, params):
    cfg = (params or {}).get('trim')
    if not cfg or data.shape[0] != params.get('naxis2'):
        return data
    return data[slice(cfg.get('top_row'), cfg.get('bottom_row')),
                slice(cfg.get('left_col'), cfg.get('right_col'))]


def central_region(data, params, half=REGION_HALF):
    img = trimmed(data, params)
    cy, cx = img.shape[0] // 2, img.shape[1] // 2
    return np.asarray(img[cy - half:cy + half, cx - half:cx + half], dtype=float)


# ----------------------------------------------------------- classification

@functools.lru_cache(maxsize=1)
def _telescope_aliases():
    try:
        from calibration.pipeutils import load_instrument_config
        return dict(load_instrument_config().get('telescope_aliases', {}))
    except Exception:
        return {}


def canonical_telescope(telescop):
    """The configured name for a TELESCOP value ('ACP->Artemis' is Artemis)."""
    name = str(telescop or '').strip()
    return _telescope_aliases().get(name, name)


def regimes_for(telescope):
    return REGIMES.get(canonical_telescope(telescope), [])


def expected_regime(telescope, night):
    """The REGIMES entry for this telescope and night, or None."""
    night = str(night)
    match = None
    for r in regimes_for(telescope):
        if r.first <= night <= r.last:
            if match is None or (match.mode == 'NORMAL' and r.mode != 'NORMAL'):
                match = r
    return match


def regime_for_mode(telescope, mode, night):
    """The telescope's entry for a mode nearest in time to the night, or None."""
    night = str(night)
    best, best_gap = None, None
    for r in regimes_for(telescope):
        if r.mode != mode:
            continue
        if r.first <= night <= r.last:
            gap = 0
        else:
            gap = min(abs(_days(night) - _days(r.first)), abs(_days(night) - _days(r.last)))
        if best is None or gap < best_gap:
            best, best_gap = r, gap
    return best


def _days(night):
    try:
        return int(np.datetime64('%s-%s-%s' % (night[:4], night[4:6], night[6:8]), 'D').astype(int))
    except (ValueError, TypeError):
        return 0


def classify(noise, level, telescope, night):
    """Readout mode of a frame from its overscan noise and level.

    The mode the dated table expects is kept unless the frame clearly says
    otherwise; any other mode needs a closer match.  'UNKNOWN' when nothing
    fits (e.g. Europa's 2017-12 bias noise of 35-146 ADU).  Without an
    overscan measurement the expected mode is returned.
    """
    regime = expected_regime(telescope, night)
    expected = regime.mode if regime else 'NORMAL'
    if noise is None or level is None or not np.isfinite(noise) or noise <= 0:
        return expected

    def fits_mode(name, tolerance):
        mode = MODES[name]
        lo, hi = mode.pedestal
        if not lo <= level <= hi:
            return None
        ratio = max(noise / mode.noise, mode.noise / noise)
        return ratio if ratio <= tolerance else None

    keep = max(MODES[expected].tolerance, EXPECTED_TOLERANCE.get(expected, EXPECTED_TOLERANCE_DEFAULT))
    if fits_mode(expected, keep) is not None:
        return expected
    used = {r.mode for r in regimes_for(telescope)}
    best, best_ratio = 'UNKNOWN', None
    for name, mode in MODES.items():
        if name == expected or not (mode.anywhere or name in used):
            continue
        ratio = fits_mode(name, mode.tolerance)
        if ratio is not None and (best_ratio is None or ratio < best_ratio):
            best, best_ratio = name, ratio
    return best


# ------------------------------------------------------------- measurements

def own_night(path, night):
    """True for a frame filed under the night itself (not borrowed)."""
    return ('/%s/' % night) in path.replace('\\', '/')


def _read(path):
    from calibration.pipeutils import open_fits_file
    with open_fits_file(path) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
        if data is None and len(hdul) > 1:
            hdr, data = hdul[1].header, hdul[1].data
        return hdr.copy(), (None if data is None else np.array(data, dtype=float))


def _frame_info(path, params, telescope, night, keep_region):
    hdr, data = _read(path)
    if data is None or np.ndim(data) != 2:
        return None
    level, noise = overscan_signature(data, params)
    info = {
        'path': path,
        't': str(hdr.get('DATE-OBS', '')),
        'filter': str(hdr.get('FILTER', '')).replace("'", '').strip(),
        'exptime': float(hdr.get('EXPTIME', 0) or 0),
        'shape': tuple(data.shape),
        'os_level': level,
        'os_noise': noise,
        'mode': classify(noise, level, telescope, night),
    }
    if keep_region:
        # float32 holds the integer ADU exactly and halves the memory of 80 flats
        region = central_region(data, params).astype(np.float32)
        info['region'] = region
        info['median'] = float(np.median(region))
        info['p999'] = float(np.percentile(region, 99.9))
    return info


def _read_list(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def mode_saturation(telescope, mode, night):
    regime = regime_for_mode(telescope, mode, night)
    if regime is not None and regime.saturation:
        return float(regime.saturation)
    if mode in MODES:
        return float(MODES[mode].saturation)
    return float(MODES['NORMAL'].saturation)


def measure_bias(frames):
    """Read noise (ADU) of each mode from consecutive bias pairs of that mode."""
    out = {}
    by_mode = {}
    for f in sorted(frames, key=lambda f: f['t']):
        by_mode.setdefault(f['mode'], []).append(f)
    for mode, fl in by_mode.items():
        rons, offsets = [], []
        for a, b in zip(fl[0::2], fl[1::2]):
            if a['shape'] != b['shape']:
                continue
            diff = a['region'].astype(float) - b['region'].astype(float)
            rons.append(tile_sigma(diff) / math.sqrt(2.0))
        for f in fl:
            if f['os_level'] is not None:
                offsets.append(f['median'] - f['os_level'])
        noises = [f['os_noise'] for f in fl if f['os_noise'] is not None]
        rons = [r for r in rons if np.isfinite(r)]
        out[mode] = {
            'n_bias': len(fl),
            'ron_adu': float(np.median(rons)) if rons else None,
            'ron_pairs': [round(r, 3) for r in rons],
            'bias_offset': float(np.median(offsets)) if offsets else None,
            'os_noise': float(np.median(noises)) if noises else None,
            'os_level': float(np.median([f['os_level'] for f in fl if f['os_level'] is not None]))
            if any(f['os_level'] is not None for f in fl) else None,
        }
    return out


def measure_gain(frames, bias_by_mode, telescope, night):
    """Gain (e-/ADU) of each mode from consecutive same-filter flat pairs."""
    out = {}
    by_key = {}
    for f in sorted(frames, key=lambda f: f['t']):
        by_key.setdefault((f['mode'], f['filter']), []).append(f)
    pairs_by_mode = {}
    for (mode, _), fl in by_key.items():
        bias = bias_by_mode.get(mode, {})
        offset = bias.get('bias_offset') or 0.0
        saturation = mode_saturation(telescope, mode, night)
        for a, b in zip(fl[:-1], fl[1:]):
            if a['shape'] != b['shape'] or a['os_level'] is None or b['os_level'] is None:
                continue
            za, zb = a['os_level'] + offset, b['os_level'] + offset
            hi = min(FLAT_MAX_LEVEL, FLAT_SAT_FRACTION * (saturation - min(za, zb)))
            m1, m2 = a['median'] - za, b['median'] - zb
            if not (FLAT_MIN_LEVEL < m1 < hi and FLAT_MIN_LEVEL < m2 < hi):
                continue
            if max(a['p999'], b['p999']) >= saturation:
                continue
            k = m1 / m2
            if abs(k - 1.0) > FLAT_MAX_RATIO:
                continue
            # read noise of the pair: the night's bias pairs of this mode, else
            # the flats' own overscan (it hardly matters at these levels)
            r = bias.get('ron_adu')
            if r is None:
                noises = [n for n in (a['os_noise'], b['os_noise']) if n is not None]
                r = float(np.median(noises)) if noises else 0.0
            s = tile_sigma((a['region'].astype(float) - za) - k * (b['region'].astype(float) - zb))
            denom = s * s - (1.0 + k * k) * r * r
            if not np.isfinite(denom) or denom <= 0:
                continue
            g = (m1 + k * k * m2) / denom
            pairs_by_mode.setdefault(mode, []).append(
                {'filter': a['filter'], 'm1': round(m1), 'm2': round(m2), 'gain': round(g, 4)})
    for mode, pairs in pairs_by_mode.items():
        gains = np.array([p['gain'] for p in pairs])
        med = float(np.median(gains))
        scatter = float(1.4826 * np.median(np.abs(gains - med)) / med) if med > 0 else float('inf')
        accepted = len(gains) >= MIN_GAIN_PAIRS and scatter <= MAX_GAIN_SCATTER
        out[mode] = {
            'gain': med,
            'gain_err': float(1.2533 * scatter * med / math.sqrt(len(gains))),
            'gain_pairs': len(gains),
            'gain_scatter': scatter,
            'gain_accepted': bool(accepted),
            'pairs': pairs,
        }
    return out


def measure_night(bias_files, flat_files, night, telescope, params, log=print):
    """The night's read noise and gain per readout mode.

    `bias_files`/`flat_files` are the T1 lists; only frames filed under the
    night itself are used.  Returns a JSON-serialisable record.
    """
    night = str(night)
    record = {
        'version': RECORD_VERSION,
        'telescope': canonical_telescope(telescope),
        'night': night,
        'andor': has_overscan(params),
        'config_gain': _float_or_none((params or {}).get('gain')),
        'config_saturation': _float_or_none((params or {}).get('saturation_threshold')),
        'expected_mode': None,
        'modes': {},
        'calib_mode': None,
        'notes': [],
    }
    if not valid_night(night):
        # without the night's label neither the dated table nor the own-night
        # test of the calibration frames can work
        log('readout: %r is not a YYYYMMDD night; nothing is measured' % night)
        record['notes'].append('night label %r is not YYYYMMDD' % night)
    regime = expected_regime(telescope, night)
    record['expected_mode'] = regime.mode if regime else ('NORMAL' if record['andor'] else 'DEFAULT')
    own_bias = [p for p in bias_files if own_night(p, night)]
    own_flats = [p for p in flat_files if own_night(p, night)]
    record['n_bias_listed'], record['n_bias_own'] = len(bias_files), len(own_bias)
    record['n_flats_listed'], record['n_flats_own'] = len(flat_files), len(own_flats)
    if len(own_bias) < len(bias_files) or len(own_flats) < len(flat_files):
        record['notes'].append('calibration frames borrowed from other nights are not measured')
    if not record['andor']:
        return record

    biases, flats = [], []
    for p in own_bias[:MAX_BIAS]:
        try:
            info = _frame_info(p, params, telescope, night, keep_region=True)
        except Exception as e:  # a corrupt frame must not stop T2
            log('readout: cannot read %s: %s' % (p, e))
            continue
        if info is not None:
            biases.append(info)
    for p in own_flats[:MAX_FLATS]:
        try:
            info = _frame_info(p, params, telescope, night, keep_region=True)
        except Exception as e:
            log('readout: cannot read %s: %s' % (p, e))
            continue
        if info is not None:
            flats.append(info)

    bias_by_mode = measure_bias(biases)
    gain_by_mode = measure_gain(flats, bias_by_mode, telescope, night)
    for mode in sorted(set(bias_by_mode) | set(gain_by_mode) | {f['mode'] for f in flats}):
        entry = {'n_flats': sum(1 for f in flats if f['mode'] == mode)}
        entry.update(bias_by_mode.get(mode, {'n_bias': 0, 'ron_adu': None}))
        entry.update(gain_by_mode.get(mode, {'gain': None, 'gain_pairs': 0, 'gain_accepted': False}))
        record['modes'][mode] = entry
        table = regime_for_mode(telescope, mode, night)
        if entry.get('gain_accepted') and table is not None and table.gain:
            ratio = entry['gain'] / table.gain
            if not 1.0 / TABLE_WARN_RATIO <= ratio <= TABLE_WARN_RATIO:
                record['notes'].append('%s gain %.3f measured, table says %.3f'
                                       % (mode, entry['gain'], table.gain))
    # the mode of the master calibrations: the biases', else the flats'
    counts = {}
    for f in biases or flats:
        counts[f['mode']] = counts.get(f['mode'], 0) + 1
    record['calib_mode'] = max(counts, key=counts.get) if counts else record['expected_mode']
    seen = sorted(counts)
    if any(m != record['expected_mode'] for m in seen):
        record['notes'].append('calibration frames in %s; table expects %s'
                               % ('+'.join(seen), record['expected_mode']))
    return record


def _float_or_none(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


# --------------------------------------------------------------- resolution

def mode_gain(record, mode):
    """(gain, source) for a mode: the night's flats ('night'), the dated table
    ('table', or 'estimate' for a gain that was never measured), the config."""
    modes = record.get('modes', {})
    entry = modes.get(mode) or {}
    if entry.get('gain_accepted') and entry.get('gain'):
        return float(entry['gain']), 'night'
    telescope, night = record.get('telescope'), record.get('night')
    if mode == 'UNKNOWN':
        # an unrecognised readout keeps the gain of the mode the table expects
        mode = record.get('expected_mode') or 'NORMAL'
        entry = modes.get(mode) or {}
        if entry.get('gain_accepted') and entry.get('gain'):
            return float(entry['gain']), 'night'
    regime = regime_for_mode(telescope, mode, night)
    if regime is not None and regime.gain:
        return float(regime.gain), ('table' if regime.measured else 'estimate')
    # a telescope the table knows, in a mode it has not used here before (a
    # telescope it does not know at all keeps its config gain)
    if regime is None and regimes_for(telescope) and mode in MODES and MODES[mode].gain:
        return float(MODES[mode].gain), 'table'
    gain = record.get('config_gain')
    return (float(gain) if gain else 1.0), 'config'


def frame_values(record, os_level, os_noise, legacy_ron_e=None, mode=None):
    """GAIN, RDNOISE (e-), raw saturation (ADU) and provenance for one frame.

    `mode` overrides the classification from the overscan.  `legacy_ron_e` is
    the pipeline's original read-noise value, used for instruments this
    module does not handle.
    """
    if not record.get('andor'):
        gain = record.get('config_gain') or 1.0
        return {
            'READMODE': 'DEFAULT',
            'GAIN': float(gain),
            'GAINSRC': 'config',
            'RDNOISE': _float_or_none(legacy_ron_e),
            'RONSRC': 'legacy',
            'SATURATE_RAW': record.get('config_saturation'),
            'OSNOISE': None,
        }
    telescope, night = record.get('telescope'), record.get('night')
    if mode is None:
        mode = classify(os_noise, os_level, telescope, night)
    gain, gain_src = mode_gain(record, mode)
    entry = record.get('modes', {}).get(mode) or {}
    # the night's bias pairs of this mode, unless this frame's own overscan
    # says its noise was different (Artemis 2023 varied within a night)
    ron_adu, ron_src = entry.get('ron_adu'), 'bias'
    if ron_adu is not None and os_noise is not None and entry.get('os_noise'):
        ratio = os_noise / entry['os_noise']
        if not 1.0 / FRAME_RON_RATIO <= ratio <= FRAME_RON_RATIO:
            ron_adu, ron_src = os_noise, 'overscan'
    if ron_adu is None and os_noise is not None:
        ron_adu, ron_src = os_noise, 'overscan'
    if ron_adu is None:
        base = mode if mode in MODES else (record.get('expected_mode') or 'NORMAL')
        ron_adu, ron_src = MODES.get(base, MODES['NORMAL']).noise, 'table'
    sat_mode = mode if mode != 'UNKNOWN' else (record.get('expected_mode') or 'NORMAL')
    return {
        'READMODE': mode,
        'GAIN': float(gain),
        'GAINSRC': gain_src,
        'RDNOISE': float(ron_adu) * float(gain),
        'RONSRC': ron_src,
        'SATURATE_RAW': mode_saturation(telescope, sat_mode, night),
        'OSNOISE': None if os_noise is None else float(os_noise),
    }


def calib_values(record, legacy_ron_e=None):
    """Values for the night's master calibrations: readoutnoise.dat, QC6 and
    the dark current, all in the mode of the night's own bias frames."""
    if not record.get('andor'):
        return frame_values(record, None, None, legacy_ron_e)
    mode = record.get('calib_mode') or record.get('expected_mode') or 'NORMAL'
    entry = record.get('modes', {}).get(mode) or {}
    return frame_values(record, entry.get('os_level'), entry.get('os_noise'), legacy_ron_e, mode=mode)


# ------------------------------------------------------------ record on disk

def record_path(reddir):
    return os.path.join(reddir, RECORD_NAME)


def write_record(record, reddir):
    path = record_path(reddir)
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(record, f, indent=1, sort_keys=True, default=_json_default)
    os.replace(tmp, path)
    return path


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(repr(o))


def read_record(reddir):
    path = record_path(reddir)
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            record = json.load(f)
    except (OSError, ValueError):
        return None
    if record.get('version') != RECORD_VERSION:
        return None
    return record


def night_lists(reddir, runname):
    """The T1 bias list and flat lists of a reduction directory."""
    bias = _read_list(os.path.join(reddir, '%s_bias.list' % runname))
    flats = []
    for fl in sorted(glob.glob(os.path.join(reddir, '%s_flat_*.list' % runname))):
        flats += _read_list(fl)
    return bias, flats


def night_from_reddir(reddir):
    """YYYYMMDD from .../output/<night>/reduction, the pipeline's layout."""
    parts = [p for p in os.path.normpath(reddir).split(os.sep) if p]
    if len(parts) >= 2 and parts[-1] == 'reduction':
        return parts[-2]
    return parts[-1] if parts else ''


def valid_night(night):
    return len(str(night)) == 8 and str(night).isdigit()


def load_or_measure(reddir, runname, params, telescope, night=None, refresh=False, write=True,
                    log=print):
    """The night's record: readout.json if T2 wrote one for this night, else
    measured now (and saved).  T2 passes refresh=True.

    A measurement that fails leaves a record with no modes, so every frame
    still gets its mode from its own overscan and its values from the table.
    """
    night = str(night or night_from_reddir(reddir))
    telescope = canonical_telescope(telescope)
    if not refresh:
        record = read_record(reddir)
        if record is not None and record.get('night') == night \
                and record.get('telescope') == telescope:
            return record
    try:
        bias, flats = night_lists(reddir, runname)
        record = measure_night(bias, flats, night, telescope, params, log=log)
    except Exception as e:  # never let this stop the calibration
        log('readout: measurement failed for %s %s: %s' % (telescope, night, e))
        record = measure_night([], [], night, telescope, params, log=log)
        record['notes'].append('measurement failed: %s' % e)
    if write:
        try:
            write_record(record, reddir)
        except OSError as e:
            log('readout: could not write %s: %s' % (record_path(reddir), e))
    return record


def summary(record):
    lines = ['%s %s: expected %s, calibration frames %s'
             % (record.get('telescope'), record.get('night'), record.get('expected_mode'),
                record.get('calib_mode'))]
    for mode, e in sorted(record.get('modes', {}).items()):
        gain, src = mode_gain(record, mode)
        lines.append('  %-8s bias %2d  RON %s ADU  flats %2d  gain %.3f (%s%s)'
                     % (mode, e.get('n_bias', 0),
                        '%.2f' % e['ron_adu'] if e.get('ron_adu') is not None else '  - ',
                        e.get('n_flats', 0), gain, src,
                        ', %d pairs' % e['gain_pairs'] if e.get('gain_pairs') else ''))
    for n in record.get('notes', []):
        lines.append('  note: ' + n)
    return '\n'.join(lines)


def header_cards(values, pedestal=0.0):
    """FITS cards for a calibrated frame.  `pedestal` is the level subtracted
    from the raw frame (overscan + master-bias median), so SATURATE is in the
    calibrated frame's units."""
    cards = [
        ('GAIN', values['GAIN'], '[e-/ADU] gain (%s)' % values['GAINSRC']),
        ('GAINSRC', values['GAINSRC'], 'night flats, table, estimate or config'),
        ('READMODE', values['READMODE'], 'readout mode inferred from the overscan'),
    ]
    if values.get('RDNOISE') is not None:
        cards += [('RDNOISE', values['RDNOISE'], '[e-] read noise (%s)' % values['RONSRC']),
                  ('RON', values['RDNOISE'], 'Read out noise (e-s)')]
    if values.get('SATURATE_RAW') is not None:
        cards += [('SATURATE', values['SATURATE_RAW'] - pedestal,
                   '[ADU] saturation in this frame (raw %.0f)' % values['SATURATE_RAW'])]
    if values.get('OSNOISE') is not None:
        cards += [('OSNOISE', values['OSNOISE'], '[ADU] overscan noise of the raw frame')]
    return cards


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('reddir')
    parser.add_argument('--night')
    parser.add_argument('--runname', default='1')
    parser.add_argument('--json', action='store_true', help='print the whole record')
    args = parser.parse_args(argv)
    from calibration.pipeutils import get_instrument_parameters, open_fits_file, frame_telescope
    bias, flats = night_lists(args.reddir, args.runname)
    first = (bias or flats or [None])[0]
    if first is None:
        print('no calibration lists in %s' % args.reddir)
        return 1
    with open_fits_file(first) as hdul:
        params = get_instrument_parameters(hdul)
        telescope = frame_telescope(hdul)
    record = measure_night(bias, flats, args.night or night_from_reddir(args.reddir), telescope, params)
    print(json.dumps(record, indent=1, default=_json_default) if args.json else summary(record))
    return 0


if __name__ == '__main__':
    sys.exit(main())
