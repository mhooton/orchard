"""
Readout-smear correction for shutterless CCD frames.

Without a mechanical shutter the array keeps integrating while the frame is
clocked out, so each charge packet collects light from every pixel it crosses
on its way to the serial register. A frame read out in N parallel shifts over
t_read seconds therefore records, at each pixel,

    raw[k] = true[k] + tau * sum_{j < k} true[j],      tau = t_read / (N * t_exp)

where k runs along the parallel (readout) axis away from the serial register.
Every source leaves a one-sided trail of tau times its own flux in each pixel
downstream of it. The model is inverted exactly, row by row:

    true[k] = raw[k] - tau * C[k],     C[k + 1] = (1 - tau) * C[k] + raw[k]

which is a first-order IIR filter, so the whole frame is one lfilter call.

The inverse is exact only where nothing upstream is clipped. A saturated pixel
has lost the charge that predicts its trail, so every pixel downstream of it in
the same row is under-corrected. saturated_trail_mask() marks those pixels so
that stars sitting on them can be flagged rather than trusted.

The mask itself is never stored. correct_frame() writes its compact form into
the frame header instead: each row holding a clipped pixel, with the position
of the first clip in readout order (SMTR0001, SMTR0002, ...; see
write_trail_header). The photometry reads it back with apertures_on_trail(),
which gives, for each star, the smallest aperture that reaches a trail pixel.

On SPECULOOS ANDOR frames as stored, the parallel axis is NAXIS1 (numpy axis 1)
and trails run toward +x. The parameters live under "smear" in
instrument_config.json. t_read = 4.33 s was measured on sky in August 2026 by
requiring 100% trail removal on Ganymede 20250408 (TOI-1227, 20 s). A
shutterless twilight flat gave 4.58 s and the Andor 1 MHz spec is 4.53 s; the
difference is unresolved and is what two shutterless flats at different
exposure times would settle.

Whether a frame was taken without the shutter comes from the CAMSHUTT header
(written by Astra on all ANDORs since 2026-07-17, but set by hand, so it can be
wrong) unless shutter_overrides.csv says otherwise for that telescope and time.
A CAMSHUTT value that cannot be read is reported and taken as shutter used.

The readout rate is not in the raw header either. Ganymede was read out in
~1.5 s, consistent with 3 MHz, for its last shutterless nights in
February-March 2026. T6 now writes the readout mode it recognises in each
frame's overscan as READMODE (calibration/readout.py), and from 2026-10-18
Astra writes the camera's own readout rate as READOUTM. So t_read comes from an
override window that sets one, else from READOUTM ("readout_rate_time" under
"smear": 1MHz 4.33 s, 3MHz 1.5 s), else from READMODE ("readout_time": NORMAL
4.33 s, FAST 1.5 s). A shutterless frame with no known readout time is left
smeared, says so in SMEARERR, and a warning is printed: a wrong t_read would
put a ramp into it instead.

A pixel counts as clipped at the frame's own saturation level, the SATURATE
card T6 writes (relative to the overscan + master-bias pedestal), which in
Ganymede's fast and 7,500-pedestal modes sits far below the instrument
config's 62,000 ADU; frames without it fall back to saturation_threshold.
"""

import csv
import functools
import os
from datetime import datetime

import numpy as np
from scipy.signal import lfilter

OVERRIDES_PATH = os.path.join(os.path.dirname(__file__), 'shutter_overrides.csv')

_warned = set()


def _warn(message):
    """Print a warning once per process: correct_frame runs once per frame."""
    if message not in _warned:
        _warned.add(message)
        print('WARNING: ' + message, flush=True)


def smear_tau(readout_time, n_shifts, exptime):
    """Fraction of a pixel's exposure signal added to each packet that crosses it."""
    if exptime <= 0:
        raise ValueError(f'exposure time must be positive, got {exptime}')
    return readout_time / (n_shifts * exptime)


def _along_axis(image, axis, direction):
    """View with the readout axis last and trails running toward higher index."""
    moved = np.moveaxis(image, axis, -1)
    return moved[..., ::-1] if direction < 0 else moved


def _from_axis(moved, axis, direction):
    if direction < 0:
        moved = moved[..., ::-1]
    return np.moveaxis(moved, -1, axis)


def add_smear(true, tau, axis=1, direction=1):
    """Forward model: what a shutterless readout records for a true image."""
    t = _along_axis(np.asarray(true, dtype=np.float64), axis, direction)
    upstream = np.cumsum(t, axis=-1) - t
    return _from_axis(t + tau * upstream, axis, direction)


def desmear(image, tau, axis=1, direction=1):
    """
    Remove readout smear from a bias- and dark-subtracted frame.

    Must be applied before flat-fielding: the smear charge was generated in the
    upstream pixels, with their sensitivity, not the sensitivity of the pixel it
    ends up in. NaN pixels are treated as zero signal when accumulating the
    trail and stay NaN in the output. A masked array keeps its mask.

    Args:
        image: 2-D frame in ADU (float or int; masked arrays are accepted)
        tau: smear fraction from smear_tau()
        axis: numpy axis along which charge is clocked out
        direction: +1 if trails run toward higher index along axis, -1 if lower

    Returns:
        float64 array (masked if the input was masked) with the smear removed
    """
    mask = np.ma.getmask(image) if np.ma.isMaskedArray(image) else None
    raw = np.array(np.ma.getdata(image), dtype=np.float64)
    bad = ~np.isfinite(raw)
    raw[bad] = 0.0

    r = _along_axis(raw, axis, direction)
    # y[k] = (1 - tau) * y[k-1] + r[k] is the summed true signal up to and
    # including k, so the upstream sum at k is y shifted by one.
    summed = lfilter([1.0], [1.0, -(1.0 - tau)], r, axis=-1)
    upstream = np.zeros_like(r)
    upstream[..., 1:] = summed[..., :-1]
    out = np.ascontiguousarray(_from_axis(r - tau * upstream, axis, direction))

    out[bad] = np.nan
    if mask is not None and mask is not np.ma.nomask:
        return np.ma.array(out, mask=mask)
    return out


def desmear_loop(image, tau, axis=1, direction=1):
    """
    Running-subtraction form of desmear(), one pixel at a time along the axis.
    Slow; kept as an independent check on the IIR form in the tests.
    """
    r = _along_axis(np.array(image, dtype=np.float64), axis, direction)
    out = np.empty_like(r)
    running = np.zeros(r.shape[:-1])
    for k in range(r.shape[-1]):
        out[..., k] = r[..., k] - tau * running
        running = running + out[..., k]
    return _from_axis(out, axis, direction)


def saturated_trail_mask(raw, saturation_threshold, axis=1, direction=1):
    """
    Pixels downstream of a saturated pixel in the same row.

    desmear() under-corrects these: the clipped pixel's lost charge was never
    recorded, so its trail is predicted too faint. A star whose aperture falls
    on such a pixel should be flagged, not trusted.

    Args:
        raw: frame in raw ADU (before overscan or bias subtraction), so the
            threshold compares like with like
        saturation_threshold: ADU at or above which a pixel counts as clipped

    Returns:
        boolean array, True where the smear correction is unreliable
    """
    clipped = _clipped(raw, saturation_threshold, axis, direction)
    downstream = np.logical_or.accumulate(clipped, axis=-1)
    downstream[..., 1:] = downstream[..., :-1].copy()
    downstream[..., 0] = False
    return _from_axis(downstream, axis, direction)


def _clipped(raw, saturation_threshold, axis, direction):
    return _along_axis(np.ma.filled(np.asarray(raw), 0) >= saturation_threshold, axis, direction)


def clipped_lines(raw, saturation_threshold, axis=1, direction=1):
    """
    The compact form of saturated_trail_mask() for a 2-D frame.

    A line is a row of pixels along the readout axis (an image row on the
    ANDORs). In line lines[k], every pixel beyond position first[k] in the
    direction the trails run is under-corrected, and nothing else is. A line
    whose only clip is its last pixel has no trail and is left out.

    Returns:
        (lines, first): int arrays of 0-based indices, lines across the
        readout axis and first along it
    """
    clipped = _clipped(raw, saturation_threshold, axis, direction)[..., :-1]
    lines = np.flatnonzero(clipped.any(axis=-1))
    first = np.argmax(clipped[lines], axis=-1)
    if direction < 0:
        first = clipped.shape[-1] - first
    return lines, first


# Header cards holding clipped_lines(). Each SMTRnnnn value is one run of
# consecutive lines, 'L:P,P,...': L is the first line of the run and each P the
# first clipped position in the next line, all 1-based FITS pixel numbers. A
# long run carries on in the next card. A card value holds at most 68 characters.
TRAIL_CARD = 'SMTR'
_CARD_VALUE_LEN = 68


def _trail_card_keys(header):
    return [k for k in header if k.startswith(TRAIL_CARD) and k[len(TRAIL_CARD):].isdigit()]


def _encode_runs(lines, first):
    """'L:P,P,...' strings for 1-based lines and positions, each fitting one card."""
    values = []
    value, last = None, None
    for line, pos in zip(lines, first):
        piece = str(int(pos))
        if value is not None and line == last + 1 and len(value) + 1 + len(piece) <= _CARD_VALUE_LEN:
            value += ',' + piece
        else:
            if value is not None:
                values.append(value)
            value = f'{int(line)}:{piece}'
        last = line
    if value is not None:
        values.append(value)
    return values


def write_trail_header(header, lines, first, axis, direction):
    """
    Record clipped_lines() in a frame header, replacing any earlier record.

    SMSATROW holds the number of lines, SMTRAXIS and SMTRDIR the geometry, and
    SMTR0001, SMTR0002, ... the lines themselves. A frame with a few saturated
    stars needs a handful of cards; even a clip in every row of an ANDOR frame
    needs only ~150 (12 kB against 16.7 MB of pixels).
    """
    for key in _trail_card_keys(header):
        del header[key]
    header['SMSATROW'] = (int(len(lines)), 'Rows with a clipped pixel: under-corrected')
    header['SMTRAXIS'] = (int(axis), 'numpy axis the smear trails run along')
    header['SMTRDIR'] = (int(direction), 'Trails run to higher (1) or lower (-1) index')
    for n, value in enumerate(_encode_runs(np.asarray(lines) + 1, np.asarray(first) + 1), start=1):
        header[f'{TRAIL_CARD}{n:04d}'] = value


def read_trail_header(header):
    """
    The clipped lines a de-smeared frame's header records.

    Returns:
        (lines, first, axis, direction), with 0-based lines and first as from
        clipped_lines(); empty arrays for a frame that was not de-smeared, which
        has no under-corrected pixels; None for a de-smeared frame whose
        header does not say (a frame corrected without the raw data at hand)
    """
    if not header.get('SMEARCOR', False):
        return np.zeros(0, int), np.zeros(0, int), 1, 1
    if 'SMTRAXIS' not in header:
        return None
    lines, first = [], []
    for key in sorted(_trail_card_keys(header)):
        start, positions = str(header[key]).split(':')
        positions = [int(p) for p in positions.split(',')]
        lines.extend(range(int(start), int(start) + len(positions)))
        first.extend(positions)
    return (np.array(lines, dtype=int) - 1, np.array(first, dtype=int) - 1,
            int(header['SMTRAXIS']), int(header['SMTRDIR']))


def trail_distance(x, y, lines, first, axis=1, direction=1):
    """
    For each star, how large an aperture can be before it reaches a trail pixel.

    An aperture of radius r centred on (x, y) overlaps an under-corrected pixel
    exactly when r exceeds the returned distance. Pixel k covers k - 0.5 to
    k + 0.5, and the aperture is taken to reach a pixel it overlaps at all, as
    the photometry weights partial pixels.

    Args:
        x, y: star positions in 0-based numpy pixel coordinates (x along axis
            1, y along axis 0), i.e. imcore's 1-based X_coordinate minus one
        lines, first, axis, direction: as from clipped_lines()

    Returns:
        float array, inf for a star with no trail pixel anywhere near
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = np.atleast_1d(np.asarray(y, dtype=float))
    if len(lines) == 0:
        return np.full(x.shape, np.inf)
    along, across = (x, y) if axis == 1 else (y, x)
    lines = np.asarray(lines, dtype=float)[None, :]
    first = np.asarray(first, dtype=float)[None, :]
    # gap across the lines between the aperture centre and the nearest edge of each line
    gap = np.maximum(np.abs(lines - across[:, None]) - 0.5, 0.0)
    # gap along the line between the centre and the edge of the first trail pixel
    if direction > 0:
        ahead = (first + 0.5) - along[:, None]
    else:
        ahead = along[:, None] - (first - 0.5)
    return np.min(np.hypot(np.maximum(ahead, 0.0), gap), axis=1)


SMEAR_NOT_ASSESSED = -1


def apertures_on_trail(header, x, y, radii):
    """
    Per star, the smallest aperture that reaches an under-corrected trail pixel.

    Args:
        header: header of the calibrated frame (written by correct_frame)
        x, y: star positions in 0-based numpy pixel coordinates
        radii: aperture radii in pixels, in increasing order

    Returns:
        int16 array: the 1-based index into radii of the smallest aperture on a
        trail, 0 if none is, or SMEAR_NOT_ASSESSED for every star when the
        header cannot say. A star is flagged at aperture j when 0 < value <= j.
    """
    x = np.atleast_1d(x)
    trail = read_trail_header(header)
    if trail is None:
        return np.full(x.shape, SMEAR_NOT_ASSESSED, dtype=np.int16)
    distance = trail_distance(x, y, *trail)
    smallest = np.searchsorted(np.asarray(radii, dtype=float), distance, side='right') + 1
    smallest[smallest > len(radii)] = 0
    return smallest.astype(np.int16)


def _parse_camshutt(value):
    """
    CAMSHUTT is a FITS bool, but it is set by hand, so accept 0/1 and strings too.

    Returns True (shutter used), False (shutterless), or None for a value it
    cannot read; the callers then assume the shutter was used.
    """
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return None if np.isnan(value) else bool(value)
    text = str(value).strip().upper()
    if text in ('T', 'TRUE', '1', 'Y', 'YES'):
        return True
    if text in ('F', 'FALSE', '0', 'N', 'NO'):
        return False
    return None


def _parse_utc(text):
    return datetime.fromisoformat(text.strip().replace('Z', ''))


def load_shutter_overrides(path=OVERRIDES_PATH):
    """Rows of shutter_overrides.csv as dicts; blank lines and # comments skipped."""
    if not os.path.exists(path):
        return []
    with open(path) as f:
        lines = [line for line in f if line.strip() and not line.lstrip().startswith('#')]
    rows = []
    for row in csv.DictReader(lines):
        readout_time = (row.get('readout_time') or '').strip()
        camshutt = _parse_camshutt(row['camshutt'])
        if camshutt is None:
            _warn(f'{path}: unrecognised camshutt {row["camshutt"]!r} for {row["telescope"].strip()} '
                  f'from {row["start_utc"].strip()}; shutter assumed used in that window')
            camshutt = True
        rows.append({
            'telescope': row['telescope'].strip(),
            'start': _parse_utc(row['start_utc']),
            'end': _parse_utc(row['end_utc']),
            'camshutt': camshutt,
            'readout_time': float(readout_time) if readout_time else None,
            'note': (row.get('note') or '').strip(),
        })
    return rows


@functools.lru_cache(maxsize=1)
def _shipped_overrides():
    """shutter_overrides.csv, read once per process rather than once per frame."""
    return tuple(load_shutter_overrides())


def shutter_state(header, overrides=None):
    """
    Was the mechanical shutter used for this frame?

    An override row for this TELESCOP whose [start_utc, end_utc) window holds
    DATE-OBS wins; then the CAMSHUTT header; with neither, the shutter is
    assumed used, which is true of almost every frame taken before July 2026.

    An override row can also give the readout time, for windows where the
    camera was read out at a different rate (nothing in the header says so).

    A CAMSHUTT value that cannot be read is reported and taken as shutter
    used, with source 'CAMSHUTT?'.

    Returns:
        (shutter_used, source, readout_time) with source one of 'override',
        'CAMSHUTT', 'CAMSHUTT?', 'default', and readout_time None unless an
        override sets it
    """
    if overrides is None:
        overrides = _shipped_overrides()
    telescope = str(header.get('TELESCOP', '')).strip()
    date_obs = header.get('DATE-OBS')
    if overrides and date_obs:
        when = _parse_utc(str(date_obs))
        for row in overrides:
            if row['telescope'] == telescope and row['start'] <= when < row['end']:
                return row['camshutt'], 'override', row.get('readout_time')
    if 'CAMSHUTT' in header:
        used = _parse_camshutt(header['CAMSHUTT'])
        if used is None:
            _warn(f'{telescope}: unrecognised CAMSHUTT {header["CAMSHUTT"]!r}; '
                  f'shutter assumed used, readout smear not removed')
            return True, 'CAMSHUTT?', None
        return used, 'CAMSHUTT', None
    return True, 'default', None


def shutterless(header, params, overrides=None):
    """Will correct_frame() treat this frame as shutterless (and has the instrument smear parameters)?"""
    if not params or not params.get('smear'):
        return False
    return not shutter_state(header, overrides)[0]


def mode_readout_time(smear, mode):
    """
    Full-frame readout time in seconds for a READMODE, from the "readout_time"
    table under "smear" in instrument_config.json; None when the mode (or the
    mode itself) is not known.
    """
    times = smear.get('readout_time') or {}
    mode = str(mode or '').strip()
    if not mode or mode not in times:
        return None
    return float(times[mode])


# Astra writes READOUTM on every Andor frame from 2026-10-18 (Astra 1.3; Peter
# Pedersen's email of that day): the camera driver's readout mode, as an index
# into this list. SPIRIT has no such card. The "readout_rate_time" table under
# "smear" gives the readout time of each rate that has been measured: 1 MHz
# 4.33 s (star trails) and 3 MHz 1.5 s (Ganymede 2026-02-27..03-01). The pixel
# clock arithmetic that reproduces both, 2048 rows x (2088 px / f + ~40 us),
# predicts ~0.94 s at 5 MHz and ~86 s at 0.05 MHz; neither has been measured,
# so neither is in the table and such frames are not de-smeared.
READOUTM_RATES = ('5MHz', '3MHz', '1MHz', '0.05MHz')


def readoutm_rate(value):
    """The readout rate a READOUTM index names, or None for anything else."""
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        index = float(value)
    except (TypeError, ValueError):
        return None
    if not index.is_integer() or not 0 <= index < len(READOUTM_RATES):
        return None
    return READOUTM_RATES[int(index)]


def readout_time(header, smear, override_time=None):
    """
    The readout time to de-smear a shutterless frame with.

    In order: an override window's own readout_time; Astra's READOUTM card
    (the camera's readout rate, from 2026-10-18), looked up in the smear
    config's "readout_rate_time"; T6's READMODE (inferred from the overscan,
    calibration/readout.py), looked up in its "readout_time". A frame whose
    READOUTM names a rate without a measured time, or cannot be read, gets
    None rather than falling back to READMODE: the header says what the camera
    did, the overscan only suggests it.

    Returns:
        (seconds or None, source): source is 'override', 'READOUTM <i> (<rate>)'
        or 'READMODE <mode>'
    """
    if override_time is not None:
        return float(override_time), 'override'
    if 'READOUTM' in header:
        value = header['READOUTM']
        rate = readoutm_rate(value)
        if rate is None:
            return None, 'READOUTM %r' % (value,)
        times = smear.get('readout_rate_time') or {}
        return (float(times[rate]) if rate in times else None), 'READOUTM %s (%s)' % (value, rate)
    mode = str(header.get('READMODE', '') or '').strip()
    return mode_readout_time(smear, mode), 'READMODE ' + (mode or '(none)')


def clip_threshold(header, params, pedestal=0.0):
    """
    The raw-ADU level at or above which a pixel counts as clipped.

    T6 writes the frame's own saturation level as SATURATE, relative to the
    pedestal it subtracts (overscan + master-bias median), so the raw level is
    SATURATE + pedestal. A frame without SATURATE falls back to the instrument
    config's saturation_threshold, compared with the raw pixels as before.

    Returns:
        (level or None, source): source is 'SATURATE' or 'config'
    """
    level = header.get('SATURATE')
    try:
        level = float(level)
    except (TypeError, ValueError):
        level = None
    if level is not None and np.isfinite(level):
        return level + float(pedestal), 'SATURATE'
    if params and params.get('saturation_threshold') is not None:
        return float(params['saturation_threshold']), 'config'
    return None, None


def correct_frame(signal, header, params, raw=None, overrides=None, pedestal=0.0):
    """
    De-smear one bias- and dark-subtracted frame if it was taken shutterless.

    Writes the decision into header either way (SMEARCOR), so every product
    records whether it was corrected. Does nothing for instruments without a
    "smear" entry in instrument_config.json. The readout time comes from
    readout_time(), so header should already carry T6's READMODE (and
    SATURATE, for the clipped rows). A shutterless frame whose readout time is
    not known is returned unchanged, with SMEARCOR = F and the reason in
    SMEARERR.

    Args:
        signal: bias/dark-subtracted, trimmed frame, not yet flat-fielded
        header: the frame's FITS header, updated in place
        params: instrument parameters from get_instrument_parameters()
        raw: trimmed raw frame, used only to record the rows with clipped
            pixels (write_trail_header) so the photometry can flag stars on them
        overrides: shutter override rows (shutter_overrides.csv if None)
        pedestal: what was subtracted from raw to make signal's zero
            (overscan + master-bias median), the zero of the SATURATE card

    Returns:
        the de-smeared frame, or signal unchanged
    """
    smear = params.get('smear') if params else None
    if not smear:
        return signal

    shutter_used, source, override_time = shutter_state(header, overrides)
    if shutter_used:
        if source == 'CAMSHUTT?':
            header['SMEARCOR'] = (False, 'Not removed: CAMSHUTT unreadable, shutter used')
        else:
            header['SMEARCOR'] = (False, 'Readout smear not removed (shutter used)')
        return signal

    t_read, t_source = readout_time(header, smear, override_time)
    if t_read is None:
        telescope = str(header.get('TELESCOP', '')).strip()
        header['SMEARCOR'] = (False, 'Shutterless, but readout smear NOT removed')
        header['SMEARSRC'] = (source, 'Where the shutter state came from')
        header['SMEARERR'] = (f'no readout time for {t_source}', 'why the smear was not removed')
        _warn(f'{telescope}: shutterless frames ({source}) with no known readout time for {t_source}; '
              f'readout smear NOT removed (add a measured time to "smear" in instrument_config.json, '
              f'or a readout_time to the override window)')
        return signal

    exptime = float(header['EXPTIME'])
    axis = smear.get('axis', 1)
    direction = smear.get('direction', 1)
    tau = smear_tau(t_read, smear['n_shifts'], exptime)
    corrected = desmear(signal, tau, axis=axis, direction=direction)

    header['SMEARCOR'] = (True, 'Readout smear removed (shutterless frame)')
    header['SMEARTRD'] = (float(t_read), f'[s] Readout time used ({t_source})')
    header['SMEARTAU'] = (float(tau), 'Smear per pixel crossed, t_read/(N t_exp)')
    header['SMEARSRC'] = (source, 'Where the shutter state came from')
    if raw is not None:
        level, level_source = clip_threshold(header, params, pedestal)
        if level is not None:
            lines, first = clipped_lines(raw, level, axis=axis, direction=direction)
            header['SMSATLEV'] = (float(level), f'[ADU] raw level counted as clipped ({level_source})')
            write_trail_header(header, lines, first, axis, direction)
    header.add_history(f'Readout smear removed: t_read {t_read} s ({t_source}), '
                       f'tau {tau:.4g}, shutter state from {source}')
    return corrected
