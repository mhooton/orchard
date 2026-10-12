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
"""

import csv
import os
from datetime import datetime

import numpy as np
from scipy.signal import lfilter

OVERRIDES_PATH = os.path.join(os.path.dirname(__file__), 'shutter_overrides.csv')


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
    clipped = _along_axis(np.ma.filled(np.asarray(raw), 0) >= saturation_threshold, axis, direction)
    downstream = np.logical_or.accumulate(clipped, axis=-1)
    downstream[..., 1:] = downstream[..., :-1].copy()
    downstream[..., 0] = False
    return _from_axis(downstream, axis, direction)


def _parse_camshutt(value):
    """CAMSHUTT is a FITS bool, but it is set by hand, so accept 0/1 and strings too."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    text = str(value).strip().upper()
    if text in ('T', 'TRUE', '1', 'Y', 'YES'):
        return True
    if text in ('F', 'FALSE', '0', 'N', 'NO'):
        return False
    raise ValueError(f'unrecognised CAMSHUTT value {value!r}')


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
        rows.append({
            'telescope': row['telescope'].strip(),
            'start': _parse_utc(row['start_utc']),
            'end': _parse_utc(row['end_utc']),
            'camshutt': _parse_camshutt(row['camshutt']),
            'note': (row.get('note') or '').strip(),
        })
    return rows


def shutter_state(header, overrides=None):
    """
    Was the mechanical shutter used for this frame?

    An override row for this TELESCOP whose [start_utc, end_utc) window holds
    DATE-OBS wins; then the CAMSHUTT header; with neither, the shutter is
    assumed used, which is true of almost every frame taken before July 2026.

    Returns:
        (shutter_used, source) with source one of 'override', 'CAMSHUTT', 'default'
    """
    if overrides is None:
        overrides = load_shutter_overrides()
    telescope = str(header.get('TELESCOP', '')).strip()
    date_obs = header.get('DATE-OBS')
    if overrides and date_obs:
        when = _parse_utc(str(date_obs))
        for row in overrides:
            if row['telescope'] == telescope and row['start'] <= when < row['end']:
                return row['camshutt'], 'override'
    if 'CAMSHUTT' in header:
        return _parse_camshutt(header['CAMSHUTT']), 'CAMSHUTT'
    return True, 'default'


def correct_frame(signal, header, params, raw=None, overrides=None):
    """
    De-smear one bias- and dark-subtracted frame if it was taken shutterless.

    Writes the decision into header either way (SMEARCOR), so every product
    records whether it was corrected. Does nothing for instruments without a
    "smear" entry in instrument_config.json.

    Args:
        signal: bias/dark-subtracted, trimmed frame, not yet flat-fielded
        header: the frame's FITS header, updated in place
        params: instrument parameters from get_instrument_parameters()
        raw: trimmed raw frame, used only to count rows with clipped pixels
        overrides: shutter override rows (loaded from disk if None)

    Returns:
        the de-smeared frame, or signal unchanged
    """
    smear = params.get('smear') if params else None
    if not smear:
        return signal

    shutter_used, source = shutter_state(header, overrides)
    if shutter_used:
        header['SMEARCOR'] = (False, 'Readout smear not removed (shutter used)')
        return signal

    exptime = float(header['EXPTIME'])
    axis = smear.get('axis', 1)
    direction = smear.get('direction', 1)
    tau = smear_tau(smear['readout_time'], smear['n_shifts'], exptime)
    corrected = desmear(signal, tau, axis=axis, direction=direction)

    header['SMEARCOR'] = (True, 'Readout smear removed (shutterless frame)')
    header['SMEARTRD'] = (float(smear['readout_time']), '[s] Readout time used for smear correction')
    header['SMEARTAU'] = (float(tau), 'Smear per pixel crossed, t_read/(N t_exp)')
    header['SMEARSRC'] = (source, 'Where the shutter state came from')
    if raw is not None and 'saturation_threshold' in params:
        trail = saturated_trail_mask(raw, params['saturation_threshold'], axis=axis, direction=direction)
        rows = int(np.any(trail, axis=axis).sum())
        header['SMSATROW'] = (rows, 'Rows with a clipped pixel: smear under-corrected')
    header.add_history(f'Readout smear removed: t_read {smear["readout_time"]} s, '
                       f'tau {tau:.4g}, shutter state from {source}')
    return corrected
