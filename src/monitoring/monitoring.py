#!/usr/bin/env python3
"""
Plot how each SPECULOOS telescope's detector behaves over time: readout noise, overscan level, dark current,
readout mode and gain, CCD temperature and photometric zero point.

Runs nightly from cron on appct (06:30 local, src/jobqueue/cron/crontab.live):
    python /opt/orchard/src/monitoring/monitoring.py
It reads every night of PipelineOutput/v2 and rewrites PipelineOutput/v2/<telescope>/monitoring/:

  monitoring_timeseries[_<camera>].csv/.txt   one row per night that has readoutnoise.dat, overscan.dat and
                                              darkcurrent.dat in its reduction directory
  readoutnoise_*.pdf, zeropoint*.pdf, darkcurrent_*.pdf
                                              zeropoint*.pdf and the zero_point column are the overscan level (ADU)
                                              from overscan.dat; the names predate the photometric zero point
  gain*.pdf, readout_changes*.txt             the gain and readout mode of the night's bias frames, from
                                              reduction/readout.json (calibration/readout.py), and the nights where
                                              the mode changes or disagrees with readout.REGIMES
  ccdtemp*.pdf, photzp*.pdf                   the median CCD-TEMP of the night's *_5_diff.fits frames, and the
                                              median PHOT_ZP of those in the camera's survey filter
  missing_data.txt                            nights with a reduction directory but not all three .dat files

Callisto's files carry _ANDOR or _SPIRIT; each night's camera comes from the measured SPIRIT spells in
jobqueue/config.py.  readout.json exists only for nights reduced since 2026-10-17 (PR #32), so older nights have
no mode or gain until they are rerun.  The _diff.fits values are kept in diff_values_cache.json beside the plots,
so a run reads only files that are new or were rewritten since the last one.

    python monitoring.py [--top-dir DIR] [--out-dir DIR] [--telescopes T ...] [--no-diff]
--out-dir writes <out-dir>/<telescope>/monitoring instead (for trial runs); the default is --top-dir.
"""

import argparse
import glob
import json
import os
import re
import sys
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.dates as mdates  # noqa: E402
import numpy as np  # noqa: E402

try:                        # imported as part of the package
    from calibration import readout
except ImportError:         # run as a bare script: put the source root on the path
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from calibration import readout
from jobqueue.config import DEFAULTS as QUEUE_DEFAULTS, camera  # noqa: E402

# Configuration
telescopes = ['Artemis', 'Callisto', 'Europa', 'Ganymede', 'Io']
top_dir = '/data/SPECULOOSPipeline/PipelineOutput/v2/'
start_date = datetime(2017, 6, 17)
end_date = datetime.now()

# the portal's aperture; CCD-TEMP is the same in every aperture's file, PHOT_ZP is not
diff_aperture = 5
# PHOT_ZP depends on the filter, so only each camera's survey filter is plotted
survey_filter = {'ANDOR': 'I+z', 'SPIRIT': 'zYJ'}
cache_name = 'diff_values_cache.json'
# <target id>_<filter>_<night>_<aperture>_diff.fits; the id may itself contain underscores
diff_name = re.compile(r'_(?P<filter>[^_]+)_(?P<night>\d{8})_(?P<ap>\d+)_diff\.fits$')

csv_columns = ['date', 'readout_noise', 'zero_point', 'dark_current', 'readout_mode', 'expected_mode',
               'modes_seen', 'gain', 'gain_source', 'bias_level', 'ccd_temp', 'phot_zp', 'n_diff']


def get_instrument(telescope, night):
    """'SPIRIT' or 'ANDOR': the camera on the telescope that night (YYYYMMDD).

    Callisto carried SPIRIT on 2022-01-11 and -15, 2022-05-10 to 2023-03-15 and from 2025-02-26, measured from
    every frame's dimensions; the queue's spirit_spells hold those dates.
    """
    return camera(QUEUE_DEFAULTS, telescope, night)


def list_nights(telescope, top, start=None, end=None):
    """The YYYYMMDD night directories under <top>/<telescope>/output, oldest first."""
    start = (start or start_date).strftime('%Y%m%d')
    end = (end or end_date).strftime('%Y%m%d')
    try:
        names = os.listdir(os.path.join(top, telescope, 'output'))
    except OSError:
        return []
    return sorted(n for n in names if readout.valid_night(n) and start <= n <= end)


def read_value_from_file(filepath):
    """Read a float value from file."""
    try:
        with open(filepath, 'r') as f:
            value = float(f.read().strip())
        return value
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def read_readout(reduction_dir, night):
    """The readout mode, gain and bias level of the night's bias frames from readout.json; {} without one.

    The values are the ones T2-T6 used for the master calibrations (readout.calib_values), so gain_source says
    whether the gain was measured that night ('night') or came from the dated table ('table', 'estimate') or the
    instrument config.  bias_level is the median raw bias level in ADU (overscan + bias offset).
    """
    record = readout.read_record(reduction_dir)
    if record is None or record.get('night') != night:
        return {}
    values = readout.calib_values(record)
    entry = record.get('modes', {}).get(values['READMODE']) or {}
    level = entry.get('os_level')
    if level is not None and entry.get('bias_offset') is not None:
        level += entry['bias_offset']
    return {
        'readout_mode': values['READMODE'],
        'expected_mode': record.get('expected_mode'),
        'modes_seen': '+'.join(sorted(record.get('modes', {}))),
        'gain': values['GAIN'],
        'gain_source': values['GAINSRC'],
        'bias_level': level,
        # readout.measure_night's notes that the night disagrees with readout.REGIMES
        'readout_notes': [n for n in record.get('notes', []) if 'table says' in n or 'table expects' in n],
    }


def read_diff_values(path):
    """Median CCD-TEMP of the frames and PHOT_ZP of one *_diff.fits; None for what it does not have."""
    import fitsio
    ccd_temp = phot_zp = None
    with fitsio.FITS(path) as f:
        header = f[0].read_header()
        if 'PHOT_ZP' in header:
            phot_zp = float(header['PHOT_ZP'])
        if 'IMAGELIST' in f:
            imagelist = f['IMAGELIST']
            columns = {c.upper(): c for c in imagelist.get_colnames()}
            if 'CCD-TEMP' in columns:
                temps = np.asarray(imagelist.read_column(columns['CCD-TEMP']), dtype=float)
                temps = temps[np.isfinite(temps)]
                if temps.size:
                    ccd_temp = float(np.median(temps))
    if phot_zp is not None and not np.isfinite(phot_zp):
        phot_zp = None
    return {'ccd_temp': ccd_temp, 'phot_zp': phot_zp}


def cached_diff_values(path, cache, used):
    """read_diff_values(path), from `cache` while the file's size and mtime are unchanged."""
    used.add(path)
    try:
        st = os.stat(path)
    except OSError:
        return {}
    key = [st.st_size, int(st.st_mtime)]
    hit = cache.get(path)
    if hit is not None and hit.get('stat') == key:
        return hit
    try:
        entry = read_diff_values(path)
    except Exception as e:  # a corrupt file must not stop the run; it is retried once it changes
        print(f"Error reading {path}: {e}")
        entry = {'ccd_temp': None, 'phot_zp': None, 'error': str(e)}
    entry['stat'] = key
    cache[path] = entry
    return entry


def night_diff_values(date_dir, instrument, cache, used):
    """The night's median CCD temperature over all its *_5_diff.fits, and median PHOT_ZP over those in the
    camera's survey filter."""
    temps, zps = [], []
    paths = sorted(glob.glob(os.path.join(date_dir, '*', '*_%d_diff.fits' % diff_aperture)))
    for path in paths:
        values = cached_diff_values(path, cache, used)
        if values.get('ccd_temp') is not None:
            temps.append(values['ccd_temp'])
        m = diff_name.search(os.path.basename(path))
        if m and m.group('filter') == survey_filter.get(instrument) and values.get('phot_zp') is not None:
            zps.append(values['phot_zp'])
    return {
        'ccd_temp': float(np.median(temps)) if temps else None,
        'phot_zp': float(np.median(zps)) if zps else None,
        'n_diff': len(paths),
    }


def load_cache(path):
    try:
        with open(path) as f:
            cache = json.load(f)
    except (OSError, ValueError):
        return {}
    return cache if isinstance(cache, dict) else {}


def save_cache(cache, path):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(cache, f, sort_keys=True)
    os.replace(tmp, path)


def collect_monitoring_data(telescope, nights, top, cache=None):
    """One row per night with all three .dat files, and the nights whose reduction directory lacks one.

    With a `cache` dict the night's _diff.fits values are added (and the cache updated; entries for files that
    no longer exist are dropped).
    """
    rows, missing_dates = [], []
    used = set()
    for night in nights:
        date_dir = os.path.join(top, telescope, 'output', night)
        reduction_dir = os.path.join(date_dir, 'reduction')

        # Skip if the night has no reduction directory
        if not os.path.isdir(reduction_dir):
            continue

        files = [os.path.join(reduction_dir, name) for name in ('readoutnoise.dat', 'overscan.dat', 'darkcurrent.dat')]
        if not all(os.path.exists(f) for f in files):
            missing_dates.append(night)
            continue

        ron_value, zp_value, dc_value = (read_value_from_file(f) for f in files)
        if ron_value is None or zp_value is None or dc_value is None:
            continue

        row = {
            'date': datetime.strptime(night, '%Y%m%d'),
            'night': night,
            'instrument': get_instrument(telescope, night),
            'readout_noise': ron_value,
            'zero_point': zp_value,
            'dark_current': dc_value,
        }
        row.update(read_readout(reduction_dir, night))
        if cache is not None:
            row.update(night_diff_values(date_dir, row['instrument'], cache, used))
        rows.append(row)

    if cache is not None:
        for path in [p for p in cache if p not in used]:
            del cache[path]
    return rows, missing_dates


def split_by_instrument(telescope, rows):
    """{file suffix: rows}: Callisto by camera ('_ANDOR', '_SPIRIT'), every other telescope as one series ('')."""
    if telescope != 'Callisto':
        return {'': rows}
    groups = {}
    for row in rows:
        groups.setdefault('_' + row['instrument'], []).append(row)
    return groups


def mode_changes(rows):
    """(night, previous record's night, previous mode, mode) wherever a night's readout mode differs from the
    mode of the last earlier night that had a readout.json."""
    changes, last = [], None
    for row in rows:
        mode = row.get('readout_mode')
        if not mode:
            continue
        if last is not None and mode != last['readout_mode']:
            changes.append((row['night'], last['night'], last['readout_mode'], mode))
        last = row
    return changes


def mark_mode_changes(ax, changes):
    """A dashed line at each readout-mode change, labelled while there are few enough to read."""
    for night, _, old, new in changes:
        x = datetime.strptime(night, '%Y%m%d')
        ax.axvline(x, color='0.4', linestyle='--', linewidth=0.8, alpha=0.7)
        if len(changes) <= 20:
            ax.text(x, 0.98, f' {old}→{new}', transform=ax.get_xaxis_transform(), rotation=90,
                    va='top', ha='right', fontsize=7, color='0.3')


def _title(telescope, instrument_suffix, what):
    if instrument_suffix:
        return f'{telescope} ({instrument_suffix.lstrip("_")}) - {what}'
    return f'{telescope} - {what}'


def plot_series(telescope, dates, values, monitoring_dir, filename, ylabel, title, instrument_suffix='',
                log=False, ylim=None, changes=()):
    """Scatter one quantity against date and save it as a PDF in monitoring_dir."""
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(dates, values, alpha=0.6, s=20)
    if log:
        ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(*ylim)
    mark_mode_changes(ax, changes)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    ax.set_xlabel('Date')
    ax.set_ylabel(ylabel)
    ax.set_title(_title(telescope, instrument_suffix, title))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close(fig)
    print(f"Saved {filename} to {output_path}")


def plot_gain(telescope, rows, monitoring_dir, instrument_suffix='', changes=()):
    """The gain of each night with a readout.json, coloured by readout mode; filled where the night's own flat
    pairs measured it, open where it came from the dated table or the config."""
    fig, ax = plt.subplots(figsize=(12, 6))
    modes = sorted({r['readout_mode'] for r in rows})
    for i, mode in enumerate(modes):
        colour = f'C{i}'
        for measured in (True, False):
            sel = [r for r in rows if r['readout_mode'] == mode and (r['gain_source'] == 'night') == measured]
            if not sel:
                continue
            label = f'{mode} ({"measured that night" if measured else "table or config"})'
            ax.scatter([r['date'] for r in sel], [r['gain'] for r in sel], s=24, label=label,
                       facecolors=colour if measured else 'none', edgecolors=colour)
    mark_mode_changes(ax, changes)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.xticks(rotation=45, ha='right')
    ax.set_xlabel('Date')
    ax.set_ylabel('Gain (e-/ADU)')
    ax.set_title(_title(telescope, instrument_suffix,
                        f'Gain and readout mode of the bias frames ({len(rows)} nights with readout.json)'))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    plt.tight_layout()

    output_path = os.path.join(monitoring_dir, f'gain{instrument_suffix}.pdf')
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close(fig)
    print(f"Saved gain plot to {output_path}")


def _fmt(value, spec=''):
    if value is None:
        return ''
    if isinstance(value, datetime):
        return value.strftime('%Y-%m-%d')
    return format(value, spec) if spec else str(value)


def save_csv(rows, monitoring_dir, instrument_suffix=''):
    """Save monitoring data to CSV and fixed-width text file.  Empty cells are values the night does not have."""
    base_filename = f'monitoring_timeseries{instrument_suffix}'
    csv_path = os.path.join(monitoring_dir, f'{base_filename}.csv')
    txt_path = os.path.join(monitoring_dir, f'{base_filename}.txt')

    with open(csv_path, 'w') as f:
        f.write(','.join(csv_columns) + '\n')
        for row in rows:
            f.write(','.join(_fmt(row.get(c)) for c in csv_columns) + '\n')
    print(f"Saved CSV to {csv_path}")

    widths = [12, 15, 15, 15, 9, 9, 16, 8, 9, 10, 9, 9, 6]
    specs = ['', '.6f', '.6f', '.6f', '', '', '', '.4f', '', '.1f', '.2f', '.4f', 'd']
    headers = ['Date', 'ReadoutNoise', 'ZeroPoint', 'DarkCurrent', 'Mode', 'Expected', 'ModesSeen', 'Gain',
               'GainSrc', 'BiasLevel', 'CCDTemp', 'PhotZP', 'NDiff']
    with open(txt_path, 'w') as f:
        f.write(' '.join(f'{h:<{w}}' for h, w in zip(headers, widths)).rstrip() + '\n')
        f.write(' '.join('=' * w for w in widths) + '\n')
        for row in rows:
            cells = [_fmt(row.get(c), s) for c, s in zip(csv_columns, specs)]
            f.write(' '.join(f'{c:<{w}}' for c, w in zip(cells, widths)).rstrip() + '\n')
    print(f"Saved fixed-width text to {txt_path}")


def save_readout_changes(rows, changes, monitoring_dir, instrument_suffix=''):
    """Readout-mode changes and the nights readout.json says disagree with readout.REGIMES."""
    with_record = [r for r in rows if r.get('readout_mode')]
    path = os.path.join(monitoring_dir, f'readout_changes{instrument_suffix}.txt')
    with open(path, 'w') as f:
        f.write(f'{len(with_record)} of {len(rows)} nights have a readout.json '
                f'(nights calibrated since 2026-10-17).\n\n')
        f.write('Readout-mode changes: the first night with a readout.json in the new mode. The switch happened '
                'after the previous record, which may be long before.\n')
        for night, previous, old, new in changes:
            f.write(f'  {night}  {old} -> {new}  (previous record {previous})\n')
        if not changes:
            f.write('  none\n')
        f.write('\nNights whose mode or gain disagrees with calibration/readout.py REGIMES:\n')
        flagged = [r for r in with_record if r.get('readout_notes')]
        for r in flagged:
            for note in r['readout_notes']:
                f.write(f'  {r["night"]}  {note}\n')
        if not flagged:
            f.write('  none\n')
    print(f"Saved readout changes to {path}")


def save_missing_dates(telescope, missing_dates, monitoring_dir):
    """Save list of dates with missing data files to text file."""
    if len(missing_dates) == 0:
        return

    missing_path = os.path.join(monitoring_dir, 'missing_data.txt')

    with open(missing_path, 'w') as f:
        f.write(
            f'Dates where reduction directory exists but readoutnoise.dat, overscan.dat, or darkcurrent.dat is missing:\n')
        f.write(f'Total: {len(missing_dates)} dates\n\n')
        for date_str in missing_dates:
            f.write(f'{date_str}\n')

    print(f"Saved missing dates list to {missing_path} ({len(missing_dates)} dates)")


def plot_instrument(telescope, rows, monitoring_dir, instrument_suffix=''):
    """Every plot and table for one camera's series."""
    dates = [r['date'] for r in rows]
    spirit = instrument_suffix == '_SPIRIT'
    changes = mode_changes(rows)

    save_csv(rows, monitoring_dir, instrument_suffix)
    ron = [r['readout_noise'] for r in rows]
    plot_series(telescope, dates, ron, monitoring_dir, f'readoutnoise{instrument_suffix}_log.pdf',
                'Readout Noise (electrons) - Log Scale', 'Readout Noise Evolution (Log Scale)',
                instrument_suffix, log=True, changes=changes)
    ylim = (0, 100) if spirit else (0, 20)
    plot_series(telescope, dates, ron, monitoring_dir, f'readoutnoise{instrument_suffix}_linear.pdf',
                'Readout Noise (electrons)',
                f'Readout Noise Evolution (Linear Scale, {ylim[0]}-{ylim[1]} electrons)',
                instrument_suffix, ylim=ylim, changes=changes)
    plot_series(telescope, dates, [r['zero_point'] for r in rows], monitoring_dir,
                f'zeropoint{instrument_suffix}.pdf', 'Overscan level (ADU)', 'Overscan Level Evolution',
                instrument_suffix, changes=changes)
    dc = [r['dark_current'] for r in rows]
    plot_series(telescope, dates, dc, monitoring_dir, f'darkcurrent{instrument_suffix}_log.pdf',
                'Dark Current (electrons/s) - Log Scale', 'Dark Current Evolution (Log Scale)',
                instrument_suffix, log=True)
    ylim = (0, 500) if spirit else (0, 1)
    plot_series(telescope, dates, dc, monitoring_dir, f'darkcurrent{instrument_suffix}_linear.pdf',
                'Dark Current (electrons/s)',
                f'Dark Current Evolution (Linear Scale, {ylim[0]}-{ylim[1]} e⁻/s)',
                instrument_suffix, ylim=ylim)

    with_record = [r for r in rows if r.get('readout_mode')]
    if with_record:
        plot_gain(telescope, with_record, monitoring_dir, instrument_suffix, changes)
        save_readout_changes(rows, changes, monitoring_dir, instrument_suffix)

    temps = [r for r in rows if r.get('ccd_temp') is not None]
    if temps:
        plot_series(telescope, [r['date'] for r in temps], [r['ccd_temp'] for r in temps], monitoring_dir,
                    f'ccdtemp{instrument_suffix}.pdf', 'CCD temperature (°C)',
                    f'CCD Temperature (median of the night\'s *_{diff_aperture}_diff.fits frames)',
                    instrument_suffix)
    zps = [r for r in rows if r.get('phot_zp') is not None]
    if zps:
        filt = survey_filter['SPIRIT' if spirit else 'ANDOR']
        plot_series(telescope, [r['date'] for r in zps], [r['phot_zp'] for r in zps], monitoring_dir,
                    f'photzp{instrument_suffix}.pdf', 'PHOT_ZP (mag)',
                    f'Photometric Zero Point ({filt}, aperture {diff_aperture}, median over the night\'s targets)',
                    instrument_suffix)
    print(f"Created plots for {telescope}{instrument_suffix.replace('_', ' ')} ({len(rows)} data points)")


def plot_telescope_monitoring(telescope, rows, missing_dates, monitoring_dir):
    """Create and save plots and tables for a single telescope."""
    os.makedirs(monitoring_dir, exist_ok=True)
    for instrument_suffix, group in sorted(split_by_instrument(telescope, rows).items()):
        if group:
            plot_instrument(telescope, group, monitoring_dir, instrument_suffix)
    if not rows:
        print(f"No data to plot for {telescope}")
    # Save missing dates file (same for both instruments on Callisto)
    save_missing_dates(telescope, missing_dates, monitoring_dir)


def main(argv=None):
    """Main function to process all telescopes."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--top-dir', default=top_dir, help='the v2 tree to read (default %(default)s)')
    parser.add_argument('--out-dir', default=None, help='write <out-dir>/<telescope>/monitoring (default --top-dir)')
    parser.add_argument('--telescopes', nargs='+', default=telescopes)
    parser.add_argument('--no-diff', action='store_true', help='skip CCD temperature and PHOT_ZP')
    args = parser.parse_args(argv)
    out_dir = args.out_dir or args.top_dir

    print(f"Generating monitoring time series plots")
    print(f"Date range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")
    print(f"Telescopes: {', '.join(args.telescopes)}")
    print()

    for telescope in args.telescopes:
        print(f"Processing {telescope}...")
        monitoring_dir = os.path.join(out_dir, telescope, 'monitoring')
        os.makedirs(monitoring_dir, exist_ok=True)
        cache_path = os.path.join(monitoring_dir, cache_name)
        cache = None if args.no_diff else load_cache(cache_path)
        nights = list_nights(telescope, args.top_dir)
        rows, missing_dates = collect_monitoring_data(telescope, nights, args.top_dir, cache)
        plot_telescope_monitoring(telescope, rows, missing_dates, monitoring_dir)
        if cache is not None:
            save_cache(cache, cache_path)
        print()


if __name__ == '__main__':
    main()
