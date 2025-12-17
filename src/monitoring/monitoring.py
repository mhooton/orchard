#!/usr/bin/env python3
"""
Plot readout noise, zero point, and dark current evolution over time for SPECULOOS telescopes.
"""

import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Configuration
telescopes = ['Artemis', 'Callisto', 'Europa', 'Ganymede', 'Io']
top_dir = '/data/SPECULOOSPipeline/PipelineOutput/v2/'
start_date = datetime(2017, 6, 17)
end_date = datetime.now()


def get_date_range(start, end):
    """Generate all dates between start and end."""
    dates = []
    current = start
    while current <= end:
        dates.append(current.strftime('%Y%m%d'))
        current += timedelta(days=1)
    return dates


def get_callisto_instrument(date):
    """
    Determine which instrument was used on Callisto for a given date.
    ANDOR: before 01/06/2022 and between 17/03/2023 and 10/03/2025
    SPIRIT: outside of those times
    """
    # ANDOR periods
    andor_period1_end = datetime(2022, 6, 1)
    andor_period2_start = datetime(2023, 3, 17)
    andor_period2_end = datetime(2025, 3, 10)

    if date < andor_period1_end:
        return 'ANDOR'
    elif andor_period2_start <= date <= andor_period2_end:
        return 'ANDOR'
    else:
        return 'SPIRIT'


def read_value_from_file(filepath):
    """Read a float value from file."""
    try:
        with open(filepath, 'r') as f:
            value = float(f.read().strip())
        return value
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def collect_monitoring_data(telescope, dates):
    """Collect readout noise, zero point, and dark current data for a telescope across all dates."""
    data_dates = []
    readoutnoise_values = []
    zeropoint_values = []
    darkcurrent_values = []
    missing_dates = []

    # For Callisto, separate by instrument
    if telescope == 'Callisto':
        andor_dates = []
        andor_readoutnoise = []
        andor_zeropoint = []
        andor_darkcurrent = []
        spirit_dates = []
        spirit_readoutnoise = []
        spirit_zeropoint = []
        spirit_darkcurrent = []

    for date_str in dates:
        reduction_dir = os.path.join(top_dir, telescope, 'output', date_str, 'reduction')

        # Skip if date directory doesn't exist
        if not os.path.exists(reduction_dir):
            continue

        # Check for required files
        readoutnoise_file = os.path.join(reduction_dir, 'readoutnoise.dat')
        overscan_file = os.path.join(reduction_dir, 'overscan.dat')
        darkcurrent_file = os.path.join(reduction_dir, 'darkcurrent.dat')

        if not os.path.exists(readoutnoise_file) or not os.path.exists(overscan_file) or not os.path.exists(
                darkcurrent_file):
            missing_dates.append(date_str)
            continue

        # Read the values
        ron_value = read_value_from_file(readoutnoise_file)
        zp_value = read_value_from_file(overscan_file)
        dc_value = read_value_from_file(darkcurrent_file)

        if ron_value is not None and zp_value is not None and dc_value is not None:
            date_obj = datetime.strptime(date_str, '%Y%m%d')

            # For Callisto, sort by instrument
            if telescope == 'Callisto':
                instrument = get_callisto_instrument(date_obj)
                if instrument == 'ANDOR':
                    andor_dates.append(date_obj)
                    andor_readoutnoise.append(ron_value)
                    andor_zeropoint.append(zp_value)
                    andor_darkcurrent.append(dc_value)
                else:
                    spirit_dates.append(date_obj)
                    spirit_readoutnoise.append(ron_value)
                    spirit_zeropoint.append(zp_value)
                    spirit_darkcurrent.append(dc_value)
            else:
                data_dates.append(date_obj)
                readoutnoise_values.append(ron_value)
                zeropoint_values.append(zp_value)
                darkcurrent_values.append(dc_value)

    # Return instrument-separated data for Callisto
    if telescope == 'Callisto':
        return {
            'ANDOR': (andor_dates, andor_readoutnoise, andor_zeropoint, andor_darkcurrent),
            'SPIRIT': (spirit_dates, spirit_readoutnoise, spirit_zeropoint, spirit_darkcurrent),
            'missing_dates': missing_dates
        }
    else:
        return data_dates, readoutnoise_values, zeropoint_values, darkcurrent_values, missing_dates


def plot_readoutnoise_log(telescope, dates, values, monitoring_dir, instrument_suffix=''):
    """Create and save a log scale plot for readout noise."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot scatter points
    ax.scatter(dates, values, alpha=0.6, s=20)

    # Set log scale
    ax.set_yscale('log')

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    # Labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('Readout Noise (electrons) - Log Scale')
    title = f'{telescope} - Readout Noise Evolution (Log Scale)'
    if instrument_suffix:
        title = f'{telescope} ({instrument_suffix}) - Readout Noise Evolution (Log Scale)'
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    filename = f'readoutnoise{instrument_suffix}_log.pdf'
    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close()

    print(f"Saved readout noise log scale plot to {output_path}")


def plot_readoutnoise_linear(telescope, dates, values, monitoring_dir, instrument_suffix=''):
    """Create and save a linear scale plot with ylim [0, 20] for readout noise."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot scatter points
    ax.scatter(dates, values, alpha=0.6, s=20)

    # Set y limits
    ax.set_ylim(0, 20)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    # Labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('Readout Noise (electrons)')
    title = f'{telescope} - Readout Noise Evolution (Linear Scale, 0-20 electrons)'
    if instrument_suffix:
        title = f'{telescope} ({instrument_suffix}) - Readout Noise Evolution (Linear Scale, 0-20 electrons)'
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    filename = f'readoutnoise{instrument_suffix}_linear.pdf'
    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close()

    print(f"Saved readout noise linear scale plot to {output_path}")


def plot_zeropoint_linear(telescope, dates, values, monitoring_dir, instrument_suffix=''):
    """Create and save a linear scale plot for zero point."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot scatter points
    ax.scatter(dates, values, alpha=0.6, s=20)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    # Labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('Zero Point (ADU)')
    title = f'{telescope} - Zero Point Evolution'
    if instrument_suffix:
        title = f'{telescope} ({instrument_suffix}) - Zero Point Evolution'
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    filename = f'zeropoint{instrument_suffix}.pdf'
    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close()

    print(f"Saved zero point plot to {output_path}")


def plot_darkcurrent_log(telescope, dates, values, monitoring_dir, instrument_suffix=''):
    """Create and save a log scale plot for dark current."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot scatter points
    ax.scatter(dates, values, alpha=0.6, s=20)

    # Set log scale
    ax.set_yscale('log')

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    # Labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('Dark Current (electrons/s) - Log Scale')
    title = f'{telescope} - Dark Current Evolution (Log Scale)'
    if instrument_suffix:
        title = f'{telescope} ({instrument_suffix}) - Dark Current Evolution (Log Scale)'
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    filename = f'darkcurrent{instrument_suffix}_log.pdf'
    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close()

    print(f"Saved dark current log scale plot to {output_path}")


def plot_darkcurrent_linear(telescope, dates, values, monitoring_dir, instrument_suffix=''):
    """Create and save a linear scale plot for dark current with instrument-specific y-limits."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot scatter points
    ax.scatter(dates, values, alpha=0.6, s=20)

    # Set y limits based on instrument
    if instrument_suffix == '_SPIRIT':
        ax.set_ylim(0, 100)
        ylim_text = '0-100'
    else:
        # ANDOR (default for all telescopes except Callisto SPIRIT)
        ax.set_ylim(0, 1)
        ylim_text = '0-1'

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.xticks(rotation=45, ha='right')

    # Labels and title
    ax.set_xlabel('Date')
    ax.set_ylabel('Dark Current (electrons/s)')
    title = f'{telescope} - Dark Current Evolution (Linear Scale, {ylim_text} e⁻/s)'
    if instrument_suffix:
        title = f'{telescope} ({instrument_suffix}) - Dark Current Evolution (Linear Scale, {ylim_text} e⁻/s)'
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    filename = f'darkcurrent{instrument_suffix}_linear.pdf'
    output_path = os.path.join(monitoring_dir, filename)
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.close()

    print(f"Saved dark current linear scale plot to {output_path}")


def save_csv(telescope, dates, readoutnoise_values, zeropoint_values, darkcurrent_values, monitoring_dir,
             instrument_suffix=''):
    """Save monitoring data to CSV and fixed-width text file."""
    base_filename = f'monitoring_timeseries{instrument_suffix}'
    csv_path = os.path.join(monitoring_dir, f'{base_filename}.csv')
    txt_path = os.path.join(monitoring_dir, f'{base_filename}.txt')

    # Save CSV
    with open(csv_path, 'w') as f:
        f.write('date,readout_noise,zero_point,dark_current\n')
        for date, ron, zp, dc in zip(dates, readoutnoise_values, zeropoint_values, darkcurrent_values):
            f.write(f'{date.strftime("%Y-%m-%d")},{ron},{zp},{dc}\n')

    print(f"Saved CSV to {csv_path}")

    # Save fixed-width text file
    with open(txt_path, 'w') as f:
        # Write header
        f.write(f"{'Date':<12} {'ReadoutNoise':<15} {'ZeroPoint':<15} {'DarkCurrent':<15}\n")
        f.write(f"{'=' * 12} {'=' * 15} {'=' * 15} {'=' * 15}\n")

        # Write data rows
        for date, ron, zp, dc in zip(dates, readoutnoise_values, zeropoint_values, darkcurrent_values):
            f.write(f"{date.strftime('%Y-%m-%d'):<12} {ron:<15.6f} {zp:<15.6f} {dc:<15.6f}\n")

    print(f"Saved fixed-width text to {txt_path}")


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


def plot_telescope_monitoring(telescope, data, missing_dates):
    """Create and save plots and CSV for a single telescope."""
    # Create monitoring directory if it doesn't exist
    monitoring_dir = os.path.join(top_dir, telescope, 'monitoring')
    os.makedirs(monitoring_dir, exist_ok=True)

    # Handle Callisto separately (instrument-specific data)
    if telescope == 'Callisto':
        andor_dates, andor_readoutnoise, andor_zeropoint, andor_darkcurrent = data['ANDOR']
        spirit_dates, spirit_readoutnoise, spirit_zeropoint, spirit_darkcurrent = data['SPIRIT']

        if len(andor_dates) > 0:
            save_csv(telescope, andor_dates, andor_readoutnoise, andor_zeropoint, andor_darkcurrent, monitoring_dir,
                     '_ANDOR')
            plot_readoutnoise_log(telescope, andor_dates, andor_readoutnoise, monitoring_dir, '_ANDOR')
            plot_readoutnoise_linear(telescope, andor_dates, andor_readoutnoise, monitoring_dir, '_ANDOR')
            plot_zeropoint_linear(telescope, andor_dates, andor_zeropoint, monitoring_dir, '_ANDOR')
            plot_darkcurrent_log(telescope, andor_dates, andor_darkcurrent, monitoring_dir, '_ANDOR')
            plot_darkcurrent_linear(telescope, andor_dates, andor_darkcurrent, monitoring_dir, '_ANDOR')
            print(f"Created ANDOR plots for {telescope} ({len(andor_dates)} data points)")
        else:
            print(f"No ANDOR data to plot for {telescope}")

        if len(spirit_dates) > 0:
            save_csv(telescope, spirit_dates, spirit_readoutnoise, spirit_zeropoint, spirit_darkcurrent, monitoring_dir,
                     '_SPIRIT')
            plot_readoutnoise_log(telescope, spirit_dates, spirit_readoutnoise, monitoring_dir, '_SPIRIT')
            plot_readoutnoise_linear(telescope, spirit_dates, spirit_readoutnoise, monitoring_dir, '_SPIRIT')
            plot_zeropoint_linear(telescope, spirit_dates, spirit_zeropoint, monitoring_dir, '_SPIRIT')
            plot_darkcurrent_log(telescope, spirit_dates, spirit_darkcurrent, monitoring_dir, '_SPIRIT')
            plot_darkcurrent_linear(telescope, spirit_dates, spirit_darkcurrent, monitoring_dir, '_SPIRIT')
            print(f"Created SPIRIT plots for {telescope} ({len(spirit_dates)} data points)")
        else:
            print(f"No SPIRIT data to plot for {telescope}")

    else:
        # Standard telescopes (all use ANDOR)
        dates, readoutnoise_values, zeropoint_values, darkcurrent_values = data
        if len(dates) == 0:
            print(f"No data to plot for {telescope}")
            return

        save_csv(telescope, dates, readoutnoise_values, zeropoint_values, darkcurrent_values, monitoring_dir)
        plot_readoutnoise_log(telescope, dates, readoutnoise_values, monitoring_dir)
        plot_readoutnoise_linear(telescope, dates, readoutnoise_values, monitoring_dir)
        plot_zeropoint_linear(telescope, dates, zeropoint_values, monitoring_dir)
        plot_darkcurrent_log(telescope, dates, darkcurrent_values, monitoring_dir)
        plot_darkcurrent_linear(telescope, dates, darkcurrent_values, monitoring_dir)
        print(f"Created plots for {telescope} ({len(dates)} data points)")

    # Save missing dates file (same for both instruments on Callisto)
    save_missing_dates(telescope, missing_dates, monitoring_dir)


def main():
    """Main function to process all telescopes."""
    print(f"Generating monitoring time series plots")
    print(f"Date range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")
    print(f"Telescopes: {', '.join(telescopes)}")
    print()

    # Generate date range
    dates = get_date_range(start_date, end_date)
    print(f"Scanning {len(dates)} dates...\n")

    # Process each telescope
    for telescope in telescopes:
        print(f"Processing {telescope}...")
        result = collect_monitoring_data(telescope, dates)

        if telescope == 'Callisto':
            # Callisto returns dict with instrument data
            missing_dates = result['missing_dates']
            plot_telescope_monitoring(telescope, result, missing_dates)
        else:
            # Standard telescopes return tuple
            data_dates, readoutnoise_values, zeropoint_values, darkcurrent_values, missing_dates = result
            plot_telescope_monitoring(telescope,
                                      (data_dates, readoutnoise_values, zeropoint_values, darkcurrent_values),
                                      missing_dates)

        print()


if __name__ == '__main__':
    main()