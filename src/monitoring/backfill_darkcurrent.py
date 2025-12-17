#!/usr/bin/env python3
"""
Backfill darkcurrent.dat files by extracting DARKCUR from science frame headers.

This script reads the DARKCUR header keyword from processed science images and
writes it to darkcurrent.dat in the reduction directory for each night.
"""

import os
import glob
from datetime import datetime, timedelta
from astropy.io import fits

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


def backfill_darkcurrent(telescope, date_str, top_dir):
    """
    Extract DARKCUR from a science frame and write to darkcurrent.dat.

    Args:
        telescope: Telescope name
        date_str: Date in YYYYMMDD format
        top_dir: Top level pipeline output directory

    Returns:
        str: Status message ('exists', 'created', 'no_frames', 'no_darkcur', 'error')
    """
    date_dir = f"{top_dir}/{telescope}/output/{date_str}"
    reduction_dir = f"{date_dir}/reduction"
    darkcurrent_file = f"{reduction_dir}/darkcurrent.dat"

    # Check if darkcurrent.dat already exists
    if os.path.exists(darkcurrent_file):
        return 'exists'

    # Check if date directory exists
    if not os.path.exists(date_dir):
        return 'no_date_dir'

    # Find target directories (not "reduction")
    try:
        all_items = os.listdir(date_dir)
        target_dirs = [d for d in all_items
                       if os.path.isdir(f"{date_dir}/{d}") and d != "reduction"]
    except Exception as e:
        return 'error'

    if not target_dirs:
        return 'no_targets'

    # Try each target until we find a science frame with DARKCUR
    for target in target_dirs:
        science_dir = f"{date_dir}/{target}/1"
        if not os.path.exists(science_dir):
            continue

        # Find first .fits file
        fits_files = glob.glob(f"{science_dir}/*.fits")
        if not fits_files:
            continue

        # Read DARKCUR from header
        try:
            with fits.open(fits_files[0]) as hdul:
                if 'DARKCUR' in hdul[0].header:
                    darkcur = hdul[0].header['DARKCUR']

                    # Create reduction directory if it doesn't exist
                    os.makedirs(reduction_dir, exist_ok=True)

                    # Write to darkcurrent.dat
                    with open(darkcurrent_file, 'w') as f:
                        f.write(str(darkcur))

                    return 'created'
        except Exception as e:
            print(f"Error reading {fits_files[0]}: {e}")
            continue

    # If we get here, we found science frames but no DARKCUR header
    return 'no_darkcur'


def main():
    """Main function to backfill all telescopes."""
    print(f"Backfilling darkcurrent.dat files")
    print(f"Date range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")
    print(f"Telescopes: {', '.join(telescopes)}")
    print()

    # Generate date range
    dates = get_date_range(start_date, end_date)
    print(f"Scanning {len(dates)} dates...\n")

    # Statistics
    stats = {
        'exists': 0,
        'created': 0,
        'no_date_dir': 0,
        'no_targets': 0,
        'no_darkcur': 0,
        'error': 0
    }

    # Process each telescope
    for telescope in telescopes:
        print(f"Processing {telescope}...")
        telescope_stats = {key: 0 for key in stats.keys()}

        for date_str in dates:
            result = backfill_darkcurrent(telescope, date_str, top_dir)
            telescope_stats[result] += 1

            if result == 'created':
                print(f"  Created: {date_str}")
            elif result == 'no_darkcur':
                print(f"  WARNING: {date_str} - Science frames exist but no DARKCUR header found")

        # Update global stats
        for key in stats.keys():
            stats[key] += telescope_stats[key]

        # Print telescope summary
        print(f"\n  {telescope} Summary:")
        print(f"    Already existed:     {telescope_stats['exists']}")
        print(f"    Newly created:       {telescope_stats['created']}")
        print(f"    No date directory:   {telescope_stats['no_date_dir']}")
        print(f"    No target dirs:      {telescope_stats['no_targets']}")
        print(f"    Missing DARKCUR:     {telescope_stats['no_darkcur']}")
        print(f"    Errors:              {telescope_stats['error']}")
        print()

    # Print overall summary
    print("=" * 60)
    print("OVERALL SUMMARY")
    print("=" * 60)
    print(f"Total dates scanned:     {len(dates) * len(telescopes)}")
    print(f"Already existed:         {stats['exists']}")
    print(f"Newly created:           {stats['created']}")
    print(f"No date directory:       {stats['no_date_dir']}")
    print(f"No target directories:   {stats['no_targets']}")
    print(f"Missing DARKCUR header:  {stats['no_darkcur']}")
    print(f"Errors:                  {stats['error']}")
    print("=" * 60)


if __name__ == '__main__':
    main()