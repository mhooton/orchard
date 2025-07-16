#!/usr/bin/env python3
"""
SPECULOOS telescope data download script
Downloads observation data from remote server and applies filtering based on flags
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from astropy.io import fits
import glob

# Configuration
REMOTE_HOST = "speculoos@appcs.ra.phy.cam.ac.uk"
REMOTE_TOPDIR = "/appct/data/SPECULOOSPipeline"
LOCAL_TOPDIR = "/Users/matthewhooton/data/SPECULOOSPipeline"


def run_rsync(source, destination, create_dirs=True):
    """Run rsync command to sync directories"""
    if create_dirs:
        os.makedirs(destination, exist_ok=True)

    cmd = [
        "rsync", "-avz", "--progress",
        f"{REMOTE_HOST}:{source}/",
        f"{destination}/"
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Warning: rsync failed for {source}")
        print(f"Error: {result.stderr}")
        return False
    return True


def format_date_for_path(date_str, format_type):
    """Convert YYYYMMDD to appropriate format"""
    if format_type == "dash":
        # YYYYMMDD -> YYYY-MM-DD
        return f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
    else:
        # Keep YYYYMMDD format
        return date_str


def download_basic_data(telescope_name, date):
    """Download plans, technical logs, and fits images"""
    date_dash = format_date_for_path(date, "dash")
    date_plain = format_date_for_path(date, "plain")

    # Define remote and local paths
    paths = [
        {
            "name": "plans",
            "remote": f"{REMOTE_TOPDIR}/Observations/{telescope_name}/schedule/Plans_by_date/{date_dash}",
            "local": f"{LOCAL_TOPDIR}/Observations/{telescope_name}/schedule/Plans_by_date/{date_dash}"
        },
        {
            "name": "technical_logs",
            "remote": f"{REMOTE_TOPDIR}/Observations/{telescope_name}/technical_logs/{date_dash}",
            "local": f"{LOCAL_TOPDIR}/Observations/{telescope_name}/technical_logs/{date_dash}"
        },
        {
            "name": "fits_images",
            "remote": f"{REMOTE_TOPDIR}/Observations/{telescope_name}/images/{date_plain}",
            "local": f"{LOCAL_TOPDIR}/Observations/{telescope_name}/images/{date_plain}"
        }
    ]

    # Download each directory
    for path_info in paths:
        print(f"\nDownloading {path_info['name']}...")
        success = run_rsync(path_info["remote"], path_info["local"])
        if not success:
            print(f"Failed to download {path_info['name']}")

    return paths[2]["local"]  # Return fits images directory


def read_fits_header(fits_file):
    """Read FITS header and return IMAGETYP and OBJECT values"""
    try:
        with fits.open(fits_file) as hdul:
            header = hdul[0].header
            imagetyp = header.get('IMAGETYP', '').strip()
            object_name = header.get('OBJECT', '').strip()
            return imagetyp, object_name
    except Exception as e:
        # More detailed error reporting
        print(f"Warning: Could not read FITS header from {os.path.basename(fits_file)}: {type(e).__name__}: {e}")
        return None, None


def filter_fits_images(fits_dir, target_name, verbose=False):
    """Remove FITS files that don't match the target"""
    if not os.path.exists(fits_dir):
        print(f"FITS directory does not exist: {fits_dir}")
        return []

    # Get all .fits files (case insensitive)
    fits_files = []
    for pattern in ["*.fits", "*.FITS", "*.fit", "*.FIT"]:
        fits_files.extend(glob.glob(os.path.join(fits_dir, pattern)))

    print(f"Found {len(fits_files)} FITS files to process")

    deleted_files = []
    processed_count = 0

    for fits_file in fits_files:
        processed_count += 1
        if verbose or processed_count % 100 == 0:
            print(f"Processing file {processed_count}/{len(fits_files)}: {os.path.basename(fits_file)}")

        imagetyp, object_name = read_fits_header(fits_file)

        if imagetyp and object_name:
            if verbose:
                print(f"  IMAGETYP={imagetyp}, OBJECT={object_name}")
            if imagetyp.lower() == "light frame" and object_name != target_name:
                if verbose:
                    print(f"  -> Will delete (OBJECT={object_name}, not {target_name})")
                print(f"Deleting {os.path.basename(fits_file)} (OBJECT={object_name}, not {target_name})")
                os.remove(fits_file)
                deleted_files.append(fits_file)
            elif verbose:
                print(f"  -> Keeping (IMAGETYP={imagetyp}, OBJECT={object_name})")
        else:
            if verbose:
                print(f"  -> Skipping (could not read headers)")

    return deleted_files


def get_unique_targets(fits_dir, verbose=False):
    """Get unique target names from Light frame FITS files"""
    if not os.path.exists(fits_dir):
        print(f"FITS directory does not exist: {fits_dir}")
        return []

    # Get all .fits files (case insensitive)
    fits_files = []
    for pattern in ["*.fits", "*.FITS", "*.fit", "*.FIT"]:
        fits_files.extend(glob.glob(os.path.join(fits_dir, pattern)))

    print(f"Processing {len(fits_files)} FITS files to find unique targets")

    targets = set()
    processed_count = 0

    for fits_file in fits_files:
        processed_count += 1
        if verbose or processed_count % 100 == 0:
            print(f"Processing file {processed_count}/{len(fits_files)}")

        imagetyp, object_name = read_fits_header(fits_file)

        if imagetyp and imagetyp.lower() == "light frame" and object_name:
            targets.add(object_name)
            if verbose:
                print(f"  Found target: {object_name}")

    return list(targets)


def download_preexisting_cat(telescope_name, target_names):
    """Download preexisting catalog files for given targets"""
    remote_stackimages = f"{REMOTE_TOPDIR}/PipelineOutput/v2/{telescope_name}/output/StackImages"
    local_stackimages = f"{LOCAL_TOPDIR}/PipelineOutput/v2/{telescope_name}/output/StackImages"

    # Create local directory
    os.makedirs(local_stackimages, exist_ok=True)

    for target in target_names:
        target_upper = target.upper()
        print(f"\nDownloading preexisting catalog files for target: {target} (searching for {target_upper}*)")

        # Use rsync with include pattern to download only matching files
        cmd = [
            "rsync", "-avz", "--progress",
            "--include", f"{target_upper}*",
            "--exclude", "*",
            f"{REMOTE_HOST}:{remote_stackimages}/",
            f"{local_stackimages}/"
        ]

        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Warning: rsync failed for StackImages target {target}")
            print(f"Error: {result.stderr}")
        else:
            print(f"Successfully downloaded files for target {target}")

    # List what was actually downloaded
    downloaded_files = glob.glob(os.path.join(local_stackimages, "*"))
    if downloaded_files:
        print(f"\nDownloaded {len(downloaded_files)} files to StackImages:")
        for f in downloaded_files[:10]:  # Show first 10
            print(f"  {os.path.basename(f)}")
        if len(downloaded_files) > 10:
            print(f"  ... and {len(downloaded_files) - 10} more files")
    else:
        print("No files were downloaded to StackImages")


def main():
    parser = argparse.ArgumentParser(description="Download SPECULOOS telescope observation data")
    parser.add_argument("telescope_name", help="Name of the telescope")
    parser.add_argument("date", help="Date in YYYYMMDD format")
    parser.add_argument("--only-target", help="Keep only FITS images for this target")
    parser.add_argument("--preexisting-cat", action="store_true",
                        help="Download preexisting catalog files")
    parser.add_argument("--no-download", action="store_true",
                        help="Skip downloads, only perform filtering on existing files")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose output")

    args = parser.parse_args()

    # Validate date format
    if len(args.date) != 8 or not args.date.isdigit():
        print("Error: Date must be in YYYYMMDD format")
        sys.exit(1)

    print(f"Processing data for {args.telescope_name} on {args.date}")
    if args.no_download:
        print("No-download mode: skipping all downloads")

    # Download basic data (unless no-download flag is set)
    if not args.no_download:
        fits_dir = download_basic_data(args.telescope_name, args.date)
    else:
        # Construct the fits directory path without downloading
        date_plain = format_date_for_path(args.date, "plain")
        fits_dir = f"{LOCAL_TOPDIR}/Observations/{args.telescope_name}/images/{date_plain}"

    # Construct the fits directory path for filtering (even when commented out)
    date_plain = format_date_for_path(args.date, "plain")
    fits_dir = f"{LOCAL_TOPDIR}/Observations/{args.telescope_name}/images/{date_plain}"

    # Apply only-target filter if specified
    if args.only_target:
        print(f"\nFiltering FITS images for target: {args.only_target}")
        deleted = filter_fits_images(fits_dir, args.only_target, args.verbose)
        print(f"Deleted {len(deleted)} files")

    # Handle preexisting-cat flag
    if args.preexisting_cat and not args.no_download:
        if args.only_target:
            # Use specified target
            target_names = [args.only_target]
        else:
            # Get all unique targets from FITS files
            print("\nExtracting unique targets from FITS files...")
            target_names = get_unique_targets(fits_dir, args.verbose)
            print(f"Found targets: {target_names}")

        if target_names:
            download_preexisting_cat(args.telescope_name, target_names)
        else:
            print("No targets found for preexisting catalog download")
    elif args.preexisting_cat and args.no_download:
        print("Skipping preexisting catalog download due to --no-download flag")

    if not args.no_download:
        print("\nDownload complete!")
    else:
        print("\nProcessing complete!")


if __name__ == "__main__":
    main()