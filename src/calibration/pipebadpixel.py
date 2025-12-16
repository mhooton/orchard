import os
import sys
import argparse
import numpy as np
from calibration.pipeutils import get_instrument_parameters, open_fits_file, create_dark_dict, create_flat_dict
from astropy.io import fits
from matplotlib import pyplot as pl

def detect_bad_pixels(dark, hot_sigma=5, cold_sigma=5):
    """
    Detect hot and cold pixels as outliers in the dark frame.

    Parameters:
    -----------
    dark : 2D array
        Bias-subtracted master dark (or shortest exposure dark)
    hot_sigma : float
        Number of standard deviations above median to flag hot pixels (default 5)
    cold_sigma : float
        Number of standard deviations below median to flag cold pixels (default 5)

    Returns:
    --------
    hot_pixel_map : 2D boolean array
        True = hot pixel
    cold_pixel_map : 2D boolean array
        True = cold pixel
    """

    # Calculate threshold
    dark_median = np.median(dark)
    dark_std = np.std(dark)

    hot_threshold = dark_median + hot_sigma * dark_std
    cold_threshold = dark_median - cold_sigma * dark_std

    # Flag pixels above/below thresholds
    hot_pixel_map = dark > hot_threshold
    cold_pixel_map = dark < cold_threshold

    return hot_pixel_map, cold_pixel_map

def detect_flatbad_pixels(master_flat, flat_threshold):
    """
    Detect flatbad pixels with low response in the master flat.

    Parameters:
    -----------
    master_flat : 2D array
        Bias-subtracted master flat field
    flat_threshold : float
        Absolute threshold for flat response (e.g., 0.1)

    Returns:
    --------
    flatbad_map : 2D boolean array
        True = flatbad pixel
    """

    # Flag pixels below threshold
    flatbad_map = master_flat < flat_threshold

    return flatbad_map

def bad_pixel_maps(inlist, caldir, outdir, darknames, flatnames, bpmname, biasname):
    """
    Find bad pixels in the dark image based on a given threshold

    Args:
        dark: Matching dark image
    """

    # Get parameters once from the first file
    with open(inlist) as infile:
        first_filename = infile.readline().strip()

    with open_fits_file(first_filename) as hdul:
        params = get_instrument_parameters(hdul)
        # exposure = hdul[0].header['exptime']
        # filter = hdul[0].header['filter']
        # filter = filter.replace("'", "")




    if 'bad_pixel_correction' in params:
        if params['bad_pixel_correction']:

            print("Creating bad pixel map")

            bpm_config = params['bad_pixel_correction']

            try:
                with open(caldir + "readoutnoise.dat", "r") as f:
                    ron = f.read()
                    print('Readout noise is ' + str(ron) + ' ADU')
            except Exception as e:
                print(f"Error finding or reading readout noise: {e}. Setting to 0.")
                ron = 0.

            dark_dict = create_dark_dict(darknames)
            flat_dict = create_flat_dict(flatnames)

            if bpm_config['preferred_flat_filter'] in flat_dict:
                filter = bpm_config['preferred_flat_filter']
            else:
                filter = list(flat_dict.keys())[0]
                print(f"Master flat with preferred filter of {filter} not found for this date. Using {filter} instead.")

            flat = flat_dict[filter]

            # dark_threshold_sigma = bpm_config.get('dark_threshold_sigma', 5)
            hot_sigma_threshold = bpm_config.get('hot_sigma_threshold', 5)
            cold_sigma_threshold = bpm_config.get('cold_sigma_threshold', 5)
            flat_threshold = bpm_config.get('flat_threshold', 0.1)
            # Read overwrite setting from config
            overwrite_existing = bpm_config.get('overwrite_existing', True)

            with open_fits_file(caldir + biasname) as hdul:
                master_bias = hdul[0].data
            effective_saturation_map = params['saturation_threshold'] - master_bias

            # Select appropriate dark (shortest exposure or combined)
            numeric_keys = [exp for exp in dark_dict.keys() if isinstance(exp, (int, float))]
            if len(numeric_keys) > 0:
                min_exposure = min(numeric_keys)
                dark_for_detection = dark_dict[min_exposure]
            elif 'combined' in dark_dict:
                dark_for_detection = dark_dict['combined']
            else:
                print("ERROR: No suitable dark found for bad pixel detection")
                return

            # Hot and cold pixels
            hot_filename = caldir + "hot_pixel_map.fits"
            cold_filename = caldir + "cold_pixel_map.fits"

            if not overwrite_existing and os.path.exists(hot_filename) and os.path.exists(cold_filename):
                print(f"Loading existing hot and cold pixel maps")
                hot_pixel_map = fits.getdata(hot_filename).astype(bool)
                cold_pixel_map = fits.getdata(cold_filename).astype(bool)
            else:
                print("Creating hot and cold pixel maps")
                hot_pixel_map, cold_pixel_map = detect_bad_pixels(dark_for_detection,
                                                                  hot_sigma=hot_sigma_threshold,
                                                                  cold_sigma=cold_sigma_threshold)
                fits.writeto(hot_filename, hot_pixel_map.astype(int), overwrite=True)
                fits.writeto(cold_filename, cold_pixel_map.astype(int), overwrite=True)

            # Flatbad pixels
            flatbad_filename = caldir + "flatbad_pixel_map.fits"
            if not overwrite_existing and os.path.exists(flatbad_filename):
                print(f"Loading existing flatbad pixel map from {flatbad_filename}")
                flatbad_pixel_map = fits.getdata(flatbad_filename).astype(bool)
            else:
                print("Creating flatbad pixel map")
                flatbad_pixel_map = detect_flatbad_pixels(flat, flat_threshold)
                fits.writeto(flatbad_filename, flatbad_pixel_map.astype(int), overwrite=True)

            # Combine all maps
            bad_pixel_map = hot_pixel_map | cold_pixel_map | flatbad_pixel_map

            # Calculate total pixels
            total_pixels = bad_pixel_map.size

            # Print statistics for each type
            print("\n" + "=" * 60)
            print("BAD PIXEL STATISTICS")
            print("=" * 60)

            maps = {
                'Hot': hot_pixel_map,
                'Cold': cold_pixel_map,
                'Flatbad': flatbad_pixel_map,
                'Total': bad_pixel_map
            }

            for name, pixel_map in maps.items():
                n_bad = np.sum(pixel_map)
                pct_bad = 100 * n_bad / total_pixels
                print(f"{name:12s}: {n_bad:6d} pixels ({pct_bad:.4f}%)")

            print("\n" + "-" * 60)
            print("UNIQUE BAD PIXELS (not flagged by other maps)")
            print("-" * 60)

            # Calculate unique pixels for each type (excluding Total)
            unique_maps = {
                'Hot': hot_pixel_map & ~(cold_pixel_map | flatbad_pixel_map),
                'Cold': cold_pixel_map & ~(hot_pixel_map | flatbad_pixel_map),
                'Flatbad': flatbad_pixel_map & ~(hot_pixel_map | cold_pixel_map)
            }

            for name, unique_map in unique_maps.items():
                n_unique = np.sum(unique_map)
                n_total_for_type = np.sum(maps[name])
                pct_of_total_pixels = 100 * n_unique / total_pixels
                pct_of_category = 100 * n_unique / n_total_for_type if n_total_for_type > 0 else 0
                print(
                    f"{name:12s}: {n_unique:6d} pixels ({pct_of_total_pixels:.4f}% of all, {pct_of_category:.2f}% of {name.lower()})")

            print("=" * 60 + "\n")

            # Save final bad pixel map
            fits.writeto(caldir + bpmname, bad_pixel_map.astype(int), overwrite=True)

        else:
            print("No bad pixel correction requested. Skipping bad pixel correction.")
    else:
        print("No bad pixel correction requested. Skipping bad pixel correction.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inlist')
    parser.add_argument('-c', '--caldir')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-d', '--darknames', nargs='+')
    parser.add_argument('-f', '--flatnames', nargs='+')
    parser.add_argument('-bpm', '--bpmname')
    parser.add_argument('-b', '--biasname')
    args = parser.parse_args()

    inlist = args.inlist
    caldir = args.caldir + '/'
    outdir = args.outdir + '/'
    darknames = [caldir + d for d in args.darknames]
    flatnames = [caldir + f for f in args.flatnames]
    bpmname = args.bpmname
    biasname = args.biasname

    bad_pixel_maps(inlist, caldir, outdir, darknames, flatnames, bpmname, biasname)

if __name__ == '__main__':
    main()

