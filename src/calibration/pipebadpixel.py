import os
import sys
import argparse
import numpy as np
from calibration.pipeutils import get_instrument_parameters, open_fits_file, create_dark_dict, create_flat_dict
from astropy.io import fits
from matplotlib import pyplot as plt


def detect_stuck_pixels(dark_dict, effective_saturation_map, threshold=1.0, max_exposure=32):
    """
    Detect stuck/inverse pixels by fitting dark current vs exposure time.

    Parameters:
    -----------
    dark_dict : dict
        Keys are exposure times (seconds), values are 2D bias-subtracted submaster dark arrays
    effective_saturation_map : 2D array
        Per-pixel saturation threshold (saturation_threshold - master_bias)
    threshold : float
        Dark current threshold in ADU/s (default 1.0)
    max_exposure : float
        Maximum exposure time to include in fitting (default 32s)

    Returns:
    --------
    stuck_map : 2D boolean array
        True = stuck pixel
    """

    # Filter dark_dict for exposure times < max_exposure
    filtered_dict = {exp: dark for exp, dark in dark_dict.items()
                     if isinstance(exp, (int, float)) and exp < max_exposure}

    # Check if we have enough exposure times
    if len(filtered_dict) < 2:
        print(
            f"WARNING: Need at least 2 exposure times < {max_exposure}s for stuck pixel detection, found {len(filtered_dict)}")
        # Get shape from any available dark, or from effective_saturation_map
        if len(filtered_dict) > 0:
            shape = list(filtered_dict.values())[0].shape
        else:
            shape = effective_saturation_map.shape
        return np.zeros(shape, dtype=bool)

    # Extract exposure times and darks
    exposure_times = np.array(sorted(filtered_dict.keys()))
    darks_adu_per_s = np.array([filtered_dict[exp] for exp in exposure_times])

    # Convert from ADU/s back to ADU by multiplying by exposure time
    darks = np.array([darks_adu_per_s[i] * exposure_times[i]
                      for i in range(len(exposure_times))])  # Shape: (n_exposures, ny, nx)

    # Get image shape
    ny, nx = darks.shape[1], darks.shape[2]

    # Initialize stuck pixel map
    stuck_map = np.zeros((ny, nx), dtype=bool)

    # Loop over each pixel
    for i in range(ny):
        for j in range(nx):
            # Extract pixel values across all exposure times
            pixel_values = darks[:, i, j]

            # Create valid mask: not NaN, not inf, not saturated
            valid_mask = (~np.isnan(pixel_values) &
                          ~np.isinf(pixel_values) &
                          (pixel_values <= effective_saturation_map[i, j]))

            valid_exposures = exposure_times[valid_mask]
            valid_values = pixel_values[valid_mask]

            # Check if we have enough valid points to fit
            if len(valid_exposures) >= 2:
                # Fit line and extract slope (dark current in ADU/s)
                coeffs = np.polyfit(valid_exposures, valid_values, deg=1)
                slope = coeffs[0]

                # Flag if dark current below threshold
                if slope < threshold:
                    stuck_map[i, j] = True
            # If < 2 valid points, leave as False (not stuck)

    return stuck_map


def detect_saturated_pixels(master_flat, effective_saturation_map):
    """
    Detect saturated pixels in the master flat.

    Parameters:
    -----------
    master_flat : 2D array
        Bias-subtracted master flat field
    effective_saturation_map : 2D array
        Per-pixel saturation threshold (saturation_threshold - master_bias)

    Returns:
    --------
    saturated_map : 2D boolean array
        True = saturated pixel
    """

    # Flag pixels where flat value exceeds the saturation threshold
    saturated_map = master_flat > effective_saturation_map

    return saturated_map


def detect_flatbad_pixels(master_flat, flat_low_value):
    """
    Detect flatbad pixels with low response in the master flat.

    Parameters:
    -----------
    master_flat : 2D array
        Bias-subtracted master flat field
    flat_low_value : float
        Threshold as fraction of flat median (e.g., 0.1)

    Returns:
    --------
    flatbad_map : 2D boolean array
        True = flatbad pixel
    """

    # Calculate the threshold
    flat_median = np.median(master_flat)
    threshold = flat_low_value * flat_median

    # Flag pixels below threshold
    flatbad_map = master_flat < threshold

    return flatbad_map


def detect_hot_pixels(dark_dict, effective_saturation_map, sigma=10):
    """
    Detect hot pixels as extreme outliers in the shortest exposure dark.

    Parameters:
    -----------
    dark_dict : dict
        Keys are exposure times (seconds), values are 2D bias-subtracted submaster dark arrays
    effective_saturation_map : 2D array
        Per-pixel saturation threshold (saturation_threshold - master_bias)
    sigma : float
        Number of standard deviations above median to flag (default 10)

    Returns:
    --------
    hot_map : 2D boolean array
        True = hot pixel
    """

    # Check if dark_dict is empty
    if len(dark_dict) == 0:
        print("WARNING: No darks available for hot pixel detection")
        return np.zeros(effective_saturation_map.shape, dtype=bool)

    # Filter to only numeric exposure times
    numeric_keys = [exp for exp in dark_dict.keys() if isinstance(exp, (int, float))]

    if len(numeric_keys) == 0:
        print("WARNING: No numeric exposure times available for hot pixel detection")
        return np.zeros(effective_saturation_map.shape, dtype=bool)

    # Find shortest exposure time
    min_exposure = min(numeric_keys)
    shortest_dark = dark_dict[min_exposure]

    # Calculate threshold
    dark_median = np.median(shortest_dark)
    dark_std = np.std(shortest_dark)
    threshold = dark_median + sigma * dark_std

    # Flag pixels above threshold
    hot_map = shortest_dark > threshold

    return hot_map


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
            flat_low_value = bpm_config.get('flat_low_value', 0.1)
            hot_sigma = bpm_config.get('hot_sigma', 10)
            stuck_threshold = bpm_config.get('stuck_threshold', -10.)
            # Read overwrite setting from config
            overwrite_existing = bpm_config.get('overwrite_existing', True)

            with open_fits_file(caldir + biasname) as hdul:
                master_bias = hdul[0].data
            effective_saturation_map = params['saturation_threshold'] - master_bias

            # Stuck pixels
            stuck_filename = caldir + "stuck_pixel_map.fits"
            if not overwrite_existing and os.path.exists(stuck_filename):
                print(f"Loading existing stuck pixel map from {stuck_filename}")
                stuck_pixel_map = fits.getdata(stuck_filename).astype(bool)
            else:
                print("Creating stuck pixel map")
                stuck_pixel_map = detect_stuck_pixels(dark_dict, effective_saturation_map, threshold=stuck_threshold)
                fits.writeto(stuck_filename, stuck_pixel_map.astype(int), overwrite=True)

            # Saturated pixels
            saturated_filename = caldir + "saturated_pixel_map.fits"
            if not overwrite_existing and os.path.exists(saturated_filename):
                print(f"Loading existing saturated pixel map from {saturated_filename}")
                saturated_pixel_map = fits.getdata(saturated_filename).astype(bool)
            else:
                print("Creating saturated pixel map")
                saturated_pixel_map = detect_saturated_pixels(flat, effective_saturation_map)
                fits.writeto(saturated_filename, saturated_pixel_map.astype(int), overwrite=True)

            # Flatbad pixels
            flatbad_filename = caldir + "flatbad_pixel_map.fits"
            if not overwrite_existing and os.path.exists(flatbad_filename):
                print(f"Loading existing flatbad pixel map from {flatbad_filename}")
                flatbad_pixel_map = fits.getdata(flatbad_filename).astype(bool)
            else:
                print("Creating flatbad pixel map")
                flatbad_pixel_map = detect_flatbad_pixels(flat, flat_low_value)
                fits.writeto(flatbad_filename, flatbad_pixel_map.astype(int), overwrite=True)

            # Hot pixels
            hot_filename = caldir + "hot_pixel_map.fits"
            if not overwrite_existing and os.path.exists(hot_filename):
                print(f"Loading existing hot pixel map from {hot_filename}")
                hot_pixel_map = fits.getdata(hot_filename).astype(bool)
            else:
                print("Creating hot pixel map")
                hot_pixel_map = detect_hot_pixels(dark_dict, effective_saturation_map, sigma=hot_sigma)
                fits.writeto(hot_filename, hot_pixel_map.astype(int), overwrite=True)

            # Noisy pixels
            noisy_filename = caldir + "noisy_pixel_map.fits"
            if os.path.exists(noisy_filename):
                print(f"Loading noisy pixel map from {noisy_filename}")
                noisy_pixel_map = fits.getdata(noisy_filename).astype(bool)
            else:
                print(f"WARNING: Noisy pixel map not found at {noisy_filename}")
                # Create empty map with appropriate shape
                noisy_pixel_map = np.zeros_like(stuck_pixel_map, dtype=bool)

            # After combining all maps
            bad_pixel_map = stuck_pixel_map | saturated_pixel_map | flatbad_pixel_map | hot_pixel_map | noisy_pixel_map

            # Calculate total pixels
            total_pixels = bad_pixel_map.size

            # Print statistics for each type
            print("\n" + "=" * 60)
            print("BAD PIXEL STATISTICS")
            print("=" * 60)

            maps = {
                'Stuck': stuck_pixel_map,
                'Saturated': saturated_pixel_map,
                'Flatbad': flatbad_pixel_map,
                'Hot': hot_pixel_map,
                'Noisy': noisy_pixel_map,
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
                'Stuck': stuck_pixel_map & ~(saturated_pixel_map | flatbad_pixel_map | hot_pixel_map | noisy_pixel_map),
                'Saturated': saturated_pixel_map & ~(
                            stuck_pixel_map | flatbad_pixel_map | hot_pixel_map | noisy_pixel_map),
                'Flatbad': flatbad_pixel_map & ~(
                            stuck_pixel_map | saturated_pixel_map | hot_pixel_map | noisy_pixel_map),
                'Hot': hot_pixel_map & ~(stuck_pixel_map | saturated_pixel_map | flatbad_pixel_map | noisy_pixel_map),
                'Noisy': noisy_pixel_map & ~(stuck_pixel_map | saturated_pixel_map | flatbad_pixel_map | hot_pixel_map)
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

    # inlist = str(sys.argv[1])
    # caldir = str(sys.argv[2]) + '/'
    # outdir = str(sys.argv[3]) + '/'
    # darknames = [caldir + d for d in str(sys.argv[4])]
    # flatnames = [caldir + f for f in str(sys.argv[5])]
    # bpmnames = str(sys.argv[6])
    bad_pixel_maps(inlist, caldir, outdir, darknames, flatnames, bpmname, biasname)

if __name__ == '__main__':
    main()

