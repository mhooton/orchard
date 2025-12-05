import math
import sys
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from calibration.pipeutils import get_instrument_parameters, apply_image_slice, open_fits_file
from reporting.QC import main as QC
import pipeutils
import fnmatch
from collections import defaultdict

def detect_noisy_pixels_single_exposure(frames, sigma=5, min_frames=5):
    """
    Detect noisy pixels for a single exposure time batch.

    Parameters:
    -----------
    frames : list of 2D arrays
        Individual dark frames at the same exposure time (bias-subtracted)
    sigma : float
        Number of standard deviations above median to flag (default 5)
    min_frames : int
        Minimum number of frames required (default 5)

    Returns:
    --------
    noisy_map : 2D boolean array
        True = noisy pixel
    """

    # Check if we have enough frames
    if len(frames) < min_frames:
        print(f"WARNING: Need at least {min_frames} frames for noisy pixel detection, found {len(frames)}")
        # Return all-False array with same shape as frames
        return np.zeros(frames[0].shape, dtype=bool)

    # Stack frames into 3D array (n_frames, ny, nx)
    frames_array = np.array(frames)

    # Calculate temporal std for each pixel
    temporal_std_map = np.std(frames_array, axis=0)

    # Calculate threshold from distribution of temporal stds
    median_std = np.median(temporal_std_map)
    std_of_stds = np.std(temporal_std_map)
    threshold = median_std + sigma * std_of_stds

    # Flag pixels above threshold
    noisy_map = temporal_std_map > threshold

    return noisy_map

def darkmaker(inlist, biasname, darkname, outdir, reportdir=None, run=None, targ=None):
    """
    Create master dark frame(s) from a list of dark images.

    Args:
        inlist (str): Path to file containing list of dark image filenames
        biasname (str): Filename of master bias
        darkname (str): Output filename for master dark
        outdir (str): Output directory for master dark(s)
        reportdir (str): Report directory for QC (optional)
        run (str): Run identifier for QC (optional)
        targ (str): Target name for QC (optional)

    Returns:
        list: List of created master dark filenames
    """

    # Get parameters once from the first file
    with open(inlist) as infile:
        first_filename = infile.readline().strip()

    with open_fits_file(first_filename) as hdul:
        params = get_instrument_parameters(hdul)
        gain = params['gain']

    # Clean up any existing sorted files
    for dsorted_fn in glob.glob(outdir + 'dsorted*'):
        os.remove(dsorted_fn)

    # Read and chunk input files (50 files per chunk)
    position = 0
    i = 1
    lines = []

    for line in open(inlist):
        lines.append(line.strip())
        fname = outdir + 'dsorted' + "{0:03d}".format(position)
        with open(fname, 'a') as f:
            f.write(line)
        if i == 50:
            i = 0
            position += 1
        i += 1

    # Create list of chunk files
    with open(outdir + 'removeindexlist.dat', 'w') as rem:
        for fn in glob.glob(outdir + 'dsorted*'):
            rem.write(fn + "\n")

    if len(lines) == 0:
        print("WARNING: No dark images!")
        return []

    # Detect instrument from first file
    first_filename = lines[0]
    with open_fits_file(first_filename) as hdulist:
        inst = pipeutils.detect_instrument(hdulist)

    print(f"Detected instrument: {inst}")

    # Load master bias
    bias_path = outdir + biasname
    with open_fits_file(bias_path) as hdulist:
        bias = hdulist[0].data

    # Group all dark files by exposure time
    exposure_groups = defaultdict(list)

    for line in open(outdir + 'removeindexlist.dat'):
        call = line.strip()
        for filename in open(call):
            filename = filename.strip()
            try:
                with open_fits_file(filename) as hdulist:
                    exposure_time = hdulist[0].header['exptime']
                    exposure_groups[exposure_time].append((call, filename))
            except Exception as e:
                print(f"Warning: Could not read exposure time from {filename}: {e}")
                continue

    # Determine which exposure times to process
    if params.get('dark_matching_exptime', False):
        groups_to_process = {exp_time: files for exp_time, files in exposure_groups.items()}
    else:
        # Combine all exposure times into one group
        all_files = []
        for files in exposure_groups.values():
            all_files.extend(files)
        groups_to_process = {'combined': all_files}

    output_files = []

    # Initialize noisy pixel map as None
    noisy_pixel_map = None

    # Read noisy pixel detection parameters
    if 'bad_pixel_correction' in params:
        bpm_config = params['bad_pixel_correction']
        noisy_sigma = bpm_config.get('noisy_threshold_sigma', 5)
        noisy_min_frames = bpm_config.get('noisy_min_frames', 5)
        print(f"Noisy pixel detection enabled: sigma={noisy_sigma}, min_frames={noisy_min_frames}")
    else:
        noisy_sigma = None
        noisy_min_frames = None
        print("Noisy pixel detection disabled (no bad_pixel_correction config)")

    # Process each group
    for group_key, file_pairs in groups_to_process.items():
        print(f"Processing {len(file_pairs)} dark files for group: {group_key}")

        # Noisy pixel detection (only if configured)
        if noisy_sigma is not None and len(file_pairs) >= noisy_min_frames:
            # Noisy pixel detection for this exposure time
            if len(file_pairs) >= noisy_min_frames:
                frames = []
                for chunk_file, filename in file_pairs:
                    try:
                        with open_fits_file(filename) as hdulist:
                            data = apply_image_slice(hdulist, 'trim', hdulist[0].data)
                            overscan_data = apply_image_slice(hdulist, 'overscan', 0.)
                            if isinstance(overscan_data, np.ndarray) and overscan_data.size > 1:
                                overscan_data = sigma_clip(overscan_data, maxiters=None)
                                overscan = np.ma.median(overscan_data)
                            else:
                                overscan = overscan_data

                            # Bias correction (no exposure scaling for noisy detection)
                            if np.shape(bias) == np.shape(data):
                                corrected = data - overscan - bias
                            else:
                                corrected = data - overscan

                            frames.append(corrected)
                    except Exception as e:
                        print(f"Warning: Could not load {filename} for noisy pixel detection: {e}")
                        continue

                # Detect noisy pixels for this exposure time
                if len(frames) > 0:
                    noisy_map_this_exp = detect_noisy_pixels_single_exposure(frames, sigma=noisy_sigma,
                                                                             min_frames=noisy_min_frames)

                    # Combine with accumulated map
                    if noisy_pixel_map is None:
                        noisy_pixel_map = noisy_map_this_exp
                    else:
                        noisy_pixel_map = noisy_pixel_map | noisy_map_this_exp

        # Group files by chunk file for this exposure time
        chunk_groups = defaultdict(list)
        for chunk_file, filename in file_pairs:
            chunk_groups[chunk_file].append(filename)

        mastermatrix = []

        # Process each chunk
        for chunk_file, filenames in chunk_groups.items():
            datamatrix = []

            for filename in filenames:
                try:
                    with open_fits_file(filename) as hdulist:
                        # obtain the overscan region
                        data = apply_image_slice(hdulist, 'trim', hdulist[0].data)
                        overscan_data = apply_image_slice(hdulist, 'overscan', 0.)
                        if isinstance(overscan_data, np.ndarray) and overscan_data.size > 1:
                            overscan_data = sigma_clip(overscan_data, maxiters=None)
                            overscan = np.ma.median(overscan_data)
                        else:
                            overscan = overscan_data  # Will be 0. if no overscan configured
                        exposure = hdulist[0].header['exptime']

                        # Bias correction then scale to ADU/s
                        if np.shape(bias) == np.shape(data):
                            corrected = (data - overscan - bias) / exposure
                        else:
                            print("WARNING: Bias and Dark dimensions do NOT match! Using UNCORRECTED dark images!")
                            print("Master Bias dimensions: ", np.shape(bias))
                            print("Master Dark dimensions: ", np.shape(data))
                            corrected = (data - overscan) / exposure

                        datamatrix.append(corrected)

                except Exception as e:
                    print(f"Warning: Could not process {filename}: {e}")
                    continue

            # Create master for this chunk if we have data
            if len(datamatrix) > 0:
                master = np.ma.median(sigma_clip(datamatrix, axis=0), axis=0)
                mastermatrix.append(master)

        if len(mastermatrix) == 0:
            print(f"WARNING: No valid dark files for group {group_key}")
            continue

        # Average all chunk masters to create final master dark
        dark = np.float32(np.mean(mastermatrix, axis=0))

        # Determine output filename
        if params.get('dark_matching_exptime', False):
            # Individual exposure time master
            base_name = (outdir + darkname).split('.')[0]
            outname = f"{base_name}_{int(round(group_key))}s.fits"
        else:
            # Combined master
            outname = outdir + darkname

        # Write master dark
        if os.path.exists(outname):
            os.remove(outname)

        phdu = pyfits.PrimaryHDU(dark)
        phdu.header.add_history(f'Master dark created from {sum(len(chunk) for chunk in chunk_groups.values())} images')
        if params.get('dark_matching_exptime', False):
            phdu.header.add_history(f'Exposure time: {group_key}s')
            phdu.header['EXPTIME'] = group_key
        else:
            phdu.header.add_history('Combined from all exposure times')
            exp_times = list(exposure_groups.keys())
            phdu.header.add_history(f'Exposure times ranged from {min(exp_times)}s to {max(exp_times)}s')
            phdu.header['EXPTIME'] = 1.0
        phdu.header['DARKUNIT'] = 'ADU/s'
        phdu.writeto(outname)

        output_files.append(outname)
        print(f"Created master dark: {outname}")

        # QC calculations if parameters provided
        if gain is not None and reportdir is not None:
            darkc = float(gain) * np.ma.median(sigma_clip(dark))

            if params.get('dark_matching_exptime', False):
                # Individual exposure time QC
                dirsplit = outdir.split("/")
                date = dirsplit[-3]
                field = targ if targ else f"exp_{int(round(group_key))}s"
                tel = dirsplit[-5]
                QC([0, reportdir + "/QC", date, field, tel, run, 'QC7', darkc])
            else:
                # Combined QC
                dirsplit = outdir.split("/")
                date = dirsplit[-3]
                field = targ if targ else "combined"
                tel = dirsplit[-5]
                QC([0, reportdir + "/QC", date, field, tel, run, 'QC7', darkc])

    #Assign final group to outname
    phdu.writeto(outdir + darkname, overwrite=True)

    # Save noisy pixel map if created
    if noisy_pixel_map is not None:
        noisy_fits = pyfits.PrimaryHDU(noisy_pixel_map.astype(int))
        noisy_fits.writeto(outdir + 'noisy_pixel_map.fits', overwrite=True)
        print("Created noisy pixel map")


    # Clean up temporary files
    for dsorted_fn in glob.glob(outdir + 'dsorted*'):
        os.remove(dsorted_fn)
    if os.path.exists(outdir + 'removeindexlist.dat'):
        os.remove(outdir + 'removeindexlist.dat')

    return output_files


def main():
    inlist = str(sys.argv[1])
    biasname = str(sys.argv[2])
    darkname = str(sys.argv[3])
    outdir = str(sys.argv[4]) + '/'
    reportdir = str(sys.argv[5]) + '/'
    run = sys.argv[6]
    targ = sys.argv[7]

    darkmaker(inlist, biasname, darkname, outdir, reportdir, run, targ)


if __name__ == '__main__':
    main()
