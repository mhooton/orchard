import argparse
import math
import sys
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from calibration.pipeutils import get_instrument_parameters, apply_image_slice, open_fits_file, frame_telescope
from calibration import readout
from utils import orchard_version
from reporting.QC import main as QC
import pipeutils
import fnmatch
from collections import defaultdict
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)

CHUNK = 50          # frames combined at once (memory)
RATE_RATIO = 5.0    # an exposure's master dark rate this far from the others' median is reported


def darkmaker(inlist, biasname, darkname, outdir, reportdir=None, run=None, targ=None, use_existing=False):
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
        telescope = frame_telescope(hdul)

    # the night's gain for the mode of its own calibration frames (T2 measured it)
    if run is not None:
        gain = readout.calib_values(readout.load_or_measure(outdir, run, params, telescope))['GAIN']

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

    # Group all dark files by exposure time.  With exposure-matched darks the group is the exposure rounded to
    # whole seconds, the key pipered.py looks a science frame's dark up by and the name of the master
    # (MasterDark_<n>s.fits): grouping by the exact value made e.g. 42.407 s and 42.408 s darks two masters of
    # the same name, the second overwriting the first.
    matching = params.get('dark_matching_exptime', False)
    exposure_groups = defaultdict(list)

    for line in open(outdir + 'removeindexlist.dat'):
        call = line.strip()
        for filename in open(call):
            filename = filename.strip()
            try:
                with open_fits_file(filename) as hdulist:
                    exposure_time = hdulist[0].header['exptime']
                    key = int(round(exposure_time)) if matching else exposure_time
                    exposure_groups[key].append((call, filename))
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
    rates = {}
    # the master MasterDark.fits becomes: the combined one, or with exposure-matched darks the highest-exposure
    # one, which pipered.py scales to a science exposure no master matches and pipeflat.py to every flat
    fallback = None

    # Process each group
    for group_key, file_pairs in groups_to_process.items():
        print(f"Processing {len(file_pairs)} dark files for group: {group_key}")

        # Chunks of at most CHUNK frames are cut within the group, as evenly as possible, and their masters
        # averaged weighted by frame count.  Chunking by position in the whole list (the dsorted files) could leave
        # one frame of a group alone in a chunk, and a plain mean of the chunk masters then gave it half the weight:
        # on Callisto 20260309 one 8 s dark with a BZERO fault made MasterDark_8s 2,079 ADU/s against ~28.
        filenames_in_group = [filename for _, filename in file_pairs]
        n_chunks = max(1, int(math.ceil(len(filenames_in_group) / float(CHUNK))))
        chunk_groups = {i: list(c) for i, c in enumerate(np.array_split(np.array(filenames_in_group, dtype=object),
                                                                          n_chunks)) if len(c)}

        mastermatrix = []
        weights = []

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
                weights.append(len(datamatrix))

        if len(mastermatrix) == 0:
            print(f"WARNING: No valid dark files for group {group_key}")
            continue

        # Average all chunk masters, weighted by their number of frames, to create the final master dark
        dark = np.float32(np.average([np.ma.getdata(m) for m in mastermatrix], axis=0, weights=weights))
        rates[group_key] = float(np.median(dark))

        # Determine output filename
        if params.get('dark_matching_exptime', False):
            # Individual exposure time master
            base_name = (outdir + darkname).split('.')[0]
            outname = f"{base_name}_{int(round(group_key))}s.fits"
        else:
            # Combined master
            outname = outdir + darkname

        # Write master dark
        if use_existing and os.path.exists(outname):
            print(f"Skipping T3: existing master dark found at {outname}")
            output_files.append(outname)
            if fallback is None or (group_key != 'combined' and group_key > fallback[0]):
                fallback = (group_key, None)
            continue
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
        orchard_version.stamp(phdu.header)
        phdu.writeto(outname)
        if fallback is None or (group_key != 'combined' and group_key > fallback[0]):
            fallback = (group_key, phdu)

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

    # an exposure's master far from the night's others is almost certainly a bad frame, not dark current
    if len(rates) > 2:
        for key, rate in rates.items():
            others = [r for k, r in rates.items() if k != key]
            ref = float(np.median(others))
            if ref > 0 and (rate > RATE_RATIO * ref or rate < ref / RATE_RATIO):
                print(f"WARNING: master dark {key} s has a median rate of {rate:.1f} ADU/s against {ref:.1f} for "
                      f"the night's other exposures: check its frames")

    # MasterDark.fits: the highest exposure time group for dark matching (the combined master otherwise)
    if fallback is not None and fallback[1] is not None:
        fallback[1].writeto(outdir + darkname, overwrite=True)
        if matching:
            print(f"MasterDark.fits is the {fallback[0]} s master (highest exposure)")

    # Clean up temporary files
    for dsorted_fn in glob.glob(outdir + 'dsorted*'):
        os.remove(dsorted_fn)
    if os.path.exists(outdir + 'removeindexlist.dat'):
        os.remove(outdir + 'removeindexlist.dat')

    return output_files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inlist')
    parser.add_argument('biasname')
    parser.add_argument('darkname')
    parser.add_argument('outdir')
    parser.add_argument('reportdir')
    parser.add_argument('run')
    parser.add_argument('targ')
    parser.add_argument('--use-existing-calibrations', action='store_true', default=False)
    args = parser.parse_args()

    darkmaker(args.inlist, args.biasname, args.darkname, args.outdir + '/', args.reportdir + '/', args.run, args.targ,
              use_existing=args.use_existing_calibrations)


if __name__ == '__main__':
    main()
