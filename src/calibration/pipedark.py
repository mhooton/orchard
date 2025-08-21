import math
import sys
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from calibration.pipeutils import extract_overscan, image_trim, open_fits_file
from reporting.QC import main as QC
import pipeutils
import fnmatch
from collections import defaultdict


def darkmaker(inlist, biasname, darkname, outdir, reportdir=None, gain=None, run=None, targ=None):
    """
    Create master dark frame(s) from a list of dark images.

    Args:
        inlist (str): Path to file containing list of dark image filenames
        biasname (str): Filename of master bias
        darkname (str): Output filename for master dark
        outdir (str): Output directory for master dark(s)
        gain (float): Gain value for QC calculations (optional)
        reportdir (str): Report directory for QC (optional)
        run (str): Run identifier for QC (optional)
        targ (str): Target name for QC (optional)

    Returns:
        list: List of created master dark filenames
    """

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
    if inst == 'spirit':
        groups_to_process = {exp_time: files for exp_time, files in exposure_groups.items()}
    else:
        # Combine all exposure times into one group
        all_files = []
        for files in exposure_groups.values():
            all_files.extend(files)
        groups_to_process = {'combined': all_files}

    output_files = []

    # Process each group
    for group_key, file_pairs in groups_to_process.items():
        print(f"Processing {len(file_pairs)} dark files for group: {group_key}")

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
                        overscan = extract_overscan(hdulist)
                        overscan = sigma_clip(overscan)
                        medoverscan = np.ma.median(overscan)
                        data = image_trim(hdulist)
                        exposure = hdulist[0].header['exptime']

                        # Bias correction then scale to ADU/s
                        if np.shape(bias) == np.shape(data):
                            corrected = (data - medoverscan - bias) / exposure
                        else:
                            print("WARNING: Bias and Dark dimensions do NOT match! Using UNCORRECTED dark images!")
                            print("Master Bias dimensions: ", np.shape(bias))
                            print("Master Dark dimensions: ", np.shape(data))
                            corrected = (data - medoverscan) / exposure

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
        if inst == 'spirit':
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
        if inst == 'spirit':
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

            if inst == 'spirit':
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
    gain = sys.argv[6]
    run = sys.argv[7]
    targ = sys.argv[8]

    darkmaker(inlist, biasname, darkname, outdir, reportdir, gain, run, targ)


if __name__ == '__main__':
    main()
