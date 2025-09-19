import math
import sys
import argparse
import os
import numpy as np
from astropy.io import fits as pyfits
from calibration.pipeutils import extract_overscan, image_trim, open_fits_file
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool
import timeit
from astropy.stats import sigma_clip
import pandas as pd
from itertools import compress
import warnings
from astropy.io.fits.verify import VerifyWarning
import pipeutils
warnings.filterwarnings('ignore', category=VerifyWarning, message=".*Invalid 'BLANK' keyword.*")

def reduce_file(filename, outdir, biasname, bias, dark_dict, flat, gain, ron, overscan, version):
    path, fname = os.path.split(filename)
    outname = outdir + 'proc' + fname.replace("'", "").replace(" ", "--").replace(".fts", ".fits")

    try:
        with open_fits_file(filename) as hdulist:
            exposure = hdulist[0].header['exptime']
            filter = hdulist[0].header['filter']
            filter = filter.replace("'", "")

            # Extract overscan and trim data
            overscan = extract_overscan(hdulist)
            overscan = sigma_clip(overscan, maxiters=None)
            overscan = np.ma.median(overscan)
            data = image_trim(hdulist)

            # Apply corrections to this science image (if there's a masterflat for this filter)
            if filter not in flat.keys():
                print("No flat found for filter: " + filter)
            else:
                # Select appropriate dark based on exposure time
                rounded_exposure = int(round(exposure))
                if rounded_exposure in dark_dict:
                    dark = dark_dict[rounded_exposure]
                    used_darkname = f"MasterDark_{rounded_exposure}s.fits"
                elif 'combined' in dark_dict:
                    dark = dark_dict['combined']
                    used_darkname = "MasterDark.fits"
                else:
                    dark = 0
                    used_darkname = "N/A"
                    print("WARNING: No dark correction applied!")

                if len(dark) == 0 or (hasattr(dark, 'shape') and np.shape(dark) != np.shape(data)):
                    dark = 0
                    print("WARNING: Dark dimensions don't match data!")
                    used_darkname = "N/A"

                if len(bias) == 0 or np.shape(bias) != np.shape(data):
                    bias = 0
                    print("WARNING: Bias dimensions don't match data!")
                    biasname = "N/A"

                # Apply corrections
                if hasattr(dark, 'shape') and dark.shape == data.shape:
                    corrected = (data - overscan - bias - (dark * exposure)) / flat[filter]
                else:
                    corrected = (data - overscan - bias) / flat[filter]

                if pipeutils.detect_instrument(hdulist) == 'spirit':
                    corrected = pipeutils.clean_bad_pixels(corrected, dark, flat[filter])

                corrected = np.float32(corrected)

                hdulist[0].data = corrected
                hdulist[0].header.set('OVERSCAN', overscan)
                hdulist[0].header.add_history('Overscan of ' + str(overscan) + ' subtracted')

                if np.median(bias) != 0 or np.ptp(bias) != 0:
                    hdulist[0].header.add_history('Bias subtracted using ' + str(biasname))
                    hdulist[0].header['RON'] = (float(ron), 'Read out noise (e-s)')

                if (hasattr(dark, 'shape') and (np.median(dark) != 0 or np.ptp(dark) != 0)) or np.isscalar(
                        dark) and dark != 0:
                    hdulist[0].header.add_history('Dark subtracted using ' + str(used_darkname))
                    if pipeutils.detect_instrument(hdulist) == 'spirit':
                        if hasattr(dark, 'shape'):
                            dark_for_err = pipeutils.clean_bad_pixels(dark, dark)
                            hdulist[0].header['DARKCUR'] = (float(gain) * np.nanmedian(dark_for_err),
                                                            'Dark current (e-s per second)')
                        else:
                            hdulist[0].header['DARKCUR'] = (0, 'Dark current (e-s per second)')
                    else:
                        if hasattr(dark, 'shape'):
                            hdulist[0].header['DARKCUR'] = (float(gain) * np.median(dark),
                                                            'Dark current (e-s per second)')
                        else:
                            hdulist[0].header['DARKCUR'] = (0, 'Dark current (e-s per second)')

                hdulist[0].header['GAIN'] = (float(gain), 'Gain used to calculate RON/Dark Cur')
                hdulist[0].header['PV'] = (version, 'Pipeline Version')

                if pipeutils.detect_instrument(hdulist) == 'andor' and hdulist[0].header.get('CRPIX1') != None:
                    hdulist[0].header['CRPIX1'] = hdulist[0].header['CRPIX1'] - 2
                    hdulist[0].header['CRPIX2'] = hdulist[0].header['CRPIX2'] - 22
                    hdulist[0].header.add_history(
                        'Subtracted 2 pixels from CRPIX1 and 22 from CRPIX2 to account for trimming of ANDOR images')

                if os.path.exists(outname):
                    os.remove(outname)
                hdulist.writeto(outname)

    except Exception as e:
        print("*** ERROR with file " + filename + ". Removing this file from analysis.***")
        print(e)
        return np.nan

    return filter


def reducer(inlist, outdir, biasname, darkname, flatname, gain, reddir, version, usebias, usedark):
    flat_dict = {}
    dark_dict = {}

    # Import master bias
    if usebias == "1":
        if os.path.exists(biasname):
            with open_fits_file(biasname) as hdulist:
                bias = hdulist[0].data
            with open(reddir + "readoutnoise.dat", "r") as f:
                ron = f.read()
                print('Readout noise is ' + str(ron) + ' ADU')
            with open(reddir + "overscan.dat", "r") as f:
                overscan = float(f.read())
                print('Overscan is ' + str(overscan) + ' ADU')
        else:
            print("No Bias Used")
            bias = []
            ron = 0
            overscan = 0
    else:
        print("No Bias Used")
        bias = []
        ron = 0
        overscan = 0

    # Import master darks
    if usedark == "1":
        for d in darkname:
            try:
                # Parse exposure time from filename
                basename = os.path.basename(d)
                if '_MasterDark.fits' in basename and basename.endswith('_MasterDark.fits'):
                    # Combined dark file
                    with open_fits_file(d) as hdulist:
                        dark_dict['combined'] = hdulist[0].data
                    print(f"Loaded combined dark: {basename}")
                elif '_MasterDark_' in basename and 's.fits' in basename:
                    # Exposure-specific dark file
                    exp_part = basename.split('_MasterDark_')[1].replace('s.fits', '')
                    exp_time = int(exp_part)
                    with open_fits_file(d) as hdulist:
                        dark_dict[exp_time] = hdulist[0].data
                    print(f"Loaded dark for {exp_time}s exposure: {basename}")
                else:
                    print(f"Warning: Could not parse dark filename: {basename}")
            except Exception as e:
                print(f"Warning: Could not load dark file {d}: {e}")
                continue

        if not dark_dict:
            print("No Dark Used")
    else:
        print("No Dark Used")

    # Import master flats
    for f in flatname:
        filt = f.split('.fits')[0].split('_')[-1]
        with open_fits_file(f) as hdulist:
            flat_dict[filt] = hdulist[0].data
    print(list(flat_dict.keys()))
    print(f"Available dark exposures: {list(dark_dict.keys())}")

    procfarr = {}

    start_time = timeit.default_timer()
    pool = ThreadPool()
    fn = partial(reduce_file, outdir=outdir, biasname=biasname, bias=bias,
                 dark_dict=dark_dict, flat=flat_dict, gain=gain, ron=ron, overscan=overscan, version=version)

    with open(inlist) as infile:
        filenames = [line.strip() for line in infile]

    # Apply reduce_file to each filename
    filters = pool.map(fn, filenames)

    idx = pd.isnull(filters)
    idx2 = [not i for i in idx]
    filters = list(compress(filters, idx2))
    filenames = list(compress(filenames, idx2))

    after_pool = timeit.default_timer()

    for f in range(len(filenames)):
        # if we have a flat for this filter
        if filters[f] in flat_dict.keys():
            if filters[f] in procfarr.keys():
                procfarr[filters[f]].append(outdir + 'proc' + os.path.basename(
                    filenames[f].replace("'", "").replace(" ", "--").replace(".fts", ".fits")))
            else:
                procfarr[filters[f]] = [outdir + 'proc' + os.path.basename(
                    filenames[f].replace("'", "").replace(" ", "--").replace(".fts", ".fits"))]

    procfarr = dict((k, v) for k, v in procfarr.items() if v)

    for k, v in procfarr.items():
        if v != []:
            dfile = outdir + k + '_processed.dat'
            f = open(dfile, 'w')
            # ignore the first image - WHY
            f.write("\n".join(v[1:]))
            f.close()

    total_time = timeit.default_timer() - start_time
    pool_time = after_pool - start_time
    file_time = timeit.default_timer() - after_pool
    print("total time for reduction = " + str(total_time / 60.) + " minutes")
    print("time for pooling = " + str(pool_time / 60.) + " minutes")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inlist')
    parser.add_argument('-b', '--biasname')
    parser.add_argument('-d', '--darkname', nargs='+')
    parser.add_argument('-f', '--flatname', nargs='+')
    parser.add_argument('-c', '--caldir')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-g', '--gain')
    parser.add_argument('-v', '--version')
    parser.add_argument('-ub', '--usebias')
    parser.add_argument('-ud', '--usedark')
    args = parser.parse_args()

    inlist = args.inlist
    caldir = args.caldir + '/'
    outdir = args.outdir + '/'
    biasname = caldir + args.biasname
    darkname = [caldir + d for d in args.darkname]
    flatname = [caldir + f for f in args.flatname]
    gain = args.gain
    version = args.version
    usebias = args.usebias
    usedark = args.usedark

    if os.path.exists(outdir + 'processed.dat'):
        os.remove(outdir + 'processed.dat')

    reducer(inlist, outdir, biasname, darkname, flatname, gain, caldir, version, usebias, usedark)


if __name__ == '__main__':
    main()