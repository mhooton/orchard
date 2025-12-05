import math
import sys
import argparse
import os
import numpy as np
from astropy.io import fits as pyfits
from calibration.pipeutils import (get_instrument_parameters, apply_image_slice, clean_bad_pixels, open_fits_file,
                                   create_dark_dict, create_flat_dict)
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool
import timeit
from astropy.stats import sigma_clip
import pandas as pd
from itertools import compress
import warnings
from astropy.io.fits.verify import VerifyWarning

from calibration.pipeutils import create_dark_dict

warnings.filterwarnings('ignore', category=VerifyWarning, message=".*Invalid 'BLANK' keyword.*")

def reduce_file(filename, outdir, params, biasname, bias, dark_dict, flat_dict, bpm, ron, version):

    path, fname = os.path.split(filename)
    outname = outdir + 'proc' + fname.replace("'", "").replace(" ", "--").replace(".fts", ".fits")

    # try:
    with open_fits_file(filename) as hdulist:
        exposure = hdulist[0].header['exptime']
        filter = hdulist[0].header['filter']
        filter = filter.replace("'", "")

        # Extract overscan and trim data using apply_image_slice
        overscan_data = apply_image_slice(hdulist, 'overscan', 0.)
        if isinstance(overscan_data, np.ndarray) and overscan_data.size > 1:
            overscan_data = sigma_clip(overscan_data, maxiters=None)
            overscan = np.ma.median(overscan_data)
        else:
            overscan = overscan_data  # Will be 0. if no overscan configured

        data = apply_image_slice(hdulist, 'trim', hdulist[0].data)

        # Get gain from params instead of parameter
        gain = params['gain']

        # Apply corrections to this science image (if there's a masterflat for this filter)
        if filter not in flat_dict.keys():
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
                corrected = (data - overscan - bias - (dark * exposure)) / flat_dict[filter]
            else:
                corrected = (data - overscan - bias) / flat_dict[filter]

            # print("BPM type:", type(bpm))
            # print("BPM dtype:", bpm.dtype)
            # print("BPM unique values:", np.unique(bpm))
            # print("BPM shape:", bpm.shape)
            # print("Image shape:", corrected.shape)
            if bpm is not None:
                # print("Corrected type:", type(corrected))
                # print("Corrected dtype:", corrected.dtype)
                # if isinstance(corrected, np.ma.MaskedArray):
                #     print("Corrected is a masked array")
                #     print("Number of masked elements:", np.ma.count_masked(corrected))
                # else:
                #     print("Corrected is NOT a masked array")
                corrected = clean_bad_pixels(corrected, bpm)

            corrected = np.float32(corrected)

            hdulist[0].data = corrected
            hdulist[0].header.set('OVERSCAN', overscan)
            hdulist[0].header.add_history('Overscan of ' + str(overscan) + ' subtracted')

            if np.median(bias) != 0 or np.ptp(bias) != 0:
                hdulist[0].header.add_history('Bias subtracted using ' + str(biasname))
                hdulist[0].header['RON'] = (float(ron), 'Read out noise (e-s)')

            # if (hasattr(dark, 'shape') and (np.median(dark) != 0 or np.ptp(dark) != 0)) or np.isscalar(
            #         dark) and dark != 0:
            #     hdulist[0].header.add_history('Dark subtracted using ' + str(used_darkname))
            #     if hasattr(dark, 'shape'):
            #         if params.get('bad_pixel_correction', False):
            #             dark_for_err = clean_bad_pixels(dark, dark, flat_dict[filter])
            #             hdulist[0].header['DARKCUR'] = (float(gain) * np.nanmedian(dark_for_err),
            #                                             'Dark current (e-s per second)')
            #         else:
            #             hdulist[0].header['DARKCUR'] = (float(gain) * np.median(dark),
            #                                             'Dark current (e-s per second)')
            #     else:
            #         hdulist[0].header['DARKCUR'] = (0, 'Dark current (e-s per second)')

            hdulist[0].header['GAIN'] = (float(gain), 'Gain used to calculate RON/Dark Cur')
            hdulist[0].header['PV'] = (version, 'Pipeline Version')

            if 'trim' in params and hdulist[0].header.get('CRPIX1') is not None:
                trim_config = params['trim']
                left_col = trim_config.get('left_col', 0)
                top_row = trim_config.get('top_row', 0)

                hdulist[0].header['CRPIX1'] = hdulist[0].header['CRPIX1'] - left_col
                hdulist[0].header['CRPIX2'] = hdulist[0].header['CRPIX2'] - top_row
                hdulist[0].header.add_history(
                    f'Subtracted {left_col} pixels from CRPIX1 and {top_row} from CRPIX2 to account for image trimming')

            if os.path.exists(outname):
                os.remove(outname)
            hdulist.writeto(outname)

    # except Exception as e:
    #     print("*** ERROR with file " + filename + ". Removing this file from analysis.***")
    #     print(e)
    #     return np.nan

    return filter


def reducer(inlist, outdir, biasname, darknames, flatnames, bpmname, reddir, version, usebias, usedark):
    # Get parameters once from the first file
    with open(inlist) as infile:
        first_filename = infile.readline().strip()

    with open_fits_file(first_filename) as hdul:
        params = get_instrument_parameters(hdul)

    # Import master bias
    if usebias == "1":
        if os.path.exists(biasname):
            with open_fits_file(biasname) as hdulist:
                bias = hdulist[0].data
            with open(reddir + "readoutnoise.dat", "r") as f:
                ron = f.read()
                print('Readout noise is ' + str(ron) + ' ADU')
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
        dark_dict = create_dark_dict(darknames)
    else:
        print("No Dark Used")
        dark_dict = {}

    print(f"Available dark exposures: {list(dark_dict.keys())}")

    # Import master flats
    flat_dict = create_flat_dict(flatnames)

    if 'bad_pixel_correction' in params:
        if os.path.exists(bpmname):
            with open_fits_file(bpmname) as hdulist:
                bpm = hdulist[0].data
                bpm = bpm.astype(bool)
            print("Bad pixel map found. Applying bad pixel correction.")
        else:
            print("Bad pixel map not found. Continuing without bad pixel correction.")
            bpm = None

    procfarr = {}

    start_time = timeit.default_timer()
    pool = ThreadPool()
    fn = partial(reduce_file, outdir=outdir, params=params, biasname=biasname, bias=bias,
                 dark_dict=dark_dict, flat_dict=flat_dict, bpm=bpm, ron=ron, version=version)

    with open(inlist) as infile:
        filenames = [line.strip() for line in infile]

    # Apply reduce_file to each filename
    filters = pool.map(fn, filenames)

    # with open(inlist) as infile:
    #     filenames = [line.strip() for line in infile]
    #
    # filters = []
    # for filename in filenames:
    #     result = reduce_file(filename, outdir=outdir, params=params, biasname=biasname,
    #                          bias=bias, dark_dict=dark_dict, flat_dict=flat_dict,
    #                          bpm=bpm, ron=ron, version=version)
    #     filters.append(result)

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
    parser.add_argument('-d', '--darknames', nargs='+')
    parser.add_argument('-f', '--flatnames', nargs='+')
    parser.add_argument('-bpm', '--bpmname')
    parser.add_argument('-c', '--caldir')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-v', '--version')
    parser.add_argument('-ub', '--usebias')
    parser.add_argument('-ud', '--usedark')
    args = parser.parse_args()

    inlist = args.inlist
    caldir = args.caldir + '/'
    outdir = args.outdir + '/'
    biasname = caldir + args.biasname
    darknames = [caldir + d for d in args.darknames]
    flatnames = [caldir + f for f in args.flatnames]
    bpmname = caldir + args.bpmname
    version = args.version
    usebias = args.usebias
    usedark = args.usedark

    if os.path.exists(outdir + 'processed.dat'):
        os.remove(outdir + 'processed.dat')

    reducer(inlist, outdir, biasname, darknames, flatnames, bpmname, caldir, version, usebias, usedark)

if __name__ == '__main__':
    main()