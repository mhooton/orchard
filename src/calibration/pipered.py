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
warnings.filterwarnings('ignore', category=VerifyWarning, message=".*Invalid 'BLANK' keyword.*")

def reduce_file(filename, outdir,biasname,darkname, bias, dark, flat, gain, ron, overscan, version):
    path, fname = os.path.split(filename)
    outname = outdir + 'proc' + fname.replace("'", "").replace(" ", "--").replace(".fts",".fits")

    try:
        with open_fits_file(filename) as hdulist:
            exposure = hdulist[0].header['exptime']
            filter = hdulist[0].header['filter']
            filter = filter.replace("'", "")

            # print(np.shape(hdulist[0].data))
            # print(filter)
            #open science image

            overscan = extract_overscan(hdulist)
            overscan = sigma_clip(overscan, maxiters=None)
            overscan = np.ma.median(overscan)
            data = image_trim(hdulist)

            # if "zYJ" in filter:
            #     data = hdulist[0].data
            #     overscan = 0
            # else:
            #     if np.shape(hdulist[0].data)[0] > 2048:
            #         overscan = extract_overscan(hdulist)
            #         # sigma clip the overscan
            #         overscan = sigma_clip(overscan, iters=None)
            #         overscan = np.ma.median(overscan)
            #         data = hdulist[0].data[22:2066, 2:2048]
            #     else:
            #         data = hdulist[0].data[2:2046, 2:2048]

            #apply corrections to this science image (if there's a masterflat for this filter)
            if filter not in flat.keys():
                print("No flat found for filter: " + filter)

            else:
                if len(dark)==0 or np.shape(dark)!=np.shape(data):
                    dark = 0
                    print("WARNING: Flat not corrected for dark!")
                    darkname = "N/A"
                if len(bias)==0 or np.shape(bias)!=np.shape(data):
                    bias = 0
                    print("WARNING: Flat not corrected for bias!")
                    biasname = "N/A"

                corrected = (data-overscan-bias-(dark*exposure))/flat[filter]

                corrected = np.float32(corrected)

                hdulist[0].data = corrected
                hdulist[0].header.set('OVERSCAN',overscan)
                hdulist[0].header.add_history('Overscan of '+str(overscan)+' subtracted')
                if np.median(bias)!=0 or np.ptp(bias)!=0:
                    hdulist[0].header.add_history('Bias subtracted using '+str(biasname))
                    hdulist[0].header['RON'] = (float(ron), 'Read out noise (e-s)')
                if np.median(dark)!=0 or np.ptp(dark)!=0:
                    hdulist[0].header.add_history('Dark subtracted using '+str(darkname))
                    hdulist[0].header['DARKCUR'] = (float(gain) * np.median(dark), 'Dark current (e-s per second)')
                # hdulist[0].header.add_history('Flat corrected using '+str(flatname))

                hdulist[0].header['GAIN'] = (float(gain), 'Gain used to calculate RON/Dark Cur')
                hdulist[0].header['PV'] = (version, 'Pipeline Version')

                if hdulist[0].header.get('CRPIX1') != None:
                    hdulist[0].header['CRPIX1'] = hdulist[0].header['CRPIX1'] - 2
                    hdulist[0].header['CRPIX2'] = hdulist[0].header['CRPIX2'] - 22
                    hdulist[0].header.add_history('Subtracted 2 pixels from CRPIX1 and 22 from CRPIX2 to account for trimming of images')
                #command = 'rm -f '+outname
                #os.system(command)
                if os.path.exists(outname):
                    os.remove(outname)
                #write new fits file to new folder
                hdulist.writeto(outname)

    except Exception as e:
        print("*** ERROR with file " + filename + ". Removing this file from analysis.***")
        print(e)
        return np.nan

    return filter #outname

def reducer(inlist,outdir, biasname, darkname, flatname, gain, reddir, version, usebias,usedark):
    flat_dict = {}
    #import master bias
    if usebias=="1":
        if os.path.exists(biasname):
            with open_fits_file(biasname) as hdulist:
                bias = hdulist[0].data
            with open(reddir + "readoutnoise.dat", "r") as f:
                ron = f.read()
                print('Readout noise is '+str(ron)+' ADU')
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
        bias=[]
        ron = 0
        overscan = 0

    #import master dark
    if usedark=="1":
        if os.path.exists(darkname):
            with open_fits_file(darkname) as hdulist:
                dark = hdulist[0].data
        else:
            print("No Dark Used")
            dark = []
    else:
        print("No Dark Used")
        dark=[]
    #import master flat
    for f in flatname:
        filt = f.split('.fits')[0].split('_')[-1]
        with open_fits_file(f) as hdulist:
            flat_dict[filt] = hdulist[0].data
    print(list(flat_dict.keys()))
    # print(np.shape(flat_dict['zYJ']))

    procfarr = {}

    start_time = timeit.default_timer()
    pool = ThreadPool()
    fn = partial(reduce_file, outdir=outdir, biasname=biasname, darkname=darkname, bias=bias, dark=dark, flat=flat_dict, gain=gain, ron=ron, overscan=overscan,version=version)
    with open(inlist) as infile:
        #loop through all science fits filenames
        filenames = [line.strip() for line in infile]

    #apply reduce_file to each filename
    filters = pool.map(fn, filenames)
    # print len(filters)
    # print len(filenames)
    idx = pd.isnull(filters)
    idx2 = [not i for i in idx]
    filters = list(compress(filters, idx2))
    filenames = list(compress(filenames, idx2))

    after_pool = timeit.default_timer()

    # for k in flat_dict.keys():
    for f in range(len(filenames)):
        # if we have a flat for this filter
        if filters[f] in flat_dict.keys():
            if filters[f] in procfarr.keys():
                procfarr[filters[f]].append(outdir + 'proc' + os.path.basename(
                    filenames[f].replace("'", "").replace(" ", "--").replace(".fts", ".fits")))
            else:
                procfarr[filters[f]] = [outdir + 'proc' + os.path.basename(
                    filenames[f].replace("'", "").replace(" ", "--").replace(".fts", ".fits"))]
        # procfarr[k] = [outdir+'proc'+os.path.basename(f).replace("'","") for f in filenames if "-"+k+'.fts' in os.path.basename(f).replace("'","")]
        # if ('-'+k in f)]
        # procfarr[k] = [outdir + 'proc' + os.path.basename(f) for f in filenames]
    procfarr = dict((k, v) for k, v in procfarr.items() if v)

    for k,v in procfarr.items():
        if v != []:
            dfile = outdir + k + '_processed.dat'
            f = open(dfile, 'w')
            # ignore the first image - WHY
            f.write("\n".join(v[1:]))
            f.close()

    total_time = timeit.default_timer()-start_time
    pool_time = after_pool-start_time
    file_time = timeit.default_timer()-after_pool
    print("total time for reduction = " + str(total_time/60.) + " minutes")
    print("time for pooling = " + str(pool_time/60.) + " minutes")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inlist')
    parser.add_argument('-b','--biasname')
    parser.add_argument('-d','--darkname')
    parser.add_argument('-f','--flatname', nargs='+')
    parser.add_argument('-c','--caldir')
    parser.add_argument('-o','--outdir')
    parser.add_argument('-g','--gain')
    parser.add_argument('-v', '--version')
    parser.add_argument('-ub', '--usebias')
    parser.add_argument('-ud', '--usedark')
    args = parser.parse_args()

    inlist = args.inlist
    # caldir is the location of thed master bias/dark/flat fits files
    caldir = args.caldir + '/'
    # outdir is the output location for the reduced images
    outdir = args.outdir + '/'
    biasname = caldir + args.biasname
    darkname = caldir + args.darkname
    # flatname = caldir+flatname
    # print args.flatname
    flatname = [caldir + f for f in args.flatname]
    gain = args.gain
    version = args.version
    usebias = args.usebias
    usedark = args.usedark

    if os.path.exists(outdir + 'processed.dat'):
        os.remove(outdir + 'processed.dat')

    reducer(inlist, outdir, biasname, darkname, flatname,gain,caldir, version, usebias, usedark)


if __name__ == '__main__':
    main()