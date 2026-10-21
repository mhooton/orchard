"""
Master Flat Field Creation

Creates normalized master flat field frames from dome or twilight flat images.
Processes flats separately per filter, applying bias and dark corrections before
combining multiple exposures into a master flat that corrects for pixel-to-pixel
sensitivity variations and vignetting.

Processing steps:
- Overscan subtraction and bias/dark correction for each flat
- Normalization of each flat by its median flux
- Sigma-clipped median combination to reject outliers (cosmic rays, etc.)
- Per-filter output organization (e.g., MasterFlat_I+z.fits, MasterFlat_zYJ.fits)

The master flats are normalized to unity mean, allowing direct division into
science frames during reduction. With --save-flat-products, each calibrated flat
(flats/proc*) and the per-pixel standard deviation and variance across the flats
(std.fts, variance.fts) are also written. Nothing in the pipeline reads them, so
they are off by default.

This is Stage T4 of the pipeline and produces the flat fields required for
science image reduction (Stage T6).
"""

import math
import sys
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from calibration.pipeutils import (apply_image_slice, open_fits_file, get_jd_header_name, get_instrument_parameters,
                                   frame_telescope)
from calibration import readout
from calibration.smear import correct_frame, shutterless
from calibration.shutter_detect import NightRecords
from utils import orchard_version
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)


def render_total_file(data, fname, nfiles):
    # passed flat_total, totalname, nflat_files
    hdu = pyfits.PrimaryHDU(data)
    hdu.header.set('nfiles', nfiles)
    # replaced clobber with overwrite
    hdu.writeto(fname, overwrite=True)


def reducer(inlist, biasname, darkname, flatname, usedark, usebias, outdir, reportdir, save_flat_products=False):
    biasname = outdir + biasname
    darkname = outdir + darkname

    if save_flat_products:
        os.system('mkdir ' + outdir + 'flats')

    # import master bias
    if usebias == "1":
        try:
            with open_fits_file(biasname) as hdulist:
                bias = hdulist[0].data
        except:
            bias = 0
    else:
        bias = 0

    # import master dark
    if usedark == "1":
        try:
            with open_fits_file(darkname) as hdulist:
                dark = hdulist[0].data
        except:
            dark = 0
    else:
        dark = 0

    # remove temp files
    # os.system('rm -f '+outdir+'datafile.dat')
    if os.path.exists(outdir + 'datafile.dat'):
        os.remove(outdir + 'datafile.dat')
    # os.system('rm -f '+outdir+'variance.fits')
    if os.path.exists(outdir + 'variance.fts'):
        os.remove(outdir + 'variance.fts')
    # os.system('rm -f '+outdir+flatname)
    if os.path.exists(outdir + flatname):
        os.remove(outdir + flatname)
    # os.system('rm -f '+outdir+'std.fits')
    if os.path.exists(outdir + 'std.fts'):
        os.remove(outdir + 'std.fts')
    # os.system('rm -f '+outdir+'expdata.dat')
    if os.path.exists(outdir + 'expdata.dat'):
        os.remove(outdir + 'expdata.dat')
    # I have added this:
    if os.path.exists(outdir + 'processed.dat'):
        os.remove(outdir + 'processed.dat')

    # the night's readout record (T2), read only if a flat was taken shutterless
    runname = os.path.basename(biasname).split('_MasterBias')[0]
    readout_record = None
    bias_median = None
    # the measured shutter state of each flat's own night: this night's from
    # T2, a borrowed night's measured now on that night's flats
    shutter_records = NightRecords(outdir, runname, refresh_borrowed=True)

    # start counter
    frameno = 1

    nflat_files = 0
    datamatrix = []
    datamatrix_dict = {}
    expfile = outdir + 'expdata.dat'

    # for each filter
    for i in inlist:
        filt = (i.split("_")[-1]).split(".")[0]
        filt = filt.replace("'", "")
        print("Create Master Flat in " + filt)
        datamatrix = []
        n_desmeared = 0
        n_left_smeared = 0
        # loop through each flat image
        for line in open(i):  # open(inlist):
            stripped = line.strip()
            with open_fits_file(stripped) as hdulist:
                # obtain the overscan region
                data = apply_image_slice(hdulist, 'trim', hdulist[0].data)
                overscan_data = apply_image_slice(hdulist, 'overscan', 0.)
                if isinstance(overscan_data, np.ndarray) and overscan_data.size > 1:
                    overscan_data = sigma_clip(overscan_data, maxiters=None)
                    overscan = np.ma.median(overscan_data)
                else:
                    overscan = overscan_data  # Will be 0. if no overscan configured
                exposure = hdulist[0].header['exptime']
                header = hdulist[0].header.copy()
                params = get_instrument_parameters(hdulist)
                # A shutterless flat is de-smeared with the readout time of its
                # readout mode, and its clipped rows judged against its own
                # saturation level, as T6 does for science frames.
                frame_readout = None
                detection = shutter_records.for_frame(stripped, frame_telescope(hdulist)) \
                    if params.get('smear') else None
                if shutterless(header, params, detection=detection):
                    if readout_record is None:
                        readout_record = readout.load_or_measure(outdir, runname, params, frame_telescope(hdulist))
                    os_level, os_noise = readout.overscan_signature(hdulist[0].data, params)
                    frame_readout = readout.frame_values(readout_record, os_level, os_noise)
                # SPECULOOS data doesn't have modified julian date, it uses UTC or julian date
                jd_header = get_jd_header_name(hdulist)
                jd = hdulist[0].header[jd_header.lower()]  # or however the header key formatting works

            # # why are we removing 40 more pixels? Just for robustness to remove the chance of bleeding?
            # median_data = np.median(data[20:-20, :])

            # write exposure times to a file called expdata.dat
            f = open(expfile, 'a')
            f.write(str(exposure) + '\n')
            f.close()

            # ensure the size of images are consistent:
            if np.shape(bias) == np.shape(data):
                corrected1 = data - overscan - bias
            else:
                print("WARNING: Bias and Flat dimensions do NOT match! Using UNCORRECTED Flat images!")
                print("Master Bias dimensions: ", np.shape(bias))
                print("Master Flat dimensions: ", np.shape(data))
                corrected1 = data - overscan

            if np.shape(dark) == np.shape(data):
                corrected1 = corrected1 - (dark * exposure)
            else:
                print("WARNING: Dark and Flat dimensions do NOT match! Using UNCORRECTED Flat images!")
                print("Master Dark dimensions: ", np.shape(dark))
                print("Master Flat dimensions: ", np.shape(data))

            # A flat taken shutterless carries a large readout ramp (43% on a
            # 7.75 s Ganymede twilight flat); remove it before normalising.
            pedestal = 0.0
            if frame_readout is not None:
                if np.shape(bias) == np.shape(data):
                    if bias_median is None:
                        bias_median = float(np.median(bias))
                    pedestal = float(overscan) + bias_median
                else:
                    pedestal = float(overscan)
                for key, value, comment in readout.header_cards(frame_readout, pedestal=pedestal):
                    header[key] = (value, comment)
            corrected1 = correct_frame(corrected1, header, params, raw=data, pedestal=pedestal, detection=detection)
            if header.get('SMEARCOR'):
                n_desmeared += 1
                print(f"De-smeared shutterless flat {os.path.basename(stripped)} "
                      f"(tau {header['SMEARTAU']:.3g}, {header.comments['SMEARTRD']}, "
                      f"shutter state from {header['SMEARSRC']})")
            elif 'SMEARERR' in header:
                n_left_smeared += 1

            nflat_files += 1

            fmean = np.ma.mean(sigma_clip(corrected1, maxiters=None))
            fstd = np.ma.std(sigma_clip(corrected1, maxiters=None))

            normalised = np.float32(np.divide(corrected1, fmean))
            path, fname = os.path.split(stripped)
            # write the normalised correct data for this flat file to a new file in /flats/
            outname = outdir + 'flats/' + 'proc' + fname  # .replace('.bz2', '')
            dfile = outdir + 'datafile.dat'
            f = open(dfile, 'a')
            # write counter, mean of flat, std dev of flat, exposure and filename of the new fits file of normalised corrected frame
            f.write(str(frameno) + " " + str(fmean) + " " + str(fstd) + " " + str(exposure) + " " + outname)
            f.close()

            datamatrix.append(normalised)

            if save_flat_products:
                # create new fits image phdu using the normalised corrected data
                phdu = pyfits.PrimaryHDU(normalised)
                phdu.header['exposure'] = exposure
                phdu.header['jd'] = jd
                if os.path.exists(outname):
                    os.remove(outname)

                # write filename of this hdulist to another file 'processed.dat'
                phdu.writeto(outname, overwrite=True)
                tfile = outdir + 'processed.dat'
                f = open(tfile, 'a')
                f.write(outname)
                f.close()

            frameno += 1

        # if dfile is not defined this will run the except - if this is a file we are not including
        # not sure this is the best way to do this since we don't use any of these variables?
        try:
            frame, means, stds = np.loadtxt(dfile, usecols=(0, 1, 2), unpack=True)
        except UnboundLocalError as err:
            if 'dfile' in str(err):
                raise RuntimeError("All flats invalid. Pipeline cannot continue"
                                   ", original error: {}".format(str(err)))

        # for filt,datamatrix in datamatrix_dict.iteritems():
        if save_flat_products:
            wholestd = np.std(datamatrix, axis=0)
            outname = outdir + 'std.fts'
            pyfits.PrimaryHDU(wholestd).writeto(outname, overwrite=True)
            print('std done')

            variance = wholestd ** 2  # 1/(wholestd*wholestd)
            outname = outdir + 'variance.fts'
            pyfits.PrimaryHDU(variance).writeto(outname, overwrite=True)
            print('var done')

        # write normalised, corrected medianed flat frame to master file
        flat = np.median(datamatrix, axis=0)
        normflat = np.divide(flat, np.mean(flat))
        ext = os.path.splitext(flatname)[1]
        flname = os.path.splitext(flatname)[0]
        outname = outdir + flname + '_' + filt + ext
        flat_hdu = pyfits.PrimaryHDU(normflat)
        flat_hdu.header['NSMEARED'] = (n_desmeared, 'Input flats de-smeared (taken shutterless)')
        if n_left_smeared:
            flat_hdu.header['NSMRLEFT'] = (n_left_smeared, 'Shutterless input flats NOT de-smeared')
            print(f"WARNING: {n_left_smeared} shutterless flats in {filt} could not be de-smeared "
                  f"(readout time unknown); the master flat carries their readout ramp")
        orchard_version.stamp(flat_hdu.header)
        flat_hdu.writeto(outname, overwrite=True)
        print('flat done')

    # write non-normalised, corrected flat frame to file 'totalname'
    # render_total_file(flat_total, totalname, nflat_files)


def main():
    # the optional flag is removed first, because the positional arguments are counted from the end
    args = sys.argv[1:]
    save_flat_products = '--save-flat-products' in args
    args = [a for a in args if a != '--save-flat-products']

    inlist = args[:-7]
    biasname = str(args[-5])
    darkname = str(args[-4])
    flatname = str(args[-3])
    usedark = str(args[-6])
    usebias = str(args[-7])
    outdir = str(args[-2]) + '/'
    reportdir = str(args[-1]) + '/'

    print(inlist)
    print(len(inlist))
    print(type(inlist))

    reducer(inlist, biasname, darkname, flatname, usedark, usebias, outdir, reportdir,
            save_flat_products=save_flat_products)


if __name__ == '__main__':
    main()