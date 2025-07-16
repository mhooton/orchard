import math
import sys
import os
import glob
import numpy as np
from astropy.io import fits as pyfits
from astropy.stats import sigma_clip
from calibration.pipeutils import extract_overscan, image_trim, open_fits_file


def render_total_file(data, fname, nfiles):
    # passed flat_total, totalname, nflat_files
    hdu = pyfits.PrimaryHDU(data)
    hdu.header.set('nfiles', nfiles)
    # replaced clobber with overwrite
    hdu.writeto(fname, overwrite=True)


def reducer(inlist, biasname, darkname, flatname, usedark, usebias, outdir, reportdir):
    biasname = outdir + biasname
    darkname = outdir + darkname

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
        # loop through each flat image
        for line in open(i):  # open(inlist):
            stripped = line.strip()
            with open_fits_file(stripped) as hdulist:
                header = hdulist[0].header
                overscan = extract_overscan(hdulist)
                overscan = sigma_clip(overscan)
                medoverscan = np.ma.median(overscan)
                data = image_trim(hdulist)

                # if filt=="zYJ":
                #     # INFRARED CAMERA - TO DO (does not have overscan, size of images are different)
                #     data = hdulist[0].data
                #     medoverscan = 0
                # else:
                #     if np.shape(hdulist[0].data)[0] > 2048:
                #         overscan = extract_overscan(hdulist)
                #         overscan = sigma_clip(overscan)
                #         medoverscan = np.ma.median(overscan)
                #         # remove the pre/overscans from the y axis (22 pixels top and bottom):
                #         # also remove 2 bad columns from LHS
                #         data = hdulist[0].data[22:2066,2:2048]
                #     else:
                #         # extract median value of the overscan
                #         with open(outdir + "overscan.dat", "r") as f:
                #             medoverscan = float(f.read())
                #         # remove the pre/overscans from the y axis (22 pixels top and bottom):
                #         # also remove 2 bad columns from LHS
                #         data = hdulist[0].data[2:2046, 2:2048]

                # exposure time:
                exposure = header['exptime']
                # SPECULOOS data doesn't have modified julian date, it uses UTC or julian date
                if header['observer'] == "Astra" or "MOANA" in header['telescop']:
                    jd = header['jd-obs']
                else:  # ACP flats do not have JD-OBS, so revert to JD (rounded to 10 decimal places rather than 7 for JD-OBS)
                    jd = header['jd']

            # why are we removing 40 more pixels? Just for robustness to remove the chance of bleeding?
            median_data = np.median(data[20:-20, :])

            # write exposure times to a file called expdata.dat
            f = open(expfile, 'a')
            f.write(str(exposure) + '\n')
            f.close()

            # ensure the size of images are consistent:
            if np.shape(bias) == np.shape(data):
                corrected1 = data - np.ma.median(medoverscan) - bias
            else:
                print("WARNING: Bias and Flat dimensions do NOT match! Using UNCORRECTED Flat images!")
                print("Master Bias dimensions: ", np.shape(bias))
                print("Master Flat dimensions: ", np.shape(data))
                corrected1 = data - np.ma.median(medoverscan)

            if np.shape(dark) == np.shape(data):
                corrected1 = corrected1 - (dark * exposure)
            else:
                print("WARNING: Dark and Flat dimensions do NOT match! Using UNCORRECTED Flat images!")
                print("Master Dark dimensions: ", np.shape(dark))
                print("Master Flat dimensions: ", np.shape(data))

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
        pyfits.PrimaryHDU(normflat).writeto(outname, overwrite=True)
        print('flat done')

    # write non-normalised, corrected flat frame to file 'totalname'
    # render_total_file(flat_total, totalname, nflat_files)


def main():
    inlist = sys.argv[1:-7]
    biasname = str(sys.argv[-5])
    darkname = str(sys.argv[-4])
    flatname = str(sys.argv[-3])
    usedark = str(sys.argv[-6])
    usebias = str(sys.argv[-7])
    outdir = str(sys.argv[-2]) + '/'
    reportdir = str(sys.argv[-1]) + '/'

    print(inlist)
    print(len(inlist))
    print(type(inlist))

    reducer(inlist, biasname, darkname, flatname, usedark, usebias, outdir, reportdir)


if __name__ == '__main__':
    main()