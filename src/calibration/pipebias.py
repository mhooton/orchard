import math
import sys
import os
import numpy as np
import glob
from calibration.pipeutils import extract_overscan, image_trim, open_fits_file
from astropy.io import fits
from astropy.stats import sigma_clip
from reporting.QC import main as QC

def biasmaker(inlist, biasname, outdir, reportdir, gain, run, targ):
    """Create master bias from list of bias images."""
    os.system('mkdir ' + outdir)
    position = 0
    i = 1
    lines = []

    for bsorted_fn in glob.glob(outdir + 'bsorted*'):
        os.remove(bsorted_fn)

    # loop through bias fits filenames and write each into a file called 'bsorted000' - I believe this is necessary to combine multiple bias frames together from different nights?
    # (or when there are more than 50 bias frames)
    for line in open(inlist):  # Changed from file() to open() for Python 3
        print(line)
        lines.append(line)
        fname = outdir + 'bsorted' + "{0:03d}".format(position)
        f = open(fname, 'a')
        f.write(line)
        f.close()
        if i == 50:
            i = 0
            position += 1
        i += 1

    # export the bsorted files into a new file 'removeindexlist.dat' - this is only necessary if there are only multiple files if there are multiple objects
    # again not sure this is necessary
    # os.system('ls '+outdir+'bsorted* >removeindexlist.dat')
    with open(outdir + 'removeindexlist.dat', 'w') as rem:
        for fn in glob.glob(outdir + ('bsorted*')):
            rem.write(fn + '\n')

    # start counter
    i = 0

    if len(lines) == 0:
        print("WARNING: No bias images!")
        ron = 0
        overscan = 0
    else:
        # loop through the bsorted filenames in removeindexlist.dat
        for line in open(outdir + 'removeindexlist.dat'):  # Changed from file() to open()

            datamatrix = []
            mastermatrix = []
            # removes '\n' from beginning and end of string
            call = line.strip('\n')
            # for the bias list file on this line, loop through each line of this file (call)
            for line_call in open(call):  # Changed from file() to open()
                # remove whitespace
                line_call = line_call.strip()
                # open the fits image file referred to in the line:
                with open_fits_file(line_call) as hdulist:
                    # obtain the overscan region
                    # print(np.shape(hdulist[0].data),len(hdulist[0].data))
                    overscan = extract_overscan(hdulist)
                    # print(np.shape(overscan))
                    # remove the pre/overscans from the y-axis (20 pixels top and bottom):
                    # hdulist[0].data is in format 2088 x 2048 therefore y by x, so is this correct?
                    # data = hdulist[0].data[0:2048,20:2068]
                    data = image_trim(hdulist)

                    # if np.shape(hdulist[0].data)[0] > 2048:
                    #     data = hdulist[0].data[22:2066,2:2048]
                    # elif np.shape(hdulist[0].data)[0] < 2000:
                    #     data = hdulist[0].data
                    # else:
                    #     data = hdulist[0].data[2:2046, 2:2048]
                    # print(np.shape(data))
                # sigma clip the overscan and
                # subtract the median of the overscan array from the bias to give the residual bias
                corrected = data - np.ma.median(sigma_clip(overscan, maxiters=None))
                datamatrix.append(corrected)

                if i == 0:
                    b0 = corrected
                if i == 1:
                    b1 = corrected
                i += 1

            print('medianing')
            # print np.shape(datamatrix)
            # find the median of the overscan corrected data and use this as the master
            # and sigma clip the datamatrix
            master = np.ma.median(sigma_clip(datamatrix, maxiters=None, axis=0), axis=0)
            mastermatrix.append(master)

        print('averaging')
        # sigma clip!!!
        # print np.shape(mastermatrix)
        bias = np.float32(np.mean(mastermatrix, axis=0))

        phdu = fits.PrimaryHDU(bias)
        outname = outdir + biasname
        # command = 'rm -f '+outname
        # os.system(command)
        if os.path.exists(outname):
            os.remove(outname)

        phdu.writeto(outname)

        for bsorted_fn in glob.glob(outdir + 'bsorted*'):
            os.remove(bsorted_fn)

        #    os.system('rm -f '+outdir+'bsorted* removeindexlist.dat')
        os.remove(outdir + 'removeindexlist.dat')

        adu_diff = np.std(b1 - b0)
        ron = adu_diff * float(gain) / (2. ** 0.5)

        with open(outdir + "readoutnoise.dat", "w") as f:
            f.write('%s' % ron)
            f.close()

        with open(outdir + "overscan.dat", "w") as f:
            f.write('%s' % np.ma.median(sigma_clip(overscan, maxiters=None)))
            f.close()

    # QC update RON
    dirsplit = outdir.split("/")
    date = dirsplit[-3]
    field = targ
    tel = dirsplit[-5]
    dirname = "/".join(dirsplit[:-4])
    QC([0, reportdir + "QC", date, field, tel, run, 'QC6', ron])
    QC([0, reportdir + "QC", date, field, tel, run, 'QC8', np.median(overscan)])


def main():
    """Main function to handle command line arguments."""
    # inlist takes in the filenames of bias fits images, biasname is the name of the file we will output the master bias to and
    # outdir is the directory to output this to
    inlist = str(sys.argv[1])
    biasname = str(sys.argv[2])
    outdir = str(sys.argv[3]) + '/'
    reportdir = str(sys.argv[4]) + '/'
    gain = sys.argv[5]
    run = sys.argv[6]
    targ = sys.argv[7]
    # version = sys.argv[8]

    biasmaker(inlist, biasname, outdir, reportdir, gain, run, targ)


if __name__ == '__main__':
    main()