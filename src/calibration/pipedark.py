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

def get_exptimes(outdir):
    """
    Get target names and exp times for matching darks
    """
    date = outdir.split("/")[-3]
    image_dir = "/SPECULOOSPipeline/Telescopes/Callisto/images/"+date+"/"

    target_names, frames, exposure_times = [], [], []
    for file in glob.glob(outdir):
        if fnmatch.fnmatch(file, "image"):
            with open(file) as f:
                frame = f.readline()
                frames.append(frame)

    for name in frames:
        with open_fits_file(image_dir+name) as hdulist:
            exposure_time = hdulist[0].header['exptime']
            target_name = hdulist[0].header['object']
        
        target_names.append(target_name)
        exposure_times.append(exposure_time)

    return target_names, exposure_times

def darkmaker(inlist, biasname, darkname, outdir, reportdir, gain, run, targ):
    biasname = outdir + biasname

    position = 0
    i = 1
    lines = []

    for dsorted_fn in glob.glob(outdir + 'dsorted*'):
        os.remove(dsorted_fn)

    # find all darks
    for line in open(inlist):
        lines.append(line)
        fname = outdir + 'dsorted' + "{0:03d}".format(position)
        f = open(fname, 'a')
        f.write(line)
        f.close()
        if i == 50:
            i = 0
            position += 1
        i += 1

    with open(outdir + 'removeindexlist.dat', 'w') as rem:
        for fn in glob.glob(outdir + ('dsorted*')):
            rem.write(fn + "\n")

    i = 1

    if len(lines) == 0:
        print("WARNING: No dark images!")
        darkc = 0
    else:
        # import master bias to correct dark images
        with open_fits_file(biasname) as hdulist:
            bias = hdulist[0].data

        for line in open(outdir + 'removeindexlist.dat'):
            datamatrix = []
            mastermatrix = []
            call = line.strip('\n')
            for line in open(call):
                line = line.strip()
                print(line)
                with open_fits_file(line) as hdulist:
                    overscan = extract_overscan(hdulist)
                    overscan = sigma_clip(overscan)
                    medoverscan = np.ma.median(overscan)
                    data = image_trim(hdulist)

                    # if np.shape(hdulist[0].data)[0] > 2048:
                    #     overscan = extract_overscan(hdulist)
                    #     overscan = sigma_clip(overscan)
                    #     medoverscan = np.ma.median(overscan)
                    #     #why is the shape of this 2048x2028?
                    #     #data = hdulist[0].data[0:2048,20:2068]
                    #     data = hdulist[0].data[22:2066,2:2048]
                    # elif np.shape(hdulist[0].data)[0] < 2000:
                    #     medoverscan = 0
                    #     data = hdulist[0].data
                    # else:
                    #     with open(outdir + "overscan.dat", "r") as f:
                    #         medoverscan = float(f.read())
                    #         # print("Extracted Overscan from median bias!: ",medoverscan)
                    #     data = hdulist[0].data[2:2046, 2:2048]

                    exposure = hdulist[0].header['exptime']
                        
                # subtract master bias image from the dark
                # sigma clip the overscan
                # print("Master Dark dimensions: ",np.shape(data))
                if np.shape(bias) == np.shape(data):
                    corrected = (data - medoverscan - bias) / exposure
                else:
                    print("WARNING: Bias and Dark dimensions do NOT match! Using UNCORRECTED dark images!")
                    print("Master Bias dimensions: ", np.shape(bias))
                    print("Master Dark dimensions: ", np.shape(data))
                    corrected = (data - medoverscan) / exposure
                # datamatrix is an array of the new corrected data with overscan and residual bias removed for every file
                datamatrix.append(corrected)
                # print np.shape(corrected)
            # shape of datamatrix: number of files x size of each file (2048x2048)
            # master = the median of the data over all files creating an image of median at every pixel
            # print np.shape(datamatrix)
            # print type(datamatrix)
            # clip_datamatrix = sigma_clip(datamatrix,iters=None,axis=0)
            # print np.shape(clip_datamatrix)
            master = np.ma.median(sigma_clip(datamatrix, axis=0), axis=0)
            # print i
            mastermatrix.append(master)
            i += 1

        print('averaging')
        # average all dark medians over all files to create a master dark image
        dark = np.float32(np.mean(mastermatrix, axis=0))

        phdu = pyfits.PrimaryHDU(dark)
        outname = outdir + darkname
        # command = 'rm -f '+outname
        # os.system(command)
        if os.path.exists(outname):
            os.remove(outname)

        phdu.writeto(outname)

        for dsorted_fn in glob.glob(outdir + 'dsorted*'):
            os.remove(dsorted_fn)

        # os.system('rm -f removeindexlist.dat '+outdir+'dsorted*')
        os.remove(outdir + 'removeindexlist.dat')

        # print np.shape(sigma_clip(dark))
        # print np.ma.median(sigma_clip(dark))
        # print np.shape(dark)
        # print gain
        # print type(gain)
        # print np.ma.median(dark)
        # print type(np.ma.median(dark))
        darkc = float(gain) * np.ma.median(sigma_clip(dark))

    # QC update dark current
    dirsplit = outdir.split("/")
    date = dirsplit[-3]
    field = targ
    tel = dirsplit[-5]
    # dirname = "/".join(dirsplit[:-4])
    QC([0, reportdir + "/QC", date, field, tel, run, 'QC7', darkc])

def matching_darkmaker(inlist, biasname, darkname, outdir, reportdir, gain, run):
    """
    Create a master dark from matching darks based on exposure time
    """
    target_list, exptimes = get_exptimes(outdir)

    # load master bias 
    with open_fits_file(os.path.join(outdir, biasname)) as hdul:
        bias = hdul[0].data

    for target, exp_time in zip(target_list, exptimes):
        for f in glob.glob(os.path.join(outdir, "dsorted*")):
            os.remove(f)

        position, i, lines = 0, 1, []
        for line in open(inlist):
            lines.append(line)
            fname = os.path.join(outdir, f"dsorted{position:03d}")
            with open(fname, "a") as f:
                f.write(line)
            if i == 50:
                i, position = 0, position + 1
            i += 1

        with open(os.path.join(outdir, "removeindexlist.dat"), "w") as rem:
            for fn in glob.glob(os.path.join(outdir, "dsorted*")):
                rem.write(fn + "\n")

        if not lines:
            print(f"WARNING: No darks found for {target}")
            darkc = 0
            continue

        mastermatrix = []

        for listfile in open(os.path.join(outdir, "removeindexlist.dat")):
            datamatrix = []
            for fname in open(listfile.strip()):
                fname = fname.strip()
                with open_fits_file(fname) as hdul:
                    data = image_trim(hdul)
                    overscan = extract_overscan(hdul)
                    medoverscan = np.ma.median(sigma_clip(overscan))
                    exposure = hdul[0].header["exptime"]

                if np.round(exposure) != np.round(exp_time):
                    continue

                if bias.shape == data.shape:
                    corrected = data - medoverscan - bias
                else:
                    print("WARNING: Bias and dark dimensions mismatch! Using uncorrected dark.")
                    corrected = (data - medoverscan) / exposure

                datamatrix.append(corrected)

            if np.ndim(datamatrix) == 3:
                master = np.ma.median(sigma_clip(datamatrix, axis=0), axis=0)
                mastermatrix.append(master)

        if not mastermatrix:
            print(f"WARNING: No matching darks for {target} ({exp_time}s)")
            darkc = 0
            continue

        dark = np.float32(np.mean(mastermatrix, axis=0))
        outname = os.path.join(outdir, f"{darkname}{target}")
        if os.path.exists(outname):
            os.remove(outname)
        pyfits.PrimaryHDU(dark).writeto(outname)

        for f in glob.glob(os.path.join(outdir, "dsorted*")):
            os.remove(f)
        os.remove(os.path.join(outdir, "removeindexlist.dat"))

        darkc = float(gain) * np.ma.median(sigma_clip(dark))

        dirsplit = outdir.split("/")
        date, tel = dirsplit[-3], dirsplit[-5]
        QC.main([0, f"{reportdir}/QC", date, target, tel, run, "QC7", darkc])


def main():
    inlist = str(sys.argv[1])
    biasname = str(sys.argv[2])
    darkname = str(sys.argv[3])
    outdir = str(sys.argv[4]) + '/'
    reportdir = str(sys.argv[5]) + '/'
    gain = sys.argv[6]
    run = sys.argv[7]
    targ = sys.argv[8]

    if pipeutils.detect_instrument(pyfits.open(biasname)) == 'spirit':
        print("Using matching darks")
        matching_darkmaker(inlist, biasname, darkname, outdir, reportdir, gain, run, targ)
    
    else:
        darkmaker(inlist, biasname, darkname, outdir, reportdir, gain, run, targ)


if __name__ == '__main__':
    main()
