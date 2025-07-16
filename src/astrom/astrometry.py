# This module makes a proper WCS-header of a FITS file.

from astropy.io import fits
import os
from multiprocessing import Pool as ThreadPool
from functools import partial
import argparse
import timeit
import astrom.find_pa
from calibration.pipeutils import detect_instrument
import fnmatch
import shutil
from astrom.twirl_speculoos import twirl_wcs

def main(args):
    import sys
    sys.stdout.flush()
    print("*************************", flush=True)
    print("Starting astrom.py", flush=True)
    print("*************************", flush=True)
    start_time = timeit.default_timer()
    print(f"Input file: {args.filelist}", flush=True)
    print(f"Extension: {args.ext}", flush=True)

    # Before the try block
    infiles = []
    print("About to process filelist...", flush=True)

    try:
        with open(args.filelist) as filel:
            # print("File opened successfully", flush=True)
            for line in filel:
                f = line.strip()
                # print(f"Found file in list: {f}", flush=True)
                print(f)
                hdulist = fits.open(f)
                # print(f"FITS file opened: {f}", flush=True)
                if 'CTYPE1' not in hdulist[0].header:
                    infiles.append(f)
                    # print(f"Added to infiles: {f}", flush=True)
                hdulist.close()  # Add this to close the file properly
    except Exception as e:
        print("Importing Astrometry failed: " + str(e), flush=True)

    print(f"Number of files in infiles: {len(infiles)}", flush=True)
    # infiles = infiles[:20]
    # try:
    #     pool = ThreadPool(int(args.nproc))
    #     fn = partial(add_astrometry,ext=args.ext)
    #     pool.map(fn, infiles)
    #     # for infile in infiles:
    #     #     add_astrometry(infile, args.ext)
    #
    # except Exception as e:
    #     print("Astrometry failed: " + str(e))

    for infile in infiles:
        print(f"Processing: {infile}")
        add_astrometry(infile, args.ext)
        print(f"Completed: {infile}")

    keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

    try:
        for f in infiles:
            new_f = f.replace(args.ext, 'new')
            f1 = f.replace(args.ext, 'axy')
            f2 = f.replace(args.ext, 'wcs')
            f3 = f.replace(args.ext, 'solved')
            f4 = f.replace(args.ext, 'corr')
            f5 = f.replace("." + args.ext, '-indx.xyls')
            f6 = f.replace(args.ext, 'match')
            f7 = f.replace(args.ext, 'rdls')
            if os.path.isfile(new_f):

                hdulist_old = fits.open(f)
                hdr_old = hdulist_old[0].header
                data = hdulist_old[0].data
                hdulist_new = fits.open(new_f)
                hdr_new  = hdulist_new[0].header
                pa, jd, solve = find_pa.pa(new_f)

                hdr_old.set('PA',pa)
                # hdr_old.set('CTYPE1',hdr_new['CTYPE1']) #'RA---TAN')
                # hdr_old.set('CTYPE2', hdr_new['CTYPE2']) #'DEC--TAN')
                # hdr_old.set('CRVAL1',hdr_new['CRVAL1'])
                # hdr_old.set('CRVAL2',hdr_new['CRVAL2'])
                # hdr_old.set('CRPIX1',hdr_new['CRPIX1'])
                # hdr_old.set('CRPIX2',hdr_new['CRPIX2'])
                # hdr_old.set('CD1_1',hdr_new['CD1_1'])
                # hdr_old.set('CD2_1',hdr_new['CD2_1'])
                # hdr_old.set('CD1_2',hdr_new['CD1_2'])
                # hdr_old.set('CD2_2',hdr_new['CD2_2'])
                for k in keywords:
                    if hdr_old.get(k) == None:
                        hdr_old.set(k,hdr_new[k])
                # if hdr_old.get('RA') == None:
                #     hdr_old.set('RA', '21 31 48.113')
                #     hdr_old.set('DEC', '48 26 29.330')
                hdr_old.set('ASTSOLVE','T')
                fits.writeto(f,data,hdr_old,overwrite=True)

                for fname in [new_f,f1,f2,f3,f4,f5,f6,f7]:
                    os.remove(fname)

        elapsed = timeit.default_timer() - start_time
        print('Total time taken for astrom: ' + str(elapsed/60.) + ' minutes')

    except Exception as e:
        print("Replacing Astrometry failed: " + str(e))
        print("on file: " + f)


def add_astrometry(f, ext):
    if fnmatch.fnmatch(f, '*.' + ext):
        file = f
        print("Add astrom.net data to header of file: " + file)

        # # Find solve-field command
        # solve_field_cmd = shutil.which('solve-field')
        # if not solve_field_cmd:
        #     raise RuntimeError("solve-field not found in PATH")

        twirl_wcs(str(file),verbose=True)

        #
        # hdulist = fits.open(file)
        # prihdr = hdulist[0].header
        #
        # # Extract RA and DEC
        # ra_val = prihdr.get('RA')
        # dec_val = prihdr.get('DEC')
        #
        # if ra_val is not None and dec_val is not None:
        #     # Process RA
        #     if " " in str(ra_val):
        #         strra = str(ra_val.replace(" ", ":"))
        #     else:
        #         strra = str(ra_val)
        #
        #     # Process DEC
        #     if " " in str(dec_val):
        #         strdec = str(dec_val.replace(" ", ":"))
        #     else:
        #         strdec = str(dec_val)
        #
        #     if detect_instrument(hdulist) == 'spirit':
        #         print("Running solve-field with spirit astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/spirit_astrometry.cfg " \
        #                      f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                      f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                      f"--scale-low 0.32 --scale-high 0.38 --scale-units arcsecperpix {file}"
        #     elif detect_instrument(hdulist) == 'andor':
        #         print("Running solve-field with andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"
        #     else:
        #         print("Instrument not spirit or andor, defaulting to andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"

        # print(f"Running: {cmd}")
        # os.system(cmd)
        # else:
        #     print(f"Warning: Missing RA and/or DEC in {file}, skipping astrom")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filelist')
    parser.add_argument('-n','--nproc')
    parser.add_argument('-e', '--ext')
    main(parser.parse_args())
