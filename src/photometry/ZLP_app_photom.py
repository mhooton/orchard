#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Aperture Photometry Execution

Orchestrates forced aperture photometry across all science frames using reference
catalogues from stacked images. Applies proper motion corrections to catalogue
positions and performs quality assessment of photometric measurements.

Args:
    catfile (str): Path to stack catalogue FITS file with Gaia cross-match
    filter (str): Filter name (e.g., 'I+z', 'zYJ')
    filelist (str): Path to file containing list of science frame paths
    outdir (str): Output directory for photometry files
    date (str): Observation date in YYYYMMDD format
    nproc (int, optional): Number of parallel processes (default: 1)
    apsize (float, optional): Aperture radius in pixels (default: 4)
    ext (str, optional): File extension 'fits' or 'fts' (default: 'fits')
    ipix (float, optional): Minimum pixels for detection (default: 6)
    s_thresh (float, optional): Detection threshold in sigma (default: 20)
    catsrc (str, optional): Catalogue for WCS solving (default: 'vizgaia3')
    norunwcs (bool, optional): Skip WCS refinement (default: False)
    verbose (bool, optional): Print detailed output (default: False)

Processing workflow:
    1. Optional WCS refinement of science frames (if norunwcs=False)
    2. Proper motion correction from Gaia epoch (J2016.0) to observation date
    3. Forced aperture photometry via casutools.imcore_list
    4. Image quality assessment (FWHM, seeing, cloud status, frame shifts)
    5. Validation and removal of failed frames

Output:
    - {image}.phot files for each science frame (FITS tables with aperture fluxes)
    - {filelist}_phot file listing successfully processed frames
    - Quality metrics in FITS headers (SEEING, FWHM, PSF_a_*, PSF_b_*, etc.)

Example:
    python ZLP_app_photom.py --catfile stack_catalogue_I+z.fits --filter I+z \\
        --filelist I+z_processed.dat --outdir /output/target/ \\
        --date 20240115 --nproc 8 --apsize 4

Notes:
    - Proper motion correction is critical for multi-year baselines
    - Aperture size should typically match FWHM of images
    - Serial mode (nproc=1) recommended for debugging
"""

import sys
import linecache
import threading
from os.path import isfile, join
from photometry.wcs_fitting import m_solve_images
from photometry.wcs_photom import m_wcs_photom
import json
import argparse
from photometry.wcs_status import wcs_succeeded
from astrom.proper_motion import generate_new_cat
import timeit
import warnings
warnings.filterwarnings("ignore", message=".*the RADECSYS keyword is deprecated, use RADESYSa*",
                        module="astropy.wcs.wcs")
warnings.filterwarnings("ignore", message=".*invalid value encountered in log10*", category=RuntimeWarning)

def main(argv):
    # if you don't provide an outlist name i'll assume you just want to add _phot to the end
    if not argv.outlist:
        argv.outlist = argv.filelist + '_phot'

    outfile = open(argv.outlist, 'w')
    outfile.close()

    # in our case we don't have a distortion map
    # dist_map = extract_dist_map(argv.dist)
    dist_map = None

    if not argv.norunwcs and argv.norunwcs != 'zYJ':
        beforewcs = timeit.default_timer()
        # SOLVE ASTROMETRY:
        # apply wcs solution to each image using the wcs reference frame from generate_reference_frame
        m_solve_images(argv.filelist, argv.outlist,
                       nproc=argv.nproc,
                       thresh=argv.s_thresh,
                       verbose=argv.verbose,
                       catsrc=argv.catsrc,
                       rcore=argv.apsize,
                       ipix=argv.ipix,
                       ext=argv.ext)
        afterwcs = timeit.default_timer()
        print("WCS SOLVING TOOK: " + str((afterwcs - beforewcs) / 60.))

    # CORRECT FOR PROPER MOTION
    gaia_date = '20160101'  # DR3 epoch is J2016.0

    # Create the PM catalogue path in the target output directory
    import os
    pm_filename = os.path.basename(argv.catfile).replace("stack_catalogue", "stack_catalogue_pm")
    cat_pm = os.path.join(argv.outdir, pm_filename)

    # Generate the PM catalogue directly in the desired location
    generate_new_cat(argv.catfile, gaia_date, cat_pm)

    # PERFORM APERTURE PHOTOMETRY:
    beforephotom = timeit.default_timer()
    m_wcs_photom(argv.filelist, argv.outlist, argv.apsize, cat_pm,
                 nproc=argv.nproc,
                 verbose=argv.verbose)
    afterphotom = timeit.default_timer()
    print("PHOTOMETRY TOOK: " + str((afterphotom - beforephotom) / 60.))

    check_wcs_succeed(argv.outlist)

    print("\n**Aperture Photometry Complete**\n")

def check_wcs_succeed(outlist):
    line_write = []

    with open(outlist) as f:
        for line in f:
            image = line.strip()
            if wcs_succeeded(image) == True:
                line_write.append(image)
            else:
                print("Removing file " + str(image))

    f = open(outlist, 'w')
    f.write(".phot\n".join(line_write) + ".phot")
    # with open(outlist,'w') as f:
    #     for line in line_write:
    #         f.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", '--catfile', help='Input catalogue', required=True)
    parser.add_argument('--filter', help='Filter', required=True)
    parser.add_argument('-f', '--filelist', help='List of files', required=True)
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('-d', '--date', help='Date of Observation', required=True)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--outlist', help="List of completed files")
    parser.add_argument('--nproc', type=int, default=1, help='Number of processors')
    parser.add_argument('--apsize', type=float, default=4, help='Aperture size')
    parser.add_argument('--ext', default='fits', help='File Extension')
    parser.add_argument('--ipix', type=float, default=6, help='Minimum number of pixels above background')
    parser.add_argument('--s_thresh',
                        type=float,
                        default=20,
                        help='Detection threshold in sigma')
    parser.add_argument('--catsrc',
                        default='vizgaia3',
                        help='Catalogue for wcs solving')
    parser.add_argument('--catpath',
                        default='catcache',
                        help='Local casutools catalogue cache')
    parser.add_argument('--norunwcs',
                        action='store_true',
                        default=False,
                        help='Do not astrometrically solve the images')

    main(parser.parse_args())