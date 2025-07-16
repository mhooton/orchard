# -*- coding: utf-8 -*-

import tempfile
from photometry.catmatch import shift_wcs_axis, apply_correct, apply_correct_old, lmq_fit, calc_seps,\
    load_wcs_from_keywords, correct_catfile
import photometry.casutools as casutools
from photometry.vector_plot import wcsf_QCheck
from photometry.wcs_status import set_wcs_status
from multiprocessing import Pool as ThreadPool
from functools import partial
import photometry.casutools as casutools
from astropy.io import fits
import os
import numpy as np
import fitsio
from astropy.io import fits as pf
from collections import namedtuple
try:
    import pickle as pickle
except ImportError:
    import pickle
import json
import psutil
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)


Catalogue = namedtuple('Catalogue', ['cat_name', 'ra_lims', 'dec_lims'])


def initialise_wcs_cache(fname, outdir, catsrc, thresh, verbose, force=False):
    print("\n**Constructing initial wcs cache**\n")
    catalogue_name = outdir+'_initial-catalogue.fits'
    #imcore - source detect on this first image and output fits catalogue to this file 'catalogue_name'
    casutools.imcore(fname, catalogue_name, threshold=thresh, verbose=verbose)
    #wcsfit - matches up objects from image with astrometric positions standards. Input first image and imcore catalogue
    print(("fname = "+str(fname)))
    print(("catalogue_name = "+str(catalogue_name)))
    print(("catsrc = "+str(catsrc)))
    print(("verbose = "+str(verbose)))
    casutools.wcsfit(fname, catalogue_name, catsrc=catsrc, verbose=verbose)


def  m_solve_images(filelist, outfile,
                   nproc=None,
                   thresh=20.0,
                   verbose=False,
                   catsrc='vizgaia2',
                   rcore =4,
                    ipix=6,
                    ext='fits'): #catsrc is the source of standard stars to be used in astrometric fit
    # type: (object, object, object, object, object, object, object, object) -> object
    # viz2mass - 2MASS standards are extracted from a VizieR distribution
    # vizgaia2 - Gaia DR2 standards are extracted from a VizieR distribution
    # wcsref is the full path to the file if catsrc is defined as a local option

    infiles = []

    with open(filelist) as infile:
        for line in infile:
            image = line.strip()
            status_checks = ['ok', 'ok']


            if all(status == 'ok' for status in status_checks):
                infiles.append(image)

    # print infiles
    outdir = os.path.dirname(os.path.realpath(infiles[0]))
    print('\n**m_solve_images:**\n')
    # print outdir
    # pass second (first doesn't have all nec. info) processed image file to initialise wcs cache
    #  by performing imcore and wcsfit on it:
    initialise_wcs_cache(infiles[0], outdir, catsrc, thresh, verbose)

    #run casu_solve on the images
    print(ext)
    fn = partial(handle_errors_in_casu_solve,
                 thresh=thresh,
                 verbose=verbose,
                 catsrc=catsrc,
                 rcore=rcore,
                 ipix=ipix,
                 ext=ext)

    pool = ThreadPool(nproc)

    return pool.map(fn, infiles)


def handle_errors_in_casu_solve(casuin, *args, **kwargs):
    '''
  Catch any exceptions that may be thrown by the image solving routine to
  prevent the pipeline from crashing
  '''
    try:
        return_value = casu_solve(casuin, *args, **kwargs)
    except Exception as err:
        print(("Exception handled in `casu_solve`: {}".format(str(err))))
        set_wcs_status(casuin, succeeded=False)
    else:
        set_wcs_status(casuin, succeeded=True)
        return return_value



def casu_solve_old(casuin,
                   thresh=20,
                   verbose=False,
                   catsrc='viz2mass',
                   catpath=None):
    hdulist = fits.open(casuin)

    cen = [[dist_map['CRPIX1'], dist_map['CRPIX2']]]

    TEL_RA = hdulist['TEL_RA']
    TEL_DEC = hdulist['TEL_DEC']

    for key in dist_map:
        print((key, dist_map[key], hdulist.get(key)))

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    catfile_name = casuin.replace('.fits', '.cat')
    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose)

    cat_names = []
    RA_lims = []
    DEC_lims = []

    catpath = (catpath if catpath is not None
            else os.path.join(os.getcwd(), 'catcache'))
    for line in open(catpath + '/index'):
        vals = line.strip('\n').split(' ')
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]), float(vals[3])]]
        DEC_lims += [[float(vals[4]), float(vals[5])]]

    n = 0

    cat_name = cat_names[n]

    with pf.open(catpath + '/' + cat_name) as catd:
        catt = catd[1].data.copy()
    cat = {'ra': catt['ra'], 'dec': catt['dec'], 'Jmag': catt['Jmag']}

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    with fitsio.FITS(catfile_name) as mycatt:
        mycat = {'Aper_flux_3': mycatt[1]['Aper_flux_3'][:]}
        my_X = mycatt[1]['x_coordinate'][:]
        my_Y = mycatt[1]['y_coordinate'][:]

    try:
        dist_map = shift_wcs_axis(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y,
                                  TEL_RA, TEL_DEC,
                                  iters=10)
        dist_map = lmq_fit(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA,
                           TEL_DEC,
                           fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])
    except IOError:
        print("Performing initial fit")
        casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)
        dist_map = shift_wcs_axis(casuin, catfile_name, thresh=thresh, iters=30)
        dist_map = lmq_fit(dist_map, mycat, cat, RA_lims, DEC_lims, my_X, my_Y, TEL_RA,
                           TEL_DEC,
                           fitlist=['RA_s', 'DEC_s', 'CD1_1', 'CD2_2', 'CD1_2', 'CD2_1'])

    apply_correct_old(dist_map, casuin, TEL_RA, TEL_DEC)

    # wcs keywords may have changed since imcore was done, so we have to update the RA and DEC values.
    correct_catfile(catfile_name, casuin, nstars=2000)

    # Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, catpath=wcsref, verbose=verbose)

    # Do QC checks. should really break this out.

    plot = True
    wcsf_QCheck(mycat, casuin, os.path.basename(casuin).strip('.fits') + '.png', cat,
                RA_lims, DEC_lims, my_X, my_Y,
                plot=plot)

    return 'ok'


def casu_solve(casuin, thresh=2, verbose=False, catsrc='vizgaia2', rcore=4, ipix=6, ext='fits'):
    # ADD debug output at the start

    catpath = os.path.join(os.getcwd(), 'catcache')

    # give the catalogue the same name as the image file it's been run on
    catfile_name = casuin.replace('.' + ext, '.cat')
    # we need an input cat for each image in order to wcsfit each image
    print(casuin)
    print(catfile_name)

    casutools.imcore(casuin, catfile_name, threshold=thresh, verbose=verbose,
                     ipix=ipix, rcore=rcore)  # changed ipix from 2 to 6, rcore was 3 - changed to match appsize???

    # Now we're ready to solve wcs
    casutools.wcsfit(casuin, catfile_name, catsrc=catsrc, verbose=verbose)

    # Come back to this testing section later - perhaps break out into separate script?
    # find frame limits from index file in catcache directory

    catalogue = compute_frame_limits(catpath)

    if catalogue is None:
        raise ValueError(f"compute_frame_limits returned None for catpath: {catpath}")

    cat = reference_catalogue_objects(catalogue, catpath, catsrc)

    if cat is None:
        raise ValueError(f"reference_catalogue_objects returned None")

    # ADD error checking for catalog file
    if not os.path.exists(catfile_name):
        raise FileNotFoundError(f"Catalog file not found: {catfile_name}")

    with fits.open(catfile_name) as mycatt:
        # ADD error checking for data extension
        if len(mycatt) < 2:
            raise ValueError(f"Catalog file {catfile_name} has no data extension")

        # extract stacked image catalogue and extra X, Y coordinates and flux of each source
        mycatt_data = mycatt[1].data

        # ADD null data check
        if mycatt_data is None:
            raise ValueError(f"No data in catalog file {catfile_name}")

        # ADD column existence checks
        required_columns = ['Aper_flux_3', 'X_coordinate', 'Y_coordinate', 'Sequence_number']
        missing_columns = [col for col in required_columns if col not in mycatt_data.columns.names]
        if missing_columns:
            raise ValueError(f"Missing columns in catalog: {missing_columns}")

        # print mycatt[1].header
        mycat = {'Aper_flux_3': mycatt_data['Aper_flux_3']}

        my_X = mycatt_data['X_coordinate']
        my_Y = mycatt_data['Y_coordinate']
        my_ID = mycatt_data['Sequence_number']

    # Do QC checks. should really break this out.
    # plot separation of each source from the ref catalogue as a quality check
    # FIXED: removed catsrc parameter and fixed filename extension
    # Remove the duplicate try/except block and fix this call:
    wcsf_QCheck(mycat, catfile_name, catsrc, casuin,
                os.path.basename(casuin).replace('.' + ext, '') + '.png',
                cat, catalogue.ra_lims[0], catalogue.dec_lims[0],
                my_X, my_Y, my_ID, plot=True)

    return 'ok'


class InitialAstrometricSolutionParser(object):
    KEYS = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2',
            'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
            'PV2_1', 'PV2_3', 'PV2_5', 'PV2_7']

    def __init__(self, filename):
        self.filename = filename

    @property
    def mimetype(self):
        fileformat_map = {
            '.fts': 'fits',
            '.json': 'json',
            '.p': 'pickle',
        }
        for key in fileformat_map:
            if self.filename.endswith(key):
                return fileformat_map[key]

        raise ValueError("Invalid file extension. Must be one of %s" % (
            list(fileformat_map.keys())))

    def parse(self):
        fn = getattr(self, 'parse_{mimetype}'.format(mimetype=self.mimetype))
        solution = fn()
        self.check_solution(solution)
        return solution

    def parse_json(self):
        with open(self.filename) as infile:
            json_data = json.load(infile)['wcs']
        return { key: json_data[key] for key in self.KEYS }

    def parse_fits(self):
        header = fits.getheader(self.filename)
        return { key: header[key] for key in self.KEYS }

    def parse_pickle(self):
        with open(self.filename) as infile:
            solution = pickle.load(infile)
        return { key: solution[key] for key in self.KEYS }

    @staticmethod
    def check_solution(solution):
        if 'CD1_1' not in solution:
            raise KeyError("Cannot find valid wcs solution in map {}".format(
                json.dumps(solution)))


def extract_dist_map(filename):
    return InitialAstrometricSolutionParser(filename).parse()


def reference_catalogue_objects(catalogue, catpath, catsrc):
    cat_name = os.path.join(catpath, catalogue.cat_name)
    # catfitsio = fitsio.FITS(cat_name)
    # print catfitsio
    # print catfitsio[1]
    # print catfitsio[1].read(rows=[:], columns=['ra','dec','__Gmag_','Source'])
    # ra = catfitsio[1]['ra'][:]
    # dec = catfitsio[1]['dec'][:]
    # gmag = catfitsio[1]['__Gmag_'][:]
    # source = catfitsio[1]['Source'][:]
    with fits.open(cat_name) as catd:
        catt = catd[1].data.copy()

    if catsrc == 'vizgaia2':
        return {'ra': catt['ra'], 'dec': catt['dec'], 'Gmag': catt['Gmag']}
    elif catsrc == 'viz2mass':
        return {'ra': catt['ra'], 'dec': catt['dec'], 'Jmag': catt['Jmag'], 'Kmag':catt['Kmag'],'Hmag':catt['Hmag']}




def compute_frame_limits(catpath):
    index_filename = os.path.join(catpath, 'index')
    cat_names = []
    RA_lims = []
    DEC_lims = []

    for line in open(index_filename):
        vals = line.strip().split()
        cat_names += [vals[0]]
        RA_lims += [[float(vals[2]), float(vals[3])]]
        DEC_lims += [[float(vals[4]), float(vals[5])]]

    return Catalogue(cat_names[0], RA_lims, DEC_lims)
