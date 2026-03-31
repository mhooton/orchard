#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
condense_photometry.py - Stage T10

Condenses per-image aperture photometry files into a single FITS output
cube for each identified target in the field of view.

Target identification is NOT performed here. Targets are read directly
from the TARGET_ROLE column in the Gaia_Crossmatch extension of the
stack catalogue, which was written by catalogue_fov.py (Stage T8).
One output file is produced per identified target (primary or secondary).

Processing workflow:
    1. Sort input .phot files by JD
    2. Filter out files with missing WCS headers
    3. Read identified targets from stack catalogue TARGET_ROLE column
    4. Build imagelist, catalogue, and flux HDUs from phot files
    5. For each identified target: write per-night output file and
       update (or create) the global multi-night output file
"""

import argparse
import logging
from astropy.io import fits
from contextlib import contextmanager
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
import sys
import collections
from astropy.table import Table
from shutil import copyfile
import fitsio
import os
from astropy.io import ascii as astropy_ascii

logging.basicConfig(level='INFO', format='%(levelname)7s %(message)s')
logger = logging.getLogger(__name__)


# this script takes all the phot files in and outputs one fits file cube


class SourceFile(object):
    '''
    Container for catalogue file
    '''

    def __init__(self, hdulist):
        self.hdulist = hdulist
        self.hdu = self.hdulist[1] #['apm-binarytable']
        self.data =  self.hdu.data
        self.header = self.hdu.header
        self.header0 = self.hdulist[0].header

    @classmethod
    def open_file(cls, filename):
        with fits.open(filename) as hdulist:
            return cls(hdulist)

    def __len__(self):
        return len(self.data)

    def get_aperture_sizes(self, indices):
        out = {}
        for index in indices:
            flux_key = 'ttype{i}'.format(i=2 * (index-1) + 20)
            error_key = 'ttype{i}'.format(i=2 * (index-1) + 21)
            out[index] = {
                'flux': self.header.comments[flux_key],
                'error': self.header.comments[error_key],
            }

        return out


class Image(SourceFile):
    '''
    Container for source image
    '''

    def __init__(self, hdulist):
        self.hdulist = hdulist
        self.header = hdulist[1].header


@contextmanager
def create_output_file(filename):
    outfile = fits.HDUList()
    phdu = fits.PrimaryHDU()
    outfile.append(phdu)
    yield outfile
    outfile.writeto(filename, overwrite=True)


class FITSImage(object):

    def __init__(self, name, nimages, napertures, data_type=np.float64):
        self.name = name
        self.data = np.zeros((napertures, nimages), dtype=data_type)

    def set_data(self, i, data):
        self.data[:, i] = data

    def hdu(self):
        return fits.ImageHDU(self.data, name=self.name)

def get_jd(fname):
    # to print out all headers:
    #print(fits.getheader(fname))
    return fname, fits.getheader(fname)['jd-obs']  #fits.getheader(fname)['mjd']

def sort_by_jd(files):
    pool = ThreadPool()
    mapping = pool.map(get_jd, files)
    return [row[0] for row in sorted(mapping, key=lambda row: row[1])]


def check_wcs_headers(source_header, filename):
    """Check if all required WCS headers are present"""
    # Critical WCS headers required for valid photometry
    required_wcs_headers = {'TC3_3', 'TC3_5', 'TC5_3', 'TC5_5', 'TCTYP3', 'TCTYP5', 'TCRPX3', 'TCRPX5', 'TCRVL3',
                            'TCRVL5'}

    missing_headers = []
    for header in required_wcs_headers:
        if header not in source_header:
            missing_headers.append(header)

    if missing_headers:
        logger.warning(f"File {filename} missing WCS headers: {missing_headers}")
        logger.warning(f"Excluding {filename} from processing (likely wcsfit failure)")
        return False
    return True

def main(args):
    if args.verbose:
        logger.setLevel('DEBUG')
    logger.debug(args)

    # print(args.filename)

    # with fits.open(args.filename) as hdulist:
    #     print(hdulist.info())

    # logger.info('Sorting images by mjd')
    # we only have JD stored in headers
    logger.info('Sorting images by jd')
    sorted_images = sort_by_jd(args.filename)

    # Filter out images with missing WCS headers
    logger.info('Checking for required WCS headers')
    valid_images = []
    excluded_images = []

    for filename in sorted_images:
        try:
            with fits.open(filename) as hdulist:
                if check_wcs_headers(hdulist[1].header, filename):
                    valid_images.append(filename)
                else:
                    excluded_images.append(filename)
        except Exception as e:
            logger.warning(f"Error reading {filename}: {e}")
            logger.warning(f"Excluding {filename} from processing")
            excluded_images.append(filename)

    if excluded_images:
        logger.warning(f"Excluded {len(excluded_images)} files due to missing WCS headers:")
        for excluded_file in excluded_images:
            logger.warning(f"  - {excluded_file}")

    if not valid_images:
        logger.error("No valid images remaining after WCS header check")
        return

    sorted_images = valid_images
    # number of images = nimages
    nimages = len(sorted_images)
    print(nimages)
    print(sorted_images[0])
    first = SourceFile.open_file(sorted_images[0])

    napertures = len(first.data)

    logger.info('{:d} images, {:d} apertures'.format(nimages, napertures))

    logger.debug('Allocating memory for catalogue')
    # Construct an ndarray that allows field access using attributes
    #('GMAG', np.float64),('GAIA_ID', '19a'),
    # cat = args.outputdir + '/StackImages/' + args.target + '_stack_catalogue_' + args.filter + "." + args.ext
    with fitsio.FITS(args.cat, 'rw') as catfile:
        gaia = catfile['Gaia_Crossmatch']
        pmra = gaia['pmra'].read()
        pmdec = gaia['pmdec'].read()
        g_rp = gaia['g_rp'].read()
        bp_rp = gaia['bp_rp'].read()
        gaia_dr2_id = gaia['gaia_dr2_id'].read()  # Keep using DR2 ID for compatibility
        gaia_dr3_id = gaia['gaia_dr3_id'].read()  # Also read DR3 ID
        parallax = gaia['parallax'].read()
        teff = gaia['teff'].read()
        gmag = gaia['gmag'].read()

    obj_ids = np.array(['SP{:06d}'.format(i) for i in np.arange(napertures)])
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(args.outputdir)))
    print(basedir)

    catalogue_data = np.recarray(napertures,
                                 dtype=[('OBJ_ID', '26a'),
                                        ('RA', np.float64),
                                        ('DEC', np.float64),
                                        ('FLUX_MEAN_' + str(args.date), np.float64),
                                        ('FLUX_MEDIAN_' + str(args.date), np.float64),
                                        ('NPTS', np.int64),
                                        ('PMRA', np.float64),
                                        ('PMDEC', np.float64),
                                        ('GAIA_DR2_ID', 'a19'),
                                        ('GAIA_DR3_ID', 'a19'),
                                        ('G_RP', np.float64),
                                        ('BP_RP', np.float64),
                                        ('PARALLAX', np.float64),
                                        ('TEFF', np.float64),
                                        ('GMAG', np.float64)])

    logger.debug('Filling catalogue data')
    catalogue_data['OBJ_ID'] = obj_ids
        #np.array(['SP{:06d}'.format(i) for i in np.arange(napertures)])
    catalogue_data['NPTS'] = np.ones(napertures) * nimages
    catalogue_data['RA'] = np.degrees(first.data['ra'])
    catalogue_data['DEC'] = np.degrees(first.data['dec'])
    catalogue_data['PMRA'] = pmra
    catalogue_data['PMDEC'] = pmdec
    catalogue_data['GAIA_DR2_ID'] = gaia_dr2_id
    catalogue_data['GAIA_DR3_ID'] = gaia_dr3_id
    catalogue_data['G_RP'] = g_rp
    catalogue_data['BP_RP'] = bp_rp
    catalogue_data['GMAG'] = gmag
    catalogue_data['TEFF'] = teff
    catalogue_data['PARALLAX'] = parallax

    # crossmatch = 'Gaia'
    # if crossmatch == '2MASS':
    #     jmag,kmag,hmag = [],[],[]
    #     for line in file(args.outputdir+'/StackImages/'+args.target+'_stack_catalogue_'+args.filter+'_2MASS.dat'):
    #         line = line.strip().split("\t")
    #         jmag.append(line[0])
    #         kmag.append(line[1])
    #         hmag.append(line[2])
    #     catalogue_data['JMAG'] = np.array(jmag)
    #     catalogue_data['KMAG'] = np.array(kmag)
    #     catalogue_data['HMAG'] = np.array(hmag)
    # elif crossmatch == 'Gaia':
    #     gmag, gid = [], []
    #     gaia_file = args.outputdir+'/StackImages/'+args.target+'_stack_catalogue_'+args.filter+'_Gaia.dat'
    #     if os.path.isfile(gaia_file):
    #         for line in file(gaia_file):
    #             line = line.strip().split("\t")
    #             gmag.append(line[0])
    #             gid.append(line[1])
    #         catalogue_data['GMAG'] = np.array(gmag)
    #         catalogue_data['GAIA_ID'] = np.array(gid)

    logger.debug('Allocating memory for imagelist')
    imagelist_from_header_data_dtype = [('airmass', np.float32),
                                        ('altitude', np.float32),
                                        ('ambtemp', np.int64),
                                        ('azimuth', np.float32),
                                        ('bjd-obs', np.float64),
                                        ('ccd-temp', np.float32),
                                        ('crowded', np.int32),
                                        ('darkcur',np.float32),
                                        ('dec', 'a16'),
                                        ('dec_move', np.float32),
                                        ('dewpoint', np.int64),
                                        ('exptime', np.float32),
                                        ('filter', 'a4'),
                                        ('filtfwhm', np.float32),
                                        ('focallen', np.float32),
                                        ('focuspos', np.int64),
                                        ('focustem', np.float32),
                                        ('frame_sn', np.float32),
                                        ('fwhm', np.float32),
                                        ('gain', np.float32),
                                        ('ha', 'a16'),
                                        ('hjd-obs', np.float64),
                                        ('humidity', np.float32),
                                        ('imagetyp', 'a16'),
                                        ('jd-obs', np.float64),
                                        ('minpix', np.float32),
                                        ('nbsize', np.int32),
                                        ('numbrms', np.int64),
                                        ('object', 'a20'),
                                        ('pa', np.float32),
                                        ('ra', 'a16'),
                                        ('ra_move', np.float32),
                                        ('rapos', np.float32),
                                        ('rcore', np.float32),
                                        ('ron',np.float32),
                                        ('seeing', np.float32),
                                        ('sky_move', np.float32),
                                        ('skylevel', np.float32),
                                        ('skynoise', np.float32),
                                        ('skytemp', np.float32),
                                        ('stdcrms', np.float32),
                                        ('telescop', 'a20'),
                                        ('threshol', np.float32),
                                        ('wcscompl', bool),
                                        ('wcsf_dec', np.float32),
                                        ('wcsf_ns', np.int64),
                                        ('wcsf_ra', np.float32),
                                        ('wcsf_rms', np.float32),
                                        ('windspd', np.int64),
                                        ('xbinning', np.float32),
                                        ('ybinning', np.float32),
                                        ]

    # Add the psf measurements
    for i in range(1, 10):
        for key in ['a', 'b', 't']:
            imagelist_from_header_data_dtype.append(
                    ('psf_{key}_{i}'.format(key=key, i=i), np.float32)
                    )

    # # Add the psf measurements - ONLY N=1 components for now
    # for key in ['a', 'b', 't']:
    #     imagelist_from_header_data_dtype.append(
    #         ('psf_{key}_1'.format(key=key), np.float32)
    #     )
    # # TODO: Add PSF_*_2 through PSF_*_9 when twirl is working properly

    imagelist_extra_dtype = [('cd1_1', np.float64),
                             ('cd1_2', np.float64),
                             ('cd2_1', np.float64),
                             ('cd2_2', np.float64),
                             ('crpix1', np.float64),
                             ('crpix2', np.float64),
                             ('crval1', np.float64),
                             ('crval2', np.float64),
                             ('ctype1', 'a8'),
                             ('ctype2', 'a8'),
                             ('hicount', np.int64),
                             ('locount', np.int64),
                             ('tmid', np.float64),
                             ('nsources', np.int64)]

    # Don't error if these keys are missing
    optional_keys = {
        'frame_sn', 'decpos', 'dec_move', 'dec_s', 'focuspos', 'rapos', 'ra_move', 'ra_s', 'pa',
        # , 'numbrms', 'stdcrms',
        'fwhm', 'seeing', 'shift', 't', 'wcscompl', 'vi_plus', 'vi_minus',
    }

    # .upper changes string to uppercase
    # make new list of all headers/data types by combining current headers with extra headers
    full_imagelist_data_type = [
        (key.upper(), typ)
        for (key, typ) in imagelist_from_header_data_dtype + imagelist_extra_dtype
    ]

    # create an nimages size array where each element contains full_imagelist_data_type
    imagelist_data = np.recarray(nimages, dtype=full_imagelist_data_type)

    logger.debug('Allocating memory for image HDUs')
    flux_hdu_indexes = list(range(1, 14))
    image = lambda name: FITSImage(name, nimages, napertures)
    image_names = ['JD', 'HJD', 'BJD', 'FLUX', 'FLUXERR', 'CCDX', 'CCDY', 'SKYBKG', 'PEAK']
    image_names.extend(['FLUX_{}'.format(i) for i in flux_hdu_indexes])
    image_names.extend(['ERROR_{}'.format(i) for i in flux_hdu_indexes])
    image_map = {name: image(name) for name in image_names}

    aperture_sizes = first.get_aperture_sizes(flux_hdu_indexes)

    # For every image file...
    logger.info('Iterating over files')
    for i, filename in enumerate(sorted_images):
        logger.info(filename)
        source = SourceFile.open_file(filename)
        image_filename = filename.replace('.phot', '')
        source_extract = SourceFile.open_file(image_filename.replace('.'+args.ext, '.cat')) #'.'+args.ext, '.cat'))

        logger.debug('Extracting header data')
        jd = source.header['jd-obs']
        bjd = source.header['bjd-obs']
        hjd = source.header['hjd-obs']
        exposure = source.header['exptime']
        # set the time of image as the mid-exposure time:
        imagelist_data['TMID'][i] = jd + ((exposure/86400)/2)

        for (key, typ) in imagelist_from_header_data_dtype:
            if key.lower() in optional_keys:
                imagelist_data[key.upper()][i] = source.header.get(key, 0)
            else:
                try:
                    imagelist_data[key.upper()][i] = source.header[key]
                except Exception as err:
                    print(filename)
                    imagelist_data[key.upper()][i] = 0
                    print("Missing Keyword!: {}".format(str(err)),
                          file=sys.stderr)

        # WCS headers that have been renamed (now guaranteed to exist)
        imagelist_data['CTYPE1'][i] = source.header['TCTYP3']
        imagelist_data['CTYPE2'][i] = source.header['TCTYP5'] #6
        imagelist_data['CRPIX1'][i] = source.header['TCRPX3']
        imagelist_data['CRPIX2'][i] = source.header['TCRPX5'] #6
        imagelist_data['CRVAL1'][i] = source.header['TCRVL3']
        imagelist_data['CRVAL2'][i] = source.header['TCRVL5'] #6
        imagelist_data['CD1_1'][i] = source.header['TC3_3']
        imagelist_data['CD1_2'][i] = source.header['TC3_5'] #6
        imagelist_data['CD2_1'][i] = source.header['TC5_3']
        imagelist_data['CD2_2'][i] = source.header['TC5_5'] #6

        imagelist_data['NSOURCES'][i] = len(source_extract)

        logger.debug('Extracting image data')
        image_map['JD'].set_data(i, jd + ((exposure/86400)/2))  #+source.data['hjd_correction'])
        image_map['HJD'].set_data(i, hjd + ((exposure/86400)/2))
        image_map['BJD'].set_data(i, bjd + ((exposure / 86400) / 2))
        image_map['FLUX'].set_data(i, source.data['Aper_flux_3'])
        image_map['FLUXERR'].set_data(i, source.data['Aper_flux_3_err'])
        image_map['CCDX'].set_data(i, source.data['X_coordinate'])
        image_map['CCDY'].set_data(i, source.data['Y_coordinate'])
        image_map['SKYBKG'].set_data(i, source.data['Sky_level'])
        image_map['PEAK'].set_data(i,source.data['Peak_height'])

        for image_key in flux_hdu_indexes:
            hdu_key = 'FLUX_{}'.format(image_key)
            error_key = 'ERROR_{}'.format(image_key)
            source_flux_key = 'Aper_flux_{}'.format(image_key)
            source_error_key = '{}_err'.format(source_flux_key)
            image_map[hdu_key].set_data(i, source.data[source_flux_key])
            image_map[error_key].set_data(i, source.data[source_error_key])

        del source

    logger.info('Post processing')
    logger.info(f'Successfully processed {nimages} images')
    if excluded_images:
        logger.info(f'Excluded {len(excluded_images)} images due to WCS failures')
    imagelist_data['LOCOUNT'] = np.zeros(nimages)
    imagelist_data['HICOUNT'] = np.zeros(nimages)
    catalogue_data['FLUX_MEAN_'+str(args.date)] = np.nanmean(image_map['FLUX'].data, axis=1)
    medflux = np.nanmedian(image_map['FLUX'].data, axis=1)
    catalogue_data['FLUX_MEDIAN_'+str(args.date)] = medflux #np.median(image_map['FLUX'].data, axis=1)
    # print(medflux)
    # output_date = args.outputdir + '/' + args.date + '/' + args.oldtarget + '/' + args.date + '_' + args.outname
    # logger.info('Rendering to %s', output_date)

    # ------------------------------------------------------------------
    # Read identified targets from stack catalogue TARGET_ROLE column.
    # Target identification was performed by catalogue_fov.py (Stage T8)
    # and written into the Gaia_Crossmatch extension. We do not repeat
    # that logic here.
    # ------------------------------------------------------------------
    logger.info("Reading identified targets from stack catalogue")
    targets = []
    try:
        with fitsio.FITS(args.cat, 'rw') as catfile:
            gaia_ext = catfile['Gaia_Crossmatch']
            cat_dr2_ids = gaia_ext['gaia_dr2_id'].read()
            cat_roles = gaia_ext['target_role'].read()
            cat_teffs = gaia_ext['teff'].read()

        seen = set()
        for dr2_id, role, teff_val in zip(cat_dr2_ids, cat_roles, cat_teffs):
            dr2_id_clean = dr2_id.strip()
            role_clean = role.strip()
            if role_clean in ('primary', 'secondary') and dr2_id_clean not in seen:
                seen.add(dr2_id_clean)
                targets.append((dr2_id_clean, role_clean, teff_val))

    except Exception as e:
        logger.error("Failed to read TARGET_ROLE from catalogue '%s': %s",
                     args.cat, e)
        return

    if not targets:
        logger.error(
            "No identified targets found in catalogue '%s' — "
            "TARGET_ROLE column may be absent or empty. "
            "Was catalogue_fov.py (Stage T8) run successfully?",
            args.cat
        )
        return

    primary = [t for t in targets if t[1] == 'primary']
    if not primary:
        logger.warning(
            "=" * 60
        )
        logger.warning(
            "WARNING: No primary target in catalogue '%s'. "
            "Processing %d secondary target(s) only.",
            args.cat, len(targets)
        )
        logger.warning(
            "=" * 60
        )
    else:
        logger.info(
            "Found %d target(s): %d primary, %d secondary",
            len(targets), len(primary), len(targets) - len(primary)
        )

    # Build a DR2-ID → OBJ_ID lookup from the catalogue data.
    # gaia_dr2_id and obj_ids are parallel arrays of the same length
    # (napertures), so we can map directly.
    dr2_to_obj_id = {}
    for idx, gid in enumerate(gaia_dr2_id):
        gid_clean = str(gid).strip()
        if gid_clean not in ('', 'nan'):
            dr2_to_obj_id[gid_clean] = obj_ids[idx]

    # Load target list for Teff fallback — only done once, not per target.
    _tlist_data = None
    if hasattr(args, 'tlist') and args.tlist and os.path.isfile(args.tlist):
        try:
            _tlist_data = astropy_ascii.read(
                args.tlist, delimiter=' ', header_start=0, data_start=1
            )
        except Exception as _e:
            logger.warning("Could not read target list for Teff fallback: %s", _e)

    def _get_teff(dr2_id, gaia_teff_val):
        """
        Return a Teff header value. Use Gaia teff_gspphot first; if that
        is missing or NaN, fall back to the 40pc target list value.
        """
        # Try Gaia value first
        try:
            val = float(gaia_teff_val)
            if not (val != val):  # NaN check without importing math
                return int(round(val))
        except (TypeError, ValueError):
            pass

        # Fallback: target list
        if _tlist_data is not None:
            try:
                ids = [str(x).strip() for x in _tlist_data['Gaia_ID,']]
                if dr2_id.strip() in ids:
                    idx = ids.index(dr2_id.strip())
                    teff_tl = _tlist_data['T_eff,'][idx]
                    return int(teff_tl)
            except Exception:
                pass

        return "N"

    for (dr2_id, role, teff_val) in targets:
        output_date = (args.outputdir + '/' + args.date + '/' +
                       args.targname + '/' + dr2_id + "_" +
                       args.filter + "_" + args.date + '_' + args.outname)
        logger.info('Rendering to %s', output_date)

        # SP_ID: the OBJ_ID (e.g. SP000047) of the aperture corresponding
        # to this target in the photometry catalogue. Identifies which row
        # of the flux arrays belongs to the science target.
        sp_id = dr2_to_obj_id.get(dr2_id.strip(), "N")

        teff_header = _get_teff(dr2_id, teff_val)

        # create a new fits image with name given by output
        with create_output_file(output_date) as hdulist:
            print('add imagelist_data')
            hdulist.append(fits.BinTableHDU(imagelist_data, name='IMAGELIST'))
            print('add catalogue_data')
            hdulist.append(fits.BinTableHDU(catalogue_data, name='CATALOGUE'))
            print('add image_map')
            for image in list(image_map.values()):
                hdulist.append(image.hdu())

            hdulist[0].header['NDATES'] = 1
            hdulist[0].header['VERSION'] = args.version
            hdulist[0].header['GAIA_DR2_ID'] = dr2_id
            hdulist[0].header['TARGET_ROLE'] = role
            hdulist[0].header['SP_ID'] = sp_id
            hdulist[0].header['TEFF'] = teff_header
            hdulist[0].header['history'] = args.date

            if args.sha:
                key_length = 72
                logger.info('Rendering pipeline sha %s', args.sha)
                short_sha = args.sha[:key_length]
                logger.debug('short sha %s', short_sha)
                hdulist[0].header['PIPESHA'] = (short_sha, 'git sha of the pipeline')

            print('add flux hdu indices')
            for index in flux_hdu_indexes:
                flux_str = aperture_sizes[index]['flux']
                err_str = aperture_sizes[index]['error']
                hdulist['FLUX_{}'.format(index)].header['radius'] = flux_str
                hdulist['ERROR_{}'.format(index)].header['radius'] = err_str

            hdulist.close()

        output = args.outputdir + '/' + dr2_id + "_" + args.filter + "_" + args.outname

        if os.path.exists(output):
            exists = True
            print("Global output file already exists!")
        else:
            exists = False

        try:
            if exists:
                f = fits.open(output)
                hdr = f[0].header
                hist = hdr['HISTORY']
                f.close()
                print("Update the global output file")
                with fits.open(output, mode='update') as hdulist_update:
                    print(hdulist_update['CATALOGUE'].header)
                    print(hdulist_update['CATALOGUE'].header.get('FLUX_MEAN_'+str(args.date)))
                    header = hdulist_update[0].header
                    if 'NDATES' in header and args.date not in hist:
                        header['NDATES'] = header['NDATES'] + 1
                        header['history'] = args.date
                        header['VERSION'] = args.version
                        header['GAIA_DR2_ID'] = dr2_id
                        header['TARGET_ROLE'] = role
                        header['SP_ID'] = sp_id
                        header['TEFF'] = teff_header

                    # extract existing column names
                    print('updating imagelist_data')
                    name = 'IMAGELIST'
                    fnames = hdulist_update[name].data.dtype.names

                    imagelist_dict = [collections.OrderedDict(list(zip(fnames, x))) for x in imagelist_data]
                    hdu1 = fits.BinTableHDU(Table(imagelist_dict))
                    hdu,to_replace = rowappend_fits_tables(hdulist_update[name], hdu1, name,False)
                    to_replace = True
                    if to_replace:
                        # at the moment we extract and add the new data to the previous hdulist and then delete the previous hdulist
                        # and append the combined one - think this can be improved
                        fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
                        hdulist_update.append(hdu)

                        print('update catalogue_data')
                        name = 'CATALOGUE'
                        import_cat = hdulist_update[name]
                        orig_cols = import_cat.columns
                        # new_cols = fits.ColDefs([fits.Column(name='FLUX_MEAN_' + str(args.date), format='D',
                        #                                      array=catalogue_data['FLUX_MEAN_' + str(args.date)]),
                        #                          fits.Column(name='FLUX_MEDIAN_' + str(args.date), format='D',
                        #                                      array=catalogue_data['FLUX_MEDIAN_' + str(args.date)])])
                        # new_cat = fits.BinTableHDU.from_columns(orig_cols + new_cols, name=name)
                        # fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
                        # hdulist_update.append(new_cat)

                        if 'FLUX_MEAN_' + str(args.date) in orig_cols.names:
                            print("This date already exists in the catalogue!")
                            print("Replace with new columns")
                            hdulist_update[name].data['FLUX_MEAN_' + str(args.date)]= catalogue_data['FLUX_MEAN_' + str(args.date)]
                            hdulist_update[name].data['FLUX_MEDIAN_' + str(args.date)]= catalogue_data['FLUX_MEDIAN_' + str(args.date)]

                        else:
                            print("This date does not exist in the catalogue")
                            new_cols = fits.ColDefs([fits.Column(name='FLUX_MEAN_' + str(args.date), format='D',
                                                                 array=catalogue_data['FLUX_MEAN_' + str(args.date)]),
                                                     fits.Column(name='FLUX_MEDIAN_' + str(args.date), format='D',
                                                                 array=catalogue_data['FLUX_MEDIAN_' + str(args.date)])])

                            new_cat = fits.BinTableHDU.from_columns(orig_cols + new_cols, name=name)
                            fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
                            hdulist_update.append(new_cat)

                        print('update image_map')

                        cntr = 0

                        to_del_night = find_dups(np.transpose(hdulist_update['JD'].data),
                                                 np.transpose(image_map['JD'].data))
                        to_del_hdu = find_dups(np.transpose(image_map['JD'].data),
                                               np.transpose(hdulist_update['JD'].data))

                        # image_map is a dictionary of hdulist name: data of [0:number of sources] where each element is an ndarray of values
                        for image in list(image_map.values()):
                            hdu = None
                            # print(image.name)
                            cntr = cntr+1
                            name = image.name
                            # print(name)
                            # hdu = colappend_fits_tables(hdulist_update[name], image.hdu(), name)
                            # fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
                            # hdulist_update.append(hdu)

                            action = decide_action(to_del_night, image.data.shape[1], False)

                            if action == "none":
                                if cntr == 1:
                                    print(str(name) + " already exists!")

                            elif action == "replace":
                                if cntr == 1:
                                    print(str(name) + " already exists! Replace with new version.")
                                newdata = []
                                for n in range(len(hdulist_update[name].data)):
                                    newrow = np.delete(hdulist_update[name].data[n], to_del_hdu)
                                    newdata.append(newrow)
                                hdulist_update[name].data = np.array(newdata)

                            elif action == "add":
                                if cntr == 1:
                                    print(str(name) +  " doesn't Exist!")

                            elif action == "append":
                                if cntr == 1:
                                    print(str(name) + " partially exists! Add missing data!")
                                newdata = []
                                for n in range(len(image.data)):
                                    newrow = np.delete(image.data[n], to_del_night)
                                    newdata.append(newrow)
                                image.data = np.array(newdata)

                            hdu = colappend_fits_tables(hdulist_update[name], image, name)
                            fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
                            hdulist_update.append(hdu)

                        hdulist_update.flush()
                    # except ValueError:
                    #     print('Already run this date!')

                # else:
                #     print("This date is already in the global output fits file")

            else:
                copyfile(output_date, output)
        except Exception as e:
            print("Global LC FAILED, " + str(e))


def find_dups(prev_arr,new_arr):
    to_del = []
    for i in range(len(new_arr)):
        if new_arr[i] in prev_arr:
            to_del.append(i)

    return to_del

def decide_action(to_del,nrows,replace):
    if len(to_del) == nrows:
        if replace:
            return "replace"
        else:
            return "none"
    elif len(to_del)==0:
        return "add"
    else:
        if replace:
            return "replace"
        else:
            return "append"


def rowappend_fits_tables(hdu1, hdu2, name,replace):
    print('row append')
    nrows1 = hdu1.data.shape[0]
    nrows2 = hdu2.data.shape[0]
    # nrows = nrows1 + nrows2
    # hdu = fits.BinTableHDU.from_columns(hdu1.columns, nrows=nrows, name=name)
    # for colname in hdu1.columns.names:
    #     hdu.data[colname][nrows1:] = hdu2.data[colname]


    to_del2 = find_dups(hdu1.data['JD-OBS'],hdu2.data['JD-OBS'])
    to_del1 = find_dups(hdu2.data['JD-OBS'], hdu1.data['JD-OBS'])

    action = decide_action(to_del2, nrows2, replace)

    if action == "none":
        return hdu1, False

    elif action == "replace":
        print("Replace existing data")
        hdu1.data = np.delete(hdu1.data, to_del1)
        nrows1 = hdu1.data.shape[0]

    elif action == "add":
        print("This date is not in the global outputfits, add it")

    elif action == "append":
        print("This data is only partially in the global outputfits")
        mask = np.ones(hdu2.data.shape[0], dtype=bool)
        mask[to_del2] = False
        result = hdu2.data[mask]
        hdu2.data = result
        # hdu2.data = np.delete(hdu2.data,to_del2)
        for colname in hdu2.columns.names:
            if np.shape(hdu2.data[colname])[0]!=hdu2.data.shape[0]:
                print("Error with shape of column " + colname)

            # print(colname,np.shape(hdu2.data[colname]),type(hdu2.data[colname]),type(hdu2.data[colname][0]))
        nrows2 = hdu2.data.shape[0]

    nrows = nrows1 + nrows2


    hdu = fits.BinTableHDU.from_columns(hdu1.columns, nrows=nrows, name=name)
    for colname in hdu1.columns.names:
        hdu.data[colname][nrows1:] = hdu2.data[colname]


    return hdu, True

def colappend_fits_tables(hdu1,hdu2,name):
    # print("col append")
    # print(hdu1.data.shape[0])
    # print(hdu2.data.shape[0])
    hdu1.data = np.append(hdu1.data,hdu2.data,axis=1)
    return hdu1

if __name__ == '__main__':
    description = 'Condense per-image photometry into output FITS cube (Stage T10)'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filename', nargs='+')
    parser.add_argument('-o', '--outputdir', required=True)
    parser.add_argument('-u', '--update', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--sha', help='Store pipeline sha', required=False)
    parser.add_argument('-d', '--date')
    parser.add_argument('-tn', '--targname', required=True,
                        help='Target name (used for output directory path)')
    parser.add_argument('-on', '--outname', required=True)
    parser.add_argument('-f', '--filter', required=True)
    parser.add_argument('-e', '--ext', required=True)
    parser.add_argument('-ve', '--version', required=True)
    parser.add_argument('-c', '--cat', required=True,
                        help='Path to stack catalogue FITS file containing '
                             'Gaia_Crossmatch extension with TARGET_ROLE column')
    parser.add_argument('-tl', '--tlist', required=False, default=None,
                        help='Path to 40pc target list (for Teff fallback)')
    main(parser.parse_args())