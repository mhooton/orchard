# -*- coding: utf-8 -*-
from astropy.io import fits as pf
import os
import sys
import linecache
import threading
from os.path import isfile, join
import multiprocessing
from multiprocessing import Pool
from functools import partial
import numpy as np
from photometry.wcs_status import wcs_succeeded
from photometry.hjd_correction import append_hjd_correction_column
from photometry.photom_quality_checks import frame_shift, m_frame_shift, cloud_check
from photometry.super_sample import call_find_fwhm
from photometry import casutools
import psutil


def m_wcs_photom(filelist, outlist, appsize, cat_file,
                 nproc=1,
                 verbose=False):

    # cat_file is the catalogue from the stacked image
    # by default outlist = filelist + '_phot'
    # this command will copy the filelist contents to outlist using a shell command
    os.system('cp ' + filelist + ' ' + outlist)

    infiles = []
    with open(filelist) as infile:
        #for each science image file
        for line in infile:
            image = line.strip()
            #check it is a valid file to use
            status_checks = ['ok', 'ok']

            # don't need to perform a WCS succeed check here as we check elsewhere

            if all(status == 'ok' for status in status_checks):
                infiles.append(image)

    # with open(outlist,'w') as outfile:
    #     for f in infiles:
    #         outfile.write(f + ".phot")


    pool = Pool(nproc)

    print("POOL M_WCS_PHOTOM")

    fn = partial(handle_errors_in_wcs_photom,
                 cat_file=cat_file,
                 appsize=appsize,
                 verbose=verbose)
    pool.map(fn, infiles)

    infiles = remove_wcs_failed(infiles)

    first_image = infiles[0] + '.phot'
    pf.setval(first_image, 'SKY_MOVE', 1, value=0)

    # input the shift of 1st image from 1st image - this is just used to initialise these headers
    # for use by the following images
    RA_shift, DEC_shift, tot_shift, RA, DEC = frame_shift(infiles[0],
                                                          infiles[0])
    pf.setval(first_image, 'RA_MOVE', 1,
              value=RA_shift,
              comment='RA shift from previous image [arcseconds]')
    pf.setval(first_image, 'DEC_MOVE', 1,
              value=DEC_shift,
              comment='Dec shift from previous image [arcseconds]')
    pf.setval(first_image, 'SKY_MOVE', 1,
              value=tot_shift,
              comment='Total movement on sky [arcseconds]')

    pf.setval(first_image, 'WCSF_RA', 1, value=RA, comment='RA center pix')
    pf.setval(first_image, 'WCSF_DEC', 1, value=DEC, comment='Dec center pix')

    # calculate the RA and DEC shifts of each image from the previous:
    indexes = np.arange(1, len(infiles))
    fn = partial(m_frame_shift, infiles)
    pool.map(fn, indexes)

def remove_wcs_failed(infiles):
    new_infiles = []

    for i in infiles:
        if wcs_succeeded(i) == True:
            new_infiles.append(i)

    # print('difference in length before and after: ' + str(len(infiles) - len(new_infiles)))

    return new_infiles

def handle_errors_in_wcs_photom(image, *args, **kwargs):
    try:
        return wcs_photom(image, *args, **kwargs)
    except Exception as err:
        print("Exception handled in wcs_photom: {}".format(str(err)),
                file=sys.stderr)

def wcs_photom(image,cat_file='nocat',conf_file='noconf',appsize=4,verbose=False):

    #check whether the wcs succeeded for the image based on if it has 'wcscompl' (=T) in its header
    # print 'wcs_photom on image: ' + str(image)
    if not wcs_succeeded(image):
        return 'failed'

    outname = image + '.phot'
    print(outname)

    # update cat_file to include PM:
    casutools.imcore_list(image, cat_file, outname, confidence_map=conf_file,rcore=appsize, noell=False,
                        verbose=True)

     #      do some quality checks
    factor = 5
    size = 11
    stars = 100

    # extract fwhm in x and y directions and theta - angle of rotation of gaussian
    fwhm_a, fwhm_b, t = call_find_fwhm(image,cat_file,factor,size,stars,tag=image)

    cloud_status = cloud_check(image)

    # get 1st extension's keyword 'SEEING':
    pixel_fwhm = pf.getval(outname,'SEEING',1)
    plate_scale = 0.35
    # convert pixel_fwhm from pixels to arcseconds
    seeing = round(plate_scale*pixel_fwhm,2)

    pf.setval(outname,'CLOUD_S',1,value=round(cloud_status,2),comment='A measure of bulk structure in the image (S/N)')
    pf.setval(outname,'FWHM',1,value=pixel_fwhm,comment='[pixels] Average FWHM')
    pf.setval(outname,'SEEING',1,value=seeing,comment='[arcseconds] Average FWHM')
    # positions = ['Top left.','Top middle.','Top right.','Middle left.','Center.','Middle right.','Bottom left.','Bottom middle.','Bottom right.']
    positions = ['Top left.', 'Top middle.', 'Top right.', 'Middle left.', 'Center.', 'Middle right.', 'Bottom left.',
                 'Bottom middle.', 'Bottom right.']

    # loop through f_1 to f_9 - this shouldn't really have integers in it - should all be auto
    for val in range(1,10):
        # print(val)
        pf.setval(outname,'PSF_a_'+str(val),1,value=fwhm_a['f_'+str(val)][0],comment='[pixels] psf long axis FWHM. '+positions[val-1])
        pf.setval(outname,'PSF_b_'+str(val),1,value=fwhm_b['f_'+str(val)][0],comment='[pixels] psf short axis FWHM. '+positions[val-1])
        pf.setval(outname,'PSF_t_'+str(val),1,value=t['f_'+str(val)][0],comment='[degrees] psf rotation angle. '+positions[val-1])

    # Compute the HJD values
    # append_hjd_correction_column(outname)

    return 'ok'
