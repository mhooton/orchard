#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Zero Level Pipeline catalog generation

Usage:
  ZLP_create_cat [options] --confmap=CONFMAP --filelist=FILELIST

Options:
  -h --help                                 Show help text
  -v, --verbose                             Print more text
  -o <OUTNAME>, --outname <OUTNAME>         Specify the name of the output catalog [default: catfile.fits]
  -s <STACKLIST>, --stacklist <STACKLIST>   The name of the file that stores the names of the images used in the stack [default: stackfilelist]
  --c_thresh <C_THRESH>                     The detection threshold to use when defining the input [default: 2]
  --s_thresh <S_THRESH>                     The detection threshold to use when WCS solving images - typically higher than when doing actual photometry [default: 20]
  -n <NPROC>, --nproc <NPROC>               Enable multithreading if you're analysing a lot of files at once
  -N <NFILES>, --nfiles <NFILES>            Maximum number of files to use in the stack
  --no-wcs                                  Do not solve each image for WCS.  However images must have a solution somehow

This is the catalog generation tool, requires a filelist input. need to work on being selective on the files used in input.

"""

#from docopt import docopt
import os
import argparse
import sys
import photometry.casutools as casutools
import photometry.gaia_dr2_test as gaia_dr2
from photometry.wcs_status import wcs_succeeded
from tempfile import NamedTemporaryFile
from calibration.pipeutils import open_fits_file
from astropy.io import fits
import numpy as np
import string
from photometry.wcs_fitting import m_solve_images, compute_frame_limits, reference_catalogue_objects
from photometry.vector_plot import wcsf_QCheck
from reporting.QC import main as QC
import utils.gaia_id_from_schedule
import utils.already_run_targs

print(os.getcwd())


def main(filelist,outname,catname,outdir,reportdir, filter, date,
         no_wcs= False,
         c_thresh= 2,
         s_thresh= 20,
         verbose= False,
         nfiles= None,
         ipix=6,
         catsrc='vizgaia2',
         rcore=4,
         ncores=1,
         ext='fits'): #argv):


    stack_filelist = os.path.dirname(filelist) + '/'+filter+'_stacked.dat'

    if verbose == True:
        print('Creating source catalogue from first {} images...'.format(nfiles))

    # Pick the first N files if argument given
    # nfiles = int(argv['--nfiles']) if argv['--nfiles'] else None
    if nfiles != 'None':
        nfiles = int(nfiles)
    else:
        nfiles = None

    count = 0
    # print filelist
    #NamedTemporaryFile() returns a file-like object that can be used temporarily but with a visible name in the file system
    with NamedTemporaryFile(mode='w+') as tmp:
        name = tmp.name
        #write the file names of the processed images to be stacked to this temporary file
        with open(filelist,'r') as infile:
            # print type(infile)
            lines = infile.readlines()
            for i, line in enumerate(lines):
                if nfiles and count >= nfiles: #(i-len(lines)/2) >= nfiles:
                    print(i)
                    print(len(lines))
                    break

                # if we only have a small number of files (<nfiles) then use them all:
                if len(lines)<nfiles:
                    print(str(i) + ': ' + line)
                    tmp.write(line)
                    count = count + 1

                elif i>= len(lines)/2:
                    # stack from images in the middle of the night as they have the lowest bkg
                    # if there are more than 1000 images then these will be in the middle of the filelist
                    # even though they're taken at the end of the night. Therefore ignore any images
                    # > 1000
                    # if int(line.split("-")[-2][1:]) < 1000:
                    # if len(lines) < 1000:
                        # print i-len(lines)
                    if 'dupe' not in line:
                        print(str(i) + ': ' + line)
                        tmp.write(line)
                        count=count+1



        #Sets the r/w pointer to start of file
        tmp.seek(0)

        if no_wcs=='False':#argv['--no-wcs']:
            #call a function to wcs solve the images in order to stack
            print('\n**wcs solving images**\n')
            m_solve_images(
                name, name,
                thresh= s_thresh,#argv['--s_thresh']
                nproc= None, #int(argv['--nproc']) if argv['--nproc'] else None,
                verbose= verbose,
                rcore=rcore,
                catsrc=catsrc,
                ipix=ipix,
                ext=ext) #argv['--verbose'])
        else:
            print('no wcs solving')

        count = 0
        with open(stack_filelist, 'w') as stacklist: #argv['--stacklist']
            for line in tmp:
                image = line.strip()
                # print image

                # with open_fits_file(image) as hdulist:
                #     field = hdulist[0].header['OBJECT']
                    # print 'TARGET:' + str(field)

                status_check = ['ok', 'ok']
                # REMOVE THIS BLOCK, IT BARELY MATTERS WHETHER OR NOT CASUTOOLS RUNS
                # if wcs_succeeded(image) == False:
                #     status_check.append('not_ok')
                #     print image + ' could not be WCS fit and will not be included in stacked image'

                if all([status == 'ok' for status in status_check]):
                    if count==1:
                        stacklist.write('\n' + image)
                    else:
                        stacklist.write(image)
                        count =1
                    # stacklist.write(image + '\n')

    outstack_name = outdir + '/' + outname
    outstackconf_name = outdir +'/'+ outname.replace('.'+ext, '_conf.'+ext)
    outcatname = outdir + '/' + catname
    print('\n**Stacked image stored as:**\n')
    print(outstack_name)
    print('\n**Stacked image catalogue stored as:**\n')
    print(outcatname)

    # stack_filelist_test = stack_filelist + '_test'
    casutools.imstack(stack_filelist,
                      confidence_map=None,
                      verbose= verbose,
                      outstack=outstack_name,
                      outconf=outstackconf_name)

    #TEST
    # median_filter_stack(outstack_name)

    trim_stack(outstack_name)
    trim_stack(outstackconf_name)
    casutools.imcore(outstack_name, outcatname, #argv['--outname'],
                     threshold= c_thresh, #argv['--c_thresh'],
                     confidence_map=outstackconf_name,
                     verbose= verbose,
                     ellfile = True) #argv['--verbose'])

    # print outcatname, outcatname+"_Gaia.log", ncores, ext, date
    perc,num_sources = gaia_dr2.crossmatch(outcatname,logfile=outcatname+"_Gaia.log", n=ncores, ext=ext, date=date)
    dirname=os.path.dirname(filelist)
    dirsplit = dirname.split("/")
    run = dirsplit[-1]
    field = dirsplit[-2]
    date = dirsplit[-3]
    tel = dirsplit[-5]
    # print dirsplit
    # print perc
    # print num_sources
    # print dirsplit[-1]
    # print dirsplit[-2]
    # print dirsplit[-3]
    # print dirsplit[-5]
    dirname= "/".join(dirsplit[:-4])

    try:
        QC([0, reportdir + "/QC", date, field, tel, run, 'QC11', perc])
        QC([0, reportdir + "/QC", date, field, tel, run, 'QC10', num_sources])
    except:
        print("ERROR: Couldn't update report")
    # crossmatchGaia(outcatname, outstack_name)

    if verbose == True: #argv['--verbose'] == True:
        print('Catalogue complete')

def trim_stack(fname):
    from scipy.ndimage import median_filter

    hdulist_old = fits.open(fname)
    hdr_old = hdulist_old[0].header
    data_old = hdulist_old[0].data
    print(np.shape(data_old))
    data_new = data_old[:,20:-20]

    # subtract 20 from the reference pixel coordinates to account for trimming
    if hdr_old.get('CRPIX1') != None:
        hdr_old.set('CRPIX1', hdr_old['CRPIX1'] - 20)
        # hdulist[0].header['CRPIX2'] = hdulist[0].header['CRPIX2']
        hdr_old.add_history('Subtracted 20 pixels from CRPIX1 to account for trimming of stack')

    fits.writeto(fname, data_new, hdr_old, overwrite=True)

    hdulist = fits.open(fname)
    print(np.shape(hdulist[0].data))

def median_filter_stack(fname):
    from scipy.ndimage import median_filter

    hdulist_old = fits.open(fname)
    hdr_old = hdulist_old[0].header
    data_old = hdulist_old[0].data
    data_new = median_filter(data_old,size=3)

    fits.writeto(fname, data_new, hdr_old, overwrite=True)


def crossmatch2MASS(catfile_name, stackimage):
    catpath = os.path.join(os.getcwd(), 'catcache')

    #Come back to this testing section later - perhaps break out into separate script?
    #find frame limits from index file in catcache directory
    catalogue = compute_frame_limits(catpath)
    #retrieve JMAG, RA and DEC of reference catalogue objects
    cat = reference_catalogue_objects(catalogue, catpath)

    with fits.open(catfile_name) as mycatt:
        # extract stacked image catalogue and extra X, Y coordinates and flux of each source
        mycatt_data = mycatt[1].data
        #print mycatt[1].header
        mycat = {'Aper_flux_3': mycatt_data['Aper_flux_3']}

        my_X = mycatt_data['x_coordinate']
        my_Y = mycatt_data['y_coordinate']
        my_ID = mycatt_data['Sequence_number']

    # Do QC checks. should really break this out.
    # plot separation of each source from the ref catalogue as a quality check
    ids, jmags, kmags, hmags = wcsf_QCheck(mycat, stackimage,
               os.path.basename(stackimage).replace('.fits', '') + '.png', cat,
               catalogue.ra_lims, catalogue.dec_lims, my_X, my_Y, my_ID,
               plot=True)

    sorttuple = [x for x in sorted(zip(ids, jmags,kmags,hmags))]
    jmag = [x[1] for x in sorttuple]
    kmag = [x[2] for x in sorttuple]
    hmag = [x[3] for x in sorttuple]
    seqnum = [x[0] for x in sorttuple]

    print("Adding JMAG/KMAG/HMAG to catalogues from crossmatch with 2MASS")

    for i in np.arange(1,len(my_ID)+1):
        if i not in seqnum:
            print("No crossmatch for SEQNUM " + str(i))
            jmag = np.insert(jmag,i-1,np.nan)
            kmag = np.insert(kmag, i - 1, np.nan)
            hmag = np.insert(hmag, i - 1, np.nan)

    file = open(catfile_name.replace('.fts','_2MASS.dat'),'w')
    for j in range(len(jmag)):
        file.write('%s\t%s\t%s\n' % (jmag[j],kmag[j],hmag[j]))
    file.close()

def crossmatchGaia(catfile_name, stackimage):

    catpath = os.path.join(os.getcwd(), 'catcache')

    # Come back to this testing section later - perhaps break out into separate script?
    # find frame limits from index file in catcache directory
    catalogue = compute_frame_limits(catpath)
    # retrieve JMAG, RA and DEC of reference catalogue objects
    cat = reference_catalogue_objects(catalogue, catpath)

    with fits.open(catfile_name) as mycatt:
        # extract stacked image catalogue and extra X, Y coordinates and flux of each source
        mycatt_data = mycatt[1].data
        # print mycatt[1].header
        mycat = {'Aper_flux_3': mycatt_data['Aper_flux_3']}

        my_X = mycatt_data['x_coordinate']
        my_Y = mycatt_data['y_coordinate']
        my_ID = mycatt_data['Sequence_number']

    # Do QC checks. should really break this out.
    # plot separation of each source from the ref catalogue as a quality check
    ids, gmags, gaia = wcsf_QCheck(mycat, stackimage,
                                           os.path.basename(stackimage).replace('.fits', '') + '.png', cat,
                                           catalogue.ra_lims, catalogue.dec_lims, my_X, my_Y, my_ID,
                                           plot=True)

    sorttuple = [x for x in sorted(zip(ids, gmags, gaia))]
    gmag = [x[1] for x in sorttuple]
    gaiaid = [str(x[2]) for x in sorttuple]
    seqnum = [x[0] for x in sorttuple]

    # print gaiaid[5]
    # print "%.19f" % gaiaid[5]
    # print str(gaiaid[5])

    print("Adding GMAG to catalogues from crossmatch with Gaia")

    for i in np.arange(1, len(my_ID) + 1):
        if i not in seqnum:
            print("No crossmatch for SEQNUM " + str(i))
            gmag = np.insert(gmag, i - 1, np.nan)
            gaiaid = np.insert(gaiaid, i - 1, np.nan)

    file = open(catfile_name.replace('.fts', '_Gaia.dat'), 'w')
    for j in range(len(gmag)):
        file.write('%s\t%s\n' % (gmag[j], gaiaid[j]))
    file.close()


def create_cat():
    #This script is only run once for each field so that for later observations of the same fields the images will not
    #be wcs solved twice

    # print os.getcwd()
    parser = argparse.ArgumentParser()
    parser.add_argument('filelist')
    parser.add_argument('outstackname')
    parser.add_argument('outcatname')
    parser.add_argument('outdir')
    parser.add_argument('reportdir')
    parser.add_argument('filter')
    parser.add_argument('date')
    parser.add_argument('no_wcs')
    parser.add_argument('c_thresh')
    parser.add_argument('s_thresh')
    parser.add_argument('verbose')
    parser.add_argument('nfiles')
    parser.add_argument('ipix')
    parser.add_argument('catsrc')
    parser.add_argument('rcore')
    parser.add_argument('ncores')
    parser.add_argument('ext')
    args = parser.parse_args()

    filelist = args.filelist  #'/data/cam217/1_image_trappist-1.list/processed.dat'
    outname = args.outstackname
    catname = args.outcatname #'catstack.fts'
    outdir = args.outdir  #'/data/cam217/create_cat_output/'
    reportdir = args.reportdir
    # no_wcs determines whether or not we solve each image for wcs (True - we don't solve)
    no_wcs = args.no_wcs #False
    c_thresh = args.c_thresh  #2
    s_thresh = args.s_thresh  #20
    # nproc = 1
    verbose = args.verbose
    nfiles = args.nfiles
    ipix = args.ipix
    catsrc = args.catsrc
    rcore= args.rcore
    ncores = args.ncores
    ext=args.ext
    filter=args.filter
    date = args.date

    main(filelist,outname,catname,outdir,reportdir,filter,date,no_wcs,c_thresh,s_thresh,verbose,nfiles, ipix,catsrc, rcore,ncores,ext) #(docopt(__doc__))


if __name__ == '__main__':
    create_cat()