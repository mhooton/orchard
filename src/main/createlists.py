#!/usr/bin/env python
"""
Image Classification and List Generation

Classifies raw FITS images from a night's observations into science, bias, dark,
and flat field categories based on FITS header keywords (IMAGETYP, EXPTIME).
Creates organized file lists for each target and calibration frame type to drive
subsequent pipeline stages.

The script handles:
- Multiple image file extensions (.fits, .fts, .fz, .bz2)
- Target name sanitization and validation
- Filter-specific flat field organization
- Comprehensive logging of classification decisions

Output structure:
- {runname}_image_{target}.list - Science images per target
- {runname}_bias.list - Bias frames
- {runname}_dark.list - Dark frames
- {runname}_flat_{filter}.list - Flat fields per filter

This is Stage T1 of the pipeline and must complete successfully before
calibration stages can proceed.
"""

import glob
import os
import time
import sys
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
from functools import partial
from datetime import datetime,timedelta
from reporting.QC import main as QC
# open_fits_file presents a Rice-compressed frame (.fz) as the plain file, its image and header in HDU 0
from calibration.pipeutils import open_fits_file
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)

DIRECTORY = os.getcwd()
DEBUG = 0

# (NAXIS1, NAXIS2) of an Andor frame as ESO delivers it. Calibration (trim, overscan, masters) assumes the
# stored MaxIm orientation, 2048 x 2088, that download/astra_transform.py produces. Frames left in the native
# shape - never rotated, or rotated twice - were calibrated with masters in the other orientation without any
# warning (Callisto 2024, Ganymede 2025-11), so a night holding one is refused rather than reduced.
NATIVE_ANDOR_SHAPE = (2088, 2048)

def get_liste(directory, root, ext):
    ### search for all image files in the directories """
    # print directory
    liste = glob.glob("%s/*.%s" % (directory, ext))
    # search for images inside folders in the image folder (such as Calibration or AutoFlat)
    add_liste = glob.glob("%s/*/*.%s" % (directory, ext))

    # if len(liste) == 0:
    #     new_ext = 'fts'
    #     liste = glob.glob("%s/*.%s" % (directory, new_ext))
    #     add_liste = glob.glob("%s/*/*.%s" % (directory, new_ext))

    if len(add_liste) != 0:
        liste = liste + add_liste

    # frames stored Rice-compressed (X.fits.fz, utils/framepack.py), as an earlier night's calibration frames may be:
    # T2-T6 read them through pipeutils.open_fits_file.  The night being reduced has been decompressed for the run
    # (utils/rawpack.py), so a frame stored both ways is its working copy and counts once, as the plain file.
    if ext != "fz":
        have = set(liste)
        liste = liste + [f for f in glob.glob("%s/*.%s.fz" % (directory, ext)) + glob.glob("%s/*/*.%s.fz" % (directory, ext))
                         if f[:-3] not in have]

    if DEBUG:
        print("| get_liste")
        print("| ",)
        print("%s/*.%s" % (directory, ext))
    liste.sort()
    return (liste)


def classify(files, logroot, runnumber):
    pool = ThreadPool()
    fn = partial(classify_file, logroot=logroot, runnumber=runnumber)
    results = pool.map(fn, files)
    keys = ['dithered', 'science', 'dark', 'bias', 'flat', 'native_andor']
    pool.close()

    return {key: sum([entry[key] for entry in results], []) for key in keys}


def classify_file(filename, logroot, runnumber):
    out = defaultdict(list)
    #print filename
    try:
        with open_fits_file(filename) as hdulist:
            #header name is different for ngts (IMGTYPE) and trappist (IMAGETYP)
            try:
                imtype = hdulist[0].header['IMAGETYP']
                naxis1 = hdulist[0].header['NAXIS1']
                naxis2 = hdulist[0].header.get('NAXIS2')
                filter = hdulist[0].header['FILTER']
                # print("normal")
            except:
                pass
                # print("abnormal")
            if (naxis1, naxis2) == NATIVE_ANDOR_SHAPE:
                out['native_andor'].append(filename)
                return out
            # print(naxis1)
            #action = hdulist[0].header['ACTION']
            try:
                #not relevant for SPECULOOS
                dither = hdulist[0].header['DITHER']
            except:
                dither = 'DISABLED'
        string = "%20s %10s %30s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), imtype, filename)
        result = write_log(logroot, runnumber, string, 2)
        # print(str.upper(imtype))
        #changed because there is inconsistency across the image types: dark frame, light image, bias frame, flat
        if 'LIGHT' in str.upper(imtype): #imtype == 'Light Frame':
            if dither == 'ENABLED':
                out['dithered'].append(filename)
            else:
                out['science'].append(filename)
        elif 'DARK' in str.upper(imtype): #imtype == 'Dark Frame':
            out['dark'].append(filename)
        elif 'BIAS' in str.upper(imtype): #imtype == 'Bias Frame':
            out['bias'].append(filename)
        elif 'FLAT' in str.upper(imtype): #imtype == 'FLAT':
            out['flat'].append(filename)
    except:
        print("Unable to open corrupted FITS file: " + str(filename))

    return out


def refuse_native_andor(files):
    """Stop T1 (and with it the run, which is under errexit) if any frame is in the native Andor shape. This
    also covers frames of earlier nights borrowed for calibration."""
    if not files:
        return
    files = sorted(files)
    sys.exit("ERROR: {} frame(s) in the native Andor orientation ({} x {}), e.g. {}. Turn them into the stored "
             "orientation with download/repair_orientation.py before reducing this night.".format(
                 len(files), NATIVE_ANDOR_SHAPE[0], NATIVE_ANDOR_SHAPE[1], files[0]))


def sort_liste(liste, logroot, runnumber):
    ### Sort images for bias, darks, flats, and scientific frames"""
    # print liste
    classification = classify(liste, logroot, runnumber)
    refuse_native_andor(classification['native_andor'])
    biaslist = classification['bias']
    darklist = classification['dark']
    flatlist = classification['flat']
    science = classification['science']
    dithered = classification['dithered']
    # print(classification)

    #print the number of bias, dark, flat and science images identified
    string = "\n"
    string = string + "%20s %8d bias images identified\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), len(biaslist))
    string = string + "%20s %8d dark images identified\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), len(darklist))
    string = string + "%20s %8d flat images identified\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), len(flatlist))
    string = string + "%20s %8d science images identified\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), len(science))
    string = string + "%20s \n" % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + "# ------------------------------------------------------\n"
    string = string + "# \n"
    result = write_log(logroot, runnumber, string, 1)
    result = write_log(logroot, runnumber, string, 2)
    return (biaslist, darklist, flatlist, science, dithered)

def sort_flatlist(liste):
    #sort flat images based on the filter used
    flatlists = []
    filters = []

    for item in liste:
        with open_fits_file(item) as hdulist:
            filter = hdulist[0].header['FILTER']
            # print filter

        if filters.count(filter) == 0:
            filters.append(filter)
            flatlists.append([])

        idx = filters.index(filter)
        flatlists[idx].append(item)

    return (filters,flatlists)


def sort_scilist(liste):
    #sort science images based on the object we're targetting
    fields = []
    scilists = []
    filters = []
    scifilt={}

    for item in liste:
        #call open_fits_file procedure from pipeutils
        with open_fits_file(item) as hdulist:
            try:
                field = hdulist[0].header['OBJECT']
            except:
                field = hdulist[0].header['FIELD']

            if '_2' in field and len(field)==13:
                field = field.replace("_2", "")

            filter = hdulist[0].header['FILTER']
            # print field
            # ra = hdulist[0].header['RA']
            # dec = hdulist[0].header['DEC']

        # field = field.upper()
        # if "SP" in field and "WASP" not in field and "TESS" not in field:
        #     field = fix_target_name(field, ra, dec)

        if fields.count(field) == 0:
            #for every object targeted add this object to the array 'fields'
            fields.append(field)
            #and add an empty entry to scilists for each new object
            scilists.append([])

        if filters.count(filter) == 0:
            filters.append(filter)

        #replace this entry in scilists with the data from hdulist
        idx = fields.index(field)
        scilists[idx].append(item)
        scifilt[field] = filter

    return (fields, scilists,filters, scifilt)


def write_liste(liste, filename, outputdir, logroot, runnumber, date):
    """ write output files """
    # print imagedir
    outdir = "%s/%s" % (outputdir, date)
    output = "%s/reduction/%s" % (outdir, filename) #changed from DIRECTORY

    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    f = open(output, 'w')
    for item in liste:
        f.write("%s\n" % item)
    f.close()
    string = "%20s List written  (%4d entries) in file %40s\n" % (
    time.strftime("%Y-%m-%d %H:%M:%S"), len(liste), filename)#output)
    write_log(logroot, runnumber, string, 1)
    write_log(logroot, runnumber, string, 2)
    return (0)


def write_log(logroot, runnumber, string, lognumber=2):
    if DEBUG:
        print("| Write_log", "| ", logroot, runnumber, string, lognumber, "|")
    if lognumber == 1:
        logtype = 'short'
    elif lognumber == 2:
        logtype = 'long'

    # f = open("%s/logfiles/%s_%s_%03d.log" % (DIRECTORY, logroot, logtype, runnumber), 'a')
    # f.write(string)
    # f.close()
    if lognumber == 1:
        print(string)
    return (0)


def check_log(logroot):
    runnumber = 1
    while os.path.exists("%s/logfiles/%s_short_%03d.log" % (DIRECTORY, logroot, runnumber)):
        runnumber += 1
    return (runnumber)


def write_logstart(directory, imageroot, ext, run, runnumber):
    if DEBUG:
        print("| write_logstart","| ",directory, imageroot, ext, run, runnumber)

    string = "-------------------------------------------------\n"
    string = string + "%20s Creating lists of images for Zero Level pipeline\n " % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + "%20s using the script createlists.py\n " % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + "%20s \n " % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + "%20s All %s*.%s files in the directory \n " % (
    time.strftime("%Y-%m-%d %H:%M:%S"), imageroot, ext)
    string = string + "%20s %s \n " % (time.strftime("%Y-%m-%d %H:%M:%S"), directory)
    string = string + "%20s will be sorted.\n " % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + "%20s \n " % time.strftime("%Y-%m-%d %H:%M:%S")
    string = string + "%20s Working Directory is %s\n " % (time.strftime("%Y-%m-%d %H:%M:%S"), DIRECTORY)
    string = string + "%20s \n " % (time.strftime("%Y-%m-%d %H:%M:%S"))
    string = string + " -----------------------------------------------\n "
    string = string + "# \n "
    result = write_log(run, runnumber, string, 1)
    result = write_log(run, runnumber, string, 2)
    return (0)

def main(args):
    # SRW: strips quotation marks so this script can be submitted
    #  as a queue job

    #for now input manually:
    # directory = sys.argv[1].replace('"','') #DIRECTORY  # sys.argv[1].replace('"', '')
    imageroot = args[1].replace('"','')
    #imageroot - folder location of fits files:
    directory = args[2] #sys.argv[2] #'example_fits\\all_fits\\'
    outputroot = args[3].replace('"','')
    reportroot = args[4].replace('"','')
    #ext - image file extension
    im = args[5] #sys.argv[3]
    ext = args[6] #sys.argv[4] #'fts'
    run = args[7] #sys.argv[5]
    date = args[8] #sys.argv[6]
    tel = args[9]

    print(directory, imageroot, im, ext, run, date, tel)

    if DEBUG:
        print("----------------------------","| Main","| ",directory, im, ext, run)
    runnumber = check_log(run)
    result = write_logstart(imageroot, im, ext, run, runnumber)
    #get a list of all fits files:
    liste = get_liste(imageroot, im, ext)

    if DEBUG:
        print("----------------------------","| Main (2)","| ",liste)

    #separate into lists based on whether the image is bias/dark/flat/science
    (biaslist, darklist, flatlist, sciencelist, ditherlist) = sort_liste(liste, run, runnumber)
    write_log(run, runnumber, "# \n", 1)
    write_log(run, runnumber, "# \n", 2)

    fnames, scilists, filters, scifilt = sort_scilist(sciencelist)
    filtnames, flatlists= sort_flatlist(flatlist)
    print("Flat filters: ",filters)
    print("Science Image filters: ", scifilt)

    dir = imageroot

    # for the IR camera don't use bias and darks temporarily until the dimension issue is sorted
    if False:
    #if 'zYJ' in scifilt.values():
        biaslist=[]
        darklist=[]
        print("WARNING: Infrared Camera, not using BIAS or DARK in calibration")
    else:

        #if there are no calibration images found on this night then look for them on previous nights
        if len(biaslist) < 2:
            new_biaslist = []
            new_date = date
            bcount = 0
            while len(new_biaslist)<2:
                bcount = bcount + 1
                if bcount < 100:
                    print("Not enough bias images found for " + str(new_date) + ". Searching for bias images taken on previous dates.")
                else:
                    break
                date_form = (datetime.strptime(new_date, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
                dir = dir.replace(str(new_date),str(date_form))
                new_date = str(date_form)
                if os.path.exists(dir):
                    print("Looking in directory: " + dir)
                    if ext == "fz":
                        liste = get_liste(dir, im, ext)
                        (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste, run,
                                                                                                                 runnumber)
                        if len(new_biaslist)==0:
                            liste = get_liste(dir, im, "fits")
                    else:
                        liste = get_liste(dir, im, ext)
                    (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste, run, runnumber)
            biaslist = new_biaslist
            # darklist = new_darklist

        dir = imageroot

        if len(darklist) == 0:
            new_darklist = []
            new_date = date
            dcount = 0
            while len(new_darklist)==0:
                dcount = dcount + 1
                if dcount < 100:
                    print("No dark images found for " + str(new_date) + ". Searching for dark images taken on previous dates.")
                else:
                    break
                date_form = (datetime.strptime(new_date, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
                dir = dir.replace(str(new_date),str(date_form))
                new_date = str(date_form)
                if os.path.exists(dir):
                    print("Looking in directory: " + dir)
                    if ext == "fz":
                        liste = get_liste(dir, im, ext)
                        (new_biaslist2, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste,
                                                                                                                  run,
                                                                                                                  runnumber)
                        if len(new_darklist)==0:
                            liste = get_liste(dir, im, "fits")
                    else:
                        liste = get_liste(dir, im, ext)
                    (new_biaslist2, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste, run, runnumber)
            # biaslist = new_biaslist
            darklist = new_darklist

    # if there are no flat images found on this night then look for them on previous nights
    # in addition if there are no flat images found for the correct filter look for these on
    # previous nights as well

    # filterlist = []
    

    
    flatdict = {}
    for i in range(len(filtnames)):
        flatdict[filtnames[i]] = flatlists[i]

    ######### TEST REMOVAL #########


    # find if there are any science images with filters which do not have flats
    # for f in filters:
    #     dir = imageroot
    #     print("Searching for flats in " + f)
    #     if f not in filtnames:
    #         # new_flatlist = []
    #         new_date = date
    #         fcount = 0
    #         while f not in flatdict.keys():
    #             fcount = fcount + 1
    #             if fcount<30:
    #                 print("No flat images found for filter " + f + " on " + str(new_date) + ". Searching for flat images taken on previous dates.\n")
    #             date_form = (datetime.strptime(new_date, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
    #             dir = dir.replace(str(new_date),str(date_form))
    #             new_date = str(date_form)
    #             if os.path.exists(dir):
    #                 print("Looking in directory: " + dir)
    #                 if ext == "fz":
    #                     liste = get_liste(dir, im, ext)
    #                     (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste,
    #                                                                                                              run,
    #                                                                                                              runnumber)
    #                     if len(new_flatlist)==0:
    #                         liste = get_liste(dir, im, "fits")
    # 
    #                 else:
    #                     liste = get_liste(dir, im, ext)
    #                 (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste, run, runnumber)
    #                 new_filtnames, new_flatlists = sort_flatlist(new_flatlist)
    #                 for i in range(len(new_filtnames)):
    #                     # print f
    #                     # print new_filtnames
    #                     if new_filtnames[i] == f:
    #                         if new_filtnames[i] not in flatdict.keys():
    #                             flatdict[new_filtnames[i]] = new_flatlists[i]
    #                         else:
    #                             flatdict[new_filtnames[i]].extend(new_flatlists[i])

    ############################

    # find if there are any science images with filters which do not have flats
    for f in filters:
        dir = imageroot
        print("Searching for flats in " + f)
        if f not in filtnames or len(flatdict[f]) < 4:
            
            new_date = date
            fcount = 0
            # Continue looping until both f is in flatdict and it has more than 3 entries
            while f not in list(flatdict.keys()) or len(flatdict[f]) < 4:
                fcount = fcount + 1
                if fcount < 30:
                    print("No sufficient flat images found for filter " + f + " on " + str(
                        new_date) + ". Searching for flat images taken on previous dates.\n")
                date_form = (datetime.strptime(new_date, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
                dir = dir.replace(str(new_date), str(date_form))
                new_date = str(date_form)
                if os.path.exists(dir):
                    print("Looking in directory: " + dir)
                    if ext == "fz":
                        liste = get_liste(dir, im, ext)
                        (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(
                            liste, run, runnumber)
                        if len(new_flatlist) == 0:
                            liste = get_liste(dir, im, "fits")
                    else:
                        liste = get_liste(dir, im, ext)

                    (new_biaslist, new_darklist, new_flatlist, new_sciencelist, new_ditherlist) = sort_liste(liste,
                                                                                                             run,
                                                                                                             runnumber)
                    new_filtnames, new_flatlists = sort_flatlist(new_flatlist)
                    for i in range(len(new_filtnames)):
                        if new_filtnames[i] == f:
                            if new_filtnames[i] not in list(flatdict.keys()):
                                flatdict[new_filtnames[i]] = new_flatlists[i]
                            else:
                                flatdict[new_filtnames[i]].extend(new_flatlists[i])
                
                
                
        # print flatdict

                if fcount > 30:
                    # if tel != "Europa":
                    print("WARNING: No flat images found for " + f + " in the past 30 days \n")
                        # break
                    # else:
                    #     print "No flat images found for " + f + " in the past 30 days \n"
                    #     print "Europa is currently not taking flat images"
                if fcount > (365*3):
                    print("WARNING: No flat images found for " + f + " in the past 3 years \n")
                    break



    #write list of bias image filenames to file "(run)_bias.list"
    print(biaslist)
    result = write_liste(biaslist, "%s_bias.list" % (run), outputroot, run, runnumber, date)
    # write list of dark image filenames to file "(run)_dark.list"
    result = write_liste(darklist, "%s_dark.list" % (run), outputroot, run, runnumber, date)

    max_flat=0
    # write list of flat image filenames to file "(run)_flat.list
    for k,v in list(flatdict.items()):
        result = write_liste(v, "%s_flat_%s.list" % (run,k), outputroot, run, runnumber, date)
        # result = write_liste(flatlist, "%s_flat.list" % (run), imageroot, run, runnumber, date)
        if len(v) > max_flat:
            max_flat = len(v)

    for i, field in enumerate(fnames):
        # if there are less than 5 science images for a field then don't create list file
        if len(scilists[i]) > 4:
            # Convert spaces to double hyphens for filename consistency
            field_filename = field.replace(' ', '--')
            # for each object write list of science image filenames to file "(run)_image_(object).list"
            result = write_liste(scilists[i], "%s_image_%s.list" % (run, field_filename), outputroot, run, runnumber,
                                 date)
        else:
            print("Less than 5 science images for " + field)

        # QC update number of images
        QC([0,reportroot + "/QC", date, field, tel, run, 'QC1', len(scilists[i])])
        QC([0,reportroot + "/QC", date, field, tel, run, 'QC2', len(biaslist)])
        QC([0,reportroot + "/QC", date, field, tel, run, 'QC3', len(darklist)])
        try:
            QC([0,reportroot + "/QC", date, field, tel, run, 'QC4', len(flatdict[scifilt[field]])])
        except:
            QC([0,reportroot + "/QC", date, field, tel, run, 'QC4', 0])
    # fnames, ditlists = sort_scilist(ditherlist)
    # for i, field in enumerate(fnames):
    #     result = write_liste(ditlists[i], "%s_dither_%s.list" % (run, field), imageroot, run, runnumber, date)

    string = "# \n"
    string = string + "%20s Finished Job\n " % time.strftime("%Y-%m-%d %H:%M:%S")
    string = string + "# -----------------------------------------------\n "
    string = string + "\n"
    write_log(run, runnumber, string, 1)
    write_log(run, runnumber, string, 2)


if __name__ == '__main__':
    main(sys.argv)
