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

Calibration frames come from the calibration register (calibration/calregister.py),
which T1 first updates with this night's own frames: for each type the nearest
compatible night on either side (calibration/calselect.py), never another camera,
never several nights combined.  The choice is logged and written to
{runname}_calib_selection.json; a night nothing can calibrate stops here.

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
from reporting.QC import main as QC
# open_fits_file presents a Rice-compressed frame (.fz) as the plain file, its image and header in HDU 0
from calibration.pipeutils import open_fits_file
from calibration import calregister, calselect
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


def classify(files, logroot, runnumber, telescope=''):
    pool = ThreadPool()
    fn = partial(classify_file, logroot=logroot, runnumber=runnumber, telescope=telescope)
    results = pool.map(fn, files)
    keys = ['dithered', 'science', 'dark', 'bias', 'flat', 'native_andor', 'meta']
    pool.close()

    return {key: sum([entry[key] for entry in results], []) for key in keys}


def classify_file(filename, logroot, runnumber, telescope=''):
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
            # what the calibration register and the calibration choice need from the header
            out['meta'].append((filename, calregister.header_fields(hdulist[0].header, telescope)))
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


def sort_liste(liste, logroot, runnumber, telescope='', with_meta=False):
    ### Sort images for bias, darks, flats, and scientific frames"""
    # print liste
    classification = classify(liste, logroot, runnumber, telescope)
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
    if with_meta:
        return (biaslist, darklist, flatlist, science, dithered), dict(classification['meta'])
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

def relative_to(path, root):
    rel = os.path.relpath(path, root)
    return rel.replace(os.sep, '/')


def science_readmodes(paths, telescope, night, n=5):
    """Readout modes of up to n science frames spread through the night (from their overscan strips)."""
    paths = sorted(paths)
    step = max(1, len(paths) // n)
    modes = []
    for p in paths[::step][:n]:
        info = calregister.read_frame(p, telescope, calregister.night_of(night) or night, stats='always')
        if not info.get('error') and info.get('readmode'):
            modes.append(info['readmode'])
    return modes


def dark_matching(telescop, instrument):
    params = calregister._instrument_config().get('telescopes', {}).get(telescop, {}).get(instrument, {})
    return bool(params.get('dark_matching_exptime', False))


def choose_calibrations(tel, date, imageroot, outputroot, run, liste, biaslist, darklist, flatlist, sciencelist,
                        meta):
    """This night's bias, dark and flat lists, from the calibration register.

    The night's own calibration frames are registered first (pixel statistics read here), so they compete with
    every other night on equal terms and win whenever they are enough.  Returns (biaslist, darklist,
    {filter: flats}).  A night without science frames keeps its own frames.  Exits (stopping the run) when
    nothing in the register can calibrate the science frames.
    """
    flatdict_own = {}
    for f in flatlist:
        flatdict_own.setdefault((meta.get(f) or {}).get('filter') or '', []).append(f)
    if not sciencelist:
        print("No science frames: the calibration lists hold this night's own frames only")
        return biaslist, darklist, flatdict_own

    path = calregister.resolve_register_path(imageroot=imageroot)
    try:
        register = calregister.open_existing(path)
    except calregister.RegisterError as e:
        sys.exit("ERROR: T1 cannot choose calibrations for {} {}: {}. Build it with "
                 "'python -m calibration.calregister build' or point {} at one ({}=1 to leave it "
                 "unchanged).".format(tel, date, e, calregister.REGISTER_ENV, calregister.READONLY_ENV))
    images_root = os.path.dirname(os.path.abspath(imageroot))
    night_path = os.path.abspath(imageroot)
    cal = biaslist + darklist + flatlist
    science = [dict(meta[f], relpath=relative_to(f, night_path)) for f in sciencelist if f in meta]
    try:
        night_row, rows = calregister.register_rows_from_t1(
            register, tel, date, night_path, [relative_to(f, night_path) for f in cal], science, len(liste))
        print("Calibration register {}: {} {} {} bias, {} dark, {} flat frames of this night{}".format(
            path, tel, date, night_row['n_bias'], night_row['n_dark'], night_row['n_flat'],
            ' (read-only: not saved)' if register.readonly else ' registered'))
        profile = calselect.science_profile(science)
        majority = [f for f in sciencelist if (meta.get(f) or {}).get('config') == (profile or {}).get('config')]
        modes = science_readmodes(majority, tel, date)
        telescop = (profile or {}).get('telescop')
        instrument = (profile or {}).get('instrument')
        record = calselect.select(register, tel, date, images_root, rows, science, modes,
                                  dark_matching=dark_matching(telescop, instrument),
                                  tonight_detector=night_row['detector'])
    finally:
        register.close()

    print("Calibrations chosen for {} {} (science: {}):".format(
        tel, date, (record.get('profile') or {}).get('config')))
    for line in calselect.describe(record):
        print(line)
    reddir = os.path.join(outputroot, date, 'reduction')
    try:
        print("Calibration choice written to " + calselect.write_record(record, reddir, run))
    except OSError as e:
        print("WARNING: could not write the calibration choice: " + str(e))
    if record['errors']:
        sys.exit("ERROR: T1 cannot calibrate {} {}: {}".format(tel, date, '; '.join(record['errors'])))

    biaslist = list(record['bias']['frames'])
    darklist = list(record['dark']['frames'])
    flatdict = {f: list(e['frames']) for f, e in record['flat'].items()}
    # this night's usable flats in filters no science frame needs are kept, as before (T5 may use them)
    usable = {os.path.join(night_path, r['relpath']) for r in rows if r['kind'] == 'flat' and r['usable']}
    for f, frames in flatdict_own.items():
        if f not in flatdict:
            mine = [x for x in frames if x in usable]
            if mine:
                flatdict[f] = mine
    return biaslist, darklist, flatdict


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
    (biaslist, darklist, flatlist, sciencelist, ditherlist), meta = sort_liste(liste, run, runnumber, tel,
                                                                               with_meta=True)
    write_log(run, runnumber, "# \n", 1)
    write_log(run, runnumber, "# \n", 2)

    fnames, scilists, filters, scifilt = sort_scilist(sciencelist)
    print("Science Image filters: ", scifilt)

    # calibrations are chosen for the science frames T6 will reduce (targets with an image list)
    reduced = [f for frames in scilists if len(frames) > 4 for f in frames]
    biaslist, darklist, flatdict = choose_calibrations(tel, date, imageroot, outputroot, run, liste, biaslist,
                                                       darklist, flatlist, reduced, meta)

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
            QC([0,reportroot + "/QC", date, field, tel, run, 'QC4', len(flatdict[scifilt[field].replace("'", "")])])
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
