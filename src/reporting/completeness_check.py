import sys
import os
from astropy.io import ascii
from astropy.table import Table,Column
# from tabulate import tabulate
import argparse
import glob
import numpy as np


def create_file(fname,dir,all_nights,version):
    dates = all_nights['date']
    tel = all_nights['tel']
    targets = all_nights['targets']
    s_images = all_nights['science_images_exist']
    c_images = all_nights['calibration_images_exist']
    images = all_nights['images_exist']
    pipeline = all_nights['pipeline_started']
    outputfits = all_nights['output_file']
    lcurves = all_nights['lightcurves']
    r = all_nights['report']
    notes = all_nights['notes']
    all_im = all_nights['all_images']

    # if the report file does not exist then initialise it
    # row = [date,target,'[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]','[o]']
    # row = []
    # create a test row to set the width of the columns when reading in the table
    test_row = ['--------','---------','--------------------------------------------------','------------','----------------','-----------','-------','----------','----------------','-----------','-------','------------------------------']
    # column headers:
    hdr = ['Night','Telescope','Targets','Images Exist','Pipeline ('+version+') Started','Calibration','Science','All Images','Pipeline ('+version+') Completed','Lightcurves','Reports','Notes']
    # column data types:
    dtypes = ['U8','U9','U50','U12','U11','U7','U16','U10','U16','U11','U7','U30']
    # Under Python 2 numpy cut longer values to these widths (e.g. a night's Targets to 50 characters, and the
    # test row's 16 and 11 dashes under Started and Calibration). astropy 5 widens the column instead, which
    # would change the report, so cut them here.
    widths = [int(d[1:]) for d in dtypes]

    # create the table, add the test row and first row to the table and write to report file
    report = Table(names=hdr,dtype=dtypes)
    report.add_row([v[:w] for v,w in zip(test_row,widths)])

    for i in range(int(len(dates))):
        row = [dates[i],tel[i]," ".join(targets[i]),images[i],pipeline[i],c_images[i],s_images[i],all_im[i]," ".join(outputfits[i])," ".join(lcurves[i])," ".join(r[i]),notes[i]]
        row = [str(v)[:w] for v,w in zip(row,widths)]
        try:
            report.add_row(row)
        except Exception as e:
            print(e)
    ascii.write(report,dir + "/" + fname,format='fixed_width',overwrite=True)

def read_images(flist):
    with open(flist) as infile:
        #loop through all science fits filenames
        filenames = [line.strip() for line in infile]
    return filenames

def num_calib_images(fdir,datestring):
    flatimages = glob.glob("%s/*_flat_*.list" % (fdir))
    biasimages = glob.glob("%s/*_bias.list" % (fdir))
    darkimages = glob.glob("%s/*_dark.list" % (fdir))

    calib_images = 0
    if len(biasimages) > 0:
        bias = read_images(biasimages[0])
        calib_images = calib_images + len([1 for b in bias if datestring in b])
    if len(darkimages) > 0:
        dark = read_images(darkimages[0])
        calib_images = calib_images + len([1 for b in dark if datestring in b])
    if len(flatimages) > 0:
        for f in flatimages:
            flat = read_images(f)
            calib_images = calib_images + len([1 for b in flat if datestring in b])

    return calib_images


def num_science_images(fdir):
    scienceimages = glob.glob("%s/*_image_*.list" % (fdir))
    numfiles = 0

    if len(scienceimages)>0:
        for i in scienceimages:
            numfiles = numfiles + len(read_images(i))

    return numfiles

def completeness_check(dir,version):
    all_nights =  {'date':[],'tel':[],'targets':[],'images_exist':[],'science_images_exist':[],'calibration_images_exist':[],'all_images':[],'pipeline_started':[],'output_file':[],'lightcurves':[],'report':[],'notes':[]}

    tels = ['Europa','Io','Callisto', 'Ganymede','Artemis']
    sdates = ['20170318','20170617','20180121','20180930','20190426']
    obsdir = "Observations"
    datdir = "PipelineOutput"
    europa_im = os.listdir(dir + "/" + obsdir + "/Europa/images/")
    io_im = os.listdir(dir + "/" + obsdir + "/Io/images/")
    callisto_im = os.listdir(dir + "/" + obsdir + "/Callisto/images/")
    ganymede_im = os.listdir(dir + "/" + obsdir + "/Ganymede/images/")
    artemis_im = os.listdir(dir + "/" + obsdir + "/Artemis/images/")

    europa_out = os.listdir(dir + "/" + datdir + "/" + version + "/Europa/output/")
    io_out = os.listdir(dir + "/" + datdir + "/" + version + "/Io/output/")
    callisto_out = os.listdir(dir + "/" + datdir + "/" + version + "/Callisto/output/")
    ganymede_out = os.listdir(dir + "/" + datdir + "/" + version + "/Ganymede/output/")
    artemis_out = os.listdir(dir + "/" + datdir + "/" + version + "/Artemis/output/")
    # print(callisto_out)

    dates_notes = {'20201008_Artemis_TOI-2267':"Guiding Issue",
                   '20180609_Europa_SP1201-2737':"Can't fit astrometry (CCD Temp = 12C)",
                   '20190904_Europa_GAIJIN-1':"Can't fit astrometry (no stars)",
                   '20181115_Europa_SP0024-0158': "Can't fit astrometry (dome closed?)",
                   '20190104_Europa_SP0216-2200':"Bad Images (moon?)",
                   '20190105_Europa_SP0216-2200': "Bad Images (moon?)",
                   '20201004_Ganymede_BD-16-52-':"Gaia DR1 ID",
                   '20201030_Ganymede_HD-14326-':"Change in defocus",
                   '20210203_Artemis':"Guiding Issue"} #_TOI-2096.01.02
    print(dates_notes.keys())

    ims = [europa_im,io_im,callisto_im,ganymede_im,artemis_im]
    outs = [europa_out,io_out,callisto_out,ganymede_out,artemis_out]


    from datetime import date, timedelta, datetime

    sdate = date(2017,3,18)
    # edate = date(2019,10,2)
    edate = date.today()

    delta = edate - sdate

    # all_nights['notes'] = [" " for r in range(delta.days + 1)]

    # print(len(all_nights['notes']))


    for i in range(delta.days + 1):
        day = sdate + timedelta(days=i)
        datestring = day.strftime('%Y%m%d')
        print(datestring)

        for t in range(len(ims)):
            if datestring in ims[t]:
                # note=0
                if os.path.exists(dir + "/" + datdir + "/" + version + "/" + tels[t] + "/reports/" + tels[t] + "_" + datestring + ".pdf"):
                    all_nights['report'].append('1')
                elif os.path.exists(dir + "/" + datdir + "/" + version + "/" + tels[t] + "/reports/" + tels[t] + "_" + datestring +  "_"  + version +".pdf"):
                    all_nights['report'].append('1')
                else:
                    all_nights['report'].append('0')
                all_nights['date'].append(datestring)
                all_nights['tel'].append(tels[t])

                all_nights['images_exist'].append(1)

                science_images = num_science_images(dir + "/" + datdir + "/" + version + "/" + tels[t] + "/output/" + datestring + "/reduction/")
                calib_images = num_calib_images(dir + "/" + datdir + "/" + version + "/" + tels[t] + "/output/" + datestring + "/reduction/",datestring)

                # scienceimages = glob.glob("%s/%s/output/%s/%s/reduction/*_image_*.list" % (dir, tels[t], version, datestring))

                if science_images>0:
                    all_nights['science_images_exist'].append(1)
                else:
                    all_nights['science_images_exist'].append(0)
                    # if calib_images>0:
                        # all_nights['notes'][i] = "no science images"

                    # all_nights['notes'].append("no science images")
                    # note=1

                if calib_images>0:
                    all_nights['calibration_images_exist'].append(1)
                else:
                    all_nights['calibration_images_exist'].append(0)


                total_processed_images = calib_images + science_images
                # raw frames stored plain or Rice-compressed (.fz, utils/framepack.py), counted once each
                night = dir + "/" + obsdir + "/" + tels[t] + "/images/" + datestring
                frames = set()
                for pattern in ("/*/*.fits", "/*/*.fts", "/*.fits", "/*.fts"):
                    for f in glob.glob(night + pattern) + glob.glob(night + pattern + ".fz"):
                        frames.add(f[:-3] if f.endswith(".fz") else f)
                total_images = len(frames)
                print(tels[t],total_images,total_processed_images)

                if (total_processed_images == total_images) or (total_processed_images == total_images-1) :
                    all_nights['all_images'].append(1)
                else:
                    all_nights['all_images'].append(0)

                if datestring in outs[t]:
                    all_nights['pipeline_started'].append(1)
                    # targets = os.listdir(dir+"/"+tels[t]+"/output/" + version + "/" +datestring)
                    targets = glob.glob("%s/%s/%s/%s/output/%s/*/*_output.fits"%(dir,datdir,version,tels[t],datestring))
                    # targets.remove('reduction')
                    if len(targets)>0:
                        sptargets = [targ_temp.split("/")[-2] for targ_temp in targets]
                        targets = [os.path.basename(targ_temp).split("_")[0] for targ_temp in targets]
                        # print(targets)
                        # print(sptargets)
                        all_nights['targets'].append(tuple(sptargets))
                        of,lc,notes = [],[],[]
                        for targ in range(len(targets)):

                            if datestring + "_" + tels[t] + "_" + sptargets[targ] in dates_notes.keys():
                                print(datestring + "_" + tels[t] + "_" + sptargets[targ])
                                print(dates_notes[datestring + "_" + tels[t] + "_" + sptargets[targ]])
                                notes.append(dates_notes[datestring + "_" + tels[t] + "_" + sptargets[targ]])
                                # note=1
                            else:
                                notes.append(' ')

                            # if targ == "Gaijin-1":
                            #     print "test"
                            # fname = "%s/*_output.fits"%(dir + "/" + tels[t] + "/output/" + version + "/" + datestring + "/" + targ)
                            fname =  glob.glob("%s/%s/%s/%s/output/%s/*/%s_*_output.fits"%(dir,datdir,version,tels[t],datestring,targets[targ]))
                            # print(fname)
                            if len(fname)>0:
                                of.append('1')

                                # lcurves = glob.glob("%s/lightcurves/*/*global.csv" % (dir + "/" + tels[t] + "/output/" + version + "/" + datestring + "/" + targ))
                                # lcurves = glob.glob("%s/*/lightcurves/%s_*_%s_*_MCMC" % (
                                #             dir + "/" + tels[t] + "/output/" + version + "/" + datestring,targ,datestring))
                                lcurves = glob.glob("%s/%s/%s/%s/output/%s/%s/*_diff.fits" % (dir,datdir,version,tels[t],datestring,sptargets[targ]))

                                if len(lcurves)>0:
                                    lc.append('1')
                                else:
                                    lc.append('0')
                            else:
                                of.append('0')
                                lc.append('0')

                        all_nights['output_file'].append(tuple(of))
                        all_nights['lightcurves'].append(tuple(lc))
                        # if note==0:
                        all_nights['notes'].append(" ".join(notes))

                        # print(len(all_nights['lightcurves']))
                        # print(len(all_nights['notes']))
                    else:
                        all_nights['targets'].append('N/A')
                        all_nights['output_file'].append('0')
                        all_nights['lightcurves'].append('0')
                        all_nights['notes'].append(' ')

                        if datestring + "_" + tels[t] in dates_notes.keys():
                            print(datestring + "_" + tels[t])
                            print(dates_notes[datestring + "_" + tels[t]])
                            notes.append(dates_notes[datestring + "_" + tels[t]])
                            # note=1
                        else:
                            notes.append(' ')

                else:
                    all_nights['targets'].append('N/A')
                    all_nights['pipeline_started'].append('0')
                    all_nights['output_file'].append('0')
                    all_nights['lightcurves'].append('0')
                    # if note == 0:
                    if science_images == 0 and calib_images>0:
                        all_nights['notes'].append('no science images')
                    else:
                        all_nights['notes'].append(' ')
            else:
                if int(datestring) > int(sdates[t]):
                    all_nights['date'].append(datestring)
                    all_nights['tel'].append(tels[t])
                    all_nights['targets'].append(' ')
                    all_nights['science_images_exist'].append(' ')
                    all_nights['calibration_images_exist'].append(' ')
                    all_nights['images_exist'].append(' ')
                    all_nights['pipeline_started'].append(' ')
                    all_nights['output_file'].append(' ')
                    all_nights['lightcurves'].append(' ')
                    all_nights['report'].append(' ')
                    all_nights['all_images'].append(' ')
                    all_nights['notes'].append(' ')
                else:
                    print("Before " + tels[t] + " started observations.")

    return all_nights

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fname", help="Filename of this report", required=True)
    parser.add_argument("-d", '--dir', help='Directory of Telescope Folders', required=True)
    parser.add_argument("-v",'--version',help='Pipeline Version', required=True)
    args = parser.parse_args()

    all_nights = completeness_check(args.dir,args.version)
    create_file(args.fname, args.dir,all_nights,args.version)