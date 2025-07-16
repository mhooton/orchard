import argparse
import glob
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
import matplotlib.gridspec as gridspec
import os
import csv
from functools import partial
from scipy.ndimage import interpolation as interp
from skimage.registration import phase_cross_correlation
import astropy.io.fits as pyfits
import astropy.visualization as pyvis
from astropy.io import ascii
import datetime
from matplotlib.font_manager import FontProperties
import fitsio
from astropy.time import Time
import math
matplotlib.use('pdf')
import warnings

def get_chronological_order(jdstart):
    order = np.argsort(jdstart)
    return order


def import_spids(fname, args):
    target_id = []
    s = (args.target).split(" ")
    if os.path.exists(fname):
        data = ascii.read(fname)
        target_list = [d.upper() for d in data['Name']]
        sp_id_list = [str(x) for x in data['SP_id']]
        for i in s:
            try:
                where_target = np.where(np.array(target_list) == i.upper())[0][0]
                targ_id = sp_id_list[where_target]
                target_id.append(targ_id)
            except:
                target_id.append('unknown')
    else:
        for i in s:
            target_id.append('unknown')
    return target_id


def import_masterframes_img(args, filt):
    reductdir = args.datdir + "/output/" + args.date + "/reduction"

    masterbias_im = glob.glob("%s/*_Master%s.*" % (reductdir, "Bias"))
    masterdark_im = glob.glob("%s/*_Master%s.*" % (reductdir, "Dark"))

    masterbias = args.datdir + "/reports/temp/masterbias.png"
    masterdark = args.datdir + "/reports/temp/masterdark.png"
    # masterbias = args.datdir + "/reports/temp_test/masterbias.png"
    # masterdark = args.datdir + "/reports/temp_test/masterdark.png"

    if (len(masterbias_im) != 0):
        with fits.open(masterbias_im[0]) as mb:
            data = mb[0].data
            plt.imshow(data, interpolation='nearest', cmap='gray')
            plt.clim(0, 10)
            plt.colorbar()
            plt.title("MasterBias")
            plt.savefig(masterbias, bbox_inches='tight')
            plt.close()
    if (len(masterdark_im) != 0):
        with fits.open(masterdark_im[0]) as md:
            data = md[0].data
            plt.imshow(data, interpolation='nearest', cmap='gray')
            plt.clim(0, 1)
            plt.colorbar()
            plt.title("MasterDark")
            plt.savefig(masterdark, bbox_inches='tight')
            plt.close()

    f_old = 'dum'
    masterflat = []

    for f in filt:
        if f != f_old:
            masterflat_name = args.datdir + "/reports/temp/masterflat_" + f + ".png"
            # masterflat_name = args.datdir + "/reports/temp_test/masterflat_" + f + ".png"
            masterflat.append(masterflat_name)
            f_old = f
            masterflat_im = glob.glob("%s/*_Master%s_%s.*" % (reductdir, "Flat", f))
            if (len(masterflat_im) != 0):
                with fits.open(masterflat_im[0]) as mf:
                    data = mf[0].data
                    plt.imshow(data, interpolation='nearest', cmap='gray')
                    plt.clim(0.95, 1.05)
                    plt.colorbar()
                    plt.title("MasterFlat " + f)
                    plt.savefig(masterflat_name, bbox_inches='tight')
                    plt.close()

    # masterflat = []
    # for f in masterflats:
    #     print f
    #     sep_ext = f.split(".")
    #     filt = sep_ext[0].split('_')[-1]
    #     with fits.open(f) as mf:
    #         data = mf[0].data
    #         plt.imshow(data, interpolation='nearest', cmap='gray')
    #         plt.clim(0.95, 1.05)
    #         plt.colorbar()
    #         plt.title("MasterFlat " + str(filt))
    #         masterflat.append(args.outdir + "/" + args.telescope + "/reports/tests/masterflat_" + filt + ".png")
    #         plt.savefig(masterflat[-1])
    #         plt.close()

    return masterbias, masterdark, masterflat


def import_night_img(args):
    ramove = []
    decmove = []
    airmass = []
    fwhm = []
    skybkg = []
    altitude = []
    jdstart = []
    jdend = []
    nimages = []
    filt = []
    ellip = []
    target_x_pos = []
    target_y_pos = []
    gaia = []
    targs = []

    label_names = ['FWHM [pix]', 'BACKGROUND [ADU/pix]', 'RA MOVE [arcsec]', 'DEC MOVE [arcsec]', 'AIRMASS',
                   'ALTITUDE [deg]']
    names = ['FWHM', 'BKG SKYLEVEL', 'RA MOVE', 'DEC MOVE', 'AIRMASS', 'ALTITUDE']

    s = (args.target).split(" ")
    # targ_count = 0
    for i in s:
        tardir = args.datdir + "/output/" + args.date + "/" + i

        outputfits = glob.glob("%s/*_output.fits" % (tardir))
        if (len(outputfits) != 0):
            for fitsfile in outputfits:
                gaia_id = fitsfile.split("/")[-1].split("_")[0]
                gaia.append(gaia_id)
                targs.append(i)

                for a in range(len(names)):
                    outfig = args.datdir + "/reports/temp/" + names[a] + "_" + gaia_id + ".png"
                    # outfig = args.datdir + "/reports/temp_test/" + names[a] + "_" + i + ".png"
                    if (a == 0):
                        fwhm.append(outfig)
                    elif (a == 1):
                        skybkg.append(outfig)
                    elif (a == 2):
                        ramove.append(outfig)
                    elif (a == 3):
                        decmove.append(outfig)
                    elif (a == 4):
                        airmass.append(outfig)
                    elif (a == 5):
                        altitude.append(outfig)

                # fitsfile = glob.glob("%s/*_output.fits" % (tardir))[0]
                with fitsio.FITS(fitsfile) as infile:
                    imagelist = infile['imagelist']
                    airmass_read = imagelist['airmass'].read()
                    altitude_read = imagelist['altitude'].read()
                    ramove_read = imagelist['ra_move'].read()
                    decmove_read = imagelist['dec_move'].read()
                    fwhm_read = imagelist['fwhm'].read()
                    skybkg_read = imagelist['skylevel'].read()
                    filt_read = imagelist['filter'].read()
                    psf_a_5_read = imagelist['psf_a_5'].read()
                    psf_b_5_read = imagelist['psf_b_5'].read()
                    psf_t_5_read = imagelist['psf_t_5'].read()
                    cat = infile['catalogue']
                    obj_id = cat['obj_id'].read()
                    jd = infile['jd'].read()
                    ccdx = infile['ccdx'].read()
                    ccdy = infile['ccdy'].read()

                # get all the dates for this target
                with fits.open(fitsfile) as f:
                    # dates = f[0].header['HISTORY']
                    try:
                        sp_id = f[0].header['SP_ID']
                        # teff_header = f[0].header['TEFF']
                        targetid = np.where(np.array(obj_id)==sp_id)[0][0]
                        target_x_pos.append(ccdx[targetid, -1])
                        target_y_pos.append(ccdy[targetid, -1])
                    except:
                        target_x_pos.append('unknown')
                        target_y_pos.append('unknown')
                        pass

                # if target_ids[targ_count] != 'unknown':
                #     targetid = int(target_ids[targ_count][2:])
                #     target_x_pos.append(ccdx[targetid, -1])
                #     target_y_pos.append(ccdy[targetid, -1])
                # else:
                #     target_x_pos.append('unknown')
                #     target_y_pos.append('unknown')
                # targ_count += 1

                ellip_indiv = []
                time_ellip_indiv = []
                theta = []
                for j in range(len(psf_a_5_read)):
                    try:
                        ell = math.sqrt(1 - np.power((psf_b_5_read[j] / psf_a_5_read[j]), 2))
                        ellip_indiv.append(ell)
                        time_ellip_indiv.append(jd[0][j])
                        theta.append(psf_t_5_read[j])
                    except:
                        ()
                ellip.append(np.median(ellip_indiv))

                inputs = [fwhm_read, skybkg_read, ramove_read, decmove_read, airmass_read, altitude_read]

                for a in range(len(inputs)):
                    if (a == 2) or (a == 3):
                        if (np.min(inputs[a]) > -2.0) and (np.max(inputs[a]) < 2.0):
                            plt.figure(dpi=300, figsize=(10, 10))
                            plt.rc('font', size=20)
                            plt.plot(jd[0], inputs[a], 'r.')
                            plt.ylim(-2, 2)
                            plt.xlabel('JD')
                            plt.ylabel(label_names[a])
                            outfig = args.datdir + "/reports/temp/" + names[a] + "_" + gaia_id + ".png"
                            # outfig = args.datdir + "/reports/temp_test/" + names[a] + "_" + i + ".png"
                            plt.savefig(outfig, bbox_inches='tight')
                            plt.close()
                        else:
                            plt.figure(dpi=300, figsize=(10, 10))
                            plt.rc('font', size=20)
                            plt.plot(jd[0], inputs[a], 'r.')
                            plt.xlabel('JD')
                            plt.ylabel(label_names[a])
                            outfig = args.datdir + "/reports/temp/" + names[a] + "_" + gaia_id + ".png"
                            # outfig = args.datdir + "/reports/temp_test/" + names[a] + "_" + i + ".png"
                            plt.savefig(outfig, bbox_inches='tight')
                            plt.close()
                    else:
                        plt.figure(dpi=300, figsize=(10, 10))
                        plt.rc('font', size=20)
                        plt.plot(jd[0], inputs[a], 'r.')
                        plt.xlabel('JD')
                        plt.ylabel(label_names[a])
                        outfig = args.datdir + "/reports/temp/" + names[a] + "_" + gaia_id + ".png"
                        # outfig = args.datdir + "/reports/temp_test/" + names[a] + "_" + i + ".png"
                        plt.savefig(outfig, bbox_inches='tight')
                        plt.close()

                mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
                outfig = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_PSF_ellipticity_" + gaia_id + ".png"
                if not os.path.exists(args.obsdir + "/technical_logs/" + mydate + "/"):
                    os.mkdir(args.obsdir + "/technical_logs/" + mydate)
                plt.figure(dpi=300, figsize=(10, 10))
                plt.rc('font', size=20)
                plt.plot(time_ellip_indiv, ellip_indiv)
                plt.xlabel('JD')
                plt.ylabel('PSF ellipticity')
                plt.savefig(outfig, bbox_inches='tight')
                plt.close()

                outfig = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_PSF_theta_" + gaia_id + ".png"
                plt.figure(dpi=300, figsize=(10, 10))
                plt.rc('font', size=20)
                plt.plot(time_ellip_indiv, theta)
                plt.xlabel('JD')
                plt.ylabel('PSF theta [deg]')
                plt.savefig(outfig, bbox_inches='tight')
                plt.close()

                jdstart.append(jd[0][0])
                jdend.append(jd[0][-1])
                nimages.append(len(jd[0]))
                filt.append(filt_read[0])

        else:
            jdstart.append(0)
            jdend.append(0)
            nimages.append(0)
            filt.append('unknown')
            # targ_count += 1
            target_x_pos.append('unknown')
            target_y_pos.append('unknown')

    return ramove, decmove, airmass, fwhm, skybkg, altitude, jdstart, jdend, nimages, filt, ellip, target_x_pos, target_y_pos, targs, gaia

def import_lcs_new(args,gaia):
    # s = (gaia).split(" ")
    dlc_arr = []
    mlc = []
    complc = []
    rmsflux = []
    newap = []

    for i in gaia:
        lcfits = glob.glob("%s/output/%s/*/%s_*_diff.fits" % (args.datdir,args.date,i))
        bestap_pwv_file = glob.glob(
            "%s/output/%s/*/lightcurves/%s_*_bestap_pwv.txt" % (args.datdir, args.date, i))
        bestap_file = glob.glob("%s/output/%s/*/lightcurves/%s_*_bestap.txt" % (args.datdir,args.date,i))
        if (len(bestap_pwv_file) != 0):
            # Read in the file
            with open(bestap_pwv_file[0], 'r') as file:
                filedata = file.read()
            bestap = filedata.split()[0]
            # ap_err = float(filedata.split()[1])
            lcplot = glob.glob("%s/output/%s/*/lightcurves/%s_*_%s_pwv.png" % (args.datdir,args.date,i,bestap))
            if (len(lcplot)!=0):
                mlc.append(lcplot[0])
            else:
                mlc.append('fake_path')
            # if (len(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestappwv))) != 0):
            #     rmsflux.append(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestappwv))[0])
            # else:
            rmsflux.append('fake_path')
        elif (len(bestap_file) != 0):
            with open(bestap_file[0], 'r') as file:
                filedata = file.read()
            bestap = filedata.split()[0]
            # ap_err = float(filedata.split()[1])
            lcplot = glob.glob("%s/output/%s/*/lightcurves/%s_*_%s.png" % (args.datdir,args.date,i,bestap))
            if (len(lcplot)!=0):
                mlc.append(lcplot[0])
            else:
                mlc.append('fake_path')
            rmsflux.append('fake_path')
        else:
            bestap = str(args.ap)
            lcplot = glob.glob("%s/output/%s/*/lightcurves/%s_*_%s.png" % (args.datdir,args.date,i,bestap))
            if (len(lcplot)!=0):
                mlc.append(lcplot[0])
            else:
                mlc.append('fake_path')
            rmsflux.append('fake_path')

        newap.append(bestap)
        # if (len(glob.glob("%s/comparison_stars/*_%s_%s_v1_field.png" % (lcdir, args.ap, args.date))) != 0):
        #     complc.append(glob.glob("%s/comparison_stars/*_%s_%s_v1_field.png" % (lcdir, args.ap, args.date))[0])
        # else:
        complc.append('fake_path')

        for a in range(3, 9):
            lcplot_pwv = glob.glob(
                "%s/output/%s/*/lightcurves/%s_*_%s_pwv.png" % (args.datdir, args.date, i, a))
            lcplot = glob.glob(
                "%s/output/%s/*/lightcurves/%s_*_%s.png" % (args.datdir, args.date, i, a))

            if (len(lcplot_pwv) != 0):
                dlc_arr.append(lcplot_pwv[0])
            elif (len(lcplot) != 0):
                dlc_arr.append(lcplot[0])
            else:
                dlc_arr.append('fake_path')

    return dlc_arr, mlc, complc, rmsflux, newap


def import_lcs(args):
    dlc_arr = []
    mlc = []
    complc = []
    rmsflux = []
    newap = []

    s = (args.target).split(" ")
    for i in s:
        lcdir = args.datdir + "/output/" + args.date + "/" + i + "/lightcurves"
        bestap_pwv_file = glob.glob("%s/*_bestap_pwv.txt" % (lcdir))
        bestap_file = glob.glob("%s/*_bestap.txt" % (lcdir))
        if (len(bestap_pwv_file) != 0):
            with open(bestap_pwv_file[0], 'r') as f:
                first_line = f.readline().split()
                bestappwv = str(first_line[0])
                newaper = bestappwv
            if (len(glob.glob("%s/*_%s_%s_pwv_3plot.png" % (lcdir, args.date, bestappwv))) != 0):
                mlc.append(glob.glob("%s/*_%s_%s_pwv_3plot.png" % (lcdir, args.date, bestappwv))[0])
            else:
                mlc.append('fake_path')
            if (len(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestappwv))) != 0):
                rmsflux.append(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestappwv))[0])
            else:
                rmsflux.append('fake_path')
        elif (len(bestap_file) != 0):
            with open(bestap_file[0], 'r') as f:
                first_line = f.readline().split()
                bestap = str(first_line[0])
                newaper = bestap
            if (len(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,bestap))) != 0):
                mlc.append(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,bestap))[0])
            else:
                mlc.append('fake_path')
            if (len(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestap))) != 0):
                rmsflux.append(glob.glob("%s/%s_rms_bin5.png" % (lcdir, bestap))[0])
            else:
                rmsflux.append('fake_path')
        else:
            newaper = str(args.ap)
            if (len(glob.glob("%s/*_%s_%s_pwv_3plot.png" % (lcdir, args.date,args.ap))) != 0):
                mlc.append(glob.glob("%s/*_%s_%s_pwv_3plot.png" % (lcdir, args.date,args.ap))[0])
            elif (len(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,args.ap))) != 0):
                mlc.append(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,args.ap))[0])
            else:
                mlc.append('fake_path')
            if (len(glob.glob("%s/%s_rms_bin5.png" % (lcdir, args.ap))) != 0):
                rmsflux.append(glob.glob("%s/%s_rms_bin5.png" % (lcdir, args.ap))[0])
            else:
                rmsflux.append('fake_path')

        newap.append(newaper)
        if (len(glob.glob("%s/comparison_stars/*_%s_%s_v1_field.png" % (lcdir, args.ap, args.date))) != 0):
            complc.append(glob.glob("%s/comparison_stars/*_%s_%s_v1_field.png" % (lcdir, args.ap, args.date))[0])
        else:
            complc.append('fake_path')

        for a in range(3, 9):
            if (len(glob.glob("%s/*_%s_%s_pwv.png" % (lcdir, args.date,a))) != 0):
                dlc_arr.append(glob.glob("%s/*_%s_%s_pwv.png" % (lcdir, args.date,a))[0])
            elif (len(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,a))) != 0):
                dlc_arr.append(glob.glob("%s/*_%s_%s.png" % (lcdir, args.date,a))[0])
            else:
                dlc_arr.append('fake_path')

    return dlc_arr, mlc, complc, rmsflux, newap


def import_weather(args):
    wea = []

    mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    bolt = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_Boltwood.png"
    wea.append(bolt)
    esoplot = args.datdir + "/reports/temp/asm.png"
    # esoplot = args.datdir + "/reports/temp_test/asm.png"
    start = (mydate + " 12:00:00")
    starttime = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    reference_date = datetime.datetime(1980, 1, 1, 00, 00)
    nr_sec = starttime - reference_date
    asm_input = str(int(float(nr_sec.total_seconds())))
    try:
        os.system(
            "wget -O " + esoplot + " -q \"http://www.eso.org/asm/api/ambientServer?print=" + asm_input + "&site=paranal&footer=true\"")
    except:
        ()
    wea.append(esoplot)
    return wea


def import_donuts(args, jdstart, jdend,gaia):
    mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    offsets = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_offsets.csv"

    HAoff = []
    DECoff = []
    # s = (args.target).split(" ")
    for i in gaia:
        HAoff_targ = args.datdir + "/reports/temp/HAoff_" + i + ".png"
        DECoff_targ = args.datdir + "/reports/temp/DECoff_" + i + ".png"
        # HAoff_targ = args.datdir + "/reports/temp_test/HAoff_" + i + ".png"
        # DECoff_targ = args.datdir + "/reports/temp_test/DECoff_" + i + ".png"
        HAoff.append(HAoff_targ)
        DECoff.append(DECoff_targ)

    jd = []
    ha = []
    dec = []

    if os.path.exists(offsets):
        with open(offsets) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader)
            next(csv_reader)
            next(csv_reader)
            for row in csv_reader:
                jd.append(row[0])
                ha.append(row[4])
                dec.append(row[5])
        csv_file.close()
        jd = np.array(jd).astype(float)
        ha = np.array(ha).astype(float)
        dec = np.array(dec).astype(float)

        targ_count = 0
        for i in gaia:
            if (jdstart[targ_count] != 0 and jdend[targ_count] != 0):
                new_jd = []
                new_ha = []
                new_dec = []
                for jd_val, ha_val, dec_val in zip(jd, ha, dec):
                    if jdstart[targ_count] <= jd_val <= jdend[targ_count]:
                        new_jd.append(jd_val)
                        new_ha.append(ha_val)
                        new_dec.append(dec_val)

                plt.figure(dpi=300, figsize=(10, 10))
                plt.rc('font', size=20)
                plt.plot(new_jd, new_ha, 'r.')
                plt.xlabel('JD')
                plt.ylabel('DONUTS HA offset [arcsec]')
                plt.savefig(HAoff[targ_count], bbox_inches='tight')
                plt.close()

                plt.figure(dpi=300, figsize=(10, 10))
                plt.rc('font', size=20)
                plt.plot(new_jd, new_dec, 'r.')
                plt.xlabel('JD')
                plt.ylabel('DONUTS DEC offset [arcsec]')
                plt.savefig(DECoff[targ_count], bbox_inches='tight')
                plt.rc('font', size=12)
                plt.close()

            targ_count += 1

    return HAoff, DECoff


def import_water_vapor(args, jdstart, jdend, gaia):
    mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    lhatpro = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_ESO_LHATPRO_log.txt"

    water_vapor = []
    # s = (args.target).split(" ")
    for i in gaia:
        water_vapor_targ = args.datdir + "/reports/temp/Water_vapor_" + i + ".png"
        # water_vapor_targ = args.datdir + "/reports/temp_test/Water_vapor_" + i + ".png"
        water_vapor.append(water_vapor_targ)

    jd = []
    wv = []

    if os.path.exists(lhatpro) and jdstart[0] >= 2458550.0:  # LHATPRO files don't have right format before 07.03.2019
        with open(lhatpro, 'r') as f:
            next(f)
            for line in f:
                part = np.array(line.split(','))
                jd.append(part[1])
                wv.append(part[7])
        jd = np.array(jd).astype(float)
        wv = np.array(wv).astype(float)

        targ_count = 0
        for i in gaia:
            if (jdstart[targ_count] != 0 and jdend[targ_count] != 0):
                new_jd = []
                new_wv = []
                for jd_val, wv_val in zip(jd, wv):
                    if jdstart[targ_count] <= jd_val <= jdend[targ_count]:
                        new_jd.append(jd_val)
                        new_wv.append(wv_val)

                plt.figure(dpi=300, figsize=(10, 10))
                plt.rc('font', size=20)
                plt.plot(new_jd, new_wv, 'r.')
                plt.xlabel('JD')
                plt.ylabel('Water Vapor [mm]')
                plt.savefig(water_vapor[targ_count], bbox_inches='tight')
                plt.close()

            targ_count += 1

    return water_vapor


def stack_and_vals(args, liste, outstack, outstack_save, outstack_save_dir, outfig_dx, outfig_dy, target_x_pos,
                   target_y_pos, target_count):
    print('Creating stack')
    datamatrix = []
    pa_vector = []
    overscan_vector = []
    ccdtemp_vector = []
    dx_vector = []
    dy_vector = []
    jd_vector = []
    ron = 'unknown'
    dark = 'unknown'
    vals = []
    j = 0
    count=0
    num_i=0
    # len(liste)
    for i in liste:
        with fits.open(i) as infile:
            data = infile[0].data
            count = count + 1

            if (count >= len(liste) / 2) and (num_i < 30):
                datamatrix.append(data.copy())
                num_i=num_i+1

            hdr = infile[0].header
            if 'OVERSCAN' in infile[0].header:
                overscan_vector.append(infile[0].header['OVERSCAN'])
            if 'PA' in infile[0].header:
                pa_vector.append(infile[0].header['PA'])
            if 'CCD-TEMP' in infile[0].header:
                ccdtemp_vector.append(infile[0].header['CCD-TEMP'])
            if j == 0 and 'RON' in infile[0].header and 'DARKCUR' in infile[0].header:
                ron = infile[0].header['RON']
                dark = infile[0].header['DARKCUR']
                j = 1
            if 'DX' in infile[0].header:
                dx_vector.append(infile[0].header['DX'])
            if 'DY' in infile[0].header:
                dy_vector.append(infile[0].header['DY'])
                if 'JD-OBS' in infile[0].header:
                    jd_vector.append(infile[0].header['JD-OBS'])
            infile.close()
    if len(overscan_vector) != 0:
        overscan = np.ma.median(overscan_vector)
    else:
        overscan = 'unknown'
    if len(pa_vector) != 0:
        pa = np.ma.median(pa_vector)
    else:
        pa = 'unknown'
    if len(ccdtemp_vector) != 0:
        ccdtemp = np.ma.median(ccdtemp_vector)
    else:
        ccdtemp = 'unknown'
    vals = [overscan, ron, dark, ccdtemp, pa]
    # print(vals)
    if (np.min(dx_vector) > -2.0) and (np.max(dx_vector) < 2.0):
        plt.figure(dpi=300, figsize=(10, 10))
        plt.rc('font', size=20)
        plt.plot(jd_vector, dx_vector, 'r.')
        plt.ylim(-2, 2)
        plt.xlabel('JD')
        plt.ylabel('DX [pix]')
        plt.savefig(outfig_dx, bbox_inches='tight')
        plt.close()
    else:
        plt.figure(dpi=300, figsize=(10, 10))
        plt.rc('font', size=20)
        plt.plot(jd_vector, dx_vector, 'r.')
        plt.xlabel('JD')
        plt.ylabel('DX [pix]')
        plt.savefig(outfig_dx, bbox_inches='tight')
        plt.close()
    if (np.min(dy_vector) > -2.0) and (np.max(dy_vector) < 2.0):
        plt.figure(dpi=300, figsize=(10, 10))
        plt.rc('font', size=20)
        plt.plot(jd_vector, dy_vector, 'r.')
        plt.ylim(-2, 2)
        plt.xlabel('JD')
        plt.ylabel('DY [pix]')
        plt.savefig(outfig_dy, bbox_inches='tight')
        plt.close()
    else:
        plt.figure(dpi=300, figsize=(10, 10))
        plt.rc('font', size=20)
        plt.plot(jd_vector, dy_vector, 'r.')
        plt.xlabel('JD')
        plt.ylabel('DY [pix]')
        plt.savefig(outfig_dy, bbox_inches='tight')
        plt.close()
    stack = np.ma.median(datamatrix, axis=0)
    stack = np.ma.filled(stack, fill_value=None)
    if os.path.exists(outstack):
        os.remove(outstack)
    if not os.path.exists(outstack_save_dir):
        os.makedirs(outstack_save_dir)
    if os.path.exists(outstack_save):  # to be commented for tests
        os.remove(outstack_save)  # to be commented for tests
    with fits.open(liste[-1]) as infile:
        head = infile[0].header
        fits.writeto(outstack, stack, head)
        fits.writeto(outstack_save, stack, head)  # to be commented for tests
    stackim = outstack[:-5] + '.png'
    if os.path.exists(stackim):
        os.remove(stackim)
    stack_data = pyfits.getdata(outstack)
    norm = pyvis.ImageNormalize(stack_data, interval=pyvis.ZScaleInterval(), stretch=pyvis.SqrtStretch())
    plt.imshow(stack_data, cmap='gray', norm=norm, origin='lower')
    if target_x_pos[target_count] != 'unknown' and target_y_pos[target_count] != 'unknown':
        plt.scatter(target_x_pos[target_count], target_y_pos[target_count], s=70, facecolors='none', edgecolors='b')
    plt.title("Stack image")
    plt.savefig(stackim, bbox_inches='tight')
    plt.close()
    # Cleaning
    for i in liste:
        if os.path.exists(i):
            os.remove(i)
    if os.path.exists(outstack):
        os.remove(outstack)
    # stackim = 'fake_path'

    return stackim, vals


def align_image(image_args):
    """
    Fixed version that takes a tuple of arguments for multiprocessing compatibility
    """
    image2, image1, datdir = image_args

    # alignment of image2 over image1
    im1 = fits.open(image1)[0].data
    im2 = fits.open(image2)[0].data
    im2head = fits.open(image2)[0].header

    # Suppress the overflow warnings from phase_cross_correlation
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message="overflow encountered in scalar multiply")

        try:
            # New API (scikit-image >= 0.18)
            shift, error, diffphase = phase_cross_correlation(im1, im2, upsample_factor=100)
            result = shift
        except Exception as e:
            print(f"Warning: Phase cross correlation failed for {image2}: {e}")
            # Fallback to no shift if registration fails
            result = [0.0, 0.0]

    # Rest of the function
    shifted_im2 = interp.shift(im2, result)
    path, filename = os.path.split(image2)
    outname = datdir + "/reports/temp/" + filename[:-5] + '_aligned.fits'

    if os.path.exists(outname):
        os.remove(outname)
    im2head.add_history('Shifts compared to ' + image1 + ' are ' + str(result))
    im2head['DX'] = (result[0], 'DX')
    im2head['DY'] = (result[1], 'DY')
    fits.writeto(outname, shifted_im2, im2head)


def align(liste, datdir):
    """
    Updated align function to pass datdir to align_image
    """
    print('Aligning images')
    # parallelization for the alignment
    import multiprocessing as mp

    print("Number of processors: ", mp.cpu_count(), ". Using ", str(min(20, mp.cpu_count())), " cores")
    pool = mp.Pool(min(20, mp.cpu_count()))  # Don't use more cores than available

    # Create argument tuples for each image (image2, image1, datdir)
    image_args = [(img, liste[-1], datdir) for img in liste]

    pool.map(align_image, image_args)
    pool.close()
    pool.join()  # Wait for all processes to complete

    liste2 = []
    for i in liste:
        path, filename = os.path.split(i)
        outname = datdir + "/reports/temp/" + filename[:-5] + '_aligned.fits'
        liste2.append(outname)
    liste2.sort()

    return liste2

def import_stack_and_vals(args, target_x_pos, target_y_pos, gaia, targs):
    """
    Updated function to pass args.datdir to align function
    """
    stacks = []
    vals = []
    dx = []
    dy = []

    target_count = 0
    for i in targs:
        print(i)
        procdir = args.datdir + "/output/" + args.date + "/" + i + "/1/"
        outstack = args.datdir + "/reports/temp/" + i + "_stack.fits"
        outstack_save = args.datdir + "/output/StackImages_nightly/" + args.date + "/" + i + "_stack.fits"
        outstack_save_dir = args.datdir + "/output/StackImages_nightly/" + args.date
        outfig_dx = args.datdir + "/reports/temp/" + "DX_" + i + ".png"
        outfig_dy = args.datdir + "/reports/temp/" + "DY_" + i + ".png"
        dx.append(outfig_dx)
        dy.append(outfig_dy)
        liste = glob.glob(procdir + 'proc*fits')
        if (len(liste) != 0):
            liste.sort()
            # Pass args.datdir to align function
            liste2 = align(liste, args.datdir)
            stackim, values = stack_and_vals(args, liste2, outstack, outstack_save, outstack_save_dir, outfig_dx,
                                             outfig_dy, target_x_pos, target_y_pos, target_count)
        else:
            stackim = 'fake_path'
            values = ['unknown', 'unknown', 'unknown', 'unknown', 'unknown']
        stacks.append(stackim)
        vals.append(values)
        target_count += 1

    return stacks, vals, dx, dy

def import_images_and_vals(args):
    wea = import_weather(args)
    ramove, decmove, airmass, fwhm, skybkg, altitude, jdstart, jdend, nimages, filt, ellip, target_x_pos, target_y_pos,targs,gaia = import_night_img(
        args)
    mb, md, mfs = import_masterframes_img(args, filt)
    HAoff, DECoff = import_donuts(args, jdstart, jdend,gaia)
    water_vapor = import_water_vapor(args, jdstart, jdend,gaia)
    dlc, mlc, complc, rmsflux, newap = import_lcs_new(args,gaia)
    stacks, vals, dx, dy = import_stack_and_vals(args, target_x_pos, target_y_pos,gaia,targs)

    imgs = (
        wea, mb, md, mfs, dlc, mlc, complc, rmsflux, ramove, decmove, airmass, fwhm, skybkg, altitude, HAoff, DECoff,
        stacks, water_vapor, dx, dy)  # if order modified, pay attention to modify clean accordingly
    # print(imgs)
    return imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia


def import_schedule(args):
    mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    infile = args.obsdir + "/schedule/Archive_night_blocks/night_blocks_" + args.telescope + "_" + mydate + ".txt"
    schedule = []
    if os.path.exists(infile):
        tab = ascii.read(infile, delimiter=" ")
        targ = tab['target']
        start = tab['start time (UTC)']
        end = tab['end time (UTC)']
        duration = tab['duration (minutes)']
        rah = tab['ra (h)']
        ram = tab['ra (m)']
        ras = tab['ra (s)']
        decd = tab['dec (d)']
        decm = tab['dec (m)']
        decs = tab['dec (s)']
        conf = tab['configuration']
        schedule = [targ, start, end, duration, rah, ram, ras, decd, decm, decs, conf]
    return schedule


def check_logs(args):
    mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    # Check for weather alert in ACP log file
    acplog = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_ACP.txt"
    weatherflag = 'unknown'
    if os.path.exists(acplog):
        weatherflag = 'N'
        with open(acplog, 'r') as f:
            for line in f:
                if 'WEATHER ALERT' in line:
                    weatherflag = 'Y'
                    print('Weather Alert during the night!')
                    break
    # Check for an error in the events log file
    eventslog = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_events.csv"
    eventsflag = 'unknown'
    if os.path.exists(eventslog):
        eventsflag = 'N'
        with open(eventslog, 'r') as f:
            for line in f:
                if 'TELESCOPE.STATUS.LIST' in line:
                    # ignore standard line
                    if 'TELESCOPE.STATUS.LIST="DRIVES|0::,SYSTEM|0::,AUXILIARY|0::,UNKNOWN|0::"' in line:
                        ()
                    # ignore Callisto error message due to weather alert
                    elif (
                            args.telescope == 'Callisto' and 'DRIVES|2:DOME[0]|2;DOME[1]|2:ERR_DeviceError|axis (0) unexpectedly changed to powered on state' in line):
                        ()
                    # ignore pressure loss error
                    elif 'Working pressure suddenly lost' in line:
                        ()
                    else:
                        eventsflag = 'Y'
                        print('Technical error during the night!')
                        break
    # Check data transfer
    logfile = args.obsdir + "/download_log.txt"
    transferflag = 'unknown'
    if os.path.exists(logfile):
        report = ascii.read(logfile, format='fixed_width')
        for i in range(len(report)):
            # if this date is in the report file
            if str(report[i]['Night']) == args.date:
                i_eso = report[i]['# ESO Archive']
                i_transferred = report[i]['# Transferred']
                i_downloaded = report[i]['# Downloaded']
                if i_eso == i_transferred == i_downloaded == 0:
                    transferflag = 'No data'
                elif i_eso == i_transferred == i_downloaded:
                    transferflag = 'OK'
                else:
                    transferflag = 'KO'
                break
    # Check emergency dome closure
    dome_closure_log = args.obsdir + "/technical_logs/log_dome_closure.txt"
    closureflag = 'N'
    if os.path.exists(dome_closure_log):
        mydate = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
        start = (mydate + " 12:00:00")
        starttime = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
        endtime = starttime + datetime.timedelta(days=1)
        jdstart = Time(starttime, scale='utc').jd
        jdend = Time(endtime, scale='utc').jd
        with open(dome_closure_log, 'r') as f:
            next(f)
            for line in f:
                part = np.array(line.split(','))
                jd = float(part[0])
                if jd >= jdstart and jd <= jdend:
                    closureflag = 'Y'
    # Check for shaking indicative of strong wind of earthquake
    shaking_log = args.obsdir + "/technical_logs/" + mydate + "/" + args.telescope + "_downtime.csv"
    shakingflag = 'unknown'
    if os.path.exists(shaking_log):
        with open(shaking_log) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader)
            for row in csv_reader:
                tracking = row[8]
        csv_file.close()
        if float(tracking) <= 1.2:
            shakingflag = 'N'
        else:
            shakingflag = 'Y'
    # Return flags
    flags = [weatherflag, eventsflag, transferflag, closureflag, shakingflag]
    return flags


def clean(imgs, args,targs):
    if os.path.exists(imgs[0][1]):
        os.remove(imgs[0][1])
    if os.path.exists(imgs[1]):
        os.remove(imgs[1])
    if os.path.exists(imgs[2]):
        os.remove(imgs[2])

    for f in range(len(imgs[3])):
        if os.path.exists(imgs[3][f]):
            os.remove(imgs[3][f])

    targ_count = 0
    # s = (args.target).split(" ")
    for i in targs:
        if os.path.exists(imgs[8][targ_count]):
            os.remove(imgs[8][targ_count])
        if os.path.exists(imgs[9][targ_count]):
            os.remove(imgs[9][targ_count])
        if os.path.exists(imgs[10][targ_count]):
            os.remove(imgs[10][targ_count])
        if os.path.exists(imgs[11][targ_count]):
            os.remove(imgs[11][targ_count])
        if os.path.exists(imgs[12][targ_count]):
            os.remove(imgs[12][targ_count])
        if os.path.exists(imgs[13][targ_count]):
            os.remove(imgs[13][targ_count])
        if os.path.exists(imgs[14][targ_count]):
            os.remove(imgs[14][targ_count])
        if os.path.exists(imgs[15][targ_count]):
            os.remove(imgs[15][targ_count])
        if os.path.exists(imgs[16][targ_count]):
            os.remove(imgs[16][targ_count])
        if os.path.exists(imgs[17][targ_count]):
            os.remove(imgs[17][targ_count])
        if os.path.exists(imgs[18][targ_count]):
            os.remove(imgs[18][targ_count])
        if os.path.exists(imgs[19][targ_count]):
            os.remove(imgs[19][targ_count])
        targ_count += 1


def make_pdf(imgs, outname, args, schedule, vals, jdstart, jdend, nimages, newap, flags, ellip, order,targs,gaia):
    print('Creating PDF')

    pdf = PdfPages(outname)
    # s = (args.target).split(" ")

    appix = ['2', '2.8', '4', '5.7', '8', '11.3', '16', '20', '24', '28', '32', '40', '48']

    # Page 1 - Schedule + Weather
    page = plt.figure(dpi=300, figsize=(11.69, 8.27))
    height_ratios = [0.2, 1, 0.2, 3.0]
    gs = gridspec.GridSpec(ncols=2, nrows=4, height_ratios=height_ratios)
    gs.update(wspace=0.01, hspace=0.01, left=0.05, bottom=0.05, right=0.95, top=0.95)

    ax = page.add_subplot(gs[0, :])
    txt2 = args.date[0:4] + '-' + args.date[4:6] + '-' + args.date[6:8]
    txt = args.telescope + ' night report: ' + txt2
    ax.text(0.5, 0.8, txt, size=16, fontweight='bold', ha="center")
    ax.axis('off')

    ax = page.add_subplot(gs[1, :])
    columns = (
        'Target(s)', 'Start time\n[UTC]', 'End time\n[UTC]', 'Duration\n[minutes]', 'RA\n[hms]', 'DEC\n[dms]', 'Config')
    colours = (
        'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue',
        'lightsteelblue')
    data = []
    if (len(schedule) == 0):
        data = [['unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown', 'unknown']]
    else:
        for i in range(len(schedule[0])):
            data.append([schedule[0][i], schedule[1][i][0:10] + '\n' + schedule[1][i][11:-4],
                         schedule[2][i][0:10] + '\n' + schedule[2][i][11:-4], \
                         "{0:.2f}".format(schedule[3][i]),
                         str(int(schedule[4][i])).zfill(2) + ' ' + str(abs(int(schedule[5][i]))).zfill(
                             2) + ' ' + "{0:.2f}".format(
                             schedule[6][i]), \
                         str(int(schedule[7][i])).zfill(2) + ' ' + str(abs(int(schedule[8][i]))).zfill(
                             2) + ' ' + "{0:.2f}".format(
                             abs(schedule[9][i])), \
                         schedule[10][i].split(',')[0][2:-1] + '\n' + schedule[10][i].split(',')[1][2:-2]])
    the_table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center', colColours=colours)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)
    if (len(schedule) == 0):
        for i in range(2):
            for j in range(7):
                the_table._cells[(i, j)].set_height(.25)
                the_table._cells[(i, j)].set_width(.15)
    else:
        for i in range(len(schedule[0]) + 1):
            for j in range(7):
                the_table._cells[(i, j)].set_height(.25)
                the_table._cells[(i, j)].set_width(.15)
    ax.axis('off')
    ax.set_title('Schedule', y=0.9, fontweight='bold', size=12)

    ax = page.add_subplot(gs[3, 0])
    if os.path.exists(imgs[0][0]):
        im = Image.open(imgs[0][0])
        ax.imshow(im, interpolation='nearest')
    else:
        ax.text(0.5, 0.5, 'No plot', size=10, ha="center")
    ax.axis('off')
    ax.set_title('Boltwood weather station', fontweight='bold', size=12)

    ax = page.add_subplot(gs[3, 1])
    if os.path.exists(imgs[0][1]):
        try:
            im = Image.open(imgs[0][1])
            ax.imshow(im, interpolation='nearest')
        except:
            ax.text(0.5, 0.5, 'No plot', size=10, ha="center")
    else:
        ax.text(0.5, 0.5, 'No plot', size=10, ha="center")
    ax.axis('off')
    ax.set_title('Paranal Astronomical Site Monitor (ASM)', fontweight='bold', size=12)

    pdf.savefig(dpi=300)
    plt.close()

    # Page 2 - Events + Health Check + Master frames
    page = plt.figure(dpi=300, figsize=(11.69, 8.27))
    height_ratios = [1, 1, 4]
    ncols = 2 + len(imgs[3])
    gs = gridspec.GridSpec(ncols=ncols, nrows=3, height_ratios=height_ratios)
    gs.update(wspace=0.02, hspace=0.02, left=0.05, bottom=0.05, right=0.95, top=0.95)

    ax = page.add_subplot(gs[0, :])
    columns = ('Weather alert', 'Error(s)', 'Emergency Closure', 'Strong shaking', 'Data transfer', 'ESO report')
    colours = (
    'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue')
    data = [[flags[0], flags[1], flags[3], flags[4], flags[2], 'OK/KO [TBD]']]
    the_table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center', colColours=colours)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)
    for i in range(2):
        for j in range(6):
            the_table._cells[(i, j)].set_height(.3)
    if (flags[0] == 'N'):
        the_table._cells[(1, 0)]._text.set_color('green')
    if (flags[0] == 'Y'):
        the_table._cells[(1, 0)]._text.set_color('red')
        the_table._cells[(1, 0)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (flags[1] == 'N'):
        the_table._cells[(1, 1)]._text.set_color('green')
    if (flags[1] == 'Y'):
        the_table._cells[(1, 1)]._text.set_color('red')
        the_table._cells[(1, 1)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (flags[3] == 'N'):
        the_table._cells[(1, 2)]._text.set_color('green')
    if (flags[3] == 'Y'):
        the_table._cells[(1, 2)]._text.set_color('red')
        the_table._cells[(1, 2)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (flags[4] == 'N'):
        the_table._cells[(1, 3)]._text.set_color('green')
    if (flags[4] == 'Y'):
        the_table._cells[(1, 3)]._text.set_color('red')
        the_table._cells[(1, 3)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (flags[2] == 'OK'):
        the_table._cells[(1, 4)]._text.set_color('green')
    if (flags[2] == 'KO'):
        the_table._cells[(1, 4)]._text.set_color('red')
        the_table._cells[(1, 4)].set_text_props(fontproperties=FontProperties(weight='bold'))
    ax.axis('off')
    ax.set_title('Special Events', y=0.9, fontweight='bold', size=12)

    ax = page.add_subplot(gs[1, :])
    columns = ('Overscan \n[ADU]', 'RON \n[el]', 'Dark current \n[el/sec/pix]', 'CCD Temp. \n[deg]', 'PA \n[deg]',
               'PSF Ellipticity')
    colours = (
    'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue', 'lightsteelblue')
    overscan = 0.
    ron = 0.
    dark = 0.
    ccdtemp = 0.
    pa = 0.
    count_overscan = 0
    count_ron = 0
    count_dark = 0
    count_ccdtemp = 0
    count_pa = 0
    for i in range(len(gaia)):
        if vals[i][0] != 'unknown':
            overscan += float(vals[i][0])
            count_overscan += 1
        if vals[i][1] != 'unknown':
            ron += float(vals[i][1])
            count_ron += 1
        if vals[i][2] != 'unknown':
            dark += float(vals[i][2])
            count_dark += 1
        if vals[i][3] != 'unknown':
            ccdtemp += float(vals[i][3])
            count_ccdtemp += 1
        if vals[i][4] != 'unknown':
            pa += float(vals[i][4])
            count_pa += 1
    if count_overscan != 0:
        overscan /= count_overscan
        overscan_print = "{0:.1f}".format(overscan)
    else:
        overscan_print = 'unknown'
    if count_ron != 0:
        ron /= count_ron
        ron_print = "{0:.2f}".format(ron)
    else:
        ron_print = 'unknown'
    if count_dark != 0:
        dark /= count_dark
        dark_print = "{0:.2f}".format(dark)
    else:
        dark_print = 'unknown'
    if count_ccdtemp != 0:
        ccdtemp /= count_ccdtemp
        ccdtemp_print = "{0:.2f}".format(ccdtemp)
    else:
        ccdtemp_print = 'unknown'
    if count_pa != 0:
        pa /= count_pa
        pa_print = "{0:.2f}".format(pa)
    else:
        pa_print = 'unknown'
    if len(ellip) != 0:
        ell = max(ellip)
        ell_print = "{0:.2f}".format(ell)
    else:
        ell_print = 'unknown'

    data = [[overscan_print, ron_print, dark_print, ccdtemp_print, pa_print, ell_print]]
    the_table = ax.table(cellText=data, colLabels=columns, loc='center', cellLoc='center', colColours=colours)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)
    for i in range(2):
        for j in range(6):
            the_table._cells[(i, j)].set_height(.3)
    if (overscan_print != 'unknown'):
        if (abs(overscan - 300.0) <= 2):
            the_table._cells[(1, 0)]._text.set_color('green')
        else:
            the_table._cells[(1, 0)]._text.set_color('red')
            the_table._cells[(1, 0)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (ron_print != 'unknown'):
        if (abs(ron - 6.5) <= 1.5):
            the_table._cells[(1, 1)]._text.set_color('green')
        else:
            the_table._cells[(1, 1)]._text.set_color('red')
            the_table._cells[(1, 1)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (dark_print != 'unknown'):
        if (abs(dark - 0.27) <= 0.1):
            the_table._cells[(1, 2)]._text.set_color('green')
        else:
            the_table._cells[(1, 2)]._text.set_color('red')
            the_table._cells[(1, 2)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (ccdtemp_print != 'unknown'):
        if (abs(ccdtemp - (-60.0)) <= 1):
            the_table._cells[(1, 3)]._text.set_color('green')
        else:
            the_table._cells[(1, 3)]._text.set_color('red')
            the_table._cells[(1, 3)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (pa_print != 'unknown'):
        if (abs(pa - 360.0) <= 2 or abs(pa - 0.0) <= 2):
            the_table._cells[(1, 4)]._text.set_color('green')
        else:
            the_table._cells[(1, 4)]._text.set_color('red')
            the_table._cells[(1, 4)].set_text_props(fontproperties=FontProperties(weight='bold'))
    if (ell_print != 'unknown'):
        if (ell <= 0.5):
            the_table._cells[(1, 5)]._text.set_color('green')
        else:
            the_table._cells[(1, 5)]._text.set_color('red')
            the_table._cells[(1, 5)].set_text_props(fontproperties=FontProperties(weight='bold'))
    ax.axis('off')
    ax.set_title('Health Check', y=0.9, fontweight='bold', size=12)

    ax = page.add_subplot(gs[2, 0])
    if os.path.exists(imgs[1]):
        im = Image.open(imgs[1])
        ax.imshow(im, interpolation='nearest')
    else:
        ax.text(0.5, 0.5, 'No Masterbias', size=10, ha="center")
    ax.axis('off')
    ax = page.add_subplot(gs[2, 1])
    if os.path.exists(imgs[2]):
        im = Image.open(imgs[2])
        ax.imshow(im, interpolation='nearest')
        ax.set_title('Master calibration frames', y=1.1, fontweight='bold', size=12)
    else:
        ax.text(0.5, 0.5, 'No Masterdark', size=10, ha="center")
        ax.set_title('Master calibration frames', y=0.95, fontweight='bold', size=12)
    ax.axis('off')

    for f in range(len(imgs[3])):
        ax = page.add_subplot(gs[2, 2 + f])
        if os.path.exists(imgs[3][f]):
            im = Image.open(imgs[3][f])
            ax.imshow(im, interpolation='nearest')
        else:
            ax.text(0.5, 0.5, 'No Masterflat', size=10, ha="center")
        ax.axis('off')

    pdf.savefig(dpi=300)
    plt.close()

    # Page 3 - Summary target observations : stacked images + differential lightcurves best aperture
    page = plt.figure(dpi=300, figsize=(11.69, 8.27))
    ncols = len(targs) + (len(targs) - 1)
    width_ratios = [1]
    for i in range(len(targs)):
        width_ratios.append(0.1)
        width_ratios.append(1)
    height_ratios = [0.1, 0.1, 0.1, 1, 1]
    gs = gridspec.GridSpec(ncols=ncols, nrows=5, height_ratios=height_ratios)
    gs.update(wspace=0.0, hspace=0.0, left=0.05, bottom=0.05, right=0.95, top=0.95)

    try:

        ax = page.add_subplot(gs[0, :])
        txt = 'Target(s) overview\n'
        ax.text(0.5, 0.8, txt, size=16, fontweight='bold', ha="center")
        ax.axis('off')

        targ_count = 0
        for i in order:
            ax = page.add_subplot(gs[1, targ_count * 2])
            txt = gaia[i] + "\n (" + targs[i] + ")"
            ax.text(0.5, 0.8, txt, size=12, fontweight='bold', ha="center")
            ax.axis('off')

            ax = page.add_subplot(gs[2, targ_count * 2])
            if (jdstart[i] != 0 and jdend[i] != 0 and nimages[i] != 0):
                tstart = Time(jdstart[i], format='jd')
                tend = Time(jdend[i], format='jd', scale='utc')
                tstart_iso = tstart.iso
                tend_iso = tend.iso
                starttime = datetime.datetime.strptime(tstart_iso, "%Y-%m-%d %H:%M:%S.%f")
                endtime = datetime.datetime.strptime(tend_iso, "%Y-%m-%d %H:%M:%S.%f")
                starttime_str = (str(starttime.hour) + ":" + str(starttime.minute).zfill(2))
                endtime_str = (str(endtime.hour) + ":" + str(endtime.minute).zfill(2))
                txt = str(nimages[i]) + ' images from ' + starttime_str + ' to ' + endtime_str + ' UT'
                ax.text(0.5, 0.8, txt, size=11, ha="center")
            else:
                txt = '? images from ? to ? UT'
                ax.text(0.5, 0.8, txt, size=11, ha="center")
            ax.axis('off')

            ax = page.add_subplot(gs[3, targ_count * 2])
            if os.path.exists(imgs[16][i]):
                im = Image.open(imgs[16][i])
                ax.imshow(im, interpolation='nearest')
            else:
                ax.text(0.5, 0.5, 'No stack image', size=10, ha="center")
            ax.axis('off')

            ax = page.add_subplot(gs[4, targ_count * 2])
            index = [i * 6, 1 + (i * 6), 2 + (i * 6), 3 + (i * 6), 4 + (i * 6),
                     5 + (i * 6)]
            if os.path.exists(imgs[4][index[int(newap[i]) - 3]]):
                im = Image.open(imgs[4][index[int(newap[i]) - 3]])  # differential LC for the best aperture
                ax.imshow(im, interpolation='nearest')
            else:
                ax.text(0.5, 0.5, 'No differential LC', size=10, ha="center")
            ax.axis('off')

            if (targ_count != 0):
                ax = page.add_subplot(gs[1:, targ_count + (targ_count - 1)])
                ax.axvline(x=0.5, color='lightsteelblue')
                ax.axis('off')

            targ_count += 1

        pdf.savefig(dpi=300)
        plt.close()

        # More details for each target
        targ_count = 0
        for i in order:

            # New page - External parameters
            page = plt.figure(dpi=300, figsize=(11.69, 8.27))
            # Airmass
            plt.subplot(2, 4, 1)
            if os.path.exists(imgs[10][i]):
                im = Image.open(imgs[10][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Airmass', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Airmass', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # FWHM
            plt.subplot(2, 4, 2)
            if os.path.exists(imgs[11][i]):
                im = Image.open(imgs[11][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     FWHM', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'FWHM', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # RA MOVE
            # plt.subplot(2, 4, 3)
            # if os.path.exists(imgs[8][targ_count]):
            #    im = Image.open(imgs[8][targ_count])
            #    plt.imshow(im, interpolation='nearest')
            #    plt.gca().set_title('     RA move', fontweight='bold', size=10)
            # else:
            #    plt.text(0.5, 0.6, 'RA move', size=10, ha="center")
            #    plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            # plt.axis('off')
            # DX
            plt.subplot(2, 4, 3)
            if os.path.exists(imgs[18][i]):
                im = Image.open(imgs[18][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     DX', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'DX', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # DONUTS HA offset
            plt.subplot(2, 4, 4)
            if os.path.exists(imgs[14][i]):
                im = Image.open(imgs[14][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     DONUTS HA offset', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'DONUTS HA offset', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # Water vapor
            plt.subplot(2, 4, 5)
            if os.path.exists(imgs[17][i]):
                im = Image.open(imgs[17][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Water vapor', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Water vapor', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # Background
            plt.subplot(2, 4, 6)
            if os.path.exists(imgs[12][i]):
                im = Image.open(imgs[12][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Background', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Background', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # DEC MOVE
            # plt.subplot(2, 4, 7)
            # if os.path.exists(imgs[9][targ_count]):
            #    im = Image.open(imgs[9][targ_count])
            #    plt.imshow(im, interpolation='nearest')
            #    plt.gca().set_title('     DEC move', fontweight='bold', size=10)
            # else:
            #    plt.text(0.5, 0.6, 'DEC move', size=10, ha="center")
            #    plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            # plt.axis('off')
            # DEC MOVE
            plt.subplot(2, 4, 7)
            if os.path.exists(imgs[19][i]):
                im = Image.open(imgs[19][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     DY', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'DY', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            # DONUTS DEC offset
            plt.subplot(2, 4, 8)
            if os.path.exists(imgs[15][i]):
                im = Image.open(imgs[15][i])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     DONUTS DEC offset', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'DONUTS DEC offset', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            plt.subplots_adjust(left=0.02, bottom=0.05, right=0.98, top=0.95, wspace=0.02, hspace=0.02)
            plt.gcf().suptitle('Target #' + str(targ_count + 1) + ' - ' + gaia[i] + ' : Detailed view - 1/4\n',
                               fontweight='bold',
                               size=14)
            pdf.savefig(dpi=300)
            plt.close()

            # New Page - lightcurves (absolute and differential) for best aperture
            plt.figure(dpi=300, figsize=(11.69, 8.27))
            if os.path.exists(imgs[5][i]):
                im = Image.open(imgs[5][i])
                plt.imshow(im, interpolation='nearest')
                plt.title(
                    'Lightcurves for aperture ' + newap[i] + ' [' + appix[int(newap[i]) - 1] + ' pix]',
                    fontweight='bold', size=12)
            else:
                plt.text(0.5, 0.65, 'Lightcurves for aperture ' + newap[i] + ' [' + appix[
                    int(newap[i]) - 1] + ' pix]', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')
            plt.gcf().suptitle('Target #' + str(targ_count + 1) + ' - ' + gaia[i] + ' : Detailed view - 2/4',
                               fontweight='bold',
                               size=14)
            pdf.savefig(dpi=300)
            plt.close()

            # New Page - differential light curves for 6 different apertures
            page = plt.figure(dpi=300, figsize=(11.69, 8.27))

            plt.subplot(2, 3, 1)
            if os.path.exists(imgs[4][i * 6]):
                im = Image.open(imgs[4][i * 6])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 3', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 3', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplot(2, 3, 2)
            if os.path.exists(imgs[4][1 + (i * 6)]):
                im = Image.open(imgs[4][1 + (i * 6)])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 4', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 4', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplot(2, 3, 3)
            if os.path.exists(imgs[4][2 + (i * 6)]):
                im = Image.open(imgs[4][2 + (i * 6)])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 5', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 5', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplot(2, 3, 4)
            if os.path.exists(imgs[4][3 + (i * 6)]):
                im = Image.open(imgs[4][3 + (i * 6)])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 6', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 6', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplot(2, 3, 5)
            if os.path.exists(imgs[4][4 + (i * 6)]):
                im = Image.open(imgs[4][4 + (i * 6)])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 7', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 7', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplot(2, 3, 6)
            if os.path.exists(imgs[4][5 + (i * 6)]):
                im = Image.open(imgs[4][5 + (i * 6)])
                plt.imshow(im, interpolation='nearest')
                plt.gca().set_title('     Aperture 8', fontweight='bold', size=10)
            else:
                plt.text(0.5, 0.6, 'Aperture 8', size=10, ha="center")
                plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            plt.axis('off')

            plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.9, wspace=0.02, hspace=0.07)
            plt.gcf().suptitle('Target #' + str(targ_count + 1) + ' - ' + gaia[i] + ' : Detailed view - 3/4',
                               fontweight='bold',
                               size=14)
            pdf.savefig(dpi=300)
            plt.close()

            # New page - rms vs flux (quality of the night) for best aperture
            # plt.figure(dpi=300, figsize=(11.69, 8.27))
            # if os.path.exists(imgs[7][i]):
            #     im = Image.open(imgs[7][i])
            #     plt.imshow(im, interpolation='nearest')
            #     plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)
            #     plt.title(
            #         'Quality of the night: binned rms [5 min] vs flux for aperture ' + newap[i] + ' [' + appix[
            #             int(newap[i]) - 1] + ' pix]', fontweight='bold', size=12)
            # else:
            #     plt.text(0.5, 0.65,
            #              'Quality of the night: binned rms [5 min] vs flux for aperture ' + newap[i] + ' [' +
            #              appix[int(newap[i]) - 1] + ' pix]', size=10, ha="center")
            #     plt.text(0.5, 0.5, 'No plot', size=10, ha="center")
            # plt.axis('off')
            # plt.gcf().suptitle('Target #' + str(targ_count + 1) + ' - ' + gaia[i] + ' : Detailed view - 4/4',
            #                    fontweight='bold',
            #                    size=14)
            # pdf.savefig(dpi=300)
            # plt.close()

            targ_count += 1
    except Exception as e:
        print(e)

    pdf.close()


def create_pdf(args):
    outname = args.datdir + "/reports/" + args.telescope + "_" + args.date + ".pdf"  # to be replaced with the next line for tests
    # outname = args.datdir + "/reports/" + args.telescope + "_" + args.date + "_test.pdf" #for tests
    # fname = args.datdir + "/target_list_ids.txt"
    if not os.path.exists(args.datdir + "/reports/temp"):
        os.mkdir(args.datdir + "/reports/temp")
    # target_ids = import_spids(fname, args)
    schedule = import_schedule(args)
    flags = check_logs(args)
    imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia = import_images_and_vals(args)
    order = get_chronological_order(jdstart)
    make_pdf(imgs, outname, args, schedule, vals, jdstart, jdend, nimages, newap, flags, ellip, order, targs, gaia)
    clean(imgs, args,targs)


if __name__ == '__main__':
    description = '''Create PDF report for the night'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--datdir', required=True)
    parser.add_argument('--obsdir', required=True)
    parser.add_argument('-d', '--date', required=True)
    parser.add_argument('-t', '--target', required=True)
    parser.add_argument('-a', '--ap', required=True)
    parser.add_argument('-tel', '--telescope', required=True)
    parser.add_argument('-v', '--version', required=True)
    args = parser.parse_args()

    create_pdf(args)