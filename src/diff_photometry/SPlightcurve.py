# coding: utf-8
from astropy.io import fits, ascii
from astropy.table import Table, Column, MaskedColumn
from contextlib import contextmanager
import numpy as np
from astropy.stats import sigma_clip as sigmaclip
import astropy.units as u
import astropy.coordinates as coord
from astropy.time import Time
import matplotlib.pyplot as plt
import water_vapour.water_vapour as wv
from water_vapour.pwvGrid import getLHATPROdata
import datetime
import csv
import os
import glob
import fitsio
from astropy.coordinates.sky_coordinate import SkyCoord
import math
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
import warnings
import matplotlib
matplotlib.use('Agg')

# Suppress specific numpy warnings that clutter output during large aperture processing
warnings.filterwarnings('ignore', message="All-NaN slice encountered")
warnings.filterwarnings('ignore', message="Mean of empty slice")
warnings.filterwarnings('ignore', message="Degrees of freedom <= 0 for slice")
warnings.filterwarnings('ignore', message="The truth value of an empty array is ambiguous.*")
warnings.filterwarnings('ignore', message=".*will ignore the 'mask' of the MaskedArray.*")
warnings.filterwarnings('ignore', message="Input data contains invalid values.*", category=UserWarning)
warnings.filterwarnings('default', message=".*", category=FutureWarning)
warnings.filterwarnings('default', message=".*", category=DeprecationWarning)
warnings.filterwarnings('ignore', message="Warning: 'partition' will ignore the 'mask' of the MaskedArray.*", category=UserWarning)

init_zp = 24.16
# gain = 1.0032


class FITSImage(object):

    def __init__(self, name, nimages, napertures, data_type=np.float64):
        self.name = name
        self.data = np.zeros((napertures, nimages), dtype=data_type)

    def set_data(self, i, data):
        self.data[:, i] = data

    def hdu(self):
        return fits.ImageHDU(self.data, name=self.name)

class Telescope:

    def __init__(self,tel_diam, tel_altitude, airmass, altitude, azimuth, ccdtemp, dec_tel, exposure, filt, focuspos,gain, humidity,
                 ra_tel,ron,dark,rcore,ra_move,dec_move,pa,fwhm,seeing,psf_a_5,psf_b_5,psf_t_5,skylevel,ambtemp,telname,targname,
                 jd,bjd,bjd_tdb,hjd,date):

        self.diameter = tel_diam
        self.airmass = airmass
        self.altitude = altitude
        self.tel_alt = tel_altitude
        self.exposure = exposure
        self.ron = ron
        self.dark = dark
        self.rcore = rcore
        self.bg = skylevel
        self.gain = gain
        self.ra_move = ra_move
        self.dec_move = dec_move
        self.fwhm = fwhm
        self.psf_a_5 = psf_a_5
        self.psf_b_5 = psf_b_5
        self.psf_t_5 = psf_t_5
        self.ambtemp = ambtemp
        self.telname = telname
        self.targname = targname
        self.azimuth = azimuth
        self.ccdtemp = ccdtemp
        self.dec = dec_tel
        self.ra = ra_tel
        self.filt = filt
        self.focuspos = focuspos
        self.humidity = humidity
        self.pa = pa
        self.seeing = seeing
        self.jd = jd
        self.bjd = bjd
        self.bjd_tdb = bjd_tdb
        self.hjd = hjd
        self.date = date
        self.calculate_ellipticity(psf_a_5,psf_b_5)

    def plot(self,x,y,xlabel,ylabel,savename):
        plt.figure(figsize=(16, 8))
        plt.plot(x, y, 'k.')
        plt.title(ylabel + " against " + xlabel)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(savename)
        plt.close()

    def calculate_ellipticity(self, psf_a,psf_b):
        ellip = []
        for a, b in zip(psf_a, psf_b):
            try:
                ellip.append(math.sqrt(1 - np.power((b / a), 2)))
            except:
                ellip.append(0.0)

        self.ellipticity = ellip

    def split_vals(val):
        from collections import Counter
        vals = list(Counter(val).keys())  # equals to list(set(words))
        num_vals = list(Counter(val).values())  # counts the elements' frequency

        return vals, num_vals

    def add_water_vapour(self, targ_jd, multilc):
        from astropy.time import Time
        t = []
        water = []
        t_orig = []
        water_orig = []
        tstart, tend = [], []

        if multilc:
            times, times2, split = split_multilc(targ_jd, targ_jd)
            for ti in times:
                tstart.append(min(ti) - 0.5)
                tend.append(max(ti) + 0.5)
        else:
            tstart = [min(targ_jd) - 0.5]
            tend = [max(targ_jd) + 0.5]

        for i, (tst, tnd) in enumerate(zip(tstart, tend)):

            tS = Time(tst, format='jd')
            tS.format = 'isot'
            tE = Time(tnd, format='jd')
            tE.format = 'isot'

            try:
                pwvdata = getLHATPROdata(tS.value, tE.value, format='jd')
            except:
                try:
                    pwvdata = getLHATPROdata(tS.value, tE.value, format='jd')
                except:
                    print("No PWV data for: " + str(tS.value) + " - " + str(tE.value))
                    pwvdata = None

            if pwvdata is not None:
                t.extend(pwvdata.index.values.tolist())
                water.extend([x[0] for x in pwvdata.values])

                # Extract original data
                if hasattr(pwvdata, 'time_original') and hasattr(pwvdata, 'pwv_original'):
                    t_orig.extend(pwvdata.time_original.tolist())
                    water_orig.extend(pwvdata.pwv_original.tolist())
                else:
                    # Fallback if original data not available
                    t_orig.extend(pwvdata.index.values.tolist())
                    water_orig.extend([x[0] for x in pwvdata.values])

        # Remove duplicates while preserving order
        # Convert to pandas for easy duplicate removal
        import pandas as pd

        if len(t_orig) > 0:
            orig_df = pd.DataFrame({'time': t_orig, 'pwv': water_orig})
            orig_df = orig_df.drop_duplicates(subset=['time']).reset_index(drop=True)
            t_orig_clean = orig_df['time'].tolist()
            water_orig_clean = orig_df['pwv'].tolist()
        else:
            t_orig_clean = []
            water_orig_clean = []

        if len(t) > 0:
            proc_df = pd.DataFrame({'time': t, 'pwv': water})
            proc_df = proc_df.drop_duplicates(subset=['time']).reset_index(drop=True)
            t_clean = proc_df['time'].tolist()
            water_clean = proc_df['pwv'].tolist()
        else:
            t_clean = []
            water_clean = []

        # Store data
        self.pwv_original = water_orig_clean
        self.pwv_t_original = t_orig_clean
        self.pwv = water_clean
        self.pwv_t = t_clean

        # Find PWV points that bracket the observation period (wider search)
        i_start = find_nearest(t_clean, targ_jd[0] - (120. / 1440.))  # 2 hours before
        i_start2 = 0  # Start from beginning of observations
        i_end = find_nearest(t_clean, targ_jd[-1] + (120. / 1440.))  # 2 hours after
        i_end2 = len(targ_jd) - 1  # End at last observation

        # Make sure we have at least 2 PWV points for interpolation
        if i_end <= i_start:
            i_start = max(0, i_start - 1)
            i_end = min(len(t_clean) - 1, i_end + 1)

        pwv_subset = [water_clean[i] for i in range(i_start, i_end + 1)]
        t_subset = [t_clean[i] for i in range(i_start, i_end + 1)]

        from scipy import interpolate
        tck = interpolate.interp1d(t_subset, pwv_subset, fill_value="extrapolate")
        targ_jd_sub = (i_start2, i_end2)
        water_interp = tck(targ_jd[i_start2:i_end2 + 1])
        self.pwv_spline = water_interp

        return targ_jd_sub

class LightCurve:
    def __init__(self, id, gaia_id,nflux, flux, jd,bjd,bjd_tdb, ra, dec,pmra,pmdec,parallax, x,y, g_rp,bp_rp,teff,gmag,ap,ap_size,filt,date,exp,gain,sat,lowflux):
        self.id = id
        self.gaia_id = gaia_id
        self.nflux = nflux
        self.flux= flux
        self.jd = jd
        self.bjd = bjd
        self.bjd_tdb = bjd_tdb
        self.ra = ra
        self.dec = dec
        self.pmra = pmra
        self.pmdec = pmdec
        self.parallax = parallax
        self.x = x
        self.y = y
        self.g_rp = g_rp
        self.bp_rp = bp_rp
        self.teff = teff
        self.gmag = gmag
        self.ap = ap
        self.ap_size = ap_size
        self.filt = filt
        self.date = date
        self.target = False
        self.pwv_correct = False
        self.bw_mask = False

        self.saturated = sat
        self.lowflux = lowflux

        self.init_weights()

        # print("GAIN = ",gain)

        self.find_izmag(exp,gain)

        # print("IZMAG = ",np.median(self.izmag))


    def init_weights(self):

        if not self.saturated and not self.lowflux:
            self.weight = 1
            self.w1 = np.nanmedian(self.flux)
            self.w2 = 1.
            self.w3 = 1.
        else:
            self.weight = 0
            self.w1 = 0
            self.w2=1.
            self.w3=1.

    def diff_photom(self,comp_lc):
        difflc = np.ma.divide(self.nflux,comp_lc.nlc)
        self.lc = difflc
        self.nlc = normalise(difflc)
        self.nlc_pwv = np.empty(len(difflc))
        self.nlc_pwv[:] = np.nan
        self.nlc_pwv_bw_mask = np.empty(len(difflc))
        self.nlc_pwv_bw_mask[:] = np.nan
        self.nlc_bw_mask = np.empty(len(difflc))
        self.nlc_bw_mask[:] = np.nan

    def set_as_target(self,targname):
        self.target = True
        self.name = targname.upper()
        self.w1=0
        self.w2=0
        self.w3=0
        self.weight = 0

    def find_izmag(self,exp,gain):
        medfluxsec = (np.nanmedian(self.flux) * gain) / exp
        izmag = init_zp - (2.5 * np.log10(medfluxsec))
        self.izmag = izmag

    def sclip(self,sigma):
        self.nlc_clip = sigmaclip(self.nlc,sigma=sigma)

    def set_weight(self,w1,w2,w3, w):
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.weight = w

    def plot_standard(self,num_comp, b, savename, lim,errs, figsize):
        self.bin_lc(b)
        plt.figure(num=None, figsize=figsize)
        plt.plot(self.jd, self.nlc, 'c.', self.bin_t, self.bin_y, 'ko')
        if errs:
            plt.errorbar(self.bin_t, self.bin_y, yerr=self.bin_e, fmt='None', ecolor='k',
                         linestyle='None', alpha=0.5,
                         capthick=2)
        plt.ylabel("Relvative Flux")
        plt.xlabel("JD")
        plt.title("Lightcurve of Target %s \n binned rms = %.5f, avg error = %.5f \n Aperture Size = %.1f pixels, Number of Comp Stars = %s" %(self.name, np.nanstd(self.bin_y),np.nanmean(self.bin_e),self.ap_size,num_comp))
        if lim:
            plt.ylim(0.98, 1.02)
        else:
            if np.nanmax(self.nlc) < 1.02 and np.nanmin(self.nlc) > 0.98:
                plt.ylim(0.98, 1.02)
        plt.savefig(savename, bbox_inches='tight')
        plt.close()

        if self.pwv_correct:
            bin_t, bin_y, bin_e = bin_data(self.jd,self.nlc_pwv,b)
            plt.figure(num=None, figsize=figsize)
            plt.plot(self.jd, self.nlc_pwv, 'c.', bin_t, bin_y, 'ko')
            plt.errorbar(bin_t, bin_y, yerr=bin_e, fmt='None', ecolor='k',
                         linestyle='None', alpha=0.5,
                         capthick=2)
            plt.ylabel("Relvative Flux")
            plt.xlabel("JD")
            plt.title(
                "PWV-Corrected Lightcurve of Target %s \n binned rms = %.5f, avg error = %.5f \n Aperture Size = %.1f pixels, Number of Comp Stars = %s" % (
                self.name, np.nanstd(bin_y), np.nanmean(bin_e), self.ap_size,num_comp))
            if lim:
                plt.ylim(0.98, 1.02)
            else:
                if np.nanmax(self.nlc_pwv) < 1.02 and np.nanmin(self.nlc_pwv) > 0.98:
                    plt.ylim(0.98, 1.02)
            # plt.savefig("_".join(savename.split("_")[:-1])+"_pwv_"+str(ap), bbox_inches='tight')
            plt.savefig(savename+"_pwv", bbox_inches='tight')
            plt.close()

    def plot_multilc(self,b,savename):
        t = self.jd
        y = self.nlc
        self.sclip(5)
        y_clip = normalise(self.nlc_clip)

        pltmlc(t,y,self.date,b,savename)
        pltmlc(t,y_clip, self.date, b, savename+"_clip")

        if self.pwv_correct and self.bw_mask:
            print("Plotting bad-weather-masked, PWV-corrected LC")
            pltmlc(t, self.nlc_pwv_bw_mask, self.date, b, savename + "_pwv_bw_mask")
            pltmlc(t, normalise(sigmaclip(self.nlc_pwv_bw_mask, 5)), self.date, b, savename + "_pwv_bw_mask_clip")

        if self.pwv_correct:
            print("Plotting PWV-corrected LC")
            pltmlc(t,self.nlc_pwv,self.date,b,savename + "_pwv")
            pltmlc(t, normalise(sigmaclip(self.nlc_pwv,5)), self.date, b, savename + "_pwv_clip")

        if self.bw_mask:
            print("Plotting bad-weather-masked LC")
            pltmlc(t, self.nlc_bw_mask, self.date, b, savename + "_bw_mask")
            pltmlc(t, normalise(sigmaclip(self.nlc_bw_mask, 5)), self.date, b, savename + "_bw_mask_clip")

    def error(self,compflux, tel,ap):
        noise, noise_bin = noise_model(tel,self.izmag,5,ap)

        medfluxsec = (np.nanmedian(compflux) * tel.gain) / tel.exposure
        comp_mag = init_zp - (2.5 * np.log10(medfluxsec))
        cnoise, cnoise_bin = noise_model(tel,comp_mag,5,ap)
        # print(noise,cnoise)
        # sigma = [d * np.sqrt((n_t / (t_lc ** 2)) + (n_c / (c_lc ** 2))) for d, n_t, t_lc, n_c, c_lc in zip(self.nlc, noise, self.nflux, cnoise, comp)]
        sigma = self.nlc * np.sqrt(noise**2 + cnoise**2)
        self.err = sigma

    def custom_error(self,errs):
        self.err = errs

    def pwv_correction(self,correct_y):
        self.nlc_pwv = correct_y
        self.pwv_correct = True

    def bin_lc(self,b):
        bin_t, bin_y, bin_e = bin_data(self.jd,self.nlc,b)
        self.bin_t = bin_t
        self.bin_e = bin_e
        self.bin_y = bin_y
        self.bin = b

    def filippazzo_teff(self,dir):
        # NEW TEFF CORRECTION
        targetlist_teffs = dir + "/tests/targetListTeffEsts.csv"
        print("Extract Teffs from targetlist: ", targetlist_teffs)
        # read in Peter's target list with updated Teffs
        gaia, teffs, teffs_ppp = [], [], []
        with open(targetlist_teffs, 'r') as tfile:
            reader = csv.reader(tfile)
            for row in reader:
                try:
                    if row[2] != "gaia":
                        gaia.append(row[2].upper())
                        teffs.append(float(row[56]))
                        teffs_ppp.append(float(row[55]))
                except Exception as e:
                    print(e)
        tfile.close()

        ind = np.where(np.array(gaia) == self.gaia_id)[0][0]
        self.teff = teffs[ind]
        print("Target's T_eff from Filippazzo: " + str(self.teff))


    def targetlist_teff(self,fname_ids):
        print("Extract Teffs from targetlist: ", fname_ids)
        # find the id of the target in the field
        # basedir = os.path.dirname(os.path.dirname(os.path.dirname(dir)))
        # print(basedir)
        # fname_ids = basedir + "/SSO_targetlist_20191104.txt"
        data = ascii.read(fname_ids, delimiter=" ")
        # sp_name = [d.upper() for d in data['Sp_ID']]
        gaia_ids = [str(x) for x in data['Gaia_ID']]
        # print(data)
        teff = data['T_eff']
        # print(teff)

        try:
            # print(self.gaia_id,type(self.gaia_id))
            i = np.where(np.array(gaia_ids)==self.gaia_id)[0][0]
            self.teff = teff[i]
            print("Target Teff = " + str(self.teff))

        except:
            print("Target not in targetlist, therefore no T_eff")

    def save_mcmc_txt(self, tel, outname):
        jd = self.jd - 2450000
        bjd = self.bjd - 2450000
        bjd_tdb = self.bjd_tdb - 2450000
        data = Table([jd, bjd,bjd_tdb, self.nlc, self.err, self.nlc_pwv, tel.ra_move, tel.dec_move, tel.fwhm, tel.psf_a_5, tel.psf_b_5, tel.bg, tel.airmass, tel.exposure],names=['TMID-2450000','BJDMID-2450000','BJDTDBMID-2450000','DIFF_FLUX','ERROR','DIFF_FLUX_PWV','RA_MOVE','DEC_MOVE','FWHM','PSF_a_5', 'PSF_b_5','SKYLEVEL', 'AIRMASS', 'EXPOSURE'])
        ascii.write(data,outname, overwrite=True)

    def save_npy(self,savename):
        print("Saving numpy files: ",savename,"_*.npy")
        np.save(savename + "_jd",self.jd)
        np.save(savename + "_bjd", self.bjd)
        np.save(savename + "_bjdtdb",self.bjd_tdb)
        np.save(savename + "_lc", self.nlc)
        np.save(savename + "_lc_pwv", self.nlc_pwv)

    def best_ap(self,fname,ap,num_comp):
        new_err = np.nanmean(self.bin_e)
        new_rms = np.nanstd(self.bin_y)

        bin_t, bin_y, bin_e = bin_data(self.jd, self.nlc_pwv, self.bin)
        new_err_pwv = np.nanmean(bin_e)
        new_rms_pwv = np.nanstd(bin_y)

        if os.path.exists(fname + "_bestap.txt"):
            # Read in the file
            with open(fname + "_bestap.txt", 'r') as file:
                filedata = file.read()

            # aperture = filedata.split()[0]
            ap_err = float(filedata.split()[1])

            if ap_err > (new_err * new_rms):
                # Replace the target string
                filedata = str(ap) + " " + str(new_err * new_rms)

                # plot target diff lc within limits
                self.plot_standard(num_comp, self.bin, fname + "_bestap", True,True,(8, 8))

            # Write the file out again
            with open(fname + "_bestap.txt", 'w') as file:
                file.write(filedata)

        else:
            filedata = str(ap) + " " + str(new_err * new_rms)
            with open(fname + "_bestap.txt", 'w') as file:
                file.write(filedata)

            # plot target diff lc within limits
            self.plot_standard(num_comp, self.bin, fname + "_bestap", True, True,(8, 8))

        if os.path.exists(fname + "_bestap_pwv.txt"):
            # Read in the file
            with open(fname + "_bestap_pwv.txt", 'r') as file:
                filedata = file.read()

            ap_err_pwv = float(filedata.split()[1])

            if ap_err_pwv > (new_err_pwv * new_rms_pwv):
                # Replace the target string
                filedata = str(ap) + " " + str(new_err_pwv * new_rms_pwv)
            # Write the file out again
            with open(fname + "_bestap_pwv.txt", 'w') as file:
                file.write(filedata)
        else:
            filedata = str(ap) + " " + str(new_err_pwv * new_rms_pwv)
            if ~np.isnan(new_err_pwv * new_rms_pwv):
                with open(fname + "_bestap_pwv.txt", 'w') as file:
                    file.write(filedata)

class ArtificialLightCurve:
    def __init__(self,lcurves):
        nflux = [l.nflux for l in lcurves]
        w = [l.weight for l in lcurves]
        self.calculate_alc(nflux,np.array(w))
        self.jd = lcurves[0].jd

    def calculate_alc(self,nflux,w):
        nmult = np.ma.multiply(nflux, w[:, None])
        lc = np.nansum(nmult, axis=0)
        self.lc = lc
        self.nlc = normalise(lc)

    def bin_lc(self,b):
        bin_t, bin_y, bin_e = bin_data(self.jd, self.nlc, b)
        self.bin_t = bin_t
        self.bin_e = bin_e
        self.bin_y = bin_y

    def mask_zeros(self):
        self.nlc = np.ma.masked_less(self.nlc, 0.00001)

def convert_to_tdb(jd,ra,dec,telname):
    try:
        ip_peg = SkyCoord(float(ra), float(dec), unit=(u.deg, u.deg), frame='icrs')
    except:
        # SkyCoord('00 42 30 +41 12 00', unit=(u.hourangle, u.deg))
        ip_peg = SkyCoord(ra,dec, unit=(u.hourangle, u.deg), frame='icrs')

    if 'ARTEMIS' in telname.upper():
        # SPECULOOS North (Tiede Obs)
        site = coord.EarthLocation.from_geodetic(height=2390, lon=-16.5097, lat=28.3)
    elif 'SAINT-EX' in telname.upper():
        site = coord.EarthLocation.from_geodetic(height=2800, lon=-115.4637, lat=31.0439)
    else:
        # paranal = coord.EarthLocation.of_site('paranal')
        # PARANAL
        site = coord.EarthLocation.from_geodetic(height=2669, lon=-70.40498688, lat=-24.62743941)
    utc_times = Time(jd, format='jd', scale='utc', location=site)
    ltt_bary = utc_times.light_travel_time(ip_peg)
    bjd_tdb = (utc_times.tdb + ltt_bary).jd

    return bjd_tdb


def pltmlc(t, y, dates, b, savename):
    ts, lcs, split = split_multilc(t, y)

    # Use the same figure size as lightcurve plots
    fig, axs = plt.subplots(1, len(ts), sharey='row', sharex=True,
                            gridspec_kw={'wspace': 0, 'hspace': 0},
                            figsize=(16, 8))  # Changed from (24, 12)

    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0.5)
    fig.tight_layout()
    errs = []

    # Plot each day next to each other with no gaps
    for i in range(len(ts)):
        tminus = [math.modf(j)[0] for j in ts[i]]
        bin_t, bin_y, bin_e = bin_data(tminus, lcs[i], b)

        if len(split) > 0:
            # Larger markers to match lightcurve plots
            axs[i].plot(tminus, lcs[i], 'c.', markersize=8)
            axs[i].plot(bin_t, bin_y, 'ko', markersize=8)
        else:
            axs.plot(tminus, lcs[i], 'c.', markersize=8)
            axs.plot(bin_t, bin_y, 'ko', markersize=8)
        errs.append(np.nanmean(bin_e))

    # Increase font sizes to match lightcurve plots
    plt.xlabel('JD', fontsize=14)
    plt.ylabel('Parameter Value', fontsize=14)
    plt.tick_params(labelsize=12)

    plt.savefig(savename, bbox_inches='tight')
    plt.close()


def plot_single_param(jd, param, param_name, savename, binning):
    """
    Plot single parameter vs time in the same style as plot_standard
    """
    bin_t, bin_y, bin_e = bin_data(jd, param, binning)

    plt.figure(num=None, figsize=(8, 8))  # Match plot_standard
    plt.plot(jd, param, 'c.', bin_t, bin_y, 'ko')  # Match plot_standard colors
    plt.ylabel(param_name)
    plt.xlabel("JD")
    plt.title(f"{param_name} vs Time")
    plt.savefig(savename, bbox_inches='tight')
    plt.close()

def find_targ_id(gaia_ids,fluxes,sp_ids,tlist):
    # find the id of the target in the field
    count = 0
    names, gaia_ids_tlist= import_targetlist(tlist)
    multitarg = False
    intarg = True
    # targetid = []
    match_gaia = []
    tid_cat, tid_tlist = [],[]

    x = list(set(gaia_ids) & set(gaia_ids_tlist))
    if len(x) == 1:
        g = np.where(np.array(gaia_ids)==x[0])[0]
        if len(g)>1:
            print("WARNING: More than one star with same Gaia ID!")
            print(fluxes[g])
            i = np.argmax(np.array(fluxes)[g])
            g = g[i]
        else:
            g = g[0]
        g_tlist = np.where(np.array(gaia_ids_tlist) == x[0])[0][0]
        match_gaia.append(x[0])
        tid_cat.append(g)
        tid_tlist.append(g_tlist)
        print("\nFound object in target list, SP ID = " + str(sp_ids[tid_cat[-1]]) + " (" + str(x[0]) + ")")

    elif len(x)>1:
        # if there's more than one object in the target list in the current field then save multiple targets
        print("WARNING: More than one SPECULOOS target in field!")
        for i in x:
            g = np.where(np.array(gaia_ids) == i)[0][0]
            g_tlist = np.where(np.array(gaia_ids_tlist) == i)[0][0]
            tid_cat.append(g)
            tid_tlist.append(g_tlist)
            match_gaia.append(i)
            print("\nFound more than one object in target list, SP ID = " + str(sp_ids[g]) + " (" + str(i) + ")")
        multitarg = True
    else:
        print("WARNING: Not in current target list!")
        intarg = False

    targ_id = [sp_ids[i] for i in tid_cat]

    return targ_id, match_gaia, multitarg, intarg

def findtarg_gaia(gaia_ids_cat, fluxes_cat, sp_ids, gaia_ids_tlist):
    # find the id of the target in the field
    # count = 0
    multitarg = False
    intarg = True
    tid_cat, tid_tlist = [],[]
    match_gaia = []
    error_gaia = False

    x = list(set(gaia_ids_cat) & set(gaia_ids_tlist))
    if len(x) == 1:
        g = np.where(np.array(gaia_ids_cat)==x[0])[0]
        if len(g)>1:
            error_gaia = True
            i = np.argmax(fluxes_cat[g])
            print(fluxes_cat[g])
            g = g[i]
        else:
            g = g[0]
        g_tlist = np.where(np.array(gaia_ids_tlist) == x[0])[0][0]
        match_gaia.append(x[0])
        tid_cat.append(g)
        tid_tlist.append(g_tlist)
        print("\nFound object in target list, SP ID = " + str(sp_ids[tid_cat[-1]]) + " (" + str(x[0]) + ")")

    elif len(x)>1:
        # if there's more than one object in the target list in the current field then save multiple targets
        print("WARNING: More than one SPECULOOS target in field!")
        for i in x:
            g = np.where(np.array(gaia_ids_cat) == i)[0][0]
            g_tlist = np.where(np.array(gaia_ids_tlist) == i)[0][0]
            tid_cat.append(g)
            tid_tlist.append(g_tlist)
            match_gaia.append(i)
        multitarg = True
    else:
        print("WARNING: Not in current target list!")
        intarg = False

    if error_gaia:
        print("WARNING: More than one star with same Gaia ID!")

    return tid_cat,tid_tlist, match_gaia, multitarg, intarg, error_gaia

def normalise(a):
    # mask_zeros = np.ma.masked_where(a==0.0,a)
    # mask_zeros = np.ma.masked_less_equal(a, 0.001)
    mask_zeros = a.copy()
    mask_zeros[(a<=0.001)]=np.nan
    return a / np.nanmedian(mask_zeros)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return idx

@contextmanager
def create_output_file(filename):
    outfile = fits.HDUList()
    phdu = fits.PrimaryHDU()
    outfile.append(phdu)
    yield outfile
    outfile.writeto(filename, overwrite=True)

def outputfits(lcurves,targ, alc,bw,tel,outname,ap,version):
    # SAVE all differential lightcurves to an outputfits file
    # MAIN TABLE: All differential LCs: NUMBER OF STARS BY TIME - N x T THIS SHOULD INCLUDE TARGET: 1 x T
    # PWV CORRECTED TABLE: N x T
    # ALC: 1 x T
    # SYSTEMATICS + BAD WEATHER FLAG + WATER VAPOUR: 1 x T (things that are the same for all stars, but change with time)
    # STAR PARAMS: 1 x N (things that are the same for all times but different per star),
    # WEIGHTS: NUMBER OF WEIGHTS, W, by STAR - W x N

    if os.path.exists(outname):
        print("WARNING: Diff.fits already exists!")
    else:
        print("Creating new differential photometry fits file...")

    # sort lightcurves by weight and add target to the 'top'
    ws = extract(lcurves,'weight')
    lcurves,ws = sort_by(lcurves,ws)
    lcurves = list(lcurves[::-1])
    lcurves.insert(0, targ)

    ntimes = len(lcurves[0].jd)
    nstars = len(lcurves)

    # CREATE LCURVE, ERROR AND PWV-CORRECTED LCURVE FITS TABLES
    image = lambda name: FITSImage(name, nstars, ntimes)
    image_names = []
    image_names.extend(['LIGHTCURVE_'+str(ap)])
    image_names.extend(['ERROR_' + str(ap)])
    image_names.extend(['PWV_LIGHTCURVE_' + str(ap)])
    image_map = {name: image(name) for name in image_names}

    try:
        pwv_table = np.recarray(len(tel.pwv_t),
                                     dtype=[('JD', np.float64),
                                            ('PWV', np.float64)])
        pwv_table['JD'] = tel.pwv_t
        pwv_table['PWV'] = tel.pwv
        pwv = True
    except:
        pwv = False

    alc_table = np.recarray(ntimes,
                                 dtype=[('ALC', np.float64)])
    alc_table['ALC'] = alc.nlc

    flag_table = np.recarray(ntimes,
                                 dtype=[('FLAG', np.int64),
                                        ('BW_FLAG', np.int64)])

    flag_table['FLAG'] = [int("0b" + str(int(b)),2) for b in tel.bw_mask]#int("0b" + str(int(tel.bw_mask)),2)
    flag_table['BW_FLAG'] = tel.bw_mask

    # FLAGS
    # IMAGE SPECIFIC
    # 1 - BAD WEATHER
    # 2 - DEFOCUSSED?
    # 3 - NO WATER VAPOUR CORRECTION
    # 4 - HIGH SKY LEVEL? (is this contained in BW?)
    # 5 - ISSUE WITH ASTROMETRY? WCS RMS

    # NIGHT SPECIFIC?
    #  - DONUTS WASN'T WORKING?
    #  - DEFOCUSSED?

    # FIELD / STAR? SPECIFIC
    #  - CROWDED FIELD / BLENDED STAR
    #  - ZP
    #  - SATURATED

    for i in range(nstars):
        nlc = np.ma.filled(lcurves[i].nlc,np.nan)
        nlc_pwv = np.ma.filled(lcurves[i].nlc_pwv, np.nan)
        lc_key = 'LIGHTCURVE_'+str(ap)
        error_key = 'ERROR_' + str(ap)
        pwv_lc_key = 'PWV_LIGHTCURVE_' + str(ap)
        image_map[lc_key].set_data(i, nlc)
        image_map[error_key].set_data(i, lcurves[i].err)
        image_map[pwv_lc_key].set_data(i, nlc_pwv)

    print("initialising catalogue of objects...")
    catalogue_data = np.recarray(nstars,
                                 dtype=[('OBJ_ID', '26a'),
                                        ('RA', np.float64),
                                        ('DEC', np.float64),
                                        ('PMRA', np.float64),
                                        ('PMDEC', np.float64),
                                        ('GAIA_DR2_ID', 'a19'),
                                        ('G_RP', np.float64),
                                        ('BP_RP', np.float64),
                                        ('PARALLAX', np.float64),
                                        ('TEFF', np.float64),
                                        ('GMAG', np.float64),
                                        ('IZMAG', np.float64),
                                        ('TARGET', np.int64)])


    # all_lcurves = lcurves.append(targ)
    # all_lcurves = sort_by(all_lcurves,extract(lcurves,'id'))
    catalogue_data['OBJ_ID'] = extract(lcurves,'id')
    catalogue_data['RA'] = extract(lcurves,'ra')
    catalogue_data['DEC'] = extract(lcurves,'dec')
    catalogue_data['PMRA'] = extract(lcurves,'pmra')
    catalogue_data['PMDEC'] = extract(lcurves,'pmdec')
    catalogue_data['GAIA_DR2_ID'] = extract(lcurves,'gaia_id')
    catalogue_data['G_RP'] = extract(lcurves,'g_rp')
    catalogue_data['BP_RP'] = extract(lcurves,'bp_rp')
    catalogue_data['GMAG'] = extract(lcurves,'gmag')
    if bw:
        izmags = [np.nanmedian(l) for l in extract(lcurves, 'izmag')]
    else:
        izmags = [np.nanmedian(l[~tel.bw_mask]) for l in extract(lcurves, 'izmag')]

    catalogue_data['IZMAG'] = izmags
    catalogue_data['TEFF'] = extract(lcurves,'teff')
    catalogue_data['PARALLAX'] = extract(lcurves,'parallax')
    catalogue_data['TARGET'] = extract(lcurves, 'target')

    print("initialising weights table...")
    weight_data = np.recarray(nstars, dtype=[('SATURATED', np.int64),
                                            ('LOWFLUX', np.int64),
                                            ('W1', np.float64),
                                            ('W2', np.float64),
                                            ('W3', np.float64),
                                            ('W_TOTAL', np.float64)])

    weight_data['W1'] = extract(lcurves,'w1')
    weight_data['W2'] = extract(lcurves,'w2')
    weight_data['W3'] = extract(lcurves,'w3')
    weight_data['SATURATED'] = extract(lcurves, 'saturated')
    weight_data['LOWFLUX'] = extract(lcurves, 'lowflux')
    weight_data['W_TOTAL'] = extract(lcurves,'weight')

    print("initialising tables of image data...")
    image_data= np.recarray(ntimes,
                                 dtype=[('AIRMASS', np.float32),
                                        ('ALTITUDE', np.float32),
                                        ('AMBTEMP', np.int64),
                                        ('AZIMUTH', np.float32),
                                        ('BJD-OBS', np.float64),
                                        ('BJD-TDB', np.float64),
                                        ('CCD-TEMP', np.float32),
                                        ('DARK', np.float32),
                                        ('DEC', 'a16'),
                                        ('DEC_MOVE', np.float32),
                                        ('DIAMETER', np.int64),
                                        ('ELLIP', np.float32),
                                        ('EXPOSURE', np.float32),
                                        ('FILTER', 'a4'),
                                        ('FOCUSPOS', np.int64),
                                        ('FWHM', np.float32),
                                        ('GAIN', np.float32),
                                        ('HJD-OBS', np.float64),
                                        ('HUMIDITY', np.float32),
                                        ('JD-OBS', np.float64),
                                        ('PSF_A_5', np.float32),
                                        ('PSF_B_5', np.float32),
                                        ('PSF_T_5', np.float32),
                                        ('OBJECT', 'a20'),
                                        ('RA', 'a16'),
                                        ('RA_MOVE', np.float32),
                                        ('RCORE', np.float32),
                                        ('RON', np.float32),
                                        ('SEEING', np.float32),
                                        ('SKYLEVEL', np.float32),
                                        ('TEL_ALT', np.float32),
                                        ('TELNAME', 'a20')])

    image_data['AIRMASS'] = tel.airmass
    image_data['ALTITUDE'] = tel.altitude
    image_data['AMBTEMP'] = tel.ambtemp
    image_data['AZIMUTH'] = tel.azimuth
    image_data['BJD-OBS'] = tel.bjd
    image_data['BJD-TDB'] = tel.bjd_tdb
    image_data['CCD-TEMP'] = tel.ccdtemp
    image_data['DARK'] = tel.dark
    image_data['DEC'] = tel.dec
    image_data['DEC_MOVE'] = tel.dec_move
    image_data['DIAMETER'] = tel.diameter
    image_data['ELLIP'] = tel.ellipticity
    image_data['EXPOSURE'] = tel.exposure
    image_data['FILTER'] = tel.filt
    image_data['FOCUSPOS'] = tel.focuspos
    image_data['FWHM'] = tel.fwhm
    image_data['GAIN'] = tel.gain
    image_data['HJD-OBS'] = tel.hjd
    image_data['HUMIDITY'] = tel.humidity
    image_data['JD-OBS'] = tel.jd
    image_data['OBJECT'] = tel.targname
    image_data['PSF_A_5'] = tel.psf_a_5
    image_data['PSF_B_5'] = tel.psf_b_5
    image_data['PSF_T_5'] = tel.psf_t_5
    # image_data['PSF_A_5'] = np.full_like(tel.psf_a_5, -1.0)  # Sentinel value
    # image_data['PSF_B_5'] = np.full_like(tel.psf_b_5, -1.0)  # Sentinel value
    # image_data['PSF_T_5'] = np.full_like(tel.psf_t_5, -1.0)  # Sentinel value
    image_data['RA'] = tel.ra
    image_data['RA_MOVE'] = tel.ra_move
    image_data['RCORE'] = tel.rcore
    image_data['RON'] = tel.ron
    image_data['SEEING'] = tel.seeing
    image_data['SKYLEVEL'] = tel.bg
    image_data['TEL_ALT'] = tel.tel_alt
    image_data['TELNAME'] = tel.telname


    # CALCULATE PHOTOMETRIC ZP
    print("Calculating Photometric Zero Point...")
    phot_zp = calculate_zp(lcurves,izmags)
    print(phot_zp)
    # print(tel.bw_mask)
    # print(bw)

    print("Create fits file with all differential photometry products...")
    # create a new fits image with name given by output
    with create_output_file(outname) as hdulist:
        # add image data and catalogue data to this image
        print('add imagelist_data')
        hdulist.append(fits.BinTableHDU(image_data, name='IMAGELIST'))
        print('add catalogue_data')
        hdulist.append(fits.BinTableHDU(catalogue_data, name='CATALOGUE'))
        print('add weights')
        hdulist.append(fits.BinTableHDU(weight_data, name='WEIGHTS_'+str(ap)))

        print('add differential lcs')
        for image in list(image_map.values()):
            hdulist.append(image.hdu())

        if pwv:
            print('add pwv_table')
            hdulist.append(fits.BinTableHDU(pwv_table, name='PWV'))

        print('add alc')
        hdulist.append(fits.BinTableHDU(alc_table, name='ALC_' + str(ap)))
        print('add flags')
        hdulist.append(fits.BinTableHDU(flag_table, name='FLAGS'))
        hdulist[0].header['GAIA_ID'] = targ.gaia_id
        hdulist[0].header['PIPE_V'] = version

        if ~np.isnan(phot_zp):
            hdulist[0].header['PHOT_ZP'] = phot_zp
        else:
            hdulist[0].header['PHOT_ZP'] = "N"

        hdulist.close()

    print("Fits file created: " + outname)

def calculate_zp(lcurves,izmags):

    gmags = np.array(extract(lcurves,'gmag'))
    g_rps = np.array(extract(lcurves, 'g_rp'))
    bp_rps = np.array(extract(lcurves, 'bp_rp'))
    xs = extract(lcurves,'x')
    ys = extract(lcurves, 'y')

    # izmags = [np.nanmedian(l.izmag[bkg_mask]) for l in lcurves]
    grmags = [(l.gmag - l.g_rp) for l in lcurves]
    izmags = np.array(izmags)
    grmags = np.array(grmags)
    xs = np.array([np.nanmedian(x) for x in xs])
    ys = np.array([np.nanmedian(y) for y in ys])

    mask = ~np.isnan(izmags) & ~np.isnan(grmags) & ~np.isinf(izmags) & (izmags < 20.0) & (gmags <= 19) \
               & (g_rps < 0.8) & (bp_rps > 0.25) \
               & (xs > 20) & (ys > 20) & (xs < 1980) & (ys < 1980)

    izmags = izmags[mask]
    grmags = grmags[mask]
    zp = grmags - (izmags - init_zp)
    # zp = GR + 2.5*log10(medfluxsec)

    # izmag = phot_zp - (2.5 * np.log10(medfluxsec))
    # medfluxsec = (np.nanmedian(self.flux) * gain) / exp
    # izmag = init_zp - (2.5 * np.log10(medfluxsec))

    return np.nanmedian(zp)


def bad_weather_mask(t,cf,boxsize,thresh):
    #     thresh = 0.08
    #     boxsize = 0.5
    l = len(t)
    boxsize = boxsize/2

    std_arr = np.zeros(l)

    for i in range(0,l):
        s = t[i]
        box_s = find_nearest(t, s-boxsize)
        box_e = find_nearest(t, s+boxsize)

        # calculate standard deviation for the current window
        try:
            cf_ext = cf[box_s:box_e+1]
            std = np.nanstd(cf_ext)
            std_arr[i] = std
            # mask any very low flux points - very cloudy / obscured obsverations
            if cf[i] < 0.2:
                std_arr[i] == 999
        except Exception as e:
            # print(e)
            std_arr[i] = 999

    mask_f = np.ma.getmask(np.ma.masked_where(std_arr > thresh, std_arr))

    return mask_f


def ap_size_pixels(ap, rcore):
    ap_pix = [0.5,1./np.sqrt(2),1,np.sqrt(2),2,2*np.sqrt(2),4,5,6,7,8,10,12]
    a_pix = []
    for r in rcore:
        a_pix.append(ap_pix[ap-1]*r)
    return a_pix

def import_targetlist(fname):
    data = ascii.read(fname)
    target = [d.upper() for d in data['Name']]
    # sp_id = [str(x) for x in data['SP_id']]
    # print(data)
    gaia = [str(x) for x in data['Gaia_ID']]
    # try:
    #     teff = data['Teff']
    # except:
    #     teff = []
    return  target,gaia

def get_all_data(lc):
    ids, nfluxs,fluxs,jds,ras,decs,xs,ys,g_rps,bp_rps,teffs,gmags,izmags,ws,w1s,w2s,w3s = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    for i in range(len(lc)):
        ids.append(lc[i].id)
        nfluxs.append(lc[i].nflux)
        fluxs.append(lc[i].flux)
        jds.append(lc[i].jd)
        ras.append(lc[i].ra)
        decs.append(lc[i].dec)
        xs.append(lc[i].x)
        ys.append(lc[i].y)
        g_rps.append(lc[i].g_rp)
        bp_rps.append(lc[i].bp_rp)
        teffs.append(lc[i].teff)
        gmags.append(lc[i].gmag)
        izmags.append(lc[i].izmag)
        ws.append(lc[i].weight)
        w1s.append(lc[i].w1)
        w2s.append(lc[i].w2)
        w3s.append(lc[i].w3)

    return ids, nfluxs,fluxs,jds,ras,decs,xs,ys,g_rps,bp_rps,teffs,gmags,izmags,ws,w1s,w2s,w3s

def check_blended_stars(lcurves,ap,tel):
    ap_size = ap_size_pixels(ap, tel.rcore)[0]
    x = extract(lcurves,'x')
    y = extract(lcurves, 'y')

    for l in range(len(lcurves)):
        dist = [np.sqrt((i-x[l])**2 + (j-y[l])**2) for i,j in zip(x,y)]
        sorted_dist = sorted(dist)
        if sorted_dist[1] < ap_size:
            print("BLENDED STAR: ",lcurves[l].id)


def init_all_lcurves(tel, id, gaia_id, nflux, flux, peak, jd, bjd, bjd_tdb, ra, dec, pmra, pmdec, parallax, x, y, g_rp,
                     bp_rp, teff, gmag, ap, filt, date, exp):
    lcurves = []
    ap_size = ap_size_pixels(ap, tel.rcore)[0]
    zero_weight = 0

    for i in range(len(id)):
        saturated = False
        lowflux = False

        # test out masking 0 flux values:
        # flux[i] = np.ma.masked_where(flux[i]<=0.1,flux[i])
        flux[i] = flux[i].copy()  # Always create a writable copy
        flux[i][(flux[i] <= 0.1)] = np.nan
        # flux[i] = np.ma.masked_less_equal(flux[i],0.1)

        if i == 0:  # Only print for first 3 stars
            print(f"Star {i}: Zero flux count: {np.count_nonzero(flux[i] <= 0)}")
            print(f"Star {i}: Length/5: {len(flux[i]) / 5.}")
            print(f"Star {i}: Lengths: {len(flux[i])}, {len(exp)}")

        # Use numpy max to avoid potential read-only issues
        peak_max = np.max(peak[i])

        if peak_max > 62000:
            # print("REMOVED SATURATED STAR: " + id[i])
            saturated = True
            zero_weight = zero_weight + 1
        else:
            # Wrap potentially problematic operations in try-catch to handle read-only arrays
            try:
                zero_count = np.count_nonzero(flux[i] <= 0)

                if zero_count > (len(flux[i]) / 5.):
                    # print("REMOVED ZERO FLUX STAR: " + id[i])
                    lowflux = True
                    zero_weight = zero_weight + 1
                elif np.nanmedian(sigmaclip(flux[i], sigma=7)) < 100:
                    # print("REMOVED LOW FLUX STAR: " + id[i])
                    lowflux = True
                    zero_weight = zero_weight + 1
                elif np.count_nonzero(np.isnan(flux[i])) == len(flux[i]):
                    lowflux = True
                    zero_weight = zero_weight + 1
                else:
                    try:
                        if np.nanmedian(flux[i] / exp) < 10:
                            # print("REMOVED LOW FLUX STAR: " + id[i])
                            lowflux = True
                            zero_weight = zero_weight + 1
                    except Exception as e:
                        print("WARNING: Flux and Exposure are different lengths?!", e)
            except Exception as e:
                # Handle read-only array errors gracefully
                print(f"WARNING: Read-only array issue for star {i}: {e}")
                # Default to marking as low flux to exclude problematic stars
                lowflux = True
                zero_weight = zero_weight + 1

        # else:
        try:
            # print(id[i], gaia_id[i], nflux[i][0], flux[i][0], jd[0], ra[i], dec[i], pmra[i],pmdec[i],parallax[i], x[i][0], y[i][0], g_rp[i],bp_rp[i], teff[i],
            #                gmag[i], ap, ap_size, filt, date, tel.exposure[0], tel.gain[0])
            lcurves.append(
                LightCurve(id[i], gaia_id[i], nflux[i], flux[i], jd, bjd, bjd_tdb, ra[i], dec[i], pmra[i], pmdec[i],
                           parallax[i], x[i], y[i], g_rp[i], bp_rp[i], teff[i],
                           gmag[i], ap, ap_size, filt, date, tel.exposure, tel.gain, saturated, lowflux))
        except Exception as e:
            print(i, e)

    print("Number of possible comparison stars: " + str(len(lcurves)))
    print("Number of possible comparison stars (removing saturated and low flux stars): " + str(
        len(lcurves) - zero_weight))

    return lcurves

def import_pwv(fname):
    data = ascii.read(fname)
    try:
        t = np.array(data['JD'])
        pwv = np.array(data['PRECIPITABLE_WATER_VAPOUR [mm]'])
    except:
        t = np.array(data['col2'])
        pwv = np.array(data['col8'])
    return  pwv, t


def save_pwv_fits(timestamps, pwv_values, filename, description, gaia_id, filt, date, telname='UNKNOWN'):
    """Save PWV data to a FITS file with proper structure"""
    import numpy as np
    from astropy.io import fits

    # Create PWV table data
    pwv_table = np.recarray(len(timestamps),
                            dtype=[('JD', np.float64),
                                   ('PWV', np.float64)])
    pwv_table['JD'] = timestamps
    pwv_table['PWV'] = pwv_values

    # Create HDU list
    hdulist = fits.HDUList()

    # Primary HDU with metadata
    phdu = fits.PrimaryHDU()
    phdu.header['GAIA_ID'] = gaia_id
    phdu.header['FILTER'] = filt
    phdu.header['DATE'] = date
    phdu.header['DESCRIP'] = description
    phdu.header['NPOINTS'] = len(timestamps)
    phdu.header['TELESCOP'] = telname
    phdu.header['PIPE_V'] = 'v2'
    hdulist.append(phdu)

    # PWV data table
    pwv_hdu = fits.BinTableHDU(pwv_table, name='PWV_DATA')
    pwv_hdu.header['TTYPE1'] = 'JD'
    pwv_hdu.header['TUNIT1'] = 'days'
    pwv_hdu.header['TCOMM1'] = 'Julian Date'
    pwv_hdu.header['TTYPE2'] = 'PWV'
    pwv_hdu.header['TUNIT2'] = 'mm'
    pwv_hdu.header['TCOMM2'] = 'Precipitable Water Vapour'
    hdulist.append(pwv_hdu)

    # Write to file
    hdulist.writeto(filename, overwrite=True)
    hdulist.close()

    print(f"PWV FITS file saved: {filename}")

def pwv_correct(targ, tel, dir, outfits, pltname, lcurves, multilc):
    # NEW: Use basedir/pwv structure
    basedir = dir  # This should be the base directory passed from pipeline
    pwv_dir = os.path.join(basedir, "pwv")

    print(f"PWV correction: Looking for grids in {pwv_dir}")

    wgrid_files = glob.glob("%s/PWV_grids/*_%s_pwvGrid_coords.npy" % (pwv_dir, targ.filt))
    wdata_files = glob.glob("%s/PWV_grids/*_%s_pwvGrid_data.npy" % (pwv_dir, targ.filt))

    if not wgrid_files or not wdata_files:
        print(f"PWV grid files not found for filter {targ.filt}")
        print(f"Checked: {pwv_dir}/PWV_grids/")
        available = glob.glob("%s/PWV_grids/*.npy" % pwv_dir)
        print(f"Available files: {[os.path.basename(f) for f in available]}")
        return targ

    wgrid = wgrid_files[0]
    wdata = wdata_files[0]
    print(f"Using wgrid: {os.path.basename(wgrid)}")
    print(f"Using wdata: {os.path.basename(wdata)}")

    # Add water vapour data
    targ_jd_sub = tel.add_water_vapour(targ.jd, multilc)
    # print(targ_jd_sub)
    # update target with teffs derived from Filippazzo et al. (2015)
    # targ.filippazzo_teff(dir)

    # get all the data from the lightcurves
    # ids, nfluxs,fluxs,jds,ras,decs,xs,ys,g_rps,bp_rps,teffs,gmags,izmags,ws,w1s,w2s,w3s = get_all_data(lcurves)
    ids = extract(lcurves,'id')
    ws = extract(lcurves, 'weight')
    nfluxs = extract(lcurves, 'nflux')
    teffs = extract(lcurves, 'teff')
    # print(targ.teff)

    # create dictionary of all objects except target for the correction
    stars_dict = wv.create_stars_dict(ids, ids, ws, nfluxs, teffs, tel.airmass)

    # # CUBIC SPLINE INTERPOLATE PWV DATA TO MATCH LC TIMES:
    # from scipy import interpolate
    # tck = interpolate.splrep(t, water, s=0.8)
    # water_interp = interpolate.splev(targ_jd, tck, der=0)

    # apply PWV correction and save it to the target
    # print(pltname)
    # print(targ_jd_sub)
    # water_interp, orig_water = (t, water)
    # water_vapour.apply_water_correction(wgrid, wdata, targ_jd, targ_lc, stars_dict, nflux,
    #                                     targ_teff, airmass, water, orig_water, jd, [date],
    #                                     pltname, False)

    # wgrid, wdata, targ_jd, targ_lc, stars_dict, nflux, targ_Teff, airmass, targ_jd_sub, water, owater, jd, dates, pltname, oldformat
    # wgrid,wdata, targ_jd, targ_lc, stars_dict,targ_Teff,airmass,targ_jd_sub,water,oldformat
    corrected_y, rms, rms_uncorr = wv.apply_water_correction(wgrid, wdata,
                                                             targ.jd, targ.nlc, stars_dict, targ.teff,
                                                             tel.airmass, targ_jd_sub, tel.pwv_spline, False)

    if corrected_y is not None:
        targ.pwv_correction(corrected_y)
        print("PWV correction successful!")

        # Save PWV files after successful correction
        save_dir = os.path.dirname(outfits)
        date_str = targ.date if hasattr(targ, 'date') else 'unknown'

        # 1. Save original LHATPRO PWV data (if available)
        if hasattr(tel, 'pwv_original') and hasattr(tel, 'pwv_t_original'):
            original_file = os.path.join(save_dir, f"{targ.gaia_id}_{targ.filt}_{date_str}_pwv_original.fits")
            save_pwv_fits(tel.pwv_t_original, tel.pwv_original, original_file,
                          "Original LHATPRO PWV data before peak removal",
                          targ.gaia_id, targ.filt, date_str, tel.telname)

        # 2. Save processed PWV data (after peak removal)
        if hasattr(tel, 'pwv') and hasattr(tel, 'pwv_t'):
            processed_file = os.path.join(save_dir, f"{targ.gaia_id}_{targ.filt}_{date_str}_pwv_processed.fits")
            save_pwv_fits(tel.pwv_t, tel.pwv, processed_file,
                          "Processed LHATPRO PWV data after peak removal",
                          targ.gaia_id, targ.filt, date_str, tel.telname)

        # 3. Save interpolated PWV data (aligned with observations)
        if hasattr(tel, 'pwv_spline') and len(tel.pwv_spline) > 1:
            # Get the observation times that correspond to the interpolated PWV
            targ_jd_sub = tel.add_water_vapour(targ.jd, multilc)  # This returns the indices
            obs_times = targ.jd[targ_jd_sub[0]:targ_jd_sub[1] + 1]

            interp_file = os.path.join(save_dir, f"{targ.gaia_id}_{targ.filt}_{date_str}_pwv_interpolated.fits")
            save_pwv_fits(obs_times, tel.pwv_spline, interp_file,
                          "PWV data interpolated to observation timestamps",
                          targ.gaia_id, targ.filt, date_str, tel.telname)

    else:
        print("PWV correction unsuccessful.")

    return targ

def targ_separation(targetra, targetdec, ra, dec):
    coord = SkyCoord(ra=targetra, dec=targetdec, unit=(u.deg, u.deg), frame='icrs')
    coord_ref = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    return coord.separation(coord_ref).radian

def distance_weight(a,sep):
    w_dist = 1. / (1 + (np.ma.divide(a * np.array(sep), np.nanmax(sep))) ** 2)
    return w_dist

def init_all_weights(lcurves):
    for i in range(len(lcurves)):
        lcurves[i].init_weights()

    lcurves = set_all_weights(lcurves)
    return lcurves

def set_all_weights(lcurves):
    w1 = [x.w1 for x in lcurves]
    w2 = [x.w2 for x in lcurves]
    w3 = [x.w3 for x in lcurves]

    w1 = w1 / np.nansum(w1)
    w2 = w2 / np.nansum(w2)
    w3 = w3 / np.nansum(w3)
    w = w1 * w2 * w3
    w = w/np.nansum(w)

    for i in range(len(lcurves)):
        lcurves[i].set_weight(w1[i],w2[i],w3[i],w[i])

    return lcurves

def extract(complc,attr):
    if attr=="w1":
       return [x.w1 for x in complc]
    elif attr=="w2":
        return [x.w2 for x in complc]
    elif attr == "w3":
        return [x.w3 for x in complc]
    elif attr == "weight":
        return [x.weight for x in complc]
    elif attr == "flux":
        return [x.flux for x in complc]
    elif attr == "nflux":
        return [x.nflux for x in complc]
    elif attr == "nlc":
        return [x.nlc for x in complc]
    elif attr == "izmag":
        return [x.izmag for x in complc]
    elif attr == "teff":
        return [x.teff for x in complc]
    elif attr == "id":
        return [x.id for x in complc]
    elif attr == "gaia_id":
        return [x.gaia_id for x in complc]
    elif attr == "ra":
        return [x.ra for x in complc]
    elif attr == "dec":
        return [x.dec for x in complc]
    elif attr == "x":
        return [j.x for j in complc]
    elif attr == "y":
        return [x.y for x in complc]
    elif attr == "g_rp":
        return [x.g_rp for x in complc]
    elif attr == "bp_rp":
        return [x.bp_rp for x in complc]
    elif attr == "gmag":
        return [x.gmag for x in complc]
    elif attr == "sep":
        return [x.sep for x in complc]
    elif attr == "parallax":
        return [x.parallax for x in complc]
    elif attr == "pmra":
        return [x.pmra for x in complc]
    elif attr == "pmdec":
        return [x.pmdec for x in complc]
    elif attr == "target":
        return [x.target for x in complc]
    elif attr == "saturated":
        return [x.saturated for x in complc]
    elif attr == "lowflux":
        return [x.lowflux for x in complc]
    else:
        print(attr + " isn't an attribute of LightCurve object.")

def split_multilc(t, x):
    start = t[0]
    t = np.array(t)
    nextdays=t[np.absolute(t-start)>0.5]
    split=[]
    while nextdays!=[]:
        start = nextdays[0]
        ind_st = np.where(t==start)[0][0]
        split.append(ind_st)
        time = t[ind_st:]
        nextdays = time[np.absolute(time-start) > 0.5]

    times=np.split(t,split)
    xs=np.split(x,split)

    return times, xs, split

def plot_weights(complc,savename):
    ids = extract(complc,'id')
    complc,mags = sort_by(complc,ids)

    w1 = extract(complc,'w1')
    w2 = extract(complc,'w2')
    w3 = extract(complc,'w3')
    weight = extract(complc,'weight')

    fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=1,ncols=4,sharex=True, figsize=(16,8))
    ax1.plot(w1, 'ko')
    ax1.set_title("W1")
    ax2.plot(w2, 'ko')
    ax2.set_title("W2")
    ax3.plot(w3, 'ko')
    ax3.set_title("W3")
    ax4.plot(weight, 'ko')
    ax4.set_title("Total W")
    ax1.set_ylabel("Weight")
    plt.savefig(savename, bbox_inches='tight')
    plt.close()

def iter_algorithm(comp,targ):
    n_it = 0
    repeat = True

    while repeat==True:
        init_w = []
        n_it = n_it + 1
        if n_it ==2:
            comp = distweight(comp, targ)

        print("ITERATION " + str(n_it))
        print("FIVE HIGHEST WEIGHTED STARS: ")
        weights = extract(comp,'weight')
        ids = extract(comp, 'id')
        # sorted(list_, key=lambda x: float('-inf') if math.isnan(x) else x)
        # sort_ids = [x for (y,x) in sorted(zip(weights,ids))]
        sort_ids, sort_weights = sort_by(ids,weights)
        print(sort_ids[-5:][::-1])
        print(sort_weights[-5:][::-1])

        total_w = 0
        count_w = 0
        for w in sort_weights[::-1]:
            total_w = total_w + w
            count_w = count_w + 1
            if total_w > 0.95:
                print("Number of comparison stars (>95% of weight): " + str(count_w))
                if count_w == 1:
                    print("Algorithm has converged on one star: ",sort_ids[-1],"\nReduce weight of this star to 0 and reinitialise all weights.")
                    print(comp[np.ma.argmax(weights)].weight)
                    for i in range(len(comp)):
                        comp[i].init_weights()
                    comp[np.ma.argmax(weights)].w3 = 0
                    comp = set_all_weights(comp)
                    print(comp[np.ma.argmax(weights)].weight)
                break


        # ensure that the algorithm doesn't converge on one star with W~1
        # comp, conv = check_convergence(comp)

        print("Calculating Artificial LC")
        alc = ArtificialLightCurve(comp)
        mask = bad_weather_mask(alc.jd, alc.nlc, 0.02, 0.08)

        # if less than half the elements are masked (bad weather) then mask these elements when
        # choosing optimal flux threshold
        # print(mask,type(mask))
        if type(mask) == np.bool_:
            # print("No bad weather")
            mask = np.zeros(len(alc.jd), dtype=bool)
        else:
            # print("Some bad weather")
            # print(mask)
            # print(np.count_nonzero(mask))
            # print(np.ma.count_masked(mask),len(mask))
            # if np.count_nonzero(mask) == len(mask):
            #     print("Bad weather for entire night!")
            #     return None, None
            if np.count_nonzero(mask) < 0.5 * len(mask):
                print("Masking frames flagged as bad weather when performing iterative algorithm.")
            else:
                print("More than half the frames are flagged as bad weather, no mask used.")
                mask = np.zeros(len(alc.jd),dtype=bool)

        nlc_nozeros = np.ma.masked_less(alc.nlc,0.00001)
        # print "NUMBER OF ZEROS: " + str(np.ma.count_masked(nlc_nozeros))

        # print "Updating weight 1 with RMS values"
        for l in range(len(comp)):
            init_w.append(comp[l].weight)
            comp[l].diff_photom(alc)

            try:
                sigma = np.nanstd(np.ma.masked_where(mask, comp[l].nlc))
            except:
                sigma = np.inf

            if comp[l].w1 != 0:
                w1 = 1/ (sigma**2)
                if ~np.isinf(w1):
                    comp[l].w1 = w1
                else:
                    comp[l].w1 = 0
                # do I need to set the w1 = 0 where w3 = 0?
        comp = set_all_weights(comp)
        # check convergence?
        w = [l.weight for l in comp]
        w_diff = np.nanmean(np.ma.absolute(np.ma.subtract(init_w,w)))


        # diff < 1e-8 threshold too high, from iteration plots 1e-5 looks more reasonable
        if w_diff < 0.00001 or n_it > 30:
            print("Iteration converged (or there's been >30 iterations)")
            repeat = False

        # targ.diff_photom(alc)
        # targ.bw_mask = True
        # targ.nlc_bw_mask = np.ma.masked_where(mask, targ.nlc)
        # targ.plot_multilc(5,"/appct/data/SPECULOOSPipeline/tests/compare_new_diffphotom/TEST_GLOBAL_LIGHTCURVE2_IT" + str(n_it))

    return alc, comp, count_w

def bin_data(jd, y, b):
    mins_jd = float(b) / 1440.
    t = np.array(jd)

    # Add validation for input arrays
    if len(t) == 0 or len(y) == 0:
        print("WARNING: Empty input arrays to bin_data")
        return np.array([]), np.array([]), np.array([])

    # Handle masked arrays and convert to regular arrays with NaN
    if hasattr(y, 'mask'):
        # Convert masked array to regular array with NaN for masked values
        y_data = np.ma.filled(y, np.nan)
    else:
        y_data = np.array(y)

    # Check if arrays have matching lengths
    if len(t) != len(y_data):
        print(f"WARNING: Array length mismatch in bin_data - t: {len(t)}, y: {len(y_data)}")
        min_len = min(len(t), len(y_data))
        t = t[:min_len]
        y_data = y_data[:min_len]

    # Check if we have any valid data points
    valid_mask = ~np.isnan(y_data)
    if not np.any(valid_mask):
        print("WARNING: No valid data points in bin_data")
        return np.array([]), np.array([]), np.array([])

    split = []
    sorted_t = t
    sorted_y = y_data

    if len(sorted_t) == 0:
        return np.array([]), np.array([]), np.array([])

    start = sorted_t[0]
    nextbin = sorted_t[np.absolute(sorted_t - start) > mins_jd]

    while len(nextbin) > 0:
        start = start + mins_jd
        ind_st = np.argmax(sorted_t > start)
        if len(split) > 0:
            if ind_st != split[-1] and ind_st > 0:
                split.append(ind_st)
                time = sorted_t[ind_st:]
        else:
            if ind_st > 0:
                split.append(ind_st)
                time = sorted_t[ind_st:]

        if len(time) > 0:
            nextbin = time[np.absolute(time - start) > mins_jd]
        else:
            break

    if len(split) == 0:
        # No splits needed, treat as single bin
        times = [sorted_t]
        ys = [sorted_y]
    else:
        times = np.split(sorted_t, split)
        ys = np.split(sorted_y, split)

    bins = []
    binned_y = []
    binned_err = []

    for i in range(len(times)):
        if len(ys[i]) > 0:
            try:
                # Filter out NaN values for this bin
                valid_mask_bin = ~np.isnan(ys[i])
                if np.any(valid_mask_bin):
                    valid_times = times[i][valid_mask_bin]
                    valid_ys = ys[i][valid_mask_bin]

                    if len(valid_ys) > 0:
                        bins.append(np.nanmean(valid_times))
                        binned_y.append(np.nanmean(valid_ys))
                        n = len(valid_ys)
                        # standard error in the median:
                        binned_err.append(1.253 * np.nanstd(valid_ys) / np.sqrt(n))
                    else:
                        bins.append(np.nan)
                        binned_y.append(np.nan)
                        binned_err.append(np.nan)
                else:
                    bins.append(np.nan)
                    binned_y.append(np.nan)
                    binned_err.append(np.nan)
            except Exception as e:
                print(f"WARNING: Error in bin_data processing bin {i}: {e}")
                bins.append(np.nan)
                binned_y.append(np.nan)
                binned_err.append(np.nan)

    # Convert to arrays and filter out NaN bins
    bins = np.array(bins)
    binned_y = np.array(binned_y)
    binned_err = np.array(binned_err)

    # Filter out bins with NaN values
    valid_bins = ~(np.isnan(bins) | np.isnan(binned_y) | np.isnan(binned_err))

    bin_t = bins[valid_bins]
    bin_y = binned_y[valid_bins]
    bin_e = binned_err[valid_bins]

    return bin_t, bin_y, bin_e

def sort_by(lcurves,l):
    # sort_lcurves = [y for (x,y) in sorted(zip(l,lcurves))]
    # sort_l = sorted(l)
    l = np.array(l)
    lcurves = np.array(lcurves)

    sort_lcurves = lcurves[np.argsort(l)]
    sort_l= l[np.argsort(l)]
    return  sort_lcurves, sort_l

def rms(lcurves):
    stddev, bin_stddev = [],[]
    for l in lcurves:
        stddev.append(np.nanstd(l.nlc))
        bin_stddev.append(np.nanmean(l.bin_e))
    return stddev, bin_stddev

def noise_model(tel,mag,bin,ap):

    # exptime = np.mean(tel.exposure)
    # airmass = np.min(tel.airmass)
    # skylevel = np.min(tel.bg)
    exptime = tel.exposure
    airmass = tel.airmass
    skylevel = tel.bg

    ape = ap_size_pixels(ap, tel.rcore)[0]

    # MODELS
    gain = tel.gain # gain in e-/ADU
    # mag = np.arange(8, 18, 0.01)  # range of mag
    zp = 24.16  # zero point estimation
    # exptime = 10.  # exposure time in seconds
    overhead = 10.5  # overhead time in seconds
    # ape = 11.3 #8.0  # aperture size in pixels (8 pix for aperture 5, 11.3 pix for aperture 6)
    ron = np.mean(tel.ron) # readout noise in electrons
    scinfac = 0.09  # scintillation factor
    # airmass = 1.4  # airmass
    altitude = tel.altitude  # telescope altitude in meters
    # skylevel = 100  # background in ADU
    dark = 0.02  # dark current in e-/s/pixel
    # dark = np.mean(tel.dark)
    diameter = tel.diameter  # telescope diameter in cm

    npbin = (bin * 60) / (exptime + overhead)

    # signal
    signal = (10 ** ((mag - zp) / -2.5)) * exptime

    # poisson noise
    source_noise = np.sqrt(signal) / signal
    source_noise_bin = source_noise / np.sqrt(npbin)

    # readout noise
    npix = np.pi * (ape ** 2)
    ronape = npix * (ron ** 2)
    ron_noise = np.sqrt(ronape) / signal
    ron_noise_bin = ron_noise / np.sqrt(npbin)

    # scintillation noise according to Dravins et al. 1998
    scinti = scinfac * (diameter ** (-2. / 3.)) * airmass ** (1.75)
    scinti = scinti * np.exp(-altitude / 8000.)
    scinti = scinti / np.sqrt(2 * exptime)
    scinti_noise = np.sqrt((scinti * signal) ** 2) / signal
    scinti_noise_bin = scinti_noise / np.sqrt(npbin)

    # scintillation noise according to Osborn et al. 2015 (Paranal C = 1.56, CQ1 = 1.27, CQ3 = 1.90
    scinfac_osb = 1.27 #1.56
    scinti_osb = scinfac_osb * np.sqrt(10 * 10 ** (-6)) * (
            (diameter / 100.) ** (-2. / 3.)) * airmass ** (1.75)
    scinti_osb = scinti_osb * np.exp(-altitude / 8000.)
    scinti_osb = scinti_osb / np.sqrt(exptime)
    scinti_noise = np.sqrt((scinti_osb * signal) ** 2) / signal
    scinti_noise_bin = scinti_noise / np.sqrt(npbin)

    # background noise
    back = skylevel * gain
    backape = npix * back
    back_noise = np.sqrt(backape) / signal
    back_noise_bin = back_noise / np.sqrt(npbin)

    # dark noise
    dark = dark * exptime
    darkape = dark * npix
    dark_noise = np.sqrt(darkape) / signal
    dark_noise_bin = dark_noise / np.sqrt(npbin)

    # total noise
    total_noise = np.sqrt(
        signal + ronape + ((scinti_osb * signal) ** 2) + backape + darkape) / signal
    total_noise_bin = total_noise / np.sqrt(npbin)
    # print medfluxarr

    if np.isinf(total_noise[0]):
        print(len(exptime),len(tel.airmass),len(skylevel),len(mag),ap)

    return total_noise, total_noise_bin

def rms_izmag(lcurves,tel,targ,bin,svname,ap):

    mag = np.arange(10, 20, 0.01)  # range of mag
    noise, noise_bin = noise_model(tel,mag,bin,ap)
    rms_all,rms_bin_all = rms(lcurves)
    rms_targ, rms_bin_targ = rms([targ])
    izmag_all = extract(lcurves,'izmag')

    for i in range(len(rms_bin_all)):
        plt.semilogy(np.nanmean(izmag_all[i]),1000 * rms_bin_all[i],'k.')

    plt.semilogy(np.nanmean(targ.izmag), 1000 * rms_bin_targ[0], 'rx')
    plt.semilogy(mag, 1000 * noise_bin, linestyle='-', color='gray')
    plt.ylabel("Fractional RMS [5-minute binning](mmag)")
    plt.xlabel("$I+z'$-mag")
    plt.savefig(svname, bbox_inches="tight")
    plt.close()

def check_convergence(lcurves):

    weights = extract(lcurves, 'weight')
    lcurves, sort_w= sort_by(lcurves,weights)

    conv = False
    # print sorted_w[-1]/sorted_w[-2]
    if sort_w[-1] > 0.6:
        print("This iteration has converged on one value: " + str(lcurves[-1].id))
        print(lcurves[-1].weight)

        # CRITERIA 1: SET THE HIGHEST WEIGHTED STAR TO 0, THEN ALLOW TO CHANGE
        # lcurves[-1].w1 = 0

        # CRITERIA 2: FIX THE MAXIMUM WEIGHT AT 60%
        # lcurves[-1].w1 = 0.6

        # CRITERIA 3: SET HIGHEST WEIGHTED STAR TO 10th PLACE
        # lcurves[-1].w1 = lcurves[-10].w1

        # CRITERIA 4: SET THE HIGHEST WEIGHTED STAR TO 0
        # lcurves[-1].w3 = 0

        conv = True
        lcurves = set_all_weights(lcurves)
        print(lcurves[-1].weight)

    return lcurves, conv

def preclean(lcurves):
    numtimes = len(lcurves[0].jd)
    numstars = len(lcurves)
    weights = extract(lcurves,'weight')
    nonzero_weight = np.count_nonzero(weights)
    # clipped = 0
    print("Time points: " + str(numtimes))
    print("Number of comparison stars (with non-zero weights): " + str(nonzero_weight))

    # print("Mask 0 values")
    nflux_copy = [np.array([]) for i in range(numstars)]

    for n in range(numstars):
        nflux = lcurves[n].nflux
        x = np.ma.masked_less_equal(nflux, 0.2)
        nflux_copy[n] = np.ma.masked_less_equal(nflux, 0.2)
    nflux_copy = np.ma.array(nflux_copy)

    nflux_copy2 = sigmaclip(nflux_copy, 3, axis=0, maxiters=None)
    sclip_mask = np.ma.getmask(nflux_copy2)
    lowflux_mask = np.ma.getmask(nflux_copy)
    mask = np.logical_or(sclip_mask,lowflux_mask)

    nflux_copy = np.ma.masked_where(mask,nflux_copy)

    all_masked = np.zeros(numtimes)
    for i in range(numtimes):
        x = np.ma.count_masked(nflux_copy[:,i])
        # print(x,len(nflux_copy))
        if x == len(nflux_copy):
            all_masked[i]=1

    # For the case where there's good data before/after a major issue, this allows us to
    # extract at least the good parts of the night:
    goodtimes = all_masked[all_masked==0]
    badtimes = all_masked[all_masked == 1]
    numgoodtimes = len(goodtimes)
    numbadtimes = len(badtimes)
    print("Number of Time points with all data masked: " + str(numbadtimes))

    for j in range(numstars):
        masked = np.ma.count_masked(nflux_copy[j])
        # print(masked,numgoodtimes,numbadtimes)
        if (masked-numbadtimes) > 0.2 * numgoodtimes:#len(nflux_copy[j]):
            if lcurves[j].w3 != 0:
                # print("Time-clipped star, weight set to 0: " + lcurves[j].id)
                lcurves[j].w3 = 0
                # clipped = clipped + 1

    lcurves = set_all_weights(lcurves)
    nonzero_weight2 = np.count_nonzero(extract(lcurves,'weight'))
    print("Preclean time-clipped stars, reduced from " + str(nonzero_weight) + " to " + str(nonzero_weight2) + " stars")

    return lcurves

def distweight(lcurves,target):
    a = np.linspace(0.5, 3, 20, endpoint=True)
    seps = extract(lcurves,'sep')
    bin_stddev = []
    # weights = extract(lcurves,'weight')
    # w1 = extract(lcurves, 'w1')
    # w2 = extract(lcurves, 'w2')
    # w3 = extract(lcurves, 'w3')
    # print(np.count_nonzero(weights),np.count_nonzero(w1),np.count_nonzero(w2),np.count_nonzero(w3))

    for p in a:
        # nonzerow, zerow = 0, 0
        wdist = distance_weight(p, seps)
        for l in range(len(lcurves)):
            # if lcurves[l].w3 != 0:
            lcurves[l].w2 = wdist[l]
                # nonzerow = nonzerow + 1
            # else:
            #     zerow = zerow + 1
        # print(nonzerow, zerow)
        lcurves = set_all_weights(lcurves)
        alc = ArtificialLightCurve(lcurves)
        target.diff_photom(alc)
        # target.bin(5)
        bin_stddev.append(np.nanstd(target.nlc))

    # nonzerow, zerow = 0, 0
    opt_fthresh = np.nanargmin(bin_stddev)
    wdist = distance_weight(a[opt_fthresh], seps)
    for l in range(len(lcurves)):
        # if lcurves[l].w3 != 0:
        lcurves[l].w2 = wdist[l]
            # nonzerow = nonzerow + 1
        # else:
        #     zerow = zerow + 1
    lcurves = set_all_weights(lcurves)
    # print(nonzerow, zerow)

    print("Optimal a = " + str(a[opt_fthresh]))
    weights = extract(lcurves,'weight')
    print(np.count_nonzero(weights))

    return lcurves


def fluxthresh(lcurves, targ):
    min_mag = np.nanmedian(targ.izmag) + 3
    mags = np.nanmedian(extract(lcurves, 'izmag'), axis=1)
    sat = extract(lcurves, 'saturated')
    weights = extract(lcurves, 'weight')
    nonzero_weights = np.count_nonzero(weights)
    min_stars = 20

    print(f"DEBUG: Starting fluxthresh with {len(lcurves)} lcurves")
    print(f"DEBUG: mags shape: {mags.shape}, type: {type(mags)}")
    print(f"DEBUG: sat length: {len(sat)}")
    print(f"DEBUG: weights length: {len(weights)}")

    # Check if arrays have matching lengths
    if len(mags) != len(lcurves) or len(sat) != len(lcurves):
        print(f"ERROR: Array length mismatch - lcurves: {len(lcurves)}, mags: {len(mags)}, sat: {len(sat)}")
        return lcurves

    # Check if we have any comparison stars left
    if len(lcurves) == 0:
        print("ERROR: No comparison stars available")
        return lcurves

    # Check if we have valid magnitudes
    valid_mag_mask = ~np.isnan(mags)
    if not np.any(valid_mag_mask):
        print("ERROR: No valid magnitudes found for comparison stars")
        return lcurves

    print(f"DEBUG: {np.sum(valid_mag_mask)} stars have valid magnitudes")

    # Replace NaN magnitudes with a very large value so they get sorted to the end
    mags_for_sorting = mags.copy()
    mags_for_sorting[np.isnan(mags_for_sorting)] = 999.0

    # Check for saturated stars that also have valid magnitudes
    valid_unsaturated = valid_mag_mask & (np.array(sat) == False)
    unsaturated_count = np.sum(valid_unsaturated)

    print(f"DEBUG: {unsaturated_count} unsaturated stars with valid magnitudes")

    if unsaturated_count < min_stars:
        print(f"WARNING: Only {unsaturated_count} unsaturated stars available, need {min_stars}")
        if unsaturated_count == 0:
            print("ERROR: No unsaturated comparison stars available")
            return lcurves
        min_stars = max(1, unsaturated_count)

    # Find the magnitude threshold for minimum stars
    unsaturated_mags = mags[valid_unsaturated]
    if len(unsaturated_mags) == 0:
        print("ERROR: No unsaturated stars with valid magnitudes")
        return lcurves

    # Sort the unsaturated magnitudes to find the threshold
    sorted_unsaturated_mags = np.sort(unsaturated_mags)
    max_mag = sorted_unsaturated_mags[min_stars - 1] if len(sorted_unsaturated_mags) >= min_stars else \
    sorted_unsaturated_mags[-1]

    print(f"Maximum Magnitude: {max_mag}")
    print(f"Minimum Magnitude: {min_mag}")

    if max_mag >= min_mag:
        print("WARNING: max_mag >= min_mag, adjusting range")
        max_mag = min_mag - 1.0

    # Create threshold range
    perc = np.linspace(max_mag, min_mag, 50, endpoint=True)

    if len(perc) == 0:
        print("ERROR: No valid magnitude range for threshold")
        return lcurves

    # Now use sort_by function to sort ALL lcurves by magnitude
    print(f"DEBUG: Before sort_by - lcurves length: {len(lcurves)}, mags_for_sorting length: {len(mags_for_sorting)}")

    try:
        lcurves_sorted, mags_sorted = sort_by(lcurves, mags_for_sorting)
        print(
            f"DEBUG: After sort_by - lcurves_sorted length: {len(lcurves_sorted)}, mags_sorted length: {len(mags_sorted)}")
    except Exception as e:
        print(f"ERROR in sort_by: {e}")
        print(f"lcurves type: {type(lcurves)}, mags_for_sorting type: {type(mags_for_sorting)}")
        return lcurves

    bin_stddev, snr = [], []
    count = 0

    for p in perc:
        # Find threshold index using sorted magnitudes
        i = np.searchsorted(mags_sorted, p, side='right')

        # Reset all weights first - use original lcurves array, not sorted
        for l in range(len(lcurves)):
            lcurves[l].w3 = 0.

        # Set weights for stars passing the threshold
        # We need to map back from sorted indices to original indices
        stars_selected = 0
        for sorted_idx in range(min(i, len(lcurves_sorted))):
            # Find the original index of this sorted lightcurve
            original_idx = None
            sorted_lc = lcurves_sorted[sorted_idx]
            for orig_idx, orig_lc in enumerate(lcurves):
                if orig_lc is sorted_lc:
                    original_idx = orig_idx
                    break

            if original_idx is not None and mags_sorted[sorted_idx] < 900:  # Valid magnitude
                lcurves[original_idx].w3 = 1.
                stars_selected += 1

        if stars_selected == 0:
            bin_stddev.append(np.inf)
            snr.append(np.inf)
            count += 1
            continue

        lcurves = set_all_weights(lcurves)
        alc = ArtificialLightCurve(lcurves)

        if count == 0:
            mask = bad_weather_mask(alc.jd, alc.nlc, 0.02, 0.08)
            if type(mask) == np.bool_:
                mask = np.zeros(len(alc.jd), dtype=bool)
            else:
                if np.count_nonzero(mask) < 0.8 * len(mask):
                    print("Masking frames flagged as bad weather when determining optimal Mag threshold.")
                else:
                    print("More than 80% of the frames are flagged as bad weather, no mask.")
                    mask = np.zeros(len(alc.jd), dtype=bool)
        count = count + 1

        # Check if alc.nlc is valid
        if len(alc.nlc) == 0:
            print("WARNING: Empty artificial lightcurve")
            bin_stddev.append(np.inf)
            snr.append(np.inf)
            continue

        mask_alc = np.ma.masked_where(mask, alc.nlc)

        try:
            bin_t, bin_y, bin_e = bin_data(alc.jd, mask_alc, 5)

            # Check if binned data is valid
            if len(bin_e) == 0 or len(bin_y) == 0:
                print("WARNING: No valid binned data")
                bin_stddev.append(np.inf)
                snr.append(np.inf)
                continue

            bin_stddev.append(np.nanmean(bin_e))
            snr.append(np.nanmean(bin_e) * np.nanstd(bin_y))
        except Exception as e:
            print(f"ERROR in binning: {e}")
            bin_stddev.append(np.inf)
            snr.append(np.inf)

    # Find optimal threshold
    if len(bin_stddev) == 0 or all(np.isinf(bin_stddev)):
        print("ERROR: No valid thresholds found")
        return lcurves

    opt_fthresh = np.nanargmin(bin_stddev)
    opt_fthresh2 = np.nanargmin(snr)

    optimal_threshold = perc[opt_fthresh]
    print(f"Optimal threshold: {optimal_threshold}")

    # Apply optimal threshold using sorted arrays
    i = np.searchsorted(mags_sorted, optimal_threshold, side='right')

    print(f"Optimal index: {i}, mags_sorted length: {len(mags_sorted)}")
    if i > 0 and i <= len(mags_sorted):
        print("Optimal Magnitude Cut-off = " + str(mags_sorted[min(i - 1, len(mags_sorted) - 1)]))

    # Reset all weights and apply final threshold
    numstars = 0
    for l in range(len(lcurves)):
        lcurves[l].w3 = 0.

    # Set weights for selected stars using sorted mapping
    for sorted_idx in range(min(i, len(lcurves_sorted))):
        # Find the original index of this sorted lightcurve
        original_idx = None
        sorted_lc = lcurves_sorted[sorted_idx]
        for orig_idx, orig_lc in enumerate(lcurves):
            if orig_lc is sorted_lc:
                original_idx = orig_idx
                break

        if original_idx is not None and mags_sorted[sorted_idx] < 900:  # Valid magnitude
            lcurves[original_idx].w3 = 1.
            numstars += 1

    lcurves = set_all_weights(lcurves)
    final_numstars = np.count_nonzero(extract(lcurves, 'weight'))
    print("Reduced from " + str(nonzero_weights) + " to " + str(final_numstars) + " stars.")

    return lcurves

def sort_on_time(time,*sortarrs):
    sorted_arrs = []
    for arr in sortarrs:
        sorted_arrs.append(np.array([j for i, j in sorted(zip(time, arr))]))
    sorted_time = np.array(sorted(time))
    return sorted_time,sorted_arrs


def resid(a,model,y):
  return a[0] + (a[1]*(np.subtract(model,1)))- (y-1)

def mod(a,model):
  return 1 + (a[0] + (a[1]*(model-1)))


def background_fix_test(lc,alc, complc,savename,lcdir):
    print("Testing common field structure fix...")
    if not os.path.exists(lcdir + "/test_background_fix/"):
        os.mkdir(lcdir + "/test_background_fix/")

    w = extract(complc,'weight')
    complc, w = sort_by(complc, w)
    complc = list(complc[::-1])
    complc.insert(0, lc)

    mags = extract(complc,'izmag')
    mags = [np.nanmedian(m) for m in mags]
    comp = extract(complc,'nlc')
    w = extract(complc, 'weight')
    sat = extract(complc,'saturated')
    targ_mag = np.nanmedian(lc.izmag)

    gauss_fact = 5
    numstars = 5
    corrected=False

    new_alc = (w[:numstars] / np.sum(w[:numstars]))[:, None] * (comp[:numstars] * alc.nlc)
    # new_alc = new_alc[lc.clean_mask == False]
    # bt, bf, be = bin_data(lc.jd,(np.sum(new_alc,axis=0)/alc.nlc),b=5,e=np.zeros(len(lc.jd)))

    structure = (np.sum(new_alc, axis=0) / alc.nlc)#[lc.clean_mask == False]
    smooth_struct = gaussian_filter1d(structure, gauss_fact)
    # smooth_struct_bin = gaussian_filter1d(bf, 3)

    plt.plot(lc.jd, (np.sum(new_alc,axis=0)/alc.nlc), '.',color="green")
    plt.plot(lc.jd, smooth_struct, color='k')
    plt.ylim(0.98, 1.02)
    plt.xlabel("JD")
    plt.ylabel("Average Relative Flux")
    plt.title(str(np.nanstd(smooth_struct)))
    plt.tight_layout()
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_nightstructure", bbox_inches='tight')
    plt.close()

    a0 = [0.0, 1.0]
    a, rsqs = [], []
    model = smooth_struct

    res_lsq = least_squares(resid, a0,
                            args=(np.array(model)[~np.isnan(np.array(lc.nlc))], np.array(lc.nlc)[~np.isnan(np.array(lc.nlc))]))
    t = mod(res_lsq.x, model)
    # y = lc.clean_f / t
    sres = np.nansum((lc.nlc - mod(res_lsq.x, model)) ** 2)
    stot = np.nansum((lc.nlc - 1) ** 2)
    # rsq = 1 - (sres / stot)
    targ_a = res_lsq.x

    for c in comp:
        c2 = c#[lc.clean_mask == False]
        res_lsq = least_squares(resid, a0, args=(np.array(model)[~np.isnan(np.array(c2))], np.array(c2)[~np.isnan(np.array(c2))]))
        t = mod(res_lsq.x, model)
        y = c2 / t
        sres = np.nansum((c2 - mod(res_lsq.x, model)) ** 2)
        stot = np.nansum((c2 - 1) ** 2)
        rsq = 1 - (sres / stot)
        rsqs.append(rsq)
        a.append(res_lsq.x)

    a0 = [i[1] for i in a]

    new_mags = np.array(mags)[np.array(sat) != True]
    new_a0 = np.array(a0)[np.array(sat) != True]

    new_mags_w_targ = np.append(new_mags, targ_mag)
    new_a0_w_targ = np.append(new_a0, targ_a[1])

    stepsize = 0.5
    bins = np.arange(12,18,step=stepsize)
    bin_m, bin_a0= [],[]
    print("Create I+z-Magnitude Bins")
    for b in range(len(bins)):
        bin_m_temp, bin_a0_temp=[],[]
        for m in range(len(new_mags_w_targ)):
            if (new_mags_w_targ[m] < 15) and (new_a0_w_targ[m]<-5):
                continue
            if new_mags_w_targ[m]>bins[b]-(0.5*stepsize) and new_mags_w_targ[m]<bins[b]+(0.5*stepsize) and new_a0_w_targ[m]<10 and new_a0_w_targ[m]>-30:
                bin_a0_temp.append(new_a0_w_targ[m])
        if len(bin_a0_temp)>0:
            bin_m.append(bins[b])
            bin_a0.append(np.nanmedian(bin_a0_temp))

    # TRY LINEAR FIT
    print("Perform linear fit for scale factor against I+z-mag")
    print(targ_mag)
    min_mag=targ_mag-2
    max_mag = min(targ_mag+1,17)

    new_a0 = np.array(new_a0)[np.array(new_mags) > min_mag]
    new_mags = np.array(new_mags)[np.array(new_mags) > min_mag]

    new_a0 = np.array(new_a0)[np.array(new_mags) < max_mag]
    new_mags = np.array(new_mags)[np.array(new_mags) < max_mag]

    new_mags = np.array(new_mags)[np.logical_and(np.array(new_a0) < 5,np.array(new_a0) > -5)]
    # new_mags = np.array(new_mags)[np.array(new_a0) > -10]
    new_a0 = np.array(new_a0)[np.logical_and(np.array(new_a0) < 5,np.array(new_a0) > -5)]
    # new_a0 = np.array(new_a0)[np.array(new_a0) < 10]

    lin = np.polyfit(np.append(new_mags,targ_mag), np.append(new_a0,targ_a[1]), 1,w=1/np.append(new_mags,targ_mag))
    yfit = np.poly1d(lin)
    x_fit = np.linspace(11, 18, 1000)

    new_bin_a0 = np.array(bin_a0)[np.array(bin_m) > min_mag]
    new_bin_m = np.array(bin_m)[np.array(bin_m) > min_mag]

    new_bin_a0 = np.array(new_bin_a0)[np.array(new_bin_m) < max_mag]
    new_bin_m = np.array(new_bin_m)[np.array(new_bin_m) < max_mag]
    new_bin_m = np.array(new_bin_m)[np.logical_and(np.array(new_bin_a0) < 5, np.array(new_bin_a0) > -15)]
    new_bin_a0 = np.array(new_bin_a0)[np.logical_and(np.array(new_bin_a0) < 5, np.array(new_bin_a0) > -15)]

    try:
        lin_bin = np.polyfit(new_bin_m, new_bin_a0, 1,w=[1/b for b in new_bin_m])
        yfit_bin = np.poly1d(lin_bin)
    except Exception as e:
        print(e)
        return lc, False

    print(np.nanstd(smooth_struct))
    if np.nanstd(smooth_struct) < 0.0002:
        return lc, corrected
    print("Structure large enough to be corrected")

    if (targ_a[1] > 5) or (targ_a[1] <-5) or (lin[0]>0):
        print("Likely bad fit")
    else:
        print("Not bad fit")

    j = find_nearest(x_fit,targ_mag)
    # ymod = mod(targ_a, model)
    # new_nlc = lc.clean_f/ymod

    print(targ_a)

    print("Generate plots")
    # plt.figure(figsize=set_size(width, fraction=1,ratio=0.8))
    plt.plot(np.array(mags)[np.array(sat) != True],np.array(a0)[np.array(sat) != True],'.',color="green")
    plt.plot(x_fit,yfit_bin(x_fit),'k')
    # plt.plot(x_fit,  m + (b * x_fit),'r')
    plt.plot(new_bin_m,new_bin_a0,'k.')
    plt.plot(targ_mag,targ_a[1],color="orange",marker="^")
    plt.ylim(-30,10)
    plt.ylabel("Scale Factor")
    plt.xlabel("I+z'-Magnitude")
    # plt.show()
    plt.tight_layout()
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_m_a", bbox_inches='tight')
    plt.close()

    # if len(lc.clean_f)!=len(lc.clean_t):
    #     print('pause')
    # if len(new_nlc)!=len(lc.clean_t):
    #     print('pause')

    y4 = mod([targ_a[0], yfit_bin(x_fit[j])], model)

    print(targ_a[1],yfit(x_fit[j]),targ_a[1]-yfit(x_fit[j]))

    # lc.correct_f = new_nlc
    new_nlc = lc.nlc/y4
    corrected=True

    bt, bf, be = bin_data(lc.jd, new_nlc, 5)
    plt.figure(num=None, figsize=(8, 8))
    plt.plot(lc.jd, new_nlc, 'c.', bt, bf, 'ko')
    plt.ylabel("Relative Flux")
    plt.xlabel("JD")
    plt.ylim(0.98, 1.02)
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename), bbox_inches='tight')
    plt.close()

    fig, (ax1, ax2) = plt.subplots(nrows=2)
    bt, bf, be = bin_data(lc.jd, lc.nlc, 5)
    ax1.plot(lc.jd, lc.nlc, 'c.')
    ax1.plot(bt, bf, 'k.')
    ax1.plot(lc.jd, y4, 'g.')
    bt, bf, be = bin_data(lc.jd, new_nlc, 5)
    ax2.plot(lc.jd, new_nlc, 'c.')
    ax2.plot(bt, bf, 'k.')
    ax1.set_ylim(0.98, 1.02)
    ax2.set_ylim(0.98, 1.02)
    # plt.show()
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_beforeafter", bbox_inches='tight')
    plt.close()

    j = 1
    for i, c in zip(a[1:numstars + 1], comp[1:numstars + 1]):
        fig, (ax1, ax2) = plt.subplots(nrows=2)
        bt, bf, be = bin_data(lc.jd, c, 5)
        ax1.plot(lc.jd, c, 'c.')
        ax1.plot(bt, bf, 'k.')
        ymod = mod(i, model)
        new_nlc_comp = c / ymod
        ax1.plot(lc.jd, ymod, 'g.')
        bt, bf, be = bin_data(lc.jd, new_nlc_comp, 5)
        ax2.plot(lc.jd, new_nlc_comp, 'c.')
        ax2.plot(bt, bf, 'k.')
        ax1.set_ylim(0.98, 1.02)
        ax2.set_ylim(0.98, 1.02)
        # plt.show()
        plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_beforeafter_comp" + str(j),
                    bbox_inches='tight')
        plt.close()
        j = j + 1


def background_fix(lc,alc, complc,savename,lcdir):
    if not os.path.exists(lcdir + "/test_background_fix/"):
        os.mkdir(lcdir + "/test_background_fix/")


    w = extract(complc,'weight')
    complc, w = sort_by(complc, w)
    complc = list(complc[::-1])
    complc.insert(0, lc)

    mags = extract(complc,'izmag')
    mags = [np.nanmedian(m) for m in mags]
    comp = extract(complc,'nlc')
    w = extract(complc, 'weight')
    sat = extract(complc,'saturated')
    print(len(complc),len(w),len(sat),len(mags))

    # print(w)
    gauss_fact = 5

    numstars = 5
    new_alc = (w[1:numstars + 1] / np.sum(w[1:numstars + 1]))[:, None] * (comp[1:numstars + 1] * alc.nlc)

    structure = np.sum(new_alc, axis=0) / alc.nlc
    smooth_struct = gaussian_filter1d(structure, gauss_fact)

    # bt, bf, be = bin_data(lc.jd, np.sum(new_alc,axis=0)/alc.nlc, 5)

    plt.plot(lc.jd, np.sum(new_alc,axis=0)/alc.nlc, 'r.')
    plt.plot(lc.jd, smooth_struct, 'b')
    # plt.plot(bt, bf, 'k.')
    plt.title(str(np.nanstd(smooth_struct)) )#+ ", " + str(np.nanstd(bf)))
    plt.ylim(0.98, 1.02)
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_nightstructure", bbox_inches='tight')
    plt.close()

    a0 = [0.0, 1.0]
    a, rsqs = [], []
    model = smooth_struct

    for c in comp:
        res_lsq = least_squares(resid, a0, args=(np.array(model)[~np.isnan(np.array(c))], np.array(c)[~np.isnan(np.array(c))]))
        t = mod(res_lsq.x, model)
        y = c / t
        sres = np.nansum((c - mod(res_lsq.x, model)) ** 2)
        stot = np.nansum((c - 1) ** 2)
        rsq = 1 - (sres / stot)
        rsqs.append(rsq)
        # print(rsq)
        a.append(res_lsq.x)

    a0 = [i[1] for i in a]

    # print(a0)

    # from scipy.interpolate import make_interp_spline,make_lsq_spline
    # from scipy.signal import medfilt, wiener

    # TRY CURVE FIT
    min_mag = 13
    max_mag = 18

    new_mags = np.array(mags)[np.array(sat) != True]
    new_a0 = np.array(a0)[np.array(sat) != True]
    print(len(sat),len(mags), len(new_mags),len(new_a0))

    new_a0 = np.array(new_a0)[np.array(new_mags) > min_mag]
    new_mags = np.array(new_mags)[np.array(new_mags) > min_mag]

    new_a0 = np.array(new_a0)[np.array(new_mags) < max_mag]
    new_mags = np.array(new_mags)[np.array(new_mags) < max_mag]

    bins = np.arange(13,18,step=0.3)
    bin_m, bin_a0,bin_m_temp, bin_a0_temp = [],[],[],[]
    for b in range(len(bins[:-1])):
        for m in range(len(new_mags)):
            if new_mags[m]>bins[b] and new_mags[m]<bins[b+1]:
                bin_a0_temp.append(new_a0[m])
        if len(bin_a0_temp)>0:
            bin_m.append((bins[b] + bins[b + 1]) / 2)
            bin_a0.append(np.nanmedian(sorted(bin_a0_temp)[:5]))

    n =2
    z = np.polyfit(new_mags, new_a0, n)
    f = np.poly1d(z)
    x_new = np.linspace(np.min(new_mags), np.max(new_mags), 1000)
    y_new = f(x_new)

    z2 = np.polyfit(bin_m, bin_a0, n)
    f2 = np.poly1d(z2)
    x_new2 = np.linspace(np.min(bin_m), np.max(bin_m), 1000)
    y_new2 = f2(x_new2)

    plt.plot(mags, a0, '.')
    plt.plot(np.array(mags)[np.array(sat)==True], np.array(a0)[np.array(sat)==True], 'r.')
    plt.plot(np.nanmedian(lc.izmag), a0[0], '.',color='orange')
    plt.plot(x_new, y_new, 'b')
    plt.plot(x_new2, y_new2, 'g')
    sort_a0 = [j for i,j in sorted(zip(new_mags,new_a0))]
    sort_mags = sorted(new_mags)
    smooth_fit = gaussian_filter1d(sort_a0, 10)
    plt.plot(sort_mags, smooth_fit, 'c')

    mags_smooth = np.linspace(min(sort_mags), max(sort_mags), 100)

    # k = 3
    # t = np.linspace(min(sort_mags), max(sort_mags), 5)
    # t = np.r_[(sort_mags[0],) * k, t, (sort_mags[-1],) * k]

    # filt_a0 = wiener(sort_a0, 5)
    # plt.plot(sort_mags, filt_a0, 'r')

    # spl = make_lsq_spline(sort_mags, sort_a0, t)
    # a0_smooth = spl(mags_smooth)
    # plt.plot(mags_smooth, a0_smooth, 'orange')

    plt.ylabel("a")
    plt.xlabel("I+z mag")
    plt.xlim(11,19)
    plt.ylim(-30,5)
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_magfit", bbox_inches='tight')
    plt.close()

    if np.nanstd(smooth_struct) < 0.0002:
        return


    j = find_nearest(x_new, np.nanmedian(lc.izmag))
    # ymod = mod([a[0][0], y_new[j]], model)
    ymod = mod(a[0], model)
    new_nlc = lc.nlc/ymod

    # for c in range(len(complc)):
    #     j = find_nearest(x_new, complc[c].izmag)
    #     ymod = mod([a[c][0], y_new[j]], model)
        # complc[c].nlc = normalise(complc[c].nlc/ymod)

    bt, bf, be = bin_data(lc.jd,new_nlc,5)
    plt.figure(num=None, figsize=(8,8))
    plt.plot(lc.jd, new_nlc, 'c.', bt, bf, 'ko')
    plt.ylabel("Relvative Flux")
    plt.xlabel("JD")
    plt.title(
        "Lightcurve of Target %s \n binned rms = %.5f, avg error = %.5f \n Aperture Size = %.1f pixels, Number of Comp Stars = %s" % (
        lc.name, np.nanstd(bf), np.nanmean(be), lc.ap_size, 0))
    plt.ylim(0.98, 1.02)
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename),bbox_inches='tight')
    plt.close()

    fig, (ax1, ax2) = plt.subplots(nrows=2)
    bt, bf, be = bin_data(lc.jd, lc.nlc, 5)
    ax1.plot(lc.jd, lc.nlc, 'c.')
    ax1.plot(bt, bf, 'k.')
    ax1.plot(lc.jd, ymod, 'g.')
    bt, bf, be = bin_data(lc.jd, new_nlc, 5)
    ax2.plot(lc.jd, new_nlc, 'c.')
    ax2.plot(bt, bf, 'k.')
    ax1.set_ylim(0.98, 1.02)
    ax2.set_ylim(0.98, 1.02)
    # plt.show()
    plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_beforeafter", bbox_inches='tight')


    j=1
    for i,c in zip(a[1:numstars+1],comp[1:numstars+1]):
        fig, (ax1, ax2) = plt.subplots(nrows=2)
        bt, bf, be = bin_data(lc.jd, c, 5)
        ax1.plot(lc.jd, c, 'c.')
        ax1.plot(bt, bf, 'k.')
        ymod = mod(i, model)
        new_nlc = c / ymod
        ax1.plot(lc.jd, ymod, 'g.')
        bt, bf, be = bin_data(lc.jd, new_nlc, 5)
        ax2.plot(lc.jd, new_nlc, 'c.')
        ax2.plot(bt, bf, 'k.')
        ax1.set_ylim(0.98, 1.02)
        ax2.set_ylim(0.98, 1.02)
        # plt.show()
        plt.savefig(lcdir + "/test_background_fix/" + os.path.basename(savename) + "_beforeafter_comp"+str(j), bbox_inches='tight')
        plt.close()
        j = j+1

#         j = find_nearest(x_new,mag)
#         print(mag,x_new[j],i,y_new[j])
#         ymod = mod([i[0],y_new[j]],model)
#         fig, (ax1,ax2) = plt.subplots(nrows=2)
#         bt, bf, be = bin_data(lc.t, c, b)
#         ax1.plot(lc.t, c, 'c.')
#         ax1.plot(bt, bf, 'k.')
#         ax1.plot(lc.t, ymod, 'g.')
#         bt, bf, be = bin_data(lc.t, c/ymod, b)
#         ax2.plot(lc.t, c/ymod,'c.')
#         ax2.plot(bt, bf, 'k.')
#         ax1.set_ylim(0.98,1.02)
#         ax2.set_ylim(0.98, 1.02)
#         plt.show()

def import_outfits(outfits,goutfits,ap,targ_gaia):
    hdu = 'FLUX_' + str(int(ap))
    print("FILES:",outfits,goutfits)

    sp_id = "N"
    teff_header = "N"
    intarg = True
    no_teff = False

    # get all the dates for this target
    with fits.open(outfits) as f:
        dates = f[0].header['HISTORY']
        print(dates)
        # g_id = f[0].header['GAIA_ID']
        try:
            sp_id = f[0].header['SP_ID']
            print("FOUND SP ID IN HEADER: " + sp_id)
            teff_header = f[0].header['TEFF']
            print("FOUND TEFF IN HEADER: " + str(teff_header))
            # Fix: properly check if TEFF is valid and set no_teff accordingly
            if str(teff_header).strip() not in ["N", "", "nan"] and float(teff_header) > 0:
                no_teff = False
                print("DEBUG: no_teff set to False because valid TEFF found")
            else:
                no_teff = True
                print("DEBUG: no_teff set to True because TEFF is invalid")
        except Exception as e:
            print(f"DEBUG: Exception reading header: {e}")
            no_teff = True

    # if there isn't an SP ID in the header (the initial gaia xmatch didn't work) then flag
    if sp_id == "N":
        intarg = False
        print("TARGET NOT IN CURRENT TARGET LIST, OR PIPELINE V2 HAS NOT BEEN RUN")

    # if there isn't a TEFF in the header then flag
    if teff_header == "N":
        no_teff = True
        print("TARGET DOESN'T HAVE A GIVEN EFFECTIVE TEMPERATURE")


    # import the outfits file and extract info
    with fitsio.FITS(outfits) as infile:

        flux = infile[hdu].read()
        # hdr = infile[0].read_header()
        imagelist = infile['imagelist']
        cat = infile['catalogue']
        ra = cat['ra'].read()
        numstars = len(ra)
        dec = cat['dec'].read()
        peak = infile['peak'].read()
        exposure = imagelist['exptime'].read()
        rcore = imagelist['rcore'].read()
        obj_id = cat['obj_id'].read()
        azimuth = imagelist['azimuth'].read()
        ccdtemp = imagelist['ccd-temp'].read()
        dec_tel = imagelist['dec'].read()
        filt = imagelist['filter'].read()
        focuspos = imagelist['focuspos'].read()
        humidity = imagelist['humidity'].read()
        ra_tel = imagelist['ra'].read()
        pa = imagelist['pa'].read()
        obj_ra = imagelist['ra'].read()[0]
        obj_dec = imagelist['dec'].read()[0]
        seeing = imagelist['seeing'].read()
        ra_move = imagelist['ra_move'].read()
        dec_move = imagelist['dec_move'].read()
        fwhm = imagelist['fwhm'].read()
        psf_a_5 = imagelist['psf_a_5'].read()
        psf_b_5 = imagelist['psf_b_5'].read()
        psf_t_5 = imagelist['psf_t_5'].read()
        skylevel = imagelist['skylevel'].read()
        airmass = imagelist['airmass'].read()
        ambtemp = imagelist['ambtemp'].read()
        altitude = imagelist['altitude'].read()
        targname = imagelist['object'].read()
        telname = imagelist['telescop'].read()[0]
        print(telname)

        if "ARTEMIS" in telname.upper():
            tel_diam = 100
            # ron = 6.328
            tel_altitude = 2482
            init_gain = 1.1
        elif "SAINT-EX" in telname.upper():
            tel_diam = 100
            tel_altitude = 2800
            init_gain = 3.48
        else:
            tel_diam = 100
            # ron = 6.328
            tel_altitude = 2390
            init_gain = 1.0032

        try:
            ron = imagelist['ron'].read()
        except:
            ron = 6.328 * np.ones(len(airmass))

        try:
            darkcur = imagelist['darkcur'].read()
        except:
            darkcur = 0.3 * np.ones(len(airmass))

        try:
            gain = imagelist['gain'].read()
        except:
            gain =  init_gain * np.ones(len(airmass))

        x = infile['ccdx'].read()
        y = infile['ccdy'].read()
        jd = infile['jd'].read()[0]
        try:
            bjd = infile['bjd'].read()[0] # BJD in headers is in UTC
            hjd = infile['hjd'].read()[0]
        except:
            bjd = imagelist['bjd-obs'].read() + ((exposure/86400)/2)
            hjd = imagelist['hjd-obs'].read() + ((exposure/86400)/2)

        bjd_tdb = convert_to_tdb(jd, obj_ra, obj_dec,telname.upper()) # get BJD in TDB
        # print(bjd_tdb)

        # get the gaia crossmatch info from the global fits file
        try:
            # print("0")
            with fitsio.FITS(outfits) as gf:
                cat = gf['catalogue']
                gaia_id = cat['gaia_dr2_id'].read()
                g_rp = cat['g_rp'].read()
                bp_rp = cat['bp_rp'].read()
                gmag = cat['gmag'].read()
                teff = cat['teff'].read()
                pmra = cat['pmra'].read()
                pmdec = cat['pmdec'].read()
                parallax = cat['parallax'].read()

        except Exception as e:
            print(e)
            print("No crossmatch with Gaia DR2")
            # g_rp, bp_rp,gmag,teff,gaia_id,pmra,pmdec,parallax = [],[],[],[],[],[],[],[]

            gaia_id = np.empty(numstars)
            gaia_id[:] = np.nan
            g_rp = np.empty(numstars)
            g_rp[:] = np.nan
            bp_rp = np.empty(numstars)
            bp_rp[:] = np.nan
            gmag = np.empty(numstars)
            gmag[:] = np.nan
            teff = np.empty(numstars)
            teff[:] = np.nan
            pmra = np.empty(numstars)
            pmra[:] = np.nan
            pmdec = np.empty(numstars)
            pmdec[:] = np.nan
            parallax = np.empty(numstars)
            parallax[:] = np.nan
            # include crossmatch here

        # if the Gaia ID has been provided by the user use that instead of SP ID in the header:
        if len(targ_gaia) >= 17:
            try:
                print("Use Gaia ID provided by user: " + str(targ_gaia))
                # print(gaia_id)
                i = np.where(np.array(gaia_id) == targ_gaia)[0]
                # print(i)
                if len(i)>1:
                    print("WARNING: More than one star with same Gaia ID!")
                    print("Select the brighter object")
                    fluxes = np.nanmedian(flux,axis=1)
                    print(fluxes[i])
                    j = np.argmax(np.array(fluxes)[i])
                    i = i[j]
                else:
                    i = i[0]
                sp_id = obj_id[i]
                intarg = True
                print("Found target in fits file: " + str(sp_id))
                # teff_header = "N"
                # no_teff = True
            except:
                intarg = False

        print("TIME SORTING...")
        # jd = [range(len(jd[0]))[::-1]]
        # print type(skylevel)

        # for var in [exposure,rcore,ra_move,dec_move,fwhm,psf_a_5,psf_b_5,skylevel,ambtemp,targname,airmass]:
        sorted_jd, [exposure, rcore, ra_move, dec_move, fwhm, ron, darkcur, psf_a_5, psf_b_5, psf_t_5, skylevel, airmass, ambtemp, \
            altitude, targname, azimuth, bjd, bjd_tdb, ccdtemp, dec_tel, filt, focuspos, gain, hjd, humidity, ra_tel, pa, \
             seeing] = sort_on_time(jd,exposure,rcore,ra_move,dec_move,fwhm,ron, darkcur, psf_a_5, psf_b_5, psf_t_5,skylevel, \
                                    airmass, ambtemp, altitude, targname, azimuth, bjd, bjd_tdb,ccdtemp, dec_tel, filt, focuspos,\
                                    gain, hjd, humidity,ra_tel, pa, seeing)


        for i in range(len(flux)):
            flux[i] = np.array([j for _, j in sorted(zip(jd, flux[i]))])

        jd = np.array(sorted_jd)
        nflux = [normalise(f) for f in flux]

        if "TRAPPIST" in targname[0].upper():
            targname[0] = "SP2306-0502"


        print("Initialising telescope object...")

        tel_args = [tel_diam, tel_altitude, airmass, altitude, azimuth,ccdtemp, dec_tel,exposure,filt,focuspos,gain,humidity,ra_tel, ron, darkcur,
                    rcore, ra_move, dec_move,pa, fwhm,seeing, psf_a_5, psf_b_5,psf_t_5, skylevel,ambtemp,telname,targname,jd,bjd,bjd_tdb,hjd,dates[0]]
        # arg_names = ['tel_diam', 'airmass', 'altitude','azimuth','ccdtemp','dec_tel', 'exposure','filt','focuspos','gain',
        #              'humidity','ra_tel', 'ron', 'rcore', 'ra_move', 'dec_move','pa', 'fwhm', 'seeing', 'psf_a_5',
        #              'psf_b_5', 'skylevel', 'ambtemp','telname','targname','jd','bjd','hjd','date']

        tel = Telescope(*tel_args)
        return intarg, sp_id, no_teff, teff_header, obj_id, gaia_id, nflux, flux, peak, jd,bjd,bjd_tdb, ra, dec, pmra,pmdec,parallax, x,y, g_rp, bp_rp,teff,gmag,tel,exposure, targname


def plot_pwv_comparison(tel, targ_jd, targ_gaia, filt, date, savename):
    """
    Plot comparison of original, processed, and interpolated PWV data
    """
    plt.figure(num=None, figsize=(8, 8))  # Match plot_standard figsize

    # Plot 1: Original LHATPRO data (if available) - 2x larger points
    if hasattr(tel, 'pwv_original') and hasattr(tel, 'pwv_t_original'):
        plt.plot(tel.pwv_t_original, tel.pwv_original, 'c.', alpha=0.6,
                 label='Original LHATPRO', markersize=12)  # 2x default size (6)

    # Plot 2: Processed PWV data - only if it exists AND is actually different from original
    plot_processed = False
    if hasattr(tel, 'pwv') and hasattr(tel, 'pwv_t'):
        if hasattr(tel, 'pwv_original') and hasattr(tel, 'pwv_t_original'):
            # Check if processed data is actually different from original
            times_different = not np.array_equal(tel.pwv_t, tel.pwv_t_original)
            values_different = not np.array_equal(tel.pwv, tel.pwv_original)

            if times_different or values_different:
                plot_processed = True
                print(f"DEBUG: Plotting processed data - data differs from original")
                print(f"  - Times different: {times_different}")
                print(f"  - Values different: {values_different}")
                print(f"  - Original length: {len(tel.pwv_original)}, Processed length: {len(tel.pwv)}")
            else:
                print(f"DEBUG: Not plotting processed data - identical to original")
                print(f"  - Length: {len(tel.pwv)} (both)")
        else:
            # If no original data, plot processed anyway
            plot_processed = True
            print(f"DEBUG: Plotting processed data - no original data for comparison")

    if plot_processed:
        plt.plot(tel.pwv_t, tel.pwv, 'r.', alpha=0.8,
                 label='Processed (peaks removed)', markersize=24)  # 4x default size (6)

    # Plot 3: Interpolated PWV data (aligned with observations) - line only, no markers
    if hasattr(tel, 'pwv_spline') and len(tel.pwv_spline) > 1:
        obs_subset = targ_jd
        plt.plot(obs_subset, tel.pwv_spline, 'k-',
                 label='Interpolated to observations')  # Just 'k-' for black line, no markers

    # Add observation time range as vertical lines
    plt.axvline(targ_jd[0], color='orange', linestyle='--', alpha=0.7,
                label='Observation start/end')
    plt.axvline(targ_jd[-1], color='orange', linestyle='--', alpha=0.7)

    plt.ylabel("PWV (mm)")
    plt.xlabel("JD")
    plt.title(f"PWV Data Comparison - {targ_gaia} {filt} {date}")
    plt.legend()

    plt.savefig(savename, bbox_inches='tight')
    plt.close()

    print(f"PWV comparison plot saved: {savename}")

def main(date, targ_gaia, ap, filt, outfits, goutfits, globallc, binning, version, outdir, lcdir, targ_teff, tlist,
             oldtlists, basedir=None):
    if basedir is not None:
        # Use the basedir passed from the pipeline
        print(f"Using basedir from pipeline: {basedir}")
    else:
        # Fallback to the old method if not provided
        if os.path.basename(outdir) == "output":
            basedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(outdir))))
        else:
            # Handle other cases
            count = 0
            newdir = outdir
            while os.path.basename(newdir) != "output" and count < 10:
                newdir = os.path.dirname(newdir)
                count = count + 1

            if count < 10:
                # Found output directory, go up to base
                basedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(newdir))))
            else:
                print("Can't find base directory for PWV")
                # Fallback: assume PWV is in same directory as outdir
                basedir = outdir

    print("BASE DIRECTORY:", basedir)
    print("PWV DIRECTORY:", os.path.join(basedir, "pwv"))
    bw = False

    # IMPORT all the objects from this field
    intarg, targ_id, no_teff, teff_header, id, gaia_id, nflux, flux,peak, jd, bjd,bjd_tdb, ra, dec,pmra,pmdec,parallax, x,y, g_rp,bp_rp,teff,gmag,tel,exp, targname = import_outfits(outfits,goutfits,ap,targ_gaia)
    # INITIALISE ALL LIGHTCURVE OBJECTS
    lcurves = init_all_lcurves(tel, id, gaia_id, nflux, flux, peak, jd, bjd,bjd_tdb, ra, dec,pmra,pmdec,parallax, x, y, g_rp, bp_rp,teff, gmag, ap, filt, date,exp)

    # GET GAIA IDS from output.fits file and find which objects are in the target list
    gaia_id = extract(lcurves,'gaia_id')
    id = extract(lcurves,'id')

    if intarg == False:
        if len(oldtlists)==0:
            oldtlists = ["/appct/data/SPECULOOSPipeline/tests/target_list_ids_201905.txt"]

        for otlist in oldtlists:
            # tlist = basedir + "/" + tel + "/target_list_ids.txt"
            # tlist = basedir + "/SSO_targetlist_20191104.txt"
            # old_tlist = basedir + "/tests/target_list_ids_201905.txt"
            print("TRY FINDING TARGET IN OLDER TARGET LIST: "+ otlist)
            try:
                targ_id, match_gaia, multitarg, intarg = find_targ_id(gaia_id,np.nanmedian(flux,axis=1),id,otlist)
            except Exception as e:
                print(e)
            if intarg == True:
                targ_gaia = match_gaia[0]
                if multitarg:
                    print("MULTIPLE TARGETS IN THIS FIELD")
                break

    # if the target is in the target list
    if intarg == True:
        # if there is only one object in this field in the target list
        # if mtarg==False:
        # for t_i in targ_id:

        # FIND AND SET TARGET OBJECT
        try:
            # targetid = np.where(np.array(gaia_id)==targ_gaia)[0]
            # print(id[targetid[0]],targ_id)
            targetid = np.where(np.array(id)==targ_id)[0][0]
        except:
            print("ISSUE WITH TARGET (SATURATED OR NO FLUX)")
            # return None
        target_lcurve = lcurves[targetid]
        target_lcurve.set_as_target(targname[0])

        pltname = lcdir + "/" + targ_gaia + "_" +filt + "_"
        orig_pltname = pltname
        if globallc == False:
            pltname = pltname + date + "_" + str(target_lcurve.ap)
        else:
            pltname = pltname + str(target_lcurve.ap)
        print("Save lightcurves: ",pltname)
        # REMOVE TARGET FROM COMPARISON LCs
        print("Remove target star")
        comp_lcurves = np.delete(lcurves,targetid)

        # INITIALISE LC WEIGHTS
        print("SET WEIGHTS")
        comp_lcurves = set_all_weights(comp_lcurves)

        # MAGNITUDE CUT-OFF based on ALC
        print("FLUX THRESH")
        comp_lcurves = fluxthresh(comp_lcurves,target_lcurve)
        # INITIAL VARIABILITY CHECK - remove LCs which behave differently to other LCs
        print("PRECLEAN LCS")
        comp_lcurves = preclean(comp_lcurves)

        # DISTANCE WEIGHT
        for l in comp_lcurves:
            sep = targ_separation(target_lcurve.ra, target_lcurve.dec, l.ra, l.dec)
            l.sep = sep
        # comp_lcurves = distweight(comp_lcurves,target_lcurve)

        # PERFORM ITERATIVE ALGORITHM using the std dev of differential LCs as new weights and return final ALC
        alc, comp_lcurves, num_comp = iter_algorithm(comp_lcurves,target_lcurve)
        if alc == None:
            return None
        # FINAL TARGET DIFFERENTIAL LC (using final ALC)
        target_lcurve.diff_photom(alc)

        weights = extract(comp_lcurves, 'weight')
        flux = extract(comp_lcurves, 'flux')
        i = np.argsort(weights)[-10:]
        # plot_weights(comp_lcurves, "_WEIGHTS")


        compflux = np.nanmedian(np.array(flux)[i], axis=0)
        comp = normalise(compflux)
        target_lcurve.error(compflux, tel,ap)

        # BAD WEATHER MASK
        mask = bad_weather_mask(alc.jd, alc.nlc, 0.02, 0.08)
        target_lcurve.bw_mask = True
        if type(mask) == np.bool_:
            mask = np.zeros(len(alc.jd), dtype=bool)
        else:
            if np.count_nonzero(mask) == len(mask):
                print("Bad weather all night")
                bw = True
                # return None
        tel.bw_mask = mask
        target_lcurve.nlc_bw_mask = np.ma.masked_where(mask, target_lcurve.nlc)
        if bw == False:
            target_lcurve.nlc_bw_mask = normalise(target_lcurve.nlc_bw_mask)

        # CALCULATE PHOTOMETRIC ZP
        # phot_zp = calculate_zp(mask, comp_lcurves)

        # GET ALL FINAL DIFF LCS
        for c in range(len(comp_lcurves)):
            comp_lcurves[c].diff_photom(alc)
            comp_lcurves[c].error(compflux, tel,ap)
            # comp_lcurves[c].bin_lc(binning)
            # comp_lcurves[c].nlc_bw_mask = normalise(np.ma.masked_where(mask, comp_lcurves[c].nlc))

        # PWV CORRECTION
        # target_lcurve.targetlist_teff(tlist)

        pwv = True

        if targ_teff is None:
            if no_teff:
                try:
                    target_lcurve.filippazzo_teff(basedir)
                except:
                    print("Target is not in the Filippazzo list from Peter")
                    print("WARNING: NO TEFF")
                    pwv = False
            else:
                print("Using Teff from header: " + str(teff_header))
                target_lcurve.teff = int(teff_header)
                print(f"Target TEFF set to: {target_lcurve.teff}K")
        else:
            print("T_eff provided by user: ", str(targ_teff))
            target_lcurve.teff = int(targ_teff)
            print(f"Target TEFF set to: {target_lcurve.teff}K")

        # Ensure we don't override the TEFF value
        print(f"Final target TEFF before PWV check: {target_lcurve.teff}K")

        # warn about Teff bounds
        if target_lcurve.teff > 3200:
            print("WARNING: Effective Temperature is higher than expected: " + str(target_lcurve.teff) + "K")
        elif target_lcurve.teff < 2000:
            print("WARNING: Effective Temperature is too low for reliable PWV correction: " + str(target_lcurve.teff) + "K")
            pwv = False

        if "ARTEMIS" in tel.telname.upper() or "SAINT-EX" in tel.telname.upper():
            pwv = False

        if not np.isnan(target_lcurve.teff) and pwv:
            print("Performing PWV Correction")
            try:
                target_lcurve = pwv_correct(target_lcurve, tel, basedir, outfits,pltname, comp_lcurves,globallc)
                target_lcurve.nlc_pwv_bw_mask = np.ma.masked_where(mask, target_lcurve.nlc_pwv)
                if bw == False: target_lcurve.nlc_pwv_bw_mask = normalise(target_lcurve.nlc_pwv_bw_mask)
            except Exception as e:
                print ("PWV Failed:",e)
        else:
            print("PWV Failed: No T_eff for target")

        outputfits(comp_lcurves, target_lcurve, alc,bw, tel, os.path.dirname(outfits) + "/" + os.path.basename(pltname) + "_diff.fits", ap,version)

        # PLOT LCs
        if not os.path.exists(lcdir):
            os.mkdir(lcdir)

        if globallc:
            # target_lcurve.split_lc()
            print("Plot Global LC")
            target_lcurve.plot_multilc(binning, pltname + "_global")
            target_lcurve.plot_standard(num_comp, binning, pltname, True, False,(16, 8))
        else:
            target_lcurve.plot_standard(num_comp,binning,pltname,True, True,(8, 8))
            target_lcurve.plot_standard(num_comp, binning, pltname + "_nolim", False, True, (8, 8))
        print("_".join(pltname.split("_")[:-1]))
        target_lcurve.best_ap("_".join(pltname.split("_")[:-1]),ap,num_comp)

        target_lcurve.sclip(5)
        print("STANDARD DEV OF LC: " + str(np.std(target_lcurve.nlc_clip)))
        print("STANDARD DEV OF PWV LC: " + str(np.std(sigmaclip(target_lcurve.nlc_pwv,5))))

        # outputfits(comp_lcurves, target_lcurve, alc, tel, pltname + ".fits",ap)
        # rms_izmag(comp_lcurves,tel,target_lcurve,binning,pltname + "_quality_of_night")

        # SAVE LC FILES
        mcmc = pltname + "_MCMC"
        target_lcurve.save_mcmc_txt(tel, mcmc)
        # target_lcurve.save_npy(pltname)

        # PLOT TELESCOPE VARIABLES AGAINST TIME:
        # arg_names = ['tel_diam', 'airmass', 'altitude','azimuth','ccdtemp','dec_tel', 'exposure','filt','focuspos','gain',
        #              'humidity','ra_tel', 'RON', 'rcore', 'ra_move', 'dec_move','pa', 'fwhm', 'seeing', 'psf_a_5',
        #              'psf_b_5', 'skylevel', 'ambtemp','telname','targname','jd','bjd','hjd','date']

        print(tel.airmass[0],tel.ccdtemp[0], tel.exposure[0], tel.fwhm[0], tel.gain[0], tel.bg[0])
        # Replace the existing pltmlc calls with:
        plot_single_param(tel.jd, tel.airmass, "Airmass", orig_pltname + "_airmass", binning)
        plot_single_param(tel.jd, tel.ccdtemp, "CCD Temperature", orig_pltname + "_ccdtemp", binning)
        plot_single_param(tel.jd, tel.exposure, "Exposure Time", orig_pltname + "_exposure", binning)
        plot_single_param(tel.jd, tel.fwhm, "FWHM", orig_pltname + "_fwhm", binning)
        plot_single_param(tel.jd, tel.gain, "Gain", orig_pltname + "_gain", binning)
        plot_single_param(tel.jd, tel.bg, "Sky Level", orig_pltname + "_background", binning)
        plot_single_param(target_lcurve.jd, target_lcurve.izmag, "I+z Magnitude", orig_pltname + "_izmag", binning)

        # Add PWV comparison plot
        if hasattr(tel, 'pwv') or hasattr(tel, 'pwv_original'):
            plot_pwv_comparison(tel, target_lcurve.jd, gaiaid, filt, date,
                                orig_pltname + "_pwv_comparison")


        # background_fix(target_lcurve, alc, comp_lcurves,pltname,lcdir)
        background_fix_test(target_lcurve, alc, comp_lcurves, pltname, lcdir)


        return target_lcurve

        # else:
        #     if more than one object in field is in target list
            # print("Still to be developed...")
            # return None
    else:
        # if none of the objects in the field are in the target list (older targets)
        return None

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('outfits')
    parser.add_argument('-f', '--filt', required=True)
    parser.add_argument('-o', '--outputdir', required=True)
    parser.add_argument('-l', '--lcdir', required=True)
    parser.add_argument('-a', '--ap', required=True)
    parser.add_argument('-g', '--gaia', required=False)
    parser.add_argument('-t', '--teff', required=False)
    parser.add_argument('-d', '--date', required=False)
    parser.add_argument('-b', '--bin', required=True)
    parser.add_argument('-v', '--version', required=True)
    parser.add_argument('-targ', '--targlist', help="Current Target list",required=True)
    parser.add_argument('-oldtarg', '--oldtarglists', help="Older Target lists to check, separated by spaces", required=False)
    parser.add_argument('--no_fits',help="Don't produce any fits files",action='store_true')
    parser.add_argument('--basedir', help="Base directory path", required=False)
    args = parser.parse_args()
    # print(args)

    version = args.version
    if args.oldtarglists is not None:
        oldtarglists = args.oldtarglists.split(" ")
    else:
        oldtarglists = []

    if args.date is None:
        globallc = True
        # goutfits = args.outputdir + "/" + args.gaiaid + "_" + args.filt + "_output.fits"
        # outfits = goutfits
        outfits = args.outfits
        gaiaid = os.path.basename(outfits).split("_")[0]
        goutfits = outfits
        outdir = args.outputdir

        if args.gaia is not None:
            gaiaid = args.gaia
        lcdir = args.lcdir + "/" + gaiaid
        if not os.path.exists(lcdir):
            os.mkdir(lcdir)

    else:
        globallc = False
        # goutfits = args.outputdir + "/" + args.gaiaid + "_" + args.filt + "_output.fits"
        # outfits = args.outputdir + "/" + args.date + "/" + args.targname  + "/" + args.gaiaid + "/" + args.gaiaid + "_" + args.filt + "_" + args.date + "_output.fits"
        outfits = args.outfits
        if len(os.path.basename(outfits).split("_")) == 4:
            gaiaid = os.path.basename(outfits).split("_")[0]
        elif len(os.path.basename(outfits).split("_")) == 5:
            gaiaid = "_".join(os.path.basename(outfits).split("_")[:2])
        else:
            gaiaid = "0000000000000000000"
        goutfits = args.outputdir + "/" + gaiaid + "_" + args.filt + "_output.fits"
        # outdir = args.outputdir + "/" + args.date + "/" + args.targname  + "/" + args.gaiaid + "/"
        outdir = args.outputdir
        lcdir = args.lcdir

        if args.gaia is not None:
            gaiaid = args.gaia

        if not os.path.exists(lcdir):
            os.mkdir(lcdir)

    try:
        lc = main(args.date,gaiaid,int(args.ap), args.filt, outfits, goutfits, globallc, args.bin,version,outdir,lcdir,args.teff,args.targlist,oldtarglists,args.basedir)

    except Exception as e:
        print(e)
