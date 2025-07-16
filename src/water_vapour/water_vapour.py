import numpy as np
from scipy.interpolate import LinearNDInterpolator
from itertools import product
import argparse
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
# import json
from astropy.io import ascii
import fitsio
from astropy.time import Time
from astropy.io import fits
import math
from astropy.stats import sigma_clip as sigmaclip
import pandas as pd
from scipy import ndimage,signal
import csv
from scipy.optimize import leastsq
import os
import glob
from water_vapour import pwvGrid

def latex_plot():
    # inspired by http://nipunbatra.github.io/2014/08/latexify/
    params = {
        'text.latex.preamble': ['\\usepackage{gensymb}','\\usepackage{siunitx}',],
        'image.origin': 'lower',
        #    'image.interpolation': 'nearest',
        #    'image.cmap': 'gray',
        #    'axes.grid': False,
        #    'savefig.dpi': 150,  # to adjust notebook inline plot size
        'axes.labelsize': 22,  # fontsize for x and y labels (was 10)
        'axes.titlesize': 22,
        'font.size': 22,  # was 10
        'legend.fontsize': 22,  # was 10
        'xtick.labelsize': 22,
        'ytick.labelsize': 22,
        'text.usetex': True,
        #    'figure.figsize': [3.39, 2.10],
        'font.family': 'serif',
    }
    plt.rcParams.update(params)

def plot_graph(x,y,z,c,point):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111, projection='3d')
    p=ax.scatter(x, y, z, c=c, cmap='viridis')
    ax.set_xlabel("PWV (mm)")
    ax.set_ylabel("Airmass")
    ax.set_zlabel("Stellar Teff (K)")
    plt.colorbar(p)

    ax.scatter(point[0],point[1],point[2],marker='x',color='r')

    # for angle in range(0, 360):
    #     ax.view_init(30, angle)
    #     plt.draw()
    #     plt.pause(.001)

    plt.show()
    # plt.draw()

def plot_4d(data,pwv_values,airmass_values,temperature_values,p):

    pwv,airmass,temp,flux = [],[],[],[]

    for d in range(len(data)):
        for e in range(len(data[d])):
            for f in range(len(data[d][e])):
                pwv.append(pwv_values[d])
                airmass.append(airmass_values[e])
                temp.append(temperature_values[f])
                flux.append(data[d][e][f])

    plot_graph(pwv, airmass, temp, flux,p)


def peak_removal(pwv_time):
    peak_widths = np.arange(1, 3)
    peak_indices = signal.find_peaks_cwt(pwv_time.pwv, peak_widths)
    peak_count = len(peak_indices)  # the number of peaks in the array

    ind_array = []
    for i in peak_indices:
        ind_array.append(i + 1)
        ind_array.append(i - 1)
        ind_array.append(i - 2)
        ind_array.append(i)

    # Use lists instead of Series.append()
    pwv_datetime_list = []
    pwv_pwv_samp = []
    for i, val in enumerate(pwv_time.index):
        if i not in ind_array:
            pwv_datetime_list.append(pwv_time.datetime.iloc[i])
            pwv_pwv_samp.append(pwv_time.pwv.iloc[i])

    # Convert to Series at the end
    pwv_datetime_samp = pd.Series(pwv_datetime_list)

    sigma = 1
    # pwv_pwv_samp = ndimage.gaussian_filter1d(pwv_pwv_samp, sigma)
    return pwv_datetime_samp, pwv_pwv_samp

def interpolator(coords, data, point):
    dims = len(point)
    indices = []
    sub_coords = []
    for j in range(dims):
        # find the indices of the values the point[j] is between
        idx = np.digitize([point[j]], coords[j])[0]
        indices += [[idx - 1, idx]]
        # find the values the point[j] is between
        try:
            sub_coords += [coords[j][indices[-1]]]
        except:
            print("ERROR")
    indices = np.array([j for j in product(*indices)])
    sub_coords = np.array([j for j in product(*sub_coords)])
    # find the relative flux for the point
    sub_data = data[list(np.swapaxes(indices, 0, 1))]
    # interpolate between the values closest to the point
    li = LinearNDInterpolator(sub_coords, sub_data, rescale=False)
    return li([point])[0]

def new_interp(coords, data, pwv, airmass, Teff, method):
    from scipy.interpolate import griddata
    Teffs = coords[..., 2][0, 0]
    Teff_lower = np.max(Teffs[Teffs <= Teff])
    Teff_upper = np.min(Teffs[Teffs >= Teff])

    if Teff_lower == Teff_upper:
        x = coords[..., 0][coords[..., 2] == Teff]  # pwv
        y = coords[..., 1][coords[..., 2] == Teff]  # airmass
        z = data[coords[..., 2] == Teff]  # effect

        zin = griddata((x, y), z, (pwv, airmass), method=method)  # interpolated value
    else:
        x_lower = coords[..., 0][coords[..., 2] == Teff_lower]  # pwv
        y_lower = coords[..., 1][coords[..., 2] == Teff_lower]  # airmass
        z_lower = data[coords[..., 2] == Teff_lower]  # effect
        zin_lower = griddata((x_lower, y_lower), z_lower, (pwv, airmass),
                             method=method)  # interpolated value lower Teff

        x_upper = coords[..., 0][coords[..., 2] == Teff_upper]  # pwv
        y_upper = coords[..., 1][coords[..., 2] == Teff_upper]  # airmass
        z_upper = data[coords[..., 2] == Teff_upper]  # effect
        zin_upper = griddata((x_upper, y_upper), z_upper, (pwv, airmass),
                             method=method)  # interpolated value upper Teff

        w_lower = (Teff_upper - Teff) / (Teff_upper - Teff_lower)  # lower weight
        w_upper = (Teff - Teff_lower) / (Teff_upper - Teff_lower)  # upper weight

        zin = w_lower * zin_lower + w_upper * zin_upper  # final interpolated value

    return zin

# new_interp(coords, data, 2.4, 1, 2600, 'cubic')

def normalise(a):
    return np.ma.divide(a, np.nanmedian(a))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def convert_to_jd(time):
    for t in range(len(time)):
        jd = Time(time[t],format='isot').jd
        time[t] = jd
    return time

def optim(lc,pwv_lc):
    p0 = (0,1)
    # res = leastsq(lambda a: (lc - fct(pwv_lc, a)), p0)
    res = leastsq(lambda a: (np.ma.divide(lc,fct(pwv_lc, a))-1), p0)
    ab = res[0]

    poly = fct(pwv_lc, ab)
    opt_lc = np.ma.divide(lc,poly)
    return opt_lc, poly,ab

def fct(x, a):
    return (a[0] + (a[1] * x))

def remove_15min_peaks(t,water):
    pwv = pd.DataFrame({"pwv":water, "datetime":t})
    smooth_t, smooth_water = peak_removal(pwv)
    return smooth_t,smooth_water

def import_outputfits(fname,ap):
    with fits.open(fname) as f:
        dates = f[0].header['HISTORY']

    with fitsio.FITS(fname) as gf:
        flux = gf['FLUX_'+ap].read()
        cat = gf['catalogue']
        obj_id = cat['obj_id'].read()
        # colour = cat['g_rp'].read()
        # gmag = cat['gmag'].read()
        teff = cat['teff'].read()
        imagelist = gf['imagelist']
        airmass = imagelist['airmass'].read()
        jd = gf['jd'].read()[0]
    nflux = [normalise(f) for f in flux]

    return nflux, teff, obj_id, airmass, jd, dates

def import_weights(wfile):
    weights = ascii.read(wfile)
    star_w = np.array(weights['total weight'])
    ids = np.array(weights['ID'])
    return star_w, ids

def create_stars_dict(ids,obj_id,star_w,nflux,teff,airmass):
    stars_dict = {}
    for w in range(len(ids)):
        i = np.where(np.array(obj_id) == ids[w])[0][0]
        if stars_dict == {}:
            stars_dict = {'id':[ids[w]], 'w':[star_w[w]], 'nflux':[nflux[i]],'teff':[teff[i]], 'airmass':airmass}
        else:
            stars_dict['id'].append(ids[w])
            stars_dict['w'].append(star_w[w])
            stars_dict['nflux'].append(nflux[i])
            stars_dict['teff'].append(teff[i])

    return stars_dict

def get_weight_temp_lists(stars_dict):
    # create weight and Temp lists for only stars with a Gaia effective Temp
    w_list,t_list = [],[]

    for i in range(len(stars_dict['teff'])):
        if not np.isnan(stars_dict['teff'][i]):
            w_list.append(stars_dict['w'][i])
            t_list.append(stars_dict['teff'][i])

    w_list = w_list / np.nansum(w_list)

    return w_list, t_list

def import_pwv(pwv_file):
    import csv
    t, water =[],[]
    with open(pwv_file, 'r') as sfile:
        reader = csv.reader(sfile)
        for row in reader:
            try:
                if len(row)==5 and row[0]!="Date time":
                    t.append(row[0])
                    water.append(float(row[4]))
            except Exception as e:
                print(e)
    sfile.close()
    return t, water

def import_pwv2(pwv_file):
    import csv
    t, water =[],[]
    with open(pwv_file, 'r') as sfile:
        reader = csv.reader(sfile)
        for row in reader:
            try:
                # if len(row)==5 and row[0]!="Date time":
                t.append(float(row[0]))
                water.append(float(row[1]))
            except Exception as e:
                print(e)
    sfile.close()

    return t, water

def import_pwv_night(pwv_file):

    try:
        data = ascii.read(pwv_file)
    except:
        print("No PWV Data")
        quit()
    try:
        t = np.array(data['JD'])
        pwv = np.array(data['PRECIPITABLE_WATER_VAPOUR [mm]'])
    except:
        t = np.array(data['col2'])
        pwv = np.array(data['col8'])
    return t, pwv

def import_target_lc(lc_file):
    import csv
    jd, lc = [], []
    with open(lc_file, 'r') as sfile:
        reader = csv.reader(sfile)
        for row in reader:
            try:
                jd.append(float(row[0]))
                try:
                    lc.append(float(row[1]))
                except:
                    lc.append(np.nan)
            except Exception as e:
                print(e)
    sfile.close()
    return jd, lc


def apply_water_correction(wgrid,wdata, targ_jd, targ_lc, stars_dict,targ_Teff,airmass,targ_jd_sub,water,oldformat):
     # get list of temperature and weights for the stars which have gaia Temps
    w_list, t_list = get_weight_temp_lists(stars_dict)

    if np.count_nonzero(w_list) > 0:
        # if #LC point >> #water points then the correction will be insufficient
        # if float(len(targ_jd))/float(len(t)) < 10000:

        # estimate ALC's temp using Gaia temperatures and weights
        ALC_Teff = np.sum(np.multiply(w_list, t_list))

        alc_water, targ_water, correction_worked = find_water_flux(targ_jd_sub,water,airmass,targ_Teff,ALC_Teff,targ_jd,wgrid,wdata,oldformat)

        if correction_worked:
            # rms,rms_uncorr,corrected_y = plot_corrections(alc_water, alc, targ_water,targ_jd,targ_lc,targ_jd_sub,water,owater,jd,dates,pltname)
            y = targ_lc
            # if sclip == True: y = sigmaclip(y, sigma=4)
            y = normalise(y)
            y2 = np.ma.divide(targ_water, alc_water)
            corrected_y = np.ma.divide(y, y2)
            rms = np.nanstd(corrected_y)
            rms_uncorr = np.nanstd(y)
            return corrected_y,rms,rms_uncorr
        else:
            print("Water Correction Failed")
            return None, None, None
    else:
        print("Water Correction Failed: Not enough stars with Teff for this field")
        return None, None, None


def find_water_flux(targ_jd_sub, water, airmass, targ_Teff, ALC_Teff, jd, wgrid, wdata, oldformat):
    coords, data = load_water_grid(wgrid, wdata)
    targ_water, alc_water, pwv_closest, pwv_t_closest = [], [], [], []
    correction_worked = False

    for j in range(len(jd)):
        if (j >= targ_jd_sub[0] and j <= targ_jd_sub[1]):
            # Convert absolute index to relative index for water array
            water_idx = j - targ_jd_sub[0]

            if water_idx < len(water) and water[water_idx] >= 0.05:
                try:
                    if oldformat:
                        # ... existing code using water[water_idx] ...
                        value_int = new_interp(coords, data, water[water_idx], airmass[j], targ_Teff, 'cubic')
                    else:
                        value_int = new_interp(coords, data, water[water_idx], airmass[j], targ_Teff, 'cubic')
                    targ_water.append(value_int)

                    # Same for ALC
                    if oldformat:
                        value_int = new_interp(coords, data, water[water_idx], airmass[j], ALC_Teff, 'cubic')
                    else:
                        value_int = new_interp(coords, data, water[water_idx], airmass[j], ALC_Teff, 'cubic')
                    alc_water.append(value_int)

                    correction_worked = True
                except Exception as e:
                    print(f"ERROR: {e}")
                    targ_water.append(np.nan)
                    alc_water.append(np.nan)
            else:
                if water_idx >= len(water):
                    print(f"PWV index {water_idx} out of bounds for water array of size {len(water)}")
                else:
                    print(f"PWV Value of {water[water_idx]} mm is too low (<0.05mm)")
                targ_water.append(np.nan)
                alc_water.append(np.nan)
        else:
            targ_water.append(np.nan)
            alc_water.append(np.nan)

    alc_water = normalise(alc_water)
    targ_water = normalise(targ_water)
    # pwv_closest = np.array(pwv_closest)
    # pwv_t_closest = np.array(pwv_t_closest)


    return alc_water, targ_water, correction_worked

def load_data(dir,filter, ap, targname,targ_teff,targ,date):
    weight_file = dir + "/" + ap + "weights.dat"
    uppath = lambda _path, n: os.sep.join(_path.split(os.sep)[:-n])

    # Find the base directory (where the pipeline is running from)
    if os.path.basename(uppath(dir, 4)) == "output":
        bname = uppath(dir, 6)  # Go up to base directory
    elif os.path.basename(uppath(dir, 5)) == "output":
        bname = uppath(dir, 7)  # Go up to base directory
    else:
        # Fallback: assume we're in a standard structure
        bname = uppath(dir, 4)

    # Add the pwv subdirectory
    pwv_base_dir = os.path.join(bname, "pwv")
    print(f"Looking for PWV files in: {pwv_base_dir}")

    if date is not None:
       lc_file =  dir + "/" + targname + "_" + ap + "_" + date + "_v1global.csv"
       pltname = dir + "/" + targname + "_" + targ + "_" + date + "_PLOT_" + ap
       outputfits = uppath(dir, 3) + "/" + date + "_" + targname.upper() + "_" + filter + "_output.fits"

       # targ_jd, targ_lc = import_target_lc(lc_file)
       # targ_lc = [j for (i, j) in sorted(zip(targ_jd, targ_lc))]
       # targ_jd = sorted(targ_jd)

       # if int(date)>20190101:
       #     print "Using Nightly PWV Data."
       #     dformat = datetime.datetime.strptime(date, "%Y%m%d").strftime("%Y-%m-%d")
       #     pwv_file = bname + "/" + tel + "/technical_logs/" + dformat + "/" + tel + "_ESO_LHATPRO_log.txt"
           # t, water = import_pwv_night(pwv_file)
       # else:
       #     t, water = import_pwv2(pwv_file)
    else:
        print("Using Archive PWV Data.")
        lc_file = dir + "/" + targname + "_" + ap + "_v1global.csv"
        pltname = dir + "/" + targname + "_" + targ + "_PLOT_" + ap
        outputfits = uppath(dir, 4) + "/" + targname.upper() + "_" + filter + "_output.fits"

    targ_jd, targ_lc = import_target_lc(lc_file)
    targ_lc = [j for (i, j) in sorted(zip(targ_jd, targ_lc))]
    targ_jd = sorted(targ_jd)

    from astropy.time import Time
    tS = Time(min(targ_jd)-0.5, format='jd')
    tS.format ='isot'
    tE = Time(max(targ_jd)+0.5, format='jd')
    tE.format = 'isot'

    # dateS = datetime.datetime.strptime(min(targ_jd), "%y%j").date()
    # dateS = dateS.strftime("%Y-%m-%d")
    # try:
    #     pwvdata = getLHATPROdata(tS.value, tE.value, peakR=1, format='jd')
    # except:
    #     pwvdata = getLHATPROdata(tS.value, tE.value, peakR=0, format='jd')
    # t=pwvdata.index.values.tolist()
    # water = [x[0] for x in pwvdata.values]
    try:
        pwvdata = pwvGrid.getLHATPROdata(tS.value, tE.value, peakR=1, format='jd')
    except:
        try:
            pwvdata = pwvGrid.getLHATPROdata(tS.value, tE.value, peakR=0, format='jd')
        except:
            print(("No PWV data for: " + str(tS.value) + " - " + str(tE.value)))
            pwvdata = None

    if pwvdata is not None:
        t = pwvdata.index.values.tolist()
        water = [x[0] for x in pwvdata.values]
        #
        # # CUBIC SPLINE INTERPOLATE PWV DATA TO MATCH LC TIMES:
        # i_start = find_nearest(targ_jd,t[0]-15./1440.)
        # i_end = find_nearest(targ_jd, t[-1]+15./1440.)
        # from scipy import interpolate
        # tck = interpolate.splrep(t, water, s=0.2)
        # targ_jd_sub = (i_start,i_end)
        # water_interp = interpolate.splev(targ_jd[i_start:i_end], tck, der=0)
        # CUBIC SPLINE INTERPOLATE PWV DATA TO MATCH LC TIMES:
        i_start = find_nearest(t, targ_jd[0] - (15. / 1440.))
        i_start2 = find_nearest(targ_jd, t[i_start])
        print((targ_jd[i_start2], i_start, t[i_start], water[i_start], i_start2))
        i_end = find_nearest(t, np.array(targ_jd)[-1] + (15. / 1440.))
        i_end2 = find_nearest(targ_jd, t[i_end])
        print((np.array(targ_jd)[i_end2], i_end, t[i_end], water[i_end], i_end2))
        from scipy import interpolate
        tck = interpolate.interp1d(t, water,fill_value="extrapolate")
        targ_jd_sub = (i_start2, i_end2)
        water_interp = tck(targ_jd[i_start2:i_end2])
        print((water_interp[0], targ_jd[i_start2]))

    oldformat = False

    import glob
    # Look in the new pwv directory structure
    wdata_files = glob.glob("%s/PWV_grids/*_%s_pwvGrid_data.npy" % (pwv_base_dir, filter))
    wgrid_files = glob.glob("%s/PWV_grids/*_%s_pwvGrid_coords.npy" % (pwv_base_dir, filter))

    if not wdata_files:
        print(f"No PWV data files found for filter {filter} in {pwv_base_dir}/PWV_grids/")
        print("Available files:")
        all_files = glob.glob("%s/PWV_grids/*.npy" % pwv_base_dir)
        for f in all_files:
            print(f"  {os.path.basename(f)}")
        return None

    if not wgrid_files:
        print(f"No PWV grid files found in {pwv_base_dir}/PWV_grids/")
        return None

    wdata = wdata_files[0]
    wgrid = wgrid_files[0]

    print(f"Using PWV data: {os.path.basename(wdata)}")
    print(f"Using PWV grid: {os.path.basename(wgrid)}")


    # if len(wdata_files)
    #     print "The up-to-date PWV data for " + filter + " doesn't exist. Try finding older data..."
    #     if os.path.exists(bname + '/water_grid_'+filter+'_data.npy'):
    #         wdata = bname + '/water_grid_'+filter+'_data.npy'
    #         print "Using older data: " + wdata
    #         oldformat=True
    #     else:
    #         print "No PWV Data available"
    #         quit()
    #
    # if not os.path.exists(wgrid):
    #     print "The up-to-date PWV grid doesn't exist. Try finding older grid..."
    #     if os.path.exists(bname + '/water_grid_coords.npy'):
    #         wgrid = bname + '/water_grid_coords.npy'
    #         print "Using older grid: " + wgrid
    #         oldformat=True
    #     else:
    #         print "No PWV Grid available"
    #         quit()

    # import downloaded water data (which has already been converted to JD)
    # t,water = import_pwv2(pwv_file)

    # isolate only the relevant water time period to save time
    print((min(t), max(t)))
    print((min(targ_jd),max(targ_jd)))
    # i_start = np.where(np.array(t) < targ_jd[0])[0][-1]
    # i_end = np.where(np.array(t)>targ_jd[-1])[0][0]
    # t = t[i_start:i_end]
    # water = water[i_start:i_end]

    # t = convert_to_jd(t)

    # remove 15 minute peaks from cone survey
    # t, water = remove_15min_peaks(t, water)
    nflux, teff, obj_id, airmass, jd, dates = import_outputfits(outputfits, ap)
    # import final star weights from diff photom
    star_w, ids = import_weights(weight_file)
    stars_dict = create_stars_dict(ids, obj_id, star_w, nflux, teff, airmass)

    # for a_i in np.linspace(0,2, 20, endpoint=True):
    #     for b_i in np.linspace(0, 2, 20, endpoint=True):
    #         water = fct(water,(a_i,b_i))
    # wgrid,wdata, targ_jd, targ_lc, stars_dict,nflux,targ_Teff,airmass,water,t,jd,dates,pltname):
    # wgrid,wdata, targ_jd, targ_lc, stars_dict,nflux,targ_Teff,airmass,water,owater,jd,dates,pltname,oldformat

    # wgrid, wdata, targ_jd, targ_lc, stars_dict, targ_Teff, airmass, targ_jd_sub, water, oldformat
    corrected_y, rms, rms_uncorr = apply_water_correction(wgrid,wdata, targ_jd,targ_lc, stars_dict, targ_teff, airmass, targ_jd_sub, water_interp, oldformat)

    mcmcfile = dir + "/" + "MCMC_text_" + str(ap)
    mcmc_update(mcmcfile,targ_jd,corrected_y)

    pltname = dir + "/" + targname + "_" + ap + "_" + date
    bin_jd, bin_y, bin_e = bin_data(jd, corrected_y, 5)
    plot_standard(jd,corrected_y,bin_jd,bin_y,"JD","Relative Flux",bin_e,
                  "PWV-Corrected Lightcurve of Target %s \n binned rms = %.5f, avg error = %.5f\nAperture Size = %.1f pixels" %(targname, np.nanstd(bin_y),np.nanmean(bin_e),ap_size_pixels(int(ap),4)[0]),
                  pltname + "_final_pwvcorrected",True)


def ap_size_pixels(ap, r):
    ap_pix = [0.5,1./np.sqrt(2),1,np.sqrt(2),2,2*np.sqrt(2),4,5,6,7,8,10,12]
    a_pix = []

    a_pix.append(ap_pix[ap-1]*r)

    return a_pix

def bin_data(t, y, b):
    mins_jd = float(b)/1440.
    # start = t[0]
    t = np.array(t)
    # nextbin=t[np.absolute(t-start)>mins_jd]
    split=[]
    # time = t
    sorted_y = [x for _,x in sorted(zip(t,y))]
    sorted_t = np.array(sorted(t))
    start = sorted_t[0]
    nextbin = sorted_t[np.absolute(sorted_t - start) > mins_jd]

    while nextbin!=[]:
        # start = nextbin[0]
        start = start + mins_jd
        # ind_st = np.where(t==start)[0][0]
        # ind_st = np.argmax(time>start)
        ind_st = np.argmax(sorted_t > start)
        if len(split) > 0:
            if ind_st != split[-1]:
                split.append(ind_st)
                # time = t[ind_st:]
                time = sorted_t[ind_st:]
        else:
            split.append(ind_st)
            # time = t[ind_st:]
            time = sorted_t[ind_st:]
        nextbin = time[np.absolute(time-start) > mins_jd]
        # nextbin = time[(time - start) > mins_jd]

    # times=np.split(t,split)
    times = np.split(sorted_t, split)
    ys = np.split(sorted_y, split)
    # ys=np.split(y,split)

    bins = np.zeros(len(times))
    binned_y = np.zeros(len(times))
    binned_err = np.zeros(len(times))

    for i in range(len(times)):
        if len(ys[i])>0:
            bins[i] = times[i][0]
            binned_y[i] = np.nanmedian(ys[i])
            n = len(ys[i])
            binned_err[i] = 1.253 * np.nanstd(ys[i]) / np.sqrt(n)

    bins = bins[binned_y != 0]
    binned_err = binned_err[binned_y != 0]
    binned_y = binned_y[binned_y != 0]

    return(bins, binned_y,binned_err)


def load_water_grid(wgrid,wdata):
    # Load I+z grid
    coords = np.load(wgrid)
    # data contains the flux change for every pwv, airmass, temp combo
    data = np.load(wdata)

    # Get axis values for pwv, airmass, temperature
    # pwv_values = np.unique(coords[..., 0])
    # airmass_values = np.unique(coords[..., 1])
    # temperature_values = np.unique(coords[..., 2])

    # return pwv_values, airmass_values, temperature_values, data
    return coords, data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dir')
    parser.add_argument('-f', '--filt', required=True)
    parser.add_argument('-a', '--ap', required=True)
    parser.add_argument('-t', '--targname', required=True)
    parser.add_argument('-teff', '--targ_teff', required=True)
    parser.add_argument('-id', '--sp_id', required=True)
    parser.add_argument('-d', '--date', required=False)
    args = parser.parse_args()
    # latex_plot()
    # filter = "I+z"
    # ap = "6"
    # targname = "SP1003-0105"
    # targ_teff = 2500  # 2511 #2300
    # targ = "SP000086"  # "SP000291"#"SP000104" #"SP000354" #'SP000104' #"SP000354" # targ = "SP000476" SP0731
    filter = args.filt
    ap = args.ap
    targname = args.targname
    targ_teff = int(args.targ_teff)
    targ = args.sp_id
    date=args.date
    load_data(args.dir,filter,ap,targname,targ_teff,targ,date)
    # filter = "i"
    # wgrid = args.dir + '/water_grid_coords.npy'
    # wdata = args.dir + '/water_grid_'+filter+'_data.npy'
    #
    # pwv_values, airmass_values, temperature_values, data = load_water_grid(wgrid,wdata)
    #
    # # flux_UCD, flux_other,f_UCD = [],[],[]
    # # water_range = np.linspace(0.5,19.5,39)
    # targstar, refstar = [],[]
    # pwv_array = np.linspace(0.5,19,30)
    # for pwv in pwv_array:
    #     # TRAPPIST type star
    #     point = (pwv, 1, 2600)  # pwv, airmass, temperature
    #     value_int = interpolator((pwv_values, airmass_values, temperature_values),
    #                              data, point)
    #     targstar.append(value_int)
    #
    #     # example comparison star
    #     point = (pwv, 1, 4600)  # pwv, airmass, temperature
    #     value_int = interpolator((pwv_values, airmass_values, temperature_values),
    #                              data, point)
    #     refstar.append(value_int)
    #
    # y = np.log10(np.divide(targstar,refstar))
    # y = [i+y[0] for i in y]
    # plt.plot(pwv_array,y)
    # plt.xlabel("PWV (mm)")
    # plt.title("Differential Flux between 2600K and 4600K stars")
    # plt.ylabel("Relative Flux")
    # plt.show()
    # plt.close()

    # point = (0.5, 1, 5000)  # pwv, airmass, temperature
    # first_val_other = interpolator((pwv_values, airmass_values, temperature_values),
    #                          data, point)


    # Interpolate between points on grid!
    # for p in water_range:
    #     point = (p, 1, 2550)  # pwv, airmass, temperature
    #     value_int = interpolator((pwv_values, airmass_values, temperature_values),
    #                              data, point)
    #     flux_UCD.append((value_int - first_val)/first_val)
    #     f_UCD.append(value_int)
    #
    # for p in water_range:
    #     point = (p, 1, 5000)  # pwv, airmass, temperature
    #     value_int = interpolator((pwv_values, airmass_values, temperature_values),
    #                              data, point)
    #     flux_other.append((value_int - first_val_other)/first_val_other)

    # plt.plot(water_range,flux_UCD,'r.')
    # plt.plot(water_range,flux_other,'k.')
    # plt.show()
    # plt.close()
    #
    # plt.plot(water_range, f_UCD, 'k.')
    # plt.show()

    # plot_4d(data,pwv_values,airmass_values,temperature_values,point)

    # print value_int