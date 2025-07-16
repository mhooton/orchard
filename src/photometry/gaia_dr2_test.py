from astroquery.gaia import Gaia
from astropy.coordinates.sky_coordinate import SkyCoord
import astropy.units as u
from astropy.units import Quantity
import argparse
import sys
import numpy as np
import timeit
from multiprocessing import Pool as ThreadPool
from functools import partial
import itertools
# import matplotlib.pyplot as plt
import math
import glob
import fitsio
import datetime as dt
from astropy.time import Time
from astropy.io import fits


# print current time to keep track of how long the querying takes to run
# print datetime.datetime.now().time()
start_time = timeit.default_timer()
# astroquery.gaia.


def crossmatch(fitsfile,
         logfile=None,
         n=8,
         ext='fits',
         date='20180101'):
    # set the output to a file called queryresults.txt
    print('\n **Crossmatching sources with Gaia DR2...** \n')
    if logfile != None:
        oldstdout = sys.stdout
        f = open(logfile, 'w')
        sys.stdout = f

    print("n: " + str(n))
    print("ext: " + ext)
    # choose a cone radius of approx. 5 arcseconds = 0.0014 degrees
    # rad_deg = 0.0028
    rad_deg = 0.0084
    outfits = False

    ##### POOL VERSION #####

    try:
        with fits.open(fitsfile) as fname:
            if 'output.' + ext in fitsfile:
                date = fname[0].header['HISTORY'][0]
            else:
                date = fname[0].header['DATE-OBS'][:10].replace("-", "")

        # extract ra and dec of all sources from the appropriate fits file
        with fitsio.FITS(fitsfile) as infile:
            if 'output.' + ext in fitsfile:
                cat = infile['catalogue']
                outfits = True
                id = cat['obj_id'].read()
            elif 'stack_catalogue' in fitsfile:
                cat = infile[1]
                id = cat['sequence_number'].read()
            ra = cat['ra'].read()#[100:120]
            dec = cat['dec'].read()#[100:120]
            # ids = cat['obj_id'].read()

        # gaia epoch is J2015.5
        gaia_epoch = 2457174.5
        obs_epoch = Time(dt.datetime.strptime(date, "%Y%m%d")).jd
        # convert delta_t into years
        delta_t = (obs_epoch - gaia_epoch) / 365.

        # show Gaia table columns:
        # table = Gaia.load_table('gaiadr2.gaia_source')
        # for column in (table.get_columns()):
        #     print(column.get_name())

        pool = ThreadPool(int(n))

        # carry out a conesearch for gaia targets around each identified object
        fn = partial(f_query_gaia,query=conesearch)
        results = pool.map(fn, zip(id,ra,dec,itertools.repeat(rad_deg),itertools.repeat(delta_t)))
        pool.close()

        crossmatch_lists = list(zip(*results))
        crossmatch_keys = ['source_id','pmra','pmdec','gmag','g_rp','bp_rp','parallax','teff_val','separation']
        crossmatch = dict([(crossmatch_keys[i],crossmatch_lists[i]) for i in range(len(crossmatch_lists))])
        print(crossmatch)

        elapsed = timeit.default_timer()-start_time
        print('Total time taken: ' + str(elapsed/60.) + ' minutes')

        # print output on terminal again
        if logfile != None:
            sys.stdout = oldstdout

        num_match = 0
        num_temp = 0

        try:
            for i in range(len(crossmatch['source_id'])):
                if crossmatch['source_id'][i] > 0:
                    num_match = num_match + 1
                if not np.ma.is_masked(crossmatch['teff_val'][i]):
                    num_temp = num_temp + 1

                # if ids[i] == "SP002336":
                #     print crossmatch['source_id'][i]

            #str(len(crossmatch['pmra']))
            # str(len(crossmatch['teff_val']))
            print("\n" + str(num_match) + "/" + str(len(ra)) + ' object matches with Gaia DR2\n')
            perc = 100 * num_match / float(len(ra))
            print("\n" + str(num_temp) + "/" + str(len(ra)) + ' object with Teff values from Gaia DR2\n')

            write_to_output(fitsfile, crossmatch,outfits)
        except Exception as e:
        #     print "Writing to fits file failed for: " + fitsfile
        #     print "Because crossmatch with Gaia DR2 failed"
            print(e)
            perc = 0
    except Exception as e:
        print("Gaia Crossmatch Failed: ")
        print(e)
        perc = 0
        ra = []


    ##### DEBUGGING VERSION #####
    #
    # # try:
    # with fits.open(fitsfile) as fname:
    #     if 'output.' + ext in fitsfile:
    #         date = fname[0].header['HISTORY'][0]
    #     else:
    #         date = fname[0].header['DATE-OBS'][:10].replace("-", "")
    #
    # # extract ra and dec of all sources from the appropriate fits file
    # with fitsio.FITS(fitsfile) as infile:
    #     if 'output.' + ext in fitsfile:
    #         cat = infile['catalogue']
    #         outfits = True
    #         id = cat['obj_id'].read()
    #     elif 'stack_catalogue' in fitsfile:
    #         cat = infile[1]
    #         id = cat['sequence_number'].read()
    #     ra = cat['ra'].read()#[100:120]
    #     dec = cat['dec'].read()#[100:120]
    #     # ids = cat['obj_id'].read()
    #
    # # gaia epoch is J2015.5
    # gaia_epoch = 2457174.5
    # obs_epoch = Time(dt.datetime.strptime(date, "%Y%m%d")).jd
    # # convert delta_t into years
    # delta_t = (obs_epoch - gaia_epoch) / 365.
    #
    # # show Gaia table columns:
    # # table = Gaia.load_table('gaiadr2.gaia_source')
    # # for column in (table.get_columns()):
    # #     print(column.get_name())
    #
    # results = [
    #     f_query_gaia((obj_id, ra, dec, rad_deg, delta_t), conesearch)
    #     for obj_id, ra, dec in zip(id, ra, dec)
    # ]
    #
    # crossmatch_lists = zip(*results)
    # crossmatch_keys = ['source_id','pmra','pmdec','gmag','g_rp','bp_rp','parallax','teff_val','separation']
    # crossmatch = dict([(crossmatch_keys[i],crossmatch_lists[i]) for i in range(len(crossmatch_lists))])
    # print(crossmatch)
    #
    # elapsed = timeit.default_timer()-start_time
    # print('Total time taken: ' + str(elapsed/60.) + ' minutes')
    #
    # # print output on terminal again
    # if logfile != None:
    #     sys.stdout = oldstdout
    #
    # num_match = 0
    # num_temp = 0
    #
    # try:
    #     for i in range(len(crossmatch['source_id'])):
    #         if crossmatch['source_id'][i] > 0:
    #             num_match = num_match + 1
    #         if not np.ma.is_masked(crossmatch['teff_val'][i]):
    #             num_temp = num_temp + 1
    #
    #         # if ids[i] == "SP002336":
    #         #     print crossmatch['source_id'][i]
    #
    #     #str(len(crossmatch['pmra']))
    #     # str(len(crossmatch['teff_val']))
    #     print("\n" + str(num_match) + "/" + str(len(ra)) + ' object matches with Gaia DR2\n')
    #     perc = 100 * num_match / float(len(ra))
    #     print("\n" + str(num_temp) + "/" + str(len(ra)) + ' object with Teff values from Gaia DR2\n')
    #
    #     write_to_output(fitsfile, crossmatch,outfits)
    # except Exception as e:
    # #     print "Writing to fits file failed for: " + fitsfile
    # #     print "Because crossmatch with Gaia DR2 failed"
    #     print(e)
    #     perc = 0

    return perc, len(ra)

def f_query_gaia(a_b,query):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return query(*a_b)

def twoMASScrossmatch():
    a=0
    return a

def conesearch(id, ra, dec,rad_deg,delta_t):
    # set coordinates using RA and DEC from source in output.fits
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    radius = Quantity(rad_deg, u.deg)

    r_sourceid = np.nan
    r_pmra = np.nan
    r_pmdec = np.nan
    r_gmag = np.nan
    r_flux = np.nan
    r_variable = np.nan
    r_g_rp = np.nan
    r_bp_rp = np.nan
    r_parallax = np.nan
    r_teff = np.nan
    r_sep = np.nan
    print(id)

    # gaia_dr2 = Gaia.load_table('gaiadr2.gaia_source')

    # cross match with gaia:
    ra = coord.ra.deg  # RA already in degrees
    dec = coord.dec.deg  # Dec already in degrees
    if rad_deg is not None:
        if hasattr(radius, 'unit'):
            radiusDeg = radius.to(u.deg).value
        else:
            # Assume radius is in arcsec if no unit specified
            radiusDeg = (radius * u.arcsec).to(u.deg).value
        query = "SELECT DISTANCE(POINT('ICRS',ra,dec), \
            POINT('ICRS',"+str(ra)+","+str(dec)+")) AS dist, * \
            FROM gaiadr2.gaia_source WHERE CONTAINS(\
            POINT('ICRS',ra,dec),\
            CIRCLE('ICRS',"+str(ra)+","+str(dec)+", "+str(radiusDeg)+"))=1 \
            ORDER BY dist ASC"

    j = Gaia.launch_job(query=query, verbose = False)
    r = j.get_results()

    # FIX ISSUE WHERE SOURCE_ID COLNAME CASE VARIES DEPENDING ON SEARCH
    if 'SOURCE_ID' in r.colnames:
        r.rename_column('SOURCE_ID', 'source_id')

    # isolate the columns we're interested in
    r_sub = r['ra', 'dec', 'pmra', 'pmdec', 'source_id', 'phot_g_mean_mag', 'bp_rp', 'g_rp', 'parallax', 'teff_val']


    matches = len(r_sub['pmra'])
    # print ra['ra']
    # print ra['dec']

    if matches == 1:

        if not np.ma.is_masked(r_sub['pmra'][0]):
            new_ra = r_sub['ra'][0] + (delta_t * (r_sub['pmra'][0] / 1000.) / 3600.)
            new_dec = r_sub['dec'][0] + (delta_t * (r_sub['pmdec'][0] / 1000.) / 3600.)
            coord_gaia = SkyCoord(new_ra, new_dec, unit=(u.deg, u.deg), frame='icrs')
            r_sep = (coord.separation(coord_gaia).radian)
        else:
            coord_gaia = SkyCoord(r_sub['ra'][0], r_sub['dec'][0], unit=(u.deg, u.deg), frame='icrs')
            r_sep = (coord.separation(coord_gaia).radian)

        print(r_sep)

        r_sourceid = r_sub['source_id'][0]
        r_pmra = r_sub['pmra'][0]
        r_pmdec = r_sub['pmdec'][0]
        r_gmag = r_sub['phot_g_mean_mag'][0]
        r_g_rp = r_sub['g_rp'][0]
        r_bp_rp = r_sub['bp_rp'][0]
        r_parallax = r_sub['parallax'][0]
        r_teff = r_sub['teff_val'][0]

    elif matches > 1:
        sep = []
        oldsep = []
        for m in range(matches):
            # new method
            if not np.ma.is_masked(r_sub['pmra'][m]):
                new_ra = r_sub['ra'][m] + (delta_t * (r_sub['pmra'][m]/1000.)/3600.)
                new_dec = r_sub['dec'][m] + (delta_t * (r_sub['pmdec'][m]/1000.)/3600.)
                coord_gaia = SkyCoord(new_ra, new_dec, unit=(u.deg, u.deg), frame='icrs')
                sep.append(coord.separation(coord_gaia).radian)
            else:
                coord_gaia = SkyCoord(r_sub['ra'][m], r_sub['dec'][m], unit=(u.deg, u.deg), frame='icrs')
                sep.append(coord.separation(coord_gaia).radian)

            # old method
            coord_gaia = SkyCoord(r_sub['ra'][m], r_sub['dec'][m], unit=(u.deg, u.deg), frame='icrs')
            oldsep.append(coord.separation(coord_gaia).radian)

        min_sep = np.argmin(sep)
        # print "NEW method: " + str(min_sep) + ", OLD method: " + str(np.argmin(oldsep))
        if min_sep != np.argmin(oldsep):
            print("MISMATCH")
            for m in range(matches):
                print(r_sub['ra'][m])
                print(delta_t * (r_sub['pmra'][m]/1000.)/3600.)
                print(r_sub['pmra'][m])
                print(type(r_sub['pmra'][m]))
                print(r_sub['dec'][m])
                print(delta_t * (r_sub['pmdec'][m] / 1000.) / 3600.)


        print(sep[min_sep])
        r_sourceid = r_sub['source_id'][min_sep]
        r_pmra = r_sub['pmra'][min_sep]
        r_pmdec = r_sub['pmdec'][min_sep]
        r_gmag = r_sub['phot_g_mean_mag'][min_sep]
        r_g_rp = r_sub['g_rp'][min_sep]
        r_bp_rp = r_sub['bp_rp'][min_sep]
        r_parallax = r_sub['parallax'][min_sep]
        r_teff = r_sub['teff_val'][min_sep]
        r_sep = sep[min_sep]

    if r_sep > 0.00001:
        print('hold')
        r_sourceid = np.nan
        r_pmra = np.nan
        r_pmdec = np.nan
        r_gmag = np.nan
        r_flux = np.nan
        r_variable = np.nan
        r_g_rp = np.nan
        r_bp_rp = np.nan
        r_parallax = np.nan
        r_teff = np.nan
        r_sep = np.nan


    return (r_sourceid,r_pmra,r_pmdec, r_gmag, r_g_rp, r_bp_rp, r_parallax, r_teff,r_sep)

def write_to_output(output, dict,outfits):

    # 'pm_RA', 'pm_DEC',
    colnames = ['GAIA_DR2_ID','parallax','Gmag','G_RP','BP_RP','Teff','pmra','pmdec']
    # 'pmra' 'pmdec'
    colnames_upper = [c.upper() for c in colnames]

    # gaia_id = dict['source_id']
    # for g in range(len(gaia_id)):
    #     if np.isnan(gaia_id[g]):
    #         gaia_id[g] = 0
    gaia_id = np.array([str(x) for x in dict['source_id']])

    cols_upper = {'GAIA_DR2_ID':gaia_id,'PMRA':dict['pmra'],'PMDEC':dict['pmdec'],
            'PARALLAX':dict['parallax'],'GMAG':dict['gmag'],'G_RP':dict['g_rp'],'BP_RP':dict['bp_rp'],'TEFF':dict['teff_val']}
    cols = {'GAIA_DR2_ID':gaia_id,'pmra':dict['pmra'],'pmdec':dict['pmdec'],
            'parallax':dict['parallax'],'Gmag':dict['gmag'],'G_RP':dict['g_rp'],'BP_RP':dict['bp_rp'],'Teff':dict['teff_val']}

    for k,v in list(cols.items()):
        cols[k] = np.array(v)

    for k,v in list(cols_upper.items()):
        cols_upper[k] = np.array(v)

    if outfits == False:
        with fits.open(output, mode='update') as hdulist_update:
            print('update gaia crossmatch data')
            name = 'Gaia_Crossmatch'
            # import_cat = hdulist_update[name]
            # orig_cols = import_cat.columns
            columns = []
            # names = orig_cols.names

            for c in colnames_upper:
                # if c not in names:
                print("Adding " + c + " from new crossmatch")
                if c == "GAIA_DR2_ID":
                    new_col = fits.ColDefs([fits.Column(name=c, format='26A', array=cols_upper[c])])
                    columns = new_col
                else:
                    new_col = fits.ColDefs([fits.Column(name=c, format='D', array=cols_upper[c])])
                    columns = columns + new_col

            new_cat = fits.BinTableHDU.from_columns(columns, name=name)
            try:
                fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
            except Exception as e:
                print(e)
            hdulist_update.append(new_cat)
            hdulist_update.flush()

            # with fitsio.FITS(output, 'rw') as infile:
            #     print "CAREFUL WITH CAPITALS"
            #     if 'Gaia_Crossmatch' in infile:
            #         # try:
            #         print "Gaia Crossmatch already exists, overwriting..."
            #
            #         x = infile['Gaia_Crossmatch']
            #         # for c in colnames:
            #         #     if c not in infile['Gaia_Crossmatch'].get_colnames():
            #                 # infile['Gaia_Crossmatch'].insert_column(c, cols[c])
            #                 # del cols[c]
            #
            #         infile['Gaia_Crossmatch'].write(data=cols.values(), names=cols.keys())

            # else:
            #     print "Writing Gaia_Crossmatch table to " + str(output)
            #     infile.write(data=cols.values(),names=cols.keys(),extname='Gaia_Crossmatch')
            # for n in colnames:
            #     if n not in orig_cols:
            #         # if the new columns don't already exist in the file then
            #         print "First time Gaia crossmatch, insert columns:"
            #         print n
            #         cat.insert_column(n,cols[n])
            #
            # for n in orig_cols:
            #     if n not in colnames:
            #         cols[n] = cat[n].read()
            #
            # cat.write(cols.values(), names=cols.keys())
    else:
        with fits.open(output, mode='update') as hdulist_update:
            print('update catalogue data')
            name = 'CATALOGUE'
            import_cat = hdulist_update[name]
            orig_cols = import_cat.columns
            columns = []
            names = orig_cols.names

            for c in names:
                # print c
                if c.upper()!="PM_RA" and c.upper()!="PM_DEC":
                    if c.upper() not in colnames_upper:
                        print("Copying " + c + " from old fits file")
                        # cols = cols + orig_cols[c]
                        if columns == []:
                            columns = fits.ColDefs([orig_cols[c]])
                        else:
                            # print type(orig_cols[c])
                            columns = columns + orig_cols[c]
                    else:
                        print("Replacing " + c + " with new Crossmatch")
                        if c.upper() == "GAIA_DR2_ID":
                            new_col = fits.ColDefs([fits.Column(name=c.upper(), format='26A', array=cols_upper[c.upper()])])
                            columns = columns + new_col
                        else:
                            new_col = fits.ColDefs([fits.Column(name=c.upper(), format='D',array=cols_upper[c.upper()])])
                            columns = columns + new_col

            for c in colnames_upper:
                if c not in names:
                    print("Adding " + c + " from new crossmatch")
                    new_col = fits.ColDefs([fits.Column(name=c.upper(), format='D', array=cols_upper[c.upper()])])
                    columns = columns + new_col

            new_cat = fits.BinTableHDU.from_columns(columns, name=name)
            fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
            hdulist_update.append(new_cat)
            hdulist_update.flush()


def colappend_fits_tables(hdu1,hdu2):
    # print("col append")
    # print(hdu1.data.shape[0])
    # print(hdu2.data.shape[0])
    hdu1.data = np.append(hdu1.data,hdu2.data,axis=1)
    return hdu1

def run_all_targets(outputdir,logfile,n,ext):
    fitslist = glob.glob("%s/*_stack_catalogue_*.%s" % (outputdir, ext))

    for f in fitslist:
        if "pm" not in f:
            print(f)
            with fits.open(f) as fname:
                # dates = fname[0].header['HISTORY'][0]
                dates = fname[0].header['DATE-OBS'][:10].replace("-","")
                print(dates)
            # outfitsname = f.split("/")[-1]
            # target = outfitsname.split("_")[0]
            # filter = outfitsname.split("_")[1]
            # only run for speculoos targets
            # if 'Sp' in target:
            p = crossmatch(f, logfile, n, ext, dates)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fitsfile')
    parser.add_argument('-l', '--logfile',required=False)
    parser.add_argument('-n','--nproc',required=True)
    parser.add_argument('-e', '--ext', required=True)
    parser.add_argument('-d','--date',required=False)
    args = parser.parse_args()

    fitsfile = args.fitsfile
    logfile = args.logfile
    n = args.nproc
    ext = args.ext
    date = args.date

    # run_all_targets(fitsfile,logfile,n,ext)

    p = crossmatch(fitsfile,logfile,n,ext,date)
    # run_all_targets(outputdir,logfile,n,ext)

if __name__ == '__main__':
    main()