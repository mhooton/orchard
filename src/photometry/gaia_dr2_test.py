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
import sqlite3
import os


# print current time to keep track of how long the querying takes to run
# print datetime.datetime.now().time()
start_time = timeit.default_timer()


# ---------------------------------------------------------------------------
# Local database connection
# ---------------------------------------------------------------------------

def _db_query(min_dec, max_dec, min_ra, max_ra):
    """
    Query the unified local Gaia database for sources within the given
    RA/dec bounding box. Returns a list of row dicts.

    A fresh SQLite connection is created and closed on each call. This is
    required for correctness under multiprocessing — SQLite connections
    cannot be safely shared across processes, so a module-level cached
    connection must not be used here.

    Inputs:
        min_dec, max_dec : float  dec limits in degrees
        min_ra,  max_ra  : float  ra  limits in degrees

    Output:
        list of dicts with keys: ra, dec, pmra, pmdec, phot_g_mean_mag,
        g_rp, bp_rp, parallax, teff_gspphot, source_id, dr2_source_id
    """
    db_path = os.environ.get('GAIADATABASEPATH',
                             '/gaia_database/gaia_dr3_unified_16jcut.db')
    conn = sqlite3.connect(db_path)

    min_dec = max(min_dec, -90.0)
    max_dec = min(max_dec, 90.0)

    arr = np.arange(np.floor(min_dec), np.ceil(max_dec) + 1, 1)
    rows = []
    for i in range(len(arr) - 1):
        shard = f"{int(arr[i])}_{int(arr[i + 1])}"
        query = (
            f"SELECT ra, dec, pmra, pmdec, phot_g_mean_mag, g_rp, bp_rp, "
            f"parallax, teff_gspphot, source_id, dr2_source_id "
            f"FROM '{shard}' "
            f"WHERE dec BETWEEN {min_dec} AND {max_dec} "
            f"AND ra BETWEEN {min_ra} AND {max_ra}"
        )
        try:
            cursor = conn.execute(query)
            cols = [d[0] for d in cursor.description]
            for row in cursor.fetchall():
                rows.append(dict(zip(cols, row)))
        except Exception as e:
            print(f"DB query error for shard {shard}: {e}")

    conn.close()
    return rows


# ---------------------------------------------------------------------------
# Crossmatch entry point
# ---------------------------------------------------------------------------

def crossmatch(fitsfile,
               logfile=None,
               n=8,
               ext='fits',
               date='20180101',
               catsrc='vizgaia3',
               use_local_db=True,
               ):
    """
    Crossmatch sources in a FITS catalogue against Gaia DR3, using either
    the local unified database (default) or the Gaia archive.

    Inputs:
        fitsfile     : str   Path to input FITS catalogue
        logfile      : str   Optional path to write log output
        n            : int   Number of parallel workers
        ext          : str   File extension ('fits')
        date         : str   Observation date YYYYMMDD
        catsrc       : str   Gaia catalogue version ('vizgaia3', 'vizgaia2')
        use_local_db : bool  If True (default), query local SQLite database.
                             If False, query the Gaia archive directly.
    """
    print('\n **Crossmatching sources with ' + str(catsrc) + '...** \n')
    if use_local_db:
        print(' ** Using local Gaia database **\n')
    else:
        print(' ** Using Gaia archive **\n')

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
            ra = cat['ra'].read()
            dec = cat['dec'].read()

        # Set Gaia epoch based on catalog version
        if catsrc == 'vizgaia3':
            gaia_epoch = 2457754.5  # J2016.0 for Gaia DR3
        elif catsrc == 'vizgaia2' or catsrc == 'vizgaia':
            gaia_epoch = 2457174.5  # J2015.5 for Gaia DR2
        else:
            # Default to DR3
            gaia_epoch = 2457754.5

        obs_epoch = Time(dt.datetime.strptime(date, "%Y%m%d")).jd
        # convert delta_t into years
        delta_t = (obs_epoch - gaia_epoch) / 365.

        pool = ThreadPool(int(n))

        # carry out a conesearch for gaia targets around each identified object
        if use_local_db:
            fn = partial(f_query_gaia, query=conesearch_local)
        else:
            fn = partial(f_query_gaia, query=conesearch_archive)

        results = pool.map(fn, zip(id, ra, dec, itertools.repeat(rad_deg), itertools.repeat(delta_t)))
        pool.close()

        crossmatch_lists = list(zip(*results))
        crossmatch_keys = ['dr3_id', 'dr2_id', 'pmra', 'pmdec', 'gmag', 'g_rp', 'bp_rp', 'parallax', 'teff_val',
                           'separation']
        crossmatch = dict([(crossmatch_keys[i], crossmatch_lists[i]) for i in range(len(crossmatch_lists))])
        print(crossmatch)

        # print output on terminal again
        if logfile != None:
            sys.stdout = oldstdout

        # Count successful matches
        teffs_valid = sum(1 for x in crossmatch['teff_val'] if x is not None and str(x) != 'nan' and not np.isnan(float(x)))
        n_dr3_targets = sum(1 for x in crossmatch['dr3_id'] if x is not None and str(x) != 'nan' and not np.isnan(float(x)))
        n_dr2_matches = sum(1 for x in crossmatch['dr2_id'] if x is not None and str(x) != 'nan' and not np.isnan(float(x)))

        print("=" * 50, file=oldstdout)
        print("GAIA CROSSMATCH REPORT", file=oldstdout)
        print("=" * 50, file=oldstdout)
        print(f"Sources identified: {len(crossmatch['dr3_id'])}", file=oldstdout)
        print(f"Sources successfully crossmatched with DR3: {n_dr3_targets}", file=oldstdout)
        print(f"DR3 objects successfully crossmatched with DR2: {n_dr2_matches}", file=oldstdout)
        print(f"Objects with valid Teff values: {teffs_valid}", file=oldstdout)
        print("*" * 50, file=oldstdout)

        perc = 100 * n_dr3_targets / float(len(ra))

        elapsed = timeit.default_timer() - start_time
        print('Total time taken: ' + str(elapsed / 60.) + ' minutes')

        # Create plot showing DR2 ID availability vs Gmag
        import matplotlib.pyplot as plt

        dr2_missing = []
        gmags = []

        for i, (dr2_id, gmag) in enumerate(zip(crossmatch['dr2_id'], crossmatch['gmag'])):
            if gmag is not None and str(gmag) != 'nan' and not np.isnan(float(gmag)):
                gmags.append(gmag)
                is_missing = (dr2_id is None or str(dr2_id) == 'nan' or np.isnan(float(dr2_id)))
                dr2_missing.append(1 if is_missing else 0)

        if gmags:
            plt.figure(figsize=(10, 6))
            plt.scatter(gmags, dr2_missing, alpha=0.6)
            plt.xlabel('Gmag')
            plt.ylabel('DR2 ID Missing (0=has DR2 ID, 1=missing DR2 ID)')
            plt.title('DR2 ID Availability vs Gmag')
            plt.grid(True, alpha=0.3)
            plt.ylim(-0.1, 1.1)

            plot_filename = fitsfile.replace('.fits', '_dr2_availability.png')
            plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"DR2 availability plot saved as: {plot_filename}")

        try:
            write_to_output(fitsfile, crossmatch, outfits)
            # Write crossmatch quality metric to Gaia_Crossmatch header
            # so that backup management can compare catalogue quality
            if not outfits:
                with fits.open(fitsfile, mode='update') as hdul:
                    ext_names = [hdu.name.upper() for hdu in hdul]
                    if 'GAIA_CROSSMATCH' in ext_names:
                        hdul['GAIA_CROSSMATCH'].header['N_MATCH'] = (
                            n_dr3_targets,
                            'Number of sources crossmatched with Gaia DR3'
                        )
                        hdul.flush()
        except Exception as e:
            print("Gaia Crossmatch Failed: ")
            print(e)
            perc = 0
            ra = []

    except Exception as e:
        print("Gaia Crossmatch Failed: ")
        print(e)
        perc = 0
        ra = []

    return perc, len(ra)


def f_query_gaia(a_b, query):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return query(*a_b)


def twoMASScrossmatch():
    a = 0
    return a


# ---------------------------------------------------------------------------
# Local database conesearch
# ---------------------------------------------------------------------------

def conesearch_local(id, ra, dec, rad_deg, delta_t):
    """
    Look up the closest Gaia DR3 source to (ra, dec) within rad_deg
    using the local unified SQLite database, applying proper motion correction.

    Inputs:
        id       : source identifier (for logging)
        ra       : float  RA in degrees
        dec      : float  Dec in degrees
        rad_deg  : float  search radius in degrees
        delta_t  : float  time baseline in years from Gaia epoch

    Output:
        tuple: (dr3_id, dr2_id, pmra, pmdec, gmag, g_rp, bp_rp,
                parallax, teff, separation)
    """
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

    r_dr3_id = np.nan
    r_dr2_id = np.nan
    r_pmra = np.nan
    r_pmdec = np.nan
    r_gmag = np.nan
    r_g_rp = np.nan
    r_bp_rp = np.nan
    r_parallax = np.nan
    r_teff = np.nan
    r_sep = np.nan

    print(id)

    rows = _db_query(dec - rad_deg, dec + rad_deg, ra - rad_deg, ra + rad_deg)

    if not rows:
        return (r_dr3_id, r_dr2_id, r_pmra, r_pmdec, r_gmag,
                r_g_rp, r_bp_rp, r_parallax, r_teff, r_sep)

    # Compute proper-motion-corrected separation for each candidate
    seps = []
    for row in rows:
        try:
            pmra = row['pmra']
            pmdec = row['pmdec']
            if pmra is not None and pmdec is not None:
                new_ra = row['ra'] + (delta_t * (pmra / 1000.) / 3600.)
                new_dec = row['dec'] + (delta_t * (pmdec / 1000.) / 3600.)
            else:
                new_ra = row['ra']
                new_dec = row['dec']
            coord_gaia = SkyCoord(new_ra, new_dec, unit=(u.deg, u.deg), frame='icrs')
            seps.append(coord.separation(coord_gaia).radian)
        except Exception:
            seps.append(np.inf)

    min_idx = int(np.argmin(seps))
    best = rows[min_idx]
    best_sep = seps[min_idx]

    print(best_sep)

    if best_sep > 0.00001:
        print('hold')
        return (r_dr3_id, r_dr2_id, r_pmra, r_pmdec, r_gmag,
                r_g_rp, r_bp_rp, r_parallax, r_teff, r_sep)

    r_dr3_id = best['source_id']
    r_dr2_id = best['dr2_source_id']
    r_pmra = best['pmra']
    r_pmdec = best['pmdec']
    r_gmag = best['phot_g_mean_mag']
    r_g_rp = best['g_rp']
    r_bp_rp = best['bp_rp']
    r_parallax = best['parallax']
    r_teff = best['teff_gspphot']
    r_sep = best_sep

    return (r_dr3_id, r_dr2_id, r_pmra, r_pmdec, r_gmag,
            r_g_rp, r_bp_rp, r_parallax, r_teff, r_sep)


# ---------------------------------------------------------------------------
# Archive conesearch (original behaviour, preserved as fallback)
# ---------------------------------------------------------------------------

def get_dr2_id_from_dr3(dr3_source_id):
    """Get DR2 source_id from DR3 source_id via the Gaia archive."""
    try:
        query = "SELECT dr2_source_id FROM gaiadr3.dr2_neighbourhood WHERE dr3_source_id = " + str(dr3_source_id)
        j = Gaia.launch_job(query=query, verbose=False)
        r = j.get_results()

        if len(r) > 0:
            return r['dr2_source_id'][0]
        else:
            return np.nan
    except Exception as e:
        print(f"Error getting DR2 ID for DR3 ID {dr3_source_id}: {e}")
        return np.nan


def conesearch_archive(id, ra, dec, rad_deg, delta_t):
    """
    Look up the closest Gaia DR3 source to (ra, dec) within rad_deg
    by querying the Gaia archive directly.

    Inputs:
        id       : source identifier (for logging)
        ra       : float  RA in degrees
        dec      : float  Dec in degrees
        rad_deg  : float  search radius in degrees
        delta_t  : float  time baseline in years from Gaia epoch

    Output:
        tuple: (dr3_id, dr2_id, pmra, pmdec, gmag, g_rp, bp_rp,
                parallax, teff, separation)
    """
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    radius = Quantity(rad_deg, u.deg)

    r_dr3_id = np.nan
    r_dr2_id = np.nan
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

    ra = coord.ra.deg
    dec = coord.dec.deg
    if rad_deg is not None:
        if hasattr(radius, 'unit'):
            radiusDeg = radius.to(u.deg).value
        else:
            radiusDeg = (radius * u.arcsec).to(u.deg).value
        query = "SELECT DISTANCE(POINT('ICRS',ra,dec), \
                POINT('ICRS'," + str(ra) + "," + str(dec) + ")) AS dist, * \
                FROM gaiadr3.gaia_source WHERE CONTAINS(\
                POINT('ICRS',ra,dec),\
                CIRCLE('ICRS'," + str(ra) + "," + str(dec) + ", " + str(radiusDeg) + "))=1 \
                ORDER BY dist ASC"

    j = Gaia.launch_job(query=query, verbose=False)
    r = j.get_results()

    # FIX ISSUE WHERE SOURCE_ID COLNAME CASE VARIES DEPENDING ON SEARCH
    if 'SOURCE_ID' in r.colnames:
        r.rename_column('SOURCE_ID', 'source_id')

    r_sub = r['ra', 'dec', 'pmra', 'pmdec', 'source_id', 'phot_g_mean_mag', 'bp_rp', 'g_rp', 'parallax', 'teff_gspphot']

    matches = len(r_sub['pmra'])

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

        r_dr3_id = r_sub['source_id'][0]
        r_dr2_id = get_dr2_id_from_dr3(r_dr3_id)
        r_pmra = r_sub['pmra'][0]
        r_pmdec = r_sub['pmdec'][0]
        r_gmag = r_sub['phot_g_mean_mag'][0]
        r_g_rp = r_sub['g_rp'][0]
        r_bp_rp = r_sub['bp_rp'][0]
        r_parallax = r_sub['parallax'][0]
        r_teff = r_sub['teff_gspphot'][0]

    elif matches > 1:
        sep = []
        oldsep = []
        for m in range(matches):
            if not np.ma.is_masked(r_sub['pmra'][m]):
                new_ra = r_sub['ra'][m] + (delta_t * (r_sub['pmra'][m] / 1000.) / 3600.)
                new_dec = r_sub['dec'][m] + (delta_t * (r_sub['pmdec'][m] / 1000.) / 3600.)
                coord_gaia = SkyCoord(new_ra, new_dec, unit=(u.deg, u.deg), frame='icrs')
                sep.append(coord.separation(coord_gaia).radian)
            else:
                coord_gaia = SkyCoord(r_sub['ra'][m], r_sub['dec'][m], unit=(u.deg, u.deg), frame='icrs')
                sep.append(coord.separation(coord_gaia).radian)

            coord_gaia = SkyCoord(r_sub['ra'][m], r_sub['dec'][m], unit=(u.deg, u.deg), frame='icrs')
            oldsep.append(coord.separation(coord_gaia).radian)

        min_sep = np.argmin(sep)
        if min_sep != np.argmin(oldsep):
            print("MISMATCH")
            for m in range(matches):
                print(r_sub['ra'][m])
                print(delta_t * (r_sub['pmra'][m] / 1000.) / 3600.)
                print(r_sub['pmra'][m])
                print(type(r_sub['pmra'][m]))
                print(r_sub['dec'][m])
                print(delta_t * (r_sub['pmdec'][m] / 1000.) / 3600.)

        print(sep[min_sep])
        r_dr3_id = r_sub['source_id'][min_sep]
        r_dr2_id = get_dr2_id_from_dr3(r_dr3_id)
        r_pmra = r_sub['pmra'][min_sep]
        r_pmdec = r_sub['pmdec'][min_sep]
        r_gmag = r_sub['phot_g_mean_mag'][min_sep]
        r_g_rp = r_sub['g_rp'][min_sep]
        r_bp_rp = r_sub['bp_rp'][min_sep]
        r_parallax = r_sub['parallax'][min_sep]
        r_teff = r_sub['teff_gspphot'][min_sep]
        r_sep = sep[min_sep]

    if r_sep > 0.00001:
        print('hold')
        r_dr3_id = np.nan
        r_dr2_id = np.nan
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

    return (r_dr3_id, r_dr2_id, r_pmra, r_pmdec, r_gmag, r_g_rp, r_bp_rp, r_parallax, r_teff, r_sep)


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def write_to_output(output, dict, outfits):
    colnames = ['GAIA_DR2_ID', 'GAIA_DR3_ID', 'parallax', 'Gmag', 'G_RP', 'BP_RP', 'Teff', 'pmra', 'pmdec']
    colnames_upper = [colname.upper() for colname in colnames]

    gaia_dr2_id = np.array([str(x) for x in dict['dr2_id']])
    gaia_dr3_id = np.array([str(x) for x in dict['dr3_id']])

    cols_upper = {'GAIA_DR2_ID': gaia_dr2_id, 'GAIA_DR3_ID': gaia_dr3_id, 'PMRA': dict['pmra'], 'PMDEC': dict['pmdec'],
                  'PARALLAX': dict['parallax'], 'GMAG': dict['gmag'], 'G_RP': dict['g_rp'], 'BP_RP': dict['bp_rp'],
                  'TEFF': dict['teff_val']}
    cols = {'GAIA_DR2_ID': gaia_dr2_id, 'GAIA_DR3_ID': gaia_dr3_id, 'pmra': dict['pmra'], 'pmdec': dict['pmdec'],
            'parallax': dict['parallax'], 'Gmag': dict['gmag'], 'G_RP': dict['g_rp'], 'BP_RP': dict['bp_rp'],
            'Teff': dict['teff_val']}

    string_cols = {'GAIA_DR2_ID', 'GAIA_DR3_ID'}

    for k, v in list(cols.items()):
        if k.upper() in string_cols:
            cols[k] = np.array(v)
        else:
            cols[k] = np.array([np.nan if x is None else x for x in v], dtype=np.float64)

    for k, v in list(cols_upper.items()):
        if k.upper() in string_cols:
            cols_upper[k] = np.array(v)
        else:
            cols_upper[k] = np.array([np.nan if x is None else x for x in v], dtype=np.float64)

    if outfits == False:
        with fits.open(output, mode='update') as hdulist_update:
            print('update gaia crossmatch data')
            name = 'Gaia_Crossmatch'
            columns = []

            for c in colnames_upper:
                print("Adding " + c + " from new crossmatch")
                if c == "GAIA_DR2_ID" or c == "GAIA_DR3_ID":
                    new_col = fits.ColDefs([fits.Column(name=c, format='26A', array=cols_upper[c])])
                    if columns == []:
                        columns = new_col
                    else:
                        columns = columns + new_col
                else:
                    new_col = fits.ColDefs([fits.Column(name=c, format='D', array=cols_upper[c])])
                    if columns == []:
                        columns = new_col
                    else:
                        columns = columns + new_col

            new_cat = fits.BinTableHDU.from_columns(columns, name=name)
            try:
                fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
            except Exception as e:
                print(e)
                print("Gaia_Crossmatch extension doesn't exist yet, so no need to delete it")
            hdulist_update.append(new_cat)
            hdulist_update.flush()

    else:
        with fits.open(output, mode='update') as hdulist_update:
            print('update catalogue data')
            name = 'CATALOGUE'
            import_cat = hdulist_update[name]
            orig_cols = import_cat.columns
            columns = []
            names = orig_cols.names

            for c in names:
                if c.upper() != "PM_RA" and c.upper() != "PM_DEC":
                    if c.upper() not in colnames_upper:
                        print("Copying " + c + " from old fits file")
                        if columns == []:
                            columns = fits.ColDefs([orig_cols[c]])
                        else:
                            columns = columns + orig_cols[c]
                    else:
                        print("Replacing " + c + " with new Crossmatch")
                        if c.upper() == "GAIA_DR2_ID":
                            new_col = fits.ColDefs([fits.Column(name=c.upper(), format='26A', array=cols_upper[c.upper()])])
                            columns = columns + new_col
                        else:
                            new_col = fits.ColDefs([fits.Column(name=c.upper(), format='D', array=cols_upper[c.upper()])])
                            columns = columns + new_col

            for c in colnames_upper:
                if c not in names:
                    print("Adding " + c + " from new crossmatch")
                    if c == "GAIA_DR2_ID" or c == "GAIA_DR3_ID":
                        new_col = fits.ColDefs([fits.Column(name=c, format='26A', array=cols_upper[c])])
                    else:
                        new_col = fits.ColDefs([fits.Column(name=c, format='D', array=cols_upper[c])])
                    columns = columns + new_col

            new_cat = fits.BinTableHDU.from_columns(columns, name=name)
            fits.hdu.hdulist.HDUList.__delitem__(hdulist_update, name)
            hdulist_update.append(new_cat)
            hdulist_update.flush()


def colappend_fits_tables(hdu1, hdu2):
    hdu1.data = np.append(hdu1.data, hdu2.data, axis=1)
    return hdu1


def run_all_targets(outputdir, logfile, n, ext):
    fitslist = glob.glob("%s/*_stack_catalogue_*.%s" % (outputdir, ext))

    for f in fitslist:
        if "pm" not in f:
            print(f)
            with fits.open(f) as fname:
                dates = fname[0].header['DATE-OBS'][:10].replace("-", "")
                print(dates)
            p = crossmatch(f, logfile, n, ext, dates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fitsfile')
    parser.add_argument('-l', '--logfile', required=False)
    parser.add_argument('-n', '--nproc', required=True)
    parser.add_argument('-e', '--ext', required=True)
    parser.add_argument('-d', '--date', required=False)
    parser.add_argument('--archive', action='store_true',
                        help='Query the Gaia archive instead of the local database')
    args = parser.parse_args()

    fitsfile = args.fitsfile
    logfile = args.logfile
    n = args.nproc
    ext = args.ext
    date = args.date
    use_local_db = not args.archive

    p = crossmatch(fitsfile, logfile, n, ext, date, use_local_db=use_local_db)


if __name__ == '__main__':
    main()