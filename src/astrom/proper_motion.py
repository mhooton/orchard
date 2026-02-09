import argparse
import fitsio
from shutil import copyfile
import math
from collections import OrderedDict
import numpy as np
import datetime
import astropy.time
import photometry.gaia_dr2_test as gaia_dr2_test

"""
Proper Motion Correction for Reference Catalogues

Updates reference catalogue star positions from Gaia DR3 epoch (J2016.0) to
observation date using Gaia proper motion measurements.

Functions:
    generate_new_cat(catfile, date, output_path)
        Main entry point for PM correction

        Args:
            catfile (str): Path to stack catalogue with Gaia_Crossmatch extension
            date (str): Target date in YYYYMMDD format
            output_path (str): Path for output PM-corrected catalogue

        Returns:
            None (writes new catalogue to output_path)

    calculate_pm(stars, cat_jd, new_jd)
        Calculate proper motion corrected positions

        Args:
            stars (dict): Star data with 'ra', 'dec', 'pmra', 'pmdec', 'parallax'
            cat_jd (float): Julian date of catalogue epoch (Gaia DR3 = 2457389.0)
            new_jd (float): Target Julian date for correction

        Returns:
            dict: Stars dictionary with 'new_ra' and 'new_dec' added

Proper motion correction:
    - Uses Gaia pmra, pmdec (mas/yr) and parallax
    - Converts PM to degrees/yr accounting for declination
    - Applies correction: new_ra = ra + (pmra/cos(dec)) * dt
    - Stars without Gaia PM data keep original positions

Example:
    generate_new_cat('SP1234567_stack_catalogue_I+z.fits', 
                     '20240115',
                     'SP1234567_stack_catalogue_pm.fits')

Typical PM values:
    - Low PM stars: <10 mas/yr → <0.1 pixel shift over 3 years
    - High PM stars: >50 mas/yr → >0.5 pixel shift over 3 years
    - Barnard's star: ~10000 mas/yr → ~100 pixel shift over 3 years

Notes:
    - Essential for accurate photometry over multi-year baselines
    - Uncorrected PM causes centroiding errors and flux loss
    - Automatically triggers Gaia cross-match if not already done
"""

# Script to calculate the Proper Motion for all stars in a Catalogue based on Gaia DR2 crossmatch
# and write a new stack catalogue with the new RA and DECs at a given time

def import_stars(cat):
    with fitsio.FITS(cat, 'rw') as catfile:
        cat = catfile[1]
        gaia = catfile['Gaia_Crossmatch']
        h = catfile[0].read_header()
        jd = h['jd-obs']
        stars = {'ra':cat['ra'].read(),'dec':cat['dec'].read(),'pmra':gaia['pmra'].read(),'pmdec':gaia['pmdec'].read(),
                 'parallax':gaia['parallax'].read(),'new_ra':[],'new_dec':[]}
    return stars, jd

def create_new_cat(old_cat,new_cat,stars):

    copy = OrderedDict()
    copyfile(old_cat,new_cat)

    with fitsio.FITS(new_cat,'rw') as cfile:
        cat = cfile[1]
        for c in cat.get_colnames():
            copy[c] = cat[c].read()

        copy['RA'] = np.array(stars['new_ra'])
        copy['DEC'] = np.array(stars['new_dec'])

        cfile[1].write(list(copy.values()),names=list(copy.keys()))
        print("Written new RA and DEC to " + str(new_cat))

def calculate_pm(stars,cat_jd,new_jd):
    for s in range(len(stars['ra'])):
        if not np.isnan(stars['pmra'][s]):

            ra = stars['ra'][s]
            dec = stars['dec'][s]

            # convert t to years:
            t = (cat_jd - float(new_jd)) / 365

            # convert from mas/yr to degrees/yr:
            pmdec = stars['pmdec'][s] / 3600000
            pmra = (stars['pmra'][s] / 3600000) / math.cos(math.radians(dec))
            # calculate new ra and dec including the PM
            new_ra = ra + (pmra * t)
            new_dec = dec + (pmdec * t)

            # can also calculate by converting to cartesian and back (gives same answer):
            # new_ra, new_dec = longmethod(stars['parallax'][s],ra,dec,stars['pmra'][s],stars['pmdec'][s],t)

            stars['new_ra'].append(new_ra)
            stars['new_dec'].append(new_dec)

        else:
            # no proper motion data from Gaia for this star
            stars['new_ra'].append(stars['ra'][s])
            stars['new_dec'].append(stars['dec'][s])

    return stars


def generate_new_cat(catfile, date, output_path):
    dt = datetime.datetime.strptime(date, "%Y%m%d")
    new_jd = astropy.time.Time(dt).jd

    try:
        stars, cat_jd = import_stars(catfile)
    except:
        print("Gaia crossmatch has not been done for this field, crossmatching now...")
        gaia_dr2_test.crossmatch(catfile, logfile=None, n=8, ext='fits', date=date)  # Add ext parameter
        stars, cat_jd = import_stars(catfile)

    stars = calculate_pm(stars, cat_jd, new_jd)

    create_new_cat(catfile, output_path, stars)

def longmethod(parallax,ra,dec,pmra,pmdec,t):
    # calculate distance in parsecs
    distance = 1. / (parallax / 1000)

    # find transverse velocities in km/sec
    vtra = (pmra / 1000) * distance * 4.740
    vtd = (pmdec / 1000) * distance * 4.740

    # find cartesian x,y,z velocities
    ra_rad = math.radians(ra)
    dec_rad = math.radians(dec)
    vx = 0 - (vtra * math.sin(ra_rad)) - (vtd * math.sin(dec_rad) * math.cos(ra_rad))
    vy = 0 + (vtra * math.cos(ra_rad)) - (vtd * math.sin(ra_rad) * math.sin(dec_rad))
    vz = 0 + (vtd * math.cos(dec_rad))

    # convert from km/yr to parsecs/yr
    vx = vx / 977780
    vy = vy / 977780
    vz = vz / 977780

    x0 = distance * math.cos(dec_rad) * math.cos(ra_rad)
    y0 = distance * math.cos(dec_rad) * math.sin(ra_rad)
    z0 = distance * math.sin(dec_rad)

    # find new positions after time t
    xt = x0 + (vx * t)
    yt = y0 + (vy * t)
    zt = z0 + (vz * t)

    dxy = math.sqrt(xt ** 2 + yt ** 2)
    new_dec = math.degrees(math.atan(zt / dxy))
    new_ra = math.degrees(math.atan2(yt, xt))

    if new_ra < 0:
        new_ra = new_ra + 360

    return new_ra, new_dec

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('catfile')
        parser.add_argument('-d','--date',required=True)
        args = parser.parse_args()

        catfile = args.catfile
        date = args.date

        generate_new_cat(catfile,date)