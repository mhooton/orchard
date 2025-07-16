from astropy.io import fits
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from scipy import ndimage
import twirl
import pandas as pd
from datetime import datetime
from typing import Optional, Tuple, Union
from astropy.units import Quantity


def gaia_radecs(
        center: Union[Tuple[float, float], SkyCoord],
        fov: Union[float, Quantity],
        limit: int = 10000,
        circular: bool = True,
        tmass: bool = False,
        dateobs: Optional[datetime] = None,
        verbose: bool = False,
) -> np.ndarray:
    """
    Query the Gaia archive to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 10000.
    circular : bool, optional
        Whether to perform a circular or a rectangular query. By default, it is set to True.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.
    verbose : bool, optional
        Whether to print verbose debugging information. By default, it is set to False.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.

    Raises
    ------
    ImportError
        If the astroquery package is not installed.

    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from twirl import gaia_radecs
    >>> center = SkyCoord(ra=10.68458, dec=41.26917, unit='deg')
    >>> fov = 0.1
    >>> radecs = gaia_radecs(center, fov)
    """
    from astroquery.gaia import Gaia

    if verbose:
        print(f"[VERBOSE GAIA] Starting Gaia query...")
        print(f"[VERBOSE GAIA] Center input: {center}")
        print(f"[VERBOSE GAIA] FOV input: {fov}")
        print(f"[VERBOSE GAIA] Limit: {limit}")
        print(f"[VERBOSE GAIA] Circular: {circular}")
        print(f"[VERBOSE GAIA] tmass: {tmass}")
        print(f"[VERBOSE GAIA] dateobs: {dateobs}")

    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
        if verbose:
            print(f"[VERBOSE GAIA] Extracted RA/DEC from SkyCoord: {ra:.6f}, {dec:.6f}")
    else:
        ra, dec = center
        if verbose:
            print(f"[VERBOSE GAIA] Using RA/DEC from tuple: {ra:.6f}, {dec:.6f}")

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg
        if verbose:
            print(f"[VERBOSE GAIA] Converted FOV to Quantity: {fov}")

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
        if verbose:
            print(f"[VERBOSE GAIA] FOV array - RA: {ra_fov:.6f}, DEC: {dec_fov:.6f}")
    else:
        ra_fov = dec_fov = fov.to(u.deg).value
        if verbose:
            print(f"[VERBOSE GAIA] FOV scalar - Both: {ra_fov:.6f}")

    radius = np.min([ra_fov, dec_fov]) / 2
    if verbose:
        print(f"[VERBOSE GAIA] Calculated radius: {radius:.6f} deg")

    # Determine which query to use
    if circular and not tmass:
        query_type = "circular, no tmass"
        query = f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec
            FROM gaiadr2.gaia_source AS gaia
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY gaia.phot_g_mean_mag
            """
    elif circular and tmass:
        query_type = "circular, with tmass"
        query = f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE 1=CONTAINS(
                POINT('ICRS', {ra}, {dec}), 
                CIRCLE('ICRS', gaia.ra, gaia.dec, {radius}))
            ORDER BY tmass.j_m
            """
    elif not circular and tmass:
        query_type = "rectangular, with tmass"
        query = f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec
            FROM gaiadr2.gaia_source AS gaia
            INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = gaia.source_id
            INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = tmass_match.tmass_oid
            WHERE gaia.ra BETWEEN {ra - ra_fov / 2} AND {ra + ra_fov / 2} AND
            gaia.dec BETWEEN {dec - dec_fov / 2} AND {dec + dec_fov / 2}
            ORDER BY tmass.j_m
            """
    else:
        query_type = "rectangular, no tmass"
        query = f"""
            SELECT top {limit} gaia.ra, gaia.dec, gaia.pmra, gaia.pmdec
            FROM gaiadr2.gaia_source AS gaia
            WHERE gaia.ra BETWEEN {ra - ra_fov / 2} AND {ra + ra_fov / 2} AND
            gaia.dec BETWEEN {dec - dec_fov / 2} AND {dec + dec_fov / 2}
            ORDER BY gaia.phot_g_mean_mag
            """

    if verbose:
        print(f"[VERBOSE GAIA] Query type: {query_type}")
        print(f"[VERBOSE GAIA] SQL Query: {query.strip()}")
        print(f"[VERBOSE GAIA] Launching Gaia job...")

    try:
        job = Gaia.launch_job(query)
        if verbose:
            print(f"[VERBOSE GAIA] Job launched successfully")

        table = job.get_results()
        if verbose:
            print(f"[VERBOSE GAIA] Query returned {len(table)} results")
            if len(table) > 0:
                print(f"[VERBOSE GAIA] First few results:")
                for i in range(min(5, len(table))):
                    print(f"[VERBOSE GAIA]   {i}: RA={table['ra'][i]:.6f}, DEC={table['dec'][i]:.6f}")
    except Exception as e:
        if verbose:
            print(f"[VERBOSE GAIA] ERROR during Gaia query: {e}")
        raise

    # add proper motion to ra and dec
    if dateobs is not None:
        if verbose:
            print(f"[VERBOSE GAIA] Applying proper motion correction for date: {dateobs}")

        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25  # type: ignore

        years = dateobs - 2015.5  # type: ignore
        if verbose:
            print(f"[VERBOSE GAIA] Years from Gaia epoch (2015.5): {years:.3f}")

        # Apply proper motion corrections
        ra_corrections = years * table["pmra"] / 1000 / 3600
        dec_corrections = years * table["pmdec"] / 1000 / 3600

        table["ra"] += ra_corrections
        table["dec"] += dec_corrections

        if verbose:
            print(f"[VERBOSE GAIA] Applied proper motion corrections")
            print(f"[VERBOSE GAIA] Max RA correction: {np.max(np.abs(ra_corrections)):.6f} deg")
            print(f"[VERBOSE GAIA] Max DEC correction: {np.max(np.abs(dec_corrections)):.6f} deg")

    result = np.array([table["ra"].value.data, table["dec"].value.data]).T
    if verbose:
        print(f"[VERBOSE GAIA] Final result shape: {result.shape}")
        print(f"[VERBOSE GAIA] Gaia query completed successfully")

    return result


def twirl_wcs(filepath, verbose=False):
    if verbose:
        print(f"[VERBOSE] Starting twirl_wcs for file: {filepath}")

    with fits.open(filepath, mode='update') as hdu:
        header = hdu[0].header
        data = hdu[0].data

        if verbose:
            print(f"[VERBOSE] FITS file opened successfully")
            print(f"[VERBOSE] Data shape: {data.shape}")
            print(f"[VERBOSE] IMAGETYP: {header.get('IMAGETYP', 'NOT_FOUND')}")

        if header['IMAGETYP'] == 'Light Frame':
            if verbose:
                print(f"[VERBOSE] Confirmed Light Frame, starting image cleaning...")

            # clean image
            sigma_clip = SigmaClip(sigma=3.0)
            bkg_estimator = MedianBackground()
            if verbose:
                print(f"[VERBOSE] Creating background estimator...")

            bkg = Background2D(data, (32, 32), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            if verbose:
                print(f"[VERBOSE] Background estimation complete")

            bkg_clean = data - bkg.background
            med_clean = ndimage.median_filter(bkg_clean, size=5, mode='mirror')
            band_corr = np.median(med_clean, axis=1).reshape(-1, 1)
            image_clean = med_clean - band_corr
            if verbose:
                print(f"[VERBOSE] Image cleaning complete, cleaned image shape: {image_clean.shape}")

            # get RA and DEC from header
            ra = header['RA']
            dec = header['DEC']
            if verbose:
                print(f"[VERBOSE] Header RA: {ra}, DEC: {dec}")
            if header['TELESCOP'] == 'Artemis':
                center = SkyCoord(ra, dec, unit=[u.hourangle, u.deg])
            else:
                center = SkyCoord(ra, dec, unit=[u.deg, u.deg])
            print("COORDINATES AFTER SKYCOORD")
            print(center.ra.deg)
            print(center.dec.deg)
            if verbose:
                print(f"[VERBOSE] Center coordinates: RA={center.ra.deg:.6f} deg, DEC={center.dec.deg:.6f} deg")

            # image fov
            shape = image_clean.shape
            plate_scale = np.arctan((header['XPIXSZ'] * 1e-6) / (header['FOCALLEN'] * 1e-3)) * (180 / np.pi)
            fovx = (1 / np.abs(np.cos(center.dec.rad))) * shape[0] * plate_scale
            fovy = shape[1] * plate_scale
            if verbose:
                print(f"[VERBOSE] Plate scale: {plate_scale * 3600:.3f} arcsec/pixel")
                print(f"[VERBOSE] FOV: X={fovx * 60:.2f} arcmin, Y={fovy * 60:.2f} arcmin")
                print(f"[VERBOSE] Search FOV: X={1.2 * fovx * 60:.2f} arcmin, Y={1.2 * fovy * 60:.2f} arcmin")

            # detect and match stars
            if verbose:
                print(f"[VERBOSE] Starting star detection...")

            stars = twirl.find_peaks(image_clean, threshold=5)
            if verbose:
                print(f"[VERBOSE] Found {len(stars)} stars in image")
                if len(stars) > 0:
                    print(f"[VERBOSE] First few stars: {stars[:min(5, len(stars))]}")

            if len(stars) < 4:
                raise Exception("Not enough stars detected for plate solve")

            # Limit number of stars to prevent memory issues
            star_limit = min(20, len(stars))  # increased from 12 but still limited
            stars = stars[0:star_limit]
            if verbose:
                print(
                    f"[VERBOSE] Using {len(stars)} stars (limited from {len(twirl.find_peaks(image_clean, threshold=5))})")

            # Limit number of Gaia stars too
            if verbose:
                print(f"[VERBOSE] Querying Gaia catalog...")
                print(f"[VERBOSE] Query center: {center}")
                print(f"[VERBOSE] Query FOV: {1.2 * np.array([fovx, fovy])}")

            gaias = gaia_radecs(center, 1.2 * np.array([fovx, fovy]), verbose=verbose)
            if verbose:
                print(f"[VERBOSE] Gaia query returned {len(gaias)} stars")
                if len(gaias) > 0:
                    print(f"[VERBOSE] First few Gaia stars: {gaias[:min(5, len(gaias))]}")

            gaia_limit = min(20, len(gaias))  # roughly 2x star_limit
            gaias = gaias[0:gaia_limit]
            if verbose:
                print(f"[VERBOSE] Using {len(gaias)} Gaia stars (limited from {len(gaia_radecs(center, 1.2 * np.array([fovx, fovy]))) if len(gaia_radecs(center, 1.2 * np.array([fovx, fovy]))) <= gaia_limit else 'original query'})")

            print(f"Using {len(stars)} image stars and {len(gaias)} Gaia stars for WCS solution...")

            if verbose:
                print(f"[VERBOSE] Starting WCS computation...")

            wcs = twirl.compute_wcs(stars, gaias)
            if verbose:
                print(f"[VERBOSE] WCS computation complete")
                print(f"[VERBOSE] WCS object: {wcs}")

            # convert gaia stars to pixel coordinates for validation
            if verbose:
                print(f"[VERBOSE] Converting Gaia stars to pixel coordinates for validation...")

            gaias_pixel = np.array(SkyCoord(gaias, unit="deg").to_pixel(wcs)).T
            if verbose:
                print(f"[VERBOSE] Converted {len(gaias_pixel)} Gaia stars to pixel coordinates")
                if len(gaias_pixel) > 0:
                    print(f"[VERBOSE] First few Gaia pixel coords: {gaias_pixel[:min(5, len(gaias_pixel))]}")

            # validate star matching
            if verbose:
                print(f"[VERBOSE] Validating star matches...")

            count = 0
            matches = []
            for x, y in gaias_pixel:
                for i, j in stars:
                    distance = np.sqrt((x - i) ** 2 + (y - j) ** 2)
                    if distance < 10:  # if within 10 pixels
                        count += 1
                        matches.append((x, y, i, j, distance))
                        break  # Only count first match for each Gaia star

            if verbose:
                print(f"[VERBOSE] Found {count} star matches within 10 pixels")
                if matches:
                    print(f"[VERBOSE] Match examples (gaia_x, gaia_y, star_x, star_y, distance):")
                    for match in matches[:min(5, len(matches))]:
                        print(f"[VERBOSE]   {match}")

            if count < 4:
                raise Exception("Plate solve failed, not enough stars matched")

            # Update header with WCS information
            if verbose:
                print(f"[VERBOSE] Updating FITS header with WCS information...")

            wcs_header = wcs.to_header()
            if verbose:
                print(f"[VERBOSE] Generated WCS header with {len(wcs_header)} keywords")
                for key in wcs_header:
                    print(f"[VERBOSE]   {key}: {wcs_header[key]}")

            header.update(wcs_header)

            # Create diagnostic plots if verbose
            if verbose:
                print(f"[VERBOSE] Creating diagnostic plots...")
                _create_diagnostic_plots(filepath, data, image_clean, stars, gaias_pixel, wcs)

            if verbose:
                print(f"[VERBOSE] WCS solution completed successfully!")
        else:
            if verbose:
                print(f"[VERBOSE] Skipping - not a Light Frame (IMAGETYP={header.get('IMAGETYP')})")


def _create_diagnostic_plots(filepath, data, image_clean, stars, gaias_pixel, wcs):
    """Create diagnostic plots showing detected stars and Gaia catalog matches"""
    import matplotlib.pyplot as plt
    from photutils.aperture import CircularAperture
    import os

    # Get base filename without extension
    base_path = os.path.splitext(filepath)[0]

    # Create plot with original data
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.imshow(data, vmin=400, vmax=700, cmap="Greys_r")

    # Plot detected stars (yellow circles, radius 10)
    if len(stars) > 0:
        detected_apertures = CircularAperture(stars, r=10.0)
        detected_apertures.plot(color='yellow', alpha=0.8, linewidth=2)

    # Plot Gaia stars (red circles, radius 5)
    if len(gaias_pixel) > 0:
        gaia_apertures = CircularAperture(gaias_pixel, r=5.0)
        gaia_apertures.plot(color='red', alpha=0.8, linewidth=2)

    ax.set_title(f'WCS Diagnostic - Original Data\n{os.path.basename(filepath)}')
    plt.tight_layout()
    plt.savefig(f"{base_path}_wcs.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Create plot with cleaned data
    fig, ax = plt.subplots(figsize=(12, 10))
    correction = np.median(data-image_clean)
    ax.imshow(image_clean, vmin=400-correction, vmax=700-correction, cmap="Greys_r")

    # Plot detected stars (yellow circles, radius 10)
    if len(stars) > 0:
        detected_apertures = CircularAperture(stars, r=10.0)
        detected_apertures.plot(color='yellow', alpha=0.8, linewidth=2)

    # Plot Gaia stars (red circles, radius 5)
    if len(gaias_pixel) > 0:
        gaia_apertures = CircularAperture(gaias_pixel, r=5.0)
        gaia_apertures.plot(color='red', alpha=0.8, linewidth=2)

    ax.set_title(f'WCS Diagnostic - Cleaned Data\n{os.path.basename(filepath)}')
    plt.tight_layout()
    plt.savefig(f"{base_path}_wcsclean.png", dpi=150, bbox_inches='tight')
    plt.close()

    print(f"[VERBOSE] Saved diagnostic plots: {base_path}_wcs.png and {base_path}_wcsclean.png")