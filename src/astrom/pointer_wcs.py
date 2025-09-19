from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Union
import astropy.units as u
import numpy as np
import pandas as pd
import twirl
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.units import Quantity
from astropy.wcs.utils import WCS, pixel_to_skycoord
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from scipy import ndimage
import sqlite3

import warnings
from astropy.utils.exceptions import AstropyWarning
# Suppress specific photutils and astropy warnings
warnings.filterwarnings("ignore", message="Input data contains invalid values.*", category=AstropyWarning)
warnings.filterwarnings("ignore", message=".*Input data contains invalid values.*", module="photutils.background.background_2d")
warnings.filterwarnings("ignore", message=".*Input data contains invalid values.*", module="astropy.stats.sigma_clipping")


def find_stars_dao(
        data: np.ndarray, threshold: float = 5.0, fwhm: float = 3.0
) -> np.ndarray:
    """
    Find stars using DAOStarFinder algorithm.

    Parameters
    ----------
    data : np.ndarray
        The 2D image data
    threshold : float, optional
        Detection threshold in units of background standard deviation
    fwhm : float, optional
        Expected FWHM of stars in pixels

    Returns
    -------
    np.ndarray
        Array of (x, y) coordinates sorted by brightness
    """
    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # Use DAOStarFinder for star detection
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std)
    sources = daofind(data - median)

    if sources is None or len(sources) == 0:
        return np.array([]).reshape(0, 2)

    # Convert to (x, y) coordinates
    coordinates = np.column_stack([sources["xcentroid"], sources["ycentroid"]])

    # Sort by flux (brightness)
    fluxes = sources["flux"]
    return coordinates[np.argsort(fluxes)[::-1]]


def remove_duplicates(detections: np.ndarray, width: int, height: int) -> np.ndarray:
    """
    Remove duplicate detections using priority-based selection.

    Parameters
    ----------
    detections : np.ndarray
        Array of detections with [x, y, scale, brightness] for each detection
    width, height : int
        Image dimensions

    Returns
    -------
    np.ndarray
        Filtered array of unique detections [x, y, scale, brightness]
    """
    if len(detections) <= 1:
        return detections

    # Calculate priority scores (prefer central position and smaller scale)
    center = np.array([width / 2, height / 2])
    positions = detections[:, :2]
    scales = detections[:, 2]

    center_distances = np.linalg.norm(positions - center, axis=1)
    center_scores = 1.0 / (1.0 + center_distances / max(width, height))
    scale_scores = 1.0 / scales  # Prefer smaller scales

    priority_scores = center_scores * scale_scores

    # Sort by priority (highest first)
    sorted_indices = np.argsort(priority_scores)[::-1]
    sorted_detections = detections[sorted_indices]

    # Remove duplicates (within 10 pixels)
    unique_detections = []
    for detection in sorted_detections:
        is_duplicate = False
        pos = detection[:2]

        for existing in unique_detections:
            if np.linalg.norm(pos - existing[:2]) < 10.0:
                is_duplicate = True
                break

        if not is_duplicate:
            unique_detections.append(detection)

    return np.array(unique_detections)


def find_stars_multiscale(
        data: np.ndarray,
        scales: list = [1, 2, 3],
        threshold: float = 5.0,
        edge_buffer: int = 15,
) -> np.ndarray:
    """
    Multi-scale star detection for noisy astronomical images.

    Parameters
    ----------
    data : np.ndarray
        The 2D image data
    scales : list, optional
        List of smoothing scales to try (in pixels)
    threshold : float, optional
        Detection threshold in units of background standard deviation
    edge_buffer : int, optional
        Minimum distance from image edges to accept a detection

    Returns
    -------
    np.ndarray
        Array of (x, y) coordinates of detected stars, sorted by brightness
    """
    all_detections = []
    height, width = data.shape

    for scale in scales:
        # Smooth the image with proper edge handling
        smoothed = ndimage.gaussian_filter(data, sigma=scale, mode="reflect")

        # Detect stars
        stars = find_stars_dao(smoothed, threshold=threshold, fwhm=scale * 2)

        if len(stars) > 0:
            # Filter out stars too close to edges
            valid_stars = []
            for star in stars:
                x, y = star
                if (
                        edge_buffer <= x < width - edge_buffer
                        and edge_buffer <= y < height - edge_buffer
                ):
                    brightness = data[int(y), int(x)]
                    valid_stars.append([x, y, scale, brightness])

            if valid_stars:
                all_detections.extend(valid_stars)

    if not all_detections:
        return np.array([]).reshape(0, 2)

    # Convert to numpy array for easier manipulation
    all_detections = np.array(all_detections)

    # Remove duplicates
    unique_stars = remove_duplicates(all_detections, width, height)

    # Sort by brightness (descending order)
    sorted_indices = np.argsort(unique_stars[:, 3])[::-1]
    sorted_stars = unique_stars[sorted_indices]

    # Return only x, y coordinates
    return sorted_stars[:, :2]


def db_query(
        db: Union[str, Path], min_dec: float, max_dec: float, min_ra: float, max_ra: float
) -> pd.DataFrame:
    """
    Queries a federated database for astronomical data within a specified range of declination and right ascension.

    Parameters
    ----------
    db : str or Path
        The path to the SQLite database file
    min_dec : float
        The minimum declination value to query (degrees)
    max_dec : float
        The maximum declination value to query (degrees)
    min_ra : float
        The minimum right ascension value to query (degrees)
    max_ra : float
        The maximum right ascension value to query (degrees)

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the queried astronomical data
    """
    conn = sqlite3.connect(db)

    if min_dec < -90:
        min_dec = -90

    if max_dec > 90:
        max_dec = 90

    # Determine the relevant shard(s) based on the query parameters.
    arr = np.arange(np.floor(min_dec), np.ceil(max_dec) + 1, 1)
    relevant_shard_ids = set()
    for i in range(len(arr) - 1):
        shard_id = f"{arr[i]:.0f}_{arr[i + 1]:.0f}"
        relevant_shard_ids.add(shard_id)

    # Execute the federated query across the relevant shard(s).
    df_total = pd.DataFrame()
    for shard_id in relevant_shard_ids:
        shard_table_name = f"{shard_id}"
        q = f"SELECT * FROM `{shard_table_name}` WHERE dec BETWEEN {min_dec} AND {max_dec} AND ra BETWEEN {min_ra} AND {max_ra}"
        df = pd.read_sql_query(q, conn)
        df_total = pd.concat([df, df_total], axis=0)

    # Close the conn and return the results.
    conn.close()
    return df_total


def gaia_db_query(
        center: Union[Tuple[float, float], SkyCoord],
        fov: Union[float, Quantity],
        limit: int = 1000,
        tmass: bool = False,
        dateobs: Optional[datetime] = None,
        db_path: Union[str, Path] = "/gaia_database/gaia_tmass_16_jm_cut.db"
) -> np.ndarray:
    """
    Query the local Gaia database to retrieve the RA-DEC coordinates of stars within a given field-of-view (FOV) centered on a given sky position.

    Parameters
    ----------
    center : tuple or astropy.coordinates.SkyCoord
        The sky coordinates of the center of the FOV. If a tuple is given, it should contain the RA and DEC in degrees.
    fov : float or astropy.units.Quantity
        The field-of-view of the FOV in degrees. If a float is given, it is assumed to be in degrees.
    limit : int, optional
        The maximum number of sources to retrieve from the Gaia archive. By default, it is set to 1000.
    tmass : bool, optional
        Whether to retrieve the 2MASS J magnitudes catelog. By default, it is set to False.
    dateobs : datetime.datetime, optional
        The date of the observation. If given, the proper motions of the sources will be taken into account. By default, it is set to None.

    Returns
    -------
    np.ndarray
        An array of shape (n, 2) containing the RA-DEC coordinates of the retrieved sources in degrees.
    """
    if isinstance(center, SkyCoord):
        ra = center.ra.deg
        dec = center.dec.deg
    else:
        ra, dec = center

    if not isinstance(fov, u.Quantity):
        fov = fov * u.deg

    if fov.ndim == 1:
        ra_fov, dec_fov = fov.to(u.deg).value
    else:
        ra_fov = fov[0].to(u.deg).value
        dec_fov = fov[1].to(u.deg).value

    min_dec = dec - dec_fov / 2
    max_dec = dec + dec_fov / 2
    min_ra = ra - ra_fov / 2
    max_ra = ra + ra_fov / 2

    table = db_query(
        db_path, min_dec, max_dec, min_ra, max_ra
    )
    if tmass:
        table = table.sort_values(by=["j_m"]).reset_index(drop=True)
    else:
        table = table.sort_values(by=["phot_g_mean_mag"]).reset_index(drop=True)

    # Fix for pandas FutureWarning - separate replace and infer_objects calls
    table = table.replace("", np.nan)
    table = table.infer_objects(copy=False)
    table.dropna(inplace=True)

    # limit number of stars
    table = table[0:limit]

    # add proper motion to ra and dec
    if dateobs is not None:
        # calculate fractional year
        dateobs = dateobs.year + (dateobs.timetuple().tm_yday - 1) / 365.25  # type: ignore

        years = dateobs - 2015.5  # type: ignore
        table["ra"] += years * table["pmra"] / 1000 / 3600
        table["dec"] += years * table["pmdec"] / 1000 / 3600

    return np.array([table["ra"].values, table["dec"].values]).T

def clean_image(data: np.ndarray) -> np.ndarray:
    """
    Clean the image by removing background and applying filters.

    Parameters
    ----------
    data : np.ndarray
        The raw image data

    Returns
    -------
    np.ndarray
        The cleaned image data
    """
    # Background subtraction
    bkg = Background2D(
        data,
        (32, 32),
        filter_size=(3, 3),
        sigma_clip=SigmaClip(sigma=3),
        bkg_estimator=MedianBackground(),
    )

    bkg_clean = data - bkg.background

    # Median filter to reduce noise
    med_clean = ndimage.median_filter(bkg_clean, size=5, mode="mirror")

    # Band correction (remove horizontal stripes)
    band_corr = np.median(med_clean, axis=1).reshape(-1, 1)
    image_clean = med_clean - band_corr

    # Clip negative values
    image_clean = np.clip(image_clean, 1, None)

    return image_clean


@dataclass
class PointingCorrection:
    """Class to store the pointing correction between the desired target center and the plating center."""
    target_ra: float
    target_dec: float
    plating_ra: float
    plating_dec: float

    @property
    def offset_ra(self):
        return self.plating_ra - self.target_ra

    @property
    def offset_dec(self):
        return self.plating_dec - self.target_dec

    @property
    def angular_separation(self) -> float:
        desired_center = SkyCoord(self.target_ra, self.target_dec, unit=[u.deg, u.deg])
        plating_center = SkyCoord(
            self.plating_ra, self.plating_dec, unit=[u.deg, u.deg]
        )
        return desired_center.separation(plating_center).deg


@dataclass
class ImageStarMapping:
    """
    A class to handle the mapping of stars detected in an image to their corresponding
    Gaia star coordinates using World Coordinate System (WCS) transformations.
    """
    wcs: WCS
    stars_in_image: np.ndarray
    gaia_stars_in_image: np.ndarray

    @classmethod
    def from_gaia_coordinates(cls, stars_in_image: np.ndarray, gaia_stars: np.ndarray):
        wcs = twirl.compute_wcs(stars_in_image, gaia_stars)
        gaia_stars_in_image = np.array(SkyCoord(gaia_stars, unit="deg").to_pixel(wcs)).T
        return cls(wcs, stars_in_image, gaia_stars_in_image)

    def get_plating_center(self, image_shape: Tuple[int, int]) -> Tuple[float, float]:
        plating_center = pixel_to_skycoord(
            image_shape[1] / 2, image_shape[0] / 2, self.wcs
        )
        return float(plating_center.ra.deg), float(plating_center.dec.deg)

    def find_gaia_match(self):
        squared_distances = np.sum(
            (
                    self.stars_in_image[:, np.newaxis]
                    - self.gaia_stars_in_image[np.newaxis, :, :]
            )
            ** 2,
            axis=-1,
        )
        match_index = np.argmin(squared_distances, axis=1)
        distance = np.sqrt(np.min(squared_distances, axis=1))

        return self.gaia_stars_in_image[match_index], distance

    def number_of_matched_stars(self, pixel_threshold: int = 10):
        distance_to_closest_star = self.find_gaia_match()[1]
        return np.sum(distance_to_closest_star < pixel_threshold)


def _create_diagnostic_plots(filepath, data, image_clean, stars, gaias_pixel, wcs):
    """Create diagnostic plots showing detected stars and Gaia catalog matches"""
    try:
        import matplotlib.pyplot as plt
        from astropy.visualization import ZScaleInterval
        import os

        # Get base filename without extension
        base_path = os.path.splitext(filepath)[0]

        # Create plot with original data
        fig, ax = plt.subplots(figsize=(10, 10))

        interval = ZScaleInterval(contrast=0.05)
        vmin, vmax = interval.get_limits(data)
        ax.imshow(data, cmap="viridis", origin="lower", vmin=vmin, vmax=vmax)

        # Plot detected stars (red x markers)
        if len(stars) > 0:
            ax.scatter(
                stars[:, 0],
                stars[:, 1],
                c="red",
                s=300,
                marker="x",
                linewidths=2,
                alpha=0.8,
                label=f"{len(stars)} detected stars",
            )

        # Plot Gaia stars (blue circles)
        if len(gaias_pixel) > 0:
            ax.scatter(
                gaias_pixel[:, 0],
                gaias_pixel[:, 1],
                s=150,
                marker="o",
                facecolors="none",
                edgecolors="blue",
                linewidths=2,
                alpha=0.8,
                label=f"{len(gaias_pixel)} Gaia stars",
            )

        ax.set_title(f'WCS Diagnostic - Original Data\n{os.path.basename(filepath)}')
        ax.set_xlabel("X pixel")
        ax.set_ylabel("Y pixel")
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{base_path}_wcs.png", dpi=150, bbox_inches='tight')
        plt.close()

        # Create plot with cleaned data
        fig, ax = plt.subplots(figsize=(10, 10))

        vmin, vmax = interval.get_limits(image_clean)
        ax.imshow(image_clean, cmap="viridis", origin="lower", vmin=vmin, vmax=vmax)

        # Plot detected stars (red x markers)
        if len(stars) > 0:
            ax.scatter(
                stars[:, 0],
                stars[:, 1],
                c="red",
                s=300,
                marker="x",
                linewidths=2,
                alpha=0.8,
                label=f"{len(stars)} detected stars",
            )

        # Plot Gaia stars (blue circles)
        if len(gaias_pixel) > 0:
            ax.scatter(
                gaias_pixel[:, 0],
                gaias_pixel[:, 1],
                s=150,
                marker="o",
                facecolors="none",
                edgecolors="blue",
                linewidths=2,
                alpha=0.8,
                label=f"{len(gaias_pixel)} Gaia stars",
            )

        ax.set_title(f'WCS Diagnostic - Cleaned Data\n{os.path.basename(filepath)}')
        ax.set_xlabel("X pixel")
        ax.set_ylabel("Y pixel")
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{base_path}_wcsclean.png", dpi=150, bbox_inches='tight')
        plt.close()

        print(f"Saved diagnostic plots: {base_path}_wcs.png and {base_path}_wcsclean.png")

    except ImportError:
        print("Matplotlib not available - skipping diagnostic plots")
    except Exception as e:
        print(f"Error creating diagnostic plots: {e}")


def clear_wcs_headers(header, wcs_keywords, verbose=False):
    """
    Remove all existing WCS headers from a FITS header.

    Parameters
    ----------
    header : astropy.io.fits.Header
        FITS header to clean
    wcs_keywords : list
        List of WCS keyword names to remove
    verbose : bool, optional
        If True, print which headers are being removed

    Returns
    -------
    int
        Number of headers removed
    """
    removed_count = 0
    for keyword in wcs_keywords:
        if keyword in header:
            if verbose:
                print(f"Removing existing WCS header: {keyword} = {header[keyword]}")
            del header[keyword]
            removed_count += 1
    return removed_count


def pointer_wcs(filepath, db_path, wcs_keywords=None, clear_existing_wcs=False, verbose=False):
    """
    Perform WCS solving on a FITS file using local Gaia database and multiscale star detection.

    This is a drop-in replacement for twirl_wcs that uses:
    - Local Gaia database instead of online queries
    - Multiscale star detection instead of simple peak finding
    - No dark frame subtraction (assumes pre-reduced images)

    Parameters
    ----------
    filepath : str
        Path to the FITS file to plate solve
    db_path : str
        Path to the Gaia database file
    wcs_keywords : list, optional
        List of WCS keywords for header operations
    clear_existing_wcs : bool, optional
        If True, remove existing WCS headers before plate solving
    verbose : bool, optional
        Whether to print verbose debugging information

    Returns
    -------
    dict
        Dictionary containing results:
        - success: bool, whether WCS solving succeeded
        - error: str or None, error message if failed
        - crpix: tuple or None, (CRPIX1, CRPIX2) if successful
        - crval: tuple or None, (CRVAL1, CRVAL2) if successful
        - sources_detected: int, number of sources detected in image
        - sources_used: int, number of sources used for WCS (after limiting)
        - gaia_queried: int, number of Gaia sources returned from query
        - gaia_used: int, number of Gaia sources used for WCS (after limiting)
        - matches: int, number of successful star matches
    """
    # Initialize return dictionary
    result = {
        'success': False,
        'error': None,
        'crpix': None,
        'crval': None,
        'sources_detected': 0,
        'sources_used': 0,
        'gaia_queried': 0,
        'gaia_used': 0,
        'matches': 0
    }

    if verbose:
        print(f"Starting pointer_wcs for file: {filepath}")

    try:
        with fits.open(filepath, mode='update') as hdu:
            header = hdu[0].header
            data = hdu[0].data.astype(float)

            if verbose:
                print(f"FITS file opened successfully")
                print(f"Data shape: {data.shape}")
                print(f"IMAGETYP: {header.get('IMAGETYP', 'NOT_FOUND')}")

            # Clear existing WCS headers if requested
            if clear_existing_wcs and wcs_keywords:
                if verbose:
                    print("Clearing any existing WCS headers...")
                removed = clear_wcs_headers(header, wcs_keywords, verbose=verbose)
                if verbose and removed > 0:
                    print(f"Removed {removed} existing WCS headers")
                elif verbose:
                    print("No existing WCS headers found to remove")

            if header['IMAGETYP'] == 'Light Frame':
                if verbose:
                    print(f"Confirmed Light Frame, starting image cleaning...")

                # Clean image (no dark frame subtraction - images are pre-reduced)
                image_clean = clean_image(data)

                if verbose:
                    print(f"Image cleaning complete, cleaned image shape: {image_clean.shape}")

                # Get RA and DEC from header
                ra = header['RA']
                dec = header['DEC']

                if verbose:
                    print(f"Header RA: {ra}, DEC: {dec}")

                # Handle coordinate format based on telescope
                if header['TELESCOP'] == 'Artemis':
                    center = SkyCoord(ra, dec, unit=[u.hourangle, u.deg])
                else:
                    center = SkyCoord(ra, dec, unit=[u.deg, u.deg])

                if verbose:
                    print("COORDINATES AFTER SKYCOORD")
                    print(center.ra.deg)
                    print(center.dec.deg)
                    print(f"Center coordinates: RA={center.ra.deg:.6f} deg, DEC={center.dec.deg:.6f} deg")

                # Calculate field of view
                shape = image_clean.shape
                plate_scale = np.arctan((header['XPIXSZ'] * 1e-6) / (header['FOCALLEN'] * 1e-3)) * (180 / np.pi)
                fovx = (1 / np.abs(np.cos(center.dec.rad))) * shape[0] * plate_scale
                fovy = shape[1] * plate_scale
                fov = np.array([fovx, fovy])

                if verbose:
                    print(f"Plate scale: {plate_scale * 3600:.3f} arcsec/pixel")
                    print(f"FOV: X={fovx * 60:.2f} arcmin, Y={fovy * 60:.2f} arcmin")

                # Extract observation date for proper motion correction
                dateobs = None
                if 'DATE-OBS' in header:
                    dateobs = pd.to_datetime(header["DATE-OBS"])
                    if verbose:
                        print(f"Using observation date: {dateobs}")

                # Detect stars using multiscale detection
                if verbose:
                    print(f"Starting multiscale star detection...")

                stars_in_image = find_stars_multiscale(
                    image_clean, scales=[1, 2, 3], threshold=7, edge_buffer=10
                )
                result['sources_detected'] = len(stars_in_image)

                if verbose:
                    print(f"Detected {len(stars_in_image)} stars in the image.")
                    print(f"Detected {len(stars_in_image)} stars in the image after multiscale detection.")

                if len(stars_in_image) < 4:
                    result['error'] = "Not enough stars detected for plate solve"
                    return result

                # Limit number of stars to prevent memory issues
                star_limit = min(16, len(stars_in_image))
                stars_in_image = stars_in_image[0:star_limit]
                result['sources_used'] = len(stars_in_image)

                if verbose:
                    print(f"Using {len(stars_in_image)} stars (limited from {result['sources_detected']})")

                # Query local Gaia database
                if verbose:
                    print(f"Querying local Gaia database...")

                gaia_stars = gaia_db_query(
                    center=(center.ra.deg, center.dec.deg),
                    fov=1.2 * fov,  # 20% margin
                    tmass=True,
                    dateobs=dateobs,
                    limit=2 * star_limit,
                    db_path=db_path
                )
                result['gaia_queried'] = len(gaia_stars)

                if verbose:
                    print(f"Found {len(gaia_stars)} Gaia stars in the field of view. with limit={2 * star_limit}")

                if len(gaia_stars) < 4:
                    result['error'] = "Not enough Gaia stars found in field"
                    return result

                # Limit Gaia stars
                gaia_limit = min(2 * star_limit, len(gaia_stars))
                gaia_stars = gaia_stars[0:gaia_limit]
                result['gaia_used'] = len(gaia_stars)

                if verbose:
                    print(f"Using {len(gaia_stars)} Gaia stars (limited from {result['gaia_queried']})")
                    print(
                        f"Using {len(stars_in_image)} image stars and {len(gaia_stars)} Gaia stars for WCS solution...")

                # Create image star mapping and compute WCS
                if verbose:
                    print(f"Starting WCS computation...")

                image_star_mapping = ImageStarMapping.from_gaia_coordinates(
                    stars_in_image, gaia_stars
                )

                if verbose:
                    print(f"WCS computation complete")

                # Validate star matching
                matches = image_star_mapping.number_of_matched_stars(pixel_threshold=10)
                result['matches'] = matches

                if verbose:
                    print(f"Number of matched stars: {matches}")

                min_matches = 4
                if matches < min_matches:
                    result['error'] = f"Plate solve failed, not enough stars matched ({matches}/{min_matches})"
                    return result

                # Get plating center
                plating_ra, plating_dec = image_star_mapping.get_plating_center(
                    image_shape=image_clean.shape
                )

                # Update header with WCS information
                if verbose:
                    print(f"Updating FITS header with WCS information...")

                wcs_header = image_star_mapping.wcs.to_header()
                if verbose:
                    print(f"Generated WCS header with {len(wcs_header)} keywords")
                    for key in wcs_header:
                        print(f"  {key}: {wcs_header[key]}")

                header.update(wcs_header)

                # Extract WCS values for return
                result['crpix'] = (header['CRPIX1'], header['CRPIX2'])
                result['crval'] = (header['CRVAL1'], header['CRVAL2'])
                result['success'] = True

                # Create diagnostic plots if verbose
                if verbose:
                    print(f"Creating diagnostic plots...")
                    _create_diagnostic_plots(filepath, data, image_clean, stars_in_image,
                                             image_star_mapping.gaia_stars_in_image, image_star_mapping.wcs)

                if verbose:
                    print(f"WCS solution completed successfully!")

            else:
                result['error'] = f"Not a Light Frame (IMAGETYP={header.get('IMAGETYP')})"
                if verbose:
                    print(f"Skipping - not a Light Frame (IMAGETYP={header.get('IMAGETYP')})")

    except Exception as e:
        result['error'] = str(e)
        if verbose:
            print(f"ERROR in pointer_wcs: {e}")

    return result


def print_wcs_headers(filepath):
    """
    Print WCS header information for a FITS file after plate solving.

    Parameters
    ----------
    filepath : str
        Path to the FITS file
    """
    print(f"Processing: {filepath}")

    # Perform plate solving
    result = pointer_wcs(filepath, verbose=True)

    if result['success']:
        print("\n" + "=" * 50)
        print("WCS SOLUTION SUCCESSFUL")
        print("=" * 50)

        # Read the updated header
        with fits.open(filepath) as hdu:
            header = hdu[0].header

        # Print key WCS parameters
        wcs_keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                        'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

        print("\nWCS Header Keywords:")
        print("-" * 30)
        for keyword in wcs_keywords:
            if keyword in header:
                print(f"{keyword:<8} = {header[keyword]}")

        print(f"\nDetection Summary:")
        print(f"  Sources detected: {result['sources_detected']}")
        print(f"  Sources used: {result['sources_used']}")
        print(f"  Gaia stars queried: {result['gaia_queried']}")
        print(f"  Gaia stars used: {result['gaia_used']}")
        print(f"  Star matches: {result['matches']}")

    else:
        print("\n" + "=" * 50)
        print("WCS SOLUTION FAILED")
        print("=" * 50)
        print(f"Error: {result['error']}")
        print(f"Sources detected: {result['sources_detected']}")
        print(f"Gaia stars queried: {result['gaia_queried']}")


if __name__ == "__main__":
    # Example usage
    import sys

    if len(sys.argv) > 1:
        print_wcs_headers(sys.argv[1])
    else:
        print("Usage: python pointer_wcs.py <fits_file>")