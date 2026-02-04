# This module makes a proper WCS-header of a FITS file.

from astropy.io import fits
import os
from functools import partial
import argparse
import timeit
import astrom.find_pa
from calibration.pipeutils import detect_instrument, get_instrument_parameters, open_fits_file
import fnmatch
import shutil
from astrom.twirl_speculoos import twirl_wcs
from astrom.pointer_wcs import pointer_wcs
def has_valid_wcs(header):
    """
    Check if header contains valid WCS data, not just WCS keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
        FITS header to check

    Returns
    -------
    bool
        True if header contains valid WCS data, False otherwise
    """
    # Check if CTYPE1 exists and has valid value
    ctype1 = header.get('CTYPE1', '')
    if ctype1 in ['nan', '', None] or str(ctype1).lower() == 'nan':
        return False

    # Check if reference values are reasonable (not zero/nan)
    crval1 = header.get('CRVAL1', 0.0)
    crval2 = header.get('CRVAL2', 0.0)
    if crval1 == 0.0 or crval2 == 0.0:
        return False

    # Check if reference pixels are reasonable (not zero)
    crpix1 = header.get('CRPIX1', 0.0)
    crpix2 = header.get('CRPIX2', 0.0)
    if crpix1 == 0.0 or crpix2 == 0.0:
        return False

    return True


def add_astrometry_with_timeout(infile, timeout, ext, db_path, raw_images, trim_offsets,
                                total_files, file_num=None):
    """
    Wrapper that runs add_astrometry with a signal-based timeout.

    Parameters
    ----------
    infile : str
        Path to image file
    timeout : int
        Timeout in seconds
    ext : str
        File extension
    db_path : str
        Path to Gaia database
    raw_images : list
        List of raw image paths
    trim_offsets : tuple
        Trim offsets for CRPIX adjustment
    total_files : int
        Total number of files being processed
    file_num : int, optional
        File number for progress display

    Returns
    -------
    dict
        Result dictionary from add_astrometry or timeout error
    """
    import signal

    class TimeoutException(Exception):
        pass

    def timeout_handler(signum, frame):
        raise TimeoutException(f"Plate solving timed out after {timeout}s")

    # Set up the timeout
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)

    try:
        result = add_astrometry(infile, ext, db_path, raw_images, trim_offsets,
                                file_num=file_num, total_files=total_files)
        signal.alarm(0)  # Cancel the alarm
        return result
    except TimeoutException as e:
        signal.alarm(0)  # Cancel the alarm
        print(f"[TIMEOUT] {infile}: {str(e)}")
        return {
            'success': False,
            'error': str(e),
            'crpix': None,
            'crval': None,
            'sources_detected': 0,
            'sources_used': 0,
            'gaia_queried': 0,
            'gaia_used': 0,
            'matches': 0
        }
    except Exception as e:
        signal.alarm(0)  # Cancel the alarm
        print(f"[ERROR] {infile}: {str(e)}")
        return {
            'success': False,
            'error': str(e),
            'crpix': None,
            'crval': None,
            'sources_detected': 0,
            'sources_used': 0,
            'gaia_queried': 0,
            'gaia_used': 0,
            'matches': 0
        }

def main(args):
    import sys
    sys.stdout.flush()
    print("*************************", flush=True)
    print("Starting plate solving", flush=True)
    print("*************************", flush=True)
    start_time = timeit.default_timer()
    print(f"Input file: {args.filelist}", flush=True)
    print(f"Extension: {args.ext}", flush=True)
    print(f"Force plate solve: {args.force_platesolve}", flush=True)

    # Before the try block
    infiles = []
    print("About to process filelist...", flush=True)

    total_files = 0
    already_solved = 0
    newly_solved = 0
    failed_solve = 0

    try:
        with open(args.filelist) as filel:
            # print("File opened successfully", flush=True)
            for line in filel:
                f = line.strip()
                total_files += 1
                # print(f"Found file in list: {f}", flush=True)
                # print(f)
                hdulist = fits.open(f)
                # print(f"FITS file opened: {f}", flush=True)

                has_wcs = has_valid_wcs(hdulist[0].header)

                if args.force_platesolve:
                    # Force mode: process all Light Frames
                    infiles.append(f)
                    # print(f"Added to infiles: {f}", flush=True)
                    if has_wcs:
                        already_solved += 1
                else:
                    # Normal mode: only process files without WCS
                    if has_wcs:
                        already_solved += 1
                    else:
                        infiles.append(f)
                        # print(f"Added to infiles: {f}", flush=True)

                hdulist.close()  # Add this to close the file properly
    except Exception as e:
        print("Importing Astrometry failed: " + str(e), flush=True)

    print(f"Total files: {total_files}", flush=True)
    print(f"Files to process: {len(infiles)}", flush=True)
    if args.force_platesolve:
        print(f"Files with existing WCS (will be overwritten): {already_solved}", flush=True)
    else:
        print(f"Files already solved (skipped): {already_solved}", flush=True)

    # Read raw image list if provided
    raw_images = []
    if args.raw_list and os.path.exists(args.raw_list):
        with open(args.raw_list) as f:
            raw_images = [line.strip() for line in f if line.strip()]

    # Get parameters once from the first file in the list
    with open(args.filelist) as filel:
        first_filename = filel.readline().strip()

    with open_fits_file(first_filename) as hdul:
        params = get_instrument_parameters(hdul, trimmed=True)  # <-- Add trimmed=True

    # Extract trim offsets once
    trim_offsets = (0, 0)  # default
    if 'trim' in params:
        trim_config = params['trim']
        trim_offsets = (trim_config.get('left_col', 0), trim_config.get('top_row', 0))

    files_to_process = len(infiles)

    if int(args.nproc) > 1:
        # Parallel processing with timeout
        from multiprocessing import Pool
        from functools import partial

        print(f"Running in parallel mode with {args.nproc} processes", flush=True)
        pool = Pool(processes=int(args.nproc))

        worker_func = partial(add_astrometry_with_timeout,
                              timeout=args.timeout,
                              ext=args.ext,
                              db_path=args.db_path,
                              raw_images=raw_images,
                              trim_offsets=trim_offsets,
                              total_files=files_to_process)

        results = pool.map(worker_func, infiles)
        pool.close()
        pool.join()

        # Count successes and failures
        for result in results:
            if result and result.get('success', False):
                newly_solved += 1
            else:
                failed_solve += 1
    else:
        # Sequential processing with timeout
        print(f"Running in serial mode (nproc={args.nproc})", flush=True)
        for i, infile in enumerate(infiles, 1):
            result = add_astrometry_with_timeout(infile,
                                                 timeout=args.timeout,
                                                 ext=args.ext,
                                                 db_path=args.db_path,
                                                 raw_images=raw_images,
                                                 trim_offsets=trim_offsets,
                                                 total_files=files_to_process,
                                                 file_num=i)

            if result and result.get('success', False):
                newly_solved += 1
            else:
                failed_solve += 1

    try:
        for f in infiles:
            new_f = f.replace(args.ext, 'new')
            f1 = f.replace(args.ext, 'axy')
            f2 = f.replace(args.ext, 'wcs')
            f3 = f.replace(args.ext, 'solved')
            f4 = f.replace(args.ext, 'corr')
            f5 = f.replace("." + args.ext, '-indx.xyls')
            f6 = f.replace(args.ext, 'match')
            f7 = f.replace(args.ext, 'rdls')
            if os.path.isfile(new_f):

                hdulist_old = fits.open(f)
                hdr_old = hdulist_old[0].header
                data = hdulist_old[0].data
                hdulist_new = fits.open(new_f)
                hdr_new = hdulist_new[0].header
                pa, jd, solve = find_pa.pa(new_f)

                hdr_old.set('PA', pa)
                # hdr_old.set('CTYPE1',hdr_new['CTYPE1']) #'RA---TAN')
                # hdr_old.set('CTYPE2', hdr_new['CTYPE2']) #'DEC--TAN')
                # hdr_old.set('CRVAL1',hdr_new['CRVAL1'])
                # hdr_old.set('CRVAL2',hdr_new['CRVAL2'])
                # hdr_old.set('CRPIX1',hdr_new['CRPIX1'])
                # hdr_old.set('CRPIX2',hdr_new['CRPIX2'])
                # hdr_old.set('CD1_1',hdr_new['CD1_1'])
                # hdr_old.set('CD2_1',hdr_new['CD2_1'])
                # hdr_old.set('CD1_2',hdr_new['CD1_2'])
                # hdr_old.set('CD2_2',hdr_new['CD2_2'])
                for k in WCS_KEYWORDS:
                    if hdr_old.get(k) == None:
                        hdr_old.set(k, hdr_new[k])
                # if hdr_old.get('RA') == None:
                #     hdr_old.set('RA', '21 31 48.113')
                #     hdr_old.set('DEC', '48 26 29.330')
                hdr_old.set('ASTSOLVE', 'T')
                fits.writeto(f, data, hdr_old, overwrite=True)

                for fname in [new_f, f1, f2, f3, f4, f5, f6, f7]:
                    os.remove(fname)

        elapsed = timeit.default_timer() - start_time
        print('Total time taken for astrom: ' + str(elapsed / 60.) + ' minutes')

    except Exception as e:
        print("Replacing Astrometry failed: " + str(e))
        print("on file: " + f)

    # Print summary
    print("\n" + "=" * 50)
    print("\n" + "=" * 50)
    print("PLATE SOLVING SUMMARY")
    print("=" * 50)
    print(f"Total images in list: {total_files}")
    if args.force_platesolve:
        print(f"Images with existing WCS (overwritten): {already_solved}")
        print(f"Images without existing WCS: {files_to_process - already_solved}")
    else:
        print(f"Images already plate solved (skipped): {already_solved}")
    print(f"Images processed: {files_to_process}")
    print(f"Successfully plate solved: {newly_solved}")
    print(f"Failed plate solving: {failed_solve}")
    if files_to_process > 0:
        print(f"Success rate: {newly_solved / files_to_process * 100:.1f}%")
    print("=" * 50)


def find_raw_image(processed_path, raw_images):
    """
    Find the corresponding raw image for a processed image.

    Parameters
    ----------
    processed_path : str
        Path to the processed image (e.g., 'proc_IMAGE123.fits')
    raw_images : list
        List of raw image paths

    Returns
    -------
    str or None
        Path to corresponding raw image, or None if not found
    """
    processed_basename = os.path.basename(processed_path)

    if processed_basename.startswith('proc'):
        raw_basename_with_ext = processed_basename[4:]  # Remove 'proc' prefix
    else:
        raw_basename_with_ext = processed_basename

    # Remove extension for matching
    raw_basename_no_ext = os.path.splitext(raw_basename_with_ext)[0]

    # Look for matching raw image (ignoring extension)
    for raw_path in raw_images:
        raw_file_basename = os.path.basename(raw_path)
        raw_file_no_ext = os.path.splitext(raw_file_basename)[0]

        if raw_file_no_ext == raw_basename_no_ext:
            return raw_path

    return None

def get_trim_offsets(instrument):
    """
    Get pixel offsets for trimming based on instrument type.

    Parameters
    ----------
    instrument : str
        Instrument name ('spirit', 'andor', 'moana')

    Returns
    -------
    tuple
        (x_offset, y_offset) - pixels removed from left and top edges
    """
    if instrument == 'spirit':
        return (0, 0)  # No trimming
    elif instrument == 'andor':
        return (2, 22)  # Left=2, Top=22 based on [22:2066, 2:2048]
    elif instrument == 'moana':
        return (1, 1)  # Based on [1:-1, 1:-1]
    else:
        return (0, 0)  # Unknown instrument, no adjustment


# Define WCS keywords once at module level
WCS_KEYWORDS = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CUNIT1', 'CUNIT2',
                'WCSAXES', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1',
                'CDELT2', 'LONPOLE', 'LATPOLE', 'RADESYS']


def copy_wcs_to_raw(processed_path, raw_path, trim_offsets):
    """
    Copy WCS headers from processed image to raw image.

    Parameters
    ----------
    processed_path : str
        Path to processed image with WCS solution
    raw_path : str
        Path to raw image to update
    """
    from calibration.pipeutils import get_instrument_parameters

    # Use module-level WCS keywords
    wcs_keywords = WCS_KEYWORDS

    # Read WCS headers from processed image
    with fits.open(processed_path) as processed_hdu:
        processed_header = processed_hdu[0].header

        # Extract WCS keywords
        wcs_headers = {}
        for keyword in wcs_keywords:
            if keyword in processed_header:
                wcs_headers[keyword] = processed_header[keyword]

    x_offset, y_offset = trim_offsets

    # Adjust CRPIX values for trimming
    if 'CRPIX1' in wcs_headers:
        wcs_headers['CRPIX1'] += x_offset
    if 'CRPIX2' in wcs_headers:
        wcs_headers['CRPIX2'] += y_offset

    # Update raw image with WCS headers
    with fits.open(raw_path, mode='update') as raw_hdu:
        raw_header = raw_hdu[0].header

        # Add WCS keywords to raw image
        for keyword, value in wcs_headers.items():
            raw_header[keyword] = value

        # # Log the adjustment for debugging
        # if x_offset != 0 or y_offset != 0:
        #     print(f"Debug: Adjusted CRPIX by ({x_offset}, {y_offset}) using config-based trim offsets")

def add_astrometry(f, ext, db_path, raw_images, trim_offsets, file_num=None, total_files=None):
    if fnmatch.fnmatch(f, '*.' + ext):
        file = f
        # print("Add astrometry.net data to header of file: " + file)

        # Attempt WCS solving
        # result = twirl_wcs(str(file), verbose=False)
        result = pointer_wcs(str(file), db_path, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True, verbose=False)

        # If successful and raw images provided, update raw image with WCS
        if result['success'] and raw_images:
            raw_file = find_raw_image(file, raw_images)
            if raw_file:
                try:
                    copy_wcs_to_raw(file, raw_file, trim_offsets)
                except Exception as e:
                    print(f"Warning: Failed to update raw image {raw_file}: {e}")

        # Create summary line
        filename = os.path.basename(file)
        if file_num is not None and total_files is not None:
            prefix = f"[{file_num}/{total_files}]"
        else:
            prefix = ""

        if result['success']:
            crpix1, crpix2 = result['crpix']
            crval1, crval2 = result['crval']
            summary = (f"{prefix} {filename} SUCCESS - "
                       f"Sources: {result['sources_detected']}→{result['sources_used']} "
                       f"Gaia: {result['gaia_queried']}→{result['gaia_used']} "
                       f"Matches: {result['matches']} "
                       f"CRPIX=({crpix1:.2f},{crpix2:.2f}) "
                       f"CRVAL=({crval1:.6f},{crval2:.6f})")
        else:
            summary = (f"{prefix} {filename} FAILED - "
                       f"Sources: {result['sources_detected']}→{result['sources_used']} "
                       f"Gaia: {result['gaia_queried']}→{result['gaia_used']} "
                       f"Matches: {result['matches']} "
                       f"Error: {result['error']}")

        print(summary)
        return result

    return None  # Add this for non-matching files

        #
        # hdulist = fits.open(file)
        # prihdr = hdulist[0].header
        #
        # # Extract RA and DEC
        # ra_val = prihdr.get('RA')
        # dec_val = prihdr.get('DEC')
        #
        # if ra_val is not None and dec_val is not None:
        #     # Process RA
        #     if " " in str(ra_val):
        #         strra = str(ra_val.replace(" ", ":"))
        #     else:
        #         strra = str(ra_val)
        #
        #     # Process DEC
        #     if " " in str(dec_val):
        #         strdec = str(dec_val.replace(" ", ":"))
        #     else:
        #         strdec = str(dec_val)
        #
        #     if detect_instrument(hdulist) == 'spirit':
        #         print("Running solve-field with spirit astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/spirit_astrometry.cfg " \
        #                      f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                      f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                      f"--scale-low 0.32 --scale-high 0.38 --scale-units arcsecperpix {file}"
        #     elif detect_instrument(hdulist) == 'andor':
        #         print("Running solve-field with andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"
        #     else:
        #         print("Instrument not spirit or andor, defaulting to andor astrometry.cfg file.")
        #         cmd = f"{solve_field_cmd} --config /opt/orchard/src/astrom/andor_astrometry.cfg " \
        #                     f"--no-plots --no-remove-lines --uniformize 0 --overwrite " \
        #                     f"--ra {strra} --dec {strdec} --radius 4 --downsample 2 " \
        #                     f"--scale-low 0.60 --scale-high 0.68 --scale-units arcsecperpix {file}"

        # print(f"Running: {cmd}")
        # os.system(cmd)
        # else:
        #     print(f"Warning: Missing RA and/or DEC in {file}, skipping astrom")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filelist')
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help='Number of parallel processes (default: 1)')
    parser.add_argument('-e', '--ext')
    parser.add_argument('-d', '--db_path')
    parser.add_argument('-r', '--raw-list', help='Path to raw image list file')
    parser.add_argument('-f', '--force-platesolve', action='store_true',
                        help='Force plate solving even for images that already have WCS')
    parser.add_argument('-t', '--timeout', type=int, default=60,
                        help='Timeout per image in seconds (default: 60)')
    main(parser.parse_args())