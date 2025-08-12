import numpy as np
import sys
import os
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import glob

def solve_astrometry(filenameold, frame_type):
    # images have to be fits not fts (otherwise astrometry solving doesn't work)
    filename = filenameold.replace('.fts', '.fits')
    filename = filename.replace('.FIT', '.fits')
    os.rename(filenameold, filename)
    filenameb = filename
    with fits.open(filename) as infile_init:
        if infile_init[0].header['IMAGETYP'] == 'Light Frame':
            if infile_init[0].header['FILTER'] != 'zYJ':
                if len(filename.split('_')) == 1:  # Only need to remove ' characters from filenames of early images
                    print(filenameold)
                    print(filename)
                    if infile_init[0].header['FILTER'][
                        0] == 'u':  # remove ' character from image name (otherwise astrometry solving doesn't work)
                        filenameb = filename[:-6] + '.fits'
                        os.rename(filename, filenameb)
                    if infile_init[0].header['FILTER'][0] == 'g':
                        filenameb = filename[:-6] + '.fits'
                        os.rename(filename, filenameb)
                    if infile_init[0].header['FILTER'][0] == 'r':
                        filenameb = filename[:-6] + '.fits'
                        os.rename(filename, filenameb)
                    if infile_init[0].header['FILTER'][0] == 'i':
                        filenameb = filename[:-6] + '.fits'
                        os.rename(filename, filenameb)
                    if infile_init[0].header['FILTER'][0] == 'z':
                        filenameb = filename[:-6] + '.fits'
                        os.rename(filename, filenameb)
                    print(filenameb)

                # Doing 3 times try/except to mitigate the fact that the current angle parsing code is not thread safe
                try:
                    c = SkyCoord(ra=infile_init[0].header['RA'], dec=infile_init[0].header['DEC'], unit=(u.deg, u.deg),
                                 frame='fk5')
                except Exception as e1:
                    print(e1)
                    try:
                        c = SkyCoord(ra=infile_init[0].header['RA'], dec=infile_init[0].header['DEC'],
                                     unit=(u.deg, u.deg), frame='fk5')
                    except Exception as e2:
                        print(e2)
                        try:
                            c = SkyCoord(ra=infile_init[0].header['RA'], dec=infile_init[0].header['DEC'],
                                         unit=(u.deg, u.deg), frame='fk5')
                        except Exception as e3:
                            print(e3)

                rep_chars = ['h', 'd', 'm']
                str1 = c.ra.to_string(unit=u.hourangle)[:-1]
                str2 = c.dec.to_string()[:-1]
                for char in rep_chars:
                    str1 = str1.replace(char, ':')
                    str2 = str2.replace(char, ':')

                # RA = infile_init[0].header['RA']
                # str1 = str(RA.replace(" ", ":"))
                # DEC = infile_init[0].header['DEC']
                # str2 = str(DEC.replace(" ", ":"))
                if frame_type == 'Light Frame':
                    os.system(
                        f"/appct/data/SPECULOOSPipeline/astrometry.net-0.73/blind/solve-field --no-plots --no-remove-lines --uniformize 0 --overwrite --ra {str1} --dec {str2} --radius 4 --depth 100 --downsample 2 --no-tweak {filenameb}")
                filenamebis = filenameb.replace('.fits', '.new')
                dum = filenameb.replace('.fits', '-indx.xyls')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.axy')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.corr')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.match')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.rdls')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.solved')
                if os.path.exists(dum):
                    os.remove(dum)
                dum = filenameb.replace('.fits', '.wcs')
                if os.path.exists(dum):
                    os.remove(dum)
                if os.path.exists(filenamebis):
                    if os.path.exists(filenameb):
                        os.remove(filenameb)
                    os.rename(filenamebis, filenameb)

def transformation_check(imgdir, down_files, type):
    transformed_files = []
    for down_file in down_files:
        if "SPECU" in down_file and ".Z" not in down_file:
            hdul = fits.open(f"{imgdir}/{down_file}.fits")
            if type == "ANDOR":
                if "ASTRAROT" in hdul[0].header and len(down_file.split('.')[-1]) == 3:
                    print(f"{down_file} is an ANDOR astra image and has already been rotated. No need to redownload.")
                    transformed_files.append(down_file)
                else:
                    print(f"{down_file} is an ANDOR astra image, but might not be rotated. Delete and redownload.")
                    os.remove(f"{imgdir}/{down_file}.fits")
            elif type == "SPIRIT":
                if "ASTRAMIR" in hdul[0].header and len(down_file.split('.')[-1]) == 3:
                    print(f"{down_file} is a SPIRIT astra image and has already been flipped. No need to redownload.")
                    transformed_files.append(down_file)
                else:
                    print(f"{down_file} is a SPIRIT astra image, but might not be flipped. Delete and redownload.")
                    os.remove(f"{imgdir}/{down_file}.fits")

    return transformed_files


def rotate_wcs_acw(header, n_rotations=1):
    """
    Update WCS headers for 90 deg  anticlockwise IMAGE rotation.
    Note: WCS transforms in opposite direction to compensate.
    """

    nx = header.get('NAXIS1', 0) * 1
    ny = header.get('NAXIS2', 0) * 1

    for rot in range(n_rotations):
        # Update CRPIX for 90 deg CCW image rotation
        if 'CRPIX1' in header and 'CRPIX2' in header:
            old_crpix1 = header['CRPIX1']
            old_crpix2 = header['CRPIX2']

            # if rot % 2 == 1:
            #     # Formula for 90 deg CCW image rotation (WCS transforms opposite)
            #     header['CRPIX1'] = ny - old_crpix2 + 1
            #     header['CRPIX2'] = old_crpix1
            # else:
            # Formula for 90 deg CCW image rotation (WCS transforms opposite)
            header['CRPIX1'] = old_crpix2
            header['CRPIX2'] = nx - old_crpix1 + 1

        # Update CD matrix for 90 deg CCW image rotation (WCS rotates CW to compensate)
        cd_keys = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
        if all(key in header for key in cd_keys):
            old_cd11 = header['CD1_1']
            old_cd12 = header['CD1_2']
            old_cd21 = header['CD2_1']
            old_cd22 = header['CD2_2']

            # 90 deg CW rotation of CD matrix (opposite to image rotation)
            header['CD1_1'] = -old_cd21
            header['CD1_2'] = -old_cd22
            header['CD2_1'] = old_cd11
            header['CD2_2'] = old_cd12

        # Alternative: Update CDELT/CROTA if present instead of CD matrix
        elif 'CDELT1' in header and 'CDELT2' in header and 'CROTA2' in header:
            header['CROTA2'] = header['CROTA2'] - 90  # Subtract 90 deg (opposite to image)

        # Update non-standard dimension keywords if present
        if 'IMAGEW' in header and 'IMAGEH' in header:
            old_imagew = header['IMAGEW']
            old_imageh = header['IMAGEH']
            header['IMAGEW'] = old_imageh
            header['IMAGEH'] = old_imagew

        tmp = ny * 1
        ny = nx * 1
        nx = tmp * 1


def mirror_wcs_horizontal(header):
    """
    Update WCS headers for horizontal mirroring (flip left-right).
    This flips the image around the vertical axis (y-axis).
    """
    nx = header.get('NAXIS1', 0)

    # Update CRPIX for horizontal flip
    if 'CRPIX1' in header:
        header['CRPIX1'] = nx - header['CRPIX1'] + 1

    # Update CD matrix for horizontal flip (negate first column)
    cd_keys = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    if all(key in header for key in cd_keys):
        header['CD1_1'] = -header['CD1_1']
        header['CD2_1'] = -header['CD2_1']

    # Alternative: Update CDELT if present
    elif 'CDELT1' in header:
        header['CDELT1'] = -header['CDELT1']


def get_output_filename(filename, dirname, n_repeats, test=False):
    """
    Generate output filename, optionally with '_test' suffix.
    """
    if test:
        name, ext = os.path.splitext(filename)
        return os.path.join(dirname, f"{name}_test{n_repeats}{ext}")
    else:
        return os.path.join(dirname, filename)


def astra_transform_unified(input_data, transform_mode='auto', n_repeats=1, completeness="all",
                            test=False, dirname=None):
    """
    Unified function to transform FITS files with proper WCS handling.

    Parameters:
    -----------
    input_data : str, list
        Can be:
        - Single directory path (string)
        - List of directory paths (for multiple dates)
        - List of filenames (without .fits extension) - requires dirname parameter
    transform_mode : str, optional
        Transformation mode:
        - 'auto': Auto-detect based on NAXIS1 (SPIRIT=mirror, ANDOR=rotate)
        - 'rotate': Force 90 deg CCW rotation
        - 'mirror': Force horizontal mirroring
        - 'rotate_wcs_only': Rotate WCS only (no image transformation)
        - 'rotate_image_only': Rotate image only (no WCS transformation)
    n_repeats : int, optional
        The number of times to repeat the transformation (for ANDOR images)
    test : bool, optional
        If True, write files with '_test' suffix to avoid overwriting originals
    dirname : str, optional
        Directory path when input_data is a list of filenames
    """

    # Determine input type and build file list
    files_to_process = []

    if isinstance(input_data, str):
        # Single directory path
        directory = input_data
        all_fits_files = glob.glob(os.path.join(directory, "*.fits"))
        fits_files = [f for f in all_fits_files if not f.endswith('test.fits')]

        for filepath in fits_files:
            files_to_process.append({
                'filepath': filepath,
                'filename': os.path.basename(filepath),
                'name_without_ext': os.path.splitext(os.path.basename(filepath))[0],
                'directory': directory
            })

    elif isinstance(input_data, list):
        # Check if it's a list of directories or filenames
        if len(input_data) > 0 and os.path.isdir(input_data[0]):
            # List of directories
            for directory in input_data:
                if not os.path.isdir(directory):
                    print(f"Warning: Directory does not exist: {directory}")
                    continue

                all_fits_files = glob.glob(os.path.join(directory, "*.fits"))
                fits_files = [f for f in all_fits_files if not f.endswith('test.fits')]

                for filepath in fits_files:
                    files_to_process.append({
                        'filepath': filepath,
                        'filename': os.path.basename(filepath),
                        'name_without_ext': os.path.splitext(os.path.basename(filepath))[0],
                        'directory': directory
                    })
        else:
            # List of filenames - requires dirname
            if dirname is None:
                print("Error: dirname parameter required when input_data is a list of filenames")
                return

            # Check observer for backward compatibility
            if len(input_data) > 0:
                first_file = os.path.join(dirname, input_data[0] + ".fits")
                if os.path.exists(first_file):
                    hdul = fits.open(first_file)
                    observer = hdul[0].header.get('OBSERVER', '')
                    hdul.close()

                    if observer != 'Astra' and transform_mode == 'auto':
                        print("Images for this night created by ACP, so no need to rotate them.")
                        return

            for filename in input_data:
                filepath = os.path.join(dirname, filename + ".fits")
                if os.path.exists(filepath):
                    files_to_process.append({
                        'filepath': filepath,
                        'filename': filename + ".fits",
                        'name_without_ext': filename,
                        'directory': dirname
                    })
    else:
        print("Error: input_data must be a string (directory) or list (directories or filenames)")
        return

    if not files_to_process:
        print("No FITS files found to process")
        return

    if completeness == 'single':
        files_to_process = [files_to_process[0]]

    print(f"Processing {len(files_to_process)} FITS files...")

    # Process each file
    for file_info in files_to_process:
        filepath = file_info['filepath']
        filename = file_info['filename']
        name_without_ext = file_info['name_without_ext']
        directory = file_info['directory']

        output_path = get_output_filename(filename, directory, n_repeats, test)

        hdul = fits.open(filepath)

        # Make copies of header and data
        header = hdul[0].header.copy()
        data = hdul[0].data.copy()

        # Close the original file
        hdul.close()

        # Get frame information
        frame_type = header.get('IMAGETYP', '')
        naxis1 = header.get('NAXIS1', 0)
        is_light_frame = (frame_type == 'Light Frame')

        # Determine transformation type
        if transform_mode == 'auto':
            # Auto-detect based on NAXIS1 (backward compatibility)
            is_spirit = (naxis1 == 1024)
            actual_transform = 'mirror' if is_spirit else 'rotate'
        else:
            actual_transform = transform_mode.lower()

        # Apply transformation
        if actual_transform == 'mirror':
            # Horizontal mirroring
            data = np.flip(data, axis=1)  # flip along axis 1 (horizontal)
            header['ASTRAMIR'] = ('Y', 'Y axis flip applied to astra images')

            # Update WCS if this is a Light Frame
            if is_light_frame:
                mirror_wcs_horizontal(header)
                header.add_history('WCS updated for horizontal mirroring')
                print(f"Mirrored Light Frame {name_without_ext} with WCS updated.")
            else:
                print(f"Mirrored calibration image {name_without_ext}.")

        elif actual_transform == 'rotate':
            # 90 deg  anticlockwise rotation
            data = np.rot90(data, k=n_repeats)  # k=1 means 90 deg CCW
            header['ASTRAROT'] = (90*n_repeats, 'Anticlockwise rotation applied to astra images')

            # Update WCS if this is a Light Frame
            if is_light_frame:
                rotate_wcs_acw(header, n_rotations=n_repeats)
                header.add_history(f'WCS updated for {n_repeats*90} deg  anticlockwise rotation')
                print(f"Rotated Light Frame {name_without_ext} {n_repeats*90} deg with WCS updated.")
            else:
                print(f"Rotated calibration image {name_without_ext} {n_repeats*90} deg.")

        elif actual_transform == 'rotate_wcs_only':
            # WCS-only rotation (no image transformation)
            if is_light_frame:
                rotate_wcs_acw(header, n_rotations=n_repeats)
                header['ASTRAROT'] = (90*n_repeats, f'WCS rotated {n_repeats*90} deg ACW (image unchanged)')
                header.add_history(f'WCS rotated {n_repeats*90} deg  anticlockwise (image data unchanged)')
                print(f"Rotated WCS only for Light Frame {name_without_ext}.")
            else:
                print(f"Skipped calibration image {name_without_ext} (no WCS to rotate).")
                continue

        elif actual_transform == 'rotate_image_only':
            # 90 deg anticlockwise rotation
            data = np.rot90(data, k=n_repeats)  # k=1 means 90 deg CCW
            print(f"Rotated image {name_without_ext} {n_repeats * 90} deg.")
        else:
            print(f"Invalid transform_mode: {actual_transform}. Use 'auto', 'rotate', 'mirror', or 'rotate_wcs_only'")
            continue

        # Create new HDU with transformed data and updated header
        new_hdu = fits.PrimaryHDU(data=data, header=header)

        # Write to output file
        new_hdu.writeto(output_path, overwrite=True)

    if test:
        print("Test mode: All files written with '_test' suffix.")
    else:
        print(f"Successfully transformed {len(files_to_process)} files.")

# Backward compatibility wrapper for download script
def astra_transform(filenames, dirname, n_repeats, test=False):
   """
   Backward compatibility wrapper for existing download script.
   """
   return astra_transform_unified(filenames, transform_mode='auto', n_repeats=n_repeats, test=test,
                                  dirname=dirname)


def main():
    """
    Command line interface for astra_transform.
    """
    if len(sys.argv) < 3:
        print("Usage: python astra_transform.py <input> <transform_mode> [test]")
        print("  input: Directory path, or comma-separated list of directories")
        print("  transform_mode: 'auto', 'rotate', 'mirror', or 'rotate_wcs_only'")
        print("  test: Optional. If present, files will be saved with '_test' suffix")
        print("")
        print("Transform modes:")
        print("  auto: Auto-detect (SPIRIT=mirror, ANDOR=rotate)")
        print("  rotate: Force 90 deg CCW rotation")
        print("  mirror: Force horizontal mirroring")
        print("  rotate_wcs_only: Rotate WCS only (no image change)")
        print("  rotate_image_only: Rotate image only (no WCS change)")
        print("")
        print("Examples:")
        print("  python astra_transform.py /path/to/fits auto")
        print("  python astra_transform.py /path/to/fits rotate test")
        print("  python astra_transform.py /path/dir1,/path/dir2,/path/dir3 mirror")
        sys.exit(1)

    input_arg = sys.argv[1]
    transform_mode = sys.argv[2]
    n_repeats = int(sys.argv[3]) # SHOULD BE SET TO 1 FOR MIRRORING
    completeness = sys.argv[4]
    test_mode = len(sys.argv) > 5 and sys.argv[5].lower() == 'test'

    # Parse input - could be single directory or comma-separated list
    if ',' in input_arg:
        input_data = [path.strip() for path in input_arg.split(',')]
    else:
        input_data = input_arg

    # Validate transform mode
    valid_modes = ['auto', 'rotate', 'mirror', 'rotate_wcs_only', 'rotate_image_only']
    if transform_mode.lower() not in valid_modes:
        print(f"Error: transform_mode must be one of: {', '.join(valid_modes)}")
        sys.exit(1)

    print(f"Input: {input_data}")
    print(f"Transform mode: {transform_mode}")
    print(f"Number of repeats: {n_repeats}")
    print(f"Test mode: {test_mode}")
    print("")

    # Run the transformation
    astra_transform_unified(input_data, transform_mode, n_repeats, completeness, test_mode)

if __name__ == "__main__":
    main()