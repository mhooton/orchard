import math
import numpy as np
import sys
import os
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import glob

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from download.raw_provenance import stamp_ingest  # noqa: E402
from utils import orchard_version  # noqa: E402

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
            if not os.path.exists(f"{imgdir}/{down_file}.fits") and os.path.exists(f"{imgdir}/{down_file}.fits.fz"):
                # stored Rice-compressed (utils/framepack.py rawpack), which happens only after the night's first
                # complete run, so after its transform; a compressed original is never deleted or fetched again here
                print(f"{down_file} is stored compressed (.fits.fz). No need to redownload.")
                transformed_files.append(down_file)
                continue
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


# The pixel transforms below move the image, so the WCS has to follow. Every one of them is linear in pixel
# offsets: an offset d' in the output image came from the offset R d' in the input, so a linear WCS matrix M
# (pixel offset -> intermediate world coordinates) becomes M R, and CRPIX moves to where the input CRPIX lands.
#
# Astra/PinPoint write the matrix as PCi_j (no CDELTi), astrometry.net and CASU as CDi_j, astropy as PCi_j with
# CDELTi = 1. The earlier code only knew CDi_j (with a formula that is right only for an unflipped matrix) and
# CDELTi/CROTA2, so on every Astra frame it moved CRPIX and left the PC matrix as it was: the stored WCS was 90
# degrees out (up to ~12 arcmin at the corners of an Andor frame). Negating CDELT1 for the mirror flipped the sky
# axis, not the pixel axis, which is only the same thing when the matrix has no rotation.
ROT90_ACW = ((0, -1), (1, 0))   # np.rot90(k=1): X' = Y, Y' = NAXIS1 + 1 - X, so dX = -dY', dY = dX'
MIRROR_X = ((-1, 0), (0, 1))    # np.flip(axis=1): X' = NAXIS1 + 1 - X

WCS_LINEAR_KEYS = ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2')
WCS_CELESTIAL_KEYS = ('WCSAXES', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CUNIT1', 'CUNIT2',
                      'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2', 'LONPOLE', 'LATPOLE', 'IMAGEW', 'IMAGEH') + WCS_LINEAR_KEYS


def linear_matrix(header):
    """
    The header's pixel -> intermediate-world matrix (degrees per pixel) and the convention it was written in,
    following wcslib's precedence: PCi_j (scaled by CDELTi, default 1) when any PCi_j is present, even if CDi_j
    is there too; then CDi_j; then CDELTi with CROTA2. Returns (None, None) when there is no matrix.
    """
    if any(k in header for k in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2')):
        # A missing diagonal PC element is 1 and a missing off-diagonal one 0 (astropy omits zeros on writing)
        pc = np.array([[header.get('PC1_1', 1.0), header.get('PC1_2', 0.0)],
                       [header.get('PC2_1', 0.0), header.get('PC2_2', 1.0)]], dtype=float)
        cdelt = np.array([header.get('CDELT1', 1.0), header.get('CDELT2', 1.0)], dtype=float)
        return np.diag(cdelt) @ pc, 'PC'
    if any(k in header for k in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2')):
        return np.array([[header.get('CD1_1', 0.0), header.get('CD1_2', 0.0)],
                         [header.get('CD2_1', 0.0), header.get('CD2_2', 0.0)]], dtype=float), 'CD'
    if 'CDELT1' in header and 'CDELT2' in header:
        rho = np.radians(header.get('CROTA2', 0.0))
        c1, c2 = float(header['CDELT1']), float(header['CDELT2'])
        return np.array([[c1 * np.cos(rho), -c2 * np.sin(rho)],
                         [c1 * np.sin(rho), c2 * np.cos(rho)]]), 'CROTA'
    return None, None


def set_linear_matrix(header, m, convention):
    """Write m back in the convention it was read in. PC keeps the existing CDELTi; CROTA becomes PC + CDELTi."""
    if convention in ('PC', 'CROTA'):
        cdelt = np.array([header.get('CDELT1', 1.0), header.get('CDELT2', 1.0)], dtype=float)
        pc = np.diag(1.0 / cdelt) @ m
        for (i, j), v in np.ndenumerate(pc):
            header['PC{}_{}'.format(i + 1, j + 1)] = float(v)
        for k in ('CROTA1', 'CROTA2'):
            if k in header:
                del header[k]
    elif convention == 'CD':
        for (i, j), v in np.ndenumerate(m):
            header['CD{}_{}'.format(i + 1, j + 1)] = float(v)


def plate_pa(m):
    """
    The PA keyword as the SPECULOOS transfer (headerfix.py) defines it: the position angle, from north through
    east, of the -Y pixel axis (towards row 1), for either parity. Recomputed so that it stays true of the
    stored orientation after a transform.
    """
    det = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]
    parity = 1.0 if det >= 0 else -1.0
    t = parity * m[0, 0] + m[1, 1]
    a = parity * m[1, 0] - m[0, 1]
    return (180.0 - math.degrees(math.atan2(a, t))) % 360.0


def drop_singular_wcs(header):
    """
    Remove a celestial WCS whose matrix is singular, and the PA computed from it. SPIRIT frames reach ESO with an
    all-zero CD block (PA = 180.0); left in place, it reads as 'already solved' to astrometry.py and as an error to
    astropy. Returns True if anything was removed.
    """
    m, _ = linear_matrix(header)
    if m is None or not str(header.get('CTYPE1', '')).startswith('RA'):
        return False
    if abs(np.linalg.det(m)) > 0:
        return False
    for k in WCS_CELESTIAL_KEYS + ('PA',):
        if k in header:
            del header[k]
    header.add_history('Singular WCS matrix removed (no sky solution); PA removed with it')
    return True


def transform_sip(header, r):
    """
    Carry SIP distortion (A/B forward, AP/BP inverse) through the pixel transform r. With pixel offsets
    d = r d', the distortion F(d) becomes r^T F(r d'), which for a signed permutation r is a relabelling of the
    coefficients: rotation swaps A and B (and their orders), the mirror flips signs.
    """
    r = np.array(r, dtype=int)
    for a, b in (('A', 'B'), ('AP', 'BP')):
        if a + '_ORDER' not in header and b + '_ORDER' not in header:
            continue
        polys = []
        for name in (a, b):
            order = int(header.get(name + '_ORDER', 0))
            coeffs = {}
            for p in range(order + 1):
                for q in range(order + 1 - p):
                    key = '{}_{}_{}'.format(name, p, q)
                    if key in header:
                        coeffs[(p, q)] = float(header[key])
                        del header[key]
            polys.append((order, coeffs))
        # P'(u', v') = P(u, v) with u = s1 d'[i1], v = s2 d'[i2], (u, v) = r (u', v')
        i1, i2 = int(np.flatnonzero(r[0])[0]), int(np.flatnonzero(r[1])[0])
        s1, s2 = int(r[0, i1]), int(r[1, i2])
        subst = []
        for order, coeffs in polys:
            out = {}
            for (p, q), c in coeffs.items():
                pq = (p, q) if (i1, i2) == (0, 1) else (q, p)
                out[pq] = out.get(pq, 0.0) + (s1 ** p) * (s2 ** q) * c
            subst.append((order, out))
        rt = r.T
        for row, name in enumerate((a, b)):
            src = int(np.nonzero(rt[row])[0][0])
            sign = int(rt[row, src])
            order, coeffs = subst[src]
            header[name + '_ORDER'] = order
            for (p, q), c in coeffs.items():
                header['{}_{}_{}'.format(name, p, q)] = sign * c


def transform_wcs(header, r, new_crpix):
    """Apply the pixel transform r (see ROT90_ACW) to the linear matrix, move CRPIX and recompute PA."""
    m, convention = linear_matrix(header)
    if 'CRPIX1' in header and 'CRPIX2' in header:
        header['CRPIX1'], header['CRPIX2'] = new_crpix(float(header['CRPIX1']), float(header['CRPIX2']))
    if m is None:
        return
    transform_sip(header, r)
    m = m @ np.array(r, dtype=float)
    set_linear_matrix(header, m, convention)
    if 'PA' in header:
        header['PA'] = (plate_pa(m), '[deg, 0-360 CCW] Position angle of plate')
    # Frames transformed before this marker existed kept Astra's PC matrix; repair_astra_wcs.py looks for it
    header['ASTRAWCS'] = (2, 'WCS matrix transformed with the image')


def rotate_wcs_acw(header, n_rotations=1):
    """
    Update the WCS for np.rot90(data, k=n_rotations), i.e. 90 deg anticlockwise per step as displayed with row 1
    at the top (the MaxIm 'Rotate 90 CCW' orientation). header must still describe the input image.
    """
    nx = header.get('NAXIS1', 0)
    ny = header.get('NAXIS2', 0)

    for rot in range(n_rotations):
        transform_wcs(header, ROT90_ACW, lambda c1, c2, nx=nx: (c2, nx - c1 + 1))

        # Update non-standard dimension keywords if present
        if 'IMAGEW' in header and 'IMAGEH' in header:
            header['IMAGEW'], header['IMAGEH'] = header['IMAGEH'], header['IMAGEW']

        nx, ny = ny, nx


def mirror_wcs_horizontal(header):
    """
    Update WCS headers for horizontal mirroring (flip left-right).
    This flips the image around the vertical axis (y-axis).
    """
    nx = header.get('NAXIS1', 0)
    transform_wcs(header, MIRROR_X, lambda c1, c2: (nx - c1 + 1, c2))


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
                            test=False, dirname=None, force=False, ingest=False):
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
    force : bool, optional
        Transform frames that already carry ASTRAROT (rotate, rotate_wcs_only) or ASTRAMIR (mirror).
        By default they are skipped, so running the transform twice cannot turn a frame twice.
    ingest : bool, optional
        True when the frames have just been downloaded or unpacked (the astra_transform() wrapper): each frame
        written is stamped as new, with ORCHVER and the card 'orchard <version> astra_transform'. False for a later
        edit of frames already in the archive (the command line): ORCHMOD and an 'orchard <version>
        astra_transform: <mode>' card. Either way the descriptive HISTORY cards below are written as before (the
        repair scripts count them) and the file is written with a fresh CHECKSUM/DATASUM.
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
    skipped = 0

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
        if is_light_frame and drop_singular_wcs(header):
            print(f"Removed the singular WCS of {name_without_ext}; astrometry.py will solve it.")

        # Determine transformation type
        if transform_mode == 'auto':
            # Auto-detect based on NAXIS1 (backward compatibility)
            is_spirit = (naxis1 == 1024)
            actual_transform = 'mirror' if is_spirit else 'rotate'
        else:
            actual_transform = transform_mode.lower()

        # A frame that already carries the flag this transform sets has been transformed: doing it again
        # turns it a further 90 deg (or flips it back). On 2025-12-09 this left 2,441 Ganymede frames
        # 180 deg from the ESO orientation.
        done_flag = {'rotate': 'ASTRAROT', 'rotate_wcs_only': 'ASTRAROT', 'mirror': 'ASTRAMIR'}.get(actual_transform)
        if done_flag and done_flag in header and not force:
            print(f"Skipped {name_without_ext}: already has {done_flag}, not transforming it again.")
            skipped += 1
            continue

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

        if ingest:
            stamp_ingest(header, 'astra_transform')
        else:
            what = actual_transform if actual_transform == 'mirror' else '{} k={}'.format(actual_transform, n_repeats)
            orchard_version.record_edit(header, 'astra_transform', what)

        # Create new HDU with transformed data and updated header
        new_hdu = fits.PrimaryHDU(data=data, header=header)

        # Write to output file. checksum=True replaces the CHECKSUM/DATASUM ESO computed for the untransformed image,
        # which astropy would otherwise write back unchanged (and wrong).
        new_hdu.writeto(output_path, overwrite=True, checksum=True)

    if test:
        print("Test mode: All files written with '_test' suffix.")
    else:
        print(f"Successfully transformed {len(files_to_process) - skipped} files"
              + (f" ({skipped} already transformed, skipped)." if skipped else "."))

# Backward compatibility wrapper for download script
def astra_transform(filenames, dirname, n_repeats, test=False):
   """
   The ingest transform, for frames just downloaded or unpacked (SSO_download.py, the queue's fetch): each frame
   it writes is stamped as new (ingest=True).
   """
   return astra_transform_unified(filenames, transform_mode='auto', n_repeats=n_repeats, test=test,
                                  dirname=dirname, ingest=True)


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

    # Run the transformation: frames already in the archive, so each one written is recorded as an edit
    astra_transform_unified(input_data, transform_mode, n_repeats, completeness, test_mode, ingest=False)

if __name__ == "__main__":
    main()