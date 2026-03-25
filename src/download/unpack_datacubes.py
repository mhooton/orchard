"""
unpack_datacubes.py — SPECULOOS Datacube Unpacker
==================================================

Detects and unpacks FITS datacubes produced by create_datacubes.py into
individual 2D FITS files, restoring original filenames from the embedded
METADATA extension.

Detection
---------
A file is considered a datacube if it contains an HDU named 'METADATA'.
This is specific to the create_datacubes.py convention and is robust to
single-frame cubes (NAXIS3=1) and any other 3D FITS files that may appear.

Memory efficiency
-----------------
Frames are extracted one at a time. The file is opened with memmap=False
because BZERO/BSCALE scaling keywords prevent astropy from memory-mapping
integer arrays. Peak RAM usage is approximately one full cube at a time.

Filename normalisation
----------------------
Filenames stored in the METADATA table may differ from the original ESO
dp_id convention (e.g. 'SPECULOOS4' instead of 'SPECU4', underscores
instead of colons in the timestamp). These are normalised on output so
that unpacked files match the expected dp_id format.

Deletion
--------
The source cube is deleted only after all frames have been successfully
written. If the number of written frames does not match NAXIS3, the cube
is retained and a warning is printed.
"""

import os
import numpy as np
from astropy.io import fits


def is_datacube(filepath):
    """
    Return True if the FITS file at filepath is a datacube produced by
    create_datacubes.py, identified by the presence of a 'METADATA' extension.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    bool
    """
    try:
        with fits.open(filepath, memmap=False) as hdul:
            return 'METADATA' in hdul
    except Exception as e:
        print(f"Warning: could not inspect {os.path.basename(filepath)}: {e}")
        return False


def _normalise_filename(filename):
    """
    Normalise a filename from the METADATA table to match the ESO dp_id
    convention used by the downloader.

    Applies two corrections:
      - 'SPECULOOS' -> 'SPECU'  (e.g. SPECULOOS4 -> SPECU4)
      - Underscores in the timestamp portion replaced with colons
        (e.g. 23_46_11 -> 23:46:11)

    Parameters
    ----------
    filename : str

    Returns
    -------
    str
    """
    filename = filename.replace('SPECULOOS', 'SPECU')
    filename = filename.replace('_', ':')
    return filename


def unpack_datacube(filepath, output_dir):
    """
    Unpack a datacube into individual 2D FITS files.

    Each frame is written to output_dir using the original filename stored
    in the METADATA extension, normalised to the ESO dp_id convention
    (e.g. SPECU4.2025-12-16T04:49:08.869.fits).

    The PRIMARY header is used as the base header for all frames, with
    per-frame values (DATE-OBS, EXPTIME, FILTER, OBJECT, AIRMASS) overwritten
    from the METADATA table for each respective frame.

    The source cube is deleted only if all frames are written successfully.

    Parameters
    ----------
    filepath   : str  — path to the datacube FITS file
    output_dir : str  — directory to write unpacked frames into

    Returns
    -------
    tuple : (list of str, int)
        - dp_ids (filenames without .fits extension) of successfully written
          frames, suitable for passing directly to astra_transform
        - total number of frames successfully written
    """
    print(f"Unpacking datacube: {os.path.basename(filepath)}")

    try:
        # memmap=False required: BZERO/BSCALE/BLANK keywords prevent astropy
        # from memory-mapping scaled integer arrays.
        hdul = fits.open(filepath, memmap=False)
    except Exception as e:
        print(f"Error opening datacube {os.path.basename(filepath)}: {e}")
        return [], 0

    try:
        if 'METADATA' not in hdul:
            print(f"Error: no METADATA extension found in {os.path.basename(filepath)}")
            return [], 0

        metadata = hdul['METADATA'].data
        n_frames_meta = len(metadata)

        primary_header = hdul[0].header.copy()
        n_frames_cube = int(primary_header.get('NAXIS3', 0))

        if n_frames_cube == 0:
            print(f"Error: NAXIS3=0 in {os.path.basename(filepath)}, nothing to unpack")
            return [], 0

        if n_frames_meta != n_frames_cube:
            print(f"Warning: METADATA has {n_frames_meta} rows but NAXIS3={n_frames_cube} "
                  f"in {os.path.basename(filepath)} — will use min of the two")

        n_frames = min(n_frames_meta, n_frames_cube)

        # Build a 2D base header from the PRIMARY header:
        # drop NAXIS3 and reset NAXIS to 2, strip the dummy pixel WCS
        # (CTYPE3, CUNIT3, CRPIX3, CRVAL3, CDELT3) added by create_datacubes.
        base_header = primary_header.copy()
        base_header['NAXIS'] = 2
        for kw in ('NAXIS3', 'CTYPE3', 'CUNIT3', 'CRPIX3', 'CRVAL3', 'CDELT3',
                   'CUBETYPE', 'CUBEKEY', 'NFRAMES', 'CREATED'):
            while kw in base_header:
                del base_header[kw]

        written_ids = []

        for i in range(n_frames):
            row = metadata[i]

            original_filename = _normalise_filename(str(row['FILENAME']).strip())
            if not original_filename.endswith('.fits'):
                original_filename = original_filename + '.fits'

            output_path = os.path.join(output_dir, original_filename)

            try:
                frame_data = np.array(hdul[0].data[i])
            except Exception as e:
                print(f"  Error reading frame {i} from {os.path.basename(filepath)}: {e}")
                continue

            # Build per-frame header from base, overwriting per-frame values
            frame_header = base_header.copy()
            frame_header['DATE-OBS'] = str(row['DATE-OBS'])
            frame_header['EXPTIME']  = float(row['EXPTIME'])
            frame_header['FILTER']   = str(row['FILTER'])
            frame_header['OBJECT']   = str(row['OBJECT'])
            frame_header['AIRMASS']  = float(row['AIRMASS'])

            try:
                fits.writeto(output_path, frame_data, frame_header, overwrite=True)
                dp_id = original_filename.replace('.fits', '')
                written_ids.append(dp_id)
                print(f"  ✓ {i + 1}/{n_frames} {original_filename}")
            except Exception as e:
                print(f"  Error writing frame {i} ({original_filename}): {e}")

    finally:
        hdul.close()

    n_written = len(written_ids)

    # Delete the cube only if all frames were written successfully
    if n_written == n_frames:
        try:
            os.remove(filepath)
            print(f"Deleted cube: {os.path.basename(filepath)}")
        except Exception as e:
            print(f"Warning: could not delete cube {os.path.basename(filepath)}: {e}")
    else:
        print(f"Warning: only {n_written}/{n_frames} frames written — "
              f"cube retained: {os.path.basename(filepath)}")

    print(f"Unpacked {n_written} frames from {os.path.basename(filepath)}")
    return written_ids, n_written