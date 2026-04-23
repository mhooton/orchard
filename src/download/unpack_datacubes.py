"""
unpack_datacubes.py — SPECULOOS Datacube Unpacker
==================================================

Detects and unpacks FITS datacubes produced by create_datacubes.py into
individual 2D FITS files, restoring original filenames and per-frame header
keywords from the embedded METADATA extension.

Detection
---------
A file is considered a datacube if it contains an HDU named 'METADATA'.
This is specific to the create_datacubes.py convention and is robust to
single-frame cubes (NAXIS3=1) and any other 3D FITS files that may appear.

Header restoration
------------------
Two cube formats are handled, distinguished by whether 'JD-OBS' is present
as a column in the METADATA table:

  New-format cubes (JD-OBS in METADATA, OLDCUBE=N):
    All METADATA columns are written directly to each frame's header using
    the column name as the keyword name, with one exception: FILENAME is
    written to ARCFILE (normalised to ESO dp_id convention) rather than
    as a FILENAME keyword.

  Old-format cubes (JD-OBS absent from METADATA, OLDCUBE=Y):
    Only the six original METADATA columns are available. The following
    keywords are derived from DATE-OBS, EXPTIME, and other header values:
      MJD-OBS, JD-OBS, JD-END, MJD-END, DATE-END, UTC  — from DATE-OBS + EXPTIME
      LST, HA                                            — from DATE-OBS + LONG-OBS + RA
      HJD-OBS, BJD-OBS                                  — from DATE-OBS + RA/DEC + observer location
      ALTITUDE, AZIMUTH                                  — from DATE-OBS + RA/DEC + observer location
    All other instantaneous sensor/telescope keywords (CCD-TEMP, WINDSPD,
    etc.) are inherited from frame 0 of the cube and cannot be recovered.
    OLDCUBE=Y is written to the header to flag this limitation.

Filename normalisation
----------------------
Filenames stored in METADATA may use 'SPECULOOS4' instead of 'SPECU4' and
underscores instead of colons in the timestamp. Both are corrected on output.

Memory efficiency
-----------------
memmap=False is required because BZERO/BSCALE scaling keywords prevent
astropy from memory-mapping integer arrays. The full cube is loaded into
RAM on open. For very large cubes this is unavoidable.

Deletion
--------
The source cube is renamed to a temporary name before unpacking begins, to
prevent any unpacked frame from overwriting it if filenames collide. It is
deleted only after all frames have been successfully written. If the number
of written frames does not match NAXIS3, the cube is renamed back to its
original name and a warning is printed.
"""

import os
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u


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
      - Underscores replaced with colons (e.g. 23_46_11 -> 23:46:11)

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


def _derive_frame_keywords(date_obs, exptime, base_header):
    """
    Derive per-frame time and pointing keywords from DATE-OBS, EXPTIME,
    and fixed header values (RA, DEC, LAT-OBS, LONG-OBS, ALT-OBS).

    Used for old-format cubes where METADATA does not contain per-frame
    engineering keywords.

    Parameters
    ----------
    date_obs    : str   — ISO 8601 datetime string (DATE-OBS value)
    exptime     : float — exposure time in seconds
    base_header : fits.Header — the cube PRIMARY header (provides RA, DEC,
                                LAT-OBS, LONG-OBS, ALT-OBS)

    Returns
    -------
    dict : keyword -> value for all derivable keywords
    """
    derived = {}

    try:
        t_start = Time(date_obs, format='isot', scale='utc')
        t_end   = t_start + exptime * u.second

        # --- Simple time arithmetic -------------------------------------------
        derived['MJD-OBS']  = round(t_start.mjd, 6)
        derived['JD-OBS']   = round(t_start.jd,  6)
        derived['JD-END']   = round(t_end.jd,    6)
        derived['MJD-END']  = round(t_end.mjd,   6)
        derived['DATE-END'] = t_end.isot[:23]

        # UTC seconds since midnight
        midnight = Time(t_start.iso[:10] + 'T00:00:00', format='isot', scale='utc')
        derived['UTC'] = round((t_start - midnight).sec, 3)

        # --- Observer location ------------------------------------------------
        lat  = base_header.get('LAT-OBS',  None)
        lon  = base_header.get('LONG-OBS', None)
        alt  = base_header.get('ALT-OBS',  None)
        ra   = base_header.get('RA',       None)
        dec  = base_header.get('DEC',      None)

        if lat is not None and lon is not None and alt is not None:
            location = EarthLocation(lat=float(lat) * u.deg,
                                     lon=float(lon) * u.deg,
                                     height=float(alt) * u.m)

            # --- LST (seconds) ------------------------------------------------
            try:
                lst_angle = t_start.sidereal_time('apparent', longitude=float(lon) * u.deg)
                derived['LST'] = round(lst_angle.to(u.deg).value * 240.0, 3)  # deg -> seconds
            except Exception as e:
                print(f"    Warning: could not derive LST: {e}")

            if ra is not None and dec is not None:
                target = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg, frame='icrs')

                # --- HA -------------------------------------------------------
                try:
                    ha_angle = lst_angle - target.ra.to(u.hourangle)
                    # Format as hms string matching original convention
                    ha_hours = ha_angle.to(u.hourangle).value
                    ha_h = int(ha_hours)
                    ha_m = int(abs(ha_hours - ha_h) * 60)
                    ha_s = abs((ha_hours - ha_h) * 60 - ha_m) * 60
                    derived['HA'] = f"{ha_h:03d} {ha_m:02d} {ha_s:08.5f}"
                except Exception as e:
                    print(f"    Warning: could not derive HA: {e}")

                # --- ALTITUDE and AZIMUTH ------------------------------------
                try:
                    altaz_frame = AltAz(obstime=t_start, location=location)
                    altaz = target.transform_to(altaz_frame)
                    derived['ALTITUDE'] = round(altaz.alt.deg, 3)
                    derived['AZIMUTH']  = round(altaz.az.deg,  3)
                except Exception as e:
                    print(f"    Warning: could not derive ALTITUDE/AZIMUTH: {e}")

                # --- HJD-OBS and BJD-OBS -------------------------------------
                try:
                    t_hjd = t_start.light_travel_time(target, 'heliocentric', location=location)
                    derived['HJD-OBS'] = round((t_start + t_hjd).jd, 6)
                except Exception as e:
                    print(f"    Warning: could not derive HJD-OBS: {e}")

                try:
                    t_bjd = t_start.light_travel_time(target, 'barycentric', location=location)
                    derived['BJD-OBS'] = round((t_start + t_bjd).jd, 6)
                except Exception as e:
                    print(f"    Warning: could not derive BJD-OBS: {e}")

    except Exception as e:
        print(f"    Warning: keyword derivation failed for {date_obs}: {e}")

    return derived


def unpack_datacube(filepath, output_dir):
    """
    Unpack a datacube into individual 2D FITS files.

    Each frame is written to output_dir using the normalised original filename
    from the METADATA extension (e.g. SPECU4.2025-12-16T04:49:08.869.fits).

    Per-frame header keywords are restored from METADATA columns where
    available. For old-format cubes (no JD-OBS column), time and pointing
    keywords are derived analytically. An OLDCUBE header card is added to
    each unpacked frame to indicate which path was taken.

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

    # Rename the cube to a temporary name before unpacking. This prevents any
    # unpacked frame from overwriting the source cube if their filenames collide
    # (both follow the SPECU#.YYYY-MM-DDTHH:MM:SS.sss.fits convention).
    tmp_filepath = filepath + '.unpacking'
    try:
        os.rename(filepath, tmp_filepath)
    except Exception as e:
        print(f"Error renaming cube to temporary name {os.path.basename(tmp_filepath)}: {e}")
        return [], 0

    try:
        # memmap=False required: BZERO/BSCALE/BLANK keywords prevent astropy
        # from memory-mapping scaled integer arrays.
        hdul = fits.open(tmp_filepath, memmap=False, do_not_scale_image_data=True)
    except Exception as e:
        print(f"Error opening datacube {os.path.basename(filepath)}: {e}")
        os.rename(tmp_filepath, filepath)
        return [], 0

    try:
        if 'METADATA' not in hdul:
            print(f"Error: no METADATA extension found in {os.path.basename(filepath)}")
            os.rename(tmp_filepath, filepath)
            return [], 0

        metadata       = hdul['METADATA'].data
        meta_colnames  = [col.name for col in hdul['METADATA'].columns]
        n_frames_meta  = len(metadata)

        primary_header = hdul[0].header.copy()
        n_frames_cube  = int(primary_header.get('NAXIS3', 0))

        if n_frames_cube == 0:
            print(f"Error: NAXIS3=0 in {os.path.basename(filepath)}, nothing to unpack")
            os.rename(tmp_filepath, filepath)
            return [], 0

        if n_frames_meta != n_frames_cube:
            print(f"Warning: METADATA has {n_frames_meta} rows but NAXIS3={n_frames_cube} "
                  f"in {os.path.basename(filepath)} — will use min of the two")

        n_frames = min(n_frames_meta, n_frames_cube)

        # Determine cube format: new (JD-OBS in METADATA) or old (derive keywords)
        is_old_cube = 'JD-OBS' not in meta_colnames
        if is_old_cube:
            print(f"  Old cube format detected (JD-OBS absent from METADATA) — "
                  f"per-frame keywords will be derived where possible")
        else:
            print(f"  New cube format detected (JD-OBS present in METADATA) — "
                  f"per-frame keywords restored from METADATA")

        # Build a 2D base header from the PRIMARY header:
        # drop NAXIS3, reset NAXIS to 2, strip cube-specific 3D pixel WCS
        # keywords added by create_datacubes (axis 3 only — axis 1/2 pixel WCS
        # keywords are left for wcsfit to overwrite with the real sky WCS).
        base_header = primary_header.copy()
        base_header['NAXIS'] = 2
        for kw in ('NAXIS3', 'CTYPE3', 'CUNIT3', 'CRPIX3', 'CRVAL3', 'CDELT3',
                   'CUBETYPE', 'CUBEKEY', 'NFRAMES', 'CREATED'):
            while kw in base_header:
                del base_header[kw]

        # Patch missing BZERO/BSCALE for BITPIX=16 cubes created before
        # create_datacubes.py was fixed to write these keywords explicitly.
        # Each is checked independently in case only one is absent.
        # BZERO=32768 + BSCALE=1 is the standard FITS uint16 encoding convention
        # for signed int16 data; without these keywords downstream readers cannot
        # recover the correct physical pixel values.
        if base_header.get('BITPIX') == 16:
            if 'BZERO' not in base_header:
                base_header['BZERO'] = (32768, 'Offset for unsigned 16-bit integer encoding')
            if 'BSCALE' not in base_header:
                base_header['BSCALE'] = (1, 'Scale factor for integer encoding')

        written_ids = []

        for i in range(n_frames):
            row = metadata[i]

            # --- Filename: normalise and map to ARCFILE ----------------------
            original_filename = _normalise_filename(str(row['FILENAME']).strip())
            if not original_filename.endswith('.fits'):
                original_filename += '.fits'

            output_path = os.path.join(output_dir, original_filename)

            # --- Extract frame data ------------------------------------------
            try:
                frame_data = np.array(hdul[0].data[i])
            except Exception as e:
                print(f"  Error reading frame {i} from {os.path.basename(filepath)}: {e}")
                continue

            # --- Build per-frame header --------------------------------------
            frame_header = base_header.copy()

            # Write ARCFILE from normalised FILENAME
            frame_header['ARCFILE'] = (original_filename,
                                       'Archive File Name')

            if is_old_cube:
                # Old format: write the six available METADATA columns directly,
                # then derive time/pointing keywords analytically.
                date_obs = str(row['DATE-OBS'])
                exptime  = float(row['EXPTIME'])

                frame_header['DATE-OBS'] = date_obs
                frame_header['EXPTIME']  = exptime
                frame_header['FILTER']   = str(row['FILTER'])
                frame_header['OBJECT']   = str(row['OBJECT'])
                frame_header['AIRMASS']  = float(row['AIRMASS'])

                derived = _derive_frame_keywords(date_obs, exptime, base_header)
                for kw, val in derived.items():
                    frame_header[kw] = val

                frame_header['OLDCUBE'] = ('Y',
                    'Old cube data, site headers do not vary')

            else:
                # New format: write all METADATA columns directly to header,
                # skipping FILENAME (handled via ARCFILE above).
                for colname in meta_colnames:
                    if colname == 'FILENAME':
                        continue
                    try:
                        val = row[colname]
                        # Convert numpy scalars to native Python types
                        if isinstance(val, (np.floating, np.integer)):
                            val = val.item()
                        elif isinstance(val, bytes):
                            val = val.decode('utf-8').strip()
                        else:
                            val = str(val).strip()
                        frame_header[colname] = val
                    except Exception as e:
                        print(f"  Warning: could not write column {colname} "
                              f"to frame {i} header: {e}")

                frame_header['OLDCUBE'] = ('N',
                    'Old cube data, site headers do not vary')

            # --- Write frame to disk -----------------------------------------
            try:
                new_hdu = fits.PrimaryHDU(data=frame_data, header=frame_header)
                new_hdu._do_not_scale_image_data = True
                fits.HDUList([new_hdu]).writeto(output_path, overwrite=True)
                dp_id = original_filename.replace('.fits', '')
                written_ids.append(dp_id)
                print(f"  ✓ {i + 1}/{n_frames} {original_filename}")
            except Exception as e:
                print(f"  Error writing frame {i} ({original_filename}): {e}")

    finally:
        hdul.close()

    n_written = len(written_ids)

    # Delete the cube only if all frames were written successfully.
    # The cube is currently named tmp_filepath; restore the original name
    # if deletion is not appropriate.
    if n_written == n_frames:
        try:
            os.remove(tmp_filepath)
            print(f"Deleted cube: {os.path.basename(filepath)}")
        except Exception as e:
            print(f"Warning: could not delete cube {os.path.basename(filepath)}: {e}")
    else:
        # Restore original filename so the cube is not silently lost
        try:
            os.rename(tmp_filepath, filepath)
        except Exception as e:
            print(f"Warning: could not restore cube filename {os.path.basename(filepath)}: {e}")
        print(f"Warning: only {n_written}/{n_frames} frames written — "
              f"cube retained: {os.path.basename(filepath)}")

    print(f"Unpacked {n_written} frames from {os.path.basename(filepath)}")
    return written_ids, n_written