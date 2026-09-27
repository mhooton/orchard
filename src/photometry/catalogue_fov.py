#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
catalogue_fov.py — Stage T8

Field-of-View Catalogue Creation

Creates a deep stacked image from a subset of science frames and generates
a source catalogue for the field of view. The catalogue is crossmatched
against Gaia DR3, and all identified targets (primary from schedule/TOI,
secondary from the 40pc target list) are recorded in the Gaia_Crossmatch
extension with their roles.

Processing workflow:
    1.  Select images from the night using integration time threshold
    2.  WCS solve selected images (unless --no-wcs)
    3.  Filter out frames with failed WCS solutions
    4.  Stack surviving images via casutools.imstack
    5.  Detect sources in stack via casutools.imcore
    6.  Crossmatch detected sources with Gaia DR3
    7.  Guard: abort if crossmatch produced no sources
    8.  Identify targets in the FOV via target_management
    9.  Write TARGET_ROLE column into the Gaia_Crossmatch extension
    10. Per identified target: validate catalogue, manage backups,
        and signal ZLP_pipeline if no targets could be identified

Exit codes:
    0 — success (primary target present, or secondary targets present)
    1 — no targets identified in FOV, or fatal error before crossmatch

Output files (named by Gaia DR2 ID where available):
    {dr2_id}_outstack_{filter}.fits         — deep stacked image
    {dr2_id}_outstack_{filter}_conf.fits    — confidence map
    {dr2_id}_stack_catalogue_{filter}.fits  — source catalogue with
                                              Gaia crossmatch and target roles

This is Stage T8 of the pipeline. It produces the reference catalogues
required for aperture photometry (Stage T9).
"""

import os
import sys
import math
import argparse
import logging
from astropy.io import fits
import numpy as np

import photometry.casutools as casutools
import photometry.gaia_dr2_test as gaia_dr2
from photometry.wcs_status import wcs_succeeded
from photometry.wcs_fitting import m_solve_images
from reporting.QC import main as QC
from calibration.pipeutils import detect_instrument
from utils.target_management import (
    identify_targets,
    load_id_aliases,
    apply_id_alias,
    get_target_from_target_list_by_name,
    catalogue_is_valid,
    update_backup,
    restore_backup,
    find_backup,
    get_output_filename,
    get_target_from_schedule,
    get_target_from_toi,
)

logging.basicConfig(level=logging.INFO, format='%(levelname)7s %(message)s')
logger = logging.getLogger(__name__)

print(os.getcwd())

# ---------------------------------------------------------------------------
# WCS header check
# ---------------------------------------------------------------------------

def _has_wcs_headers(image_path):
    """
    Return True if the image has the minimum WCS headers needed to
    attempt astrometric solving (CRVAL1 and CRVAL2 present and non-zero).

    Inputs:
        image_path : str

    Output:
        bool
    """
    try:
        hdr = fits.getheader(image_path)
        return (hdr.get('CRVAL1') not in (None, 0) and
                hdr.get('CRVAL2') not in (None, 0))
    except Exception as e:
        logger.warning("Could not read header for '%s': %s", image_path, e)
        return False


# ---------------------------------------------------------------------------
# Image selection
# ---------------------------------------------------------------------------

def select_images(all_images):
    """
    Select the subset of images to include in the stack.

    The target number of images N is determined by the integration time
    threshold: N = max(ceil(1200 / EXPTIME), 40). A minimum of 40 images
    and a maximum of 400 images
    is always targeted regardless of exposure time. If the total number of
    available images is less than N, all images are used as the target.

    Selection proceeds in two passes:
        Pass 1: Start at the midpoint of the night, work forward (later
                in the night), adding images that pass the WCS header
                check until N is reached or the end of the list.
        Pass 2: If N is not yet reached, start at midpoint-1, work
                backwards (earlier in the night), adding further
                passing images until N is reached or the beginning
                of the list.

    If N is still not reached after both passes, all passing images
    are used and a warning is logged stating the actual total
    integration time achieved.

    Inputs:
        all_images : list of str  All candidate image paths, chronologically
                                  sorted

    Output:
        list of str  Selected image paths in the order they were added
    """
    if not all_images:
        return []

    # Determine target N from EXPTIME of first image
    exptime = None
    for candidate in all_images:
        try:
            exptime = fits.getheader(candidate).get('EXPTIME')
            if exptime and float(exptime) > 0:
                exptime = float(exptime)
                break
        except Exception:
            continue

    if exptime is None:
        logger.warning(
            "Could not read EXPTIME from any image — using all images"
        )
        n_target = len(all_images)
    else:
        n_target = min(max(math.ceil(1200.0 / exptime), 40), 400, len(all_images))
        logger.info(
            "EXPTIME=%.1fs → target stack size N=%d (threshold 1200s, minimum 40, maximum 400)",
            exptime, n_target
        )

    midpoint = len(all_images) // 2
    selected = []

    # Pass 1: midpoint → end
    for image in all_images[midpoint:]:
        if len(selected) >= n_target:
            break
        if _has_wcs_headers(image):
            selected.append(image)
        else:
            logger.warning(
                "Skipping '%s' — missing WCS headers", image
            )

    # Pass 2: midpoint-1 → start (working backwards)
    if len(selected) < n_target:
        for image in reversed(all_images[:midpoint]):
            if len(selected) >= n_target:
                break
            if _has_wcs_headers(image):
                selected.append(image)
            else:
                logger.warning(
                    "Skipping '%s' — missing WCS headers", image
                )

    actual_exptime = len(selected) * (exptime or 0)
    if len(selected) < n_target:
        logger.warning(
            "Integration time threshold not met: selected %d image(s), "
            "total integration %.0fs (target 1200s, minimum 40 images, maximum 400)",
            len(selected), actual_exptime
        )
    else:
        logger.info(
            "Selected %d image(s) for stack (total integration ~%.0fs)",
            len(selected), actual_exptime
        )

    return selected


# ---------------------------------------------------------------------------
# Gaia_Crossmatch extension — target role column
# ---------------------------------------------------------------------------

def write_target_roles(catalogue_path, targets):
    """
    Write a TARGET_ROLE column into the Gaia_Crossmatch extension of the
    stack catalogue FITS file.

    Each row receives one of:
        'primary'   — the scheduled (or TOI) target
        'secondary' — a 40pc target list star in the FOV
        ''          — not an identified target

    Inputs:
        catalogue_path : str   Path to the stack catalogue FITS file
        targets        : list  of (DR2_ID, role, Teff) from identify_targets
    """
    target_role_map = {dr2_id: role for (dr2_id, role, _) in targets}

    try:
        with fits.open(catalogue_path, mode='update') as hdul:
            ext_names = [hdu.name.upper() for hdu in hdul]

            if 'GAIA_CROSSMATCH' not in ext_names:
                logger.warning(
                    "Cannot write TARGET_ROLE — Gaia_Crossmatch extension "
                    "absent from '%s'", catalogue_path
                )
                return

            gaia_ext = hdul['GAIA_CROSSMATCH']
            dr2_ids = gaia_ext.data['GAIA_DR2_ID']
            if 'GAIA_DR3_ID' in gaia_ext.columns.names:
                dr3_ids = gaia_ext.data['GAIA_DR3_ID']
            else:
                dr3_ids = [''] * len(dr2_ids)
            # Targets may be identified by a DR2 or a DR3 ID; a row takes
            # the role of whichever of its IDs is in the target map.
            roles = np.array([
                target_role_map.get(str(d2).strip(), '')
                or target_role_map.get(str(d3).strip(), '')
                for d2, d3 in zip(dr2_ids, dr3_ids)
            ])

            existing_cols = gaia_ext.columns
            existing_header = gaia_ext.header
            new_col = fits.ColDefs([
                fits.Column(name='TARGET_ROLE', format='10A', array=roles)
            ])

            if 'TARGET_ROLE' in existing_cols.names:
                keep = fits.ColDefs([
                    c for c in existing_cols if c.name != 'TARGET_ROLE'
                ])
                new_table = fits.BinTableHDU.from_columns(
                    keep + new_col, name='GAIA_CROSSMATCH'
                )
            else:
                new_table = fits.BinTableHDU.from_columns(
                    existing_cols + new_col, name='GAIA_CROSSMATCH'
                )

            # Preserve non-structural header keywords from the original
            # extension (e.g. N_MATCH written by gaia_dr2_test.crossmatch).
            # from_columns only carries column-definition keywords, so any
            # additional keywords must be copied explicitly.
            skip = set(new_table.header.keys()) | {'XTENSION', 'BITPIX',
                       'NAXIS', 'NAXIS1', 'NAXIS2', 'PCOUNT', 'GCOUNT',
                       'TFIELDS', 'EXTNAME', 'END', ''}
            for key, value, comment in existing_header.cards:
                if key not in skip and not key.startswith('TTYPE') \
                        and not key.startswith('TFORM') \
                        and not key.startswith('TUNIT'):
                    new_table.header.set(key, value, comment)

            del hdul['GAIA_CROSSMATCH']
            hdul.append(new_table)
            hdul.flush()

    except Exception as e:
        logger.error(
            "Failed to write TARGET_ROLE to '%s': %s", catalogue_path, e
        )
        return

    logger.info(
        "TARGET_ROLE column written for %d identified target(s)",
        len(targets)
    )


def read_catalogue_info(path):
    """
    Read what target identification needs from a stack catalogue:
    both Gaia ID columns of Gaia_Crossmatch, plus RA/Dec (degrees) and
    isophotal flux of every detection from the imcore table.
    """
    with fits.open(path) as hdul:
        g = hdul['GAIA_CROSSMATCH'].data
        dr2 = [str(x) for x in g['GAIA_DR2_ID']]
        dr3 = ([str(x) for x in g['GAIA_DR3_ID']]
               if 'GAIA_DR3_ID' in g.columns.names else None)
        coords = None
        flux = None
        try:
            apm = hdul[1].data
            names = [n.upper() for n in apm.columns.names]
            ra = np.asarray(apm[apm.columns.names[names.index('RA')]], dtype=float)
            dec = np.asarray(apm[apm.columns.names[names.index('DEC')]], dtype=float)
            if np.nanmax(np.abs(ra)) <= 2 * np.pi + 1e-6:   # radians
                ra, dec = np.degrees(ra), np.degrees(dec)
            coords = np.column_stack([ra, dec])
            if 'ISOPHOTAL_FLUX' in names:
                flux = np.asarray(apm[apm.columns.names[names.index('ISOPHOTAL_FLUX')]], dtype=float)
            if len(coords) != len(dr2):
                coords, flux = None, None
        except Exception as e:
            logger.warning("Could not read detection positions from '%s': %s", path, e)
    return {'dr2': dr2, 'dr3': dr3, 'coords': coords, 'flux': flux}


def inject_primary_identity(catalogue_path, row, gaia_id, db_row=None):
    """
    Write the schedule's Gaia ID (and, if the star is in the local
    database, its proper motion, parallax, magnitudes and Teff) into the
    Gaia_Crossmatch row that the coordinate fallback identified as the
    primary target.  Only blank/NaN cells are filled; nothing is overridden.
    The extension header records the injection.
    """
    try:
        with fits.open(catalogue_path, mode='update') as hdul:
            ext = hdul['GAIA_CROSSMATCH']
            g = ext.data
            cols = list(g.columns.names)
            id_col = 'GAIA_DR3_ID' if 'GAIA_DR3_ID' in cols else 'GAIA_DR2_ID'
            if str(g[id_col][row]).strip().lower() in ('', 'nan', 'none'):
                g[id_col][row] = str(gaia_id)
            if db_row is not None:
                for col, key in (('PMRA', 'pmra'), ('PMDEC', 'pmdec'),
                                 ('PARALLAX', 'parallax'),
                                 ('GMAG', 'phot_g_mean_mag'), ('G_RP', 'g_rp'),
                                 ('BP_RP', 'bp_rp'), ('TEFF', 'teff_gspphot')):
                    val = db_row.get(key)
                    if col in cols and val is not None:
                        try:
                            if not np.isfinite(float(g[col][row])):
                                g[col][row] = float(val)
                        except (ValueError, TypeError):
                            pass
                if 'GAIA_DR2_ID' in cols and db_row.get('dr2_source_id') \
                        and str(g['GAIA_DR2_ID'][row]).strip().lower() in ('', 'nan', 'none'):
                    g['GAIA_DR2_ID'][row] = str(db_row['dr2_source_id'])
            ext.header['PRIMINJ'] = (str(gaia_id),
                                     'primary ID injected by coordinate identification')
            ext.header['PRIMIROW'] = (int(row), 'catalogue row of injected primary')
            hdul.flush()
        logger.info("Injected Gaia ID %s into row %d of '%s'", gaia_id, row, catalogue_path)
        return True
    except Exception as e:
        logger.error("Failed to inject primary identity into '%s': %s", catalogue_path, e)
        return False


def record_primary_provenance(catalogue_path, info):
    """Record how the primary target was identified in the Gaia_Crossmatch header."""
    try:
        with fits.open(catalogue_path, mode='update') as hdul:
            h = hdul['GAIA_CROSSMATCH'].header
            h['PRIM_ST'] = (str(info.get('status', ''))[:60], 'primary identification status')
            h['PRIM_HOW'] = (str(info.get('primary_method') or 'none')[:20], 'how the primary was identified')
            h['PRIM_ID'] = (str(info.get('primary_id') or '')[:30], 'primary Gaia ID used')
            if info.get('schedule_id'):
                h['SCHED_ID'] = (str(info['schedule_id'])[:30], 'Gaia ID from schedule (after aliases)')
            hdul.flush()
    except Exception as e:
        logger.warning("Could not record primary provenance in '%s': %s", catalogue_path, e)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(filelist, outdir, backupcatdir, reportdir, filter, date,
         obsdir, targname, target_list_path, toi_table_path,
         telescope='unknown',
         no_wcs=False,
         c_thresh=2,
         s_thresh=20,
         verbose=False,
         ipix=6,
         catsrc='vizgaia3',
         rcore=4,
         ncores=1,
         ext='fits'):
    """
    Main entry point for FOV catalogue creation.

    Returns:
        int  Exit code: 0 for success, 1 if no targets found or fatal error.
    """
    ncores = int(ncores)

    # Change to data directory so catcache is written in the correct place
    os.chdir(filelist.split('output')[0])

    # Detect instrument using pipeutils.detect_instrument, which matches on
    # image dimensions against the instrument config. Images are trimmed at
    # this point in the pipeline so trimmed=True is required.
    instrument = 'unknown'
    try:
        with open(filelist, 'r') as _f:
            _first = [l.strip() for l in _f if l.strip()][0]
        with fits.open(_first) as _hdul:
            instrument = detect_instrument(_hdul, trimmed=True)
        logger.info("Instrument detected: '%s'", instrument)
    except Exception as _e:
        logger.warning("Could not detect instrument: %s", _e)

    stack_filelist = os.path.dirname(filelist) + '/' + filter + '_stacked.dat'

    # ------------------------------------------------------------------
    # Read full image list
    # ------------------------------------------------------------------
    with open(filelist, 'r') as f:
        all_images = [l.strip() for l in f if l.strip()]

    if not all_images:
        logger.error("Image list '%s' is empty — aborting", filelist)
        return 1

    # ------------------------------------------------------------------
    # Step 1: identify primary target early for output file naming.
    # This is done before stacking so that output files can be named
    # correctly from the outset, with no temporary names or renames.
    # identify_targets() is still called after the crossmatch to validate
    # presence in the catalogue and identify secondary targets.
    # ------------------------------------------------------------------
    aliases = load_id_aliases()
    primary_gaia_id = get_target_from_schedule(obsdir, date, targname,
                                               aliases=aliases)
    if primary_gaia_id is None:
        # no plan file for tonight: resolve the name via the 40 pc list
        _by_name = get_target_from_target_list_by_name(targname, target_list_path)
        if _by_name is not None:
            primary_gaia_id = apply_id_alias(_by_name[0], targname, date, aliases)[0]
    if primary_gaia_id is None and 'toi' in targname.lower():
        primary_gaia_id = get_target_from_toi(targname, toi_table_path)

    if primary_gaia_id is not None:
        logger.info(
            "Primary target identified early for naming: '%s'", primary_gaia_id
        )
    else:
        logger.warning(
            "Could not identify primary target before stacking — "
            "output files will use target name '%s' as identifier", targname
        )

    # Construct final output filenames immediately.
    # Convention: {gaia_id}_{telescope}_{instrument}_{date}_{suffix}_{filter}.{ext}
    stack_name = get_output_filename(
        primary_gaia_id, telescope, instrument, date, filter, 'outstack', ext
    )
    conf_name = stack_name.replace(f'.{ext}', f'_conf.{ext}')
    cat_name = get_output_filename(
        primary_gaia_id, telescope, instrument, date, filter,
        'stack_catalogue', ext
    )

    stack_path = os.path.join(outdir, stack_name)
    conf_path = os.path.join(outdir, conf_name)
    cat_path = os.path.join(outdir, cat_name)

    # Write sidecar file recording the expected catalogue filename for
    # this night. T9 and condense read this to identify tonight's catalogue
    # and prefer any other catalogue in the directory (i.e. a backup).
    sidecar_path = os.path.join(outdir, '.expected_stack_catalogue')
    try:
        with open(sidecar_path, 'w') as _sc:
            _sc.write(cat_name + '\n')
        logger.info("Sidecar written: '%s' → %s", sidecar_path, cat_name)
    except Exception as _e:
        logger.warning("Could not write sidecar file: %s", _e)

    # ------------------------------------------------------------------
    # Step 2: image selection
    # ------------------------------------------------------------------
    selected_images = select_images(all_images)

    if not selected_images:
        logger.error("No images passed WCS header check — aborting")
        return 1

    # ------------------------------------------------------------------
    # Step 3: WCS solve selected images
    # ------------------------------------------------------------------
    tmp_name = stack_filelist + '.tmp'
    try:
        with open(tmp_name, 'w') as tmp:
            for image in selected_images:
                tmp.write(image + '\n')

        if not no_wcs:
            logger.info("WCS solving %d selected image(s)", len(selected_images))
            m_solve_images(
                tmp_name, tmp_name,
                thresh=s_thresh,
                nproc=ncores,
                verbose=verbose,
                rcore=rcore,
                catsrc=catsrc,
                ipix=ipix,
                ext=ext,
            )
        else:
            logger.info("Skipping WCS solving (no_wcs=True)")

        # ------------------------------------------------------------------
        # Step 4: filter out failed WCS solutions
        # ------------------------------------------------------------------
        passed = []
        for image in selected_images:
            if wcs_succeeded(image):
                passed.append(image)
            else:
                logger.warning(
                    "WCS solve failed for '%s' — excluding from stack", image
                )

    finally:
        if os.path.exists(tmp_name):
            os.remove(tmp_name)

    if not passed:
        logger.error("No images survived WCS filtering — aborting")
        return 1

    logger.info("%d image(s) will be included in the stack", len(passed))

    with open(stack_filelist, 'w') as f:
        f.write('\n'.join(passed))

    # ------------------------------------------------------------------
    # Step 5: stack images — written directly to final output paths
    # ------------------------------------------------------------------
    logger.info("Stacking images → %s", stack_path)
    casutools.imstack(
        stack_filelist,
        confidence_map=None,
        verbose=verbose,
        outstack=stack_path,
        outconf=conf_path,
    )

    # ------------------------------------------------------------------
    # Step 6: source detection
    # ------------------------------------------------------------------
    logger.info("Detecting sources in stack → %s", cat_path)
    casutools.imcore(
        stack_path,
        cat_path,
        threshold=c_thresh,
        confidence_map=conf_path,
        verbose=verbose,
        ellfile=True,
    )

    # ------------------------------------------------------------------
    # Step 7: Gaia DR3 crossmatch
    # ------------------------------------------------------------------
    logger.info("Crossmatching detected sources with Gaia DR3")
    perc, num_sources = gaia_dr2.crossmatch(
        cat_path,
        logfile=cat_path + '_Gaia.log',
        n=ncores,
        ext=ext,
        date=date,
        catsrc=catsrc,
    )

    # ------------------------------------------------------------------
    # Step 8: guard — abort if crossmatch produced nothing
    # ------------------------------------------------------------------
    if num_sources == 0:
        logger.error(
            "Gaia crossmatch produced no sources for '%s' — aborting",
            targname
        )
        return 1

    logger.info(
        "Crossmatch complete: %d source(s), %.1f%% matched with DR3",
        num_sources, perc
    )

    # ------------------------------------------------------------------
    # Step 9: full target identification (validates presence in catalogue,
    # identifies secondary targets, writes TARGET_ROLE column)
    # ------------------------------------------------------------------
    logger.info("Identifying targets in FOV for '%s'", targname)

    # Parse directory structure for QC reporting (needed throughout)
    dirname = os.path.dirname(filelist)
    dirsplit = dirname.split('/')
    run = dirsplit[-1]
    field = dirsplit[-2]
    dir_date = dirsplit[-3]
    tel = dirsplit[-5]

    try:
        cat_info = read_catalogue_info(cat_path)
    except Exception as e:
        logger.error(
            "Could not read Gaia IDs from catalogue '%s': %s", cat_path, e
        )
        return 1

    ident_info = {}
    targets = identify_targets(
        obsdir=obsdir,
        date=date,
        targname=targname,
        catalogue_dr2_ids=cat_info['dr2'],
        target_list_path=target_list_path,
        toi_table_path=toi_table_path,
        catalogue_dr3_ids=cat_info['dr3'],
        catalogue_coords=cat_info['coords'],
        catalogue_fluxes=cat_info['flux'],
        aliases=aliases,
        info=ident_info,
    )
    logger.info("Primary identification status: %s (method %s)",
                ident_info.get('status'), ident_info.get('primary_method'))

    # Coordinate fallback found the star but the crossmatch had no ID for
    # it: write the identity into the catalogue so every later stage
    # (TARGET_ROLE, condense, light curves) can key on it as usual.
    if ident_info.get('inject'):
        inj = ident_info['inject']
        inject_primary_identity(cat_path, inj['row'], inj['gaia_id'],
                                inj.get('db_row'))

    primary_targets = [t for t in targets if t[1] == 'primary']

    # ------------------------------------------------------------------
    # If primary target not found in tonight's catalogue, attempt to
    # restore a backup. The backup lands in the run directory with its
    # own filename (different date/telescope), so tonight's failed stack
    # and catalogue are preserved alongside it for diagnostic purposes.
    # The backup catalogue is NOT copied to StackImages — only successful
    # nightly stacks are written there.
    # ------------------------------------------------------------------
    using_backup = False
    active_cat_path = cat_path
    active_stack_path = stack_path

    if not primary_targets and primary_gaia_id is not None:
        logger.warning("=" * 60)
        logger.warning(
            "PRIMARY TARGET %s NOT IN TONIGHT'S CATALOGUE — "
            "searching for backup", primary_gaia_id
        )
        logger.warning("=" * 60)

        backup_stack, backup_cat = find_backup(
            backupcatdir, primary_gaia_id, instrument, date
        )

        if backup_stack is not None:
            # Restore backup into run directory preserving its original
            # filename (which encodes its own date/telescope). Tonight's
            # failed files remain untouched.
            backup_stack_dest = os.path.join(
                outdir, os.path.basename(backup_stack)
            )
            backup_cat_dest = os.path.join(
                outdir, os.path.basename(backup_cat)
            )
            try:
                import shutil as _shutil
                _shutil.copy2(backup_stack, backup_stack_dest)
                _shutil.copy2(backup_cat, backup_cat_dest)
                logger.info(
                    "Backup catalogue restored to run directory: %s",
                    os.path.basename(backup_cat_dest)
                )
                active_cat_path = backup_cat_dest
                active_stack_path = backup_stack_dest
                using_backup = True
            except Exception as e:
                logger.error(
                    "Failed to copy backup to run directory: %s", e
                )
                logger.error(
                    "No valid catalogue available for '%s' — aborting",
                    targname
                )
                return 1

            # Re-identify targets using the restored catalogue
            try:
                backup_info = read_catalogue_info(active_cat_path)
            except Exception as e:
                logger.error(
                    "Could not read Gaia IDs from backup catalogue '%s': %s",
                    active_cat_path, e
                )
                return 1
            # No coordinate fallback on a backup: its positions belong to
            # another night's epoch.
            ident_info = {}
            targets = identify_targets(
                obsdir=obsdir,
                date=date,
                targname=targname,
                catalogue_dr2_ids=backup_info['dr2'],
                target_list_path=target_list_path,
                toi_table_path=toi_table_path,
                catalogue_dr3_ids=backup_info['dr3'],
                aliases=aliases,
                info=ident_info,
            )
            primary_targets = [t for t in targets if t[1] == 'primary']

            if not primary_targets:
                logger.error(
                    "Primary target %s still not found in backup catalogue — "
                    "aborting", primary_gaia_id
                )
                return 1

            logger.info(
                "Primary target identified in backup catalogue — "
                "proceeding with backup"
            )

        else:
            logger.error(
                "No backup available for primary target %s — aborting",
                primary_gaia_id
            )
            return 1

    secondary_targets = [t for t in targets if t[1] == 'secondary']

    if not targets:
        logger.error(
            "NO TARGETS IDENTIFIED IN FOV FOR '%s' — cannot produce "
            "any light curves. Aborting.", targname
        )
        return 1

    # ------------------------------------------------------------------
    # Step 10: write TARGET_ROLE column into the active catalogue
    # ------------------------------------------------------------------
    write_target_roles(active_cat_path, targets)
    record_primary_provenance(active_cat_path, ident_info)

    # ------------------------------------------------------------------
    # Step 11: backup management
    # Only copy to StackImages if tonight's stack was successful
    # (primary target crossmatched). Restored backups are never
    # re-written to StackImages.
    # ------------------------------------------------------------------
    if not using_backup:
        logger.info(
            "Tonight's stack valid for primary target — copying to backup"
        )
        update_backup(backupcatdir, active_stack_path, active_cat_path)
    else:
        logger.info(
            "Using backup catalogue — tonight's stack not copied to StackImages"
        )

    # QC reporting
    try:
        QC([0, reportdir + '/QC', dir_date, field, tel, run, 'QC11', perc])
        QC([0, reportdir + '/QC', dir_date, field, tel, run, 'QC10', num_sources])
    except Exception as e:
        logger.warning("Could not update QC report: %s", e)

    if verbose:
        logger.info("catalogue_fov complete for '%s'", targname)

    return 0


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def create_cat():
    parser = argparse.ArgumentParser(
        description='Create FOV stack image and source catalogue (Stage T8)'
    )
    parser.add_argument('filelist')
    parser.add_argument('outdir')
    parser.add_argument('stackimagesdir',
                        help='Shared StackImages backup directory '
                             '(PipelineOutput/vN/StackImages)')
    parser.add_argument('reportdir')
    parser.add_argument('filter')
    parser.add_argument('date')
    parser.add_argument('obsdir',
                        help='Observatory base directory (for schedule lookup)')
    parser.add_argument('targname',
                        help='Target name as passed by scheduler')
    parser.add_argument('telescope',
                        help='Telescope name e.g. Ganymede')
    parser.add_argument('target_list_path',
                        help='Path to 40pc target list')
    parser.add_argument('toi_table_path',
                        help='Path to TOI/Gaia ID lookup table')
    parser.add_argument('no_wcs')
    parser.add_argument('c_thresh')
    parser.add_argument('s_thresh')
    parser.add_argument('verbose')
    parser.add_argument('ipix')
    parser.add_argument('catsrc')
    parser.add_argument('rcore')
    parser.add_argument('ncores')
    parser.add_argument('ext')
    args = parser.parse_args()

    exit_code = main(
        filelist=args.filelist,
        outdir=args.outdir,
        backupcatdir=args.stackimagesdir,
        reportdir=args.reportdir,
        filter=args.filter,
        date=args.date,
        obsdir=args.obsdir,
        targname=args.targname,
        telescope=args.telescope,
        target_list_path=args.target_list_path,
        toi_table_path=args.toi_table_path,
        no_wcs=args.no_wcs == 'True',
        c_thresh=args.c_thresh,
        s_thresh=args.s_thresh,
        verbose=args.verbose == 'True',
        ipix=args.ipix,
        catsrc=args.catsrc,
        rcore=args.rcore,
        ncores=args.ncores,
        ext=args.ext,
    )

    sys.exit(exit_code)


if __name__ == '__main__':
    create_cat()