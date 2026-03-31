#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Target Management

Provides target identification, catalogue validation, backup management,
and output file naming for use by catalogue_fov.py and condense_photometry.py.

Target identification follows a strict hierarchy:
    1. Schedule plan file (primary target, authoritative)
    2. TOI lookup table (primary target, fallback if schedule fails)
    3. 40pc target list crossmatch (secondary targets, always runs)

No side effects on import. All state is passed explicitly via arguments.
"""

import re
import logging
import shutil
import glob
import os
from pathlib import Path
from datetime import datetime

import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table

from utils import gaia_id_from_schedule

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Target identification
# ---------------------------------------------------------------------------

def get_target_from_schedule(obsdir, date, targname):
    """
    Attempt to identify the primary target's Gaia DR2 ID from the
    observatory schedule plan file.

    Inputs:
        obsdir   : str  Observatory base directory
        date     : str  Observation date in YYYYMMDD format
        targname : str  Target name as passed by the scheduler

    Output:
        str  Gaia DR2 ID if found, else None
    """
    targname_clean = targname.replace('--', ' ')

    plan_files = gaia_id_from_schedule.find_plan(obsdir, date, targname_clean)

    if not plan_files:
        logger.warning(
            "No plan file found for target '%s' on date %s", targname, date
        )
        return None

    dr2_id = gaia_id_from_schedule.read_file(plan_files[0], silent=True)

    if dr2_id is None:
        logger.warning(
            "Plan file found for '%s' but contains no Gaia DR2 ID", targname
        )
        return None

    logger.info(
        "Gaia DR2 ID '%s' extracted from schedule plan for '%s'",
        dr2_id, targname
    )
    return dr2_id.strip()


def get_target_from_toi(targname, toi_table_path):
    """
    Attempt to identify the primary target's Gaia DR2 ID from the TOI
    lookup table, using the TOI number parsed from the target name.
    Planet designators (e.g. '.01', 'b') are stripped automatically —
    the Gaia ID is for the host star, not the planet.

    Inputs:
        targname       : str  Target name (expected to contain 'TOI')
        toi_table_path : str  Path to the TOI/Gaia ID lookup CSV

    Output:
        str  Gaia DR2 ID if found, else None
    """
    match = re.search(r'toi[-_\s]?(\d+)', targname, re.IGNORECASE)

    if match is None:
        logger.warning(
            "Could not parse TOI number from target name '%s'", targname
        )
        return None

    toi_no = int(match.group(1))
    logger.info("Parsed TOI number %d from target name '%s'", toi_no, targname)

    try:
        toi_table = Table.read(toi_table_path, format='ascii.csv')
    except Exception as e:
        logger.error("Failed to read TOI table at '%s': %s", toi_table_path, e)
        return None

    row_idx = np.where(toi_table['TOI'] == toi_no)[0]

    if len(row_idx) == 0:
        logger.warning("TOI %d not found in TOI table", toi_no)
        return None

    dr2_id = str(int(toi_table['GAIA'][row_idx[0]]))
    logger.info("TOI %d resolved to Gaia DR2 ID '%s'", toi_no, dr2_id)
    return dr2_id


def get_targets_from_target_list(catalogue_dr2_ids, target_list_path):
    """
    Crossmatch Gaia DR2 IDs from the FOV catalogue against the 40pc
    target list. SP names are ignored entirely — matching is by DR2 ID only.

    Inputs:
        catalogue_dr2_ids : list of str  Cleaned DR2 IDs from Gaia_Crossmatch
        target_list_path  : str          Path to the 40pc target list

    Output:
        list of (DR2_ID: str, Teff: int or None) tuples
        Empty list if no matches or on read failure.
    """
    try:
        data = ascii.read(target_list_path, delimiter=' ',
                          header_start=0, data_start=1)
    except Exception as e:
        logger.error(
            "Failed to read target list at '%s': %s", target_list_path, e
        )
        return []

    try:
        target_list_ids = [str(x).strip() for x in data['Gaia_ID,']]
    except KeyError:
        logger.error(
            "Column 'Gaia_ID,' not found in target list '%s'", target_list_path
        )
        return []

    catalogue_id_set = set(catalogue_dr2_ids)
    target_list_id_set = set(target_list_ids)
    intersection = catalogue_id_set & target_list_id_set

    results = []
    for dr2_id in intersection:
        row_idx = target_list_ids.index(dr2_id)
        try:
            teff_val = data['T_eff,'][row_idx]
            teff = int(teff_val) if not np.ma.is_masked(teff_val) else None
        except (KeyError, ValueError):
            teff = None
        results.append((dr2_id, teff))

    logger.info(
        "Target list crossmatch found %d target(s) in FOV", len(results)
    )
    return results


def identify_targets(obsdir, date, targname, catalogue_dr2_ids,
                     target_list_path, toi_table_path):
    """
    Orchestrate the full target identification process for a given FOV.

    Identification hierarchy:
        1. Schedule plan file  → primary target
        2. TOI lookup table    → primary target (fallback if schedule fails)
        3. 40pc target list    → secondary targets (always runs)

    Inputs:
        obsdir            : str       Observatory base directory
        date              : str       Observation date YYYYMMDD
        targname          : str       Target name as passed by scheduler
        catalogue_dr2_ids : list[str] Raw DR2 IDs from Gaia_Crossmatch extension
        target_list_path  : str       Path to 40pc target list
        toi_table_path    : str       Path to TOI/Gaia ID lookup table

    Output:
        list of (DR2_ID: str, role: str, Teff: int or None)
        role is 'primary' or 'secondary'
        Empty list if catalogue contains no usable DR2 IDs.
    """
    # Guard: empty catalogue
    clean_ids = [x.strip() for x in catalogue_dr2_ids
                 if x.strip() not in ('', 'nan')]

    if not clean_ids:
        logger.error(
            "Catalogue contains no usable DR2 IDs — cannot identify targets"
        )
        return []

    match_gaia = []   # list of (DR2_ID, role, Teff)
    primary_found = False

    # ------------------------------------------------------------------
    # Step 1: schedule
    # ------------------------------------------------------------------
    dr2 = get_target_from_schedule(obsdir, date, targname)

    if dr2 is not None:
        if dr2 in clean_ids:
            match_gaia.append((dr2, 'primary', None))
            primary_found = True
            logger.info(
                "Primary target identified from schedule: %s", dr2
            )
        else:
            logger.warning(
                "Schedule DR2 ID '%s' not found in FOV catalogue — "
                "target may not have been observed", dr2
            )

    # ------------------------------------------------------------------
    # Step 2: TOI fallback (only if schedule failed)
    # ------------------------------------------------------------------
    if not primary_found:
        logger.warning(
            "Primary target could not be identified from schedule for '%s'",
            targname
        )

        if 'toi' in targname.lower():
            dr2 = get_target_from_toi(targname, toi_table_path)

            if dr2 is not None:
                if dr2 in clean_ids:
                    match_gaia.append((dr2, 'primary', None))
                    primary_found = True
                    logger.info(
                        "Primary target identified from TOI table: %s", dr2
                    )
                else:
                    logger.warning(
                        "TOI DR2 ID '%s' not found in FOV catalogue", dr2
                    )
            else:
                logger.warning(
                    "Could not identify primary target from TOI table "
                    "for '%s'", targname
                )
        else:
            logger.warning(
                "Target '%s' is not a TOI — no further primary "
                "identification possible", targname
            )

    if not primary_found:
        logger.warning(
            "PRIMARY TARGET COULD NOT BE IDENTIFIED BY ANY METHOD for '%s'",
            targname
        )

    # ------------------------------------------------------------------
    # Step 3: target list crossmatch (always runs)
    # ------------------------------------------------------------------
    additional = get_targets_from_target_list(clean_ids, target_list_path)

    existing_ids = [x[0] for x in match_gaia]
    for (dr2_id, teff) in additional:
        if dr2_id not in existing_ids:
            match_gaia.append((dr2_id, 'secondary', teff))
            existing_ids.append(dr2_id)

    # Backfill Teff for primary target if it appears in the target list
    for i, (dr2_id, role, teff) in enumerate(match_gaia):
        if role == 'primary' and teff is None:
            tlist_match = [t for (d, t) in additional if d == dr2_id]
            if tlist_match:
                match_gaia[i] = (dr2_id, role, tlist_match[0])

    logger.info(
        "Identified %d target(s) in FOV in total (%d primary, %d secondary)",
        len(match_gaia),
        sum(1 for x in match_gaia if x[1] == 'primary'),
        sum(1 for x in match_gaia if x[1] == 'secondary')
    )

    return match_gaia


# ---------------------------------------------------------------------------
# Catalogue validation
# ---------------------------------------------------------------------------

def catalogue_is_valid(catalogue_path, target_dr2_id):
    """
    Check whether a stack catalogue contains a valid Gaia crossmatch
    entry for a specific target.

    A catalogue is considered valid for a target if:
        - The Gaia_Crossmatch FITS extension exists
        - The target's DR2 ID appears in the GAIA_DR2_ID column
          (excluding 'nan' entries)

    Inputs:
        catalogue_path : str  Path to the stack catalogue FITS file
        target_dr2_id  : str  Gaia DR2 ID of the target to check for

    Output:
        bool
    """
    target_dr2_id = target_dr2_id.strip()

    try:
        with fits.open(catalogue_path) as hdul:
            ext_names = [hdu.name.upper() for hdu in hdul]

            if 'GAIA_CROSSMATCH' not in ext_names:
                logger.warning(
                    "No Gaia_Crossmatch extension in '%s'", catalogue_path
                )
                return False

            ids = hdul['GAIA_CROSSMATCH'].data['GAIA_DR2_ID']

    except Exception as e:
        logger.warning(
            "Could not open catalogue '%s': %s", catalogue_path, e
        )
        return False

    clean_ids = [x.strip() for x in ids if x.strip() not in ('', 'nan')]

    if target_dr2_id in clean_ids:
        logger.info(
            "Target '%s' found in catalogue '%s'",
            target_dr2_id, catalogue_path
        )
        return True

    logger.warning(
        "Target '%s' not found in catalogue '%s'",
        target_dr2_id, catalogue_path
    )
    return False


# ---------------------------------------------------------------------------
# Backup management
# ---------------------------------------------------------------------------

def _parse_date_from_filename(filename):
    """
    Extract the 8-digit observation date from a stack catalogue filename
    of the form {gaia_id}_{telescope}_{instrument}_{date}_stack_catalogue_{filter}.fits
    Returns None if no date can be parsed.
    """
    m = re.search(r'_(\d{8})_stack_catalogue_', filename)
    if m:
        return m.group(1)
    return None


def find_backup(backup_dir, primary_gaia_id, instrument, obs_date):
    """
    Search for the best backup stack catalogue for a given primary target
    and instrument, selecting the one whose date is closest to obs_date.

    The search key is (primary_gaia_id, instrument) — telescope and filter
    are ignored for matching. Any candidate catalogue for this primary target
    and instrument is valid regardless of which telescope or filter produced it.

    Inputs:
        backup_dir       : str  Path to shared StackImages backup directory
        primary_gaia_id  : str  Gaia DR2 ID of the primary (scheduled) target
        instrument       : str  Instrument name e.g. 'andor', 'spirit'
        obs_date         : str  Observation date YYYYMMDD

    Output:
        (stack_path: str, cat_path: str) if found, else (None, None)
    """
    backup_dir_path = Path(backup_dir)

    # Glob for all catalogues matching (primary_gaia_id, instrument)
    pattern = f'{primary_gaia_id}_*_{instrument}_*_stack_catalogue_*.fits'
    cat_candidates = list(backup_dir_path.glob(pattern))

    if not cat_candidates:
        logger.warning(
            "No backup catalogue found for primary target '%s', "
            "instrument '%s' in '%s'",
            primary_gaia_id, instrument, backup_dir
        )
        return None, None

    # Parse dates and find closest to obs_date
    try:
        obs_dt = datetime.strptime(obs_date, '%Y%m%d')
    except ValueError:
        logger.warning("Could not parse obs_date '%s' as YYYYMMDD", obs_date)
        return None, None

    best_cat = None
    best_dist = None

    for cat_path in cat_candidates:
        date_str = _parse_date_from_filename(cat_path.name)
        if date_str is None:
            continue
        try:
            cat_dt = datetime.strptime(date_str, '%Y%m%d')
        except ValueError:
            continue
        dist = abs((obs_dt - cat_dt).days)
        if best_dist is None or dist < best_dist:
            # Check corresponding stack exists
            stack_name = cat_path.name.replace('stack_catalogue', 'outstack')
            stack_path = backup_dir_path / stack_name
            if stack_path.exists():
                best_dist = dist
                best_cat = cat_path

    if best_cat is None:
        logger.warning(
            "No complete backup pair (stack + catalogue) found for "
            "primary target '%s', instrument '%s'",
            primary_gaia_id, instrument
        )
        return None, None

    best_stack = backup_dir_path / best_cat.name.replace(
        'stack_catalogue', 'outstack'
    )
    logger.info(
        "Backup found for primary target '%s', instrument '%s': "
        "date %s (distance %d days from %s)",
        primary_gaia_id, instrument,
        _parse_date_from_filename(best_cat.name),
        best_dist, obs_date
    )
    return str(best_stack), str(best_cat)


def restore_backup(backup_dir, primary_gaia_id, instrument, obs_date,
                   outstack_path, outcat_path):
    """
    Restore the closest backup stack image and catalogue for a given
    primary target and instrument.

    Inputs:
        backup_dir       : str  Path to shared StackImages backup directory
        primary_gaia_id  : str  Gaia DR2 ID of the primary target
        instrument       : str  Instrument name
        obs_date         : str  Observation date YYYYMMDD
        outstack_path    : str  Destination path for stack image
        outcat_path      : str  Destination path for catalogue

    Output:
        bool  True if restoration successful, False otherwise
    """
    stack_path, cat_path = find_backup(
        backup_dir, primary_gaia_id, instrument, obs_date
    )

    if stack_path is None:
        logger.warning(
            "No backup available for primary target '%s', "
            "instrument '%s' — cannot restore",
            primary_gaia_id, instrument
        )
        return False

    try:
        shutil.copy2(stack_path, outstack_path)
        shutil.copy2(cat_path, outcat_path)
        logger.info(
            "Backup restored for primary target '%s' from '%s'",
            primary_gaia_id, stack_path
        )
        return True

    except Exception as e:
        logger.error(
            "Failed to restore backup for '%s': %s", primary_gaia_id, e
        )
        return False


def update_backup(backup_dir, outstack_path, outcat_path):
    """
    Copy the current stack image and catalogue to the shared StackImages
    backup directory. The filename already encodes the primary Gaia ID,
    telescope, instrument, date and filter, so no additional metadata is
    needed. Each night's successful stack is preserved independently —
    filenames are unique per date so nothing is overwritten.

    Inputs:
        backup_dir    : str  Path to shared StackImages backup directory
        outstack_path : str  Source stack image path
        outcat_path   : str  Source catalogue path
    """
    backup_dir_path = Path(backup_dir)
    backup_dir_path.mkdir(parents=True, exist_ok=True)

    dest_stack = backup_dir_path / Path(outstack_path).name
    dest_cat = backup_dir_path / Path(outcat_path).name

    try:
        shutil.copy2(outstack_path, dest_stack)
        shutil.copy2(outcat_path, dest_cat)
        logger.info(
            "Backup written: '%s' and '%s'",
            dest_stack.name, dest_cat.name
        )
    except Exception as e:
        logger.error(
            "Failed to copy files to backup directory: %s", e
        )


# ---------------------------------------------------------------------------
# Output file naming
# ---------------------------------------------------------------------------

def get_output_filename(primary_gaia_id, telescope, instrument, date,
                        filter, suffix, ext):
    """
    Construct a stack output filename encoding all parameters needed
    for the backup search system.

    Convention: {gaia_id}_{telescope}_{instrument}_{date}_{suffix}_{filter}.{ext}

    Matching key: (gaia_id, instrument) — telescope, date and filter
    provide uniqueness and are used for date-proximity search.

    Inputs:
        primary_gaia_id : str or None  Gaia DR2 ID of primary target
        telescope       : str          Telescope name e.g. 'Ganymede'
        instrument      : str          Instrument name e.g. 'andor'
        date            : str          Observation date YYYYMMDD
        filter          : str          Filter name e.g. 'r', 'I+z'
        suffix          : str          e.g. 'outstack', 'stack_catalogue'
        ext             : str          e.g. 'fits'

    Output:
        str  Filename (not a full path)
    """
    if primary_gaia_id is not None:
        return (f"{primary_gaia_id}_{telescope}_{instrument}_{date}"
                f"_{suffix}_{filter}.{ext}")

    logger.warning(
        "Primary DR2 ID unavailable — output files will lack Gaia ID "
        "in filename and will not be found by the backup system"
    )
    return f"unknown_{telescope}_{instrument}_{date}_{suffix}_{filter}.{ext}"