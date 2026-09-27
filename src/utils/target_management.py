#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Target Management

Provides target identification, catalogue validation, backup management,
and output file naming for use by catalogue_fov.py and condense_photometry.py.

Target identification follows a strict hierarchy:
    1. Schedule plan file (primary target, authoritative), with an alias
       table applied for plan files known to carry a wrong Gaia ID
    1b. Coordinate fallback: if the schedule ID is not in the field
       catalogue, the plan-file position is propagated to the night with
       the star's proper motion and the nearest detection is taken as
       primary, with the schedule ID injected into the catalogue
    2. TOI lookup table (primary target, fallback if schedule fails)
    3. 40pc target list crossmatch (secondary targets, always runs)

IDs from the schedule, the TOI table and the target list may be Gaia DR2
or DR3 source IDs; catalogue rows are matched on either column.  When the
schedule ID matches nothing, the local Gaia database is consulted so the
log distinguishes "target not in our database" (a build problem, see
PRIMARY_NOT_IN_DB) from "target not in tonight's field" (PRIMARY_NOT_IN_FOV).

No side effects on import. All state is passed explicitly via arguments.
"""

import re
import csv
import math
import sqlite3
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
from utils.target_list import clean_id as _clean_target_id, load_target_list

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Target identification
# ---------------------------------------------------------------------------

DEFAULT_ALIAS_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'calibration', 'target_id_aliases.csv')

DEFAULT_DB_PATH = '/gaia_database/gaia_dr3_unified_16jcut.db'

PLAN_COORD_EPOCH = 2000.0   # plan-file coordinates are J2000
GAIA_DR3_EPOCH = 2016.0


def _norm_name(name):
    return (name or '').replace('--', ' ').replace('_', ' ').strip().upper()


def load_id_aliases(path=None):
    """
    Read calibration/target_id_aliases.csv: rows of
    target, schedule_id, use_id, from_date, to_date, note.
    Lines starting with '#' are comments.  Returns a list of dicts;
    an unreadable file yields an empty list (aliases are optional).
    """
    path = path or DEFAULT_ALIAS_PATH
    aliases = []
    try:
        with open(path, newline='') as f:
            lines = [l for l in f if not l.lstrip().startswith('#') and l.strip()]
        for row in csv.DictReader(lines, fieldnames=[
                'target', 'schedule_id', 'use_id', 'from_date', 'to_date', 'note']):
            if not row.get('schedule_id') or not row.get('use_id'):
                continue
            aliases.append({k: (v or '').strip() for k, v in row.items()})
    except OSError:
        return []
    return aliases


def apply_id_alias(gaia_id, targname, date, aliases):
    """
    Replace gaia_id by the alias table's use_id when a row matches on
    schedule_id, on target name (blank = any) and on date range
    (blank bounds = unbounded).  Returns (id, alias_row_or_None).
    """
    if gaia_id is None or not aliases:
        return gaia_id, None
    for a in aliases:
        if a['schedule_id'] != str(gaia_id).strip():
            continue
        if a['target'] and _norm_name(a['target']) != _norm_name(targname):
            continue
        if a['from_date'] and date < a['from_date']:
            continue
        if a['to_date'] and date > a['to_date']:
            continue
        logger.warning(
            "ID ALIAS APPLIED for '%s' on %s: schedule says %s, using %s (%s)",
            targname, date, gaia_id, a['use_id'], a.get('note', ''))
        return a['use_id'], a
    return gaia_id, None


def get_coords_from_schedule(obsdir, date, targname):
    """
    Return the (ra_deg, dec_deg) J2000 position from the plan file for
    this target and night, or None.
    """
    targname_clean = targname.replace('--', ' ')
    plan_files = gaia_id_from_schedule.find_plan(obsdir, date, targname_clean)
    if not plan_files:
        return None
    return gaia_id_from_schedule.read_coords(plan_files[0])


def lookup_gaia_id_in_db(gaia_id, dec_hint=None, db_path=None):
    """
    Look a Gaia ID up directly in the local database (band table chosen
    from dec_hint, plus the hand-curated custom_sources table).  Matches
    on source_id (DR3) or dr2_source_id.  Returns a dict of the row, or
    None if absent.  Without dec_hint only custom_sources is searched,
    because an unsharded scan of the database is far too slow.
    """
    if gaia_id is None:
        return None
    db_path = db_path or os.environ.get('GAIADATABASEPATH', DEFAULT_DB_PATH)
    if not os.path.exists(db_path):
        logger.warning("Gaia database '%s' not found — cannot check ID %s",
                       db_path, gaia_id)
        return None
    gid = str(gaia_id).strip()
    tables = ['custom_sources']
    if dec_hint is not None and math.isfinite(dec_hint):
        lo = int(math.floor(dec_hint))
        lo = max(-90, min(89, lo))
        tables.insert(0, '%d_%d' % (lo, lo + 1))
    cols = ('ra', 'dec', 'pmra', 'pmdec', 'phot_g_mean_mag', 'g_rp', 'bp_rp',
            'parallax', 'teff_gspphot', 'source_id', 'dr2_source_id')
    try:
        conn = sqlite3.connect(db_path, timeout=600)
    except sqlite3.Error as e:
        logger.warning("Could not open Gaia database '%s': %s", db_path, e)
        return None
    try:
        for t in tables:
            try:
                cur = conn.execute(
                    "SELECT %s FROM '%s' WHERE source_id = ? OR dr2_source_id = ? LIMIT 1"
                    % (', '.join(cols), t), (gid, gid))
                row = cur.fetchone()
            except sqlite3.OperationalError:
                continue   # table absent (old database) — try the next one
            if row is not None:
                d = dict(zip(cols, row))
                d['table'] = t
                return d
    finally:
        conn.close()
    return None


def _clean_id(value):
    """
    Normalise a Gaia ID, or '' when it is absent.

    Delegates to utils.target_list so that every absent spelling is handled
    in one place — including '--', which is what astropy prints for a
    masked cell and therefore what an empty field in a CSV integer column
    becomes.
    """
    return _clean_target_id(value)


def _preferred_id(gaia_ids, catalogue_ids=None):
    """
    Which of a target-list row's identifiers to carry forward.

    A row may hold both a DR2 and a DR3 identifier.  Prefer one the field
    catalogue actually contains, so that the later crossmatch succeeds;
    failing that take the first, which is DR2 before DR3 and so keeps
    output filenames what they have always been for a star whose DR2
    identifier still exists.
    """
    if not gaia_ids:
        return ''
    if catalogue_ids:
        for gid in gaia_ids:
            if gid in catalogue_ids:
                return gid
    return gaia_ids[0]


def _propagate(ra, dec, pmra, pmdec, from_epoch, to_epoch):
    """Linear proper-motion propagation in degrees; pm in mas/yr."""
    if pmra is None or pmdec is None or not (math.isfinite(pmra) and math.isfinite(pmdec)):
        return ra, dec
    dt = to_epoch - from_epoch
    return (ra + pmra * 1e-3 * dt / 3600.0 / math.cos(math.radians(dec)),
            dec + pmdec * 1e-3 * dt / 3600.0)


def _sep_arcsec(ra1, dec1, ra2, dec2):
    return 3600.0 * math.hypot((ra1 - ra2) * math.cos(math.radians(dec1)), dec1 - dec2)


def _date_to_epoch(date):
    d = datetime.strptime(date, '%Y%m%d')
    return d.year + (d.timetuple().tm_yday - 1) / 365.25


def get_target_from_schedule(obsdir, date, targname, aliases=None):
    """
    Attempt to identify the primary target's Gaia ID from the
    observatory schedule plan file, then apply the ID alias table.

    Inputs:
        obsdir   : str  Observatory base directory
        date     : str  Observation date in YYYYMMDD format
        targname : str  Target name as passed by the scheduler
        aliases  : list from load_id_aliases(); None loads the default file

    Output:
        str  Gaia ID (DR2 or DR3, as written by the scheduler) if found, else None
    """
    if aliases is None:
        aliases = load_id_aliases()
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
        "Gaia ID '%s' extracted from schedule plan for '%s'",
        dr2_id, targname
    )
    dr2_id, _ = apply_id_alias(dr2_id.strip(), targname, date, aliases)
    return dr2_id


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


def get_target_from_target_list_by_name(targname, target_list_path,
                                        catalogue_ids=None):
    """
    Resolve a target NAME (the scheduler's Sp_ID) to a Gaia ID and a J2000
    position via the 40 pc target list.  Used when there is no schedule
    plan file for the night, which is the case for roughly half of the
    observed nights of some targets.

    Duplicated Sp_IDs (resolved binaries sharing a coordinate-derived name)
    are disambiguated by preferring a row whose ID is in `catalogue_ids`;
    otherwise the first row is used and a warning is logged.

    Output: (gaia_id: str, ra_deg, dec_deg, teff) with None for unknown
            values, or None if the name is absent or has no usable ID.
    """
    tlist = load_target_list(target_list_path, logger)
    if tlist is None or not tlist.has('Sp_ID'):
        return None
    want = _norm_name(targname)
    rows = [i for i in range(len(tlist)) if _norm_name(tlist.sp_id(i)) == want]
    # A row with no usable ID cannot become a target, whatever its name.
    rows = [i for i in rows if tlist.ids(i)]
    if not rows:
        return None
    if len(rows) > 1:
        in_fov = [i for i in rows if catalogue_ids and
                  any(g in catalogue_ids for g in tlist.ids(i))]
        if len(in_fov) == 1:
            rows = in_fov
        else:
            logger.warning(
                "Target name '%s' matches %d target-list rows with different "
                "IDs; using the first (%s)", targname, len(rows),
                tlist.ids(rows[0])[0])
    i = rows[0]
    gaia_id = _preferred_id(tlist.ids(i), catalogue_ids)
    ra, dec = tlist.coords(i)
    teff = tlist.teff(i)
    logger.info("Target '%s' resolved via target list to Gaia ID %s", targname, gaia_id)
    return gaia_id, ra, dec, teff


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
    tlist = load_target_list(target_list_path, logger)
    if tlist is None:
        return []

    catalogue = set(_clean_id(x) for x in catalogue_dr2_ids)
    catalogue.discard('')

    # One result per target-list *row*, not per matching ID: a row carrying
    # both a DR2 and a DR3 identifier would otherwise match twice through
    # the two catalogue ID columns and be observed as two separate targets.
    results = []
    for row in range(len(tlist)):
        hits = [g for g in tlist.ids(row) if g in catalogue]
        if hits:
            # Report the ID the catalogue matched on, preferring DR2 so
            # that output filenames stay what they have always been.
            results.append((hits[0], tlist.teff(row)))

    logger.info(
        "Target list crossmatch found %d target(s) in FOV", len(results)
    )
    return results


def _coordinate_fallback(gaia_id, targname, date, plan_coords, dbrow,
                         cat_coords, cat_fluxes, dr2_clean, dr3_clean,
                         radius, info):
    """
    Locate the primary target by position when its ID matched no catalogue
    row.  The expected position is the database position (J2016) plus
    proper motion when the star is in the database, else the plan-file
    J2000 position with no proper motion (and a warning, because a fast
    star will then be missed).  Accepts the nearest detection within
    `radius` arcsec provided it carries no conflicting Gaia ID.
    Returns the catalogue row index, or None.
    """
    epoch = _date_to_epoch(date)
    if dbrow is not None and dbrow.get('ra') is not None:
        ra0, dec0, e0 = float(dbrow['ra']), float(dbrow['dec']), GAIA_DR3_EPOCH
        pmra, pmdec = dbrow.get('pmra'), dbrow.get('pmdec')
        src = 'database position (J2016) + proper motion'
    elif plan_coords is not None:
        ra0, dec0, e0 = plan_coords[0], plan_coords[1], PLAN_COORD_EPOCH
        pmra, pmdec = None, None
        src = 'plan-file position (J2000)'
    else:
        logger.warning(
            "Coordinate fallback impossible for '%s': no plan-file "
            "coordinates and not in database", targname)
        return None
    pm_known = (pmra is not None and pmdec is not None
                and math.isfinite(float(pmra)) and math.isfinite(float(pmdec)))
    if not pm_known:
        logger.warning(
            "Coordinate fallback for '%s' has NO proper motion: a "
            "fast-moving star may be missed or mismatched", targname)
    ra_exp, dec_exp = _propagate(ra0, dec0,
                                 float(pmra) if pm_known else None,
                                 float(pmdec) if pm_known else None,
                                 e0, epoch)
    coords = np.asarray(cat_coords, dtype=float)
    seps = np.array([_sep_arcsec(r, d, ra_exp, dec_exp) for r, d in coords])
    if not len(seps) or not np.isfinite(seps).any():
        return None
    i = int(np.nanargmin(seps))
    if not seps[i] <= radius:
        logger.warning(
            "COORD_NO_MATCH for '%s': nearest detection to expected "
            "position (%.5f, %.5f) is %.1f arcsec away (limit %.1f)",
            targname, ra_exp, dec_exp, seps[i], radius)
        info['status'] += '+COORD_NO_MATCH'
        return None
    existing = dr3_clean[i] or dr2_clean[i]
    if existing and existing != str(gaia_id).strip():
        logger.error(
            "COORD_CONFLICT for '%s': detection at the expected position "
            "already carries Gaia ID %s, not %s — refusing to override",
            targname, existing, gaia_id)
        info['status'] += '+COORD_CONFLICT'
        return None
    rank = ''
    if cat_fluxes is not None and len(cat_fluxes) == len(coords):
        f = np.nan_to_num(np.asarray(cat_fluxes, dtype=float), nan=-1.0)
        rank = ' (brightness rank %d of %d)' % (
            int(np.sum(f > f[i])) + 1, len(f))
    logger.warning(
        "PRIMARY IDENTIFIED BY COORDINATES for '%s': row %d, %.2f arcsec "
        "from expected position, using %s%s; injecting ID %s",
        targname, i, seps[i], src, rank, gaia_id)
    info['inject'] = {'row': i, 'gaia_id': str(gaia_id).strip(),
                      'db_row': dbrow, 'separation_arcsec': float(seps[i])}
    return i


def identify_targets(obsdir, date, targname, catalogue_dr2_ids,
                     target_list_path, toi_table_path,
                     catalogue_dr3_ids=None, catalogue_coords=None,
                     catalogue_fluxes=None, aliases=None, db_path=None,
                     info=None, coord_match_radius=2.0):
    """
    Orchestrate the full target identification process for a given FOV.

    Identification hierarchy:
        1.  Schedule plan file (with ID aliases)      → primary target
        1a. Target name looked up in the 40 pc list if there is no plan
            file (ID and J2000 position)               → primary target
        1b. Coordinate fallback if the schedule ID matched no row and
            catalogue_coords were supplied              → primary target
        2.  TOI lookup table (if schedule failed)       → primary target
        3.  40pc target list                            → secondary targets

    Inputs:
        obsdir            : str       Observatory base directory
        date              : str       Observation date YYYYMMDD
        targname          : str       Target name as passed by scheduler
        catalogue_dr2_ids : list[str] GAIA_DR2_ID column of Gaia_Crossmatch
        target_list_path  : str       Path to 40pc target list
        toi_table_path    : str       Path to TOI/Gaia ID lookup table
        catalogue_dr3_ids : list[str] GAIA_DR3_ID column (same length), optional
        catalogue_coords  : (N,2) array of RA, Dec in degrees, optional;
                            enables the coordinate fallback
        catalogue_fluxes  : (N,) array of fluxes, optional (for logging)
        aliases           : list from load_id_aliases(); None loads default
        db_path           : local Gaia database; None uses $GAIADATABASEPATH
        info              : dict, optional; filled with diagnostics:
                            status, primary_method, primary_id, primary_row,
                            schedule_id, db_row, inject
        coord_match_radius: arcsec, acceptance radius for the fallback

    Output:
        list of (GAIA_ID: str, role: str, Teff: int or None)
        role is 'primary' or 'secondary'
        Empty list if the catalogue contains no usable IDs and no coords.
    """
    if info is None:
        info = {}
    info.update({'status': 'NO_PRIMARY', 'primary_method': None,
                 'primary_id': None, 'primary_row': None,
                 'schedule_id': None, 'db_row': None, 'inject': None})
    if aliases is None:
        aliases = load_id_aliases()

    dr2_clean = [_clean_id(x) for x in catalogue_dr2_ids]
    if catalogue_dr3_ids is not None and len(catalogue_dr3_ids) == len(dr2_clean):
        dr3_clean = [_clean_id(x) for x in catalogue_dr3_ids]
    else:
        dr3_clean = [''] * len(dr2_clean)
    id_rows = {}
    for i, (a, b) in enumerate(zip(dr2_clean, dr3_clean)):
        for v in (a, b):
            if v and v not in id_rows:
                id_rows[v] = i
    clean_ids = list(id_rows)

    if not clean_ids and catalogue_coords is None:
        logger.error(
            "Catalogue contains no usable Gaia IDs — cannot identify targets"
        )
        return []

    match_gaia = []   # list of (ID, role, Teff)
    primary_found = False

    # ------------------------------------------------------------------
    # Step 1: schedule
    # ------------------------------------------------------------------
    dr2 = get_target_from_schedule(obsdir, date, targname, aliases=aliases)
    info['schedule_id'] = dr2
    id_source = 'schedule'
    list_coords = None
    list_teff = None

    # ------------------------------------------------------------------
    # Step 1a: no plan file for this night — resolve the target NAME via
    # the 40 pc target list instead (roughly half of some targets' nights
    # have no plan file on disk).
    # ------------------------------------------------------------------
    if dr2 is None:
        by_name = get_target_from_target_list_by_name(
            targname, target_list_path, catalogue_ids=set(clean_ids))
        if by_name is not None:
            dr2, lra, ldec, list_teff = by_name
            dr2, _ = apply_id_alias(dr2, targname, date, aliases)
            id_source = 'target_list_name'
            info['schedule_id'] = dr2
            if lra is not None and ldec is not None:
                list_coords = (lra, ldec)

    if dr2 is not None:
        if dr2 in id_rows:
            match_gaia.append((dr2, 'primary', None))
            primary_found = True
            info.update(status='OK', primary_method=id_source,
                        primary_id=dr2, primary_row=id_rows[dr2])
            logger.info(
                "Primary target identified from %s: %s", id_source, dr2
            )
        else:
            plan_coords = get_coords_from_schedule(obsdir, date, targname) \
                or list_coords
            dec_hint = None
            if plan_coords is not None:
                dec_hint = plan_coords[1]
            elif catalogue_coords is not None and len(catalogue_coords):
                dec_hint = float(np.nanmedian(np.asarray(catalogue_coords, dtype=float)[:, 1]))
            dbrow = lookup_gaia_id_in_db(dr2, dec_hint=dec_hint, db_path=db_path)
            info['db_row'] = dbrow
            if dbrow is None:
                info['status'] = 'PRIMARY_NOT_IN_DB'
                logger.error(
                    "PRIMARY_NOT_IN_DB: schedule ID %s for '%s' is not in "
                    "the local Gaia database, so it can never be crossmatched "
                    "on any night. Fix the database build (gaia-tmass-sqlite) "
                    "or add the object to calibration/supplementary_sources.csv.",
                    dr2, targname
                )
            else:
                info['status'] = 'PRIMARY_NOT_IN_FOV'
                logger.warning(
                    "PRIMARY_NOT_IN_FOV: schedule ID %s for '%s' is in the "
                    "local Gaia database (G=%s, table %s) but not in tonight's "
                    "catalogue — target may not have been observed, or the "
                    "crossmatch missed it",
                    dr2, targname, dbrow.get('phot_g_mean_mag'), dbrow.get('table')
                )
            # --------------------------------------------------------------
            # Step 1b: coordinate fallback
            # --------------------------------------------------------------
            if catalogue_coords is not None and len(catalogue_coords):
                row = _coordinate_fallback(
                    dr2, targname, date, plan_coords, dbrow,
                    catalogue_coords, catalogue_fluxes, dr2_clean, dr3_clean,
                    coord_match_radius, info)
                if row is not None:
                    match_gaia.append((dr2, 'primary', None))
                    primary_found = True
                    info.update(status='OK_COORDINATES',
                                primary_method='coordinates',
                                primary_id=dr2, primary_row=row)
                    info['id_source'] = id_source

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
                if dr2 in id_rows:
                    match_gaia.append((dr2, 'primary', None))
                    primary_found = True
                    info.update(status='OK', primary_method='toi',
                                primary_id=dr2, primary_row=id_rows[dr2])
                    logger.info(
                        "Primary target identified from TOI table: %s", dr2
                    )
                else:
                    logger.warning(
                        "TOI ID '%s' not found in FOV catalogue", dr2
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
            "PRIMARY TARGET COULD NOT BE IDENTIFIED BY ANY METHOD for '%s' "
            "(status %s)", targname, info['status']
        )

    # ------------------------------------------------------------------
    # Step 3: target list crossmatch (always runs; matches DR2 or DR3 IDs)
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
            elif list_teff is not None:
                match_gaia[i] = (dr2_id, role, list_teff)

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

def catalogue_is_valid(catalogue_path, target_dr2_id, require_pm=False):
    """
    Check whether a stack catalogue contains a valid Gaia crossmatch
    entry for a specific target.

    A catalogue is considered valid for a target if:
        - The Gaia_Crossmatch FITS extension exists
        - The target's ID appears in the GAIA_DR2_ID or GAIA_DR3_ID column
        - If require_pm: that row has finite PMRA and PMDEC, so the
          catalogue can be propagated to another night's epoch.  A
          catalogue reused across nights for a primary target without
          proper motion drifts off the star (5"/yr for Wolf 359).

    Inputs:
        catalogue_path : str  Path to the stack catalogue FITS file
        target_dr2_id  : str  Gaia ID (DR2 or DR3) of the target to check for
        require_pm     : bool

    Output:
        bool
    """
    target = _clean_id(target_dr2_id)
    try:
        with fits.open(catalogue_path) as hdul:
            ext_names = [hdu.name.upper() for hdu in hdul]
            if 'GAIA_CROSSMATCH' not in ext_names:
                logger.warning(
                    "No Gaia_Crossmatch extension in '%s'", catalogue_path
                )
                return False
            g = hdul['GAIA_CROSSMATCH'].data
            cols = [c.upper() for c in g.columns.names]
            rows = set()
            for col in ('GAIA_DR2_ID', 'GAIA_DR3_ID'):
                if col in cols:
                    for i, v in enumerate(g[col]):
                        if _clean_id(v) == target:
                            rows.add(i)
            if not rows:
                logger.warning(
                    "Target '%s' not found in catalogue '%s'",
                    target, catalogue_path
                )
                return False
            if require_pm:
                has_pm = False
                for i in rows:
                    try:
                        pmra = float(g['PMRA'][i])
                        pmdec = float(g['PMDEC'][i])
                    except (KeyError, ValueError, TypeError):
                        continue
                    if math.isfinite(pmra) and math.isfinite(pmdec):
                        has_pm = True
                if not has_pm:
                    logger.warning(
                        "Target '%s' is in catalogue '%s' but has no proper "
                        "motion — catalogue cannot be reused on another night",
                        target, catalogue_path
                    )
                    return False
    except Exception as e:
        logger.warning(
            "Could not open catalogue '%s': %s", catalogue_path, e
        )
        return False
    logger.info(
        "Target '%s' found in catalogue '%s'", target, catalogue_path
    )
    return True


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


def find_backup(backup_dir, primary_gaia_id, instrument, obs_date,
                max_age_without_pm_days=30):
    """
    Search for the best backup stack catalogue for a given primary target
    and instrument, selecting the closest date to obs_date whose catalogue
    actually contains the target (by DR2 or DR3 ID).

    The search key is (primary_gaia_id, instrument) — telescope and filter
    are ignored for matching.  A candidate more than max_age_without_pm_days
    from obs_date is accepted only if the target has a proper motion in
    it, because the catalogue positions are propagated to the night and a
    row without proper motion stays where it was.

    Inputs:
        backup_dir       : str  Path to shared StackImages backup directory
        primary_gaia_id  : str  Gaia ID of the primary (scheduled) target
        instrument       : str  Instrument name e.g. 'andor', 'spirit'
        obs_date         : str  Observation date YYYYMMDD
        max_age_without_pm_days : int

    Output:
        (stack_path: str, cat_path: str) if found, else (None, None)
    """
    backup_dir_path = Path(backup_dir)
    pattern = f'{primary_gaia_id}_*_{instrument}_*_stack_catalogue_*.fits'
    cat_candidates = list(backup_dir_path.glob(pattern))

    if not cat_candidates:
        logger.warning(
            "No backup catalogue found for primary target '%s', "
            "instrument '%s' in '%s'",
            primary_gaia_id, instrument, backup_dir
        )
        return None, None

    try:
        obs_dt = datetime.strptime(obs_date, '%Y%m%d')
    except ValueError:
        logger.warning("Could not parse obs_date '%s' as YYYYMMDD", obs_date)
        return None, None

    ranked = []
    for cat_path in cat_candidates:
        date_str = _parse_date_from_filename(cat_path.name)
        if date_str is None:
            continue
        try:
            cat_dt = datetime.strptime(date_str, '%Y%m%d')
        except ValueError:
            continue
        stack_path = backup_dir_path / cat_path.name.replace(
            'stack_catalogue', 'outstack')
        if stack_path.exists():
            ranked.append((abs((obs_dt - cat_dt).days), cat_path, stack_path))
    ranked.sort(key=lambda x: x[0])

    for dist, cat_path, stack_path in ranked:
        if catalogue_is_valid(str(cat_path), str(primary_gaia_id),
                              require_pm=(dist > max_age_without_pm_days)):
            logger.info(
                "Backup found for primary target '%s', instrument '%s': "
                "%s (distance %d days from %s)",
                primary_gaia_id, instrument, cat_path.name, dist, obs_date
            )
            return str(stack_path), str(cat_path)
        logger.warning(
            "Backup candidate %s rejected for '%s' (%d days from %s)",
            cat_path.name, primary_gaia_id, dist, obs_date
        )

    logger.warning(
        "No usable backup pair (stack + catalogue containing the target) "
        "found for primary target '%s', instrument '%s'",
        primary_gaia_id, instrument
    )
    return None, None


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