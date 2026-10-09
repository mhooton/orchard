#!/usr/bin/env python3
"""
promote_check.py
================

Read-only prototype of the promotion check that would replace T12's sweep.

Why
---
T12 (perform_version_migration in main/ZLP_pipeline.sh) moves any v3 target
directory holding a *4_diff.fits over the v2 one, the tree the portal reads,
without looking at either.  The 2026-10-09 audit of that sweep found 13
target-nights where it replaced a better light curve with a worse one, and
93 good light curves stranded in v3 because the run died before T12.

What
----
For every target-night it is given, this compares the candidate (normally
the v3 target directory) with the incumbent at the same path in v2 and
gives one verdict, with the reasons and the numbers behind them.  Nothing
is moved, copied or deleted.  The only files written are the output CSV and,
if asked for, a cache of raw-header counts.

Checks on the candidate alone
    A1  an aperture-5 light curve (*_5_diff.fits) exists for the target;
        the portal ingests nothing else.  A directory with none for any star
        is rejected; one with aperture 5 only for other stars is held
    A2  finite target points in LIGHTCURVE_5 (target = CATALOGUE TARGET==1)
    A3  the light curve's Gaia ID (DR3ID / DR2ID / GAIA_ID header, file
        prefix, target row) against line 4 of the schedule plan file
        Observations/<tel>/schedule/Plans_by_date/YYYY-MM-DD/Obj_<target>.txt,
        with calibration/target_id_aliases.csv applied and, where the IDs
        differ, the DR2/DR3 counterpart taken from the local Gaia database.
        Nights without a plan file fall back to the 40 pc target list by
        name, as identify_targets does
    A4  frames in the light curve against the raw science frames of that
        target and filter on disk, classified from the raw headers the way
        createlists.py does it (IMAGETYP containing LIGHT, OBJECT = target)
    A5  point-to-point scatter, and its ratio to the pipeline's own per-point
        error: a curve far below its own error bar is degenerate (the March
        2026 Callisto reruns have p2p ~4e-5 against errors of ~2e-2)

Checks against the incumbent, when v2 holds a usable aperture-5 light curve
    R1  no fewer finite target points
    R2  p2p scatter no more than 20% worse
    R3  identity not moved away from the plan's ID
R1 and R2 also cover every other star v2 serves from the directory (secondary
40 pc targets get their own light curves there), because T12 replaces the
whole directory: in the audit, Europa 2025-08-09 Sp0102-6322 lost 37% in p2p
on the secondary while the primary changed by 6%.

Verdicts
    REJECT   unusable on its own (no aperture 5 at all, no finite target
             point, degenerate), or worse than v2 beyond argument: v2 has
             the right star and the candidate a different one, or v2 has an
             aperture-5 light curve of the target and the candidate none
    HOLD     usable, but loses to v2 on R1/R2, or its identity or
             completeness cannot be confirmed: for a person to decide
    PROMOTE  passed everything

Scatter is measured exactly as the sweep audit's lc_compare.py did, so the
numbers can be checked against lc_compare.csv: every finite target point,
normalised by its median, p2p = 1.4826 * median|successive difference| / sqrt(2).

Usage (inside the orchard-server container)
-----
    python utils/promote_check.py scan --out verdicts.csv [--tel Europa] [--date 20251103]
    python utils/promote_check.py pairs pairs.csv --out verdicts.csv
    python utils/promote_check.py redecide verdicts.csv --out new.csv [--max-scatter-ratio 1.3 ...]

pairs.csv has columns tel,date,target,candidate,incumbent[,label], where
candidate and incumbent are target directories (incumbent may be blank).
redecide applies new thresholds to an earlier output without reopening any
FITS file.
"""

import argparse
import contextlib
import csv
import glob
import io
import json
import logging
import math
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
from astropy.io import fits

from utils import gaia_id_from_schedule
from utils import target_management

BASE = '/data/SPECULOOSPipeline'
APERTURE = 5                        # the portal ingests only *_5_diff.fits
RAW_EXTS = ('fits', 'fts', 'fz', 'fit')

THRESHOLDS = {
    'min_valid_points': 20,         # A2: fewer finite target points holds (none at all rejects)
    'min_valid_fraction': 0.5,      # A2: finite / frames in the light curve below this holds
    'min_frame_fraction': 0.8,      # A4: frames in the light curve / raw science frames below this holds
    'min_noise_ratio': 0.05,        # A5: p2p / median pipeline error below this rejects (degenerate)
    'max_p2p': 0.05,                # A5: p2p above this (relative flux) holds
    'min_points_ratio': 1.0,        # R1: candidate / incumbent finite points below this holds
    'max_scatter_ratio': 1.2,       # R2: candidate / incumbent p2p above this holds
}

COLUMNS = [
    'tel', 'date', 'target', 'verdict', 'reasons', 'notes', 'label',
    'expected_id', 'expected_src',
    'c_dir', 'c_stars', 'c_ap5_stars', 'c_primary', 'c_pick', 'c_id_match', 'c_filter', 'c_aps', 'c_ap5',
    'c_file', 'c_pipe_v', 'c_frames', 'c_valid', 'c_p2p', 'c_rstd', 'c_err', 'c_noise_ratio', 'c_bw',
    'raw_frames', 'c_frame_frac',
    'i_dir', 'i_stars', 'i_primary', 'i_pick', 'i_id_match', 'i_filter', 'i_aps', 'i_ap5', 'i_file',
    'i_pipe_v', 'i_frames', 'i_valid', 'i_p2p', 'i_rstd', 'i_err', 'i_noise_ratio',
    'same_star', 'points_ratio', 'scatter_ratio', 'others',
]
NUMERIC = {'c_frames', 'c_valid', 'c_p2p', 'c_rstd', 'c_err', 'c_noise_ratio', 'c_bw', 'raw_frames',
           'c_frame_frac', 'i_frames', 'i_valid', 'i_p2p', 'i_rstd', 'i_err', 'i_noise_ratio',
           'points_ratio', 'scatter_ratio'}
VERDICT_RANK = {'PROMOTE': 0, 'HOLD': 1, 'REJECT': 2}


# ---------------------------------------------------------------------------
# One light curve
# ---------------------------------------------------------------------------

def _clean_id(value):
    v = str(value).strip() if value is not None else ''
    return '' if v.upper() in ('', 'N', 'NAN', 'NONE', '0') else v


def p2p_scatter(y):
    """1.4826 * median |successive difference| / sqrt(2) of a normalised curve (NaN below 10 points)."""
    if y.size < 10:
        return np.nan
    return 1.4826 * np.median(np.abs(np.diff(y))) / np.sqrt(2)


def parse_diff_name(path):
    """(gaia_id, filter, aperture) from '<id>_<filter>_<date>_<ap>_diff.fits', or None."""
    parts = os.path.basename(path)[:-len('_diff.fits')].rsplit('_', 3)
    if len(parts) != 4 or not parts[3].isdigit():
        return None
    return parts[0], parts[1], int(parts[3])


def read_light_curve(path):
    """
    Target light-curve metrics from one *_diff.fits; never raises.

    The target column is the CATALOGUE row with TARGET==1 (SPlightcurve puts
    the target first), falling back to the row whose Gaia ID is the file
    prefix.  'ids' collects every Gaia ID the file attaches to the target.
    """
    gid, _, ap = parse_diff_name(path)
    m = {'file': os.path.basename(path), 'ids': {gid}, 'pipe_v': '', 'frames': 0, 'valid': 0,
         'p2p': np.nan, 'rstd': np.nan, 'err': np.nan, 'noise_ratio': np.nan, 'bw': 0, 'error': ''}
    try:
        with fits.open(path, memmap=True) as h:
            hdus = [x.name for x in h]
            hdr = h[0].header
            m['pipe_v'] = str(hdr.get('PIPE_V', '')).strip()
            m['ids'] |= {_clean_id(hdr.get(k)) for k in ('DR3ID', 'DR2ID', 'GAIA_ID')}
            cat = h['CATALOGUE'].data
            cols = cat.columns.names
            ti = None
            if 'TARGET' in cols:
                w = np.where(np.asarray(cat['TARGET']) == 1)[0]
                if len(w):
                    ti = int(w[0])
            if ti is None:
                for c in ('GAIA_DR3_ID', 'GAIA_DR2_ID', 'OBJ_ID'):
                    if c in cols:
                        w = np.where(np.char.strip(np.asarray(cat[c]).astype(str)) == gid)[0]
                        if len(w):
                            ti = int(w[0])
                            break
            if ti is None:
                m['error'] = 'target row not found'
                return m
            for c in ('GAIA_DR3_ID', 'GAIA_DR2_ID'):
                if c in cols:
                    m['ids'].add(_clean_id(cat[c][ti]))

            def target_column(name):
                data = h[name].data
                if data.shape[1] != len(cat) and data.shape[0] == len(cat):
                    return np.array(data[ti, :], dtype=float)
                return np.array(data[:, ti], dtype=float)

            y = target_column('LIGHTCURVE_%d' % ap)
            finite = np.isfinite(y)
            m['frames'] = int(y.size)
            m['valid'] = int(finite.sum())
            if m['valid']:
                med = np.median(y[finite])
                if med > 0:
                    yn = y[finite] / med
                    m['p2p'] = float(p2p_scatter(yn))
                    m['rstd'] = float(1.4826 * np.median(np.abs(yn - np.median(yn))))
                    if 'ERROR_%d' % ap in hdus:
                        e = target_column('ERROR_%d' % ap)[finite]
                        e = e[np.isfinite(e)]
                        if e.size:
                            m['err'] = float(np.median(e) / med)
                    if np.isfinite(m['p2p']) and m['err'] > 0:
                        m['noise_ratio'] = m['p2p'] / m['err']
            if 'FLAGS' in hdus and 'BW_FLAG' in h['FLAGS'].columns.names:
                m['bw'] = int(np.sum(np.asarray(h['FLAGS'].data['BW_FLAG']) > 0))
    except Exception as ex:          # report, never stop
        m['error'] = '%s: %s' % (type(ex).__name__, str(ex)[:80])
    m['ids'].discard('')
    return m


# ---------------------------------------------------------------------------
# One target directory
# ---------------------------------------------------------------------------

def survey(target_dir):
    """
    Every star with a light curve in a target directory.

    Returns {gaia_id: {'ids': set, 'filters': {filter: {'aps': [..], 'lc': metrics or None}}}};
    'lc' is the aperture-5 light curve.  Secondary 40 pc targets get their own
    files beside the primary's, so a directory can hold several stars.
    """
    stars = {}
    if not target_dir or not os.path.isdir(target_dir):
        return stars
    found = defaultdict(lambda: defaultdict(dict))
    for path in glob.glob(os.path.join(target_dir, '*_diff.fits')):
        parsed = parse_diff_name(path)
        if parsed:
            gid, flt, ap = parsed
            found[gid][flt][ap] = path
    for gid, by_filter in found.items():
        star = {'ids': {gid}, 'filters': {}}
        for flt, aps in by_filter.items():
            lc = read_light_curve(aps[APERTURE]) if APERTURE in aps else None
            if lc:
                star['ids'] |= lc['ids']
            star['filters'][flt] = {'aps': sorted(aps), 'lc': lc}
        if not has_ap5(star):
            # no light curve to open, so take the IDs from any aperture's header
            any_path = sorted(p for aps in by_filter.values() for p in aps.values())[0]
            try:
                hdr = fits.getheader(any_path, 0)
                star['ids'] |= {_clean_id(hdr.get(k)) for k in ('DR3ID', 'DR2ID', 'GAIA_ID')}
            except Exception:
                pass
            star['ids'].discard('')
        stars[gid] = star
    return stars


def has_ap5(star):
    return any(f['lc'] is not None for f in star['filters'].values())


def pick_primary(stars, expected):
    """
    The scheduled target among a directory's stars: (gaia_id, how).

    how = 'plan' when its IDs include the plan's, 'single' when it is the
    only star with light curves, 'ambiguous' when several stars have them
    and none is the plan's (the one with most finite points is taken),
    '' for an empty directory.
    """
    if not stars:
        return '', ''
    hits = [g for g, s in stars.items() if s['ids'] & expected]
    if hits:
        return sorted(hits, key=lambda g: (not has_ap5(stars[g]), g))[0], 'plan'
    if len(stars) == 1:
        return next(iter(stars)), 'single'

    def best(g):
        lcs = [f['lc'] for f in stars[g]['filters'].values() if f['lc']]
        return max([lc['valid'] for lc in lcs] or [-1])
    return sorted(stars, key=lambda g: (-best(g), g))[0], 'ambiguous'


# ---------------------------------------------------------------------------
# Expected identity
# ---------------------------------------------------------------------------

class Identity(object):
    """The Gaia ID a target-night should carry, and its other-release counterpart (cached)."""

    def __init__(self, base, target_list=None):
        self.base = base
        # the pipeline's own default: TARGET_LIST, else <basedir>/ml_40pc.txt
        self.target_list = target_list or os.environ.get('TARGET_LIST') or os.path.join(base, 'ml_40pc.txt')
        self.aliases = target_management.load_id_aliases()
        self.plans = {}
        self.by_name = {}
        self.crosswalk = {}

    def expected(self, tel, date, target, catalogue_ids=()):
        """
        (ID or '', where it came from, Dec or None).  Line 4 of the plan file
        with the alias table applied; with no plan file or no ID in it (about
        half the nights of some targets), the 40 pc target list by name, as
        identify_targets does it.  catalogue_ids picks between list rows that
        share a name (resolved binaries).
        """
        key = (tel, date, target)
        if key not in self.plans:
            obsdir = os.path.join(self.base, 'Observations', tel)
            files = gaia_id_from_schedule.find_plan(obsdir, date, target.replace('--', ' '))
            gid, src, dec = '', 'no plan file', None
            if files:
                with contextlib.redirect_stdout(io.StringIO()):
                    raw_id = gaia_id_from_schedule.read_file(files[0], silent=True)
                if raw_id:
                    gid, alias = target_management.apply_id_alias(raw_id.strip(), target, date, self.aliases)
                    src = 'plan, alias %s' % raw_id.strip() if alias else 'plan'
                    coords = gaia_id_from_schedule.read_coords(files[0])
                    dec = coords[1] if coords else None
                else:
                    src = 'no ID in plan file'
            self.plans[key] = (gid, src, dec)
        gid, src, dec = self.plans[key]
        if gid:
            return gid, src, dec
        name_key = (target, tuple(sorted(catalogue_ids)))
        if name_key not in self.by_name:
            hit = target_management.get_target_from_target_list_by_name(
                target, self.target_list, set(catalogue_ids) or None)
            self.by_name[name_key] = (hit[0], hit[2]) if hit else ('', None)
        list_id, list_dec = self.by_name[name_key]
        if list_id:
            return list_id, src + ', target list', list_dec
        return '', src + ', not in target list', None

    def counterpart(self, gid, dec):
        """{DR3, DR2} IDs of gid from the local Gaia DB (empty if absent); ~3 s, so only on demand."""
        if gid not in self.crosswalk:
            row = target_management.lookup_gaia_id_in_db(gid, dec_hint=dec)
            ids = {_clean_id(row.get('source_id')), _clean_id(row.get('dr2_source_id'))} if row else set()
            ids.discard('')
            self.crosswalk[gid] = ids
        return self.crosswalk[gid]


# ---------------------------------------------------------------------------
# Raw science frames
# ---------------------------------------------------------------------------

def _read_cards(f, keys):
    """Values of the wanted keys from one FITS header; (dict, reached END)."""
    out = {}
    while True:
        block = f.read(2880)
        if len(block) < 2880:
            return out, False
        for i in range(0, 2880, 80):
            card = block[i:i + 80].decode('ascii', 'replace')
            key = card[:8].strip()
            if key == 'END':
                return out, True
            if key in keys and card[8:10] == '= ':
                v = card[10:].strip()
                if v.startswith("'"):
                    v = v[1:]
                    v = v[:v.find("'")] if "'" in v else v
                else:
                    v = v.split('/')[0]
                out[key] = v.strip()


RAW_KEYS = {'IMAGETYP', 'OBJECT', 'FIELD', 'FILTER', 'DITHER'}


def raw_header(path):
    """IMAGETYP, OBJECT, FILTER, DITHER of a raw frame (the image HDU's for .fz)."""
    with open(path, 'rb') as f:
        h, ok = _read_cards(f, RAW_KEYS)
        if path.endswith('.fz') and ok:
            h.update(_read_cards(f, RAW_KEYS)[0])
    return h


def _raw_header_or_none(path):
    try:
        return raw_header(path)
    except OSError:
        return None


def _norm_target(name):
    return name.replace('--', ' ').strip().upper()


def _norm_filter(name):
    return name.replace("'", '').replace(' ', '').lower()


class RawFrames(object):
    """
    Raw science frames per target and filter on a night, classified from the
    headers as createlists.py does it: IMAGETYP containing LIGHT and no
    dithering, OBJECT (or FIELD) as the target, '_2' stripped from 13-character
    names.  Frames are counted per file extension and the largest count is
    used, since some nights hold the same frames as both .fits and .fz.
    Counts are cached per night in a JSON file when one is given.
    """

    def __init__(self, base, cache_path=None, threads=8):
        self.base = base
        self.cache_path = cache_path
        self.threads = threads
        self.cache = {}
        if cache_path and os.path.exists(cache_path):
            with open(cache_path) as f:
                self.cache = json.load(f)

    def night(self, tel, date):
        key = '%s/%s' % (tel, date)
        if key not in self.cache:
            imgdir = os.path.join(self.base, 'Observations', tel, 'images', date)
            counts = None
            if os.path.isdir(imgdir):
                counts = defaultdict(Counter)
                paths = glob.glob(os.path.join(imgdir, '*.*')) + glob.glob(os.path.join(imgdir, '*', '*.*'))
                paths = [p for p in paths if p.rsplit('.', 1)[-1] in RAW_EXTS]
                pool = ThreadPool(self.threads)       # header reads are disk-latency bound
                headers = pool.map(_raw_header_or_none, paths)
                pool.close()
                for p, h in zip(paths, headers):
                    ext = p.rsplit('.', 1)[-1]
                    if h is None or 'LIGHT' not in h.get('IMAGETYP', '').upper() or h.get('DITHER', '') == 'ENABLED':
                        continue
                    field = h.get('OBJECT', h.get('FIELD', ''))
                    if '_2' in field and len(field) == 13:
                        field = field.replace('_2', '')
                    counts[ext]['%s\t%s' % (_norm_target(field), _norm_filter(h.get('FILTER', '')))] += 1
                counts = {e: dict(c) for e, c in counts.items()}
            self.cache[key] = counts
            if self.cache_path:
                with open(self.cache_path + '.tmp', 'w') as f:
                    json.dump(self.cache, f)
                os.replace(self.cache_path + '.tmp', self.cache_path)
        return self.cache[key]

    def count(self, tel, date, target, flt):
        """Raw science frames of the target in this filter, or None if the night has no raw directory."""
        counts = self.night(tel, date)
        if counts is None:
            return None
        want = _norm_target(target)
        best = 0
        for by_key in counts.values():
            mine = {k.split('\t')[1]: n for k, n in by_key.items() if k.split('\t')[0] == want}
            # filter names differ between headers and file names (r' vs r), so
            # only insist on the filter when the target was observed in several
            n = mine.get(_norm_filter(flt), 0) if len(mine) > 1 else sum(mine.values())
            best = max(best, n)
        return best


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def _finite(x):
    try:
        return math.isfinite(x)
    except TypeError:
        return False


def incumbent_usable(r, th):
    """Why the incumbent's light curve cannot serve as a reference ('' if it can)."""
    if not r['i_dir']:
        return 'no v2 directory'
    if not r['i_ap5']:
        return 'v2 has no aperture-5 light curve'
    if r['i_valid'] < th['min_valid_points']:
        return 'v2 light curve has %d finite points' % r['i_valid']
    if _finite(r['i_noise_ratio']) and r['i_noise_ratio'] < th['min_noise_ratio']:
        return 'v2 light curve is degenerate (p2p %.2g = %.3f x its error)' % (r['i_p2p'], r['i_noise_ratio'])
    return ''


def decide(r, th=THRESHOLDS):
    """Verdict, reasons and notes for one measured row; reads nothing from disk."""
    reject, hold, notes = [], [], []
    if r['c_error']:
        notes.append('candidate read error: %s' % r['c_error'])

    # A1-A5: the candidate alone
    if not r['c_ap5']:
        aps = ' (apertures %s)' % r['c_aps'] if r['c_aps'] else ''
        others_ap5 = [g for g in r['c_ap5_stars'].split() if g != r['c_primary']]
        if others_ap5:
            hold.append('A1 no aperture-5 light curve of the target%s, only of %s' % (aps, ' '.join(others_ap5)))
        else:
            reject.append('A1 no aperture-5 light curve' + aps)
    else:
        if r['c_valid'] == 0:
            reject.append('A2 no finite target point in %d frames' % r['c_frames'])
        elif r['c_valid'] < th['min_valid_points']:
            hold.append('A2 only %d finite target points' % r['c_valid'])
        elif r['c_valid'] < th['min_valid_fraction'] * r['c_frames']:
            hold.append('A2 only %d of %d target points finite' % (r['c_valid'], r['c_frames']))
        if _finite(r['c_noise_ratio']) and r['c_noise_ratio'] < th['min_noise_ratio']:
            reject.append('A5 degenerate: p2p %.2g is %.3f x the pipeline error' % (r['c_p2p'], r['c_noise_ratio']))
        elif _finite(r['c_p2p']) and r['c_p2p'] > th['max_p2p']:
            hold.append('A5 p2p %.3f above %.3f' % (r['c_p2p'], th['max_p2p']))
        if _finite(r['raw_frames']) and r['raw_frames'] > 0:
            if r['c_frames'] < th['min_frame_fraction'] * r['raw_frames']:
                hold.append('A4 light curve has %d of %d raw science frames (%.0f%%)'
                            % (r['c_frames'], r['raw_frames'], 100.0 * r['c_frames'] / r['raw_frames']))
        else:
            notes.append('A4 no raw science frames found for the target')
    if r['c_id_match'] == 'no':
        hold.append('A3 light curve is %s, expected %s (%s)' % (r['c_primary'], r['expected_id'], r['expected_src']))
    elif r['c_id_match'] == 'unknown':
        notes.append('A3 %s, identity unchecked' % r['expected_src'])
    if r['c_pick'] == 'ambiguous':
        hold.append('A3 several stars have light curves and %s (%s)'
                    % ('none is the expected ID' if r['expected_id'] else 'no expected ID to choose', r['c_stars']))

    # R1-R3: against the incumbent
    why_not = incumbent_usable(r, th)
    if why_not:
        notes.append('no comparison: %s' % why_not)
    else:
        comparable = r['same_star'] == 'yes'
        if r['expected_id']:
            if r['i_id_match'] == 'yes' and r['c_id_match'] != 'yes':
                reject.append('R3 identity moved away from the plan ID: v2 is %s, candidate %s'
                              % (r['i_primary'], r['c_primary']))
            elif r['c_id_match'] == 'yes' and r['i_id_match'] != 'yes':
                notes.append('R3 identity corrected: v2 was %s' % r['i_primary'])
                comparable = False
        elif not comparable:
            hold.append('R3 star changed (v2 %s, candidate %s) with no plan ID to arbitrate'
                        % (r['i_primary'], r['c_primary']))
        if comparable and not r['c_ap5']:
            reject.append('R1 v2 has an aperture-5 light curve of the target (%d points), the candidate none'
                          % r['i_valid'])
        elif comparable:
            if r['c_valid'] < th['min_points_ratio'] * r['i_valid']:
                hold.append('R1 %d finite points against %d in v2' % (r['c_valid'], r['i_valid']))
            if _finite(r['c_p2p']) and _finite(r['i_p2p']) and r['c_p2p'] > th['max_scatter_ratio'] * r['i_p2p']:
                hold.append('R2 p2p %.4f against %.4f in v2 (x%.2f)' % (r['c_p2p'], r['i_p2p'], r['c_p2p'] / r['i_p2p']))
    # R1-R2 for the other stars v2 serves from this directory, which T12 would replace too
    for entry in r['others'].split():
        gid, flt, c_valid, i_valid, c_p2p, i_p2p, i_ratio = entry.split('|')
        c_valid, i_valid = int(c_valid), int(i_valid)
        c_p2p, i_p2p, i_ratio = (float(x) if x else np.nan for x in (c_p2p, i_p2p, i_ratio))
        if i_valid < th['min_valid_points'] or (_finite(i_ratio) and i_ratio < th['min_noise_ratio']):
            continue
        if c_valid < th['min_points_ratio'] * i_valid:
            hold.append('R1 %s %s: %d finite points against %d in v2' % (gid, flt, c_valid, i_valid))
        elif _finite(c_p2p) and _finite(i_p2p) and c_p2p > th['max_scatter_ratio'] * i_p2p:
            hold.append('R2 %s %s: p2p %.4f against %.4f in v2 (x%.2f)' % (gid, flt, c_p2p, i_p2p, c_p2p / i_p2p))

    verdict = 'REJECT' if reject else 'HOLD' if hold else 'PROMOTE'
    return verdict, reject + hold, notes


# ---------------------------------------------------------------------------
# One target-night
# ---------------------------------------------------------------------------

def _side(prefix, target_dir, stars, gid, how, flt, expected):
    """Columns describing one side (candidate 'c_' or incumbent 'i_') for one filter."""
    star = stars.get(gid)
    lc = None
    flt_used = ''
    id_match = ''
    if star:
        filters = star['filters']
        flt_used = flt if flt in filters else sorted(filters, key=lambda f: filters[f]['lc'] is None)[0]
        lc = filters[flt_used]['lc']
        id_match = 'unknown' if not expected else 'yes' if star['ids'] & expected else 'no'
    aps = sorted({a for f in star['filters'].values() for a in f['aps']}) if star else []
    r = {
        prefix + 'dir': target_dir if target_dir and os.path.isdir(target_dir) else '',
        prefix + 'stars': ' '.join(sorted(stars)),
        prefix + 'primary': gid,
        prefix + 'pick': how,
        prefix + 'id_match': id_match,
        prefix + 'filter': flt_used,
        prefix + 'aps': ' '.join(str(a) for a in aps),
        prefix + 'ap5': bool(lc),
        prefix + 'file': lc['file'] if lc else '',
        prefix + 'pipe_v': lc['pipe_v'] if lc else '',
    }
    for k in ('frames', 'valid', 'p2p', 'rstd', 'err', 'noise_ratio'):
        r[prefix + k] = lc[k] if lc else (0 if k in ('frames', 'valid') else np.nan)
    if prefix == 'c_':
        r['c_bw'] = lc['bw'] if lc else 0
        r['c_error'] = lc['error'] if lc else ''
    return r


def assess(tel, date, target, cand_dir, inc_dir, identity, raw, th=THRESHOLDS, label=''):
    """Measure and decide one target-night; one row per filter of the candidate's primary star."""
    cand = survey(cand_dir)
    inc = survey(inc_dir)
    seen = set().union(*[s['ids'] for s in list(cand.values()) + list(inc.values())])
    plan_id, plan_src, plan_dec = identity.expected(tel, date, target, seen)
    expected = {plan_id} if plan_id else set()
    c_gid, c_how = pick_primary(cand, expected)
    i_gid, i_how = pick_primary(inc, expected)
    if expected and ((cand and c_how != 'plan') or (inc and i_how != 'plan')):
        # the plan may name the other Gaia release's ID for the same star
        expected |= identity.counterpart(plan_id, plan_dec)
        c_gid, c_how = pick_primary(cand, expected)
        i_gid, i_how = pick_primary(inc, expected)

    # every other star v2 serves an aperture-5 light curve for, with the candidate's
    # light curve of the same star: gid|filter|c_valid|i_valid|c_p2p|i_p2p|i_noise_ratio
    others = []
    for gid, star in sorted(inc.items()):
        if gid == i_gid:
            continue
        twin = next((s for s in cand.values() if s['ids'] & star['ids']), None)
        for flt, f in sorted(star['filters'].items()):
            if f['lc']:
                c_lc = twin['filters'].get(flt, {}).get('lc') if twin else None
                others.append('|'.join([gid, flt, str(c_lc['valid'] if c_lc else 0), str(f['lc']['valid']),
                                        _fmt(c_lc['p2p']) if c_lc else '', _fmt(f['lc']['p2p']),
                                        _fmt(f['lc']['noise_ratio'])]))

    filters = sorted(cand[c_gid]['filters']) if c_gid else ['']
    rows = []
    for flt in filters:
        r = {'tel': tel, 'date': date, 'target': target, 'label': label,
             'expected_id': ' '.join(sorted(expected)), 'expected_src': plan_src}
        r.update(_side('c_', cand_dir, cand, c_gid, c_how, flt, expected))
        r.update(_side('i_', inc_dir, inc, i_gid, i_how, flt, expected))
        r['raw_frames'] = raw.count(tel, date, target, r['c_filter'] or flt) if raw else np.nan
        if r['raw_frames'] is None:
            r['raw_frames'] = np.nan
        r['c_frame_frac'] = r['c_frames'] / r['raw_frames'] if _finite(r['raw_frames']) and r['raw_frames'] > 0 \
            else np.nan
        same = c_gid and i_gid and (cand[c_gid]['ids'] & inc[i_gid]['ids'])
        r['same_star'] = 'yes' if same else ('no' if c_gid and i_gid else '')
        r['points_ratio'] = r['c_valid'] / r['i_valid'] if r['i_valid'] else np.nan
        r['scatter_ratio'] = r['c_p2p'] / r['i_p2p'] if _finite(r['i_p2p']) and r['i_p2p'] > 0 else np.nan
        r['others'] = ' '.join(others) if not rows else ''      # judged once per target-night
        r['c_ap5_stars'] = ' '.join(sorted(g for g, s in cand.items() if has_ap5(s)))
        r['verdict'], reasons, notes = decide(r, th)
        r['reasons'] = '; '.join(reasons)
        r['notes'] = '; '.join(notes)
        rows.append(r)
    return combine(rows)


def combine(rows):
    """One row per target-night: the worst filter's, with the other filters' reasons appended."""
    worst = max(rows, key=lambda r: VERDICT_RANK[r['verdict']])
    if len(rows) > 1:
        worst = dict(worst)
        worst['reasons'] = '; '.join('[%s] %s' % (r['c_filter'], r['reasons']) for r in rows if r['reasons'])
        worst['notes'] = '; '.join('[%s] %s' % (r['c_filter'], r['notes']) for r in rows if r['notes'])
        worst['c_filter'] = ' '.join(r['c_filter'] for r in rows)
        worst['others'] = rows[0]['others']
    return worst


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def v3_target_nights(base, tels=None, dates=None):
    """(tel, date, target, v3 dir, v2 dir) for every v3 target directory holding a *_diff.fits."""
    po = os.path.join(base, 'PipelineOutput')
    for tel in sorted(tels or os.listdir(os.path.join(po, 'v3'))):
        outdir = os.path.join(po, 'v3', tel, 'output')
        if not os.path.isdir(outdir):
            continue
        for date in sorted(os.listdir(outdir)):
            if not re.match(r'^\d{8}$', date) or (dates and date not in dates):
                continue
            for target in sorted(os.listdir(os.path.join(outdir, date))):
                d = os.path.join(outdir, date, target)
                if target != 'reduction' and os.path.isdir(d) and glob.glob(os.path.join(d, '*_diff.fits')):
                    yield tel, date, target, d, os.path.join(po, 'v2', tel, 'output', date, target), ''


def read_pairs(path):
    with open(path, newline='') as f:
        for p in csv.DictReader(f):
            yield p['tel'], p['date'], p['target'], p['candidate'], p.get('incumbent', ''), p.get('label', '')


def _fmt(v):
    if isinstance(v, bool):
        return 'yes' if v else 'no'
    if isinstance(v, float):
        return '' if not math.isfinite(v) else '%.6g' % v
    return v


def write_rows(rows, path):
    """Write rows as they arrive, flushing each, so a long run leaves a usable partial file."""
    done = []
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS + ['c_error'], extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(v) for k, v in r.items()})
            f.flush()
            done.append(r)
    return done


def read_rows(path):
    """Rows of an earlier output, typed for decide()."""
    rows = []
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            for k in NUMERIC:
                r[k] = float(r[k]) if r[k] != '' else np.nan
            for k in ('c_frames', 'c_valid', 'i_frames', 'i_valid', 'c_bw'):
                r[k] = int(r[k]) if _finite(r[k]) else 0
            r['c_ap5'] = r['c_ap5'] == 'yes'
            r['i_ap5'] = r['i_ap5'] == 'yes'
            rows.append(r)
    return rows


def summarise(rows, out=sys.stdout):
    by = Counter(r['verdict'] for r in rows)
    print('%d target-nights: %s' % (len(rows), ', '.join('%s %d' % (v, by[v]) for v in VERDICT_RANK)), file=out)
    labels = sorted({r.get('label', '') for r in rows} - {''})
    for lab in labels:
        sub = Counter(r['verdict'] for r in rows if r.get('label') == lab)
        print('  %-40s %s' % (lab, '  '.join('%s %3d' % (v, sub[v]) for v in VERDICT_RANK)), file=out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['scan', 'pairs', 'redecide'])
    p.add_argument('input', nargs='?', help='pairs CSV (pairs) or an earlier output CSV (redecide)')
    p.add_argument('--out', required=True, help='output CSV, one row per target-night')
    p.add_argument('--base', default=BASE)
    p.add_argument('--tel', action='append', help='scan: only this telescope (repeatable)')
    p.add_argument('--date', action='append', help='scan: only this YYYYMMDD (repeatable)')
    p.add_argument('--raw-cache', help='JSON cache of raw-header counts per night')
    p.add_argument('--no-raw', action='store_true', help='skip the raw-frame check (A4)')
    for k, v in THRESHOLDS.items():
        p.add_argument('--' + k.replace('_', '-'), type=type(v), default=v)
    a = p.parse_args(argv)
    th = {k: getattr(a, k) for k in THRESHOLDS}
    logging.getLogger(target_management.__name__).setLevel(logging.ERROR)

    if a.mode == 'redecide':
        earlier = read_rows(a.input)      # read in full first: --out may be the same file

        def judged():
            for r in earlier:
                r['verdict'], reasons, notes = decide(r, th)
                r['reasons'], r['notes'] = '; '.join(reasons), '; '.join(notes)
                yield r
    else:
        identity = Identity(a.base)
        raw = None if a.no_raw else RawFrames(a.base, a.raw_cache)
        work = v3_target_nights(a.base, a.tel, a.date) if a.mode == 'scan' else read_pairs(a.input)

        def judged():
            for n, (tel, date, target, cand, inc, label) in enumerate(work, 1):
                r = assess(tel, date, target, cand, inc, identity, raw, th, label)
                print('%4d %-8s %s %-16s %s' % (n, tel, date, target, r['verdict']), file=sys.stderr)
                yield r
    summarise(write_rows(judged(), a.out))


if __name__ == '__main__':
    main()
