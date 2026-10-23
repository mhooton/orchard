"""
Choosing a night's calibration frames from the calibration register (calregister.py) for T1.

Each type of calibration comes from one night, the nearest to the science night on either side (the past wins a
tie), among the nights whose frames are compatible with the science frames:

  - same telescope and same configuration epoch: the same camera run (a run ends wherever the detector -
    unbinned pixel size, else shape - changes from one night to the next, e.g. ANDOR -> SPIRIT -> ANDOR on
    Callisto) and the same side of every hand-maintained break in calib_epoch_breaks.csv (camera remounts);
  - the same frame configuration: shape, binning, window and orientation (MaxIm's FLIPSTAT, Astra's
    ASTRAROT/ASTRAMIR) - a master of another shape would be left out by T2-T6 with only a warning;
  - bias and darks: the same readout mode (calibration/readout.py), and darks a detector temperature within
    DARK_TEMP_TOLERANCE of the science frames';
  - enough usable frames (calregister.flag_night): MIN_BIAS (per instrument), MIN_DARK, MIN_FLAT per filter.

Darks for an instrument with dark_matching_exptime (SPIRIT) are chosen per science exposure time (rounded to
seconds, as pipered.py matches them), MIN_DARK frames at that exposure; one night that covers every exposure is
used for all of them if it is within COVERING_SLACK_DAYS of each exposure's own best night.

Nothing is combined across nights.  When nothing matches, the night fails (Selection.errors: no bias, no dark at
all, no flat for any science filter); a partial gap is flagged (Selection.flags: a filter without flats, whose
frames T6 will not reduce; a SPIRIT exposure without darks, which T6 corrects with the highest-exposure master
scaled; darks outside the temperature tolerance on an ANDOR night; a night more than FAR_DAYS away).  Every
choice is logged with |delta days|
and written to <reduction>/<run>_calib_selection.json, which the PDF report shows.
"""
import csv
import json
import os

import numpy as np

from calibration import calregister as cr

# bias frames a night needs: Andor nights take 6-14 (the overscan removes the level, the master is a structure
# template, so a neighbour's 6-14 beat a night's own 2-4); SPIRIT nights of 2025-26 take 1-3, and its bias, the
# zero level of a camera without overscan, cancels against exposure-matched darks anyway (the old threshold, 2)
MIN_BIAS = {'spirit': 2}
MIN_BIAS_DEFAULT = 5
MIN_DARK = 3
MIN_FLAT = 4
COVERING_SLACK_DAYS = 3
# a choice further away than this is used but flagged (Callisto 20250407's J flats are 208 days later)
FAR_DAYS = 30
DARK_TEMP_TOLERANCE = {'spirit': 2.0}
DARK_TEMP_TOLERANCE_DEFAULT = 5.0
# instruments whose darks may come from outside the temperature tolerance (flagged) when nothing is inside it
DARK_TEMP_RELAX = ('andor',)
BREAKS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'calib_epoch_breaks.csv')
SELECTION_NAME = '{run}_calib_selection.json'


class SelectionError(RuntimeError):
    pass


# ------------------------------------------------------------------------------------- epochs

def load_breaks(path=BREAKS_FILE):
    """{telescope: [(first_night, applies_to)]} from the hand-maintained table (comment lines start with #)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as f:
        rows = csv.DictReader(line for line in f if line.strip() and not line.lstrip().startswith('#'))
        for r in rows:
            applies = (r.get('applies_to') or 'all').strip()
            out.setdefault(r['telescope'].strip(), []).append((r['first_night'].strip(), applies))
    for v in out.values():
        v.sort()
    return out


class Epochs:
    """Camera runs of one telescope from its nights' detectors, plus the manual breaks."""

    def __init__(self, night_detectors, breaks=()):
        nights = sorted((cr.night_of(n), d) for n, d in night_detectors if d and cr.night_of(n))
        self.starts = []          # (first night of the run, detector)
        for night, det in nights:
            if not self.starts or self.starts[-1][1] != det:
                self.starts.append((night, det))
        self.breaks = list(breaks)
        self._cache = {}

    def run(self, night):
        night = cr.night_of(night)
        idx = None
        for i, (first, _) in enumerate(self.starts):
            if first <= night:
                idx = i
            else:
                break
        return idx

    def epoch(self, night, kind):
        key = (night, kind)
        if key not in self._cache:
            n = cr.night_of(night)
            n_breaks = sum(1 for first, applies in self.breaks if first <= n and applies in ('all', kind))
            self._cache[key] = (self.run(n), n_breaks)
        return self._cache[key]

    def describe(self, night, kind):
        r = self.run(night)
        first = self.starts[r][0] if r is not None else None
        det = self.starts[r][1] if r is not None else None
        last_break = max([b for b, a in self.breaks if b <= cr.night_of(night) and a in ('all', kind)] or [None],
                         key=lambda x: x or '')
        return {'detector': det, 'run_start': first, 'last_break': last_break}


# ------------------------------------------------------------------------------------- the science frames

def science_profile(science, readmodes=None):
    """What the calibrations must match: the science frames' majority configuration, its detector, filters,
    exposure keys, median CCD temperature, instrument and readout mode."""
    good = [s for s in science if not s.get('error') and s.get('config')]
    if not good:
        return None
    counts = {}
    for s in good:
        counts[s['config']] = counts.get(s['config'], 0) + 1
    config = max(sorted(counts), key=lambda c: counts[c])
    mine = [s for s in good if s['config'] == config]
    temps = [s['ccd_temp'] for s in mine if s.get('ccd_temp') is not None]
    modes = {}
    for m in (readmodes or []):
        modes[m] = modes.get(m, 0) + 1
    readmode = max(sorted(modes), key=lambda m: modes[m]) if modes else None
    filters = {}
    for s in mine:
        filters[s.get('filter') or ''] = filters.get(s.get('filter') or '', 0) + 1
    exposures = {}
    for s in mine:
        e = cr.exposure_key(s.get('exptime'))
        if e is not None:
            exposures[e] = exposures.get(e, 0) + 1
    return {
        'config': config, 'detector': mine[0].get('detector'), 'instrument': mine[0].get('instrument'),
        'telescop': mine[0].get('telescop'), 'n': len(mine), 'filters': filters, 'exposures': exposures,
        'ccd_temp': float(np.median(temps)) if temps else None, 'readmode': readmode,
        'readmodes': modes, 'other_configs': {c: n for c, n in counts.items() if c != config},
    }


# ------------------------------------------------------------------------------------- choosing

def _delta(night, science_night):
    return cr.night_days(night) - cr.night_days(science_night)


def _order(delta, n, night_dir):
    """Nearest first; the past before the future at the same distance; more frames; then the name."""
    return (abs(delta), 0 if delta <= 0 else 1, -n, night_dir)


def group_by_night(rows):
    out = {}
    for r in rows:
        out.setdefault(r['night_dir'], []).append(r)
    return out


def median_temp(rows):
    t = [r['ccd_temp'] for r in rows if r.get('ccd_temp') is not None]
    return float(np.median(t)) if t else None


class Chooser:
    """Selection for one science night.  `rows` maps kind -> register rows of the telescope (already limited to
    the science configuration), with the science night's own rows from T1 in place of any registered ones."""

    def __init__(self, telescope, night_dir, profile, rows, epochs, resolve, dark_matching=False, log=print):
        self.telescope = telescope
        self.night_dir = night_dir
        self.night = cr.night_of(night_dir)
        self.profile = profile
        self.rows = rows
        self.epochs = epochs
        self.resolve = resolve
        self.dark_matching = dark_matching
        self.log = log
        self.flags = []
        self.errors = []
        self.stale = set()
        self._by_night = {}

    def _groups(self, kind):
        """{night_dir: rows} of a kind, built once."""
        if kind not in self._by_night:
            g = {}
            for r in self.rows.get(kind, []):
                g.setdefault(r['night_dir'], []).append(r)
            self._by_night[kind] = g
        return self._by_night[kind]

    def nights_in_order(self, kind):
        """Night directories of the science night's epoch for this kind, nearest first, the past before the
        future at the same distance."""
        want = self.epochs.epoch(self.night, kind)
        out = []
        for nd, rows in self._groups(kind).items():
            if nd != self.night_dir and not cr.NIGHT_DIR.match(nd):
                continue  # a working copy such as 20251124_rotated is only used for itself
            if self.epochs.epoch(rows[0]['night'], kind) != want:
                continue
            d = _delta(rows[0]['night'], self.night)
            out.append((abs(d), 0 if d <= 0 else 1, nd))
        return [nd for _, _, nd in sorted(out)]

    def usable(self, kind, nd, extra=None):
        """The night's rows of this kind that are compatible with the science frames."""
        out = []
        for r in self._groups(kind).get(nd, []):
            if not r.get('usable') or r.get('config') != self.profile['config']:
                continue
            if kind in ('bias', 'dark') and self.profile.get('readmode') and r.get('readmode') != self.profile[
                    'readmode']:
                continue
            if extra is not None and not extra(r):
                continue
            out.append(r)
        return out

    def candidates(self, kind, extra=None):
        """{night_dir: [usable rows]} compatible with the science frames for this kind, nearest night first."""
        out = {}
        for nd in self.nights_in_order(kind):
            members = self.usable(kind, nd, extra)
            if members:
                out[nd] = members
        return out

    def _accept(self, nd, members, minimum, check):
        if len(members) < minimum:
            return None
        if check is not None and not check(nd, members):
            return None
        paths = [self.resolve(m) for m in members]
        ok = [(m, p) for m, p in zip(members, paths) if p]
        if len(ok) < minimum:
            self.stale.add(nd)
            return None
        return nd, ok

    def nearest(self, kind, minimum, extra=None, check=None):
        """(night_dir, [(row, path)]) of the nearest night with `minimum` compatible frames on disk, or None."""
        for nd in self.nights_in_order(kind):
            got = self._accept(nd, self.usable(kind, nd, extra), minimum, check)
            if got is not None:
                return got
        return None

    def pick(self, groups, minimum, check=None):
        """(night_dir, rows) of the best group with at least `minimum` frames that resolve on disk, or None."""
        order = sorted(groups.items(), key=lambda kv: _order(_delta(kv[1][0]['night'], self.night),
                                                              len(kv[1]), kv[0]))
        for nd, members in order:
            got = self._accept(nd, members, minimum, check)
            if got is not None:
                return got
        return None

    def entry(self, nd, ok, **extra):
        rows = [m for m, _ in ok]
        e = {'night': nd, 'delta_days': _delta(rows[0]['night'], self.night), 'n': len(ok),
             'own_night': nd == self.night_dir, 'frames': [p for _, p in ok]}
        e.update(extra)
        return e

    # ---- the types

    def min_bias(self):
        return MIN_BIAS.get(self.profile.get('instrument'), MIN_BIAS_DEFAULT)

    def bias(self):
        got = self.nearest('bias', self.min_bias())
        if got is None:
            self.errors.append('no night in this epoch has {} usable {} bias frames in readout mode {}'.format(
                self.min_bias(), self.profile['config'], self.profile.get('readmode')))
            return None
        nd, ok = got
        return self.entry(nd, ok, readmode=ok[0][0].get('readmode'))

    def _temp_ok(self, tol):
        want = self.profile.get('ccd_temp')

        def check(nd, members):
            if want is None or tol is None:
                return True
            t = median_temp(members)
            return t is None or abs(t - want) <= tol
        return check

    def dark_tolerance(self):
        return DARK_TEMP_TOLERANCE.get(self.profile.get('instrument'), DARK_TEMP_TOLERANCE_DEFAULT)

    def darks(self):
        if self.dark_matching:
            return self.darks_by_exposure()
        tol = self.dark_tolerance()
        got = self.nearest('dark', MIN_DARK, check=self._temp_ok(tol))
        relaxed = False
        if got is None and self.profile.get('instrument') in DARK_TEMP_RELAX:
            got = self.nearest('dark', MIN_DARK)
            relaxed = got is not None
        if got is None:
            self.errors.append('no night in this epoch has {} usable {} dark frames in readout mode {}'.format(
                MIN_DARK, self.profile['config'], self.profile.get('readmode')))
            return None
        nd, ok = got
        e = self.entry(nd, ok, ccd_temp=median_temp([m for m, _ in ok]))
        if relaxed:
            self.flags.append('darks from {} at {:.1f} C, outside {:.0f} C of the science frames ({:.1f} C)'.format(
                nd, e['ccd_temp'], tol, self.profile['ccd_temp']))
        return e

    def darks_by_exposure(self):
        tol = self.dark_tolerance()
        exposures = sorted(self.profile.get('exposures') or {})
        temp_ok = self._temp_ok(tol)
        groups = self.candidates('dark')
        # {night_dir: {exposure: rows}}
        by_night = {}
        for nd, members in groups.items():
            for m in members:
                by_night.setdefault(nd, {}).setdefault(cr.exposure_key(m.get('exptime')), []).append(m)
        best = {}
        for e in exposures:
            got = self.pick({nd: v[e] for nd, v in by_night.items() if e in v}, MIN_DARK, temp_ok)
            if got is not None:
                best[e] = got
        # one night with darks at every exposure that has any (an exposure no night has does not count)
        covering = None
        matched = sorted(best)
        cover_groups = {nd: sum((v[e] for e in matched), []) for nd, v in by_night.items()
                        if len(matched) > 1 and all(len(v.get(e, [])) >= MIN_DARK for e in matched)}

        def covers(nd, members):
            if not temp_ok(nd, members):
                return False
            return all(sum(1 for m in members if cr.exposure_key(m.get('exptime')) == e and self.resolve(m))
                       >= MIN_DARK for e in matched)
        got = self.pick(cover_groups, MIN_DARK * len(matched), covers) if cover_groups else None
        if got is not None:
            nd, ok = got
            d_cover = abs(_delta(ok[0][0]['night'], self.night))
            if all(d_cover <= abs(_delta(b[1][0][0]['night'], self.night)) + COVERING_SLACK_DAYS
                   for b in best.values()):
                covering = nd
                for e in matched:
                    best[e] = (nd, [(m, p) for m, p in ok if cr.exposure_key(m.get('exptime')) == e])
        frames, by_exp = [], {}
        for e in exposures:
            if e not in best:
                continue
            nd, ok = best[e]
            by_exp[str(e)] = {'night': nd, 'delta_days': _delta(ok[0][0]['night'], self.night), 'n': len(ok),
                              'ccd_temp': median_temp([m for m, _ in ok])}
            frames += [p for _, p in ok]
        missing = [e for e in exposures if e not in best]
        if missing and not best:
            # no exposure matched: the nearest temperature-matched night with any darks, for the scaled fallback
            got = self.pick(groups, MIN_DARK, temp_ok)
            if got is None:
                self.errors.append('no night in this epoch has {} usable {} dark frames within {:.0f} C of the '
                                   'science frames'.format(MIN_DARK, self.profile['config'], tol))
                return None
            nd, ok = got
            frames = [p for _, p in ok]
            by_exp['any'] = {'night': nd, 'delta_days': _delta(ok[0][0]['night'], self.night), 'n': len(ok),
                             'ccd_temp': median_temp([m for m, _ in ok]),
                             'exposures': sorted({cr.exposure_key(m.get('exptime')) for m, _ in ok})}
        # the primary night's darks at the exposures no science frame needs are kept too, as the old lists kept
        # all of a night's darks: MasterDark.fits, which pipeflat.py scales to every flat and pipered.py to an
        # unmatched exposure, stays the night's longest (best-measured) master.  An exposure already taken from
        # another night is not added, so no master combines two nights.
        extra = []
        if best:
            primary = covering or best[max(sorted(best), key=lambda e: self.profile['exposures'][e])][0]
            for e2, members in sorted(by_night.get(primary, {}).items()):
                if e2 in best:
                    continue
                paths = [self.resolve(m) for m in members]
                extra += [p for p in paths if p]
                if any(paths):
                    by_exp.setdefault('extra', {'night': primary, 'exposures': []})['exposures'].append(e2)
            frames += extra
        if missing:
            self.flags.append('no darks at {} s within {:.0f} C in this epoch: T6 subtracts the highest-exposure '
                              'master dark scaled to those exposures'.format(
                                  ', '.join(str(e) for e in missing), tol))
        return {'by_exposure': by_exp, 'covering_night': covering, 'missing_exposures': missing,
                'frames': sorted(set(frames)), 'n': len(set(frames))}

    def flats(self):
        out = {}
        filters = sorted(self.profile.get('filters') or {})
        for f in filters:
            got = self.nearest('flat', MIN_FLAT, extra=lambda r, f=f: (r.get('filter') or '') == f)
            if got is None:
                self.flags.append('no night in this epoch has {} usable {} flats in {}: its {} science frame(s) '
                                  'will not be reduced'.format(MIN_FLAT, self.profile['config'], f or "''",
                                                               self.profile['filters'][f]))
                continue
            nd, ok = got
            out[f] = self.entry(nd, ok, readmodes=sorted({m.get('readmode') or '' for m, _ in ok}))
        if filters and not out:
            self.errors.append('no science filter ({}) has flats in this epoch'.format(', '.join(filters)))
        return out


def resolver(images_root):
    """Path of a register row's frame on disk: its recorded name, else the same stem plain or Rice-compressed
    (raw frames are packed to X.fits.fz after a night's run since v5.0.0, and T2-T6 read either), or None."""
    def resolve(row):
        base = os.path.join(images_root, row['night_dir'])
        p = os.path.join(base, row['relpath'])
        if os.path.exists(p):
            return p
        for s in ('.fits', '.fts', '.fit', '.fits.fz', '.fts.fz', '.fit.fz'):
            q = os.path.join(base, row['stem'] + s)
            if os.path.exists(q):
                return q
        return None
    return resolve


def select(register, telescope, night_dir, images_root, tonight_rows, science, science_readmodes,
           dark_matching=False, tonight_detector=None, breaks=None, log=print):
    """Choose the night's bias, dark and flat frames.  Returns the selection record (see the module docstring);
    record['errors'] non-empty means the night cannot be calibrated."""
    profile = science_profile(science, science_readmodes)
    record = {'telescope': telescope, 'night': night_dir, 'register': getattr(register, 'path', None),
              'rules_version': cr.RULES_VERSION,
              'minimum': {'bias': MIN_BIAS.get((science_profile(science) or {}).get('instrument'), MIN_BIAS_DEFAULT),
                          'dark': MIN_DARK, 'flat': MIN_FLAT},
              'profile': profile, 'flags': [], 'errors': [], 'stale_nights': []}
    if profile is None:
        record['errors'].append('no readable science frames')
        return record
    if profile['other_configs']:
        record['flags'].append('science frames in other configurations get no calibrations: {}'.format(
            ', '.join('{} ({})'.format(c, n) for c, n in sorted(profile['other_configs'].items()))))
    if len(profile.get('readmodes') or {}) > 1:
        record['flags'].append('science frames in several readout modes ({}): bias and darks match {}'.format(
            ', '.join(sorted(profile['readmodes'])), profile['readmode']))
    nights = [(n['night_dir'], n['detector']) for n in (register.nights(telescope) if register else [])
              if n['night_dir'] != night_dir]
    nights.append((night_dir, tonight_detector or profile.get('detector')))
    if breaks is None:
        breaks = load_breaks().get(telescope, [])
    epochs = Epochs(nights, breaks)
    record['epoch'] = {k: epochs.describe(night_dir, k) for k in cr.KINDS}
    rows = {}
    for kind in cr.KINDS:
        reg_rows = []
        if register is not None:
            reg_rows = [dict(r) for r in register.conn.execute(
                'SELECT * FROM frames WHERE telescope = ? AND kind = ? AND config = ? AND usable = 1 '
                'AND night_dir != ?', (telescope, kind, profile['config'], night_dir))]
        rows[kind] = reg_rows + [r for r in tonight_rows if r.get('kind') == kind]
    chooser = Chooser(telescope, night_dir, profile, rows, epochs, resolver(images_root), dark_matching, log)
    record['bias'] = chooser.bias()
    record['dark'] = chooser.darks()
    record['flat'] = chooser.flats()
    record['flags'] += chooser.flags
    record['flags'] += far_flags(record)
    record['errors'] += chooser.errors
    record['stale_nights'] = sorted(chooser.stale)
    if chooser.stale:
        record['flags'].append('register rows of {} no longer match the files on disk (re-register them with '
                               'calregister night)'.format(', '.join(sorted(chooser.stale))))
    return record


def far_flags(record):
    """A flag for every chosen night more than FAR_DAYS away."""
    out = []

    def check(label, e):
        if e and 'delta_days' in e and abs(e['delta_days']) > FAR_DAYS:
            out.append('{} from {} ({:+d} d), more than {} days away'.format(label, e['night'], e['delta_days'],
                                                                         FAR_DAYS))
    check('bias', record.get('bias'))
    d = record.get('dark') or {}
    if 'by_exposure' in d:
        for k, e in sorted(d['by_exposure'].items()):
            if k != 'extra':
                check('dark {}{}'.format(k, '' if k == 'any' else ' s'), e)
    else:
        check('dark', d)
    for f, e in sorted((record.get('flat') or {}).items()):
        check('flat ' + f, e)
    return out


def describe(record):
    """Log lines: every chosen night with |delta days|."""
    lines = []

    def one(label, e):
        if not e:
            return
        lines.append('  {:<14s} {} ({}{} d) {} frames{}'.format(
            label, e['night'], '+' if e['delta_days'] > 0 else '', e['delta_days'], e['n'],
            '' if e.get('own_night', e['delta_days'] == 0) else '  [borrowed, |delta| = {} d]'.format(
                abs(e['delta_days']))))
    one('bias', record.get('bias'))
    d = record.get('dark')
    if d and 'by_exposure' in d:
        for k, e in sorted(d['by_exposure'].items(), key=lambda kv: (kv[0] == 'any', kv[0].zfill(6))):
            if k == 'extra':
                continue
            one('dark {}{}'.format(k, '' if k == 'any' else ' s'), e)
        if d.get('covering_night'):
            lines.append('  dark: one night covers every exposure: {}'.format(d['covering_night']))
        if 'extra' in d['by_exposure']:
            x = d['by_exposure']['extra']
            lines.append('  dark: also {}\'s darks at {} s (no science frame needs them; the longest becomes '
                         'MasterDark.fits)'.format(x['night'], ', '.join(str(e) for e in x['exposures'])))
    else:
        one('dark', d)
    for f, e in sorted((record.get('flat') or {}).items()):
        one('flat {}'.format(f), e)
    for f in record.get('flags', []):
        lines.append('  WARNING: ' + f)
    for e in record.get('errors', []):
        lines.append('  ERROR: ' + e)
    return lines


def write_record(record, reddir, run):
    path = os.path.join(reddir, SELECTION_NAME.format(run=run))
    slim = json.loads(json.dumps(record, default=str))
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(slim, f, indent=1, sort_keys=True)
    os.replace(tmp, path)
    return path


def read_record(reddir, run):
    path = os.path.join(reddir, SELECTION_NAME.format(run=run))
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def report_line(record):
    """One line for the PDF report: where each master came from."""
    if not record:
        return None
    parts = []

    def fmt(label, e):
        if not e:
            return '{} none'.format(label)
        d = e['delta_days']
        return '{} {}{}'.format(label, e['night'], '' if d == 0 else ' ({:+d} d)'.format(d))
    parts.append(fmt('bias', record.get('bias')))
    d = record.get('dark')
    if d and 'by_exposure' in d:
        chosen = [e for k, e in d['by_exposure'].items() if k != 'extra']
        nights = sorted({e['night'] for e in chosen})
        ds = [e['delta_days'] for e in chosen]
        worst = max(ds, key=abs) if ds else 0
        parts.append('dark {}{}'.format('/'.join(nights) or 'none',
                                        '' if worst == 0 else ' (max {:+d} d)'.format(worst)))
    else:
        parts.append(fmt('dark', d))
    for f, e in sorted((record.get('flat') or {}).items()):
        parts.append(fmt('flat ' + f, e))
    line = 'Calibrations: ' + '; '.join(parts)
    n_flags = len(record.get('flags') or [])
    if n_flags:
        line += '  [{} warning{}]'.format(n_flags, '' if n_flags == 1 else 's')
    return line
