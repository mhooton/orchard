"""
What the calibration register would have chosen for nights already reduced, against what T1's old night-by-night
search put in their lists (read-only: reads PipelineOutput/<v>/<Tel>/output/<night>/reduction/*_{bias,dark,flat_*}
.list and the register, writes a report).

For each night and type the old lists give the night(s) the frames came from; the register gives the night
calselect.py chooses for the same science frames.  The science frames are described from the register (the
configuration, detector and temperature of the night's sampled science frames; the readout mode of its own
calibration frames, else the dated table) and, for SPIRIT, ESO's per-night science exposure times (eso-science).
A filter counts as a science filter when its old list borrowed frames (the old search only borrowed for science
filters) or holds at least MIN_FLAT of the night's own flats; a filter with fewer own flats and nothing borrowed
was not a science filter, and keeps its own flats in both schemes.

    python -m calibration.calcompare eso-science --out spirit_science.csv      (Callisto SPIRIT nights, from ESO)
    python -m calibration.calcompare compare --out compare_dir [--science spirit_science.csv]
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

from calibration import calregister as cr
from calibration import calselect as cs

LIST_NAME = re.compile(r'^(?P<run>[^_]+)_(?P<kind>bias|dark|flat_(?P<filter>.+))\.list$')


def source(path):
    """(night directory, path in it) of a frame named in an old list, whatever prefix the list used."""
    p = path.strip().replace('\\', '/')
    i = p.find('/images/')
    if i < 0:
        return None, None
    rest = p[i + len('/images/'):]
    if '/' not in rest:
        return None, None
    nd, rel = rest.split('/', 1)
    return nd, rel


def read_lists(reddir):
    """{'bias': [(night_dir, rel)], 'dark': [...], 'flat': {filter: [...]}} of one reduction directory."""
    out = {'bias': None, 'dark': None, 'flat': {}}
    try:
        names = os.listdir(reddir)
    except OSError:
        return None
    found = False
    for name in sorted(names):
        m = LIST_NAME.match(name)
        if not m or m.group('run') != '1':
            continue
        found = True
        try:
            with open(os.path.join(reddir, name)) as f:
                frames = [source(line) for line in f if line.strip()]
        except OSError:
            continue
        frames = [x for x in frames if x[0]]
        if m.group('kind') in ('bias', 'dark'):
            out[m.group('kind')] = frames
        else:
            out['flat'][m.group('filter')] = frames
    return out if found else None


class TelescopeView:
    """Register rows of one telescope, grouped for fast lookups."""

    def __init__(self, register, telescope):
        self.telescope = telescope
        self.nights = {n['night_dir']: n for n in register.nights(telescope)}
        self.epochs = cs.Epochs([(n['night_dir'], n['detector']) for n in self.nights.values()],
                                cs.load_breaks().get(telescope, []))
        self.by_key = {}
        self.by_config = collections.defaultdict(lambda: collections.defaultdict(list))
        self.instrument = {}
        for r in register.frames(telescope):
            self.by_key[(r['night_dir'], r['stem'])] = r
            if r['usable']:
                self.by_config[r['config']][r['kind']].append(r)
            self.instrument.setdefault(r['config'], r.get('instrument'))
            self.instrument.setdefault(('telescop', r['config']), r.get('telescop'))

    def row(self, night_dir, rel):
        return self.by_key.get((night_dir, cr.strip_suffix(rel)))

    def profile(self, night_dir, science_exposures=None):
        n = self.nights.get(night_dir)
        if not n:
            return None
        sci = json.loads(n.get('science') or '{}')
        if not sci:
            return None
        config = max(sorted(sci), key=lambda c: sci[c]['n'])
        c = sci[config]
        own = [r for r in self.by_key.values() if r['night_dir'] == night_dir and r['config'] == config]
        modes = collections.Counter(r.get('readmode') for r in own if r['kind'] == 'bias' and r.get('readmode'))
        if not modes:
            modes = collections.Counter(r.get('readmode') for r in own if r.get('readmode'))
        instrument = self.instrument.get(config)
        if modes:
            readmode = modes.most_common(1)[0][0]
        elif instrument == 'andor':
            from calibration import readout
            reg = readout.expected_regime(self.instrument.get(('telescop', config)), cr.night_of(night_dir))
            readmode = reg.mode if reg else 'NORMAL'
        else:
            readmode = 'DEFAULT'
        exposures = {}
        if science_exposures is not None:
            exposures = {int(e): k for e, k in science_exposures.items()}
        return {'config': config, 'detector': c.get('detector'), 'instrument': instrument,
                'telescop': self.instrument.get(('telescop', config)), 'n': c['n'], 'filters': {},
                'exposures': exposures, 'ccd_temp': c.get('ccd_temp'), 'readmode': readmode,
                'readmodes': dict(modes), 'other_configs': {}}


def dark_matching(view, config):
    telescop = view.instrument.get(('telescop', config))
    instrument = view.instrument.get(config)
    params = cr._instrument_config().get('telescopes', {}).get(telescop, {}).get(instrument, {})
    return bool(params.get('dark_matching_exptime', False))


def nights_of(frames):
    return collections.Counter(nd for nd, _ in frames)


def classify(old_nights, night_dir, new, delta_old=None):
    """Outcome of one type: same / different (with why)."""
    if new is None:
        return 'new: none' if old_nights else 'none either'
    if not old_nights:
        return 'old: none'
    if len(old_nights) > 1:
        return 'old: several nights'
    (old,), = [list(old_nights)]
    if old == new['night']:
        return 'same'
    d_new = abs(new['delta_days'])
    d_old = abs(cr.night_days(old) - cr.night_days(night_dir))
    if old == night_dir:
        return 'old own night, new borrows'
    if new['night'] == night_dir:
        return 'old borrowed, new own night'
    if d_new < d_old:
        return 'new nearer'
    if d_new > d_old:
        return 'new farther'
    return 'other side, same distance'


def compare_night(view, night_dir, lists, science_exposures=None):
    profile = view.profile(night_dir, science_exposures)
    out = {'telescope': view.telescope, 'night': night_dir, 'profile': None, 'types': {}, 'cross': [],
           'unregistered': 0, 'old_unusable': 0}
    if profile is None:
        out['skip'] = 'no science frames registered'
        return out
    out['profile'] = {k: profile[k] for k in ('config', 'detector', 'instrument', 'readmode', 'ccd_temp')}
    # science filters, as the old lists show them
    for f, frames in lists['flat'].items():
        own = [x for x in frames if x[0] == night_dir]
        if len(own) < len(frames) or len(own) >= cs.MIN_FLAT:
            profile['filters'][f] = 1
    rows = {k: list(v) for k, v in view.by_config[profile['config']].items()}
    matching = dark_matching(view, profile['config'])
    ch = cs.Chooser(view.telescope, night_dir, profile, rows, view.epochs,
                    lambda r: '{}/{}'.format(r['night_dir'], r['relpath']), matching, log=lambda *a: None)
    new = {'bias': ch.bias(), 'dark': ch.darks(), 'flat': ch.flats()}
    out['flags'], out['errors'] = ch.flags, ch.errors
    # cross-camera and unusable frames in the old lists
    old_all = [('bias', None, x) for x in (lists['bias'] or [])] + [('dark', None, x) for x in (lists['dark'] or [])]
    for f, frames in lists['flat'].items():
        old_all += [('flat', f, x) for x in frames]
    cross = collections.Counter()
    for kind, f, (nd, rel) in old_all:
        r = view.row(nd, rel)
        if r is None:
            out['unregistered'] += 1
            continue
        if not r['usable']:
            out['old_unusable'] += 1
        if r.get('detector') != profile['detector'] or r.get('config') != profile['config']:
            cross[(kind, f, nd, r.get('detector'), r.get('config'))] += 1
    out['cross'] = [{'kind': k, 'filter': f, 'from_night': nd, 'detector': det, 'config': cfg, 'n': n,
                     'other_camera': det != profile['detector']}
                    for (k, f, nd, det, cfg), n in sorted(cross.items(), key=lambda kv: str(kv[0]))]
    # outcomes
    out['types']['bias'] = {'old': dict(nights_of(lists['bias'] or [])), 'new': _brief(new['bias']),
                            'outcome': classify(nights_of(lists['bias'] or []), night_dir, new['bias'])}
    if matching:
        old_by_exp = collections.defaultdict(list)
        for nd, rel in lists['dark'] or []:
            r = view.row(nd, rel)
            if r is not None:
                old_by_exp[cr.exposure_key(r.get('exptime'))].append((nd, rel))
        per = {}
        newd = new['dark'] or {}
        for e in sorted(profile['exposures']):
            ne = (newd.get('by_exposure') or {}).get(str(e))
            if ne is not None:
                ne = dict(ne, own_night=ne['night'] == night_dir)
            old_e = old_by_exp.get(e, [])
            if old_e:
                old_state = 'matched'
            elif lists['dark']:
                old_state = 'combined fallback'
            else:
                old_state = 'no dark'
            per[str(e)] = {'old': dict(nights_of(old_e)), 'old_state': old_state, 'new': _brief(ne),
                           'new_state': 'matched' if ne else ('combined fallback' if newd else 'none'),
                           'science_frames': profile['exposures'][e],
                           'outcome': classify(nights_of(old_e), night_dir, ne) if old_e or ne else 'none either'}
        out['types']['dark'] = {'by_exposure': per, 'covering_night': newd.get('covering_night'),
                                'old': dict(nights_of(lists['dark'] or []))}
    else:
        out['types']['dark'] = {'old': dict(nights_of(lists['dark'] or [])), 'new': _brief(new['dark']),
                                'outcome': classify(nights_of(lists['dark'] or []), night_dir, new['dark'])}
    out['types']['flat'] = {}
    for f in sorted(set(lists['flat']) | set(profile['filters'])):
        if f not in profile['filters']:
            continue
        old_f = nights_of(lists['flat'].get(f, []))
        ne = new['flat'].get(f)
        out['types']['flat'][f] = {'old': dict(old_f), 'new': _brief(ne), 'outcome': classify(old_f, night_dir, ne)}
    return out


def _brief(e):
    if not e:
        return None
    return {'night': e['night'], 'delta_days': e['delta_days'], 'n': e['n']}


def compare(register, basedir, telescopes, versions=('v2', 'v3'), science=None, only=None, log=print):
    results = []
    for tel in telescopes:
        view = TelescopeView(register, tel)
        for v in versions:
            outdir = os.path.join(basedir, 'PipelineOutput', v, tel, 'output')
            try:
                nights = sorted(d for d in os.listdir(outdir) if cr.NIGHT_DIR.match(d) and (not only or d in only))
            except OSError:
                continue
            n = 0
            for nd in nights:
                lists = read_lists(os.path.join(outdir, nd, 'reduction'))
                if lists is None:
                    continue
                res = compare_night(view, nd, lists, (science or {}).get((tel, nd)))
                res['version'] = v
                results.append(res)
                n += 1
            log('{} {}: {} nights compared'.format(tel, v, n))
    return results


def load_science(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r['telescope'], r['night'])
            try:
                e = cr.exposure_key(float(r['exposure']))
            except ValueError:
                continue  # ESO rows without an exposure time
            out.setdefault(key, {})
            out[key][e] = out[key].get(e, 0) + int(r['n'])
    return out


def eso_science(register, out_path, env_path, nights=None, log=print):
    """Per-night science exposure times of the Callisto SPIRIT nights, from ESO (one grouped query a night)."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from jobqueue.eso import EsoArchive
    from jobqueue.util import night_bounds, to_mjd
    archive = EsoArchive(env_path, log=log)
    if nights is None:
        nights = [n['night_dir'] for n in register.nights('Callisto') if n.get('detector') == '12.0um'
                  and cr.NIGHT_DIR.match(n['night_dir'])]
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['telescope', 'night', 'object', 'filter', 'exposure', 'n'])
        for i, nd in enumerate(nights):
            a, b = (to_mjd(t) for t in night_bounds(nd))
            q = ("SELECT object, filter_path, exposure, COUNT(*) AS n FROM dbo.raw WHERE prog_id = '60.A-9009(D)' "
                 "AND dp_type = 'OBJECT' AND mjd_obs >= {:.6f} AND mjd_obs < {:.6f} "
                 "GROUP BY object, filter_path, exposure").format(a, b)
            for r in archive.tap(q):
                w.writerow(['Callisto', nd, r['object'], r['filter_path'], r['exposure'], r['n']])
            if i % 50 == 0:
                log('{}/{} nights'.format(i, len(nights)))
    log('wrote {}'.format(out_path))


def summarise(results):
    """Counts for the write-up."""
    s = {'nights': collections.Counter(), 'bias': collections.Counter(), 'dark': collections.Counter(),
         'flat': collections.Counter(), 'spirit_dark_old': collections.Counter(),
         'spirit_dark_old_frames': collections.Counter(), 'spirit_dark_new': collections.Counter(),
         'cross_camera': [], 'cross_config': [], 'errors': collections.Counter(), 'flags': collections.Counter()}
    for r in results:
        key = (r['version'], r['telescope'])
        s['nights'][key] += 1
        if r.get('skip'):
            s['nights'][key + ('skipped',)] += 1
            continue
        for e in r.get('errors') or []:
            s['errors'][re.sub(r'\d+', '#', e)[:120]] += 1
        for fl in r.get('flags') or []:
            s['flags'][re.sub(r'\d+', '#', fl)[:120]] += 1
        t = r['types']
        s['bias'][(r['version'], r['telescope'], t['bias']['outcome'])] += 1
        if 'by_exposure' in t['dark']:
            for e, x in t['dark']['by_exposure'].items():
                s['dark'][(r['version'], r['telescope'], 'spirit ' + x['outcome'])] += 1
                s['spirit_dark_old'][(r['version'], x['old_state'])] += 1
                s['spirit_dark_old_frames'][(r['version'], x['old_state'])] += x['science_frames']
                s['spirit_dark_new'][(r['version'], x['new_state'])] += 1
        else:
            s['dark'][(r['version'], r['telescope'], t['dark']['outcome'])] += 1
        for f, x in t['flat'].items():
            s['flat'][(r['version'], r['telescope'], x['outcome'])] += 1
        for c in r['cross']:
            (s['cross_camera'] if c['other_camera'] else s['cross_config']).append(
                dict(c, version=r['version'], telescope=r['telescope'], night=r['night'],
                     science=r['profile']['config']))
    return s


def write_report(results, outdir):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'compare.jsonl'), 'w') as f:
        for r in results:
            f.write(json.dumps(r, sort_keys=True, default=str) + '\n')
    s = summarise(results)
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump({k: ({' | '.join(map(str, kk)) if isinstance(kk, tuple) else kk: vv for kk, vv in v.items()}
                       if isinstance(v, collections.Counter) else v) for k, v in s.items()},
                  f, indent=1, sort_keys=True, default=str)
    # the nights whose calibrations would change, one row per night and type
    with open(os.path.join(outdir, 'changes.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['version', 'telescope', 'night', 'type', 'outcome', 'old_nights', 'new_night', 'new_delta_days'])
        for r in results:
            if r.get('skip'):
                continue
            rows = [('bias', r['types']['bias'])]
            d = r['types']['dark']
            rows += [('dark {} s'.format(e), x) for e, x in sorted(d.get('by_exposure', {}).items())] \
                if 'by_exposure' in d else [('dark', d)]
            rows += [('flat ' + f, x) for f, x in sorted(r['types']['flat'].items())]
            for name, x in rows:
                if x['outcome'] in ('same', 'none either'):
                    continue
                w.writerow([r['version'], r['telescope'], r['night'], name, x['outcome'],
                            ' '.join('{}x{}'.format(k, v) for k, v in sorted(x['old'].items())),
                            (x['new'] or {}).get('night', ''), (x['new'] or {}).get('delta_days', '')])
    return s


def main(argv=None):
    p = argparse.ArgumentParser(prog='python -m calibration.calcompare')
    p.add_argument('--register', required=True)
    p.add_argument('--basedir', default='/data/SPECULOOSPipeline')
    sub = p.add_subparsers(dest='cmd', required=True)
    e = sub.add_parser('eso-science')
    e.add_argument('--out', required=True)
    e.add_argument('--env', default='/opt/orchard/src/reporting/.env')
    e.add_argument('--nights-from', help='CSV with night and shapes columns (SPIRIT nights: shapes 1024x1280)')
    c = sub.add_parser('compare')
    c.add_argument('--out', required=True)
    c.add_argument('--science')
    c.add_argument('--telescopes', default='Io,Europa,Ganymede,Callisto,Artemis')
    c.add_argument('--nights', help='comma-separated night directories (default: all)')
    a = p.parse_args(argv)
    log = lambda m: print(m, flush=True)
    if a.cmd == 'eso-science' and a.nights_from:
        with open(a.nights_from) as f:
            nights = [r['night'] for r in csv.DictReader(f) if '1024x1280' in r.get('shapes', '')
                      or '1280x1024' in r.get('shapes', '')]
        eso_science(None, a.out, a.env, nights=nights, log=log)
        return 0
    reg = cr.Register(a.register, readonly=True)
    try:
        if a.cmd == 'eso-science':
            eso_science(reg, a.out, a.env, log=log)
        else:
            res = compare(reg, a.basedir, a.telescopes.split(','), science=load_science(a.science),
                          only=set(a.nights.split(',')) if a.nights else None, log=log)
            s = write_report(res, a.out)
            log(json.dumps({k: len(v) if isinstance(v, list) else sum(v.values()) for k, v in s.items()}))
    finally:
        reg.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
