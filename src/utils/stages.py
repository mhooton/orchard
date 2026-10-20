#!/usr/bin/env python3
"""
stages.py
=========

The stage map: what each stage of main/ZLP_pipeline.sh reads and writes, and so which stages a fix has to rerun.

A fix to a late stage need not cost a rerun from raw. When the PWV correction failed (the LHATPRO fetch, during
appct's DNS outage of 2026-10-18..20), T10 made light curves without it and T11 put them in the night's PDF; T1-T9
were fine. So "PWV failed" means: rerun from T10, which runs T10 and T11 on what T1-T9 left. `cause_stage('pwv')`
says so, `rerun_plan(10)` lists what such a rerun reads and rewrites, and utils/restage.py brings exactly those
inputs back into v3 for a night that promotion has already moved into v2 (docs/rerun-from-stage.md).

Paths are patterns relative to a scope:

    raw        Observations/<TEL>/images/<NIGHT>/
    obs        Observations/<TEL>/                  (schedule/, technical_logs/)
    reduction  PipelineOutput/v3/<TEL>/output/<NIGHT>/reduction/
    target     PipelineOutput/v3/<TEL>/output/<NIGHT>/<TARGET>/
    telescope  PipelineOutput/v3/<TEL>/              (output/output.fits, reports/, logs/)
    data       the data root                        (target list, PWV grids, Teff list, v2/StackImages)

A pattern matches path components one by one: `*` never crosses a '/', `**` matches any number of components, and
{night}, {target} and {day} (YYYY-MM-DD) are filled in. A read says which stage made the file (made_by), so the
inputs of a rerun from stage N are the reads of stages N..11 that stages before N made. A file a later stage edits
in place (the WCS that T7, T8 and T9 write into the calibrated frames) cannot be shared with v2 by a hard link.

    python utils/stages.py                   the map
    python utils/stages.py --cause pwv       which stages a PWV fix reruns, and what that rerun reads and writes
    python utils/stages.py --from 10         the same for a rerun from T10
    python utils/stages.py --markdown        the tables in docs/rerun-from-stage.md
"""
import argparse
import datetime as dt
import fnmatch
import re
import sys
from collections import OrderedDict, namedtuple

VERSION = 1

RAW, OBS, REDUCTION, TARGET, TELESCOPE, DATA = 'raw', 'obs', 'reduction', 'target', 'telescope', 'data'
SCOPES = OrderedDict([
    (RAW, 'Observations/<TEL>/images/<NIGHT>/'),
    (OBS, 'Observations/<TEL>/'),
    (REDUCTION, 'PipelineOutput/v3/<TEL>/output/<NIGHT>/reduction/'),
    (TARGET, 'PipelineOutput/v3/<TEL>/output/<NIGHT>/<TARGET>/'),
    (TELESCOPE, 'PipelineOutput/v3/<TEL>/'),
    (DATA, '<BASEDIR>/'),
])

LAST = 11           # T12 (the v3 -> v2 sweep) is replaced by the promotion step: reruns run --no_T12
RESTAGE = (8, 9, 10, 11)  # stages utils/restage.py can rerun a promoted night from (restage_inputs: why not before)

Read = namedtuple('Read', 'scope pattern made_by what')
Write = namedtuple('Write', 'scope pattern what')
Stage = namedtuple('Stage', 'n name function script reads writes edits note')


def R(scope, pattern, made_by=None, what=''):
    return Read(scope, pattern, made_by, what)


def W(scope, pattern, what=''):
    return Write(scope, pattern, what)


FRAMES = '1/proc*.fits'      # calibrated frames: made by T6, WCS written into them by T7, T8 and T9
FRAME_LISTS = ('1/*_processed.dat',)   # name the frames by absolute path: restaged as copies pointed at the copies

STAGES = OrderedDict((s.n, s) for s in [
    Stage(1, 'lists', 'create_input_lists', 'main/createlists.py',
          [R(RAW, '*', what='every frame header')],
          [W(REDUCTION, '1_bias.list'), W(REDUCTION, '1_dark.list'), W(REDUCTION, '1_flat_*.list'),
           W(REDUCTION, '1_image_*.list', 'one per target: the targets of the night')], [], ''),
    Stage(2, 'master bias', 'create_master_bias', 'calibration/pipebias.py',
          [R(REDUCTION, '1_bias.list', 1), R(RAW, '*', what='bias frames')],
          [W(REDUCTION, '1_MasterBias.fits'), W(REDUCTION, 'overscan.dat'), W(REDUCTION, 'readout.json'),
           W(REDUCTION, 'readoutnoise.dat')], [], ''),
    Stage(3, 'master dark', 'create_master_dark', 'calibration/pipedark.py',
          [R(REDUCTION, '1_dark.list', 1), R(REDUCTION, '1_MasterBias.fits', 2), R(RAW, '*', what='dark frames')],
          [W(REDUCTION, '1_MasterDark*.fits')], [], ''),
    Stage(4, 'master flat', 'create_master_flat', 'calibration/pipeflat.py',
          [R(REDUCTION, '1_flat_*.list', 1), R(REDUCTION, '1_MasterBias.fits', 2),
           R(REDUCTION, '1_MasterDark*.fits', 3), R(RAW, '*', what='flat frames')],
          [W(REDUCTION, '1_MasterFlat_*.fits'), W(REDUCTION, 'datafile.dat'), W(REDUCTION, 'expdata.dat'),
           W(REDUCTION, 'flats/**', 'only with --save-flat-products')], [], ''),
    Stage(5, 'bad pixel map', 'create_bad_pixel_map', 'calibration/pipebadpixel.py',
          [R(REDUCTION, '1_image_{target}.list', 1), R(REDUCTION, '1_Master*.fits', 4)],
          [W(REDUCTION, '1_BadPixelMap.fits'), W(REDUCTION, 'hot_pixel_map.fits'), W(REDUCTION, 'cold_pixel_map.fits'),
           W(REDUCTION, 'flatbad_pixel_map.fits'), W(REDUCTION, 'darkcurrent.dat')], [], ''),
    Stage(6, 'calibrate', 'reduce_science_images', 'calibration/pipered.py',
          [R(REDUCTION, '1_image_{target}.list', 1), R(REDUCTION, '1_Master*.fits', 4),
           R(REDUCTION, '1_BadPixelMap.fits', 5), R(REDUCTION, '*.dat', 5), R(REDUCTION, 'readout.json', 2),
           R(RAW, '*', what='science frames')],
          [W(TARGET, FRAMES, 'calibrated frames'), W(TARGET, '1/*_processed.dat', 'frame list per filter')], [],
          'a target with nothing calibrated skips T7-T10'),
    Stage(7, 'astrometry', 'check_astrometry', 'astrom/astrometry.py',
          [R(TARGET, '1/*_processed.dat', 6), R(TARGET, FRAMES, 6), R(REDUCTION, '1_image_{target}.list', 1),
           R(RAW, '*', what='headers of the frames the image list names'),
           R(DATA, 'PipelineOutput/wcs_cache/{tel}/{night}.ecsv', what="earlier runs' solutions"),
           R(DATA, 'Gaia database', what='$GAIADATABASEPATH')],
          [W(DATA, 'PipelineOutput/wcs_cache/{tel}/{night}.ecsv', 'the WCS cache (PR #64): raw frames are no '
                                                                  'longer written')],
          [R(TARGET, FRAMES, 6, 'WCS')], ''),
    Stage(8, 'stack', 'create_stack_image', 'photometry/catalogue_fov.py',
          [R(TARGET, '1/*_processed.dat', 6), R(TARGET, FRAMES, 6), R(OBS, 'schedule/**', what='plan files'),
           R(DATA, 'ml_40pc_v2.csv', what='the 40 pc target list'),
           R(DATA, 'PipelineOutput/v2/StackImages/**', what='backup catalogues')],
          [W(TARGET, '1/*_stacked.dat'), W(TARGET, '1/*.cat', 'of the stacked frames'),
           W(TARGET, '1_initial-catalogue.fits'), W(TARGET, '*_outstack_*.fits'),
           W(TARGET, '*_stack_catalogue_*'), W(TARGET, '.expected_stack_catalogue'),
           W(TARGET, '.active_stack_catalogue_*'), W(DATA, 'PipelineOutput/v2/StackImages/**', 'backup catalogues')],
          [R(TARGET, FRAMES, 6, 'WCS of the stacked frames (wcsfit)')],
          'a target without a catalogue skips T9-T10'),
    Stage(9, 'photometry', 'perform_aperture_photometry',
          'photometry/ZLP_app_photom.py, condense/condense_photometry.py',
          [R(TARGET, '1/*_processed.dat', 6), R(TARGET, FRAMES, 6), R(TARGET, '.active_stack_catalogue_*', 8),
           R(TARGET, '.expected_stack_catalogue', 8), R(TARGET, '*_stack_catalogue_*.fits', 8),
           R(OBS, 'schedule/**', what="plan files (the target's Gaia ID)"),
           R(DATA, 'ml_40pc_v2.csv', what='the 40 pc target list'),
           R(DATA, 'Gaia database', what="wcsfit's reference stars ($GAIADATABASEPATH)")],
          [W(TARGET, '1/*.cat'), W(TARGET, '1/*.phot'), W(TARGET, '1/*.ell'), W(TARGET, '1/*.png', 'PSF plots'),
           W(TARGET, '1/*_processed.dat_phot'), W(TARGET, '1_initial-catalogue.fits'),
           W(TARGET, '*_stack_catalogue_pm_*.fits'), W(TARGET, '*_{night}_output.fits', 'the night\'s photometry'),
           W(TELESCOPE, 'output/StackImages/already_run.txt', 'targets seen per telescope')],
          [R(TARGET, FRAMES, 6, 'WCS: wcsfit, the stale PC matrix removed, wcscompl'),
           R(TARGET, '*_stack_catalogue_*.fits', 8, "opened read-write by condense (fitsio 'rw')"),
           R(TELESCOPE, 'output/*_output.fits', 9, 'the multi-night file per star: the night replaced')],
          'condense reads TARGET_ROLE, which T8 writes since 2026-03-31: older catalogues need a rerun from T8'),
    Stage(10, 'light curves', 'perform_differential_photometry', 'diff_photometry/SPlightcurve.py',
          [R(TARGET, '*_{night}_output.fits', 9), R(REDUCTION, '1_image_{target}.list', 1, 'its first raw frame'),
           R(RAW, '*', what='the header of the first frame of the image list'),
           R(DATA, 'ml_40pc_v2.csv', what='the 40 pc target list'), R(DATA, 'pwv/PWV_grids/*', what='PWV grids'),
           R(DATA, 'tests/targetListTeffEsts.csv', what='Teff list (Filippazzo)'),
           R(DATA, 'ESO LHATPRO', what='PWV over the network: archive.eso.org')],
          [W(TARGET, '*_diff.fits', 'apertures 3-8'), W(TARGET, '*_pwv_original.fits'),
           W(TARGET, '*_pwv_processed.fits'), W(TARGET, '*_pwv_interpolated.fits'),
           W(TARGET, 'lightcurves/**', 'MCMC files, plots, best aperture')], [],
          'the PWV correction, the Teff it needs and the comparison stars are all here'),
    Stage(11, 'PDF report', 'pdf_report', 'reporting/pdf_report_catriona.py',
          [R(REDUCTION, '1_Master*.fits', 4), R(TARGET, '*_output.fits', 9), R(TARGET, '*_diff.fits', 10),
           R(TARGET, 'lightcurves/**', 10), R(TARGET, FRAMES, 6, 'the stack panel'),
           R(OBS, 'technical_logs/{day}/*'), R(OBS, 'technical_logs/log_dome_closure.txt'),
           R(OBS, 'schedule/Archive_night_blocks/*')],
          [W(TELESCOPE, 'reports/{tel}_{night}.pdf'), W(TELESCOPE, 'output/StackImages_nightly/{night}/*'),
           W(OBS, 'technical_logs/{day}/*_PSF_*.png')], [],
          'then frame packing (16-bit) of the frames, when the run made or touched them'),
    Stage(12, 'sweep', 'perform_version_migration', 'main/ZLP_pipeline.sh',
          [], [], [], 'v3 -> v2; replaced by the promotion step (utils/promote.py): queue runs use --no_T12'),
])

# What has to be rerun when something is fixed: the first stage it affects. Everything after it reruns too, because
# each stage reads what the ones before it wrote.
CAUSES = OrderedDict([
    ('pwv', (10, 'the PWV correction: the LHATPRO fetch, the PWV grids, the correction itself')),
    ('teff', (10, 'the target Teff the PWV correction needs (header TEFF, the Filippazzo list)')),
    ('comparison-stars', (10, 'comparison-star weights and the artificial light curve')),
    ('light-curves', (10, 'anything else in SPlightcurve.py: apertures 3-8, MCMC files, plots, best aperture')),
    ('pdf', (11, 'the night report only')),
    ('stale-pc', (9, 'forced apertures placed with a stale PC matrix over the fitted CD (fixed 2026-10-17, PR #30)')),
    ('photometry', (9, 'aperture photometry, condensing, the night\'s *_output.fits')),
    ('stack', (8, 'the stack, its catalogue, the Gaia crossmatch, which stars are targets')),
    ('astrometry', (7, 'plate solving of the calibrated frames')),
    ('calibration', (6, 'calibration of the science frames (pipered.py)')),
    ('masters', (2, 'master bias, dark or flat, the bad pixel map')),
    ('lists', (1, 'which frame is what (createlists.py)')),
])

FLAG = re.compile(r'^--(no|only|from)_T(\d+)$')


# ---------------------------------------------------------------------------
# Stages a run ran
# ---------------------------------------------------------------------------

def stages_from_argv(argv):
    """The stages a ZLP_pipeline.sh invocation runs, as ZLP_pipeline.sh decides them: every one, less each --no_TN;
    only N with --only_TN; none before N with --from_TN. Returns a sorted list."""
    on, only, first = set(range(1, 13)), None, None
    for a in argv:
        m = FLAG.match(a)
        if not m:
            continue
        kind, n = m.group(1), int(m.group(2))
        if not 1 <= n <= 12:
            continue
        if kind == 'no':
            on.discard(n)
        elif kind == 'only':
            only = n
        else:
            first = n
    if only is not None:
        on = {only}
    if first is not None:
        on -= set(range(1, first))
    return sorted(on)


def parse_stages_line(line):
    """The stages from ZLP_pipeline.sh's 'Stages run: T10 T11' log line, or None if the line is not one."""
    if not line.startswith('Stages run:'):
        return None
    return sorted(int(x[1:]) for x in line[len('Stages run:'):].split() if re.match(r'^T\d+$', x))


def first_stage(stages):
    return min(stages) if stages else None


# ---------------------------------------------------------------------------
# Patterns
# ---------------------------------------------------------------------------

def fill(pattern, tel='*', night='*', target='*'):
    day = '*'
    if night != '*':
        day = dt.datetime.strptime(night, '%Y%m%d').strftime('%Y-%m-%d')
    return pattern.replace('{tel}', tel).replace('{night}', night).replace('{target}', target).replace('{day}', day)


def _match(parts, pats):
    if not pats:
        return not parts
    if pats[0] == '**':
        return any(_match(parts[k:], pats[1:]) for k in range(len(parts) + 1))
    return bool(parts) and fnmatch.fnmatchcase(parts[0], pats[0]) and _match(parts[1:], pats[1:])


def matches(rel, pattern, **kw):
    """Whether a relative path matches a pattern ('*' within a component, '**' across them)."""
    return _match(rel.replace('\\', '/').split('/'), fill(pattern, **kw).split('/'))


def writer(rel, scope=TARGET, **kw):
    """The last stage that writes a file at this path in the scope (T9 rewrites what T8's pattern also covers), or
    None if no stage writes such a file."""
    found = [s.n for s in STAGES.values() if any(w.scope == scope and matches(rel, w.pattern, **kw) for w in s.writes)]
    return max(found) if found else None


def made_before(rel, first, **kw):
    """Whether a file in a target directory is a product of the stages before `first`: an input of a rerun from
    `first`, which must predate it. False for files no stage writes."""
    n = writer(rel, **kw)
    return n is not None and n < first


# ---------------------------------------------------------------------------
# What a rerun from stage N reads and rewrites
# ---------------------------------------------------------------------------

def rerun_stages(first):
    return [n for n in STAGES if first <= n <= LAST]


def restage_inputs(first):
    """What a rerun from `first` needs in the night's reduction and the target directory, made by the stages before
    it: [(scope, pattern, made_by, link)], one per (scope, pattern). link: no stage of the rerun edits the file in
    place, so v2's own file can be shared by a hard link; otherwise it must be a copy. From T10 on, every input can be
    shared. From T9 back, the calibrated frames get new WCS headers (wcsfit writes into the file), so a rerun needs
    its own copy of every frame. A file a stage of the rerun writes is never an input, whatever pattern it matches
    (T9's *_stack_catalogue_pm_*.fits matches what T9 reads of T8's): restage.py leaves those out.
    From T7 back a rerun is not restaged yet: T7 wrote each frame's WCS into the raw frame its image list names until
    PR #64's WCS cache replaced that (2026-10-20), so T7 is the next stage to try; from T6 back the frames
    themselves are new, so a run from raw costs the same."""
    edited = [(e.scope, e.pattern) for n in rerun_stages(first) for e in STAGES[n].edits]
    out, seen = [], set()
    for n in rerun_stages(first):
        for r in STAGES[n].reads:
            if r.scope not in (REDUCTION, TARGET) or r.made_by is None or r.made_by >= first:
                continue
            key = (r.scope, r.pattern)
            if key in seen:
                continue
            seen.add(key)
            out.append((r.scope, r.pattern, r.made_by, key not in edited))
    return out


def rerun_outputs(first, scope=TARGET):
    """Patterns of what a rerun from `first` writes in a scope: the files it replaces."""
    return [w.pattern for n in rerun_stages(first) for w in STAGES[n].writes if w.scope == scope]


def is_output(rel, first, **kw):
    return any(matches(rel, p, **kw) for p in rerun_outputs(first))


def edited_by_rerun(rel, first, scope=TARGET, **kw):
    """Whether a stage of a rerun from `first` edits this file in place (an input it must have its own copy of)."""
    return any(e.scope == scope and matches(rel, e.pattern, **kw) for n in rerun_stages(first) for e in STAGES[n].edits)


def cause_stage(cause):
    key = cause.strip().lower()
    if key not in CAUSES:
        raise KeyError('unknown cause {!r}: one of {}'.format(cause, ', '.join(CAUSES)))
    return CAUSES[key][0]


def rerun_plan(first):
    """What a rerun from `first` runs, reads (made by earlier stages) and rewrites."""
    if first not in STAGES or first > LAST:
        raise ValueError('a rerun starts at T1..T{}'.format(LAST))
    return {
        'from': first,
        'stages': rerun_stages(first),
        'flags': ['--from_T{}'.format(first), '--no_T12'],
        'restage': first in RESTAGE,
        'inputs': [{'scope': s, 'pattern': p, 'made_by': m, 'link': ln} for s, p, m, ln in restage_inputs(first)],
        'outputs': [{'scope': w.scope, 'pattern': w.pattern, 'stage': n}
                    for n in rerun_stages(first) for w in STAGES[n].writes],
    }


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _items(items):
    return '; '.join('{} {}{}'.format(i.scope, i.pattern, ' ({})'.format(i.what) if i.what else '') for i in items)


def markdown():
    """The tables of docs/rerun-from-stage.md: stages, causes, and what each restageable rerun reads."""
    lines = ['| Stage | ZLP function | Script | Reads | Writes | Edits in place |', '|---|---|---|---|---|---|']
    for s in STAGES.values():
        reads = '<br>'.join('`{}` {}{}'.format(r.scope, r.pattern, ' (T{})'.format(r.made_by) if r.made_by else '')
                            for r in s.reads) or '-'
        writes = '<br>'.join('`{}` {}'.format(w.scope, w.pattern) for w in s.writes) or '-'
        edits = '<br>'.join('`{}` {}'.format(e.scope, e.pattern) for e in s.edits) or '-'
        lines.append('| T{} {} | `{}` | `{}` | {} | {} | {} |'.format(s.n, s.name, s.function, s.script, reads, writes,
                                                                     edits))
    lines += ['', '| Cause | Rerun from | Stages | What |', '|---|---|---|---|']
    for cause, (n, what) in CAUSES.items():
        lines.append('| `{}` | T{} | {} | {} |'.format(cause, n, ' '.join('T{}'.format(k) for k in rerun_stages(n)),
                                                       what))
    for n in RESTAGE:
        lines += ['', '| Rerun from T{}: input | Made by | Restaged as |'.format(n), '|---|---|---|']
        for scope, pattern, made_by, link in restage_inputs(n):
            how = 'copy, each path pointed at the restaged frame' if pattern in FRAME_LISTS else \
                'hard link' if link else 'copy'
            lines.append('| `{}` {} | T{} | {} |'.format(scope, pattern, made_by, how))
    return '\n'.join(lines) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser(description='The stage map of ZLP_pipeline.sh: what each stage reads and writes.')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--cause', help='what was fixed: ' + ', '.join(CAUSES))
    g.add_argument('--from', dest='first', type=lambda x: int(x.upper().lstrip('T')), help='first stage, e.g. T10')
    g.add_argument('--markdown', action='store_true', help='the tables of docs/rerun-from-stage.md')
    a = ap.parse_args(argv)
    if a.markdown:
        sys.stdout.write(markdown())
        return 0
    if a.cause or a.first:
        first = cause_stage(a.cause) if a.cause else a.first
        p = rerun_plan(first)
        if a.cause:
            print('{}: {}'.format(a.cause, CAUSES[a.cause.strip().lower()][1]))
        print('rerun from T{}: stages {}; ZLP_pipeline.sh {}'.format(
            first, ' '.join('T{}'.format(n) for n in p['stages']), ' '.join(p['flags'])))
        print('a night already in v2: {}'.format(
            'restage it first (utils/restage.py stage --from {}), then promote as an overlay'.format(first)
            if p['restage'] else 'rerun the whole night: a restage from T{} is not supported'.format(first)))
        print('reads, made by earlier stages:')
        for i in p['inputs']:
            print('  {:9} {:28} T{:<3} {}'.format(i['scope'], i['pattern'], i['made_by'],
                                                  'hard link' if i['link'] else 'copy'))
        print('rewrites:')
        for o in p['outputs']:
            print('  {:9} {:40} T{}'.format(o['scope'], o['pattern'], o['stage']))
        return 0
    for s in STAGES.values():
        print('T{:<3} {:14} {}  ({})'.format(s.n, s.name, s.function, s.script))
        print('     reads:  {}'.format(_items(s.reads) or '-'))
        print('     writes: {}'.format(_items(s.writes) or '-'))
        if s.edits:
            print('     edits:  {}'.format(_items(s.edits)))
        if s.note:
            print('     note:   {}'.format(s.note))
    return 0


if __name__ == '__main__':
    sys.exit(main())
