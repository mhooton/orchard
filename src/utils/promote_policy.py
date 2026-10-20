"""
promote_policy.py
=================

The approval policy of the promotion step (utils/promote.py): what to do with a candidate the checks held or
rejected, so that a person sees only what the evidence cannot settle. promote.py asks it about every pending HOLD
or REJECT when it runs with --apply, and journals each answer with the rule that gave it
(decided_by = policy:v<VERSION>:<rule>).

The philosophy (Matthew Hooton, 2026-10-19): do not judge what should and should not be on the portal; process
everything we can and put it there, ideally with the tools that show users when the conditions, not the pipeline,
make data poor. So a candidate is withheld only when it is broken (no aperture-5 curve, no finite point, a
degenerate curve), when v2 already serves something better, when its directory mixes runs, or when its identity is
in doubt. Where the pipeline lost frames (plate-solve failures on frames full of stars) the curve is still published
when v2 has nothing better, and the decision carries `rerun`, so the night is rerun once the cause is fixed at
source.

Outcomes
    promote    publish the candidate as `accept` would; the files the run did not write are set aside first
    reject     close the verdict; v2 keeps what it has
    reprocess  not publishable as it is (a directory that mixes runs' served products, or a replacement that lost
               v2's PWV correction or is worse than v2 because the pipeline lost frames): pending, for a rerun
    person     the evidence is mixed (identity, no known run, better and worse than v2 at once): pending
Any outcome can carry `rerun`: the processing, not the data, is at fault.

The run's log says why each frame missing from a curve dropped out: a plate solve that found too few stars is the
sky; one rejected on a frame with plenty of stars is the pipeline. Bad-weather flags on the curve's own frames say the
same for noise. Measured on the 122 verdicts pending on 2026-10-19 (mh_scratch/promote_v3_20261019/triage.jsonl): S1
was almost always clutter, and the partial curves had lost their frames to "Plate scale sanity check failed" on
frames with 130-3,100 sources.
"""
import os
import re
from collections import Counter

VERSION = 1

THRESHOLDS = {
    'stale_rerun_fraction': 0.9,     # S1: this share of a directory older than its run means the run did not make it
    'equal_points': 0.95,            # against v2: within 5% on points ...
    'equal_scatter': 1.05,           # ... and on scatter counts as the same
    'worse_points': 0.9,             # against v2: fewer points than this, and no less scatter, is worse
    'worse_scatter': 1.2,            # against v2: more scatter than this, and no more points, is worse
    'processing_loss_fraction': 0.1,  # raw frames lost to the pipeline above this share: rerun the night
    'processing_loss_min': 5,        # ... and at least this many
    'stars_for_processing': 20,      # a solve that failed on a frame with this many sources was the pipeline's
    'weather_fraction': 0.5,         # this share of the curve's frames flagged bad weather: the night, not the pipeline
    'short_night_fraction': 0.9,     # a night with nearly every frame in the curve, however few, is complete
}

PROMOTE, REJECT, REPROCESS, PERSON = 'promote', 'reject', 'reprocess', 'person'

CODES = re.compile(r'\b(S1|S2|P1|O1|O2|O3|A1|A2|A3|A4|A5|R1|R2|R3)\b')
DIFF = re.compile(r'_(\d+)_diff\.fits$')
MCMC = re.compile(r'_(\d+)_MCMC$')
# astrometry.py's per-frame summary: '[i/n] <frame> SUCCESS - Sources: a→b Gaia: c→d Matches: m CRPIX=... CRVAL=...'
# or '... FAILED - Sources: a→b Gaia: c→d Matches: m Error: <why>'. The rest runs to the next summary on the line, if
# parallel workers printed two together, or to its end.
SUMMARY = re.compile(r'(?P<name>proc\S+?\.fits(?:\.fz)?) (?P<state>SUCCESS|FAILED) - Sources: (?P<src>\d+|\?)\D+?'
                     r'(?:\d+|\?) Gaia: (?:\d+|\?)\D+?(?:\d+|\?) Matches: (?:\d+|\?)(?P<rest>.*?)'
                     r'(?=\s*\[\d+/\d+\]\s+proc|\s+proc\S+?\.fits(?:\.fz)? (?:SUCCESS|FAILED) - |$)')
TARGET_LINE = re.compile(r'^TARGET: (\S+)')
# 'Sources trailed': T7 (astrom.pointer_wcs) found the brightest sources stretched all one way, as when the telescope
# stopped tracking (Ganymede 20250716 and 20250825), and no solution verified. A rerun cannot recover those frames.
DATA_ERRORS = ('Not enough stars detected', 'too few stars', 'Sources trailed')


class Decision(object):
    def __init__(self, outcome, rule, why, rerun=False):
        self.outcome, self.rule, self.why, self.rerun = outcome, rule, why, bool(rerun)

    @property
    def by(self):
        return 'policy:v{}:{}'.format(VERSION, self.rule)

    def as_dict(self):
        return {'outcome': self.outcome, 'rule': self.rule, 'why': self.why, 'version': VERSION, 'rerun': self.rerun}

    def __repr__(self):
        return 'Decision({!r}, {!r}, {!r}, rerun={!r})'.format(self.outcome, self.rule, self.why, self.rerun)


def codes(reasons):
    """The check codes a verdict's reasons name (A1-A5, R1-R3, S1, S2)."""
    return set(CODES.findall(reasons or ''))


# ---------------------------------------------------------------------------
# S1: what predates the run
# ---------------------------------------------------------------------------

def stale_kind(rel):
    """What a file the run did not write is. Served products: 'diff-ap5' (what the portal reads), 'mcmc-ap5',
    'output' (a Gaia-named *_output.fits). Clutter: 'diff' and 'mcmc' of other apertures, 'output-old-name' (the
    SP-named *_output.fits of the old naming), 'run-dir' (frames, per-frame files and frame lists), 'plot', 'other'."""
    parts = rel.split(os.sep)
    name = parts[-1]
    if len(parts) > 1 and parts[0].isdigit():
        return 'run-dir'
    m = DIFF.search(name)
    if m:
        return 'diff-ap5' if m.group(1) == '5' else 'diff'
    m = MCMC.search(name)
    if m:
        return 'mcmc-ap5' if m.group(1) == '5' else 'mcmc'
    if name.endswith('_output.fits'):
        return 'output' if name[:1].isdigit() else 'output-old-name'
    if name.lower().endswith(('.png', '.pdf', '.jpg', '.jpeg', '.gif')):
        return 'plot'
    return 'other'


SERVED = ('diff-ap5', 'mcmc-ap5', 'output')


def stale_summary(stale):
    return dict(Counter(stale_kind(rel) for rel in stale or ()))


# ---------------------------------------------------------------------------
# Why frames dropped out: the run's log
# ---------------------------------------------------------------------------

def frame_losses(log_path, target, th=THRESHOLDS):
    """Plate-solve outcomes of the target's frames in a ZLP_pipeline.sh log: the frames whose last attempt failed,
    as 'data' (too few stars on the frame: the sky) or 'processing' (a solve that failed or was rejected on a frame
    with plenty of sources: the pipeline). None if the log cannot be read."""
    if not log_path or not os.path.isfile(log_path):
        return None
    want = (target or '').upper()
    last, current = {}, None
    try:
        with open(log_path, 'rb') as f:
            for raw in f:
                line = raw.decode('utf-8', 'replace').rstrip('\r\n')
                m = TARGET_LINE.match(line)
                if m:
                    current = m.group(1).upper()
                    continue
                if current != want or ' - Sources: ' not in line:
                    continue
                for s in SUMMARY.finditer(line):        # parallel workers sometimes share a line
                    rest = s.group('rest')
                    error = rest.split(' Error: ', 1)[1].strip() if ' Error: ' in rest else ''
                    last[s.group('name')] = (s.group('state'), s.group('src'), error)
    except OSError:
        return None
    out = {'solved': 0, 'data': 0, 'processing': 0, 'causes': Counter()}
    for state, src, error in last.values():
        if state == 'SUCCESS':
            out['solved'] += 1
            continue
        few = src.isdigit() and int(src) < th['stars_for_processing']
        kind = 'data' if few or any(e in error for e in DATA_ERRORS) else 'processing'
        out[kind] += 1
        out['causes'][(kind, error.split(':')[0][:60] or 'unknown')] += 1
    out['causes'] = {'{}: {}'.format(*k): n for k, n in out['causes'].most_common()}
    return out


# ---------------------------------------------------------------------------
# The decision
# ---------------------------------------------------------------------------

def _num(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if x == x else None        # NaN -> None


def decide(verdict, metrics, reasons, incumbent, stale, files, losses, th=THRESHOLDS):
    """The policy's answer for one HOLD or REJECT.

    verdict    'HOLD' or 'REJECT' (a PROMOTE needs no policy)
    metrics    promote_check.assess() row for the candidate (c_*, i_*, raw_frames, points_ratio, scatter_ratio)
    reasons    the verdict's reasons ('; '-joined, each starting with its check code: A1-A5 and R1-R3 from the
               checker, S1, S2, P1 and O1-O3 (restaged reruns) from promote.py)
    incumbent  whether v2 has a directory for the target-night
    stale      relative paths of the candidate's files that predate its science run (None if no run is known)
    files      how many files the candidate holds
    losses     frame_losses() of the run's log for the target, or None
    """
    th = dict(THRESHOLDS, **(th or {}))
    cs = codes(reasons)
    checks = cs - {'S1'}
    m = metrics or {}
    valid = _num(m.get('c_valid')) or 0
    frames = _num(m.get('c_frames')) or 0
    raw = _num(m.get('raw_frames'))
    frac = _num(m.get('c_frame_frac'))
    p2p, noise = _num(m.get('c_p2p')), _num(m.get('c_noise_ratio'))
    pr, sr = _num(m.get('points_ratio')), _num(m.get('scatter_ratio'))
    bw = _num(m.get('c_bw')) or 0
    usable_v2 = bool(incumbent) and (_num(m.get('i_valid')) or 0) >= 20

    lost = (losses or {}).get('processing', 0)
    by_pipeline = bool(raw) and lost >= max(th['processing_loss_min'], th['processing_loss_fraction'] * raw)
    lost_text = '{} of {:.0f} raw frames lost to plate-solve failures on frames with stars ({})'.format(
        lost, raw or 0, ', '.join('{} {}'.format(n, c) for c, n in sorted((losses or {}).get('causes', {}).items(),
                                                                          key=lambda kv: -kv[1])[:3]))
    rerun_text = '; rerun once the cause is fixed: ' + lost_text if by_pipeline else ''

    if 'S2' in cs:
        return Decision(PERSON, 'no-run', 'no pipeline run of the night covering the target was found, so its '
                        'files cannot be checked against the run')
    # a rerun from a late stage on a night restaged from v2 (utils/restage.py): its inputs must still be what was
    # restaged, and v2 what it was then
    if 'O1' in cs or 'O3' in cs:
        return Decision(PERSON, 'restage-inputs', "the rerun's inputs are not what was restaged (or a file of it is "
                        "another release's): look before anything moves; restage.py show lists them")
    if 'O2' in cs:
        return Decision(REPROCESS, 'v2-changed-since-restage', 'v2 changed after the night was restaged: clear the '
                        'restage, restage again and rerun', rerun=True)
    if 'P1' in cs:
        return Decision(REPROCESS, 'pwv-lost', "the run made no PWV correction and v2's curve has one: v2 stays "
                        'until a rerun gets PWV', rerun=True)
    if 'S1' in cs:
        if stale is None:
            return Decision(PERSON, 'stale-unknown', 'S1 held it, but the run is not known to recheck its files')
        kinds = stale_summary(stale)
        if files and len(stale) >= th['stale_rerun_fraction'] * files:
            return Decision(REPROCESS, 'run-did-not-write', '{} of the {} files predate the latest science run: it '
                            'did not make this directory'.format(len(stale), files), rerun=True)
        served = sorted(k for k in kinds if k in SERVED)
        if served:
            return Decision(REPROCESS, 'stale-product', 'files the portal serves ({}) predate the latest science '
                            'run, beside its new ones'.format(', '.join(served)), rerun=True)
        if not checks:
            # the checks themselves passed, against v2 too: S1 was the only objection, and it is clutter
            return Decision(PROMOTE, 'clutter-only', 'only files the portal does not serve predate the run ({}); '
                            'they are set aside'.format(', '.join('{} {}'.format(n, k) for k, n in
                                                                  sorted(kinds.items()))) + rerun_text,
                            rerun=by_pipeline)

    # identity: never settled by numbers
    if 'R3' in checks:
        return Decision(REJECT, 'v2-right-star', "v2 is the plan's star and the candidate is not, or v2 has the "
                        "target's aperture 5 and the candidate none")
    if 'A3' in checks and 'A1' not in checks:
        return Decision(PERSON, 'identity', 'the curve is not of the star the plan or the target table names')

    # broken: nothing to publish
    if 'A1' in checks:
        return Decision(REJECT, 'no-aperture-5', 'no aperture-5 light curve' +
                        ('; v2 keeps its own' if usable_v2 else ''), rerun=by_pipeline)
    if valid == 0:
        return Decision(REJECT, 'no-points', 'no finite target point' + rerun_text, rerun=by_pipeline)
    if noise is not None and noise < th.get('min_noise_ratio', 0.05):
        return Decision(REJECT, 'degenerate', 'degenerate curve: p2p {:.2g} of the pipeline error'.format(noise),
                        rerun=by_pipeline)

    # against what v2 serves: never replace it with something worse
    if usable_v2 and pr is not None:
        ratios = 'x{:.2f} the points of v2, x{} its scatter'.format(pr, '{:.2f}'.format(sr) if sr is not None else '-')
        if pr >= 1.0 and (sr is None or sr <= 1.0):
            return Decision(PROMOTE, 'better-than-v2', ratios + rerun_text, rerun=by_pipeline)
        if pr >= th['equal_points'] and sr is not None and sr <= th['equal_scatter']:
            return Decision(PROMOTE, 'same-as-v2', 'within {:.0%} of v2: '.format(1 - th['equal_points']) + ratios +
                            rerun_text, rerun=by_pipeline)
        worse = (pr < th['worse_points'] and (sr is None or sr >= 1.0)) or \
            (sr is not None and sr > th['worse_scatter'] and pr <= 1.0)
        if worse and by_pipeline:
            return Decision(REPROCESS, 'worse-than-v2-pipeline', 'v2 is better and stays (' + ratios + ')' +
                            rerun_text, rerun=True)
        if worse:
            return Decision(REJECT, 'worse-than-v2', 'v2 is better and stays: ' + ratios)
        return Decision(PERSON, 'mixed-against-v2', 'mixed against v2: ' + ratios + rerun_text, rerun=by_pipeline)

    # nothing usable in v2: publish every light curve there is, saying what made it poor
    if by_pipeline:
        return Decision(PROMOTE, 'partial-pipeline', 'published as it is (nothing usable in v2); ' + lost_text +
                        ': rerun once plate solving is fixed', rerun=True)
    weather = bool(frames) and bw >= th['weather_fraction'] * frames
    weather_text = '; {:.0%} of its frames flagged bad weather'.format(bw / frames) if frames and bw else ''
    if 'A5' in checks:
        return Decision(PROMOTE, 'noisy-night-data' if weather else 'noisy',
                        'noisy (p2p {}, x{} the pipeline error){}'.format(
                            '{:.3f}'.format(p2p) if p2p is not None else '?',
                            '{:.1f}'.format(noise) if noise is not None else '?', weather_text))
    if 'A2' in checks:
        complete = frac is not None and frac >= th['short_night_fraction'] and \
            valid >= th['short_night_fraction'] * frames
        return Decision(PROMOTE, 'short-night' if complete else 'few-points', '{:.0f} points on {} of the '
                        'frames'.format(valid, '{:.0%}'.format(frac) if frac is not None else '?') + weather_text)
    if 'A4' in checks:
        missing = max(0.0, (raw or 0) - frames)
        by_data = (losses or {}).get('data', 0)
        sky = losses is not None and missing > 0 and by_data >= 0.5 * missing
        return Decision(PROMOTE, 'partial-night-data' if sky else 'partial',
                        'a partial night ({:.0%} of the frames){}'.format(frac or 0, ': {} of the {:.0f} missing '
                        'frames had too few stars to solve (the sky)'.format(by_data, missing) if sky else
                        ', not explained by the log' if losses is None or missing == 0 else
                        ': {} of the {:.0f} missing frames had too few stars to solve'.format(by_data, missing)) +
                        weather_text)
    return Decision(PROMOTE, 'published', 'checks {} with nothing usable in v2'.format(' '.join(sorted(checks))))
