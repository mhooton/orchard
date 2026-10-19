"""
promote_policy.py
=================

The approval policy of the promotion step (utils/promote.py): what to do with a candidate the checks held or
rejected, so that a person sees only what the evidence cannot settle. promote.py asks it about every pending HOLD
or REJECT when it runs with --apply, and journals each answer with the rule that gave it
(decided_by = policy:v<VERSION>:<rule>).

Outcomes
    promote    publish the candidate as `accept` would; the files the run did not write are set aside first
    reject     close the verdict; v2 keeps what it has
    reprocess  the processing, not the data, is at fault: frames lost to plate-solve failures on frames full of
               stars, or a directory the latest science run did not write. The verdict stays pending and is
               listed for a rerun, not for a person
    person     the evidence is mixed; the verdict stays pending for a person

The deciding question for a partial or noisy curve is whether the data or the processing made it so: noisy is fine
when the target is right and the data are genuinely noisy, not when the pipeline processed them poorly (Matthew
Hooton, 2026-10-19). The run's log says why each missing frame dropped out: a plate solve that found too few stars is
the sky; a solution rejected on a frame with plenty of stars is the pipeline. Bad-weather flags on the curve's own
frames say the same for noise.

Measured on the 122 verdicts pending on 2026-10-19 (mh_scratch/promote_v3_20261019/triage.jsonl): S1 was
almost always clutter (plots, per-frame sidecars, backup stacks, SP-named _output.fits left from the old naming);
the partial curves had lost their frames to "Plate scale sanity check failed" on frames with 130-3,100 sources.
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
    'processing_loss_fraction': 0.1,  # raw frames lost to the pipeline above this share: reprocess
    'processing_loss_min': 5,        # ... and at least this many
    'stars_for_processing': 20,      # a solve that failed on a frame with this many sources was the pipeline's
    'weather_fraction': 0.5,         # this share of the curve's frames flagged bad weather: the night, not the pipeline
    'short_night_fraction': 0.9,     # a short night with nearly every frame in the curve ...
    'short_night_points': 10,        # ... and at least this many points is published
    'sound_noise_low': 0.3,          # a sound curve: p2p / pipeline error between these
    'sound_noise_high': 3.0,
}

PROMOTE, REJECT, REPROCESS, PERSON = 'promote', 'reject', 'reprocess', 'person'

CODES = re.compile(r'\b(S1|S2|P1|A1|A2|A3|A4|A5|R1|R2|R3)\b')
DIFF = re.compile(r'_(\d+)_diff\.fits$')
MCMC = re.compile(r'_(\d+)_MCMC$')
# astrometry.py's per-frame summary: '[i/n] <frame> SUCCESS - Sources: a→b Gaia: c→d Matches: m CRPIX=... CRVAL=...'
# or '... FAILED - Sources: a→b Gaia: c→d Matches: m Error: <why>'. The rest runs to the next summary on the line, if
# parallel workers printed two together, or to its end.
SUMMARY = re.compile(r'(?P<name>proc\S+?\.fits(?:\.fz)?) (?P<state>SUCCESS|FAILED) - Sources: (?P<src>\d+|\?)\D+?'
                     r'(?:\d+|\?) Gaia: (?:\d+|\?)\D+?(?:\d+|\?) Matches: (?:\d+|\?)(?P<rest>.*?)'
                     r'(?=\s*\[\d+/\d+\]\s+proc|\s+proc\S+?\.fits(?:\.fz)? (?:SUCCESS|FAILED) - |$)')
TARGET_LINE = re.compile(r'^TARGET: (\S+)')
DATA_ERRORS = ('Not enough stars detected', 'too few stars')


class Decision(object):
    def __init__(self, outcome, rule, why):
        self.outcome, self.rule, self.why = outcome, rule, why

    @property
    def by(self):
        return 'policy:v{}:{}'.format(VERSION, self.rule)

    def as_dict(self):
        return {'outcome': self.outcome, 'rule': self.rule, 'why': self.why, 'version': VERSION}

    def __repr__(self):
        return 'Decision({!r}, {!r}, {!r})'.format(self.outcome, self.rule, self.why)


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
               checker, S1, S2 and P1 from promote.py)
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

    if 'S2' in cs:
        return Decision(PERSON, 'no-run', 'no pipeline run of the night covering the target was found, so its '
                        'files cannot be checked against the run')
    if 'P1' in cs:
        return Decision(REPROCESS, 'pwv-lost', 'the run made no PWV correction and v2\'s curve has one: v2 stays '
                        'until a rerun gets PWV')
    if 'S1' in cs:
        if stale is None:
            return Decision(PERSON, 'stale-unknown', 'S1 held it, but the run is not known to recheck its files')
        kinds = stale_summary(stale)
        if files and len(stale) >= th['stale_rerun_fraction'] * files:
            return Decision(REPROCESS, 'run-did-not-write', '{} of the {} files predate the latest science run: it '
                            'did not make this directory'.format(len(stale), files))
        served = sorted(k for k in kinds if k in SERVED)
        if served:
            return Decision(REPROCESS, 'stale-product', 'files the portal serves ({}) predate the latest science '
                            'run, beside its new ones'.format(', '.join(served)))
        if not checks:
            # the checks themselves passed, against v2 too: S1 was the only objection, and it is clutter
            return Decision(PROMOTE, 'clutter-only', 'only files the portal does not serve predate the run ({}); '
                            'they are set aside'.format(', '.join('{} {}'.format(n, k) for k, n in
                                                                  sorted(kinds.items()))))

    # identity: never settled by numbers
    if 'R3' in checks:
        return Decision(REJECT, 'v2-right-star', 'v2 is the plan\'s star and the candidate is not, or v2 has the '
                        'target\'s aperture 5 and the candidate none')
    if 'A3' in checks and 'A1' not in checks:
        return Decision(PERSON, 'identity', 'the curve is not of the star the plan or the target table names')

    # nothing to publish
    if 'A1' in checks:
        return Decision(REJECT, 'no-aperture-5', 'no aperture-5 light curve' +
                        ('; v2 keeps its own' if usable_v2 else ''))
    if valid == 0:
        return Decision(REJECT, 'no-points', 'no finite target point')
    if noise is not None and noise < th.get('min_noise_ratio', 0.05):
        return Decision(REJECT, 'degenerate', 'degenerate curve: p2p {:.2g} of the pipeline error'.format(noise))

    lost = (losses or {}).get('processing', 0)
    by_pipeline = bool(raw) and lost >= max(th['processing_loss_min'], th['processing_loss_fraction'] * raw)
    lost_text = '{} of {:.0f} raw frames lost to plate-solve failures on frames with stars ({})'.format(
        lost, raw or 0, ', '.join('{} {}'.format(n, c) for c, n in sorted((losses or {}).get('causes', {}).items(),
                                                                          key=lambda kv: -kv[1])[:3]))

    # against what v2 serves
    if usable_v2 and pr is not None:
        if pr >= 1.0 and (sr is None or sr <= 1.0):
            return Decision(PROMOTE, 'better-than-v2', 'x{:.2f} the points of v2 and x{} its scatter'.format(
                pr, '{:.2f}'.format(sr) if sr is not None else '-'))
        if pr >= th['equal_points'] and sr is not None and sr <= th['equal_scatter']:
            return Decision(PROMOTE, 'same-as-v2', 'within {:.0%} of v2: x{:.2f} the points, x{:.2f} the scatter'.format(
                1 - th['equal_points'], pr, sr))
        worse = (pr < th['worse_points'] and (sr is None or sr >= 1.0)) or \
            (sr is not None and sr > th['worse_scatter'] and pr <= 1.0)
        if worse and by_pipeline:
            return Decision(REPROCESS, 'worse-than-v2-pipeline', 'v2 is better (x{:.2f} points, x{} scatter) and '
                            'stays; {}'.format(pr, '{:.2f}'.format(sr) if sr is not None else '-', lost_text))
        if worse:
            return Decision(REJECT, 'worse-than-v2', 'v2 is better and stays: x{:.2f} the points, x{} the '
                            'scatter'.format(pr, '{:.2f}'.format(sr) if sr is not None else '-'))
        return Decision(PERSON, 'mixed-against-v2', 'mixed against v2: x{:.2f} the points, x{} the scatter'.format(
            pr, '{:.2f}'.format(sr) if sr is not None else '-'))

    # nothing usable in v2: publish what the data allow, reprocess what the pipeline lost
    if by_pipeline:
        return Decision(REPROCESS, 'frames-lost-to-pipeline', lost_text)
    sound = p2p is not None and p2p <= th.get('max_p2p', 0.05) and noise is not None and \
        th['sound_noise_low'] <= noise <= th['sound_noise_high']
    weather = bool(frames) and bw >= th['weather_fraction'] * frames
    if checks <= {'A2'} and frac is not None and frac >= th['short_night_fraction'] and \
            valid >= th['short_night_points'] and valid >= th['short_night_fraction'] * frames:
        return Decision(PROMOTE, 'short-night', 'a short night: {:.0f} points on {:.0%} of its frames'.format(
            valid, frac))
    if valid < th.get('min_valid_points', 20) and (frac is None or frac < 0.5):
        return Decision(REJECT, 'too-few-points', '{:.0f} points on {} of the frames'.format(
            valid, '{:.0%}'.format(frac) if frac is not None else '?'))
    if checks <= {'A4'} and sound:
        missing = max(0.0, (raw or 0) - frames)
        by_data = (losses or {}).get('data', 0)
        if losses is not None and missing > 0 and by_data >= 0.5 * missing:
            return Decision(PROMOTE, 'partial-night-data', 'a partial night ({:.0%} of the frames): {} of the {:.0f} '
                            'missing frames had too few stars to solve (the sky, not the pipeline)'.format(
                                frac or 0, by_data, missing))
        return Decision(PERSON, 'partial-night-unexplained', 'a partial night ({:.0%} of the frames) the log does '
                        'not explain'.format(frac or 0))
    if checks <= {'A4', 'A5'} and 'A5' in checks:
        if weather and (frac is None or frac >= 0.5):
            return Decision(PROMOTE, 'noisy-night-data', 'noisy (p2p {}) on a night with {:.0%} of its frames '
                            'flagged bad weather'.format('{:.3f}'.format(p2p) if p2p is not None else '?',
                                                         bw / frames))
        if p2p is not None and p2p > 0.25:
            return Decision(REJECT, 'noise-not-signal', 'p2p {:.2f}: noise, not a light curve'.format(p2p))
        return Decision(PERSON, 'noisy-clear-night', 'noisy (p2p {}, x{} the pipeline error) without a bad-weather '
                        'night to explain it'.format('{:.3f}'.format(p2p) if p2p is not None else '?',
                                                     '{:.1f}'.format(noise) if noise is not None else '?'))
    return Decision(PERSON, 'unsettled', 'checks {} with nothing usable in v2'.format(' '.join(sorted(checks))))
