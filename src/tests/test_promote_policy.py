"""
Tests for the approval policy (utils/promote_policy.py): one case per kind of verdict found among the 122 pending on
2026-10-19 (mh_scratch/promote_v3_20261019/triage.jsonl), with the metrics of a real one, and the kinds of stale file.
"""
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import promote_policy as pp  # noqa: E402

SKY = {'solved': 50, 'data': 61, 'processing': 0, 'causes': {'data: Not enough stars detected for plate solve': 61}}
PIPELINE = {'solved': 80, 'data': 0, 'processing': 162,
            'causes': {'processing: Plate scale sanity check failed': 162}}


def m(**kw):
    """Metrics of a sound candidate with nothing in v2, changed by kw."""
    base = {'c_valid': 120, 'c_frames': 120, 'raw_frames': 121, 'c_frame_frac': 0.99, 'c_p2p': 0.004,
            'c_noise_ratio': 1.0, 'c_bw': 0, 'i_valid': 0, 'points_ratio': None, 'scatter_ratio': None}
    base.update(kw)
    return base


V2 = dict(i_valid=500)


CASES = [
    # (name, verdict, reasons, metrics, incumbent, stale, files, losses, outcome, rule)
    ('S1 stale plots only (v35)', 'HOLD', 'S1 1 of its 836 files predate the run', m(), False,
     ['lightcurves/x_5.png'], 836, None, pp.PROMOTE, 'clutter-only'),
    ('S1 old-name output and sidecars, v2 usable (v113)', 'HOLD', 'S1 44 of its 1599 files predate the run',
     m(**V2, points_ratio=1.0, scatter_ratio=1.05), True,
     ['SP1804-1616_I+z_20250724_output.fits', '1/procA.fits.ell', 'SP1804-1616_outstack_I+z.fits'], 1599, None,
     pp.PROMOTE, 'clutter-only'),
    ('S1 a stale aperture-8 curve is clutter (v191)', 'HOLD', 'S1 30 of its 6382 files predate the run', m(), False,
     ['26_zYJ_20250914_8_diff.fits', 'lightcurves/26_zYJ_20250914_8_MCMC'], 6382, None, pp.PROMOTE, 'clutter-only'),
    ('S1 a stale aperture-5 curve beside the new one', 'HOLD', 'S1 2 of its 900 files predate the run', m(), False,
     ['26_zYJ_20250914_5_diff.fits'], 900, None, pp.REPROCESS, 'stale-product'),
    ('S1 the run wrote none of the directory (v276)', 'HOLD', 'S1 the run wrote none of its 2051 files', m(), False,
     ['f%d' % k for k in range(2051)], 2051, None, pp.REPROCESS, 'run-did-not-write'),
    ('S2 no run', 'HOLD', 'S2 no pipeline run of the night covering X was found', m(), False, None, 10, None,
     pp.PERSON, 'no-run'),
    ('P1 a replacement without the PWV correction v2 has', 'HOLD', "P1 the candidate has no PWV correction (x has an "
     "empty PWV table or curve) and v2's curve has one", m(**V2, points_ratio=1.0, scatter_ratio=1.0), True, [], 50,
     None, pp.REPROCESS, 'pwv-lost'),
    ('A1 no aperture 5 (v189)', 'REJECT', 'A1 no aperture-5 light curve (apertures 3 4)', m(c_valid=0), False, [], 50,
     None, pp.REJECT, 'no-aperture-5'),
    ('A1 with v2 serving the target (v415)', 'REJECT', 'A1 no aperture-5 light curve (apertures 3); R1 v2 has an '
     'aperture-5 light curve of the target (998 points), the candidate none', m(c_valid=0, i_valid=998), True, [], 50,
     None, pp.REJECT, 'no-aperture-5'),
    ('A3 another star', 'HOLD', 'A3 light curve is 1, expected 2 (plan)', m(), False, [], 50, None, pp.PERSON,
     'identity'),
    ('no finite point (v68)', 'REJECT', 'A2 no finite target point in 1 frames; A4 light curve has 1 of 255 raw',
     m(c_valid=0, c_frames=1, raw_frames=255, c_frame_frac=0.0), True, [], 50, None, pp.REJECT, 'no-points'),
    ('degenerate (v132)', 'REJECT', 'A5 degenerate: p2p 4.9e-05 is 0.015 x the pipeline error',
     m(c_p2p=0.00005, c_noise_ratio=0.015), True, [], 50, None, pp.REJECT, 'degenerate'),
    ('noisier than v2 (v80)', 'HOLD', 'R2 p2p 0.0097 against 0.0068 in v2 (x1.42)',
     m(**V2, points_ratio=0.99, scatter_ratio=1.42), True, [], 50, None, pp.REJECT, 'worse-than-v2'),
    ('worse than v2 because the pipeline lost frames (v84)', 'HOLD', 'A4 light curve has 44 of 558 raw science '
     'frames (8%); R1 44 finite points against 557 in v2; R2 p2p x1.98', m(**dict(V2, c_valid=44, c_frames=44,
     raw_frames=558, c_frame_frac=0.08, points_ratio=0.08, scatter_ratio=1.98)), True, [], 50,
     dict(PIPELINE, processing=483), pp.REPROCESS, 'worse-than-v2-pipeline'),
    ('better than v2 though partial (v182)', 'HOLD', 'A4 light curve has 238 of 298 raw science frames (80%)',
     m(**dict(V2, c_valid=238, c_frames=238, raw_frames=298, c_frame_frac=0.799, points_ratio=1.35,
              scatter_ratio=0.30)), True, [], 50, None, pp.PROMOTE, 'better-than-v2'),
    ('within 5% of v2 (v123)', 'HOLD', 'R1 96 finite points against 99 in v2',
     m(**V2, points_ratio=0.97, scatter_ratio=0.96), True, [], 50, None, pp.PROMOTE, 'same-as-v2'),
    ('fewer points, more scatter, not by much (v131)', 'HOLD', 'R1 444 finite points against 483 in v2',
     m(**V2, points_ratio=0.92, scatter_ratio=1.16), True, [], 50, None, pp.PERSON, 'mixed-against-v2'),
    ('fewer points but cleaner (v333)', 'HOLD', 'A4 light curve has 29 of 44 raw science frames (66%); R1 29 finite',
     m(**V2, points_ratio=0.67, scatter_ratio=0.83), True, [], 50, None, pp.PERSON, 'mixed-against-v2'),
    ('partial: the pipeline lost the frames (v120)', 'HOLD', 'A4 light curve has 80 of 243 raw science frames (33%)',
     m(c_valid=80, c_frames=80, raw_frames=243, c_frame_frac=0.33), False, [], 50, PIPELINE, pp.REPROCESS,
     'frames-lost-to-pipeline'),
    ('partial: the sky lost the frames (v219)', 'HOLD', 'A4 light curve has 60 of 121 raw science frames (50%)',
     m(c_valid=60, c_frames=60, raw_frames=121, c_frame_frac=0.5), False, [], 50, SKY, pp.PROMOTE,
     'partial-night-data'),
    ('partial: no log to say why', 'HOLD', 'A4 light curve has 60 of 121 raw science frames (50%)',
     m(c_valid=60, c_frames=60, raw_frames=121, c_frame_frac=0.5), False, [], 50, None, pp.PERSON,
     'partial-night-unexplained'),
    ('noisy on a cloudy night (v72)', 'HOLD', 'A5 p2p 0.076 above 0.050',
     m(c_valid=588, c_frames=588, raw_frames=625, c_frame_frac=0.94, c_p2p=0.076, c_noise_ratio=3.58, c_bw=493),
     False, [], 50, {'solved': 588, 'data': 34, 'processing': 1, 'causes': {}}, pp.PROMOTE, 'noisy-night-data'),
    ('noisy on a clear night', 'HOLD', 'A5 p2p 0.076 above 0.050',
     m(c_p2p=0.076, c_noise_ratio=3.58, c_bw=0), False, [], 50, {'solved': 120, 'data': 0, 'processing': 1,
                                                                 'causes': {}}, pp.PERSON, 'noisy-clear-night'),
    ('noise, not a light curve, on a clear night', 'HOLD', 'A4 light curve has 25 of 278 raw (9%); A5 p2p 2.784',
     m(c_valid=152, c_frames=152, raw_frames=272, c_frame_frac=0.56, c_p2p=2.78, c_noise_ratio=367.0), False, [],
     50, {'solved': 152, 'data': 0, 'processing': 2, 'causes': {}}, pp.REJECT, 'noise-not-signal'),
    ('a short complete night (v420)', 'HOLD', 'A2 only 19 finite target points',
     m(c_valid=19, c_frames=19, raw_frames=20, c_frame_frac=0.95), False, [], 50, None, pp.PROMOTE, 'short-night'),
    ('a few points of a long night', 'HOLD', 'A2 only 14 finite target points; A4 light curve has 14 of 290 (5%)',
     m(c_valid=14, c_frames=14, raw_frames=290, c_frame_frac=0.05), False, [], 50,
     {'solved': 14, 'data': 3, 'processing': 2, 'causes': {}}, pp.REJECT, 'too-few-points'),
]


@pytest.mark.parametrize('case', CASES, ids=[c[0] for c in CASES])
def test_decide(case):
    name, verdict, reasons, metrics, incumbent, stale, files, losses, outcome, rule = case
    d = pp.decide(verdict, metrics, reasons, incumbent, stale, files, losses)
    assert (d.outcome, d.rule) == (outcome, rule), d.why
    assert d.by == 'policy:v{}:{}'.format(pp.VERSION, rule) and d.why


def test_stale_kinds():
    assert pp.stale_kind('123_I+z_20250101_5_diff.fits') == 'diff-ap5'
    assert pp.stale_kind('123_I+z_20250101_8_diff.fits') == 'diff'
    assert pp.stale_kind('lightcurves/123_I+z_20250101_5_MCMC') == 'mcmc-ap5'
    assert pp.stale_kind('lightcurves/123_I+z_20250101_5.png') == 'plot'
    assert pp.stale_kind('123_I+z_20250101_output.fits') == 'output'
    assert pp.stale_kind('SP1804-1616_I+z_20250724_output.fits') == 'output-old-name'
    assert pp.stale_kind('1/procSPECU2.2025-08-12T01:35:25.790.fits.ell') == 'run-dir'
    assert pp.stale_kind('6525046188759705728_Ganymede_andor_20230928_outstack_I+z.fits') == 'other'


def test_frame_losses_reads_the_targets_own_section_and_the_last_attempt(tmp_path):
    log = tmp_path / 'x.log'
    log.write_text('\n'.join([
        'TARGET: SpA',
        '[1/3] procA.fits FAILED - Sources: 250→16 Gaia: 51→32 Matches: 0 Error: Plate scale sanity check '
        'failed: solved=25.674 arcsec/px, expected=0.309 arcsec/px, deviation=8198.0%',
        '[2/3] procB.fits FAILED - Sources: 4→0 Gaia: 0→0 Matches: 0 Error: Not enough stars detected for '
        'plate solve: 4 [3/3] procC.fits SUCCESS - Sources: 90→16 Gaia: 80→32 Matches: 13 CRPIX=(1,2) '
        'CRVAL=(3,4)',
        '[1/3] procA.fits FAILED - Sources: 240→16 Gaia: 51→32 Matches: 0 Error: Plate solve failed, twirl '
        'found no 4-star asterism shared by x',
        'TARGET: SpB',
        '[1/1] procD.fits FAILED - Sources: 300→16 Gaia: 51→32 Matches: 0 Error: Plate scale sanity',
    ]) + '\n', encoding='utf-8')
    a = pp.frame_losses(str(log), 'SpA')
    assert (a['solved'], a['data'], a['processing']) == (1, 1, 1)
    assert pp.frame_losses(str(log), 'SpB')['processing'] == 1
    assert pp.frame_losses(str(tmp_path / 'missing.log'), 'SpA') is None
