"""
Unit tests for T7's second pass (astrom/astrometry.py): which failed frames are tried again, and from which frame's
solution.

A frame that failed with 5 or more sources, not trailed and not timed out, is tried again from the WCS of the frame
nearest to it in time that has one: solved this run, or applied from the WCS cache (every listed frame not in
to_solve). Callisto 20261012 (SPIRIT) is the case: faint frames of a sparse field with 4-5 sources on catalogue
stars, among frames that solved.

Needs the pipeline's own environment (astrometry.py imports twirl): run inside orchard-server.
Run:  pytest src/tests/test_neighbour_pass.py -q
"""
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pytest.importorskip('twirl')
from astrom import astrometry as am      # noqa: E402

WCS = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRVAL1': 150.0, 'CRVAL2': -30.0, 'CRPIX1': 512.0,
       'CRPIX2': 640.0, 'CD1_1': -8.6e-5, 'CD1_2': 0.0, 'CD2_1': 0.0, 'CD2_2': 8.6e-5}


def header(minute, solved, second=0, tag=None):
    h = {'DATE-OBS': '2026-10-13T03:%02d:%02d.500' % (minute, second)}
    if solved:
        h.update(WCS)
        h['CRVAL1'] = 150.0 + (tag or 0) * 1e-3     # tells the neighbours apart
    return h


def failed(sources, error='No verified solution (best 4 of 5 sources on catalogue stars); twirl: ...'):
    return dict(success=False, sources_detected=sources, error=error)


OK = dict(success=True, sources_detected=200)


def test_worth_a_second_pass():
    assert am.worth_a_second_pass(failed(5))
    assert not am.worth_a_second_pass(failed(4))
    assert not am.worth_a_second_pass(failed(None))
    assert not am.worth_a_second_pass(OK)
    assert not am.worth_a_second_pass(None)
    assert not am.worth_a_second_pass(failed(300, 'Sources trailed (axis ratio 3.1, all one way), the telescope ...'))
    assert not am.worth_a_second_pass(failed(300, 'Plate solving timed out after 60s during star_matching'))


def test_obs_seconds():
    a = am._obs_seconds({'DATE-OBS': '2026-10-13T03:10:00.500'})
    b = am._obs_seconds({'DATE-OBS': '2026-10-13T03:11:00'})
    assert b - a == pytest.approx(59.5)
    assert am._obs_seconds({}) is None
    assert am._obs_seconds({'DATE-OBS': 'nonsense'}) is None


def test_each_failure_gets_the_nearest_frame_with_a_solution():
    frames = ['f%d' % i for i in range(6)]
    headers = {'f0': header(0, True, tag=0), 'f1': header(1, False), 'f2': header(2, True, tag=2),
               'f3': header(9, False), 'f4': header(10, True, tag=4), 'f5': header(11, False)}
    results = [OK, failed(6), OK, failed(5), OK, failed(3)]
    jobs = am.neighbour_jobs(frames, results, frames, read_header=headers.get)
    picked = {i: cards['CRVAL1'] for i, cards in jobs}
    assert sorted(picked) == [1, 3]                       # f5 has too few sources
    assert picked[1] == pytest.approx(150.0)              # f0 and f2 tie at one minute: the first listed wins
    assert picked[3] == pytest.approx(150.004)            # f4, a minute after
    assert set(dict(jobs)[1]) <= set(am.WCS_KEYWORDS)


def test_frames_from_the_cache_are_neighbours_too():
    # on a rerun the frames that solved before come from the WCS cache and only the failures are solved again
    infiles = ['c0', 'c1', 'x2', 'c3']
    to_solve = ['x2']
    headers = {'c0': header(0, True, tag=0), 'c1': header(5, True, tag=1), 'x2': header(6, False),
               'c3': header(30, True, tag=3)}
    jobs = am.neighbour_jobs(to_solve, [failed(7)], infiles, read_header=headers.get)
    assert len(jobs) == 1 and jobs[0][0] == 0
    assert jobs[0][1]['CRVAL1'] == pytest.approx(150.001)


def test_no_neighbour_too_far_away_or_without_a_solution():
    frames = ['a', 'b', 'c']
    headers = {'a': header(0, False), 'b': header(1, False), 'c': header(2, False)}
    assert am.neighbour_jobs(frames, [failed(6), failed(6), failed(6)], frames, read_header=headers.get) == []
    # more than NEIGHBOUR_MAX_GAP away
    far = {'a': {'DATE-OBS': '2026-10-13T00:00:00', **WCS}, 'b': {'DATE-OBS': '2026-10-13T05:00:00'}}
    assert am.neighbour_jobs(['a', 'b'], [OK, failed(6)], ['a', 'b'], read_header=far.get) == []
