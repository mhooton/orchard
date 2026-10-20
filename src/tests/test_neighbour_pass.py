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
from astropy.io import fits              # noqa: E402

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



# ---------------------------------------------------------------------------------------------------------------
# Through astrometry.main, with the WCS cache (astrom/wcs_cache.py): the second pass takes its neighbours from frames
# applied from the cache as well as frames solved this run, and what it solves is cached like any other solution.
# The solves are real; test_wcs_cache's `solver` fakes only the image cleaning, detection and Gaia query.
# ---------------------------------------------------------------------------------------------------------------

import numpy as np                                                      # noqa: E402
import astropy.wcs                                                       # noqa: E402
from astrom import wcs_cache                                            # noqa: E402
from test_wcs_cache import Night, STARS, _truth, solver                 # noqa: E402,F401


def _hard_frame(monkeypatch, path):
    """Make `path` a frame no lone solve verifies, as the faint sparse-field frames of Callisto 20261012 are: alone,
    it fails with its sources counted; given a neighbour, the real solve runs."""
    real = am.pointer_wcs

    def pointer_wcs(filepath, *args, **kwargs):
        if os.path.basename(filepath) == os.path.basename(path) and not kwargs.get('neighbour'):
            return dict(success=False, crpix=None, crval=None, sources_detected=12, sources_used=12, gaia_queried=12,
                        gaia_used=12, matches=3, method=None, stage='star_matching',
                        error='No verified solution (best 3 of 12 sources on catalogue stars); twirl: ...')
        return real(filepath, *args, **kwargs)
    monkeypatch.setattr(am, 'pointer_wcs', pointer_wcs)


def _solved_right(path, center):
    """The frame's WCS puts the detections on the catalogue stars they are (test_wcs_cache._truth), within 0.1"."""
    got = astropy.wcs.WCS(fits.getheader(path)).pixel_to_world(STARS[:, 0], STARS[:, 1])
    want = _truth(center).pixel_to_world(STARS[:, 0], STARS[:, 1])
    return np.max(got.separation(want).arcsec) < 0.1


def _row_is_the_solution(row, path):
    """The frame's cache row holds the cards the solve wrote into it."""
    cached, header = wcs_cache.cards_header(row['cards']), fits.getheader(path)
    return 'CRVAL1' in cached and all(header.get(k) == v for k, v in cached.items())


def test_the_second_pass_solves_from_a_frame_solved_this_run_and_caches_it(tmp_path, solver, monkeypatch, capsys):
    night = Night(tmp_path, n=3)
    hard = night.procs[1]
    _hard_frame(monkeypatch, hard)
    night.t7()
    out = capsys.readouterr().out
    assert 'procIo_001.fits FAILED' in out                     # the first pass
    assert 'Second pass: 1 frames tried again' in out and 'Second pass: 1 of 1 solved' in out
    assert ' procIo_001.fits SUCCESS' in out and 'Method: neighbour' in out
    assert 'Successfully plate solved: 3' in out and 'Failed plate solving: 0' in out
    assert _solved_right(hard, solver[-1])                      # its query is the last one: frames 0 and 2 went first
    rows = night.rows()
    assert set(rows) == {'Io_000.fits', 'Io_001.fits', 'Io_002.fits'}
    assert _row_is_the_solution(rows['Io_001.fits'], hard)
    assert rows['Io_001.fits']['solver_epoch'] == am.SOLVER_EPOCH


def test_the_second_pass_takes_neighbours_from_the_cache_and_caches_what_it_solves(tmp_path, solver, monkeypatch,
                                                                                    capsys):
    night = Night(tmp_path, n=3)
    hard = night.procs[1]
    _hard_frame(monkeypatch, hard)
    # a run without the second pass (v4.1.0): the hard frame fails and gets no row
    real_jobs = am.neighbour_jobs
    monkeypatch.setattr(am, 'neighbour_jobs', lambda *a, **k: [])
    night.t7()
    assert set(night.rows()) == {'Io_000.fits', 'Io_002.fits'}
    # the rerun: frames 0 and 2 come from the cache and only the hard frame is solved, so its only neighbours are
    # frames the cache applied
    monkeypatch.setattr(am, 'neighbour_jobs', real_jobs)
    capsys.readouterr()
    asked = len(solver)
    night.t7()
    out = capsys.readouterr().out
    assert 'WCS cache: 2 frames from the cache, 1 to solve' in out
    assert 'Second pass: 1 frames tried again' in out and 'Second pass: 1 of 1 solved' in out
    assert ' procIo_001.fits SUCCESS' in out and 'Method: neighbour' in out
    assert len(solver) == asked + 1                             # one solve: the hard frame's second pass
    assert _solved_right(hard, solver[-1])
    rows = night.rows()
    assert set(rows) == {'Io_000.fits', 'Io_001.fits', 'Io_002.fits'}
    assert _row_is_the_solution(rows['Io_001.fits'], hard)
    # and the next run solves nothing
    capsys.readouterr()
    night.t7()
    assert 'WCS cache: 3 frames from the cache, 0 to solve' in capsys.readouterr().out
    assert len(solver) == asked + 1
