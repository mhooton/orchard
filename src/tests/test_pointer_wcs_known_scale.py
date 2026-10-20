"""
End-to-end tests of T7's plate solve (pointer_wcs) with the known-scale solver, on rendered frames.

Each test renders a small ANDOR-like frame of Gaussian stars at the positions a known WCS gives a synthetic
catalogue, writes that catalogue into a one-shard stand-in for the local Gaia/2MASS database, and runs pointer_wcs
on the frame. Covers:
  - a frame solves, and the WCS written into the header is the true one
  - an inherited WCS that is wrong (Astra's unturned matrix) is replaced; a right one is refined and kept
  - a frame of stars unrelated to the catalogue fails, with no WCS left behind and the reason given
  - star trails (a telescope that stopped tracking) are named as such
  - the old twirl-only solve is still there behind known_scale_solve=False

Needs the pipeline's own environment (twirl, photutils): run inside orchard-server.
Run:  pytest src/tests/test_pointer_wcs_known_scale.py -q
"""
import math
import os
import sqlite3
import sys

import numpy as np
import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pytest.importorskip('twirl')
pytest.importorskip('photutils')
from astropy.io import fits                     # noqa: E402
from astropy.wcs import WCS                     # noqa: E402

from astrom import known_scale as ks            # noqa: E402
from astrom import pointer_wcs as pw            # noqa: E402

SCALE = 13.5e-6 / 8.0 * 206264.806 / 3600.0     # XPIXSZ 13.5 um, FOCALLEN 8000 mm, degrees per pixel
SHAPE = (600, 640)
CENTRE = (150.0, -30.5)
TRUTH_OFFSET = (40.0, -25.0)                    # field centre from the header pointing, pixels
WCS_KEYWORDS = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                'CUNIT1', 'CUNIT2', 'WCSAXES', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2', 'LONPOLE',
                'LATPOLE', 'RADESYS']


def truth_wcs(theta=181.0, parity=1):
    rows, cols = SHAPE
    return ks.linear_wcs(CENTRE[0], CENTRE[1], theta, parity, SCALE,
                         (cols / 2.0 - TRUTH_OFFSET[0], rows / 2.0 - TRUTH_OFFSET[1]))


def catalogue(seed=0, n=260):
    rng = np.random.default_rng(seed)
    half = 1.1 * ks.CATALOGUE_BOX * np.array([SHAPE[1] / math.cos(math.radians(CENTRE[1])), SHAPE[0]]) * SCALE / 2
    ra = CENTRE[0] + rng.uniform(-half[0], half[0], n)
    dec = CENTRE[1] + rng.uniform(-half[1], half[1], n)
    j = 16.0 - 6.0 * rng.random(n) ** 2
    return ra, dec, j


def make_db(tmp_path, ra, dec, j):
    path = os.path.join(str(tmp_path), 'gaia.db')
    conn = sqlite3.connect(path)
    conn.execute('CREATE TABLE `-31_-30` (ra, dec, pmra, pmdec, phot_g_mean_mag, j_m)')
    conn.executemany('INSERT INTO `-31_-30` VALUES (?,?,?,?,?,?)',
                     [(float(a), float(d), 0.0, 0.0, float(m) + 1.0, float(m)) for a, d, m in zip(ra, dec, j)])
    conn.commit()
    conn.close()
    return path


def render(points, mags, trail=None, seed=1):
    rng = np.random.default_rng(seed)
    rows, cols = SHAPE
    img = 1000.0 + rng.normal(0.0, 10.0, SHAPE)
    yy, xx = np.mgrid[-12:13, -12:13]
    for (x, y), m in zip(points, mags):
        peak = 4000.0 * 10 ** (-0.4 * (m - 11.0))
        steps = [(0.0, 0.0)] if trail is None else [(s * trail[0], s * trail[1]) for s in np.linspace(-1, 1, 61)]
        for dx, dy in steps:
            cx, cy = x + dx, y + dy
            ix, iy = int(round(cx)), int(round(cy))
            if ix < 12 or iy < 12 or ix >= cols - 12 or iy >= rows - 12:
                continue
            stamp = peak / len(steps) * np.exp(-((xx + ix - cx) ** 2 + (yy + iy - cy) ** 2) / (2 * 1.8 ** 2))
            img[iy - 12:iy + 13, ix - 12:ix + 13] += stamp
    return img.astype(np.float32)


def write_frame(tmp_path, img, extra=None, name='procTEST.fits'):
    d = tmp_path / 'PipelineOutput' / 'v3' / 'Io' / 'output' / '20260101' / 'Sp0000-0000' / '1'
    os.makedirs(str(d), exist_ok=True)
    path = str(d / name)
    h = fits.PrimaryHDU(data=img)
    h.header.update(TELESCOP='SPECULOOS-IO', IMAGETYP='Light Frame', RA=CENTRE[0], DEC=CENTRE[1], XPIXSZ=13.5,
                    FOCALLEN=8000.0, **{'DATE-OBS': '2026-01-01T03:00:00'})
    if extra is not None:
        h.header.update(extra)
    h.writeto(path)
    return path


def frame_and_db(tmp_path, theta=181.0, parity=1, trail=None, unrelated=False, inherited=None):
    ra, dec, j = catalogue()
    db = make_db(tmp_path, ra, dec, j)
    truth = truth_wcs(theta, parity)
    x, y = truth.all_world2pix(ra, dec, 0)
    pts = np.column_stack([x, y])
    if unrelated:
        rng = np.random.default_rng(7)
        pts = np.column_stack([rng.uniform(15, SHAPE[1] - 15, len(pts)), rng.uniform(15, SHAPE[0] - 15, len(pts))])
    img = render(pts, j, trail=trail)
    extra = inherited.to_header() if inherited is not None else None
    return write_frame(tmp_path, img, extra), db, truth


def header_error_arcsec(path, truth):
    w = WCS(fits.getheader(path), naxis=2)
    rows, cols = SHAPE
    xs = np.array([0, cols - 1, 0, cols - 1, cols / 2.0])
    ys = np.array([0, 0, rows - 1, rows - 1, rows / 2.0])
    return float(np.max(w.pixel_to_world(xs, ys).separation(truth.pixel_to_world(xs, ys)).arcsec))


def test_a_frame_solves_and_the_header_holds_the_true_wcs(tmp_path):
    path, db, truth = frame_and_db(tmp_path)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert result['success'], result
    assert result['method'] in ('pairs', 'grid')
    assert result['matches'] >= ks.VERIFY_MIN_MATCHES
    assert header_error_arcsec(path, truth) < 0.5


@pytest.mark.parametrize('theta,parity', [(0.4, 1), (269.0, -1)])
def test_other_orientations(tmp_path, theta, parity):
    path, db, truth = frame_and_db(tmp_path, theta=theta, parity=parity)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert result['success'], result
    assert header_error_arcsec(path, truth) < 0.5


def test_a_wrong_inherited_wcs_is_replaced(tmp_path):
    rows, cols = SHAPE
    wrong = ks.linear_wcs(CENTRE[0], CENTRE[1], 91.0, 1, SCALE, (cols / 2.0, rows / 2.0))
    path, db, truth = frame_and_db(tmp_path, inherited=wrong)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert result['success'] and result['method'] != 'prior', result
    assert header_error_arcsec(path, truth) < 0.5


def test_a_right_inherited_wcs_is_refined_and_kept(tmp_path):
    rows, cols = SHAPE
    near = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.2, 1, SCALE,
                         (cols / 2.0 - TRUTH_OFFSET[0] + 3.0, rows / 2.0 - TRUTH_OFFSET[1] - 2.0))
    path, db, truth = frame_and_db(tmp_path, inherited=near)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert result['success'] and result['method'] == 'prior', result
    assert header_error_arcsec(path, truth) < 0.5


def test_unrelated_stars_fail_and_leave_no_wcs(tmp_path):
    rows, cols = SHAPE
    inherited = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.0, 1, SCALE, (cols / 2.0, rows / 2.0))
    path, db, _ = frame_and_db(tmp_path, unrelated=True, inherited=inherited)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True, verify_twirl=True)
    assert not result['success'], result
    assert result['error'].startswith('No verified solution'), result['error']
    assert 'twirl:' in result['error']
    assert 'CTYPE1' not in fits.getheader(path)


def test_trailed_frames_are_named(tmp_path):
    # a telescope that stopped tracking: 300-px trails, which no solver verifies
    path, db, _ = frame_and_db(tmp_path, trail=(150.0, 40.0))
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert not result['success'], result
    assert result['error'].startswith('Sources trailed'), result['error']
    assert 'CTYPE1' not in fits.getheader(path)


def test_twirl_only_solve_is_still_available(tmp_path):
    path, db, truth = frame_and_db(tmp_path)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True,
                            known_scale_solve=False, verify_twirl=False)
    if result['success']:
        assert result['method'] == 'twirl'
        assert header_error_arcsec(path, truth) < 2.0
    else:
        assert result['error']


def test_a_known_scale_crash_falls_back_to_twirl(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError('simulated')
    monkeypatch.setattr(pw.known_scale, 'solve', boom)
    path, db, truth = frame_and_db(tmp_path)
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    # twirl ran: either it solved, verified, or it failed with its own reason
    assert result['method'] in (None, 'twirl'), result
    if result['success']:
        assert header_error_arcsec(path, truth) < 2.0
    else:
        assert 'simulated' not in (result['error'] or '')


def test_a_faint_frame_solves_from_its_neighbours_solution(tmp_path, monkeypatch):
    # Callisto 20261012's faint SPIRIT frames: 4 sources on catalogue stars and a spurious one cannot verify alone,
    # but do from the solution of a frame taken minutes earlier. Detection is pinned to the five rendered sources,
    # as on those frames ("Sources: 5→5"); the image is real, so their shapes are measured.
    ra, dec, j = catalogue(n=45)
    db = make_db(tmp_path, ra, dec, j)
    truth = truth_wcs()
    x, y = truth.all_world2pix(ra, dec, 0)
    rows, cols = SHAPE
    inside = (x > 20) & (x < cols - 20) & (y > 20) & (y < rows - 20)
    order = np.argsort(j[inside])[:4]
    pts = np.column_stack([x[inside][order], y[inside][order]])
    pts = np.vstack([pts, [[cols * 0.31, rows * 0.77]]])
    path = write_frame(tmp_path, render(pts, np.full(len(pts), 12.0)))
    monkeypatch.setattr(pw, 'find_stars_multiscale', lambda *a, **k: pts.copy())
    alone = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True)
    assert not alone['success'], alone
    nb = truth.deepcopy()
    nb.wcs.crpix = nb.wcs.crpix + np.array([12.0, -8.0])
    cards = {k: v for k, v in nb.to_header().items() if k in WCS_KEYWORDS}
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True, neighbour=cards)
    assert result['success'] and result['method'] == 'neighbour', result
    assert result['matches'] >= 4
    assert header_error_arcsec(path, truth) < 1.0
    # a neighbour 300 px away is no help, and the usual failure comes back
    far = truth.deepcopy()
    far.wcs.crpix = far.wcs.crpix + np.array([300.0, 0.0])
    cards = {k: v for k, v in far.to_header().items() if k in WCS_KEYWORDS}
    result = pw.pointer_wcs(path, db, wcs_keywords=WCS_KEYWORDS, clear_existing_wcs=True, neighbour=cards)
    assert not result['success'] and result['error'].startswith('No verified solution'), result


def test_inside_query_box_matches_the_query():
    radecs = np.array([[150.0, -30.5], [150.05, -30.5], [150.2, -30.5], [150.0, -30.6]])
    mask = pw.inside_query_box(radecs, CENTRE, (0.12, 0.12))
    assert mask.tolist() == [True, True, False, False]
    assert pw.inside_query_box(np.zeros((0, 2)), CENTRE, (1, 1)).shape == (0,)


def test_inherited_wcs_reads_only_a_sky_wcs():
    h = fits.Header()
    assert pw.inherited_wcs(h) is None
    h.update(truth_wcs().to_header())
    assert pw.inherited_wcs(h) is not None
    h['CRVAL1'] = 0.0
    assert pw.inherited_wcs(h) is None
