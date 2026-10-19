"""
Unit tests for astrom.known_scale: T7's plate solve with the plate scale held fixed, and the verification every
solution has to pass.

The fields are synthetic: a catalogue scattered over the solve's box, projected through a known WCS (any rotation,
either parity, the header pointing a few arcmin off the field centre), with the detections' brightness order
scrambled, a fifth of the stars missing and spurious sources mixed in. Covers:
  - the solve recovers the true WCS, sparse and crowded, at any rotation and parity
  - nothing is ever accepted for detections unrelated to the catalogue, or for star trails
  - an inherited WCS is used when it is right and ignored when it is wrong
  - verify()'s counts and chance probability, and the acceptance thresholds
  - the pieces: pair rotation votes, offset votes, the Poisson tail, elongation

Run:  pytest src/tests/test_known_scale.py -q
"""
import math
import os
import sys

import numpy as np
import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from astrom import known_scale as ks        # noqa: E402

SCALE = 0.348 / 3600.0          # degrees per pixel (SSO ANDOR)
SHAPE = (2044, 2046)            # rows, columns of a trimmed ANDOR frame
CENTRE = (150.0, -30.0)         # header pointing


def field(seed, n_cat=400, theta=181.0, parity=1, offset_px=(150.0, -90.0), keep=0.8, spurious=20, jitter=0.3,
          mag_scatter=0.3, shuffle_top=0, shape=SHAPE, centre=CENTRE):
    """A synthetic frame: (sources brightest first, catalogue ra, dec brightest first, inner-box mask, true WCS)."""
    rng = np.random.default_rng(seed)
    rows, cols = shape
    fov = np.array([cols * SCALE / math.cos(math.radians(centre[1])), rows * SCALE])
    # catalogue over the solve's box, in standard coordinates about the header pointing
    half = ks.CATALOGUE_BOX * np.array([cols, rows]) * SCALE / 2
    xi = rng.uniform(-half[0], half[0], n_cat)
    eta = rng.uniform(-half[1], half[1], n_cat)
    ra, dec = _from_standard(xi, eta, *centre)
    j = np.sort(16.0 - 7.0 * rng.random(n_cat) ** 2)
    # the true field centre sits offset_px from the header pointing
    truth = ks.linear_wcs(centre[0], centre[1], theta, parity, SCALE,
                          (cols / 2.0 - offset_px[0], rows / 2.0 - offset_px[1]))
    x, y = truth.all_world2pix(ra, dec, 0)
    inside = (x > 10) & (x < cols - 10) & (y > 10) & (y < rows - 10)
    mags = j[inside] + rng.normal(0, mag_scatter, inside.sum())       # filter and saturation scramble the order
    pts = np.column_stack([x[inside], y[inside]]) + rng.normal(0, jitter, (inside.sum(), 2))
    pts = pts[np.argsort(mags)]
    pts = pts[rng.random(len(pts)) < keep]                              # blends, edges, missed stars
    if shuffle_top:
        pts[:shuffle_top] = pts[:shuffle_top][rng.permutation(shuffle_top)]
    junk = np.column_stack([rng.uniform(10, cols - 10, spurious), rng.uniform(10, rows - 10, spurious)])
    for k, p in zip(rng.integers(0, len(pts) + 1, spurious), junk):     # spurious sources at random ranks
        pts = np.insert(pts, min(k, len(pts)), p, axis=0)
    inner = (np.abs(ra - centre[0]) <= ks.INNER_BOX * fov[0] / 2) & (np.abs(dec - centre[1]) <= ks.INNER_BOX * fov[1] / 2)
    return pts, ra, dec, inner, truth


def _from_standard(xi, eta, ra0, dec0):
    xi, eta = np.radians(xi), np.radians(eta)
    ra0, dec0 = math.radians(ra0), math.radians(dec0)
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        dec = np.arcsin(np.cos(c) * math.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * math.cos(dec0) / rho, 0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * math.cos(dec0) * np.cos(c) - eta * math.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360.0, np.degrees(dec)


def max_error_arcsec(wcs, truth, shape=SHAPE):
    rows, cols = shape
    xs = np.array([0, cols - 1, 0, cols - 1, cols / 2.0])
    ys = np.array([0, 0, rows - 1, rows - 1, rows / 2.0])
    return float(np.max(wcs.pixel_to_world(xs, ys).separation(truth.pixel_to_world(xs, ys)).arcsec))


# ---------------------------------------------------------------------------------------------------------------
# The solve
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('theta,parity', [(181.0, 1), (0.4, 1), (37.0, 1), (270.2, -1), (359.7, -1)])
def test_recovers_the_true_wcs_at_any_rotation_and_parity(theta, parity):
    src, ra, dec, inner, truth = field(1, theta=theta, parity=parity)
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert out['success'], out
    assert out['method'] in ('pairs', 'grid')
    assert max_error_arcsec(out['wcs'], truth) < 0.3
    assert ks.verified(out['verification'])
    assert out['matches'] >= 40


def test_header_pointing_far_from_the_field_centre():
    # 450 px (2.6 arcmin) off: beyond the 1.2-field box that used to be queried, inside the 2-field one
    src, ra, dec, inner, truth = field(2, offset_px=(450.0, 300.0))
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert out['success'], out
    assert max_error_arcsec(out['wcs'], truth) < 0.3


def test_crowded_field_with_saturated_brightest_stars():
    # a bulge-like field: 1,500 stars in the frame, the 40 brightest in random order, as saturation leaves them
    src, ra, dec, inner, truth = field(3, n_cat=6000, shuffle_top=40, mag_scatter=0.6, spurious=100)
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert out['success'], out
    assert max_error_arcsec(out['wcs'], truth) < 0.3


def test_sparse_field():
    src, ra, dec, inner, truth = field(4, n_cat=60, spurious=3)
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert out['success'], out
    assert max_error_arcsec(out['wcs'], truth) < 0.5


def test_seven_star_field():
    # a sparse SPIRIT field (Callisto 20250416 Sp1153+0659): seven sources, seven catalogue stars in the frame. A full
    # TAN fit to that few stars came out 2-3% anisotropic from noise alone; the similarity fit has none.
    shape = (1280, 1024)
    hits = 0
    for seed in range(50, 56):
        src, ra, dec, inner, truth = field(seed, n_cat=30, keep=1.0, spurious=0, shape=shape, theta=179.5)
        if not 6 <= len(src) <= 9:
            continue
        out = ks.solve(src, ra, dec, CENTRE, SCALE, shape, inner=inner)
        assert out['success'], (seed, len(src), out)
        assert max_error_arcsec(out['wcs'], truth, shape) < 1.0
        hits += 1
    assert hits >= 2


def test_spirit_shaped_frame():
    shape = (1280, 1024)
    src, ra, dec, inner, truth = field(5, n_cat=150, shape=shape, theta=179.5)
    out = ks.solve(src, ra, dec, CENTRE, SCALE, shape, inner=inner)
    assert out['success'], out
    assert max_error_arcsec(out['wcs'], truth, shape) < 0.3


@pytest.mark.parametrize('seed', range(5))
def test_detections_unrelated_to_the_catalogue_are_never_accepted(seed):
    _, ra, dec, inner, _ = field(10 + seed, n_cat=800)
    rng = np.random.default_rng(100 + seed)
    noise = np.column_stack([rng.uniform(10, SHAPE[1] - 10, 300), rng.uniform(10, SHAPE[0] - 10, 300)])
    out = ks.solve(noise, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert not out['success'], out
    assert out['wcs'] is None
    assert out['error'] == 'no verified solution'


def test_star_trails_are_never_accepted():
    # a telescope that stopped tracking: each star a 300-px trail, detected as beads along it
    src, ra, dec, inner, truth = field(20, n_cat=800)
    rng = np.random.default_rng(21)
    beads = []
    for x, y in src[:60]:
        for s in np.linspace(0.0, 300.0, 6):
            beads.append((x + s * 0.96 + rng.normal(0, 1.5), y + s * 0.28 + rng.normal(0, 1.5)))
    beads = np.array([b for b in beads if 10 < b[0] < SHAPE[1] - 10 and 10 < b[1] < SHAPE[0] - 10])
    beads = beads[rng.permutation(len(beads))]
    out = ks.solve(beads, ra, dec, CENTRE, SCALE, SHAPE, inner=inner)
    assert not out['success'], out


def test_too_few_sources_or_catalogue_stars():
    src, ra, dec, inner, _ = field(6)
    assert ks.solve(src[:4], ra, dec, CENTRE, SCALE, SHAPE)['error'].startswith('too few sources')
    assert ks.solve(src, ra[:4], dec[:4], CENTRE, SCALE, SHAPE)['error'].startswith('too few catalogue')


def test_a_right_prior_is_refined_and_kept():
    src, ra, dec, inner, truth = field(7)
    rows, cols = SHAPE
    prior = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.3, 1, SCALE, (cols / 2.0 - 150.0 + 4.0, rows / 2.0 + 90.0 - 3.0))
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner, prior=prior)
    assert out['success'] and out['method'] == 'prior', out
    assert max_error_arcsec(out['wcs'], truth) < 0.3


def test_a_prior_sixty_pixels_out_is_recentred_and_kept():
    src, ra, dec, inner, truth = field(12)
    rows, cols = SHAPE
    prior = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.0, 1, SCALE, (cols / 2.0 - 150.0 - 60.0, rows / 2.0 + 90.0 + 45.0))
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner, prior=prior)
    assert out['success'] and out['method'] == 'prior', out
    assert max_error_arcsec(out['wcs'], truth) < 0.3


def test_a_wrong_prior_is_ignored():
    # an Astra WCS left unturned by the old astra_transform: 90 degrees out
    src, ra, dec, inner, truth = field(8)
    rows, cols = SHAPE
    prior = ks.linear_wcs(CENTRE[0], CENTRE[1], 91.0, 1, SCALE, (cols / 2.0, rows / 2.0))
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner, prior=prior)
    assert out['success'] and out['method'] in ('pairs', 'grid'), out
    assert max_error_arcsec(out['wcs'], truth) < 0.3


def test_a_prior_at_the_wrong_scale_is_not_tried():
    src, ra, dec, inner, truth = field(9)
    rows, cols = SHAPE
    prior = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.0, 1, SCALE * 30, (cols / 2.0, rows / 2.0))
    out = ks.solve(src, ra, dec, CENTRE, SCALE, SHAPE, inner=inner, prior=prior)
    assert out['success'] and out['method'] != 'prior'


# ---------------------------------------------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------------------------------------------

def test_verify_counts_the_true_wcs_as_verified():
    src, ra, dec, inner, truth = field(30)
    v = ks.verify(truth, src, ra, dec, SHAPE)
    assert v['tested'] == min(ks.VERIFY_SOURCES, len(src))
    assert v['matched'] >= 0.7 * v['tested']     # all but the spurious sources among the brightest
    assert v['in_frame'] > 50
    assert v['expected'] < 1.0
    assert v['log10_chance'] < -50
    assert ks.verified(v)


def test_verify_rejects_a_wcs_off_by_fifty_pixels():
    src, ra, dec, inner, truth = field(31)
    rows, cols = SHAPE
    shifted = ks.linear_wcs(CENTRE[0], CENTRE[1], 181.0, 1, SCALE, (cols / 2.0 - 150.0 + 50.0, rows / 2.0 + 90.0))
    v = ks.verify(shifted, src, ra, dec, SHAPE)
    assert not ks.verified(v)


def test_verified_thresholds():
    base = dict(matched=20, tested=100, in_frame=60, expected=0.1, log10_chance=-30.0)
    assert ks.verified(base)
    assert not ks.verified(dict(base, matched=7, log10_chance=-12.0))             # too few
    assert not ks.verified(dict(base, log10_chance=-9.8))                         # chance too likely
    assert not ks.verified(dict(base, matched=12, in_frame=600))                 # 12 of 100 possible: 12%
    assert ks.verified(dict(base, matched=8, tested=20, in_frame=12, log10_chance=-15.0))   # sparse: 8 of 12
    assert not ks.verified(dict(base, expected=None))
    assert not ks.verified(None)


def test_sparse_fields_verify_on_five_nearly_complete_matches():
    # Callisto 20251030 Sp0050-2008: 14 sources, 5 catalogue stars in the frame, all 5 matched
    sparse = dict(matched=5, tested=14, in_frame=5, expected=0.0015, log10_chance=-16.2)
    assert ks.verified(sparse)
    assert ks.verified(dict(sparse, tested=6, in_frame=15))             # 5 of 6 sources (Sp2151-1141)
    assert not ks.verified(dict(sparse, matched=4, in_frame=4))         # four is never enough
    assert not ks.verified(dict(sparse, in_frame=7))                    # 5 of 7: 71%
    # Artemis 20240211 after the pointing walked off: 5 of 92 sources on 35 in-frame stars, chance 10^-10.4
    assert not ks.verified(dict(matched=5, tested=92, in_frame=35, expected=0.08, log10_chance=-10.4), trail=0.1)
    assert not ks.verified(dict(sparse, log10_chance=-8.0))             # the chance test still applies


def test_round_sources_verify_on_a_smaller_share():
    # Artemis 20220331: 12 of the brightest 100 on 56 in-frame catalogue stars (21%), at chance 10^-40
    artemis = dict(matched=12, tested=100, in_frame=56, expected=0.004, log10_chance=-40.0)
    assert not ks.verified(artemis)                  # trail unknown: the full share
    assert ks.verified(artemis, trail=0.3)           # round sources
    assert not ks.verified(artemis, trail=0.76)      # Ganymede 20250716's trails
    # trailed frames' best wrong candidates: 8-19% of the possible matches, far beyond chance
    trailed = dict(matched=19, tested=100, in_frame=939, expected=0.6, log10_chance=-21.1)
    assert not ks.verified(trailed, trail=0.77)
    assert not ks.verified(dict(artemis, matched=8), trail=0.3)      # 14%: below even the smaller share


@pytest.mark.parametrize('k,lam', [(1, 0.5), (5, 0.1), (8, 3.0), (30, 2.0), (3, 7.0)])
def test_poisson_tail_matches_the_direct_sum(k, lam):
    direct = sum(math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1)) for i in range(k, k + 300))
    assert ks._log10_poisson_tail(k, lam) == pytest.approx(math.log10(direct), abs=1e-6)


def test_poisson_tail_far_beyond_chance_stays_finite():
    v = ks._log10_poisson_tail(100, 0.03)
    assert -400 < v < -200


# ---------------------------------------------------------------------------------------------------------------
# The pieces
# ---------------------------------------------------------------------------------------------------------------

def test_offset_vote_finds_the_translation():
    rng = np.random.default_rng(40)
    cat = rng.uniform(-1500, 1500, (300, 2))
    src = cat[:80] + np.array([317.0, -42.0]) + rng.normal(0, 0.3, (80, 2))
    t, votes = ks.offset_vote(src, cat)
    assert t == pytest.approx([317.0, -42.0], abs=0.5)
    assert votes >= 75
    # limited to a window that excludes the true offset, the vote cannot find it
    t2, votes2 = ks.offset_vote(src, cat, limit=200.0)
    assert votes2 < 10
    assert ks.offset_vote(src, cat + 5000.0, limit=200.0) == (None, 0)


def test_rotation_votes_find_the_angle_modulo_180():
    rng = np.random.default_rng(41)
    cat = rng.uniform(-1000, 1000, (90, 2))
    src = (cat[:30] @ ks._rotation(213.0).T) + 900 + rng.normal(0, 0.3, (30, 2))
    peaks = ks.rotation_votes(src, cat)
    assert peaks and peaks[0] == pytest.approx(213.0 - 180.0, abs=0.3)


@pytest.mark.parametrize('parity', [1, -1])
def test_fit_similarity_recovers_a_similarity(parity):
    truth = ks.linear_wcs(CENTRE[0], CENTRE[1], 47.0, parity, SCALE, (900.0, 1100.0))
    rng = np.random.default_rng(42)
    px = rng.uniform(0, 2000, (8, 2))
    sky = truth.pixel_to_world(px[:, 0], px[:, 1])
    w = ks.fit_similarity(px, sky.ra.deg, sky.dec.deg, parity)
    assert ks.parity_of(w) == parity == ks.parity_of(truth)
    back = np.column_stack(w.world_to_pixel(sky))
    assert back == pytest.approx(px, abs=1e-3)
    assert ks.geometry(w)[1] == pytest.approx(0.0, abs=1e-9)


def test_linear_wcs_round_trip():
    w = ks.linear_wcs(10.0, 20.0, 33.0, -1, SCALE, (500.0, 600.0))
    x, y = w.all_world2pix([10.0], [20.0], 0)
    assert (x[0], y[0]) == pytest.approx((500.0, 600.0), abs=1e-6)
    assert ks.solved_scale(w) == pytest.approx(SCALE)


def _render(shape, points, trail=None, sigma=1.8, peak=300.0, seed=0):
    rng = np.random.default_rng(seed)
    img = rng.normal(0.0, 3.0, shape)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    for x, y in points:
        steps = [(0.0, 0.0)] if trail is None else [(s * trail[0], s * trail[1]) for s in np.linspace(-1, 1, 81)]
        for dx, dy in steps:
            img += peak / len(steps) * 10 * np.exp(-((xx - x - dx) ** 2 + (yy - y - dy) ** 2) / (2 * sigma ** 2))
    return img


def test_trail_ellipticity_round_stars_trails_and_crowding():
    pts = [(40 + 60 * i, 40 + 50 * k) for i in range(5) for k in range(5)]
    round_ = ks.trail_ellipticity(_render((300, 340), pts), np.array(pts, float))
    trails = ks.trail_ellipticity(_render((300, 340), pts, trail=(25.0, 6.0)), np.array(pts, float))
    assert round_ < 0.1
    assert trails > 0.5
    assert ks.axis_ratio(trails) > 1.7
    # crowding: each star with a close neighbour in a random direction stretches it, but not all one way
    rng = np.random.default_rng(3)
    pairs = []
    for x, y in pts:
        a = rng.uniform(0, 2 * np.pi)
        pairs += [(x, y), (x + 5 * np.cos(a), y + 5 * np.sin(a))]
    crowded = ks.trail_ellipticity(_render((300, 340), pairs), np.array(pts, float))
    assert crowded < 0.3
    assert ks.trail_ellipticity(np.zeros((50, 50)), np.array([[25.0, 25.0]])) is None
