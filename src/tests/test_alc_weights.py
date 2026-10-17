"""
Unit tests for SPlightcurve's 'stable' comparison-star weighting.

The legacy chain (fluxthresh -> iter_algorithm -> top-10 cap) picked a different ensemble when its
inputs changed by ~0.3% (smear A/B, Ganymede 20260223). These tests pin the properties the
replacement is meant to have:
  - leave-one-out inverse-variance weights track each star's real noise, and no star can run away
    with the weight by dominating the ALC
  - the weights are a continuous function of the photometry: a small flux change gives a small
    weight change, and star order does not matter
  - saturated, low-flux, too-faint and variable stars get no weight; the distance kernel is smooth
"""
import os
import sys

import numpy as np
import pytest

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)


@pytest.fixture(scope='module')
def spl():
    pytest.importorskip('fitsio')
    pytest.importorskip('matplotlib')
    from diff_photometry import SPlightcurve
    return SPlightcurve


def field(n_frames=300, seed=1):
    """A synthetic night: 60 stars sharing a transparency curve, each with its own white noise.

    Star 0 is the target. Returns flux (star, frame), per-star sigma, RA, Dec."""
    rng = np.random.RandomState(seed)
    t = np.linspace(0, 1, n_frames)
    transparency = 1 - 0.15 * np.exp(-((t - 0.6) / 0.05) ** 2) + 0.02 * np.sin(9 * t)
    n = 60
    base = np.geomspace(2e5, 4e3, n)
    sigma = 0.002 * np.sqrt(2e5 / base)                     # photon-like: fainter is noisier
    flux = base[:, None] * transparency[None, :] * (1 + sigma[:, None] * rng.randn(n, n_frames))
    ra = 150.0 + rng.uniform(-0.08, 0.08, n)
    dec = -30.0 + rng.uniform(-0.08, 0.08, n)
    ra[0], dec[0] = 150.0, -30.0
    return flux, sigma, ra, dec


def make_lcurves(spl, flux, ra, dec, peak_max=None):
    class Tel:
        rcore = np.full(flux.shape[1], 4.0)
        exposure = np.full(flux.shape[1], 20.0)
        gain = np.full(flux.shape[1], 1.0)

    n, nf = flux.shape
    jd = 2461095.5 + np.arange(nf) / 1440.
    nflux = [spl.normalise(f) for f in flux]
    peak = np.full((n, nf), 1e4)
    if peak_max is not None:
        for i, p in peak_max.items():
            peak[i] = p
    # frames without a SATURATE level: the peak alone against the instrument config's level
    saturated = spl.saturated_stars(peak, None, None, None, 62000.)
    nan = [np.nan] * n
    zeros = [np.zeros(nf)] * n
    return spl.init_all_lcurves(Tel(), ['SP%03d' % i for i in range(n)], nan, nan, nflux, [f.copy() for f in flux],
                                saturated, jd, jd, jd, list(ra), list(dec), nan, nan, nan, zeros, zeros, nan, nan, nan,
                                nan, 5, 'I+z', '20260223', np.full(nf, 20.0))


def stable(spl, flux, ra, dec, **kw):
    lcurves = make_lcurves(spl, flux, ra, dec, **kw)
    target = lcurves[0]
    target.set_as_target('SP000')
    comp = spl.set_all_weights(np.delete(lcurves, 0))
    comp, num_comp = spl.stable_weights(comp, target)
    w = np.array(spl.extract(comp, 'weight'), dtype=float)
    return comp, w, num_comp


def test_loo_weights_follow_the_noise(spl):
    flux, sigma, _, _ = field()
    comp_flux, comp_sigma = flux[1:], sigma[1:]
    nflux = comp_flux / np.median(comp_flux, axis=1, keepdims=True)
    s, w = spl.loo_inverse_variance(nflux, np.ones(len(nflux)), np.median(comp_flux, axis=1))
    # measured scatter ~ the injected noise (plus the small residual of the rest of the ensemble)
    assert np.allclose(s, comp_sigma, rtol=0.15)
    expected = comp_sigma ** -2 / np.sum(comp_sigma ** -2)
    assert np.allclose(w, expected, rtol=0.3, atol=1e-4)
    assert abs(w.sum() - 1) < 1e-12


def test_no_runaway_onto_one_star(spl):
    """The legacy iteration measured each star against an ALC containing itself. With 10 stars, one
    of them twice as quiet as the rest (inverse-variance share 0.31), that iteration hands the quiet
    star 100% of the weight; in production it then hit the reset that re-admits every star."""
    rng = np.random.RandomState(1)
    t = np.linspace(0, 1, 300)
    sigma = np.full(10, 0.004)
    sigma[0] = 0.002
    nflux = (1 - 0.15 * np.exp(-((t - 0.6) / 0.05) ** 2))[None, :] * (1 + sigma[:, None] * rng.randn(10, 300))
    nflux /= np.median(nflux, axis=1, keepdims=True)

    def self_included(w):
        for _ in range(50):
            alc = w @ nflux
            d = nflux / (alc / np.median(alc))[None, :]
            w = np.std(d / np.median(d, axis=1, keepdims=True), axis=1) ** -2
            w = w / w.sum()
            if w.max() > 0.99:          # one star is the ALC: its own scatter -> 0, so stop before 1/0
                break
        return w

    assert self_included(np.full(10, 0.1))[0] > 0.99          # the legacy failure, reproduced
    _, w = spl.loo_inverse_variance(nflux, np.ones(10), np.ones(10), max_iter=50, tol=0)
    expected = sigma ** -2 / np.sum(sigma ** -2)
    assert 0.5 * expected[0] < w[0] < 1.2 * expected[0]
    assert 1 / np.sum(w ** 2) > 0.8 * (1 / np.sum(expected ** 2))


def test_prior_zero_means_zero_weight(spl):
    flux, _, _, _ = field()
    nflux = flux[1:] / np.median(flux[1:], axis=1, keepdims=True)
    prior = np.ones(len(nflux))
    prior[[0, 5, 9]] = 0
    _, w = spl.loo_inverse_variance(nflux, prior, np.ones(len(nflux)))
    assert np.all(w[[0, 5, 9]] == 0)
    assert np.all(w[prior > 0] > 0)


def test_small_input_change_small_weight_change(spl):
    """The smear A/B in miniature: a 0.3% flux change on half the stars, plus noise at the 1e-4 level."""
    flux, _, ra, dec = field()
    _, w0, _ = stable(spl, flux, ra, dec)
    rng = np.random.RandomState(7)
    pert = flux.copy()
    pert[::2] *= 1.003
    pert *= 1 + 1e-4 * rng.randn(*pert.shape)
    _, w1, _ = stable(spl, pert, ra, dec)
    assert np.sum(np.minimum(w0, w1)) > 0.97
    assert np.max(np.abs(w0 - w1)) < 0.01


def test_star_order_does_not_matter(spl):
    flux, _, ra, dec = field()
    _, w0, _ = stable(spl, flux, ra, dec)
    order = np.concatenate([[0], 1 + np.random.RandomState(3).permutation(len(flux) - 1)])
    _, w1, _ = stable(spl, flux[order], ra[order], dec[order])
    back = np.empty_like(w1)
    back[order[1:] - 1] = w1
    assert np.allclose(back, w0, rtol=1e-6, atol=1e-12)


def test_unusable_stars_get_no_weight(spl):
    flux, _, ra, dec = field()
    flux = flux.copy()
    flux[7] *= 1 + 0.3 * (np.arange(flux.shape[1]) % 2)        # wildly variable: preclean drops it
    flux[11] *= np.median(flux[0]) / np.median(flux[11]) * 10 ** (-0.4 * 3.5)   # 3.5 mag below the target
    comp, w, num_comp = stable(spl, flux, ra, dec, peak_max={3: 65000})   # star 3 saturates
    ids = spl.extract(comp, 'id')
    for sp in ('SP003', 'SP007', 'SP011'):
        assert w[ids.index(sp)] == 0
    assert abs(w.sum() - 1) < 1e-9
    assert 1 <= num_comp <= np.count_nonzero(w)


def test_distance_kernel(spl):
    sep = np.radians(np.array([0., 1., 2., 4.]) / 60.)
    k = spl.distance_kernel(sep, 2.0)
    assert np.allclose(k, [1., 0.8, 0.5, 0.2])
    assert np.all(spl.distance_kernel(sep, None) == 1)
    assert np.all(np.diff(spl.distance_kernel(np.radians(np.linspace(0, 0.2, 50)), 2.0)) < 0)
