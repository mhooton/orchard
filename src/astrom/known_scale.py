"""
known_scale.py
==============

T7's plate solve when the plate scale is known, and the test every solution has to pass.

Why: twirl (``pointer_wcs``'s first solver until 2026-10) matches 4-star asterisms with a *free* affine transform,
so a solution at any scale can win. When few of the 16 brightest detections are among the 32 brightest catalogue
stars, a transform that squeezes the whole catalogue into a ~40 px clump around a handful of detections outscores
the true one, and the plate-scale sanity check then rejects the frame ("Plate scale sanity check failed: solved=38
arcsec/px ... Matches: 0"). Measured over every pipeline log on appct (2026-10-19), that was 39,765 frames whose
last attempt failed on frames with 20 or more sources, the rate rising with field density from 0.15-0.42% below 50
catalogue stars in the box to 18% above 4,000. In bulge fields twirl's plausible-scale solutions are often wrong too.

How: the plate scale is known to ~1% (header or instrument config) and the pixel grid is square, so the sky-to-pixel
map is a rotation, a parity and a translation. The rotation comes from star pairs: a detection pair and a catalogue
pair of the same length vote for the angle between them. Each candidate rotation (both parities, and the opposite
direction, since a pair has no sign) is scored by voting on the offset between every detection and every catalogue
star; the true one gives a sharp peak. In crowded fields the pair votes drown, so a 1-degree rotation grid scored the
same way is the fallback. The winner is refined by iterative cross-matching: a TAN fit (astropy's
fit_wcs_from_points, as twirl's own last step) with 10 or more stars, a similarity below that.

``verified`` is what keeps a wrong solution out, whichever solver made it. Of the brightest 100 detections, 6 or more
must sit within 3 px of a catalogue star, at most 10^-10 likely by chance for the catalogue density in the frame, and
they must be a quarter of what could match (15% on round sources). In a sparse field, 5 that are 80% of what could
match also verify. Measured on 820 frames of 29 nights (mh_scratch/platesolve_partial, 2026-10-19): of 620 frames
from affected nights production made 281 correct solutions and 32 wrong ones, this solve 495 correct and none wrong,
losing none of production's correct ones; on 8 control nights both made 200 of 200. Trails are what the share of
matches guards against: on frames from a telescope that stopped tracking the best wrong candidates matched 8-19% of
what could match, at chance 10^-10 to 10^-21.

Everything here is pure: no file is read or written.
"""
import math

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points, proj_plane_pixel_scales
from scipy.spatial import cKDTree

# Catalogue box for the solve, in fields of view. The header pointing can be the scheduled target rather than the
# field centre, a few arcmin out in some eras; a box twice the field still holds the whole field for offsets up to
# half a field.
CATALOGUE_BOX = 2.0
# The box the pair votes use: the production query box, where catalogue and image overlap most.
INNER_BOX = 1.2

# The test a solution must pass (see verify / verified).
VERIFY_RADIUS = 3.0          # px between a detection and its catalogue star
VERIFY_SOURCES = 100         # brightest detections tested
VERIFY_MIN_MATCHES = 6
VERIFY_MIN_FRACTION = 0.25   # of min(detections tested, catalogue stars in the frame)
VERIFY_MAX_LOG10_CHANCE = -10.0
# On round sources a smaller share will do. Artemis I+z frames match a median 29% of their in-frame J < 16 stars
# (down to 21%) on correct solutions. Star trails are why the share exists: a trail holds each star's position at
# some instant, so on Ganymede 20250716's trailed frames the best wrong candidates reached 8-19% of matches at
# chance 10^-10 to 10^-21; their sources measure 0.73-0.78 on trail_ellipticity.
VERIFY_MIN_FRACTION_ROUND = 0.15
ROUND_BELOW = 0.5            # trail_ellipticity below which a frame's sources count as round
# A field with only 5 or 6 catalogue stars (to J = 16) cannot give 6 matches. Five verify there if they are nearly all
# that could match: two matches fix a similarity, and three more coincidences among 6 stars are negligible. SPIRIT
# has 36,093 frames solved on exactly 5 matches; on the sparse ones checked (Callisto 20251030 Sp0050-2008, 20260918
# Sp2151-1141, 20261008 Sp0218-0617) the old solution and this solve's best candidate agreed, 5 of 5. Four never do.
SPARSE_MIN_MATCHES = 5
SPARSE_MIN_COMPLETENESS = 0.8   # of min(detections tested, catalogue stars in the frame)
# A frame taken close in time to one that solved (solve_from_neighbour) needs fewer matches: the neighbour's solution
# fixes rotation, scale and pointing to within a few pixels, so 3 sources on catalogue stars within this shift of it
# cannot be chance (the offset vote searches a 100 px square). Faint frames of sparse SPIRIT fields have only 5-7
# sources: on Callisto 20261012 v4.1.0 refused 63 frames whose old photometry was normal. Their consecutive frames,
# 6-18 s apart, put 3 of their 5 sources on catalogue stars; twirl's "4 of 5" was its own quad, fitted.
NEIGHBOUR_MIN_MATCHES = 3
NEIGHBOUR_MAX_LOG10_CHANCE = -8.0
NEIGHBOUR_SHIFT = 50.0          # px between this frame's solution and its neighbour's, at the frame centre
NEIGHBOUR_ROTATION = 0.5        # degrees
NEIGHBOUR_SCALE = 0.01          # fractional

# Pair votes
PAIR_SOURCES = 30
PAIR_CATALOGUE_INNER = 90
PAIR_CATALOGUE_OUTER = 120
PAIR_PEAKS = 3
PAIR_LENGTH_TOLERANCE = (2.0, 0.012)   # px, and fraction of the pair length: the scale is known to ~1%
# Offset votes
OFFSET_SOURCES = 100
OFFSET_CATALOGUE = 600
OFFSET_BIN = 6.0
# How far an inherited WCS may be off and still be recentred before refining, px
PRIOR_SEARCH = 200.0
# Rotation grid (fallback)
GRID_STEP = 1.0
GRID_SOURCES = 60
GRID_CATALOGUE = 300
GRID_BIN = 10.0
# A solution must be what a telescope gives: the expected scale (the header's is within 1% of every measured solve,
# the instrument config's within 2%) and square pixels. A fit bent to a handful of chance matches rarely is.
SCALE_TOLERANCE = 0.05
MAX_ANISOTROPY = 0.02        # ratio of the two pixel-scale axes, less one (TAN fits; similarities have none)


def standard_coordinates(ra, dec, ra0, dec0):
    """Gnomonic standard coordinates (xi east, eta north), in degrees, of (ra, dec) about (ra0, dec0)."""
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    ra0, dec0 = math.radians(ra0), math.radians(dec0)
    cosc = math.sin(dec0) * np.sin(dec) + math.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (math.cos(dec0) * np.sin(dec) - math.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi), np.degrees(eta)


def _rotation(theta_deg):
    t = math.radians(theta_deg)
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def _linear_map(theta_deg, parity, scale_deg):
    """Matrix taking standard coordinates (degrees) to pixel offsets."""
    return _rotation(theta_deg) @ np.diag([float(parity), 1.0]) / scale_deg


def _pairs(points):
    i, j = np.triu_indices(len(points), 1)
    v = points[j] - points[i]
    return np.hypot(v[:, 0], v[:, 1]), np.arctan2(v[:, 1], v[:, 0])


def rotation_votes(sources, catalogue_px, n_peaks=PAIR_PEAKS, bin_deg=1.0, tolerance=PAIR_LENGTH_TOLERANCE):
    """Candidate rotations taking catalogue pairs (pixel units, parity applied, unrotated) onto source pairs.

    Every source pair votes, for each catalogue pair of the same length, for the angle between the two. A pair has
    no direction, so angles are folded into [0, 180) and the caller tries both. Returns up to n_peaks angles
    (degrees), strongest first."""
    if len(sources) < 2 or len(catalogue_px) < 2:
        return []
    sl, sa = _pairs(sources)
    cl, ca = _pairs(catalogue_px)
    order = np.argsort(cl)
    cl, ca = cl[order], ca[order]
    tol = tolerance[0] + tolerance[1] * sl
    lo = np.searchsorted(cl, sl - tol)
    hi = np.searchsorted(cl, sl + tol)
    counts = hi - lo
    total = int(counts.sum())
    if total == 0:
        return []
    si = np.repeat(np.arange(len(sl)), counts)
    ci = np.repeat(lo, counts) + (np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts))
    angle = np.degrees(sa[si] - ca[ci]) % 180.0
    nbins = int(round(180.0 / bin_deg))
    hist = np.bincount(np.minimum((angle / bin_deg).astype(int), nbins - 1), minlength=nbins)
    smooth = hist + np.roll(hist, 1) + np.roll(hist, -1)   # a peak on a bin edge is not split
    peaks, taken = [], np.zeros(nbins, bool)
    for idx in np.argsort(smooth)[::-1]:
        if len(peaks) >= n_peaks or smooth[idx] == 0:
            break
        if taken[idx]:
            continue
        centre = (idx + 0.5) * bin_deg
        delta = (angle - centre + 90.0) % 180.0 - 90.0
        near = np.abs(delta) < 1.5 * bin_deg
        peaks.append((centre + float(np.median(delta[near]))) % 180.0 if near.any() else centre)
        taken[(idx + np.arange(-3, 4)) % nbins] = True
    return peaks


def offset_vote(sources, catalogue_px, bin_px=OFFSET_BIN, limit=None):
    """The translation t with sources ~ catalogue_px + t, by voting over every source-catalogue pair (only offsets
    within `limit` px on both axes, if given).

    Returns (t, votes): votes is the count in the best 2x2-bin cell."""
    if len(sources) == 0 or len(catalogue_px) == 0:
        return None, 0
    off = (sources[:, None, :] - catalogue_px[None, :, :]).reshape(-1, 2)
    if limit is not None:
        off = off[(np.abs(off) < limit).all(1)]
        if len(off) == 0:
            return None, 0
    lo = off.min(0) - bin_px
    nx = int((off[:, 0].max() - lo[0]) // bin_px) + 2
    ny = int((off[:, 1].max() - lo[1]) // bin_px) + 2
    ix = ((off[:, 0] - lo[0]) // bin_px).astype(np.int64)
    iy = ((off[:, 1] - lo[1]) // bin_px).astype(np.int64)
    hist = np.bincount(ix * ny + iy, minlength=nx * ny).reshape(nx, ny)
    box = hist[:-1, :-1] + hist[1:, :-1] + hist[:-1, 1:] + hist[1:, 1:]
    kx, ky = np.unravel_index(int(np.argmax(box)), box.shape)
    sel = (ix >= kx) & (ix <= kx + 1) & (iy >= ky) & (iy <= ky + 1)
    return np.median(off[sel], axis=0), int(box[kx, ky])


def linear_wcs(ra0, dec0, theta_deg, parity, scale_deg, t):
    """TAN WCS with pixel = A @ (xi, eta) + t about (ra0, dec0); pixels 0-based, as the detections are."""
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [t[0] + 1.0, t[1] + 1.0]
    w.wcs.cd = np.linalg.inv(_linear_map(theta_deg, parity, scale_deg))
    return w


def _unique_matches(sources, catalogue_px, radius):
    """Pairs (source, catalogue star) within radius, at most one source per catalogue star (the closest)."""
    d, j = cKDTree(catalogue_px).query(sources)
    i = np.flatnonzero(d < radius)
    if len(i) == 0:
        return i, i
    j = j[i]
    order = np.argsort(d[i])
    _, first = np.unique(j[order], return_index=True)
    keep = order[first]
    return i[keep], j[keep]


def _recentred(wcs, sources, ra, dec, search=PRIOR_SEARCH):
    """A copy of `wcs` shifted by the offset most source-catalogue pairs agree on, if one stands out within
    `search` px: a starting WCS a few tens of pixels out still refines."""
    px = _project(wcs, ra[:OFFSET_CATALOGUE], dec[:OFFSET_CATALOGUE])
    t, votes = offset_vote(sources[:OFFSET_SOURCES], px, limit=search)
    out = wcs.deepcopy()
    if t is not None and votes >= 4:
        out.wcs.crpix = out.wcs.crpix + t
    return out


def _project(wcs, ra, dec):
    x, y = wcs.all_world2pix(ra, dec, 0)
    px = np.column_stack([x, y])
    px[~np.isfinite(px).all(1)] = 1e9
    return px


# Below this many matched stars a full TAN fit (six parameters) is too loose: on SPIRIT frames with seven stars it
# gave 2-3% anisotropy from noise alone. Fewer matches get a similarity (scale, rotation, shift; square pixels).
FULL_FIT_MATCHES = 10


def parity_of(wcs):
    """+1 or -1: the sign of the pixel-to-sky matrix's determinant, in the convention of _linear_map."""
    return 1 if np.linalg.det(wcs.pixel_scale_matrix) > 0 else -1


def fit_similarity(pixels, ra, dec, parity):
    """Least-squares WCS with pixel = s R (parity * xi, eta) + t about the matched stars' mean position."""
    ra0, dec0 = float(np.mean(ra)), float(np.mean(dec))
    xi, eta = standard_coordinates(ra, dec, ra0, dec0)
    u1, u2 = parity * xi, eta
    n = len(pixels)
    design = np.zeros((2 * n, 4))
    design[0::2] = np.column_stack([u1, -u2, np.ones(n), np.zeros(n)])
    design[1::2] = np.column_stack([u2, u1, np.zeros(n), np.ones(n)])
    (a, b, tx, ty), *_ = np.linalg.lstsq(design, np.asarray(pixels, float).reshape(-1), rcond=None)
    matrix = np.array([[a, -b], [b, a]]) @ np.diag([float(parity), 1.0])
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [tx + 1.0, ty + 1.0]
    w.wcs.cd = np.linalg.inv(matrix)
    return w


def _fit(src, sky, ra, dec, i, j, parity):
    if len(i) >= FULL_FIT_MATCHES:
        return fit_wcs_from_points(src[i].T, sky[j], projection='TAN')
    return fit_similarity(src[i], ra[j], dec[j], parity)


def refine(sources, ra, dec, wcs, radii=(10.0, 4.0, 3.0), final_radius=2.0, max_sources=200):
    """Iterative cross-match and fit from a starting WCS: a TAN fit (astropy's fit_wcs_from_points, twirl's own last
    step) with FULL_FIT_MATCHES or more stars, a similarity below. A last pass at final_radius is kept only if it
    loses no more than a tenth of the stars. Returns (wcs, matches) or (None, matches)."""
    sky = SkyCoord(ra, dec, unit='deg')
    src = sources[:max_sources]
    parity = parity_of(wcs)
    n = 0
    for r in radii:
        i, j = _unique_matches(src, _project(wcs, ra, dec), r)
        n = len(i)
        if n < 4:
            return None, n
        try:
            wcs = _fit(src, sky, ra, dec, i, j, parity)
        except Exception:
            return None, n
    i, j = _unique_matches(src, _project(wcs, ra, dec), final_radius)
    if len(i) >= max(4, 0.9 * n):
        try:
            wcs, n = _fit(src, sky, ra, dec, i, j, parity), len(i)
        except Exception:
            pass
    return wcs, n


def _log10_poisson_tail(k, lam):
    """log10 P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 0.0
    if lam <= 0:
        return -300.0
    term = -lam + k * math.log(lam) - math.lgamma(k + 1)
    terms, i = [term], k
    while i < k + 2000:
        i += 1
        term += math.log(lam) - math.log(i)
        terms.append(term)
        if term < terms[0] - 40:
            break
    top = max(terms)
    return min(0.0, (top + math.log(sum(math.exp(t - top) for t in terms))) / math.log(10))


def verify(wcs, sources, ra, dec, shape, radius=VERIFY_RADIUS, n_sources=VERIFY_SOURCES):
    """How well a WCS puts catalogue stars on the brightest detections.

    Returns a dict: matched (of the `tested` brightest detections, those within `radius` px of a catalogue star),
    tested, in_frame (catalogue stars inside the frame), expected (matches chance alone would give for that
    catalogue density) and log10_chance (log10 of the chance of `matched` or more)."""
    out = dict(matched=0, tested=0, in_frame=0, expected=None, log10_chance=0.0)
    if wcs is None or len(sources) == 0 or len(ra) == 0:
        return out
    h, w = shape
    px = _project(wcs, ra, dec)
    inside = (px[:, 0] > -radius) & (px[:, 0] < w + radius) & (px[:, 1] > -radius) & (px[:, 1] < h + radius)
    out['in_frame'] = int(inside.sum())
    if not out['in_frame']:
        return out
    src = sources[:n_sources]
    d, _ = cKDTree(px[inside]).query(src)
    out['tested'] = len(src)
    out['matched'] = int(np.sum(d < radius))
    expected = len(src) * out['in_frame'] * math.pi * radius ** 2 / float(w * h)
    out['expected'] = expected
    out['log10_chance'] = _log10_poisson_tail(out['matched'], expected)
    return out


def verified(v, trail=None, min_matches=VERIFY_MIN_MATCHES, max_log10_chance=VERIFY_MAX_LOG10_CHANCE):
    """The acceptance test for a solution, from verify()'s numbers and, if known, the frame's trail_ellipticity."""
    if not v or v.get('expected') is None:
        return False
    if v['log10_chance'] > max_log10_chance:
        return False
    possible = max(1, min(v['tested'], v['in_frame']))
    if v['matched'] >= SPARSE_MIN_MATCHES and v['matched'] >= SPARSE_MIN_COMPLETENESS * possible:
        return True             # a sparse field, nearly every possible match made
    round_sources = trail is not None and trail < ROUND_BELOW
    fraction = VERIFY_MIN_FRACTION_ROUND if round_sources else VERIFY_MIN_FRACTION
    return v['matched'] >= min_matches and v['matched'] >= fraction * possible


def solved_scale(wcs):
    """Mean pixel scale of a WCS, degrees per pixel."""
    return float(np.mean(proj_plane_pixel_scales(wcs)))


def geometry(wcs):
    """(mean pixel scale in degrees, anisotropy): the singular values of the pixel-to-sky matrix, which a skew or
    unequal axes would split."""
    sv = np.linalg.svd(wcs.pixel_scale_matrix, compute_uv=False)
    return float(math.sqrt(sv[0] * sv[1])), float(sv[0] / sv[1] - 1.0)


def _scale_ok(wcs, scale_deg, tolerance=SCALE_TOLERANCE, max_anisotropy=MAX_ANISOTROPY):
    scale, anisotropy = geometry(wcs)
    return abs(scale - scale_deg) / scale_deg <= tolerance and anisotropy <= max_anisotropy


def solve(sources, ra, dec, centre, scale_deg, shape, inner=None, prior=None, grid=True, trail=None):
    """Plate-solve with the plate scale known.

    sources   detections (x, y), 0-based pixels, brightest first
    ra, dec   catalogue stars in a box CATALOGUE_BOX fields across about `centre`, brightest first
    centre    (ra, dec) of the header pointing, degrees
    scale_deg expected plate scale, degrees per pixel
    shape     image shape (rows, columns)
    inner     boolean mask of the catalogue stars inside the INNER_BOX box (all if None)
    prior     a WCS to try first (e.g. the frame's inherited one): refined and verified like any other
    grid      allow the rotation-grid fallback
    trail     the frame's trail_ellipticity, if measured: round sources may verify on a smaller share of matches

    Returns a dict: success, wcs (or None), method ('prior', 'pairs', 'grid'), matches, tested, verification
    (verify()'s dict for the accepted or best rejected solution), candidates (how many were refined).
    """
    out = dict(success=False, wcs=None, method=None, matches=0, tested=0, verification=None, candidates=0,
               error=None)
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    if len(sources) < 5:
        out['error'] = 'too few sources ({})'.format(len(sources))
        return out
    if len(ra) < 5:
        out['error'] = 'too few catalogue stars ({})'.format(len(ra))
        return out
    inner = np.ones(len(ra), bool) if inner is None else np.asarray(inner, bool)
    ra0, dec0 = centre
    std = np.column_stack(standard_coordinates(ra, dec, ra0, dec0))
    best = dict(v=None)

    def try_wcs(wcs0, method):
        out['candidates'] += 1
        wcs, _ = refine(sources, ra, dec, wcs0)
        if wcs is None or not _scale_ok(wcs, scale_deg):
            return False
        v = verify(wcs, sources, ra, dec, shape)
        if best['v'] is None or v['matched'] > best['v']['matched']:
            best['v'] = v
        if verified(v, trail):
            out.update(success=True, wcs=wcs, method=method, matches=v['matched'], tested=v['tested'],
                       verification=v)
            return True
        return False

    def try_rotations(rotations, method, n_sources, n_catalogue, bin_px, keep=3, catalogue=None):
        cat = std if catalogue is None else catalogue
        scored = []
        for theta, parity in rotations:
            catpx = cat[:n_catalogue] @ _linear_map(theta, parity, scale_deg).T
            t, votes = offset_vote(sources[:n_sources], catpx, bin_px)
            if t is not None:
                scored.append((votes, theta, parity, t))
        scored.sort(key=lambda s: -s[0])
        return scored[:keep]

    if prior is not None:
        try:
            if _scale_ok(prior, scale_deg) and try_wcs(_recentred(prior, sources, ra, dec), 'prior'):
                return out
        except Exception:
            pass

    # pair votes: rotation candidates from the inner box (best overlap with the image) and the whole box
    rotations = []
    for cat, n_cat in ((std[inner], PAIR_CATALOGUE_INNER), (std, PAIR_CATALOGUE_OUTER)):
        for parity in (1, -1):
            catpx = cat[:n_cat] / scale_deg * np.array([parity, 1.0])
            for angle in rotation_votes(sources[:PAIR_SOURCES], catpx):
                rotations += [(angle, parity), ((angle + 180.0) % 360.0, parity)]
    for votes, theta, parity, t in try_rotations(rotations, 'pairs', OFFSET_SOURCES, OFFSET_CATALOGUE, OFFSET_BIN):
        if votes >= 4 and try_wcs(linear_wcs(ra0, dec0, theta, parity, scale_deg, t), 'pairs'):
            return out

    if grid:
        # a coarse rotation grid, scored on the inner box with bigger bins, then refined in 0.25-degree steps
        coarse = [(theta, parity) for parity in (1, -1) for theta in np.arange(0.0, 360.0, GRID_STEP)]
        top = try_rotations(coarse, 'grid', GRID_SOURCES, GRID_CATALOGUE, GRID_BIN, keep=4, catalogue=std[inner])
        fine = [((theta + d) % 360.0, parity) for _, theta, parity, _ in top for d in (-0.5, -0.25, 0.0, 0.25, 0.5)]
        for votes, theta, parity, t in try_rotations(fine, 'grid', OFFSET_SOURCES, OFFSET_CATALOGUE, OFFSET_BIN):
            if votes >= 4 and try_wcs(linear_wcs(ra0, dec0, theta, parity, scale_deg, t), 'grid'):
                return out

    out['verification'] = best['v']
    out['error'] = 'no verified solution'
    return out


def _orientation(wcs, shape):
    """(sky position of the frame centre, position angle in degrees of the frame's +y axis, scale in degrees per
    pixel) of a WCS."""
    h, w = shape
    x, y = (w - 1) / 2.0, (h - 1) / 2.0
    c, up = wcs.pixel_to_world([x, x], [y, y + 100.0])
    return c, float(c.position_angle(up).deg), solved_scale(wcs)


def agrees(wcs, neighbour, shape, max_shift=NEIGHBOUR_SHIFT, max_rotation=NEIGHBOUR_ROTATION,
           max_scale=NEIGHBOUR_SCALE):
    """Whether a solution sits where its neighbour's does: centre within max_shift px, the frame turned no more than
    max_rotation degrees, the same scale to max_scale."""
    c1, pa1, s1 = _orientation(wcs, shape)
    c2, pa2, s2 = _orientation(neighbour, shape)
    shift = c1.separation(c2).deg / s2
    turn = abs((pa1 - pa2 + 180.0) % 360.0 - 180.0)
    return shift <= max_shift and turn <= max_rotation and abs(s1 / s2 - 1.0) <= max_scale


def solve_from_neighbour(sources, ra, dec, neighbour, scale_deg, shape, trail=None):
    """Solve a frame from the solution of a frame taken close to it in time (astrometry.py's second pass).

    The neighbour's WCS is recentred by an offset vote within NEIGHBOUR_SHIFT px and refined (or, with too few
    matches to refit, kept as recentred). The result is accepted if it verifies on its own, or, when the sources are
    round, if NEIGHBOUR_MIN_MATCHES or more sit on catalogue stars at chance 10^NEIGHBOUR_MAX_LOG10_CHANCE or less and
    the solution agrees with the neighbour's. Trailed frames get no such allowance: a
    trail passes through every star's position at some instant. Returns a dict as solve() does (method
    'neighbour')."""
    out = dict(success=False, wcs=None, method=None, matches=0, tested=0, verification=None, candidates=0,
               error=None)
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    if len(sources) <= NEIGHBOUR_MIN_MATCHES or len(ra) <= NEIGHBOUR_MIN_MATCHES:
        out['error'] = 'too few sources or catalogue stars'
        return out
    if not _scale_ok(neighbour, scale_deg):
        out['error'] = "the neighbour's solution is not at the expected scale"
        return out
    out['candidates'] = 1
    start = _recentred(neighbour, sources, ra, dec, search=NEIGHBOUR_SHIFT)
    wcs, _ = refine(sources, ra, dec, start)
    if wcs is None:
        wcs = start             # too few matches to refit: the neighbour's solution, moved by the offset vote
    if not _scale_ok(wcs, scale_deg):
        out['error'] = "nothing matched near the neighbour's solution"
        return out
    v = verify(wcs, sources, ra, dec, shape)
    out['verification'] = v
    round_sources = trail is not None and trail < ROUND_BELOW
    if verified(v, trail) or (round_sources and v['matched'] >= NEIGHBOUR_MIN_MATCHES and
                              v['log10_chance'] <= NEIGHBOUR_MAX_LOG10_CHANCE and agrees(wcs, neighbour, shape)):
        out.update(success=True, wcs=wcs, method='neighbour', matches=v['matched'], tested=v['tested'])
    else:
        out['error'] = "no verified solution near the neighbour's"
    return out


def trail_ellipticity(image, sources, n_sources=30, half=12):
    """How far the brightest sources are stretched in one common direction: the size of the median ellipticity
    vector (e1, e2) of their second moments, in (2*half+1)-pixel cutouts of the cleaned image.

    The trails of a telescope that stopped tracking all point the same way and give 0.5-0.9. Round stars give
    under 0.1. Crowding stretches sources too, but in random directions that cancel in the median. None if nothing
    could be measured."""
    h, w = image.shape
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    e1, e2 = [], []
    for x, y in sources[:n_sources]:
        xi, yi = int(round(x)), int(round(y))
        if xi < half or yi < half or xi >= w - half or yi >= h - half:
            continue
        c = image[yi - half:yi + half + 1, xi - half:xi + half + 1].astype(float)
        c = c - np.median(c)
        c[c < 0] = 0
        tot = c.sum()
        if tot <= 0:
            continue
        mx, my = (c * xx).sum() / tot, (c * yy).sum() / tot
        ixx = (c * (xx - mx) ** 2).sum() / tot
        iyy = (c * (yy - my) ** 2).sum() / tot
        ixy = (c * (xx - mx) * (yy - my)).sum() / tot
        if ixx + iyy <= 0:
            continue
        e1.append((ixx - iyy) / (ixx + iyy))
        e2.append(2 * ixy / (ixx + iyy))
    if not e1:
        return None
    return float(math.hypot(np.median(e1), np.median(e2)))


def axis_ratio(ellipticity):
    """Major over minor axis for a second-moment ellipticity |e| = (a^2 - b^2) / (a^2 + b^2)."""
    return math.sqrt((1 + ellipticity) / (1 - ellipticity)) if ellipticity < 1 else float('inf')
