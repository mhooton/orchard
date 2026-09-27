"""
Unit tests for the plate-solve failure that reported

    FAILED - ... Error: 'NoneType' object has no attribute 'cpdis1'

twirl.compute_wcs returns None when it cannot match an asterism, and that None
used to be handed straight to astropy, which asks _has_distortion(wcs) before
doing anything else and looks for `cpdis1` first. Covers:
  - the mechanism, so the regression is recognisable if it comes back
  - the guard and the wording of the failure it raises instead
  - the too-few-stars gate, which fires before twirl is called
  - the (RA, Dec) ordering of the catalogue query box
  - gaia_db_query's tmass path: empty-string columns, limit, proper motion

Run:  pytest src/tests/test_platesolve_wcs_guard.py -q
"""
import os
import sqlite3
import sys

import numpy as np
import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from astrom import pointer_wcs as pw              # noqa: E402
from astrom.pointer_wcs import (                  # noqa: E402
    TWIRL_ASTERISM,
    ImageStarMapping,
    gaia_db_query,
    image_field_of_view,
)

# Eight image sources and eight catalogue stars: enough that the asterism gate
# passes and the twirl call is reached.
STARS = np.array([[100.0, 100.0], [200.0, 150.0], [300.0, 400.0], [450.0, 120.0],
                  [500.0, 600.0], [640.0, 210.0], [700.0, 800.0], [820.0, 330.0]])
GAIAS = np.array([[32.30, -3.10], [32.31, -3.11], [32.32, -3.09], [32.33, -3.12],
                  [32.34, -3.08], [32.35, -3.13], [32.36, -3.07], [32.37, -3.14]])


def test_none_wcs_in_astropy_is_the_cpdis1_error():
    """The mechanism behind the reported message, pinned down."""
    from astropy.coordinates import SkyCoord
    with pytest.raises(AttributeError) as e:
        SkyCoord(GAIAS, unit="deg").to_pixel(None)
    assert "cpdis1" in str(e.value)


def test_unmatched_asterism_reports_the_inputs(monkeypatch):
    monkeypatch.setattr(pw.twirl, "compute_wcs", lambda *a, **k: None)
    with pytest.raises(ValueError) as e:
        ImageStarMapping.from_gaia_coordinates(STARS, GAIAS)
    message = str(e.value)
    assert "cpdis1" not in message
    assert "no %d-star asterism" % TWIRL_ASTERISM in message
    assert "8 image sources" in message and "8 catalogue stars" in message


def test_twirl_raising_is_reported_as_a_failed_fit(monkeypatch):
    """scipy's "`x0` is infeasible." comes out of twirl's own least-squares fit."""
    def boom(*a, **k):
        raise ValueError("`x0` is infeasible.")
    monkeypatch.setattr(pw.twirl, "compute_wcs", boom)
    with pytest.raises(ValueError) as e:
        ImageStarMapping.from_gaia_coordinates(STARS, GAIAS)
    message = str(e.value)
    assert "could not fit a WCS" in message
    assert "x0` is infeasible" in message
    assert "8 image sources" in message


@pytest.mark.parametrize("n_stars,n_gaia", [(4, 8), (8, 4), (4, 4), (1, 20)])
def test_too_few_stars_never_reaches_twirl(monkeypatch, n_stars, n_gaia):
    called = []
    monkeypatch.setattr(pw.twirl, "compute_wcs",
                        lambda *a, **k: called.append(1))
    stars = np.tile(STARS, (10, 1))[:n_stars]
    gaias = np.tile(GAIAS, (10, 1))[:n_gaia]
    with pytest.raises(ValueError) as e:
        ImageStarMapping.from_gaia_coordinates(stars, gaias)
    assert "too few stars" in str(e.value)
    assert "more than %d of each" % TWIRL_ASTERISM in str(e.value)
    assert called == []


def test_more_than_the_asterism_is_allowed_through(monkeypatch):
    """Five and five is the smallest pair that twirl is still asked to try."""
    called = []

    def record(stars, gaias, **kwargs):
        called.append((len(stars), len(gaias), kwargs.get("asterism")))
        return None
    monkeypatch.setattr(pw.twirl, "compute_wcs", record)
    with pytest.raises(ValueError):
        ImageStarMapping.from_gaia_coordinates(STARS[:5], GAIAS[:5])
    assert called == [(5, 5, TWIRL_ASTERISM)]


def test_field_of_view_maps_ra_to_columns():
    """gaia_db_query reads fov as (RA extent, Dec extent)."""
    plate_scale = 0.3094 / 3600.0
    # SPIRIT: 1280 rows x 1024 columns
    fov = image_field_of_view((1280, 1024), plate_scale, np.radians(-3.10272))
    ra_extent, dec_extent = fov
    assert dec_extent == pytest.approx(1280 * plate_scale)
    assert ra_extent == pytest.approx(1024 * plate_scale / np.cos(np.radians(-3.10272)))
    # the long axis of the detector must be the long axis on the sky
    assert dec_extent > ra_extent
    # a square detector cannot tell the two apart
    square = image_field_of_view((2048, 2048), plate_scale, 0.0)
    assert square[0] == pytest.approx(square[1])


def make_db(tmp_path, rows):
    """A one-shard stand-in for the local Gaia/2MASS database."""
    path = os.path.join(str(tmp_path), "gaia.db")
    conn = sqlite3.connect(path)
    conn.execute('CREATE TABLE `-4_-3` (ra, dec, pmra, pmdec, phot_g_mean_mag, j_m)')
    conn.executemany('INSERT INTO `-4_-3` VALUES (?,?,?,?,?,?)', rows)
    conn.commit()
    conn.close()
    return path


def test_tmass_query_is_ordered_limited_and_epoch_corrected(tmp_path):
    import pandas as pd
    # j_m ascending is 11.0, 12.0, 13.0, 14.0; the 14.0 row has empty-string
    # proper motions, as some shards store a missing value.
    rows = [
        (32.400, -3.500, 0.0, 0.0, 15.0, 13.0),
        (32.401, -3.501, 3600000.0, 0.0, 15.0, 12.0),   # 1 deg/yr in RA
        (32.402, -3.502, 0.0, 3600000.0, 15.0, 11.0),   # 1 deg/yr in Dec
        (32.403, -3.503, '', '', 15.0, 14.0),
    ]
    db = make_db(tmp_path, rows)
    fov = np.array([1.0, 1.0])
    dateobs = pd.to_datetime("2026-01-01T00:00:00")     # 10.0 yr after J2016.0
    out = gaia_db_query(center=(32.4, -3.5), fov=fov, tmass=True,
                        dateobs=dateobs, limit=3, db_path=db)

    assert out.dtype == np.float64
    # the 2MASS path returns the whole box; the caller does the cutting, so the
    # plate-solve log can still report how rich the field is
    assert len(out) == 4
    # brightest in J first
    assert out[0][1] == pytest.approx(-3.502 + 10.0, abs=2e-3)   # pmdec star
    assert out[1][0] == pytest.approx(32.401 + 10.0 / np.cos(np.radians(-3.501)),
                                     abs=2e-3)                   # pmra star
    assert out[2] == pytest.approx([32.400, -3.500])              # no motion
    assert out[3] == pytest.approx([32.403, -3.503])              # empty pm -> zero


def test_limit_applies_to_the_magnitude_ordered_path(tmp_path):
    rows = [(32.400 + i / 1000.0, -3.500, 0.0, 0.0, 15.0 + i, 13.0 + i) for i in range(6)]
    db = make_db(tmp_path, rows)
    out = gaia_db_query(center=(32.4, -3.5), fov=np.array([1.0, 1.0]), tmass=False,
                        dateobs=None, limit=3, db_path=db)
    assert len(out) == 3
    assert out[0] == pytest.approx([32.400, -3.500])


def test_tmass_query_drops_rows_without_usable_coordinates(tmp_path):
    rows = [
        (32.400, -3.500, 0.0, 0.0, 15.0, 13.0),
        ('', -3.501, 0.0, 0.0, 15.0, 12.0),
        (32.402, '', 0.0, 0.0, 15.0, 11.0),
        (32.403, -3.503, 0.0, 0.0, 15.0, ''),
    ]
    db = make_db(tmp_path, rows)
    out = gaia_db_query(center=(32.4, -3.5), fov=np.array([1.0, 1.0]), tmass=True,
                        dateobs=None, limit=10, db_path=db)
    assert out.dtype == np.float64
    assert len(out) == 1
    assert out[0] == pytest.approx([32.400, -3.500])


def test_empty_box_returns_no_stars(tmp_path):
    db = make_db(tmp_path, [(40.0, -3.5, 0.0, 0.0, 15.0, 13.0)])
    out = gaia_db_query(center=(32.4, -3.5), fov=np.array([0.01, 0.01]), tmass=True,
                        dateobs=None, limit=10, db_path=db)
    assert out.shape == (0, 2)
