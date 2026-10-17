"""
Unit tests for the stale PC matrix that mis-centred T9's forced apertures.

Frames from Astra 0.3.2 onwards (Callisto from 2025-07, Ganymede from 2025-11,
Europa and Io from 2026-01) carry PCi_j/CDELTi next to the CDi_j that wcsfit
writes. wcslib takes PC over CD, so imcore_list, imstack and astropy used the
stale matrix, not the fitted one. Covers:
  - the precedence itself, so the mechanism is recognisable if it comes back
  - drop_stale_pc_matrix: removes PC/CDELT only when CD is there to replace it
  - headers without PC, or without CD, are left alone
  - a 16-bit frame keeps its BZERO and its data bytes
  - casu_solve drops the PC straight after a successful wcsfit

Run:  pytest src/tests/test_stale_pc_matrix.py -q
"""
import os
import sys
import warnings

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

warnings.filterwarnings('ignore', category=FITSFixedWarning)

# wcsfit's CD from a Ganymede 20260223 frame, and a PC 0.8% larger in scale and
# 0.3 deg rotated from it: the size of the stale matrices found on those nights.
CD = {'CD1_1': -9.52774491690303e-05, 'CD1_2': -1.53142482517418e-05,
      'CD2_1': 1.53337328670726e-05, 'CD2_2': -9.52715702818309e-05}
PC = {'PC1_1': -9.5947e-05, 'PC1_2': -1.6043e-05,
      'PC2_1': 1.5893e-05, 'PC2_2': -9.6034e-05}
BASE = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRVAL1': 163.623177091396,
        'CRVAL2': -54.5083842720856, 'CRPIX1': 1008.5248785684, 'CRPIX2': 837.833781521929}
CORNERS = np.array([[1.0, 1.0], [2046.0, 1.0], [1.0, 2044.0], [2046.0, 2044.0]])


def header(pc=True, cd=True):
    h = fits.Header()
    h['NAXIS'] = 2
    h['NAXIS1'] = 2046
    h['NAXIS2'] = 2044
    h.update(BASE)
    if pc:
        h.update(PC)
        h['CDELT1'] = 1.0
        h['CDELT2'] = 1.0
    if cd:
        h.update(CD)
    return h


def sky(h):
    return WCS(h).all_pix2world(CORNERS, 1)


def write(tmp_path, h, data=None, name='frame.fits'):
    path = str(tmp_path / name)
    if data is None:
        data = np.zeros((h['NAXIS2'], h['NAXIS1']), dtype=np.float32)
    fits.PrimaryHDU(data=data, header=h).writeto(path)
    return path


def test_wcslib_takes_pc_over_cd():
    """The mechanism: with both present, the header means the PC solution."""
    both, pc_only, cd_only = sky(header()), sky(header(cd=False)), sky(header(pc=False))
    assert np.allclose(both, pc_only, atol=1e-10, rtol=0)
    # ... which puts a corner star several pixels from where the fitted CD has it
    gap = WCS(header(pc=False)).all_world2pix(both, 1) - CORNERS
    assert np.max(np.hypot(gap[:, 0], gap[:, 1])) > 5


def test_drop_leaves_the_fitted_cd(tmp_path):
    from photometry.wcs_fitting import drop_stale_pc_matrix
    path = write(tmp_path, header())
    written = fits.getheader(path)
    assert drop_stale_pc_matrix(path) is True
    h = fits.getheader(path)
    for key in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2'):
        assert key not in h
    for key in list(CD) + list(BASE):
        assert h[key] == written[key]
    assert np.allclose(sky(h), sky(header(pc=False)), atol=1e-10, rtol=0)
    assert any('PC/CDELT removed' in str(c) for c in h['HISTORY'])


@pytest.mark.parametrize('pc,cd', [(False, True), (True, False)])
def test_headers_without_both_are_untouched(tmp_path, pc, cd):
    """No PC: nothing to do. No CD (wcsfit never wrote one): PC is all there is."""
    from photometry.wcs_fitting import drop_stale_pc_matrix
    path = write(tmp_path, header(pc=pc, cd=cd))
    before = fits.getheader(path).tostring()
    assert drop_stale_pc_matrix(path) is False
    assert fits.getheader(path).tostring() == before


def test_sixteen_bit_frame_keeps_bzero_and_data(tmp_path):
    from photometry.wcs_fitting import drop_stale_pc_matrix
    h = header()
    raw = np.arange(2044 * 2046, dtype=np.int64).reshape(2044, 2046) % 60000 - 1000
    hdu = fits.PrimaryHDU(data=(raw - 31000).astype(np.int16), header=h)
    hdu.header['BZERO'] = 31000
    hdu.header['BSCALE'] = 1
    path = str(tmp_path / 'frame16.fits')
    hdu.writeto(path)
    with fits.open(path, do_not_scale_image_data=True) as hdul:
        stored = hdul[0].data.copy()
    assert drop_stale_pc_matrix(path) is True
    with fits.open(path, do_not_scale_image_data=True) as hdul:
        assert hdul[0].header['BITPIX'] == 16
        assert hdul[0].header['BZERO'] == 31000
        assert np.array_equal(hdul[0].data, stored)
    assert np.array_equal(fits.getdata(path), raw)


def test_casu_solve_drops_pc_after_wcsfit(tmp_path, monkeypatch):
    """wcsfit writes CD and keeps the PC cards; casu_solve must not leave them."""
    from photometry import wcs_fitting
    path = write(tmp_path, header(cd=False))

    def fake_wcsfit(infile, incat, catsrc=None, verbose=None, site=None, wcsref=None):
        with fits.open(infile, mode='update') as hdul:
            hdul[0].header.update(CD)

    monkeypatch.setattr(wcs_fitting.casutools, 'imcore', lambda *a, **k: None)
    monkeypatch.setattr(wcs_fitting, 'make_localfits_table', lambda *a, **k: None)
    monkeypatch.setattr(wcs_fitting.casutools, 'wcsfit', fake_wcsfit)
    monkeypatch.setattr(wcs_fitting, 'wcsfit_qc', lambda *a, **k: None)
    assert wcs_fitting.casu_solve(path, ext='fits') == 'ok'
    h = fits.getheader(path)
    assert 'PC1_1' not in h and 'CDELT1' not in h
    assert np.allclose(sky(h), sky(header(pc=False)), atol=1e-10, rtol=0)
