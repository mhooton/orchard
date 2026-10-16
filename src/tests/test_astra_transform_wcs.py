"""
Unit tests for the WCS handling in download/astra_transform.py.

Astra frames are rotated (Andor) or mirrored (SPIRIT) on download so that they match the MaxIm/ACP orientation of
the earlier data. The image moved, so the WCS has to follow it. Measured on 2026-10-16: the transform moved CRPIX
but left Astra's PCi_j matrix untouched, so every Astra light frame in the archive carried a WCS rotated by 90
degrees about CRPIX (Europa 2026-09-01: 30 arcsec out at the centre, 742 arcsec at a corner, against an
independent Gaia plate solve). Covers:
  - a sky position stays on the same star through the rotation and the mirror, for every matrix convention
    (Astra's PC without CDELT, astrometry.net's CD with either parity, astropy's PC with CDELT, CDELT/CROTA2)
  - the measured Europa frame: the corrected header gives the plate-solved orientation
  - PA is recomputed for the stored orientation, with the transfer script's definition
  - an all-zero CD block (SPIRIT frames as delivered) is removed so astrometry.py solves the frame

Run:  pytest src/tests/test_astra_transform_wcs.py -q
"""
import math
import os
import sys

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from download.astra_transform import (            # noqa: E402
    astra_transform_unified,
    drop_singular_wcs,
    linear_matrix,
    plate_pa,
)

NX, NY = 9, 6   # native NAXIS1, NAXIS2: small and non-square, so a swapped axis cannot pass by accident

# Europa 2026-09-01T01:03:29, as Astra wrote it (native orientation, before the transform)
EUROPA_PC = [[-9.3859668511035e-06, -9.5989180007196e-05], [9.6034430595804e-05, -9.4006699007457e-06]]
EUROPA_CRPIX_NATIVE = (2089 - 984.2046769341, 1015.5990367466)
EUROPA_CRVAL = (292.43961447225, 11.035372553426)


def base_header(nx=NX, ny=NY):
    h = fits.Header()
    h['IMAGETYP'] = 'Light Frame'
    h['OBSERVER'] = 'Astra'
    h['CTYPE1'] = 'RA---TAN'
    h['CTYPE2'] = 'DEC--TAN'
    h['CRVAL1'] = 150.0
    h['CRVAL2'] = -30.0
    h['CRPIX1'] = 3.3
    h['CRPIX2'] = 2.2
    return h


def header_pc(m):
    h = base_header()
    for (i, j), v in np.ndenumerate(np.asarray(m)):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    return h


def header_cd(m):
    h = base_header()
    for (i, j), v in np.ndenumerate(np.asarray(m)):
        h['CD{}_{}'.format(i + 1, j + 1)] = v
    return h


def header_pc_cdelt(m, cdelt=(2.0e-5, 3.0e-5)):
    h = base_header()
    pc = np.diag(1 / np.array(cdelt)) @ np.asarray(m)
    for (i, j), v in np.ndenumerate(pc):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    h['CDELT1'], h['CDELT2'] = cdelt
    return h


def header_crota(c1=-9.6e-5, c2=9.6e-5, rho=7.0):
    h = base_header()
    h['CDELT1'], h['CDELT2'], h['CROTA2'] = c1, c2, rho
    return h


UNFLIPPED = [[-9.6e-5, 1.2e-5], [-1.1e-5, -9.5e-5]]   # det > 0, as the SPECULOOS Andors are
FLIPPED = [[-9.6e-5, 1.2e-5], [1.1e-5, 9.5e-5]]       # det < 0, where the old CD formula was 180 deg out

def header_sip():
    """astrometry.net with tweak: RA---TAN-SIP, CD, forward A/B and inverse AP/BP polynomials."""
    h = header_cd(UNFLIPPED)
    h['CTYPE1'], h['CTYPE2'] = 'RA---TAN-SIP', 'DEC--TAN-SIP'
    h['A_ORDER'] = h['B_ORDER'] = 2
    h['A_2_0'], h['A_1_1'], h['A_0_2'] = 2e-3, -1e-3, 5e-4
    h['B_2_0'], h['B_1_1'], h['B_0_2'] = -7e-4, 1.5e-3, 3e-3
    h['AP_ORDER'] = h['BP_ORDER'] = 2
    h['AP_2_0'], h['AP_0_2'] = -2e-3, -5e-4
    h['BP_1_1'], h['BP_0_2'] = -1.5e-3, -3e-3
    return h


HEADERS = {
    'sip': header_sip,
    'astra_pc': lambda: header_pc(UNFLIPPED),
    'cd_unflipped': lambda: header_cd(UNFLIPPED),
    'cd_flipped': lambda: header_cd(FLIPPED),
    'pc_cdelt': lambda: header_pc_cdelt(UNFLIPPED),
    'pc_cdelt_flipped': lambda: header_pc_cdelt(FLIPPED),
    'crota': lambda: header_crota(),
}


def write_frame(tmp_path, header, name='frame'):
    data = np.arange(NX * NY, dtype=np.int32).reshape(NY, NX) + 1000   # every pixel distinct
    fits.PrimaryHDU(data=data, header=header).writeto(tmp_path / (name + '.fits'))
    return data


def sky_of_native_pixels(header):
    w = WCS(header)
    ys, xs = np.mgrid[0:NY, 0:NX]
    return w.pixel_to_world(xs.ravel(), ys.ravel())


@pytest.mark.parametrize('mode', ['rotate', 'mirror'])
@pytest.mark.parametrize('kind', sorted(HEADERS))
def test_each_pixel_keeps_its_sky_position(tmp_path, kind, mode):
    header = HEADERS[kind]()
    data = write_frame(tmp_path, header)
    before = sky_of_native_pixels(header)
    astra_transform_unified(['frame'], transform_mode=mode, n_repeats=1, dirname=str(tmp_path))
    with fits.open(tmp_path / 'frame.fits') as h:
        out, hdr = h[0].data, h[0].header
    w = WCS(hdr)
    worst = 0.0
    for k, value in enumerate(data.ravel()):
        y, x = np.argwhere(out == value)[0]
        worst = max(worst, w.pixel_to_world(x, y).separation(before[k]).arcsec)
    assert worst < 1e-6, f'{kind}/{mode}: a pixel moved {worst:.3g} arcsec on the sky'
    if kind == 'sip':
        # the inverse polynomials must follow too: world -> pixel lands back on the same pixel
        xs, ys = np.meshgrid(np.arange(out.shape[1]), np.arange(out.shape[0]))
        back = w.world_to_pixel(w.pixel_to_world(xs.ravel(), ys.ravel()))
        assert np.allclose(back[0], xs.ravel(), atol=1e-3) and np.allclose(back[1], ys.ravel(), atol=1e-3)
        # astropy inverts SIP iteratively, so check AP/BP directly: they were set to -A/-B term by term and the
        # same relabelling must keep that true
        n_checked = 0
        for inv, fwd in (('AP', 'A'), ('BP', 'B')):
            for key in [k for k in hdr if k.startswith(inv + '_') and k != inv + '_ORDER']:
                assert hdr[key] == pytest.approx(-hdr[fwd + key[len(inv):]])
                n_checked += 1
        assert n_checked == 4


def test_measured_europa_frame_gets_the_plate_solved_orientation(tmp_path):
    """Plate-solved against Gaia on 2026-10-16: +Y at PA 174.40, +X at PA 264.40, parity +1, 0.3474 arcsec/px."""
    h = fits.Header()
    h['IMAGETYP'] = 'Light Frame'
    h['CTYPE1'], h['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    h['CRVAL1'], h['CRVAL2'] = EUROPA_CRVAL
    h['CRPIX1'], h['CRPIX2'] = EUROPA_CRPIX_NATIVE
    for (i, j), v in np.ndenumerate(np.array(EUROPA_PC)):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    h['PA'] = (84.41224818505546, '[deg, 0-360 CCW] Position angle of plate')
    fits.PrimaryHDU(data=np.zeros((2048, 2088), dtype=np.int16), header=h).writeto(tmp_path / 'europa.fits')
    astra_transform_unified(['europa'], transform_mode='rotate', n_repeats=1, dirname=str(tmp_path))
    hdr = fits.getheader(tmp_path / 'europa.fits')
    m = WCS(hdr).pixel_scale_matrix
    pa_plus_y = math.degrees(math.atan2(m[0, 1], m[1, 1])) % 360
    pa_plus_x = math.degrees(math.atan2(m[0, 0], m[1, 0])) % 360
    assert (hdr['NAXIS1'], hdr['NAXIS2']) == (2048, 2088)
    assert abs(pa_plus_y - 174.40) < 0.05
    assert abs(pa_plus_x - 264.40) < 0.05
    assert np.linalg.det(m) > 0
    assert abs(math.sqrt(abs(np.linalg.det(m))) * 3600 - 0.3474) < 0.0005
    # The stored PA describes the stored orientation (the PA of -Y, as the transfer script defines it; that
    # formula averages the two axes, so it differs from the -Y column by the matrix's slight skew, 0.006 deg)
    assert abs(hdr['PA'] - (pa_plus_y + 180) % 360) < 0.01
    assert abs(hdr['PA'] - (84.41224818505546 - 90) % 360) < 1e-6


def test_plate_pa_matches_the_transfer_script_values():
    """headerfix.py gave 359.297 (MaxIm/astrometry.net CD) and 84.412 (Astra PC) for these two frames."""
    cd_2019 = np.array([[-9.63886583146e-05, 1.18259776625e-06], [-1.18259776625e-06, -9.63886583146e-05]])
    assert abs(plate_pa(cd_2019) - 359.2970701895605) < 1e-9
    assert abs(plate_pa(np.array(EUROPA_PC)) - 84.41224818505546) < 1e-9


@pytest.mark.parametrize('theta', [0.0, 5.6, 93.0, 181.5, 266.0])
@pytest.mark.parametrize('parity', [1, -1])
def test_plate_pa_is_the_position_angle_of_minus_y(theta, parity):
    """Exactly so for a skew-free matrix of either parity; for a skewed one it is the mean of the two axes."""
    c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
    m = 9.6e-5 * (np.array([[c, -s], [s, c]]) if parity > 0 else np.array([[-c, s], [s, c]]))
    minus_y = -m[:, 1]   # IWC (east, north) of a step towards row 1
    expected = math.degrees(math.atan2(minus_y[0], minus_y[1])) % 360
    assert abs((plate_pa(m) - expected + 180) % 360 - 180) < 1e-9


def test_singular_cd_block_is_removed():
    h = header_cd([[0.0, 0.0], [0.0, 0.0]])
    h['PA'] = 180.0
    assert drop_singular_wcs(h)
    assert 'CTYPE1' not in h and 'CD1_1' not in h and 'PA' not in h and 'CRVAL1' not in h
    assert h['IMAGETYP'] == 'Light Frame'


def test_valid_wcs_is_not_removed():
    h = header_pc(UNFLIPPED)
    assert not drop_singular_wcs(h)
    assert linear_matrix(h)[1] == 'PC'


def test_pc_takes_precedence_over_cd_as_in_wcslib():
    h = header_cd([[0.0, 0.0], [0.0, 0.0]])
    for (i, j), v in np.ndenumerate(np.array(UNFLIPPED)):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    m, convention = linear_matrix(h)
    assert convention == 'PC'
    assert np.allclose(m, WCS(h).pixel_scale_matrix)


# --- repair_astra_wcs: frames transformed by the old code ------------------------------------------------------

from download.repair_astra_wcs import diagnose, repair        # noqa: E402


def old_transform_output():
    """What the old astra_transform wrote for the Europa frame: CRPIX moved, Astra's PC left as it was."""
    h = fits.Header()
    h['NAXIS'] = 2
    h['NAXIS1'], h['NAXIS2'] = 2048, 2088
    h['IMAGETYP'] = 'Light Frame'
    h['CTYPE1'], h['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    h['CRVAL1'], h['CRVAL2'] = EUROPA_CRVAL
    h['CRPIX1'], h['CRPIX2'] = 1015.5990367466, 984.2046769341
    for (i, j), v in np.ndenumerate(np.array(EUROPA_PC)):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    h['PA'] = (84.41224818505546, '[deg, 0-360 CCW] Position angle of plate')
    h['ASTRAROT'] = (90, 'Anticlockwise rotation applied to astra images')
    h.add_history('WCS updated for 90 deg  anticlockwise rotation')
    return h


def test_repair_finds_and_fixes_the_old_rotation(tmp_path):
    h = old_transform_output()
    found = dict(diagnose(h))
    assert set(found) == {'rot_pc', 'stale_pa'}
    repair(h, list(found.items()))
    # Same WCS as the fixed transform produces from the native header
    native = fits.Header()
    native['IMAGETYP'] = 'Light Frame'
    native['CTYPE1'], native['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    native['CRVAL1'], native['CRVAL2'] = EUROPA_CRVAL
    native['CRPIX1'], native['CRPIX2'] = EUROPA_CRPIX_NATIVE
    for (i, j), v in np.ndenumerate(np.array(EUROPA_PC)):
        native['PC{}_{}'.format(i + 1, j + 1)] = v
    native['PA'] = 84.41224818505546
    fits.PrimaryHDU(data=np.zeros((2048, 2088), dtype=np.int16), header=native).writeto(tmp_path / 'n.fits')
    astra_transform_unified(['n'], transform_mode='rotate', n_repeats=1, dirname=str(tmp_path))
    good = fits.getheader(tmp_path / 'n.fits')
    xs, ys = np.array([0, 2047, 0, 2047, 1023]), np.array([0, 0, 2087, 2087, 1043])
    sep = WCS(h).pixel_to_world(xs, ys).separation(WCS(good).pixel_to_world(xs, ys)).arcsec
    assert sep.max() < 1e-6
    assert abs(h['PA'] - good['PA']) < 1e-9
    assert diagnose(h) == []          # idempotent


def test_repair_leaves_an_orchard_solution_alone_but_fixes_its_pa():
    h = old_transform_output()
    # orchard's plate solution, made on the stored image and written by astropy (PC with CDELT = 1); the PA is
    # still the transfer script's, for the native orientation
    stored = np.array(EUROPA_PC) @ np.array([[0, -1], [1, 0]])
    for (i, j), v in np.ndenumerate(stored):
        h['PC{}_{}'.format(i + 1, j + 1)] = v
    h['CDELT1'], h['CDELT2'] = 1.0, 1.0
    found = dict(diagnose(h))
    assert set(found) == {'stale_pa'}
    before = WCS(h).pixel_scale_matrix.copy()
    repair(h, list(found.items()))
    assert np.allclose(WCS(h).pixel_scale_matrix, before)


def test_repair_removes_a_singular_wcs():
    h = old_transform_output()
    for k in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'):
        del h[k]
    for k in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
        h[k] = 0.0
    h['PA'] = 180.0
    found = diagnose(h)
    assert [k for k, _ in found] == ['singular']
    repair(h, found)
    assert 'CTYPE1' not in h and 'PA' not in h and h['ASTRAWCS'] == 2


def test_new_transform_output_is_not_repaired(tmp_path):
    header = header_pc(UNFLIPPED)
    header['PA'] = 10.0
    write_frame(tmp_path, header)
    astra_transform_unified(['frame'], transform_mode='rotate', n_repeats=1, dirname=str(tmp_path))
    assert diagnose(fits.getheader(tmp_path / 'frame.fits')) == []


def test_repair_finds_a_wcs_left_entirely_native(tmp_path):
    """Callisto 2024: the image rotated, ASTRAROT set, no transform HISTORY, CRPIX and CD as ESO delivered them."""
    native = header_cd(UNFLIPPED)
    native['PA'] = (plate_pa(np.array(UNFLIPPED)), '[deg, 0-360 CCW] Position angle of plate')
    data = write_frame(tmp_path, native.copy(), name='good')
    astra_transform_unified(['good'], transform_mode='rotate', n_repeats=1, dirname=str(tmp_path))
    good = fits.getheader(tmp_path / 'good.fits')
    old = native.copy()
    old['NAXIS'] = 2
    old['NAXIS1'], old['NAXIS2'] = NY, NX          # the data were rotated, the WCS cards were not
    old['ASTRAROT'] = (90, 'Anticlockwise rotation applied to astra images')
    found = dict(diagnose(old))
    assert set(found) == {'unrot', 'stale_pa'}
    repair(old, list(found.items()))
    ys, xs = np.mgrid[0:NX, 0:NY]
    sep = WCS(old).pixel_to_world(xs.ravel(), ys.ravel()).separation(
        WCS(good).pixel_to_world(xs.ravel(), ys.ravel())).arcsec
    assert sep.max() < 1e-6
    assert abs(old['PA'] - good['PA']) < 1e-9
    assert diagnose(old) == []


def test_a_native_shaped_image_with_astrarot_is_reported_not_fixed():
    h = old_transform_output()
    h['NAXIS1'], h['NAXIS2'] = 2088, 2048
    found = diagnose(h)
    assert [k for k, _ in found] == ['image_not_rotated']
    before = h.copy()
    repair(h, found)
    assert h == before


def test_a_native_matrix_under_transform_history_is_reported_not_fixed():
    """Ganymede 2025-07-28: two transform HISTORY lines, CD matrix and PA both still native."""
    h = old_transform_output()
    for k in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'):
        del h[k]
    for (i, j), v in np.ndenumerate(np.array(EUROPA_PC)):
        h['CD{}_{}'.format(i + 1, j + 1)] = v
    h.add_history('WCS updated for 90 deg counterclockwise rotation')
    found = diagnose(h)
    assert [k for k, _ in found] == ['ambiguous']
    before = h.copy()
    repair(h, found)
    assert h == before


def mirrored_spirit(native_m, history=None):
    h = base_header()
    h['NAXIS'] = 2
    h['NAXIS1'], h['NAXIS2'] = NX, NY
    for (i, j), v in np.ndenumerate(np.asarray(native_m)):
        h['CD{}_{}'.format(i + 1, j + 1)] = v
    h['PA'] = (plate_pa(np.asarray(native_m)), '[deg, 0-360 CCW] Position angle of plate')
    h['ASTRAMIR'] = ('Y', 'Y axis flip applied to astra images')
    if history:
        h.add_history(history)
    return h


def test_mirror_keeps_pa_so_a_correctly_mirrored_frame_is_not_flagged(tmp_path):
    native = header_cd(FLIPPED)          # SPIRIT's ESO parity is -1
    native['PA'] = (plate_pa(np.array(FLIPPED)), '[deg, 0-360 CCW] Position angle of plate')
    write_frame(tmp_path, native)
    astra_transform_unified(['frame'], transform_mode='mirror', n_repeats=1, dirname=str(tmp_path))
    out = fits.getheader(tmp_path / 'frame.fits')
    assert abs(out['PA'] - native['PA']) < 1e-9      # the angle of -Y does not change under a left-right mirror
    del out['ASTRAWCS']                               # as if written by the old code, which got CD mirrors right
    assert [k for k, _ in diagnose(out)] == []


def test_unmirrored_spirit_matrix_without_history_is_fixed(tmp_path):
    h = mirrored_spirit(FLIPPED)
    found = dict(diagnose(h))
    assert 'unmir' in found
    good = header_cd(FLIPPED)
    good['PA'] = h['PA']
    write_frame(tmp_path, good, name='g')
    astra_transform_unified(['g'], transform_mode='mirror', n_repeats=1, dirname=str(tmp_path))
    g = fits.getheader(tmp_path / 'g.fits')
    repair(h, list(found.items()))
    ys, xs = np.mgrid[0:NY, 0:NX]
    sep = WCS(h).pixel_to_world(xs.ravel(), ys.ravel()).separation(WCS(g).pixel_to_world(xs.ravel(), ys.ravel())).arcsec
    assert sep.max() < 1e-6


def test_native_parity_under_mirror_history_is_ambiguous():
    h = mirrored_spirit(FLIPPED, history='WCS updated for horizontal mirroring')
    h['CDELT1'] = h['CDELT2'] = 1.0       # keep it out of the mir_pc branch
    assert [k for k, _ in diagnose(h)] == ['ambiguous']


# --- writing the repair to disk: header bytes only, backed up, reversible ----------------------------------------

import tarfile                                                   # noqa: E402
from download.repair_astra_wcs import apply_to_file, read_header_bytes, restore   # noqa: E402


def raw_like_frame(path, filler_cards=0):
    """A small uint16 frame stored like the raw archive (BITPIX 16, BZERO 32768) with the old-transform WCS."""
    h = old_transform_output()
    for k in ('NAXIS', 'NAXIS1', 'NAXIS2'):
        del h[k]
    for i in range(filler_cards):
        h['FILL{:04d}'.format(i)] = (i, 'padding to fill the header block')
    data = (np.arange(30 * 20, dtype=np.uint16).reshape(30, 20) * 97) % 65535
    fits.PrimaryHDU(data=data, header=h).writeto(path)
    os.utime(path, (1_600_000_000, 1_600_000_000))
    return open(path, 'rb').read()


def data_bytes(path):
    raw, _ = read_header_bytes(path)
    return open(path, 'rb').read()[len(raw):]


@pytest.mark.parametrize('filler', [0, 'overflow'])
def test_apply_writes_header_only_and_restore_undoes_it(tmp_path, filler):
    path = str(tmp_path / 'frame.fits')
    if filler == 'overflow':
        # fill the header so that the two cards the repair adds need one more 2880-byte block
        raw = raw_like_frame(path)
        ncards = len(read_header_bytes(path)[1])
        n_fill = (36 - (ncards + 1) % 36) % 36
        os.remove(path)
        raw = raw_like_frame(path, filler_cards=n_fill)
    else:
        raw = raw_like_frame(path)
    old_header_len = len(read_header_bytes(path)[0])
    old_data = data_bytes(path)
    with tarfile.open(str(tmp_path / 'backup.tar'), 'a') as t:
        found = apply_to_file(path, t)
    assert {k for k, _ in found} == {'rot_pc', 'stale_pa'}
    new_raw, header = read_header_bytes(path)
    assert header['ASTRAWCS'] == 2
    assert data_bytes(path) == old_data                          # pixels untouched, byte for byte
    assert (len(new_raw) > old_header_len) == (filler == 'overflow')
    assert os.stat(path).st_mtime == 1_600_000_000                # timestamps kept
    with fits.open(path) as h:                                    # still a valid FITS file that astropy reads
        assert h[0].data.dtype.kind in 'ui' and h[0].data[3, 4] == (3 * 20 + 4) * 97 % 65535
    assert restore(str(tmp_path / 'backup.tar')) == 1
    assert open(path, 'rb').read() == raw                         # back to the original, byte for byte


def test_apply_leaves_reported_only_frames_alone(tmp_path):
    path = str(tmp_path / 'frame.fits')
    h = old_transform_output()
    for k in ('NAXIS', 'NAXIS1', 'NAXIS2'):
        del h[k]
    fits.PrimaryHDU(data=np.zeros((2048, 2088), dtype=np.uint16), header=h).writeto(path)   # native shape
    before = open(path, 'rb').read()
    with tarfile.open(str(tmp_path / 'b.tar'), 'a') as t:
        found = apply_to_file(path, t)
    assert [k for k, _ in found] == ['image_not_rotated']
    assert open(path, 'rb').read() == before
