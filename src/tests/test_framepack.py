"""
Unit tests for utils/framepack.py: 16-bit calibrated frames and lossless Rice
compression, done in place without risking the original.

Covers in particular:
  - the BZERO window: anchored at the faintest pixel, capped at 32,767, never
    0 or 32,768, and the wide-frame floor that only clips pixels CASU already
    treats as bad (Artemis 2021-04-01 and 2025-04-01 values)
  - rounding half away from zero, and NaN/inf stored as BLANK and read back
    as NaN
  - compress/decompress reproducing pixels, scaling and header cards exactly,
    for a 16-bit processed frame and for a raw frame (BZERO 32768)
  - in-place conversion keeping mode and mtime, leaving no hidden temporary
    files, and leaving the original untouched when anything fails

Run:  pytest src/tests/test_framepack.py -q
"""
import glob
import os
import sys
import warnings

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import framepack as fp                  # noqa: E402

warnings.simplefilter('ignore', category=fits.verify.VerifyWarning)


def _header():
    h = fits.Header()
    h['OBJECT'] = 'Sp1428+3310'
    h['EXPTIME'] = (30.0, 'exposure time')
    h['JD-OBS'] = 2460764.6123456789
    h['CRVAL1'] = 217.07208333333335
    h['CD1_1'] = -9.6123456789012e-05
    h['IMAGETYP'] = 'Light Frame'
    h['FILTER'] = 'I+z'
    h['ORIGFILE'] = 'SPECU2.' + 'x' * 80        # a long string, written with CONTINUE cards
    h['ASTSOLVE'] = True
    h['COMMENT'] = 'a comment card'
    h['HISTORY'] = 'an earlier history card'
    h['CHECKSUM'] = 'stale0000000000'            # stale, as copied from raw frames
    h['DATASUM'] = '12345'
    return h


def _float_frame(path, shape=(96, 128), seed=1, extra=None):
    rng = np.random.default_rng(seed)
    data = (1000 + 30 * rng.standard_normal(shape)).astype(np.float32)
    data[40:43, 50:53] += 20000                 # a star
    data[10, 10] = -40.3                        # a negative pixel
    if extra:
        extra(data)
    fits.PrimaryHDU(data=data, header=_header()).writeto(path)
    return data


def _raw_frame(path, shape=(96, 128), seed=2):
    rng = np.random.default_rng(seed)
    data = np.clip(300 + 1000 + 30 * rng.standard_normal(shape), 0, 65535).astype(np.uint16)
    data[5, 5] = 65535
    fits.PrimaryHDU(data=data, header=_header()).writeto(path)   # astropy stores uint16 as BZERO 32768
    return data


# --------------------------------------------------------------------------- the window

@pytest.mark.parametrize('lo, hi, expected', [
    (-9.8, 64146.8, 32757),         # a typical Ganymede frame: window -10 .. 65,524
    (9.7, 65184.6, 32767),          # faintest pixel above zero: capped, window 0 .. 65,534
    (-23.6, 67149.4, 32743),        # Artemis 2025-04-01: midpoint rule clipped the whole sky here
    (-1016.2, 11556.0, 31750),      # Callisto SPIRIT 2026
    (-48671.4, 64301.8, 31535),     # Artemis 2021-04-01: spans 113k ADU, floor raised to -1,232
])
def test_choose_bzero(lo, hi, expected):
    assert fp.choose_bzero(lo, hi) == expected


def test_window_covers_the_frame_or_clips_only_bad_or_saturated_pixels():
    for lo, hi in [(-9.8, 64146.8), (-48671.4, 64301.8), (-23.6, 67149.4)]:
        b = fp.choose_bzero(lo, hi)
        assert b + fp.STORED_MIN <= max(np.floor(lo), fp.CASU_BAD_BELOW - 1)   # faint end kept, or already bad
        assert fp.MIN_BZERO <= b <= fp.MAX_BZERO


def test_bzero_is_never_zero_or_32768():
    assert fp.choose_bzero(-32767.0, 100.0) == -1
    assert fp.choose_bzero(1.0, 60000.0) != 32768
    assert fp.choose_bzero(5000.0, 9000.0) == 32767


def test_a_frame_that_cannot_be_stored_is_refused():
    with pytest.raises(fp.FramePackError):
        fp.choose_bzero(-70000.0, -60000.0)


# --------------------------------------------------------------------------- quantising

def test_rounding_is_half_away_from_zero():
    x = np.array([[2.5, -2.5, 0.49, -0.5, 1.5, 100.0]])
    stored, bzero, n_clipped, n_blank = fp.quantise(x)
    assert (stored.astype(np.int64) + bzero).tolist() == [[3, -3, 0, -1, 2, 100]]
    assert (n_clipped, n_blank) == (0, 0)


def test_non_finite_pixels_become_blank():
    x = np.array([[1.0, np.nan, np.inf, -np.inf, 5.0]])
    stored, bzero, n_clipped, n_blank = fp.quantise(x)
    assert n_blank == 3
    assert stored[0, 1:4].tolist() == [fp.BLANK] * 3


def test_pixels_beyond_the_window_are_clipped_and_counted():
    x = np.array([[0.0, 100.0, 70000.0]])
    stored, bzero, n_clipped, _ = fp.quantise(x)
    assert bzero == 32767 and n_clipped == 1
    assert int(stored[0, 2]) == fp.STORED_MAX


# --------------------------------------------------------------------------- files

def test_int16_file_keeps_pixels_within_half_an_adu_and_every_card(tmp_path):
    src, dst = str(tmp_path / 'proc.fits'), str(tmp_path / 'out.fits')
    data = _float_frame(src, extra=lambda d: d.__setitem__((0, 0), np.nan))
    bzero, n_clipped, n_blank = fp.to_int16(src, dst)
    fp.check_int16(src, dst, bzero, n_clipped, n_blank)
    h = fits.getheader(dst)
    assert (h['BITPIX'], h['BZERO'], h['BSCALE'], h['BLANK']) == (16, bzero, 1.0, fp.BLANK)
    out = fits.getdata(dst)
    assert out.dtype == np.float32                      # scaled on reading, so safe for arithmetic
    assert np.isnan(out[0, 0]) and n_blank == 1
    finite = np.isfinite(data)
    assert np.nanmax(np.abs(out[finite] - data[finite])) <= 0.5
    assert fp._cards(fits.getheader(src)) == fp._cards(h)
    assert h['ORIGFILE'] == 'SPECU2.' + 'x' * 80
    assert os.path.getsize(dst) < 0.6 * os.path.getsize(src)


@pytest.mark.parametrize('kind', ['processed', 'raw'])
def test_compress_and_decompress_are_exact(tmp_path, kind):
    src = str(tmp_path / 'frame.fits')
    if kind == 'processed':
        _float_frame(str(tmp_path / 'float.fits'))
        fp.to_int16(str(tmp_path / 'float.fits'), src)
    else:
        _raw_frame(src)
    fz, back = src + '.fz', str(tmp_path / 'back.fits')
    fp.compress(src, fz)
    fp.check_lossless(src, fz)
    with fits.open(fz, disable_image_compression=True) as h:
        z = h[1].header
        assert (z['ZCMPTYPE'], z['ZBITPIX'], 'ZQUANTIZ' in z) == ('RICE_1', 16, False)
    assert os.path.getsize(fz) < os.path.getsize(src)
    fp.decompress(fz, back)
    fp.check_lossless(fz, back)
    fp.check_lossless(src, back)
    assert fits.getheader(back)['BZERO'] == fits.getheader(src)['BZERO']


def test_compress_refuses_a_float_frame(tmp_path):
    src = str(tmp_path / 'proc.fits')
    _float_frame(src)
    with pytest.raises(fp.FramePackError, match='only 16-bit'):
        fp.compress(src, src + '.fz')


def test_lossless_check_catches_a_changed_pixel_and_a_lost_bzero(tmp_path):
    src = str(tmp_path / 'raw.fits')
    data = _raw_frame(src)
    other = str(tmp_path / 'other.fits')
    changed = data.copy()
    changed[1, 1] += 1
    no_sums = fp._strip(_header(), ('CHECKSUM', 'DATASUM'))     # so the pixel check is what fires
    fits.PrimaryHDU(data=changed, header=no_sums).writeto(other)
    with pytest.raises(fp.FramePackError, match='physical pixel values differ'):
        fp.check_lossless(src, other)
    # the trap found on 16 Oct: the same stored integers written without BZERO
    with fits.open(src, do_not_scale_image_data=True) as h:
        fits.PrimaryHDU(data=np.array(h[0].data),
                        header=fp._strip(h[0].header, fp.SCALING + ('CHECKSUM', 'DATASUM')),
                        do_not_scale_image_data=True).writeto(str(tmp_path / 'nobzero.fits'))
    with pytest.raises(fp.FramePackError, match='physical pixel values differ'):
        fp.check_lossless(src, str(tmp_path / 'nobzero.fits'))


def test_t7_header_write_keeps_a_16_bit_frame_16_bit(tmp_path):
    """pointer_wcs.py's pattern: header updated in update mode, pixels read through a separate open.
    Touching the update-mode HDU's data instead makes astropy rewrite the frame as float32 on close."""
    src, f = str(tmp_path / 'proc.fits'), str(tmp_path / 'int16.fits')
    _float_frame(src)
    fp.to_int16(src, f)
    before = fits.getdata(f)
    with fits.open(f, mode='update') as hdu:
        header = hdu[0].header
        data = fits.getdata(f).astype(float)
        header['CRVAL1'] = 123.456
    assert data.shape == before.shape
    h = fits.getheader(f)
    assert (h['BITPIX'], h['CRVAL1']) == (16, 123.456)
    np.testing.assert_array_equal(fits.getdata(f), before)
    with fits.open(f, mode='update') as hdu:      # the old pattern, for contrast
        _ = hdu[0].data.astype(float)
        hdu[0].header['CRVAL2'] = 1.0
    assert fits.getheader(f)['BITPIX'] == -32


# --------------------------------------------------------------------------- in place

def test_pack_in_place_renames_and_keeps_mode_and_mtime(tmp_path):
    src = str(tmp_path / 'procSPECU3.2026-10-13T23:43:22.000.fits')
    _float_frame(src)
    os.chmod(src, 0o640)
    os.utime(src, (1_700_000_000, 1_700_000_000))
    rep = fp.convert_in_place(src, 'pack')
    assert rep.status == 'done', rep.reason
    assert not os.path.exists(src) and rep.final == src + '.fz'
    st = os.stat(rep.final)
    assert (st.st_mode & 0o777, int(st.st_mtime)) == (0o640, 1_700_000_000)
    assert glob.glob(str(tmp_path / '.*')) == []
    assert rep.bzero and rep.n_clipped == '0' and rep.bytes_out < rep.bytes_in


def test_decompress_in_place_restores_the_16_bit_frame(tmp_path):
    src = str(tmp_path / 'proc.fits')
    _float_frame(src)
    kept = str(tmp_path / 'int16.fits')
    fp.to_int16(src, kept)
    assert fp.convert_in_place(src, 'pack').status == 'done'
    rep = fp.convert_in_place(src + '.fz', 'decompress')
    assert rep.status == 'done', rep.reason
    assert os.path.exists(src) and not os.path.exists(src + '.fz')
    np.testing.assert_array_equal(fits.getdata(src), fits.getdata(kept))


def test_pack_compresses_a_frame_that_is_already_16_bit(tmp_path):
    f = str(tmp_path / 'proc.fits')
    _float_frame(str(tmp_path / 'float.fits'))
    fp.to_int16(str(tmp_path / 'float.fits'), f)
    rep = fp.convert_in_place(f, 'pack')
    assert rep.status == 'done' and rep.bzero == ''


def test_int16_skips_a_frame_that_is_already_16_bit(tmp_path):
    f = str(tmp_path / 'raw.fits')
    _raw_frame(f)
    rep = fp.convert_in_place(f, 'int16')
    assert rep.status == 'skipped' and 'already' in rep.reason


def test_a_failure_leaves_the_original_and_no_temporary_files(tmp_path, monkeypatch):
    src = str(tmp_path / 'proc.fits')
    _float_frame(src)
    before = open(src, 'rb').read()

    def broken(a, b):
        raise fp.FramePackError('header cards differ')
    monkeypatch.setattr(fp, 'check_lossless', broken)
    rep = fp.convert_in_place(src, 'pack')
    assert rep.status == 'failed'
    assert open(src, 'rb').read() == before
    assert sorted(os.listdir(tmp_path)) == ['proc.fits']


def test_an_existing_target_is_not_overwritten(tmp_path):
    src = str(tmp_path / 'raw.fits')
    _raw_frame(src)
    open(src + '.fz', 'w').close()
    rep = fp.convert_in_place(src, 'compress')
    assert rep.status == 'skipped' and os.path.getsize(src + '.fz') == 0


def test_wrong_suffixes_are_refused(tmp_path):
    f = str(tmp_path / 'proc.fits.phot')
    open(f, 'w').close()
    assert fp.convert_in_place(f, 'compress').status == 'failed'
    assert fp.convert_in_place(f, 'decompress').status == 'failed'


def test_directories_expand_to_frames_only(tmp_path):
    for name in ('procA.fits', 'procB.fits', 'procA.fits.phot', 'procA.cat', 'procA.fits.fz'):
        open(tmp_path / name, 'w').close()
    assert [os.path.basename(f) for f in fp.expand([str(tmp_path)], 'pack')] == ['procA.fits', 'procB.fits']
    assert [os.path.basename(f) for f in fp.expand([str(tmp_path)], 'decompress')] == ['procA.fits.fz']


def test_command_line_writes_a_manifest(tmp_path):
    d = tmp_path / '1'
    d.mkdir()
    for i in range(3):
        _float_frame(str(d / f'proc{i}.fits'), seed=i)
    manifest = str(tmp_path / 'm.csv')
    assert fp.main(['pack', str(d), '--manifest', manifest]) == 0
    rows = open(manifest).read().splitlines()
    assert len(rows) == 4 and rows[0].startswith('path,op,final,status')
    assert sorted(os.listdir(d)) == [f'proc{i}.fits.fz' for i in range(3)]


# --------------------------------------------------------------------------- byte-exact raw frames and readers

def _eso_raw(path, shape=(96, 128), seed=3):
    """A raw frame laid out as ESO serves it: uint16 as BZERO 32768, a CHECKSUM, cards in ESO's order."""
    rng = np.random.default_rng(seed)
    data = np.clip(1300 + 30 * rng.standard_normal(shape), 0, 65535).astype(np.uint16)
    h = _header()
    del h['CHECKSUM'], h['DATASUM']
    hdu = fits.PrimaryHDU(data=data, header=h)
    hdu.writeto(path, checksum=True)
    return data


def test_rawpack_rebuilds_the_original_byte_for_byte(tmp_path):
    src = str(tmp_path / 'SPECU2.2026-10-16T23:08:07.000.fits')
    data = _eso_raw(src)
    os.utime(src, ns=(1_700_000_000_000_000_000, 1_700_000_000_123_456_789))
    original = open(src, 'rb').read()
    rep = fp.convert_in_place(src, 'rawpack')
    assert rep.status == 'done', rep.reason
    assert not os.path.exists(src) and os.path.exists(src + '.fz')
    assert fp.plain_bytes(src + '.fz') == original
    info = fp.original_info(src + '.fz')
    assert info['ORIGSIZE'] == len(original) and info['ORIGMTIM'] == 1_700_000_000_123_456_789
    rep = fp.convert_in_place(src + '.fz', 'decompress')
    assert rep.status == 'done', rep.reason
    assert open(src, 'rb').read() == original
    assert os.stat(src).st_mtime_ns == 1_700_000_000_123_456_789
    with fits.open(src, checksum=True) as h:          # the original's own CHECKSUM verifies again
        assert np.array_equal(h[0].data, data)


def test_rawpack_refuses_what_it_cannot_rebuild_exactly(tmp_path):
    f32 = str(tmp_path / 'float.fits')
    fits.PrimaryHDU(data=np.ones((8, 8), np.float32)).writeto(f32)
    assert fp.convert_in_place(f32, 'rawpack').status == 'failed'
    two = str(tmp_path / 'two.fits')
    fits.HDUList([fits.PrimaryHDU(data=np.ones((8, 8), np.int16)), fits.ImageHDU(np.ones(4, np.int16))]).writeto(two)
    rep = fp.convert_in_place(two, 'rawpack')
    assert rep.status == 'failed' and 'after the image' in rep.reason
    assert os.path.exists(f32) and os.path.exists(two)
    assert not glob.glob(str(tmp_path / '.*'))


def test_a_corrupted_exact_file_is_refused_on_reading(tmp_path):
    src = str(tmp_path / 'r.fits')
    _eso_raw(src)
    fp.convert_in_place(src, 'rawpack')
    with fits.open(src + '.fz', mode='update') as h:           # tamper with the kept header
        h['ORIGHDR'].header['ORIGSHA'] = '0' * 64
    with pytest.raises(fp.FramePackError):
        fp.plain_bytes(src + '.fz')


@pytest.mark.parametrize('kind', ['raw', 'exact', 'processed'])
def test_open_frame_presents_a_compressed_frame_as_the_plain_file(tmp_path, kind):
    src = str(tmp_path / ('procX.fits' if kind == 'processed' else 'X.fits'))
    if kind == 'processed':
        _float_frame(src)
        fp.convert_in_place(src, 'int16')
    elif kind == 'raw':
        _raw_frame(src)
    else:
        _eso_raw(src)
    with fits.open(src) as h:
        want, want_header = h[0].data.copy(), h[0].header.copy()
    fp.convert_in_place(src, 'rawpack' if kind == 'exact' else 'compress')
    for name in (src, src + '.fz'):                    # the plain name is followed to the .fz
        with fp.open_frame(name) as h:
            assert h.filename() == src + '.fz'
            assert len(h) == 1 and h[0].data.dtype == want.dtype
            assert np.array_equal(h[0].data, want, equal_nan=True)
            assert h[0].header['OBJECT'] == want_header['OBJECT'] and h[0].header['NAXIS2'] == want.shape[0]
        assert fp.getheader(name)['EXPTIME'] == 30.0


def test_open_frame_opens_a_plain_file_as_it_is(tmp_path):
    src = str(tmp_path / 'X.fits')
    data = _raw_frame(src)
    with fp.open_frame(src) as h:
        assert h.filename() == src and np.array_equal(h[0].data, data)
    assert fp.stored_path(str(tmp_path / 'missing.fits')) == str(tmp_path / 'missing.fits')
