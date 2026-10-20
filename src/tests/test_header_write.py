"""
T7's header writes to calibrated frames (astrom/header_write.py), used by pointer_wcs after a solve and by
wcs_cache.apply on a cache hit. Covers:
  - CHECKSUM and DATASUM verify after every write: astropy's check and the ones' complement sum of the whole file
  - a DATASUM is trusted, and the data not read, only when the old CHECKSUM agrees with it; otherwise it is summed
    from the data bytes. That covers a stale CHECKSUM, and T6's frames, which carry the raw frame's cards
  - the data bytes never change, so a 16-bit (BZERO) frame stays 16-bit
  - a header that grows or shrinks by a 2880-byte block rewrites the file, keeping its mode; a failed rewrite leaves
    the frame as it was and no temporary file
  - an unchanged header is not written; a frame without checksum cards gets none
  - the header is the one astropy's update mode writes, but for the checksum cards' timestamps

T7's own runs (a solve, then a cache hit) on frames written as T6 writes them are in test_wcs_cache.py.

Run:  pytest src/tests/test_header_write.py -q
"""
import os
import shutil
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from astrom import header_write as hw                    # noqa: E402
from download import raw_provenance as rp                # noqa: E402

SHAPE = (40, 30)
BLOCK = 2880
CARDS_PER_BLOCK = BLOCK // 80


def verifies(path):
    """CHECKSUM and DATASUM present and valid: astropy's check, and the ones' complement sum of the whole file by
    hand (a single HDU)."""
    with fits.open(path, do_not_scale_image_data=True) as h:
        ok = h[0].verify_datasum() == 1 and h[0].verify_checksum() == 1
    with open(path, 'rb') as f:
        return ok and rp.ones_sum(f.read()) == rp.NEGATIVE_ZERO


def data_bytes(path):
    raw, _ = rp.read_header_bytes(path)
    with open(path, 'rb') as f:
        return f.read()[len(raw):]


def raw_frame(path, extra_cards=0):
    """What ESO delivers: uint16 pixels (BITPIX 16, BZERO 32768) with CHECKSUM/DATASUM that verify."""
    data = (1000 + np.arange(np.prod(SHAPE)) % 97).reshape(SHAPE).astype(np.uint16)
    h = fits.Header()
    h.update(IMAGETYP='Light Frame', EXPTIME=30.0, OBJECT='Sp0000+0000')
    for k in range(extra_cards):
        h['KEY{}'.format(k)] = k
    fits.PrimaryHDU(data=data, header=h).writeto(path, checksum=True)
    return path


def t6_frame(path, raw):
    """What T6 writes (calibration/pipered.py): the raw frame opened, its header edited and its pixels replaced by
    float32 ones, written with writeto(checksum=False). The raw frame's CHECKSUM/DATASUM come along, stale."""
    with fits.open(raw) as h:
        h[0].header['OVERSCAN'] = 300.0
        h[0].header.add_history('Overscan of 300.0 subtracted')
        h[0].data = (h[0].data - 300.0).astype(np.float32)
        h.writeto(path)
    return path


def stale(path):
    """Change one card's value in the header bytes and leave CHECKSUM as it was (DATASUM stays right)."""
    raw, header = rp.read_header_bytes(path)
    header['EXPTIME'] = 31.0
    new = header.tostring(padding=True).encode('ascii')
    assert len(new) == len(raw)
    with open(path, 'r+b') as f:
        f.write(new)
    return path


def no_data_sum(*args):
    raise AssertionError('the data were summed although the DATASUM could be trusted')


@pytest.fixture
def sums(monkeypatch):
    """The data_sum calls header_write makes."""
    calls = []
    monkeypatch.setattr(hw, 'data_sum', lambda *a: calls.append(a) or rp.data_sum(*a))
    return calls


def edit(path, **cards):
    with hw.updating(path) as header:
        header.update(cards)


# ------------------------------------------------------------------------------------------------ trust


def test_a_datasum_its_checksum_agrees_with_is_trusted(tmp_path, monkeypatch):
    f = raw_frame(str(tmp_path / 'f.fits'))
    before = data_bytes(f)
    monkeypatch.setattr(hw, 'data_sum', no_data_sum)
    edit(f, CRVAL1=313.7, CRVAL2=9.17)
    assert verifies(f) and data_bytes(f) == before
    assert fits.getheader(f)['CRVAL1'] == 313.7


def test_a_stale_checksum_means_the_data_are_summed(tmp_path, sums):
    f = stale(raw_frame(str(tmp_path / 'f.fits')))
    assert not verifies(f)
    before = data_bytes(f)
    edit(f, CRVAL1=313.7)
    assert len(sums) == 1 and verifies(f) and data_bytes(f) == before


def test_t6s_frames_are_summed_and_their_datasum_corrected(tmp_path, sums):
    raw = raw_frame(str(tmp_path / 'raw.fits'))
    f = t6_frame(str(tmp_path / 'proc.fits'), raw)
    h = fits.getheader(f)
    assert h['BITPIX'] == -32 and h['DATASUM'] == fits.getheader(raw)['DATASUM']      # the raw frame's, stale
    assert not verifies(f)
    before = data_bytes(f)
    edit(f, CRVAL1=313.7)
    assert len(sums) == 1 and verifies(f) and data_bytes(f) == before
    assert fits.getheader(f)['DATASUM'] != fits.getheader(raw)['DATASUM']
    del sums[:]
    edit(f, CRVAL1=313.8)                   # T7's second write to the frame trusts the first one's DATASUM
    assert sums == [] and verifies(f)


@pytest.mark.parametrize('make_stale', [False, True])
def test_a_16_bit_frame_stays_16_bit(tmp_path, make_stale):
    f = raw_frame(str(tmp_path / 'f.fits'))
    if make_stale:
        stale(f)
    before, pixels = data_bytes(f), fits.getdata(f)
    edit(f, CRVAL1=313.7)
    h = fits.getheader(f)
    assert (h['BITPIX'], h['BZERO'], h['BSCALE']) == (16, 32768, 1)
    assert data_bytes(f) == before and verifies(f)
    with fits.open(f) as hdul:
        assert hdul[0].data.dtype == np.uint16
        np.testing.assert_array_equal(hdul[0].data, pixels)


# ------------------------------------------------------------------------------------------------ blocks


def _cards_free(path):
    """Card slots left before the header needs another block."""
    _, header = rp.read_header_bytes(path)
    return CARDS_PER_BLOCK - (len(header) + 1) % CARDS_PER_BLOCK          # +1: the END card


def test_a_header_that_grows_a_block_rewrites_the_file(tmp_path):
    f = raw_frame(str(tmp_path / 'f.fits'))
    os.chmod(f, 0o640)
    size, before, inode = os.path.getsize(f), data_bytes(f), os.stat(f).st_ino
    edit(f, **{'NEW{}'.format(k): k for k in range(_cards_free(f) + 1)})
    assert os.path.getsize(f) == size + BLOCK
    assert data_bytes(f) == before and verifies(f)
    assert os.stat(f).st_mode & 0o7777 == 0o640 and os.stat(f).st_ino != inode
    assert os.listdir(str(tmp_path)) == ['f.fits']


def test_a_header_that_shrinks_a_block_rewrites_the_file(tmp_path):
    f = raw_frame(str(tmp_path / 'f.fits'), extra_cards=CARDS_PER_BLOCK)
    size, before = os.path.getsize(f), data_bytes(f)
    with hw.updating(f) as header:
        for k in range(CARDS_PER_BLOCK):
            del header['KEY{}'.format(k)]
    assert os.path.getsize(f) == size - BLOCK
    assert data_bytes(f) == before and verifies(f)


def test_a_header_that_keeps_its_blocks_is_written_in_place(tmp_path):
    f = raw_frame(str(tmp_path / 'f.fits'))
    size, inode = os.path.getsize(f), os.stat(f).st_ino
    edit(f, CRVAL1=313.7)
    assert (os.path.getsize(f), os.stat(f).st_ino) == (size, inode) and verifies(f)


def test_a_failed_rewrite_leaves_the_frame_and_no_temporary_file(tmp_path, monkeypatch):
    f = raw_frame(str(tmp_path / 'f.fits'))
    with open(f, 'rb') as fh:
        before = fh.read()

    def broken(*args, **kwargs):
        raise OSError('disk full')
    monkeypatch.setattr(hw.shutil, 'copyfileobj', broken)
    with pytest.raises(OSError, match='disk full'):
        edit(f, **{'NEW{}'.format(k): k for k in range(_cards_free(f) + 1)})
    with open(f, 'rb') as fh:
        assert fh.read() == before
    assert os.listdir(str(tmp_path)) == ['f.fits']


# ------------------------------------------------------------------------------------------------ when and what


def test_an_unchanged_header_is_not_written(tmp_path):
    f = raw_frame(str(tmp_path / 'f.fits'))
    with open(f, 'rb') as fh:
        before = fh.read()
    st = os.stat(f)
    with hw.updating(f) as header:
        header['EXPTIME'] = header['EXPTIME']
    with open(f, 'rb') as fh:
        assert fh.read() == before
    assert os.stat(f).st_mtime_ns == st.st_mtime_ns


def test_a_frame_without_checksums_gets_none(tmp_path, monkeypatch):
    f = str(tmp_path / 'f.fits')
    fits.PrimaryHDU(data=np.zeros(SHAPE, np.float32)).writeto(f)
    before = data_bytes(f)
    monkeypatch.setattr(hw, 'data_sum', no_data_sum)
    edit(f, CRVAL1=313.7)
    h = fits.getheader(f)
    assert 'CHECKSUM' not in h and 'DATASUM' not in h and h['CRVAL1'] == 313.7
    assert data_bytes(f) == before


def test_the_header_is_written_when_the_block_raises(tmp_path):
    """As astropy's update mode did: pointer_wcs clears an inherited WCS before a solve that may fail."""
    f = raw_frame(str(tmp_path / 'f.fits'))
    with pytest.raises(RuntimeError):
        with hw.updating(f) as header:
            header['CRVAL1'] = 313.7
            raise RuntimeError('solve failed')
    assert fits.getheader(f)['CRVAL1'] == 313.7 and verifies(f)


def _cards(path):
    """The header as (keyword, value) pairs, but for the checksum cards' values, which carry a timestamp."""
    h = fits.getheader(path)
    return [(c.keyword, c.value) for c in h.cards if c.keyword != 'CHECKSUM'], h['DATASUM']


@pytest.mark.parametrize('make_stale', [False, True])
def test_the_header_is_the_one_astropys_update_mode_writes(tmp_path, make_stale):
    a = raw_frame(str(tmp_path / 'a.fits'))
    if make_stale:
        stale(a)
    b = str(tmp_path / 'b.fits')
    shutil.copy(a, b)
    new = {'CTYPE1': 'RA---TAN', 'CRVAL1': 313.7, 'CRPIX1': 15.5}
    with fits.open(a, mode='update', do_not_scale_image_data=True) as hdul:
        del hdul[0].header['OBJECT']
        hdul[0].header.update(new)
    with hw.updating(b) as header:
        del header['OBJECT']
        header.update(new)
    assert _cards(a) == _cards(b)
    assert data_bytes(a) == data_bytes(b) and verifies(a) and verifies(b)
