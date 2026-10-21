"""
Unit tests for utils/rawpack.py: a compressed raw night is decompressed for a pipeline run and stored compressed
again afterwards, without ever being left without a stored copy of a frame.

Covers:
  - unpack: working copies byte for byte (byte-exact .fz) with the original mtime, the .fz hidden from T1's globs
  - repack: the .fz back, the working copy gone, nothing else left behind
  - a run killed after unpack, or half-way through repack, is taken over by the next unpack/repack
  - a working copy the run changed is compressed again and the earlier .fz kept
  - .fz files without the kept header (the 2019-20 nights, fpack 2022) and frames one directory down
  - mixed nights (plain and .fz frames) and a plain file that differs from its .fz are left alone

Run:  pytest src/tests/test_rawpack.py -q
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

from utils import framepack as fp     # noqa: E402
from utils import rawpack as rp       # noqa: E402

warnings.simplefilter('ignore', category=fits.verify.VerifyWarning)

MTIME = 1_700_000_000_123_456_789


def _raw(path, seed=0):
    rng = np.random.default_rng(seed)
    data = np.clip(1300 + 30 * rng.standard_normal((40, 48)), 0, 65535).astype(np.uint16)
    h = fits.Header()
    h['IMAGETYP'] = 'Light Frame'
    h['EXPTIME'] = (20.0, 'exposure time')
    fits.PrimaryHDU(data=data, header=h).writeto(path, checksum=True)
    os.utime(path, ns=(MTIME, MTIME))
    return open(path, 'rb').read()


def _night(tmp_path, n=3, op='rawpack'):
    night = tmp_path / 'images' / '20261016'
    night.mkdir(parents=True)
    originals = {}
    for k in range(n):
        p = str(night / 'SPECU2.2026-10-16T23:0{}:00.000.fits'.format(k))
        originals[p] = _raw(p, k)
        assert fp.convert_in_place(p, op).status == 'done'
    return str(night), originals


def _visible(night):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(night, '*')))


def test_unpack_then_repack_restores_the_compressed_night(tmp_path):
    night, originals = _night(tmp_path)
    n, failures = rp.unpack(night, workers=2)
    assert n == 3 and not failures
    for p, raw in originals.items():
        assert open(p, 'rb').read() == raw and os.stat(p).st_mtime_ns == MTIME
        assert not os.path.exists(p + '.fz')
    assert all(name.endswith('.fits') for name in _visible(night))      # what T1 sees
    assert rp.status(night) == {'compressed': 0, 'unpacked': 3, 'plain': 3}

    counts, messages = rp.repack(night, workers=2)
    assert counts == {'repacked': 3} and not messages
    assert all(name.endswith('.fits.fz') for name in _visible(night))
    assert not os.path.exists(os.path.join(night, rp.STATE))
    for p, raw in originals.items():
        assert fp.plain_bytes(p + '.fz') == raw


def test_unpack_and_repack_repeat_safely(tmp_path):
    night, originals = _night(tmp_path)
    rp.unpack(night)
    assert rp.unpack(night) == (3, [])          # a second unpack (a killed run's successor) takes the copies over
    counts, _ = rp.repack(night)
    assert counts == {'repacked': 3}
    assert rp.repack(night) == ({}, [])         # nothing left to do
    assert rp.status(night) == {'compressed': 3, 'unpacked': 0, 'plain': 0}


def test_a_lost_working_copy_is_made_again_and_a_half_repack_is_finished(tmp_path):
    night, originals = _night(tmp_path)
    rp.unpack(night)
    first, second, third = sorted(originals)
    os.remove(first)                                           # a working copy lost after a kill
    os.replace(os.path.join(night, rp.STATE, os.path.basename(second) + '.fz'), second + '.fz')   # repack cut short
    n, failures = rp.unpack(night)
    assert n == 3 and not failures
    assert open(first, 'rb').read() == originals[first]
    counts, _ = rp.repack(night)
    assert counts == {'repacked': 3}
    assert rp.status(night)['compressed'] == 3


def test_a_working_copy_the_run_changed_is_stored_and_the_old_fz_kept(tmp_path):
    night, originals = _night(tmp_path, n=2)
    rp.unpack(night)
    p = sorted(originals)[0]
    with fits.open(p, mode='update') as h:
        h[0].header['CRVAL1'] = 12.5
    counts, messages = rp.repack(night)
    assert counts == {'recompressed': 1, 'repacked': 1}
    assert 'changed by the run' in messages[0]
    assert fp.getheader(p + '.fz')['CRVAL1'] == 12.5
    kept = glob.glob(os.path.join(night, rp.STATE, os.path.basename(p) + '.fz.*'))
    assert len(kept) == 1 and fp.plain_bytes(kept[0]) == originals[p]
    assert not os.path.exists(p)


def test_fz_without_the_kept_header_and_frames_one_directory_down(tmp_path):
    night, originals = _night(tmp_path, n=2, op='compress')     # like the 2019-20 nights compressed with fpack
    sub = os.path.join(night, 'calib')
    os.mkdir(sub)
    q = os.path.join(sub, 'bias.fits')
    _raw(q, 9)
    fp.convert_in_place(q, 'rawpack')
    n, failures = rp.unpack(night)
    assert n == 3 and not failures
    for p in originals:
        kept = os.path.join(night, rp.STATE, os.path.basename(p) + '.fz')
        with fits.open(p, checksum=True) as a, fp.open_frame(kept) as b:
            assert np.array_equal(a[0].data, b[0].data)
            # the working copy is the .fz's image as stored, with no edit of its own recorded on top
            assert len(a) == 1 and a[0].header.get('ORCHMOD') == b[0].header.get('ORCHMOD')
    assert os.path.exists(q)
    counts, _ = rp.repack(night)
    assert counts == {'repacked': 3}
    assert os.path.exists(q + '.fz') and not os.path.exists(q)


def test_mixed_nights_and_conflicts(tmp_path):
    night, originals = _night(tmp_path, n=2)
    plain = os.path.join(night, 'SPECU2.2026-10-16T23:59:00.000.fits')
    raw = _raw(plain, 7)                                     # a frame that was never compressed
    clash = sorted(originals)[0]
    _raw(clash, 99)                                          # a plain copy that differs from its .fz
    n, failures = rp.unpack(night)
    assert n == 1 and len(failures) == 1 and 'differ' in failures[0]
    assert os.path.exists(clash + '.fz') and open(plain, 'rb').read() == raw
    counts, _ = rp.repack(night)
    assert counts == {'repacked': 1}
    assert open(plain, 'rb').read() == raw                   # untouched throughout
