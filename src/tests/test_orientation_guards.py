"""
Tests for the Andor orientation guards (audit of 2026-10-16).

Calibration assumes Andor frames in the stored MaxIm orientation, 2048 x 2088. 7,524 raw frames were found in the
native 2088 x 2048 shape: Callisto frames never rotated (np.rot90 k=1 turns them) and Ganymede frames rotated
twice, 180 deg from ESO (k=3), because astra_transform turned frames that already carried ASTRAROT. Covers:
  - download/repair_orientation.py turns both kinds into the stored frame exactly, keeps BITPIX/BZERO, removes the
    WCS, keeps the original in quarantine, writes a manifest, and refuses a frame it has already turned
  - astra_transform skips frames that already carry ASTRAROT / ASTRAMIR unless forced
  - createlists (T1) refuses a night holding a native-shaped frame

Run:  pytest src/tests/test_orientation_guards.py -q
"""
import csv
import hashlib
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from download import repair_orientation as ro                # noqa: E402
from download.astra_transform import astra_transform_unified  # noqa: E402

STORED = (2088, 2048)    # numpy shape of a stored-orientation Andor frame


def stored_frame(seed=0):
    rng = np.random.default_rng(seed)
    return rng.integers(250, 60000, size=STORED, dtype=np.uint16)


def write(path, data, **cards):
    h = fits.Header()
    h['IMAGETYP'] = cards.pop('IMAGETYP', 'Light Frame')
    h['OBSERVER'] = 'Astra'
    h['FILTER'] = 'I+z'
    for k, v in cards.items():
        h[k] = v
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fits.PrimaryHDU(data=data, header=h).writeto(path)   # uint16 -> BITPIX 16, BZERO 32768
    return path


def wcs_cards():
    return dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CRVAL1=150.0, CRVAL2=-30.0, CRPIX1=1044.0, CRPIX2=1024.0,
                PC1_1=-9.6e-5, PC1_2=0.0, PC2_1=0.0, PC2_2=9.6e-5, IMAGEW=2088, IMAGEH=2048)


def md5(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


def frames_csv(tmp_path, rows):
    p = tmp_path / 'frames.csv'
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['file', 'k'])
        w.writerows(rows)
    return str(p)


@pytest.mark.parametrize('k', [1, 3])
def test_repair_turns_a_native_frame_into_the_stored_one(tmp_path, k):
    s = stored_frame(k)
    native = np.rot90(s, -k)                    # what is on disk: rot90(native, k) gives the stored frame
    night = tmp_path / 'Observations' / 'Callisto' / 'images' / '20240408'
    path = write(str(night / 'SPECU4.x.fits'), native, **wcs_cards())
    before = md5(path)
    q, man = tmp_path / 'quarantine', tmp_path / 'manifest.csv'
    rc = ro.main([frames_csv(tmp_path, [(path, f'k{k}')]), '--apply', '--quarantine', str(q), '--manifest', str(man)])
    assert rc == 0
    with fits.open(path) as hl:
        assert np.array_equal(hl[0].data, s)
        h = hl[0].header
    raw = fits.getheader(path)
    assert (raw['NAXIS1'], raw['NAXIS2']) == (2048, 2088)
    assert raw['BITPIX'] == 16 and raw['BZERO'] == 32768
    assert h['ASTRAROT'] == 90
    assert not any(c in h for c in ('CTYPE1', 'CRVAL1', 'CRPIX1', 'PC1_1'))
    assert (h['IMAGEW'], h['IMAGEH']) == (2048, 2088)
    assert any('orient repair: rot90 k={}'.format(k) in str(x) for x in h['HISTORY'])
    kept = q / 'Callisto' / '20240408' / 'SPECU4.x.fits'
    assert md5(kept) == before
    rows = list(csv.DictReader(open(man)))
    assert len(rows) == 1 and rows[0]['old_md5'] == before and rows[0]['new_md5'] == md5(path)
    assert not [n for n in os.listdir(night) if n.endswith('.tmp')]
    # a second run refuses the frame: it is no longer in the native shape
    after = md5(path)
    rc = ro.main([frames_csv(tmp_path, [(path, f'k{k}')]), '--apply', '--quarantine', str(q), '--manifest', str(man)])
    assert rc == 1 and md5(path) == after


def test_report_only_changes_nothing(tmp_path):
    path = write(str(tmp_path / 'Observations' / 'Ganymede' / 'images' / '20251129' / 'f.fits'), np.rot90(stored_frame(), -3))
    before = md5(path)
    assert ro.main([frames_csv(tmp_path, [(path, 3)])]) == 0
    assert md5(path) == before


def test_refuses_k2_and_stored_frames(tmp_path):
    native = write(str(tmp_path / 'a.fits'), np.rot90(stored_frame(), -1))
    stored = write(str(tmp_path / 'b.fits'), stored_frame())
    assert ro.check(native, 2) is not None
    assert ro.check(stored, 1) is not None
    assert ro.check(native, 1) is None


def spirit_frame(seed=0):
    rng = np.random.default_rng(seed)
    return rng.integers(250, 60000, size=(1280, 1024), dtype=np.uint16)


def mirrored_twice(path, data, imagetyp='Light Frame', cards=2):
    write(path, data, IMAGETYP=imagetyp, ASTRAMIR='Y', CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CRVAL1=150.0,
          CRVAL2=-30.0, CRPIX1=512.0, CRPIX2=640.0, PC1_1=-8.6e-5, PC1_2=0.0, PC2_1=0.0, PC2_2=8.6e-5)
    with fits.open(path, mode='update') as hl:
        for _ in range(cards):
            hl[0].header.add_history('WCS updated for horizontal mirroring')
    return path


def test_flip_mirrors_a_doubly_mirrored_spirit_frame_back(tmp_path):
    s = spirit_frame(4)
    night = tmp_path / 'Observations' / 'Callisto' / 'images' / '20251222'
    path = mirrored_twice(str(night / 'SPECU4.s.fits'), np.flip(s, axis=1))
    before = md5(path)
    q, man = tmp_path / 'quarantine', tmp_path / 'manifest.csv'
    args = [frames_csv(tmp_path, [(path, 'flip')]), '--apply', '--quarantine', str(q), '--manifest', str(man)]
    assert ro.main(args) == 0
    with fits.open(path) as hl:
        assert np.array_equal(hl[0].data, s)
        h = hl[0].header
    raw = fits.getheader(path)
    assert (raw['NAXIS1'], raw['NAXIS2']) == (1024, 1280)
    assert raw['BITPIX'] == 16 and raw['BZERO'] == 32768
    assert h['ASTRAMIR'] == 'Y' and 'ASTRAROT' not in h
    assert not any(c in h for c in ('CTYPE1', 'CRVAL1', 'PC1_1'))
    assert any('orient repair: fliplr' in str(x) for x in h['HISTORY'])
    assert md5(q / 'Callisto' / '20251222' / 'SPECU4.s.fits') == before
    assert list(csv.DictReader(open(man)))[0]['k'] == 'flip'
    # a second run refuses it, so it cannot be flipped back again
    after = md5(path)
    assert ro.main(args) == 1 and md5(path) == after


def test_flip_refuses_a_light_frame_mirrored_once_and_andor_frames(tmp_path):
    once = mirrored_twice(str(tmp_path / 'once.fits'), spirit_frame(), cards=1)
    assert 'mirroring card' in ro.check(once, 'flip')
    andor = write(str(tmp_path / 'andor.fits'), stored_frame(), ASTRAMIR='Y')
    assert 'not a SPIRIT frame' in ro.check(andor, 'flip')
    dark = mirrored_twice(str(tmp_path / 'dark.fits'), spirit_frame(), imagetyp='Dark Frame', cards=0)
    assert ro.check(dark, 'flip') is None       # calibration frames carry no mirroring card


def test_transform_does_not_turn_a_frame_twice(tmp_path):
    native = np.rot90(stored_frame(), -1)
    write(str(tmp_path / 'f.fits'), native, IMAGETYP='Dark Frame')
    astra_transform_unified(['f'], transform_mode='auto', n_repeats=1, dirname=str(tmp_path))
    once = fits.getdata(tmp_path / 'f.fits')
    assert once.shape == STORED and fits.getheader(tmp_path / 'f.fits')['ASTRAROT'] == 90
    before = md5(tmp_path / 'f.fits')
    astra_transform_unified(['f'], transform_mode='auto', n_repeats=1, dirname=str(tmp_path))
    assert md5(tmp_path / 'f.fits') == before
    astra_transform_unified(['f'], transform_mode='rotate', n_repeats=1, dirname=str(tmp_path), force=True)
    assert fits.getdata(tmp_path / 'f.fits').shape == STORED[::-1]


def test_transform_does_not_mirror_a_frame_twice(tmp_path):
    rng = np.random.default_rng(5)
    write(str(tmp_path / 's.fits'), rng.integers(0, 60000, size=(1280, 1024), dtype=np.uint16), IMAGETYP='Dark Frame')
    astra_transform_unified(['s'], transform_mode='auto', n_repeats=1, dirname=str(tmp_path))
    before = md5(tmp_path / 's.fits')
    astra_transform_unified(['s'], transform_mode='auto', n_repeats=1, dirname=str(tmp_path))
    assert md5(tmp_path / 's.fits') == before


def test_createlists_refuses_a_native_frame(tmp_path):
    pytest.importorskip('tabulate')            # createlists imports reporting.QC
    from main import createlists
    good = write(str(tmp_path / 'good.fits'), stored_frame(), IMAGETYP='Bias Frame')
    assert createlists.sort_liste([good], 'run', 1)[0] == [good]
    bad = write(str(tmp_path / 'bad.fits'), np.rot90(stored_frame(1), -1), IMAGETYP='Bias Frame')
    with pytest.raises(SystemExit) as e:
        createlists.sort_liste([good, bad], 'run', 1)
    assert 'native Andor orientation' in str(e.value) and 'bad.fits' in str(e.value)
