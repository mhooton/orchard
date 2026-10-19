"""
Unit tests for the BZERO decision in download/unpack_datacubes.py.

SPIRIT cubes made before create_datacubes.py wrote BZERO come in two kinds: counts shifted by -32768 into the FITS
unsigned 16-bit convention, and the counts themselves. Forcing BZERO = 32768 onto every such cube (a9c4f4a) put every
pixel of Callisto 2026-03-10..25 32768 ADU too high (stored minimum 265, read as 33033). Covers:
  - a shifted cube gets BZERO = 32768, an unshifted one BZERO = 0
  - a header that already has BZERO, or data that are not 16-bit, is left to itself
  - a frame written with the chosen BZERO reads back as the counts

Run:  pytest src/tests/test_unpack_bzero.py -q
"""
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from download.unpack_datacubes import missing_bzero, write_frame  # noqa: E402

COUNTS = np.array([[265, 3052], [16128, 900]], dtype=np.int32)    # SPIRIT counts, 0..16383


def test_shifted_cube_gets_32768():
    stored = (COUNTS - 32768).astype(np.int16)
    assert missing_bzero(fits.Header({'BITPIX': 16}), stored) == 32768


def test_unshifted_cube_gets_zero():
    stored = COUNTS.astype(np.int16)
    assert missing_bzero(fits.Header({'BITPIX': 16}), stored) == 0


def test_existing_bzero_or_other_bitpix_left_alone():
    stored = COUNTS.astype(np.int16)
    assert missing_bzero(fits.Header({'BITPIX': 16, 'BZERO': 32768}), stored) is None
    assert missing_bzero(fits.Header({'BITPIX': 32}), COUNTS) is None


def test_written_frame_reads_back_as_counts(tmp_path):
    for stored in ((COUNTS - 32768).astype(np.int16), COUNTS.astype(np.int16)):
        bzero = missing_bzero(fits.Header({'BITPIX': 16}), stored)
        path = str(tmp_path / 'frame_{}.fits'.format(bzero))
        write_frame(path, stored, fits.Header([('BZERO', bzero), ('BSCALE', 1)]))   # as unpack_datacube writes it
        assert np.array_equal(fits.getdata(path).astype(np.int32), COUNTS)
        with fits.open(path, do_not_scale_image_data=True) as h:
            assert np.array_equal(h[0].data, stored) and h[0].header['BZERO'] == bzero
