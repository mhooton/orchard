"""condense_photometry.sort_by_jd leaves out the .phot files it cannot read.

imcore_list died (SIGSEGV) on seven dawn frames of Callisto 20181025 and left two zero-byte .phot files; reading
their headers stopped the condense step and, under errexit, the rest of the night (job 4473, twice)."""
import os
import sys

import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pytest.importorskip('fitsio')
from condense import condense_photometry as cp      # noqa: E402


def phot(path, jd):
    h = fits.Header()
    h['JD-OBS'] = jd
    fits.PrimaryHDU(header=h).writeto(path)
    return path


def test_an_empty_phot_file_is_left_out_of_the_night(tmp_path, caplog):
    late = phot(str(tmp_path / 'procB.fits.phot'), 2458417.90)
    early = phot(str(tmp_path / 'procA.fits.phot'), 2458417.60)
    empty = str(tmp_path / 'procC.fits.phot')
    open(empty, 'w').close()
    with caplog.at_level('WARNING'):
        assert cp.sort_by_jd([late, empty, early]) == [early, late]
    assert 'Excluding ' + empty in caplog.text
    assert cp.sort_by_jd([empty]) == []
