"""
Unit tests for telescope name aliases in the instrument configuration.

Artemis frames up to 2022-06-15 carry TELESCOP = 'ACP->Artemis'; later ones 'Artemis'.
Before the alias, the older frames matched no configured telescope, so
get_instrument_parameters fell back to generic defaults (no trim, no overscan, no
plate scale) and T8 stopped with KeyError: 'arcsec_per_pixel' in wcs_status.py.

Run:  pytest src/tests/test_telescope_aliases.py -q
"""
import os
import sys

import numpy as np
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import pipeutils as pu            # noqa: E402


def _hdul(telescop, naxis2):
    h = fits.PrimaryHDU(data=np.zeros((naxis2, 8), dtype=np.uint16))
    h.header['TELESCOP'] = telescop
    return fits.HDUList([h])


def test_every_alias_points_at_a_configured_telescope():
    config = pu.load_instrument_config()
    for old, new in config['telescope_aliases'].items():
        assert new in config['telescopes'], old


def test_resolve_telescope():
    config = pu.load_instrument_config()
    assert pu.resolve_telescope('ACP->Artemis', config) == 'Artemis'
    assert pu.resolve_telescope('Artemis ', config) == 'Artemis'
    assert pu.resolve_telescope('SPECULOOS-IO', config) == 'SPECULOOS-IO'
    assert pu.resolve_telescope(None, config) == ''


def test_acp_era_artemis_raw_frame_gets_artemis_parameters():
    old, new = _hdul('ACP->Artemis', 2088), _hdul('Artemis', 2088)
    assert pu.detect_instrument(old) == 'andor'
    p = pu.get_instrument_parameters(old)
    assert p == pu.get_instrument_parameters(new)
    assert p['arcsec_per_pixel'] == 0.34 and 'trim' in p and 'overscan' in p


def test_acp_era_artemis_calibrated_frame_gets_artemis_parameters():
    # calibrated frames are trimmed to rows 22..2066, and keep the raw TELESCOP card
    old = _hdul('ACP->Artemis', 2044)
    assert pu.detect_instrument(old, trimmed=True) == 'andor'
    assert pu.get_instrument_parameters(old, trimmed=True)['arcsec_per_pixel'] == 0.34


def test_an_unknown_telescope_still_falls_back():
    p = pu.get_instrument_parameters(_hdul('Somewhere-Else', 2088))
    assert 'arcsec_per_pixel' not in p and p['saturation_threshold'] == 62000
