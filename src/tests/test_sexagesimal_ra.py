"""
Unit tests for reading a sexagesimal RA card as hours in T7's plate solve.

ACP wrote RA as hours, minutes and seconds ('20 54 45.00'). Artemis frames still do (the
instrument config says hourangle), and SSO and SAINT-EX frames did until the transfer began
rewriting RA in degrees in autumn 2018. pointer_wcs took the unit from the config alone, so
an ACP-era SSO frame without a PinPoint WCS was centred at 20.9 deg instead of 313.7 deg:
every frame of Io 20180614 (1,673, none solved by PinPoint) queried Gaia in the wrong part
of the sky and failed with 0 matches.

Run:  pytest src/tests/test_sexagesimal_ra.py -q
"""
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

pw = pytest.importorskip('astrom.pointer_wcs')


@pytest.mark.parametrize('ra, expected', [
    ('20 54 45.00', True),
    ('20:54:45.0', True),
    (' 05 16 53.59 ', True),
    (313.6875, False),
    ('313.6875', False),
    ('', False),
    (None, False),
])
def test_ra_is_sexagesimal(ra, expected):
    assert pw.ra_is_sexagesimal(ra) is expected


def _centre_of(path, monkeypatch):
    """The field centre pointer_wcs hands to the Gaia query for this frame."""
    asked = []

    def fake_query(center, **kwargs):
        asked.append(center)
        raise RuntimeError('stop after the query')

    stars = np.array([[20.0 + 15 * k, 30.0 + 11 * k] for k in range(8)])
    monkeypatch.setattr(pw, 'clean_image', lambda data: data)
    monkeypatch.setattr(pw, 'find_stars_multiscale', lambda *a, **k: stars)
    monkeypatch.setattr(pw, 'gaia_db_query', fake_query)
    result = pw.pointer_wcs(path, db_path='/nonexistent.db')
    assert asked, f'the Gaia query was never reached: {result}'
    return asked[0]


def _frame(path, telescop, ra, dec):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = fits.PrimaryHDU(data=np.zeros((200, 200), dtype=np.float32))
    h.header.update(TELESCOP=telescop, IMAGETYP='Light Frame', RA=ra, DEC=dec,
                    XPIXSZ=13.5, FOCALLEN=8000.0, **{'DATE-OBS': '2018-06-15T03:00:00'})
    h.writeto(str(path))
    return str(path)


def test_acp_era_sso_frame_is_centred_in_hours(tmp_path, monkeypatch):
    path = _frame(tmp_path / 'PipelineOutput' / 'v3' / 'Io' / 'output' / '20180614' / 'WASP-176b' / '1' /
                  'procWASP-176b-S001-R001-C001-I+z.fits', 'ACP->Io', '20 54 45.00', '+09 10 29.0')
    ra, dec = _centre_of(path, monkeypatch)
    assert ra == pytest.approx(313.6875) and dec == pytest.approx(9.174722, abs=1e-5)


def test_sso_frame_in_degrees_is_unchanged(tmp_path, monkeypatch):
    path = _frame(tmp_path / 'PipelineOutput' / 'v3' / 'Io' / 'output' / '20181113' / 'T' / '1' / 'procT.fits',
                  'SPECULOOS-IO', 313.6875, 9.174722)
    assert _centre_of(path, monkeypatch)[0] == pytest.approx(313.6875)


def test_artemis_frame_is_unchanged(tmp_path, monkeypatch):
    path = _frame(tmp_path / 'PipelineOutput' / 'v3' / 'Artemis' / 'output' / '20211101' / 'T' / '1' / 'procT.fits',
                  'ACP->Artemis', '05 16 53.59', '+56 40 21.0')
    assert _centre_of(path, monkeypatch)[0] == pytest.approx(79.223292, abs=1e-5)
