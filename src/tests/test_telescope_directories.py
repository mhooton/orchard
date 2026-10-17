"""
Unit tests for telling the telescope from the frame's directory when TELESCOP names none.

ACP/MaxIm-era frames carry no configured telescope name: MaxIm's bias, dark and flat frames
have TELESCOP = '' (Io, Europa and Callisto until September 2018, Ganymede 2019-04-05/06,
Artemis until June 2022, SAINT-EX 2019-2021) and ACP's science frames 'ACP->IO', 'ACP->Io',
'ACP->Europa', 'ACP->M06', 'ACP->Callisto', 'ACP->Ganymede', 'Io' or 'Callisto'.  Without a
match get_instrument_parameters fell back to generic defaults, so T2-T6 neither trimmed nor
subtracted the overscan and T8 stopped with KeyError: 'arcsec_per_pixel' in wcs_status.py.
Since 'ACP->Artemis' was aliased, ACP-era Artemis science frames were trimmed while the
masters made from their '' calibration frames were not.

Run:  pytest src/tests/test_telescope_directories.py -q
"""
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import pipeutils as pu            # noqa: E402

RAW = '/data/SPECULOOSPipeline/Observations/{}/images/20180614/{}'
CALIBRATED = '/data/SPECULOOSPipeline/PipelineOutput/v3/{}/output/20180614/Sp1424-3514/1/{}'


@pytest.fixture(scope='module')
def config():
    return pu.load_instrument_config()


def _write(path, telescop, shape, **cards):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = fits.PrimaryHDU(data=np.zeros(shape, dtype=np.float32))
    h.header['TELESCOP'] = telescop
    h.header['INSTRUME'] = 'Andor Tech'
    for k, v in cards.items():
        h.header[k] = v
    h.writeto(str(path))
    return str(path)


def test_every_directory_points_at_a_configured_telescope(config):
    for directory, name in config['telescope_directories'].items():
        assert name in config['telescopes'], directory


@pytest.mark.parametrize('telescop, directory, expected', [
    ('', 'Io', 'SPECULOOS-IO'),
    ('ACP->IO', 'Io', 'SPECULOOS-IO'),
    ('ACP->Io', 'Io', 'SPECULOOS-IO'),
    ('Io', 'Io', 'SPECULOOS-IO'),
    ('', 'Europa', 'SPECULOOS-EUROPA'),
    ('ACP->M06', 'Europa', 'SPECULOOS-EUROPA'),
    ('ACP->Europa', 'Europa', 'SPECULOOS-EUROPA'),
    ('ACP->Ganymede', 'Ganymede', 'SPECULOOS-GANYMEDE'),
    ('Callisto', 'Callisto', 'SPECULOOS-CALLISTO'),
    ('', 'Artemis', 'Artemis'),
    ('ACP->SAINT-EX', 'Saint-Ex', 'SAINT-EX'),
])
def test_unconfigured_names_take_the_telescope_directory(config, telescop, directory, expected):
    assert pu.resolve_telescope(telescop, config, RAW.format(directory, 'Bias-S001-R001-C001-B1.fts')) == expected
    assert pu.resolve_telescope(telescop, config, CALIBRATED.format(directory, 'procX.fits')) == expected


def test_a_configured_name_always_beats_the_path(config):
    io = RAW.format('Io', 'x.fts')
    assert pu.resolve_telescope('SPECULOOS-EUROPA', config, io) == 'SPECULOOS-EUROPA'
    assert pu.resolve_telescope('ACP->Artemis', config, io) == 'Artemis'


def test_without_a_known_directory_nothing_changes(config):
    assert pu.resolve_telescope('', config) == ''
    assert pu.resolve_telescope('ACP->Io', config, None) == 'ACP->Io'
    assert pu.resolve_telescope('', config, '/tmp/some/where/x.fts') == ''
    # a directory with no entry, and the stale pre-2024 layout
    assert pu.resolve_telescope('', config, RAW.format('Moana-ES', 'x.fits')) == ''
    assert pu.resolve_telescope('', config, '/import/appct/SPECULOOSPipeline/Telescopes/Io/images/20180614/x.fts') == ''
    # a stack outside the per-telescope tree
    assert pu.resolve_telescope('', config, '/data/SPECULOOSPipeline/PipelineOutput/v2/StackImages/x.fits') == ''


def test_the_directory_nearest_the_frame_decides(config):
    # e.g. a scratch replay root that is itself inside someone's Observations tree
    path = '/x/Observations/Io/images/scratch/replay/Observations/Callisto/images/20180615/a.fts'
    assert pu.resolve_telescope('', config, path) == 'SPECULOOS-CALLISTO'
    assert pu.resolve_telescope('', config, 'Observations/Europa/images/20170815/a.fts') == 'SPECULOOS-EUROPA'


def test_acp_era_raw_calibration_frame_is_trimmed_and_overscan_corrected(tmp_path):
    path = _write(tmp_path / 'Observations' / 'Io' / 'images' / '20180614' / 'AutoFlat' / 'AutoFlat-Dusk-I+z-Bin1-001.fts',
                  '', (2088, 2048), IMAGETYP='FLAT')
    with fits.open(path) as hdul:
        assert pu.detect_instrument(hdul) == 'andor'
        p = pu.get_instrument_parameters(hdul)
        assert p == pu.load_instrument_config()['telescopes']['SPECULOOS-IO']['andor']
        assert pu.image_trim(hdul).shape == (2044, 2046)
        assert pu.extract_overscan(hdul).shape == (15, 2044)


def test_acp_era_science_frame_matches_its_calibration_frames(tmp_path):
    night = tmp_path / 'Observations' / 'Callisto' / 'images' / '20180615'
    light = _write(night / 'Sp1424-3514-S001-R001-C001-I+z.fts', 'ACP->Callisto', (2088, 2048), IMAGETYP='Light Frame')
    bias = _write(night / 'Bias-S001-R001-C001-B1.fts', '', (2088, 2048), IMAGETYP='Bias Frame')
    with fits.open(light) as lh, fits.open(bias) as bh:
        assert pu.get_instrument_parameters(lh) == pu.get_instrument_parameters(bh)
        assert pu.image_trim(lh).shape == pu.image_trim(bh).shape == (2044, 2046)


def test_artemis_calibration_frames_match_artemis_science_frames(tmp_path):
    # 'ACP->Artemis' science frames resolved through the alias; their '' masters did not
    night = tmp_path / 'Observations' / 'Artemis' / 'images' / '20210401'
    light = _write(night / 'Sp1122+3755-S001-R001-C001-I+z.fts', 'ACP->Artemis', (2088, 2048))
    dark = _write(night / 'Calibration' / 'Dark-S001-R001-C001-B1.fts', '', (2088, 2048))
    with fits.open(light) as lh, fits.open(dark) as dh:
        assert pu.image_trim(lh).shape == pu.image_trim(dh).shape == (2044, 2046)
        assert pu.get_instrument_parameters(dh)['gain'] == 1.1


def test_callisto_directory_still_tells_its_instruments_apart(tmp_path):
    night = tmp_path / 'Observations' / 'Callisto' / 'images' / '20180615'
    andor = _write(night / 'a.fts', '', (2088, 2048))
    spirit = _write(night / 's.fits', '', (1280, 1024))
    with fits.open(andor) as a, fits.open(spirit) as s:
        assert pu.detect_instrument(a) == 'andor'
        assert pu.detect_instrument(s) == 'spirit'


def test_calibrated_acp_era_frame_passes_the_t8_wcs_check(tmp_path):
    # T6 keeps the raw TELESCOP card; wcs_succeeded used to raise KeyError: 'arcsec_per_pixel'
    pytest.importorskip('fitsio')
    from photometry.wcs_status import wcs_succeeded
    scale, theta = 0.3450 / 3600.0, np.radians(0.3)
    path = _write(tmp_path / 'PipelineOutput' / 'v3' / 'Io' / 'output' / '20180614' / 'Sp1424-3514' / '1' /
                  'procSp1424-3514-S001-R001-C001-I+z.fits', 'ACP->Io', (2044, 2046),
                  CD1_1=-scale * np.cos(theta), CD1_2=scale * np.sin(theta),
                  CD2_1=scale * np.sin(theta), CD2_2=scale * np.cos(theta))
    with fits.open(path) as hdul:
        assert pu.get_instrument_parameters(hdul, trimmed=True)['arcsec_per_pixel'] == 0.34
    assert wcs_succeeded(path)


def test_frames_built_in_memory_still_fall_back():
    hdul = fits.HDUList([fits.PrimaryHDU(data=np.zeros((2088, 8), dtype=np.uint16))])
    hdul[0].header['TELESCOP'] = ''
    assert pu.frame_path(hdul) is None
    p = pu.get_instrument_parameters(hdul)
    assert 'trim' not in p and 'arcsec_per_pixel' not in p
