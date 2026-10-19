"""
Tests for raw-frame provenance (docs/versioning.md, "Raw frames"). Covers:
  - ingest: astra_transform() on synthetic Astra Andor and SPIRIT frames, lights and calibrations, and unpack_datacube
    on a small cube stamp each frame with ORCHVER and one 'orchard <version> <script>' card, set no ORCHMOD, keep the
    descriptive HISTORY cards as they were, and replace ESO's CHECKSUM/DATASUM (stale once the image is turned) with
    ones that verify
  - later edits: the astra_transform command line, repair_astra_wcs and repair_orientation set ORCHMOD and their own
    card, keep ORCHVER and their marker cards, and leave checksums valid; repair_astra_wcs recomputes CHECKSUM from
    the frame's DATASUM when the old CHECKSUM vouches for it, and from the data bytes when it does not
  - the repair classifiers (repair_astra_wcs.diagnose, repair_orientation.check) answer the same with the provenance
    cards as without them
  - every step run twice changes nothing the second time
  - the hand-made checksums agree with astropy's

Run:  pytest src/tests/test_raw_provenance.py -q
"""
import hashlib
import os
import shutil
import sys
import tarfile

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from download import astra_transform as at               # noqa: E402
from download import raw_provenance as rp                # noqa: E402
from download import repair_astra_wcs as raw_wcs         # noqa: E402
from download import repair_orientation as ro            # noqa: E402
from download.astra_transform import plate_pa            # noqa: E402
from download.unpack_datacubes import unpack_datacube    # noqa: E402
from utils import orchard_version as ov                  # noqa: E402

V = ov.version()
INGEST = 'orchard {} astra_transform'.format(V)
UNPACK = 'orchard {} unpack_datacubes'.format(V)

ANDOR_SMALL = (6, 9)          # numpy shape (NAXIS2, NAXIS1): a small stand-in for a native Andor frame
ANDOR_NATIVE = (2048, 2088)   # the real one, for repair_orientation's shape checks
SPIRIT_SMALL = (8, 1024)      # NAXIS1 = 1024 is how astra_transform tells SPIRIT
SPIRIT = (1280, 1024)

ROT_CARD = 'WCS updated for 90 deg  anticlockwise rotation'
MIR_CARD = 'WCS updated for horizontal mirroring'


# --- helpers -------------------------------------------------------------------------------------------------------

def orchard_cards(header):
    return [str(h) for h in header.get('HISTORY', []) if str(h).startswith('orchard ')]


def other_cards(header):
    return [str(h) for h in header.get('HISTORY', []) if not str(h).startswith('orchard ')]


def verifies(path):
    """DATASUM and CHECKSUM both present and valid: astropy's check, and the ones' complement sum of the whole file
    (a single HDU) by hand."""
    with fits.open(path, do_not_scale_image_data=True) as h:
        ok = h[0].verify_datasum() == 1 and h[0].verify_checksum() == 1
    with open(path, 'rb') as f:
        return ok and rp.ones_sum(f.read()) == rp.NEGATIVE_ZERO


def md5(path):
    with open(path, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def astra_header(imagetyp='Light Frame', shape=ANDOR_SMALL):
    h = fits.Header()
    h['IMAGETYP'] = imagetyp
    h['OBSERVER'] = 'Astra'
    h['FILTER'] = 'I+z'
    if imagetyp == 'Light Frame':
        ny, nx = shape
        for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CRVAL1=150.0, CRVAL2=-30.0, CRPIX1=nx / 2,
                         CRPIX2=ny / 2, PC1_1=-9.4e-6, PC1_2=-9.6e-5, PC2_1=9.6e-5, PC2_2=-9.4e-6).items():
            h[k] = v
    return h


def eso_frame(path, shape=ANDOR_SMALL, imagetyp='Light Frame', seed=0):
    """A raw Astra frame as ESO delivers it: uint16 (BITPIX 16, BZERO 32768) with ESO's CHECKSUM/DATASUM."""
    data = np.random.default_rng(seed).integers(200, 60000, size=shape, dtype=np.uint16)
    fits.PrimaryHDU(data=data, header=astra_header(imagetyp, shape)).writeto(str(path), checksum=True)
    assert verifies(str(path))
    return data


def strip_provenance(path):
    """A copy of path without ORCHVER/ORCHREL/ORCHMOD and the 'orchard ...' cards: the frame as it would be without
    this PR's records."""
    out = str(path) + '.bare.fits'
    shutil.copy(str(path), out)
    with fits.open(out, mode='update', do_not_scale_image_data=True) as hl:
        h = hl[0].header
        for k in ('ORCHVER', 'ORCHREL', 'ORCHMOD'):
            h.remove(k, ignore_missing=True)
        for i in reversed([i for i, c in enumerate(h.cards)
                           if c.keyword == 'HISTORY' and str(c.value).startswith('orchard ')]):
            del h[i]
    return out


# --- ingest --------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('instrument, imagetyp', [('andor', 'Light Frame'), ('andor', 'Dark Frame'),
                                                  ('spirit', 'Light Frame'), ('spirit', 'Flat Field')])
def test_ingest_transform_stamps_the_frame_and_its_checksums_verify(tmp_path, instrument, imagetyp):
    path = tmp_path / 'frame.fits'
    data = eso_frame(path, ANDOR_SMALL if instrument == 'andor' else SPIRIT_SMALL, imagetyp)
    at.astra_transform(['frame'], str(tmp_path), 1)
    with fits.open(path) as hl:
        h, out = hl[0].header, hl[0].data
    assert np.array_equal(out, np.rot90(data) if instrument == 'andor' else np.flip(data, axis=1))
    assert ('ASTRAROT' if instrument == 'andor' else 'ASTRAMIR') in h        # calibration frames too
    assert h['ORCHVER'] == V and 'ORCHMOD' not in h
    assert orchard_cards(h) == [INGEST]
    expected = ([ROT_CARD] if instrument == 'andor' else [MIR_CARD]) if imagetyp == 'Light Frame' else []
    assert other_cards(h) == expected                                       # the descriptive cards, as before
    assert verifies(str(path))


def test_ingest_transform_twice_changes_nothing(tmp_path):
    for name, imagetyp in (('light', 'Light Frame'), ('dark', 'Dark Frame')):
        eso_frame(tmp_path / (name + '.fits'), imagetyp=imagetyp)
    at.astra_transform(['light', 'dark'], str(tmp_path), 1)
    before = {n: md5(tmp_path / (n + '.fits')) for n in ('light', 'dark')}
    at.astra_transform(['light', 'dark'], str(tmp_path), 1)
    assert {n: md5(tmp_path / (n + '.fits')) for n in ('light', 'dark')} == before


def write_cube(path, n=3, with_bzero=True, seed=1):
    """A small SPIRIT datacube as create_datacubes.py writes it (METADATA table, new format), with ESO's checksums.
    Without BZERO it is one of the early cubes that hold counts - 32768 and no BZERO keyword."""
    counts = np.random.default_rng(seed).integers(200, 16000, size=(n,) + SPIRIT_SMALL)
    h = astra_header('Light Frame', SPIRIT_SMALL)
    stored = counts.astype(np.uint16) if with_bzero else (counts - 32768).astype(np.int16)
    names = ['SPECULOOS4.2026-10-16T02_00_{:02d}.000'.format(i) for i in range(n)]
    meta = fits.BinTableHDU.from_columns([
        fits.Column(name='FILENAME', format='40A', array=np.array(names)),
        fits.Column(name='JD-OBS', format='D', array=2461329.5 + np.arange(n) / 86400.0),
        fits.Column(name='DATE-OBS', format='23A', array=np.array(['2026-10-16T02:00:{:02d}.000'.format(i)
                                                                    for i in range(n)])),
        fits.Column(name='EXPTIME', format='D', array=np.full(n, 10.0)),
    ], name='METADATA')
    fits.HDUList([fits.PrimaryHDU(data=stored, header=h), meta]).writeto(str(path), checksum=True)
    return counts, ['SPECU4.2026-10-16T02:00:{:02d}.000'.format(i) for i in range(n)]


@pytest.mark.parametrize('with_bzero', [True, False])
def test_unpacked_frames_are_stamped_and_verify(tmp_path, with_bzero):
    cube = tmp_path / 'SPECU4.cube.fits'
    counts, names = write_cube(cube, with_bzero=with_bzero)
    out = tmp_path / 'frames'
    out.mkdir()
    ids, n = unpack_datacube(str(cube), str(out))
    assert sorted(ids) == names and n == 3 and not cube.exists()
    for i, name in enumerate(names):
        path = str(out / (name + '.fits'))
        with fits.open(path, do_not_scale_image_data=True) as hl:
            h = hl[0].header
            assert h['BITPIX'] == 16 and h['BZERO'] == 32768 and h['BSCALE'] == 1
            assert h['ORCHVER'] == V and 'ORCHMOD' not in h and orchard_cards(h) == [UNPACK]
            assert h['JD-OBS'] == pytest.approx(2461329.5 + i / 86400.0) and h['OLDCUBE'] == 'N'
        assert np.array_equal(fits.getdata(path).astype(np.int64), counts[i])
        assert verifies(path)
    # then the Astra transform, as the fetch and SSO_download run it on the unpacked frames
    at.astra_transform(names, str(out), 1)
    for i, name in enumerate(names):
        path = str(out / (name + '.fits'))
        h = fits.getheader(path)
        assert h['ORCHVER'] == V and h['ASTRAMIR'] == 'Y' and 'ORCHMOD' not in h
        assert orchard_cards(h) == [UNPACK, INGEST] and other_cards(h) == [MIR_CARD]
        assert np.array_equal(fits.getdata(path).astype(np.int64), np.flip(counts[i], axis=1))
        assert verifies(path)


def test_unpacking_a_cube_twice_gives_the_same_frames(tmp_path):
    for d in ('a', 'b'):
        (tmp_path / d).mkdir()
    _, names = write_cube(tmp_path / 'a' / 'cube.fits')
    shutil.copy(str(tmp_path / 'a' / 'cube.fits'), str(tmp_path / 'b' / 'cube.fits'))
    for d in ('a', 'b'):
        unpack_datacube(str(tmp_path / d / 'cube.fits'), str(tmp_path / d))
    for name in names:
        (ra, ha), (rb, hb) = (raw_wcs.read_header_bytes(str(tmp_path / d / (name + '.fits'))) for d in ('a', 'b'))
        assert open(tmp_path / 'a' / (name + '.fits'), 'rb').read()[len(ra):] == \
            open(tmp_path / 'b' / (name + '.fits'), 'rb').read()[len(rb):]
        # identical apart from when the checksums were computed (their comments, and so CHECKSUM's value)
        assert [c.image for c in ha.cards if c.keyword not in ('CHECKSUM', 'DATASUM')] == \
            [c.image for c in hb.cards if c.keyword not in ('CHECKSUM', 'DATASUM')]
        assert ha['DATASUM'] == hb['DATASUM']
        assert orchard_cards(ha) == [UNPACK]


def test_stamp_ingest_adds_its_card_once():
    h = fits.Header()
    h['ORCHMOD'] = 'x'
    rp.stamp_ingest(h, 'astra_transform')
    rp.stamp_ingest(h, 'astra_transform')
    assert orchard_cards(h) == [INGEST] and h['ORCHVER'] == V and 'ORCHMOD' not in h


# --- later edits ---------------------------------------------------------------------------------------------------

def test_command_line_transform_is_recorded_as_an_edit(tmp_path, monkeypatch):
    """A night fetched without the transform (SSO_download --skip-transform), turned later from the command line."""
    data = eso_frame(tmp_path / 'frame.fits')
    monkeypatch.setattr(sys, 'argv', ['astra_transform.py', str(tmp_path), 'rotate', '1', 'all'])
    at.main()
    path = str(tmp_path / 'frame.fits')
    h = fits.getheader(path)
    assert 'ORCHVER' not in h and h['ORCHMOD'] == V
    assert orchard_cards(h) == ['orchard {} astra_transform: rotate k=1'.format(V)] and other_cards(h) == [ROT_CARD]
    assert np.array_equal(fits.getdata(path), np.rot90(data)) and verifies(path)
    before = md5(path)
    at.main()                                                     # skipped: already carries ASTRAROT
    assert md5(path) == before


def test_a_forced_second_transform_keeps_orchver_and_adds_an_edit(tmp_path):
    eso_frame(tmp_path / 'frame.fits', SPIRIT_SMALL)
    at.astra_transform(['frame'], str(tmp_path), 1)
    at.astra_transform_unified(['frame'], transform_mode='mirror', dirname=str(tmp_path), force=True)
    h = fits.getheader(tmp_path / 'frame.fits')
    assert h['ORCHVER'] == V and h['ORCHMOD'] == V
    assert orchard_cards(h) == [INGEST, 'orchard {} astra_transform: mirror'.format(V)]
    assert other_cards(h) == [MIR_CARD, MIR_CARD]                 # repair_orientation counts these: still two
    assert verifies(str(tmp_path / 'frame.fits'))


@pytest.mark.parametrize('case', ['andor_k3', 'spirit_flip'])
def test_repair_orientation_records_the_edit_and_verifies(tmp_path, case):
    """The frames repair_orientation exists for, made the way they were made (a second, forced transform), and now
    stamped at ingest: the repair keeps ORCHVER, adds ORCHMOD and its card next to its own markers, and verifies."""
    night = tmp_path / 'Observations' / 'Ganymede' / 'images' / '20261016'
    night.mkdir(parents=True)
    andor = case == 'andor_k3'
    data = eso_frame(night / 'f.fits', ANDOR_NATIVE if andor else SPIRIT)
    at.astra_transform(['f'], str(night), 1)
    at.astra_transform_unified(['f'], transform_mode='rotate' if andor else 'mirror', dirname=str(night), force=True)
    path, op = str(night / 'f.fits'), (3 if andor else 'flip')
    assert ro.check(path, op) is None
    ro.repair(path, op, str(tmp_path / 'quarantine'), today='2026-10-19')
    with fits.open(path) as hl:
        h, out = hl[0].header, hl[0].data
    assert np.array_equal(out, np.rot90(data) if andor else np.flip(data, axis=1))
    assert h['ORCHVER'] == V and h['ORCHMOD'] == V
    assert orchard_cards(h) == [INGEST, 'orchard {} astra_transform: {}'.format(V, 'rotate k=1' if andor else 'mirror'),
                                'orchard {} repair_orientation: {}'.format(V, 'rot90 k=3' if andor else 'fliplr')]
    markers = [c for c in other_cards(h) if ro.REPAIR_CARD in c]
    assert markers == (['2026-10-19 orient repair: rot90 k=3; was rotated twice'] if andor else
                       ['2026-10-19 orient repair: fliplr; was mirrored twice']) + \
        ['2026-10-19 orient repair: celestial WCS removed; re-solve']
    assert verifies(path)
    assert ro.check(path, op) == 'refused: already repaired'     # a second run repairs nothing


def old_transform_header():
    """What the old astra_transform wrote for a Europa frame: rotated, CRPIX moved, Astra's PC left native (rot_pc)."""
    h = fits.Header()
    h['IMAGETYP'] = 'Light Frame'
    h['OBSERVER'] = 'Astra'
    h['CTYPE1'], h['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
    h['CRVAL1'], h['CRVAL2'] = 292.43961447225, 11.035372553426
    h['CRPIX1'], h['CRPIX2'] = 10.5, 15.5
    h['PC1_1'], h['PC1_2'] = -9.3859668511035e-06, -9.5989180007196e-05
    h['PC2_1'], h['PC2_2'] = 9.6034430595804e-05, -9.4006699007457e-06
    h['PA'] = (84.41224818505546, '[deg, 0-360 CCW] Position angle of plate')
    h['ASTRAROT'] = (90, 'Anticlockwise rotation applied to astra images')
    h.add_history(ROT_CARD)
    return h


def old_transformed_frame(path, checksums, filler=0):
    """A frame as the pre-ASTRAWCS astra_transform left it. checksums: 'valid' (as after T7 wrote a WCS into it
    through astropy's update mode), 'stale' (ESO's cards for the unturned image, as astra_transform kept them), or
    'none'."""
    h = old_transform_header()
    for i in range(filler):
        h['FILL{:04d}'.format(i)] = (i, 'padding to fill the header block')
    data = (np.arange(30 * 20, dtype=np.uint16).reshape(30, 20) * 97) % 65535
    if checksums == 'stale':
        eso = str(path) + '.eso'
        fits.PrimaryHDU(data=np.rot90(data, -1), header=astra_header()).writeto(eso, checksum=True)
        for k in ('CHECKSUM', 'DATASUM'):
            h[k] = fits.getheader(eso)[k]
        os.remove(eso)
    fits.PrimaryHDU(data=data, header=h).writeto(str(path), checksum=checksums == 'valid')
    os.utime(str(path), (1_600_000_000, 1_600_000_000))
    assert verifies(str(path)) == (checksums == 'valid')


@pytest.mark.parametrize('checksums', ['valid', 'stale', 'none', 'valid_overflow'])
def test_repair_astra_wcs_records_the_edit_and_leaves_checksums_valid(tmp_path, monkeypatch, checksums):
    path = str(tmp_path / 'frame.fits')
    filler = 0
    if checksums == 'valid_overflow':
        # fill the header block so that the cards the repair adds need another 2880-byte block
        old_transformed_frame(path, 'valid')
        n = len(raw_wcs.read_header_bytes(path)[1])
        os.remove(path)
        filler, checksums = (36 - (n + 1) % 36) % 36, 'valid'
    old_transformed_frame(path, checksums, filler)
    original = open(path, 'rb').read()
    old_len = len(raw_wcs.read_header_bytes(path)[0])
    old_data = original[old_len:]
    reads = []
    monkeypatch.setattr(raw_wcs, 'data_sum', lambda *a: reads.append(a) or rp.data_sum(*a))
    with tarfile.open(str(tmp_path / 'backup.tar'), 'a') as t:
        found = raw_wcs.apply_to_file(path, t)
    assert {k for k, _ in found} == {'rot_pc', 'stale_pa'}
    assert len(reads) == (0 if checksums == 'valid' else 1)      # the DATASUM is trusted only when CHECKSUM agrees
    new_raw, h = raw_wcs.read_header_bytes(path)
    assert open(path, 'rb').read()[len(new_raw):] == old_data    # pixels untouched, byte for byte
    assert (len(new_raw) > old_len) == (filler > 0)
    assert os.stat(path).st_mtime == 1_600_000_000
    assert h['ASTRAWCS'] == 2 and h['ORCHMOD'] == V and 'ORCHVER' not in h
    assert other_cards(h) == [ROT_CARD, 'repair_astra_wcs: rot_pc, stale_pa']    # the marker, as before
    assert orchard_cards(h) == ['orchard {} repair_astra_wcs: rot_pc, stale_pa'.format(V)]
    assert verifies(path)
    # a second run finds nothing to do and writes nothing
    after = open(path, 'rb').read()
    with tarfile.open(str(tmp_path / 'backup.tar'), 'a') as t:
        assert raw_wcs.apply_to_file(path, t) == []
    assert open(path, 'rb').read() == after
    assert raw_wcs.restore(str(tmp_path / 'backup.tar')) == 1 and open(path, 'rb').read() == original


# --- the repair classifiers ignore the provenance records ----------------------------------------------------------

def _rot_unrot():
    h = old_transform_header()
    del h['HISTORY']
    for k in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'):
        del h[k]
    m = np.array([[-9.6e-5, 1.2e-5], [-1.1e-5, -9.5e-5]])
    for (i, j), v in np.ndenumerate(m):
        h['CD{}_{}'.format(i + 1, j + 1)] = v
    h['PA'] = plate_pa(m)
    return h


def _ambiguous():
    h = old_transform_header()
    m = np.array([[h['PC1_1'], h['PC1_2']], [h['PC2_1'], h['PC2_2']]])
    h['CDELT1'] = h['CDELT2'] = 1.0
    h['PA'] = plate_pa(m)
    return h


def _mir_pc():
    h = old_transform_header()
    del h['ASTRAROT'], h['HISTORY']
    h['ASTRAMIR'] = 'Y'
    h.add_history(MIR_CARD)
    return h


def _unmir():
    h = _mir_pc()
    del h['HISTORY']
    for k, v in (('PC1_1', -1.1e-4), ('PC1_2', 0.0), ('PC2_1', 0.0), ('PC2_2', 1.1e-4)):
        h[k] = v
    return h


def _singular():
    h = old_transform_header()
    for k in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'):
        h[k] = 0.0
    return h


def _not_rotated():
    h = old_transform_header()
    h['NAXIS'], h['NAXIS1'], h['NAXIS2'] = 2, 2088, 2048
    return h


DIAGNOSE_CASES = {'rot_pc': old_transform_header, 'unrot': _rot_unrot, 'ambiguous': _ambiguous, 'mir_pc': _mir_pc,
                  'unmir': _unmir, 'singular': _singular, 'image_not_rotated': _not_rotated}


@pytest.mark.parametrize('case', sorted(DIAGNOSE_CASES))
def test_repair_astra_wcs_classifies_stamped_frames_as_before(case):
    bare = DIAGNOSE_CASES[case]()
    expected = raw_wcs.diagnose(bare)
    assert case in [k for k, _ in expected]
    stamped = bare.copy()
    rp.stamp_ingest(stamped, 'unpack_datacubes')
    rp.stamp_ingest(stamped, 'astra_transform')
    ov.record_edit(stamped, 'astra_transform', 'mirror')
    ov.record_edit(stamped, 'repair_orientation', 'fliplr')
    assert raw_wcs.diagnose(stamped) == expected


def test_new_ingest_output_needs_no_wcs_repair(tmp_path):
    eso_frame(tmp_path / 'frame.fits')
    at.astra_transform(['frame'], str(tmp_path), 1)
    assert raw_wcs.diagnose(fits.getheader(tmp_path / 'frame.fits')) == []


def test_repair_orientation_classifies_stamped_frames_as_before(tmp_path):
    """check() for every operation, on frames made by the stamping code, against the same frames without the
    records: mirrored once, mirrored twice, never mirrored (unpacked, transform skipped), rotated once and twice."""
    frames = {}
    for name, shape, imagetyp in (('spirit_once', SPIRIT, 'Light Frame'), ('spirit_twice', SPIRIT, 'Light Frame'),
                                  ('spirit_dark', SPIRIT, 'Dark Frame'), ('andor_once', ANDOR_NATIVE, 'Light Frame'),
                                  ('andor_twice', ANDOR_NATIVE, 'Light Frame')):
        eso_frame(tmp_path / (name + '.fits'), shape, imagetyp)
        at.astra_transform([name], str(tmp_path), 1)
        if name.endswith('twice'):
            at.astra_transform_unified([name], transform_mode='rotate' if name.startswith('andor') else 'mirror',
                                       dirname=str(tmp_path), force=True)
        frames[name] = str(tmp_path / (name + '.fits'))
    never = str(tmp_path / 'never.fits')
    fits.PrimaryHDU(data=np.zeros(SPIRIT, np.uint16),
                    header=rp.stamp_ingest(astra_header('Light Frame', SPIRIT), 'unpack_datacubes')).writeto(never)
    frames['spirit_never'] = never
    accepted = {}
    for name, path in frames.items():
        bare = strip_provenance(path)
        assert orchard_cards(fits.getheader(path)) and not orchard_cards(fits.getheader(bare))
        results = {op: ro.check(path, op) for op in (1, 3, 'flip', 'mirror')}
        assert results == {op: ro.check(bare, op) for op in (1, 3, 'flip', 'mirror')}, name
        accepted[name] = sorted(str(op) for op, why in results.items() if why is None)
    # a calibration frame has no mirroring card to count, so check() accepts a flip of any mirrored one
    assert accepted == {'spirit_once': [], 'spirit_twice': ['flip'], 'spirit_dark': ['flip'], 'andor_once': [],
                        'andor_twice': ['1', '3'], 'spirit_never': ['mirror']}


# --- the checksum helpers ------------------------------------------------------------------------------------------

def test_hand_made_checksums_agree_with_astropy(tmp_path):
    rng = np.random.default_rng(3)
    path = str(tmp_path / 'c.fits')
    for t in range(40):
        data = rng.integers(-32768, 32767, size=(rng.integers(1, 30), rng.integers(1, 30))).astype(np.int16)
        h = fits.Header()
        h['VALUE'] = int(rng.integers(0, 10 ** 9))
        h['TEXT'] = 'case {}'.format(t)
        fits.PrimaryHDU(data=data, header=h).writeto(path, overwrite=True, checksum=True)
        raw, header = raw_wcs.read_header_bytes(path)
        assert rp.data_sum(path, len(raw), header) == int(header['DATASUM'])
        assert rp.claimed_datasum(raw, header) == int(header['DATASUM'])
        blank = header.copy()
        blank['CHECKSUM'] = '0' * 16
        assert rp.encode_checksum(rp.ones_sum(blank.tostring(padding=True).encode('ascii'),
                                              int(header['DATASUM']))) == header['CHECKSUM']
        edited = header.copy()
        edited['NEW'] = t
        assert rp.claimed_datasum(edited.tostring(padding=True).encode('ascii'), edited) is None
        rp.set_checksum(edited, int(header['DATASUM']))
        with open(path, 'rb') as f:
            body = f.read()[len(raw):]
        with open(path, 'wb') as f:
            f.write(edited.tostring(padding=True).encode('ascii') + body)
        assert verifies(path)
