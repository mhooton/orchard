"""
The pipeline's readers of frames that may be stored Rice-compressed (.fz, utils/framepack.py), and ZLP_pipeline.sh's
handling of a compressed raw night (docs/frame-compression.md).

Covers:
  - pipeutils.open_fits_file / proc_name on .fz frames (T2-T6, smear_check), and a plain list entry followed to its .fz
  - T1's get_liste lists compressed frames once (an earlier night's calibrations), plain twin first
  - stages/restage: procX.fits.fz is the same product as procX.fits; restage reads its cards from HDU 1
  - the T11 report stacks proc*.fits.fz, and names aligned copies as for plain frames
  - telescope_faults sees compressed raw frames, and a plain frame wins over a .fz of the same frame
  - the job queue's header reader and suffixes (Artemis .fts.fz)
  - SSO_download's transformation check keeps a compressed frame
  - ZLP: find_ext never picks fz over plain frames; unpack_raw_night / repack_raw_night round trip on a real night
    layout; unpack_frames for a rerun from T7-T9

Run:  pytest src/tests/test_compressed_frame_readers.py -q
"""
import importlib.util
import os
import re
import subprocess
import sys
import types
import warnings

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import framepack as fp      # noqa: E402
from utils import rawpack as rp        # noqa: E402
from utils import restage, stages      # noqa: E402

warnings.simplefilter('ignore', category=fits.verify.VerifyWarning)
ZLP = os.path.join(SRC, 'main', 'ZLP_pipeline.sh')


def _raw(path, seed=0, **cards):
    rng = np.random.default_rng(seed)
    data = np.clip(1300 + 30 * rng.standard_normal((40, 48)), 0, 65535).astype(np.uint16)
    h = fits.Header()
    h['TELESCOP'] = 'Io'
    h['IMAGETYP'] = 'Light Frame'
    h['FILTER'] = 'I+z'
    h['OBJECT'] = 'Sp0000+0000'
    h['EXPTIME'] = 20.0
    h['DATE-OBS'] = '2026-10-16T23:08:07.000'
    for k, v in cards.items():
        h[k] = v
    fits.PrimaryHDU(data=data, header=h).writeto(path, checksum=True)
    return data


def _packed(path, seed=0, **cards):
    data = _raw(path, seed, **cards)
    assert fp.convert_in_place(path, 'rawpack').status == 'done'
    return data


# --------------------------------------------------------------------------- T1-T6

def test_open_fits_file_reads_a_compressed_raw_frame_and_follows_a_plain_name(tmp_path):
    from calibration import pipeutils
    p = str(tmp_path / 'Observations' / 'Io' / 'images' / '20261016' / 'SPECU1.2026-10-16T23:08:07.000.fits')
    os.makedirs(os.path.dirname(p))
    data = _packed(p)
    for name in (p, p + '.fz'):
        with pipeutils.open_fits_file(name) as h:
            assert h[0].data.dtype == np.uint16 and np.array_equal(h[0].data, data)
            assert h[0].header['IMAGETYP'] == 'Light Frame'
            assert pipeutils.frame_path(h) == p + '.fz'
    assert pipeutils.proc_name(p + '.fz') == 'procSPECU1.2026-10-16T23:08:07.000.fits'
    assert pipeutils.proc_name('/x/Sp1-S001-R001-C001-i.fts.fz') == 'procSp1-S001-R001-C001-i.fits'


def test_get_liste_lists_compressed_frames_once(tmp_path):
    pytest.importorskip('tabulate')            # createlists imports reporting.QC
    from main import createlists
    night = tmp_path / '20261015'
    (night / 'Calibration').mkdir(parents=True)
    a, b, c = str(night / 'a.fits'), str(night / 'b.fits'), str(night / 'Calibration' / 'bias.fits')
    _raw(a)
    _packed(b)
    _packed(c)
    _raw(c)                                    # a working copy beside its .fz: counted once, as the plain file
    got = createlists.get_liste(str(night), None, 'fits')
    assert sorted(got) == sorted([a, b + '.fz', c])
    assert createlists.sort_flatlist([b + '.fz'])[0] == ['I+z']


# --------------------------------------------------------------------------- stages, restage

@pytest.mark.parametrize('rel, stage', [('1/procX.fits', 6), ('1/procX.fits.fz', 6), ('1/procX.fits.phot', 9)])
def test_a_compressed_frame_is_the_same_product(rel, stage):
    assert stages.writer(rel) == stage
    assert stages.made_before(rel, 10) and not stages.made_before(rel, 1)


def test_restage_reads_the_cards_of_a_compressed_frame(tmp_path):
    p = str(tmp_path / 'procX.fits')
    _raw(p, ORCHVER='4.3.4', ORCHREL='20261021T104419Z-7e4521e')
    fp.convert_in_place(p, 'compress')
    assert p.endswith(restage.FITS[:3]) and (p + '.fz').lower().endswith(restage.FITS)
    cards = restage.fits_cards(p + '.fz')
    assert cards['orchver'] == '4.3.4' and cards['orchrel'] == '20261021T104419Z-7e4521e' and 'orchmod' in cards


# --------------------------------------------------------------------------- T11

def _report():
    for name, attrs in (('fitsio', {'FITS': object}), ('skimage', {}),
                        ('skimage.registration', {'phase_cross_correlation': None})):
        try:
            __import__(name)
        except ImportError:
            mod = types.ModuleType(name)
            for k, v in attrs.items():
                setattr(mod, k, v)
            sys.modules[name] = mod
    spec = importlib.util.spec_from_file_location('pdf_report_for_fz',
                                                  os.path.join(SRC, 'reporting', 'pdf_report_catriona.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_report_stacks_compressed_frames(tmp_path):
    rep = _report()
    one = tmp_path / '1'
    one.mkdir()
    a, b, c = str(one / 'procA.fits'), str(one / 'procB.fits'), str(one / 'procC.fits')
    _raw(a, **{'JD-OBS': 2461330.1})
    _packed(b, **{'JD-OBS': 2461330.2})
    _packed(c, **{'JD-OBS': 2461330.3})
    _raw(c, **{'JD-OBS': 2461330.3})           # a rerun's new frame beside the older .fz
    frames = rep.processed_frames(str(one) + '/')
    assert frames == [a, b + '.fz', c]
    assert rep.frame_jds(frames)[b + '.fz'] == 2461330.2
    assert rep.aligned_name('/t', b + '.fz') == rep.aligned_name('/t', b) == '/t/procB_aligned.fits'
    with rep.open_frame(b + '.fz') as h:
        assert h[0].data.shape == (40, 48)


# --------------------------------------------------------------------------- telescope faults, queue, download

def test_telescope_faults_sees_compressed_raw_and_prefers_plain_processed(tmp_path):
    from utils import telescope_faults as tf
    raw = tmp_path / 'images' / '20261016'
    raw.mkdir(parents=True)
    _packed(str(raw / 'SPECU1.2026-10-16T23:08:07.000.fits'))
    _raw(str(raw / 'SPECU1.2026-10-16T23:09:07.000.fits'), DATE_OBS='x')
    frames = tf.raw_frames(str(raw))
    assert sorted(f.name for f in frames) == ['SPECU1.2026-10-16T23:08:07.000.fits.fz',
                                              'SPECU1.2026-10-16T23:09:07.000.fits']
    assert frames[0].light and frames[0].object == 'Sp0000+0000' and frames[0].exptime == 20.0
    out = tmp_path / 'output' / '20261016' / 'Sp0000+0000' / '1'
    out.mkdir(parents=True)
    old = str(out / 'procSPECU1.2026-10-16T23:08:07.000.fits')
    _packed(old)
    _raw(old)
    assert tf.processed_frames([str(tmp_path / 'output' / '20261016')]) == {'SPECU1.2026-10-16T23:08:07.000': old}
    assert tf.frame_key('procSp1-S001-R001-C001-i.fts.fz') == 'Sp1-S001-R001-C001-i'


def test_the_queue_reads_compressed_headers_and_artemis_suffixes(tmp_path):
    from jobqueue import nights
    p = str(tmp_path / 'Sp1056+0700-S001-R001-C001-r.fts')
    _packed(p)
    assert nights.list_frames(str(tmp_path)) == ['Sp1056+0700-S001-R001-C001-r.fts.fz']
    h = nights.read_header(p + '.fz')
    assert h['IMAGETYP'] == 'Light Frame' and h['NAXIS1'] == 48 and h['NAXIS2'] == 40
    assert nights.read_header(p + '.fz', max_blocks=12)['DATE-OBS'] == '2026-10-16T23:08:07.000'


def test_sso_download_keeps_a_compressed_frame(tmp_path):
    pytest.importorskip('scipy')
    from download.astra_transform import transformation_check
    _packed(str(tmp_path / 'SPECU4.2026-10-16T23:08:07.000.fits'))
    assert transformation_check(str(tmp_path), ['SPECU4.2026-10-16T23:08:07.000'], 'SPIRIT') == \
        ['SPECU4.2026-10-16T23:08:07.000']
    assert os.path.exists(str(tmp_path / 'SPECU4.2026-10-16T23:08:07.000.fits.fz'))


# --------------------------------------------------------------------------- ZLP_pipeline.sh

def _function(name):
    text = open(ZLP).read()
    start = text.index(name + '(){')
    return text[start:text.index('\n}\n', start) + 3]


def _bash(script, env=None):
    full = dict(os.environ)
    full.update(env or {})
    return subprocess.run(['bash', '-c', 'set -o nounset; set -o pipefail\n' + script], env=full,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)


def test_the_default_frame_format_is_rice():
    assert re.search(r'^FRAME_FORMAT="rice"', open(ZLP).read(), re.M)


def test_find_ext_takes_plain_frames_over_a_leftover_fz(tmp_path):
    night = tmp_path / '20261016'
    night.mkdir()
    for name in ('SPECU1.a.fits', 'SPECU1.b.fits', 'SPECU1.z.fits.fz', 'Data_Download.txt'):
        (night / name).write_bytes(b'')
    out = _bash(_function('find_ext') + 'IMGDIRS={}; find_ext'.format(night))
    assert out.stdout.strip() == 'fits', out.stdout
    for name in ('SPECU1.a.fits', 'SPECU1.b.fits'):
        (night / name).unlink()
    assert _bash(_function('find_ext') + 'IMGDIRS={}; find_ext'.format(night)).stdout.strip() == 'fz'


def test_zlp_unpacks_a_compressed_night_for_the_run_and_repacks_it(tmp_path):
    night = tmp_path / 'images' / '20261016'
    night.mkdir(parents=True)
    originals = {}
    for k in range(3):
        p = str(night / 'SPECU1.2026-10-16T23:0{}:00.000.fits'.format(k))
        _raw(p, k)
        originals[p] = open(p, 'rb').read()
        fp.convert_in_place(p, 'rawpack')
    funcs = _function('unpack_raw_night') + _function('repack_raw_night')
    common = 'SCRIPTDIR={} CORES=2 IMGDIRS={} T1=1 T2=1 T3=1 T4=1 T5=1 T6=1 T7=1\n'.format(SRC, night)
    env = {'PYTHONPATH': SRC}

    out = _bash(funcs + common + 'unpack_raw_night; ls {}'.format(night), env)
    assert out.returncode == 0, out.stdout
    assert 'rawpack unpack exit code 0' in out.stdout
    for p, raw in originals.items():
        assert open(p, 'rb').read() == raw

    out = _bash(funcs + common + 'repack_raw_night {}'.format(night), env)
    assert out.returncode == 0 and 'rawpack repack exit code 0' in out.stdout, out.stdout
    assert rp.status(str(night)) == {'compressed': 3, 'unpacked': 0, 'plain': 0}

    # a run without T1-T7 (from T10) leaves the night compressed; a night with nothing compressed is left alone
    out = _bash(funcs + common.replace('T7=1', 'T7=0').replace('T1=1', 'T1=0').replace('T2=1', 'T2=0')
                .replace('T3=1', 'T3=0').replace('T4=1', 'T4=0').replace('T5=1', 'T5=0').replace('T6=1', 'T6=0')
                + 'unpack_raw_night', env)
    assert 'START raw unpack' not in out.stdout
    assert rp.status(str(night))['compressed'] == 3


def test_zlp_repacks_when_the_run_is_terminated(tmp_path):
    night = tmp_path / 'images' / '20261016'
    night.mkdir(parents=True)
    p = str(night / 'SPECU1.2026-10-16T23:00:00.000.fits')
    _packed(p)
    text = open(ZLP).read()
    traps = '\n'.join(line for line in text.splitlines() if line.startswith(('RAW_NIGHT=""', 'trap ')))
    assert traps.count('trap ') == 3
    script = (_function('unpack_raw_night') + _function('repack_raw_night') + traps + '\n'
              + 'SCRIPTDIR={} CORES=1 IMGDIRS={} T1=1 T2=0 T3=0 T4=0 T5=0 T6=0 T7=0\n'.format(SRC, night)
              + 'RAW_NIGHT=$IMGDIRS; unpack_raw_night; kill -TERM $$; sleep 5; echo not reached')
    out = _bash(script, {'PYTHONPATH': SRC})
    assert out.returncode == 143 and 'not reached' not in out.stdout, out.stdout
    assert rp.status(str(night)) == {'compressed': 1, 'unpacked': 0, 'plain': 0}


def test_zlp_decompresses_packed_frames_for_a_rerun_from_t7_to_t9(tmp_path):
    one = tmp_path / 'output' / '20261016' / 'Sp0000+0000' / '1'
    one.mkdir(parents=True)
    p = str(one / 'procX.fits')
    _raw(p)
    fp.convert_in_place(p, 'compress')
    logs = tmp_path / 'logs'
    logs.mkdir()
    common = 'SCRIPTDIR={} CORES=1 DATDIR={} DATE=20261016 RUNNAME=1 '.format(SRC, tmp_path)
    out = _bash(_function('unpack_frames') + common + 'T6=1 T7=1 T8=1 T9=1; unpack_frames {}'.format(one),
                {'PYTHONPATH': SRC})
    assert 'START frame unpacking' not in out.stdout and os.path.exists(p + '.fz')     # T6 remakes them
    out = _bash(_function('unpack_frames') + common + 'T6=0 T7=0 T8=0 T9=1; unpack_frames {}'.format(one),
                {'PYTHONPATH': SRC})
    assert 'framepack decompress exit code 0' in out.stdout, out.stdout
    assert os.path.exists(p) and not os.path.exists(p + '.fz')
    assert (logs / 'framepack_20261016_1.csv').exists()
