"""Regression tests for the nightly PDF report (T11).

Covers the two faults that stopped reports being produced in production:

  * a night with no reduced targets made page 3 ask GridSpec for -1 columns,
    which aborted the report and - because ZLP_pipeline.sh runs under errexit -
    the stages after it as well;
  * align() claimed 20 worker processes regardless of the pipeline's core
    budget, and each worker let OpenBLAS size its own thread pool from the host
    core count, so concurrent nights exhausted RLIMIT_NPROC and wedged.

Run with `pytest src/tests/test_pdf_report.py`, or directly with python.
"""

import multiprocessing
import os
import re
import subprocess
import sys
import types

import pytest

# PDF_REPORT_PATH lets these run against a deployed copy, e.g. the one inside
# the orchard-server container at /opt/orchard/src/reporting/.
REPORT = os.environ.get(
    'PDF_REPORT_PATH',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'reporting', 'pdf_report_catriona.py'))


def _stub(name, **attrs):
    """Minimal stand-in so the module imports without the pipeline's C deps."""
    mod = types.ModuleType(name)
    for key, value in attrs.items():
        setattr(mod, key, value)
    return mod


def _load_report_module():
    import importlib.util

    for name, attrs in (('fitsio', {'FITS': object}),
                        ('skimage', {}),
                        ('skimage.registration', {'phase_cross_correlation': None})):
        try:
            __import__(name)
        except ImportError:
            sys.modules[name] = _stub(name, **attrs)

    spec = importlib.util.spec_from_file_location('pdf_report_catriona', REPORT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rep = _load_report_module()


class _Args(object):
    def __init__(self, datdir, target='Sp0000+0000'):
        self.datdir = datdir
        self.obsdir = datdir
        self.date = '20260924'
        self.telescope = 'Io'
        self.target = target
        self.ap = '5'
        self.version = 'v3'
        self.nproc = 1
        self.align_timeout = 60.0


def _synthetic_inputs(ntargs):
    """The arguments make_pdf would get for a night with ntargs reduced targets.

    Every image path is deliberately missing, so the report falls back to its
    'No plot' panels and the test exercises layout rather than plotting.
    """
    per_target = lambda prefix: ['%s_%d.png' % (prefix, i) for i in range(ntargs)]
    imgs = (
        ['boltwood.png', 'asm.png'],                       # 0  weather
        'masterbias.png',                                  # 1
        'masterdark.png',                                  # 2
        ['masterflat_I+z.png'],                            # 3
        ['dlc_%d.png' % i for i in range(6 * ntargs)],     # 4  six apertures each
        per_target('mlc'),                                 # 5
        per_target('complc'),                              # 6
        per_target('rmsflux'),                             # 7
        per_target('ramove'),                              # 8
        per_target('decmove'),                             # 9
        per_target('airmass'),                             # 10
        per_target('fwhm'),                                # 11
        per_target('skybkg'),                              # 12
        per_target('altitude'),                            # 13
        per_target('haoff'),                               # 14
        per_target('decoff'),                              # 15
        per_target('stack'),                               # 16
        per_target('watervapor'),                          # 17
        per_target('dx'),                                  # 18
        per_target('dy'),                                  # 19
    )
    vals = [['unknown'] * 5 for _ in range(ntargs)]
    jdstart = [2460000.5 + i for i in range(ntargs)]
    jdend = [2460000.6 + i for i in range(ntargs)]
    nimages = [100 + i for i in range(ntargs)]
    newap = ['5'] * ntargs
    ellip = [0.2] * ntargs
    targs = ['Sp%04d+0000' % i for i in range(ntargs)]
    gaia = ['%019d' % (i + 1) for i in range(ntargs)]
    flags = ['N', 'N', 'OK', 'N', 'N']
    return imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia, flags


def _page_count(path):
    with open(path, 'rb') as handle:
        data = handle.read()
    assert data.rstrip().endswith(b'%%EOF'), '%s was never closed' % path
    return len(re.findall(br'/Type\s*/Page[^s]', data))


# --------------------------------------------------------------------- fault 1

@pytest.mark.parametrize('ntargs', [0, 1, 2, 3])
def test_make_pdf_survives_any_number_of_targets(tmp_path, ntargs):
    """Zero reduced targets used to raise ValueError from GridSpec(ncols=-1)."""
    datdir = str(tmp_path)
    os.makedirs(os.path.join(datdir, 'reports', 'temp'))
    outname = os.path.join(datdir, 'reports', 'Io_20260924.pdf')

    imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia, flags = _synthetic_inputs(ntargs)
    order = rep.get_chronological_order(jdstart)

    rep.make_pdf(imgs, outname, _Args(datdir), [], vals, jdstart, jdend, nimages,
                 newap, flags, ellip, order, targs, gaia)

    # Two fixed pages, the target overview, then three detail pages per target.
    assert _page_count(outname) == 3 + 3 * ntargs


def test_import_night_img_keeps_per_target_lists_aligned(tmp_path):
    """A target with no *_output.fits must not push jdstart out of step.

    It used to append a placeholder to jdstart but not to targs/gaia, so
    argsort(jdstart) indexed past the end of gaia and page 3 - along with every
    detail page - was silently dropped by the catch-all around it.
    """
    datdir = str(tmp_path)
    for name in ('Sp1300+1912', 'Sp1735-2745'):
        os.makedirs(os.path.join(datdir, 'output', '20260924', name))

    args = _Args(datdir, target='Sp1300+1912 Sp1735-2745')
    (ramove, decmove, airmass, fwhm, skybkg, altitude, jdstart, jdend, nimages,
     filt, ellip, target_x_pos, target_y_pos, targs, gaia) = rep.import_night_img(args)

    for name, values in (('jdstart', jdstart), ('jdend', jdend), ('nimages', nimages),
                         ('filt', filt), ('ellip', ellip), ('target_x_pos', target_x_pos),
                         ('target_y_pos', target_y_pos), ('ramove', ramove),
                         ('decmove', decmove), ('airmass', airmass), ('fwhm', fwhm),
                         ('skybkg', skybkg), ('altitude', altitude)):
        assert len(values) == len(targs), '%s is out of step with targs' % name
    assert targs == [] and gaia == []
    assert list(rep.get_chronological_order(jdstart)) == []


def test_import_water_vapor_handles_no_targets(tmp_path):
    """jdstart[0] was read unconditionally once an LHATPRO log was present.

    It only stayed hidden because the old placeholder left a 0 sitting in
    jdstart; with the lists kept aligned there is nothing to index.
    """
    datdir = str(tmp_path)
    logdir = os.path.join(datdir, 'technical_logs', '2026-09-24')
    os.makedirs(logdir)
    with open(os.path.join(logdir, 'Io_ESO_LHATPRO_log.txt'), 'w') as handle:
        handle.write('header\n')
        handle.write(','.join(['x'] * 8) + '\n')

    assert rep.import_water_vapor(_Args(datdir), [], [], []) == []


# --------------------------------------------------------------------- fault 2

def test_thread_pools_are_capped_at_import():
    for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        assert os.environ.get(var) == '1', '%s was left unset' % var


class _RecordingPool(object):
    """Stands in for multiprocessing.Pool; records the size it was asked for."""

    sizes = []

    def __init__(self, size, timeout_error=False):
        _RecordingPool.sizes.append(size)
        self._timeout_error = timeout_error
        self.terminated = False

    def map_async(self, func, iterable):
        pool = self

        class _Result(object):
            def get(self, timeout=None):
                if pool._timeout_error:
                    raise multiprocessing.TimeoutError()
                return [None for _ in iterable]

        return _Result()

    def close(self):
        pass

    def terminate(self):
        self.terminated = True

    def join(self):
        pass


@pytest.fixture
def fake_pool(monkeypatch):
    _RecordingPool.sizes = []
    pools = []

    def factory(size):
        pool = _RecordingPool(size, timeout_error=factory.timeout_error)
        pools.append(pool)
        return pool

    factory.timeout_error = False
    monkeypatch.setattr(multiprocessing, 'Pool', factory)
    return factory, pools


@pytest.mark.parametrize('nproc,nimages,expected', [
    (1, 50, 1),      # the default: one worker, whatever the host looks like
    (20, 50, 20),    # the pipeline's core budget is honoured
    (20, 4, 4),      # never more workers than images
    (0, 50, 1),      # nonsense input still yields a usable pool
])
def test_align_pool_size_follows_nproc(fake_pool, tmp_path, nproc, nimages, expected):
    _factory, _pools = fake_pool
    liste = [str(tmp_path / ('proc%03d.fits' % i)) for i in range(nimages)]
    rep.align(liste, str(tmp_path), nproc=nproc, timeout=60)
    assert _RecordingPool.sizes == [min(expected, multiprocessing.cpu_count())]


def test_align_gives_up_instead_of_hanging(fake_pool, tmp_path):
    factory, pools = fake_pool
    factory.timeout_error = True
    liste = [str(tmp_path / ('proc%03d.fits' % i)) for i in range(4)]

    with pytest.raises(RuntimeError, match='timed out'):
        rep.align(liste, str(tmp_path), nproc=2, timeout=0.1)
    assert pools[0].terminated, 'the wedged pool was left running'


def test_stack_failure_costs_only_the_stack(monkeypatch, tmp_path):
    """A failed alignment must not take the whole report with it."""
    datdir = str(tmp_path)
    procdir = os.path.join(datdir, 'output', '20260924', 'Sp1300+1912', '1')
    os.makedirs(procdir)
    open(os.path.join(procdir, 'proc001.fits'), 'w').close()

    def boom(*a, **kw):
        raise RuntimeError('image alignment timed out after 0.1 s')

    monkeypatch.setattr(rep, 'align', boom)
    stacks, vals, dx, dy = rep.import_stack_and_vals(
        _Args(datdir), ['unknown'], ['unknown'], ['12345'], ['Sp1300+1912'])

    assert stacks == ['fake_path']
    # overscan, RON, dark current, CCD temperature, PA, gain, readout mode
    assert vals == [['unknown'] * 7]


# ------------------------------------------------- readout mode (readout.py)

def test_stack_and_vals_reads_gain_and_readout_mode(tmp_path):
    from astropy.io import fits
    datdir = str(tmp_path)
    os.makedirs(os.path.join(datdir, 'reports', 'temp'))
    liste = []
    for i, mode in enumerate(['GAIN3.8', 'GAIN3.8', 'PED2900']):
        path = os.path.join(datdir, 'proc%03d.fits' % i)
        hdu = fits.PrimaryHDU([[1.0, 2.0], [3.0, 4.0]])
        hdu.header.update({'RON': 9.3, 'DARKCUR': 0.27, 'GAIN': 3.83, 'READMODE': mode,
                           'OVERSCAN': 300.0, 'PA': 0.5, 'CCD-TEMP': -60.0, 'DX': 0.1, 'DY': 0.1,
                           'JD-OBS': 2461000.5 + i})
        hdu.writeto(path)
        liste.append(path)
    out = os.path.join(datdir, 'reports', 'temp', 'x.fits')
    _, vals = rep.stack_and_vals(_Args(datdir), liste, out, out + '_save', datdir, out + 'dx.png',
                                 out + 'dy.png', [1.0], [1.0], 0)
    assert vals[1] == 9.3 and vals[5] == 3.83
    assert vals[6] == 'GAIN3.8+PED2900'


@pytest.mark.parametrize('ron, gain, modes, expected', [
    (6.3, 1.08, ['NORMAL'], (True, True)),
    (9.3, 1.03, ['NORMAL'], (True, False)),          # a noisy normal night
    (9.3, 3.83, ['GAIN3.8'], (False, True)),         # the mode, not the noise
    (16.0, 3.83, ['GAIN3.8'], (False, False)),
    (6.3, None, [], (None, True)),                   # frames from before readout.py
])
def test_readout_flags(ron, gain, modes, expected):
    if 'GAIN3.8' in modes and not rep.READOUT_MODES:
        pytest.skip('calibration.readout is not importable here')
    assert rep.readout_flags(ron, gain, modes) == expected


def test_make_pdf_with_readout_values(tmp_path):
    datdir = str(tmp_path)
    os.makedirs(os.path.join(datdir, 'reports', 'temp'))
    outname = os.path.join(datdir, 'reports', 'Io_20260924.pdf')
    imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia, flags = _synthetic_inputs(2)
    vals = [[300.0, 9.3, 0.27, -60.0, 0.4, 3.83, 'GAIN3.8'], ['unknown'] * 5]
    order = rep.get_chronological_order(jdstart)
    rep.make_pdf(imgs, outname, _Args(datdir), [], vals, jdstart, jdend, nimages,
                 newap, flags, ellip, order, targs, gaia)
    assert _page_count(outname) == 3 + 3 * 2
    if rep.orchard_version:                                 # the version is on the first page and in the metadata
        assert rep.version_line().startswith('orchard ')
        with open(outname, 'rb') as f:
            assert b'/Creator (orchard ' in f.read()


# ------------------------------------------------ a space in a target's name

# A space in a target's name is '--' in its directory (createlists.py).  On Ganymede 20250716 ZLP turned
# TOI-4336--c back into 'TOI-4336 c' and passed '--target "Sp2135-3153 TOI-4336 c"'; the report split that
# into three targets, read neither 'TOI-4336' nor 'c', and its three stars with photometry were missing from
# every page after the first.

ZLP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'main', 'ZLP_pipeline.sh')


def _pdf_report_argv(target):
    """The report's command line from ZLP_pipeline.sh pdf_report(), for a night whose TARGET is target."""
    text = open(ZLP).read()
    m = re.search(r'^pdf_report\(\)\{\n.*?^\}\n', text, re.S | re.M)
    assert m, 'pdf_report not found in ZLP_pipeline.sh'
    script = ('set -o nounset\n'
              'SCRIPTDIR=/s DATDIR=/d OBSDIR=/o DATE=20250716 TEL=Ganymede VERSION=v3 CORES=1 PDF_REPORT_TIMEOUT=1\n'
              'timeout(){ printf "%%s\\n" "$@"; }\n'      # stands in for the report: prints its arguments
              '%s'
              'TARGET="%s"\n'                             # as find_target_names leaves it: one string
              'pdf_report\n') % (m.group(0), target)
    out = subprocess.run(['bash', '-c', script], stdout=subprocess.PIPE, check=True,
                         universal_newlines=True).stdout
    return out.splitlines()


def _report_target(target):
    """The --target ZLP gives the report for a night whose TARGET is target."""
    argv = _pdf_report_argv(target)
    return argv[argv.index('--target') + 1]


def test_zlp_passes_the_report_the_target_directories():
    target = _report_target('Sp2135-3153 TOI-4336--c')
    assert target == 'Sp2135-3153 TOI-4336--c'
    assert rep.target_dirs(_Args('/d', target=target)) == ['Sp2135-3153', 'TOI-4336--c']


def test_target_dirs_lists_each_directory_once():
    assert rep.target_dirs(_Args('/d', target='TOI-4336 TOI-4336--c  TOI-4336')) == ['TOI-4336', 'TOI-4336--c']
    assert rep.display_name('TOI-4336--c') == 'TOI-4336 c'
    assert rep.display_name('G--261-6--b') == 'G 261-6 b'


class _FitsioRead(object):
    def __init__(self, data):
        self.data = data

    def read(self):
        return self.data

    def __getitem__(self, column):
        return _FitsioRead(self.data[column])


class _AstropyAsFitsio(object):
    """Enough of fitsio.FITS for import_night_img, where fitsio is not installed."""

    def __init__(self, path):
        from astropy.io import fits
        with fits.open(path, memmap=False) as hdul:
            self.hdus = {h.name.lower(): h.data for h in hdul}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __getitem__(self, name):
        return _FitsioRead(self.hdus[name])


def _write_output_fits(path, nframes=4):
    """The HDUs of a <gaia>_<filter>_<night>_output.fits that import_night_img reads."""
    import numpy as np
    from astropy.io import fits
    cols = [fits.Column(name=name, format='D', array=np.full(nframes, value))
            for name, value in (('airmass', 1.2), ('altitude', 60.0), ('ra_move', 0.1), ('dec_move', -0.1),
                                ('fwhm', 3.0), ('skylevel', 500.0), ('psf_a_5', 2.0), ('psf_b_5', 1.8),
                                ('psf_t_5', 10.0))]
    cols.append(fits.Column(name='filter', format='8A', array=['r'] * nframes))
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns(cols, name='IMAGELIST'),
        fits.BinTableHDU.from_columns([fits.Column(name='obj_id', format='26A', array=['SP000001'])],
                                      name='CATALOGUE'),
        fits.ImageHDU(2460873.6 + 0.01 * np.arange(nframes)[None, :], name='JD'),
        fits.ImageHDU(np.ones((1, nframes)), name='CCDX'),
        fits.ImageHDU(np.ones((1, nframes)), name='CCDY'),
    ]).writeto(path)


def test_import_night_img_reads_a_target_whose_name_has_a_space(monkeypatch, tmp_path):
    datdir = str(tmp_path)
    night = os.path.join(datdir, 'output', '20260924')
    os.makedirs(os.path.join(night, 'TOI-4336', '1'))      # an earlier run's, with no photometry, as that night
    stars = {'Sp2135-3153': ['6784738092073342848'],
             'TOI-4336--c': ['6113245033656232448', '6113245033659187200']}
    for target, gaia_ids in stars.items():
        os.makedirs(os.path.join(night, target))
        for gaia_id in gaia_ids:
            _write_output_fits(os.path.join(night, target, '%s_r_20260924_output.fits' % gaia_id))
    os.makedirs(os.path.join(datdir, 'reports', 'temp', '20260924'))
    os.makedirs(os.path.join(datdir, 'technical_logs'))     # the PSF plots go in technical_logs/<day>/
    if rep.fitsio.FITS is object:                          # the import-time stand-in
        monkeypatch.setattr(rep.fitsio, 'FITS', _AstropyAsFitsio)
    monkeypatch.setattr(rep.plt, 'savefig', lambda *a, **kw: None)

    out = rep.import_night_img(_Args(datdir, target=_report_target('Sp2135-3153 TOI-4336--c')))
    targs, gaia, jdstart = out[-2], out[-1], out[6]

    assert targs == ['Sp2135-3153', 'TOI-4336--c', 'TOI-4336--c']
    assert sorted(gaia) == sorted(stars['Sp2135-3153'] + stars['TOI-4336--c'])
    assert len(jdstart) == len(gaia)


def test_overview_labels_a_target_by_its_name(monkeypatch, tmp_path):
    import matplotlib.axes
    datdir = str(tmp_path)
    os.makedirs(os.path.join(datdir, 'reports', 'temp'))
    imgs, vals, jdstart, jdend, nimages, newap, ellip, targs, gaia, flags = _synthetic_inputs(1)
    labels = []
    text = matplotlib.axes.Axes.text

    def record(self, x, y, s, *a, **kw):
        labels.append(s)
        return text(self, x, y, s, *a, **kw)

    monkeypatch.setattr(matplotlib.axes.Axes, 'text', record)
    rep.make_pdf(imgs, os.path.join(datdir, 'reports', 'Io_20260924.pdf'), _Args(datdir), [], vals, jdstart,
                 jdend, nimages, newap, flags, ellip, rep.get_chronological_order(jdstart), ['TOI-4336--c'], gaia)
    assert gaia[0] + '\n (TOI-4336 c)' in labels


# ----------------------------------------- two nights of a telescope at once

def test_report_temp_files_are_per_night():
    """masterbias.png, asm.png, ... have fixed names: a shared reports/temp let two nights reported at once
    overwrite and delete each other's figures."""
    import inspect
    a, b = _Args('/d'), _Args('/d')
    b.date = '20260925'
    assert rep.temp_dir(a) == os.path.join('/d', 'reports', 'temp', '20260924')
    assert rep.temp_dir(a) != rep.temp_dir(b)
    code = [l for l in inspect.getsource(rep).splitlines() if not l.strip().startswith('#')]
    assert not [l for l in code if '/reports/temp' in l]


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
