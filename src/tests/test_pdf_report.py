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
    assert vals == [['unknown'] * 5]


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
