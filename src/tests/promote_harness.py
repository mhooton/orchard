#!/usr/bin/env python3
"""
promote_harness.py
==================

Drives utils/promote.py through its command line on toy data roots (tests/promote_toys.py) and prints what it
did, for the four cases T12 gets wrong:

  1. repeated supersession of one night (T12 nests on the second and fails on the third)
  2. a single-target run on a multi-target night (T12 swaps the whole night's reduction)
  3. stale files in a reused v3 directory (T12 sweeps whatever is there)
  4. an incumbent with aperture 5 against a candidate without it (T12's *4_diff gate lets it through)

    python tests/promote_harness.py [--keep DIR]

Everything happens under a temporary directory (or DIR); nights are in 2001, so no real pipeline run can hold
them.  Run it from src/ (inside orchard-server: PYTHONPATH=<src> python tests/promote_harness.py).
"""
import argparse
import datetime as dt
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, SRC)

from promote_toys import Toy, gaia_id, tree_lines  # noqa: E402

TEL = 'Europa'


def now():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


class Harness(object):
    def __init__(self, root):
        self.root = root
        self.failures = 0

    def case(self, name):
        base = os.path.join(self.root, name, 'data')
        os.makedirs(base)
        self.base, self.queue = base, os.path.join(self.root, name, 'queue')
        self.toy = Toy(base)
        self.env = dict(os.environ, ORCHARD_QUEUE_ROOT=self.queue,
                        GAIADATABASEPATH=os.path.join(self.root, 'no-gaia.db'),
                        PYTHONPATH=SRC + os.pathsep + os.environ.get('PYTHONPATH', ''), PYTHONDONTWRITEBYTECODE='1')
        self.env.pop('TARGET_LIST', None)
        self.env.pop('ORCHARD_QUEUE_CONFIG', None)
        self.env.pop('JOBQUEUE_JOB_ID', None)

    def promote(self, *args):
        cmd = [sys.executable, os.path.join(SRC, 'utils', 'promote.py')] + list(args)
        shown = ' '.join(a if ' ' not in a else '"{}"'.format(a) for a in args)
        print('\n$ promote.py {}'.format(shown))
        # every path explicit, so nothing can fall back to the production queue root
        paths = ['--base', self.base, '--db', os.path.join(self.queue, 'promote.sqlite'),
                 '--queue-db', os.path.join(self.queue, 'queue.sqlite')]
        p = subprocess.run(cmd + paths, env=self.env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True)
        out = p.stdout.replace(self.base + '/', '').replace(self.base, '<base>')
        for line in out.rstrip('\n').splitlines():
            print('  ' + line)
        print('  [exit {}]'.format(p.returncode))
        return p.returncode

    def tree(self, title, top='PipelineOutput'):
        print('\n{}:'.format(title))
        for line in tree_lines(self.base, top):
            print('    ' + line)

    def check(self, ok, text):
        print('CHECK {}: {}'.format('ok' if ok else 'FAILED', text))
        if not ok:
            self.failures += 1

    def path(self, *p):
        return os.path.join(self.base, 'PipelineOutput', *p)


def case_repeated_supersession(h):
    print('=' * 100)
    print('1. Repeated supersession of one night: three reruns of Europa 20010101 (two targets)')
    print('   T12 moves v2/<night>/<target> to v2_old/<night>/<target>: the second time it lands inside the first')
    print('   copy, the third time mv fails ("File exists") and errexit kills the run.')
    h.case('1_repeated')
    a, b = 'Sp0451+0305', 'Sp2209-3707'
    for k in range(4):
        h.toy.night(TEL, '20010101', [a, b], now() - dt.timedelta(minutes=30))
        print('\n--- run {} of the night is in v3'.format(k + 1))
        h.promote('run', '--apply', '--telescope', TEL, '--night', '20010101')
        time.sleep(1.1)
    h.tree('PipelineOutput/superseded after four runs', os.path.join('PipelineOutput', 'superseded'))
    sup = h.path('superseded', TEL, '20010101', a)
    stamps = sorted(os.listdir(sup))
    h.check(len(stamps) == 3, 'three superseded copies of {}, one per later run: {}'.format(a, ', '.join(stamps)))
    h.check(all(a not in os.listdir(os.path.join(sup, s)) for s in stamps), 'none nested inside another')
    h.check(all(sorted(os.listdir(os.path.join(sup, s, '1'))) == ['I+z_processed.dat', 'I+z_stacked.dat']
                for s in stamps), 'their frames are gone, the frame lists and products kept')
    h.check(len(os.listdir(h.path('superseded', TEL, '20010101', 'reduction'))) == 3,
            'the night reduction: installed by run 1, replaced (old one kept) by runs 2-4, which covered every target')
    h.promote('list', '--all')


def case_single_target(h):
    print('=' * 100)
    print('2. A single-target run on a multi-target night: Europa 20010102, targets A and B in v2; B is rerun alone')
    print("   T12 replaces the night's whole v2 reduction, so A's products end up beside B's rerun masters.")
    h.case('2_single_target')
    a, b = 'Sp0451+0305', 'Sp2209-3707'
    for t in (a, b):
        h.toy.raw(TEL, '20010102', t, 121)
        h.toy.plan(TEL, '20010102', t, gaia_id(t))
        h.toy.target_dir('v2', TEL, '20010102', t)
    h.toy.reduction('v2', TEL, '20010102', [a, b], tag='old')
    h.toy.log('v2', TEL, '20010102', now() - dt.timedelta(days=3), [a, b])
    a_ino = os.stat(h.path('v2', TEL, 'output', '20010102', a)).st_ino
    h.toy.log('v3', TEL, '20010102', now() - dt.timedelta(minutes=20), [b], specified=True)
    h.toy.target_dir('v3', TEL, '20010102', b)
    h.toy.reduction('v3', TEL, '20010102', [a, b], tag='new')
    h.toy.pdf('v3', TEL, '20010102')
    h.promote('run', '--telescope', TEL, '--night', '20010102')
    h.promote('run', '--apply', '--telescope', TEL, '--night', '20010102')
    h.promote('show', '1')
    from astropy.io import fits
    tag = fits.getheader(h.path('v2', TEL, 'output', '20010102', 'reduction', '1_MasterBias.fits'))['TAG']
    h.check(tag == 'old', "v2's reduction is still the full night's ({!r})".format(tag))
    h.check(os.stat(h.path('v2', TEL, 'output', '20010102', a)).st_ino == a_ino, 'A is untouched in v2')
    h.check(os.path.isdir(h.path('v2', TEL, 'output', '20010102', b)), "B's rerun is in v2")
    h.tree('PipelineOutput after the promotion')


def case_stale_files(h):
    print('=' * 100)
    print('3. Stale files in a reused v3 directory: Europa 20010103, an earlier run left apertures 3-8; the rerun')
    print('   rewrites 3-5 only, so the old aperture-8 light curve sits beside the new ones.')
    h.case('3_stale')
    a = 'Sp0451+0305'
    h.toy.raw(TEL, '20010103', a, 121)
    h.toy.plan(TEL, '20010103', a, gaia_id(a))
    h.toy.target_dir('v3', TEL, '20010103', a, aps=(3, 4, 5, 8))
    time.sleep(2.2)
    start = now()
    time.sleep(1.1)
    h.toy.log('v3', TEL, '20010103', start, [a])
    h.toy.target_dir('v3', TEL, '20010103', a, aps=(3, 4, 5))
    h.toy.reduction('v3', TEL, '20010103', [a])
    h.toy.pdf('v3', TEL, '20010103')
    h.promote('run', '--apply', '--telescope', TEL, '--night', '20010103')
    h.check(os.path.isdir(h.path('v3', TEL, 'output', '20010103', a)) and
            not os.path.exists(h.path('v2', TEL, 'output', '20010103', a)), 'HOLD: nothing moved')
    h.promote('run', '--apply', '--telescope', TEL, '--night', '20010103')
    h.promote('list')
    h.promote('accept', '1', '--reason', 'aperture 8 is not used')
    h.promote('accept', '1', '--reason', 'aperture 8 is not used', '--apply', '--by', 'harness')
    h.check(os.path.isdir(h.path('v2', TEL, 'output', '20010103', a)), 'accepted: in v2')
    h.promote('show', '1')


def case_aperture_5(h):
    print('=' * 100)
    print('4. An incumbent with aperture 5 against a candidate without it: Europa 20010104. v2 serves a 120-point')
    print('   aperture-5 curve; the rerun reached apertures 3-4 only. T12 checks *4_diff.fits, so it would sweep it')
    print('   and the portal (which reads *_5_diff.fits only) would lose the night.')
    h.case('4_aperture_5')
    a = 'Sp0451+0305'
    h.toy.raw(TEL, '20010104', a, 121)
    h.toy.plan(TEL, '20010104', a, gaia_id(a))
    h.toy.target_dir('v2', TEL, '20010104', a, aps=(3, 4, 5))
    h.toy.log('v3', TEL, '20010104', now() - dt.timedelta(minutes=20), [a])
    h.toy.target_dir('v3', TEL, '20010104', a, aps=(3, 4))
    h.toy.reduction('v3', TEL, '20010104', [a])
    h.promote('run', '--apply', '--telescope', TEL, '--night', '20010104')
    diffs = sorted(n for n in os.listdir(h.path('v2', TEL, 'output', '20010104', a)) if n.endswith('_diff.fits'))
    h.check(any(n.endswith('_5_diff.fits') for n in diffs), 'v2 still serves aperture 5: ' + ', '.join(diffs))
    h.promote('reject', '1', '--reason', 'keep the aperture-5 curve in v2', '--by', 'harness')
    h.promote('run', '--apply', '--telescope', TEL, '--night', '20010104')
    h.promote('list', '--all')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--keep', help='work here and keep it (default: a temporary directory, removed afterwards)')
    a = p.parse_args()
    root = a.keep or tempfile.mkdtemp(prefix='promote_harness_')
    os.makedirs(root, exist_ok=True)
    print('promote harness, {} UTC, python {}, code {}'.format(now(), sys.version.split()[0], SRC))
    h = Harness(root)
    try:
        for case in (case_repeated_supersession, case_single_target, case_stale_files, case_aperture_5):
            case(h)
    finally:
        if not a.keep:
            shutil.rmtree(root, ignore_errors=True)
    print('=' * 100)
    print('{} check(s) failed'.format(h.failures) if h.failures else 'all checks passed')
    return 1 if h.failures else 0


if __name__ == '__main__':
    sys.exit(main())
