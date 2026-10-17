"""
Two pipeline runs on one telescope at once (two nights): the shared state that had to go.

  - wcsfit's catcache. wcsfit reuses a cached copy of a reference table in
    preference to the table itself, and the cache was one directory per
    telescope. Each call now runs in a temporary directory of its own, so a
    rerun always solves against the reference it is given, and nothing is
    left in (or wiped from) the telescope directory.
  - reports/Pipeline_Report and reports/QC: one table per telescope, read,
    changed and rewritten by every stage. Now under a lock and replaced
    atomically, so concurrent writers lose no rows.
  - output/<gaia_id>_<filter>_output.fits, one per star gathering every
    night, and output/StackImages/already_run.txt, the telescope's run
    targets: updated under the same kind of lock. Two runs interleaving on
    already_run.txt left a doubled header and another target's row
    (side-by-side test on two Io nights, 2026-10-17).

Run:  pytest src/tests/test_concurrent_nights.py -q
"""
import os
import stat
import subprocess
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)


# ---------------------------------------------------------------- wcsfit


FAKE_WCSFIT = """#!/bin/sh
# stands in for CASU wcsfit: caches into ./catcache and leaves a match file, as the real one does
pwd >> "$WCSFIT_LOG"
echo "$@" >> "$WCSFIT_LOG"
mkdir -p catcache
echo "cch_00000001 $6 0 1 0 1" >> catcache/index
: > match_frame.fits
"""


@pytest.fixture
def fake_wcsfit(tmp_path, monkeypatch):
    bindir = tmp_path / 'bin'
    bindir.mkdir()
    exe = bindir / 'wcsfit'
    exe.write_text(FAKE_WCSFIT)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / 'wcsfit.log'
    monkeypatch.setenv('PATH', '{}{}{}'.format(bindir, os.pathsep, os.environ['PATH']))
    monkeypatch.setenv('WCSFIT_LOG', str(log))
    return log


def test_wcsfit_runs_in_a_temporary_directory_of_its_own(tmp_path, monkeypatch, fake_wcsfit):
    from photometry import casutools
    datdir = tmp_path / 'PipelineOutput' / 'v3' / 'Callisto'
    datdir.mkdir(parents=True)
    monkeypatch.chdir(datdir)                       # the pipeline runs from the telescope directory

    for _ in range(2):
        casutools.wcsfit('frame.fits', 'frame.cat', 'vizgaia3', verbose=False, wcsref='frame_localfits_ref.fits')

    assert sorted(os.listdir(datdir)) == []         # no catcache, no match files
    lines = fake_wcsfit.read_text().splitlines()
    cwds, argvs = lines[0::2], lines[1::2]
    assert len(cwds) == 2 and cwds[0] != cwds[1]    # a fresh directory per call ...
    assert not any(os.path.exists(c) for c in cwds)  # ... removed afterwards
    for argv in argvs:
        args = argv.split()
        assert args[:2] == [str(datdir / 'frame.fits'), str(datdir / 'frame.cat')]
        assert args[-1] == str(datdir / 'frame_localfits_ref.fits')


def test_wcsfit_directory_removed_when_wcsfit_fails(tmp_path, monkeypatch, fake_wcsfit):
    from photometry import casutools
    exe = tmp_path / 'bin' / 'wcsfit'
    exe.write_text(FAKE_WCSFIT + 'exit 3\n')
    monkeypatch.chdir(tmp_path)
    with pytest.raises(subprocess.CalledProcessError):
        casutools.wcsfit('frame.fits', 'frame.cat', 'vizgaia3', verbose=False, wcsref='ref.fits')
    cwd = fake_wcsfit.read_text().splitlines()[0]
    assert not os.path.exists(cwd)


# ---------------------------------------------------------------- report tables


WRITERS = 6
ROWS_PER_WRITER = 8


def _env():
    env = dict(os.environ)
    env['PYTHONPATH'] = SRC + os.pathsep + env.get('PYTHONPATH', '')
    return env


def _run_writers(script):
    procs = [subprocess.Popen([sys.executable, '-c', script.format(w=w, n=ROWS_PER_WRITER)], env=_env(),
                              stderr=subprocess.PIPE) for w in range(WRITERS)]
    for p in procs:
        _, err = p.communicate(timeout=300)
        assert p.returncode == 0, err.decode()


def test_pipeline_report_keeps_every_row_with_concurrent_writers(tmp_path):
    pytest.importorskip('tabulate')
    from astropy.io import ascii
    fname = str(tmp_path / 'Pipeline_Report')
    # each writer is one night of the telescope, reporting stages for its targets as ZLP_pipeline.sh does
    _run_writers("""
from reporting.report import main
for i in range({n}):
    target = 'SP{w}%03d' % i
    for t in (1, 6):
        main(['report.py', FNAME, '2025030{w}', target, 'Callisto', '1', str(t)])
""".replace('FNAME', repr(fname)))
    table = ascii.read(fname, format='fixed_width')
    rows = [r for r in table if str(r['Night']).startswith('2025')]
    assert len(rows) == WRITERS * ROWS_PER_WRITER
    assert all(str(r['T6']) == '[x]' and str(r['T7']) == '[o]' for r in rows)
    assert not [f for f in os.listdir(tmp_path) if '.tmp' in f]


def test_qc_table_keeps_every_row_with_concurrent_writers(tmp_path):
    pytest.importorskip('tabulate')
    from astropy.io import ascii
    fname = str(tmp_path / 'QC')
    _run_writers("""
from reporting.QC import main
for i in range({n}):
    target = 'SP{w}%03d' % i
    main([0, FNAME, '2025030{w}', target, 'Callisto', '1', 'QC10', 500])
    main([0, FNAME, '2025030{w}', target, 'Callisto', '1', 'QC11', 90])
""".replace('FNAME', repr(fname)))
    table = ascii.read(fname, format='fixed_width')
    rows = [r for r in table if str(r['Night']).startswith('2025')]
    assert len(rows) == WRITERS * ROWS_PER_WRITER
    assert all(str(r['QC10']) == '[x]' and str(r['QC11']) == '[x]' for r in rows)


def test_lock_excludes_other_processes(tmp_path):
    import time
    from utils.shared_files import acquire
    fname = str(tmp_path / '123_I+z_output.fits')
    holder = subprocess.Popen([sys.executable, '-c', """
import os, sys, time
from utils.shared_files import acquire
fd = acquire(sys.argv[1])
print('held', flush=True)
time.sleep(1.5)
os.close(fd)
""", fname], env=_env(), stdout=subprocess.PIPE)
    assert holder.stdout.readline().strip() == b'held'
    t0 = time.monotonic()
    fd = acquire(fname)
    waited = time.monotonic() - t0
    os.close(fd)
    holder.wait(timeout=30)
    assert waited > 1.0


def test_already_run_list_waits_for_the_lock(tmp_path, capsys):
    pytest.importorskip('fitsio')
    import time
    from utils.already_run_targs import already_run_targs
    outdir = tmp_path / 'output'
    (outdir / 'StackImages').mkdir(parents=True)
    runlist = outdir / 'StackImages' / 'already_run.txt'
    runlist.write_text('Gaia_ID SPECULOOS_ID\n0 SP0000+0000\n1 SP0001+0001\n')
    holder = subprocess.Popen([sys.executable, '-c', """
import os, sys, time
from utils.shared_files import acquire
fd = acquire(sys.argv[1])
print('held', flush=True)
time.sleep(1.5)
os.close(fd)
""", str(runlist)], env=_env(), stdout=subprocess.PIPE)
    assert holder.stdout.readline().strip() == b'held'
    t0 = time.monotonic()
    already_run_targs('unused', str(outdir), 'Sp0002+0002', '2')
    waited = time.monotonic() - t0
    holder.wait(timeout=30)
    assert waited > 1.0
    assert capsys.readouterr().out.strip() == 'SP0002+0002'     # a newly added target comes back upper case
    assert 'SP0002+0002' in runlist.read_text()
