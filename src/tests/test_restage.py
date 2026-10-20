"""
Tests for reruns from a late stage (docs/rerun-from-stage.md): utils/restage.py puts a promoted night's inputs back
into v3, a rerun from T10 writes new light curves beside them, and the promotion step (utils/promote.py) overlays v2
with only what the rerun wrote. Toy data roots (tests/promote_toys.py), nights in 2001, /proc pointed at an empty
directory: nothing here can touch a real night.
"""
import datetime as dt
import glob
import json
import os
import sys
import time
import types

import numpy as np
import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from astropy.io import fits  # noqa: E402

from promote_toys import FILTER, Toy, gaia_id, snapshot, write_diff  # noqa: E402
from utils import promote, restage  # noqa: E402
from utils.promote_journal import Journal  # noqa: E402

TEL, NIGHT = 'Europa', '20010101'
A, B = 'Sp0451+0305', 'Sp2209-3707'


def now():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'queue'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    monkeypatch.delenv('JOBQUEUE_JOB_ID', raising=False)
    monkeypatch.setenv('PYTHONPATH', SRC + os.pathsep + os.environ.get('PYTHONPATH', ''))
    proc = tmp_path / 'proc'
    proc.mkdir()
    monkeypatch.setattr(promote, 'PROC', str(proc))
    base = tmp_path / 'data'
    base.mkdir()
    e = types.SimpleNamespace(base=str(base), toy=Toy(str(base)), tmp=tmp_path,
                              db=str(tmp_path / 'queue' / 'promote.sqlite'),
                              qdb=str(tmp_path / 'queue' / 'queue.sqlite'))
    e.v2 = lambda *p: os.path.join(e.base, 'PipelineOutput', 'v2', TEL, *p)
    e.v3 = lambda *p: os.path.join(e.base, 'PipelineOutput', 'v3', TEL, *p)
    e.sup = lambda *p: os.path.join(e.base, 'PipelineOutput', 'superseded', TEL, *p)
    return e


def promoted(env, targets=(A,), aps=(3, 4, 5)):
    """A night as promotion leaves it: every target's directory, the reduction, the log and the PDF in v2; in v3 the
    reduction only."""
    old = now() - dt.timedelta(days=1)
    for t in targets:
        env.toy.raw(TEL, NIGHT, t, 121)
        env.toy.plan(TEL, NIGHT, t, gaia_id(t))
        env.toy.target_dir('v2', TEL, NIGHT, t, aps=aps)
    env.toy.reduction('v2', TEL, NIGHT, targets, tag='v2')
    env.toy.reduction('v3', TEL, NIGHT, targets, tag='v2')
    env.toy.log('v2', TEL, NIGHT, old, list(targets))
    env.toy.pdf('v2', TEL, NIGHT)


def stage(env, *args, first='10'):
    return restage.main(['stage', '--base', env.base, '--queue-db', env.qdb, '--telescope', TEL, '--night', NIGHT,
                         '--from', first] + list(args))


def promote_night(env, *extra):
    return promote.main(['run', '--apply', '--telescope', TEL, '--night', NIGHT, '--base', env.base, '--db', env.db,
                         '--queue-db', env.qdb] + list(extra))


def with_pwv(path, ap, frames=120):
    """Add what T10 writes when the PWV correction works: a PWV table and a finite PWV_LIGHTCURVE_<ap>."""
    with fits.open(path, mode='append') as h:
        h.append(fits.BinTableHDU.from_columns([fits.Column(name='JD', format='D', array=np.arange(5.0)),
                                                fits.Column(name='PWV', format='D', array=np.full(5, 2.1))],
                                               name='PWV'))
        h.append(fits.ImageHDU(np.ones((1, frames)), name='PWV_LIGHTCURVE_%d' % ap))


def start_run(env, target, stages, specified=True):
    time.sleep(2.2)                     # run starts are logged to the second, and S1 allows one second
    start = now()
    time.sleep(1.1)
    env.toy.log('v3', TEL, NIGHT, start, [target], specified=specified, stages=stages)
    return start


def rerun_t10(env, target, aps=(3, 4, 5), pwv=True, pdf=True, specified=True):
    """What ZLP_pipeline.sh --from_T10 --no_T12 writes into a restaged directory: a log, new curves (the same
    differential flux, now with PWV), the light-curve files, PWV files, and the night's PDF. Returns the start."""
    start = start_run(env, target, [10, 11], specified)
    write_t10(env, target, aps, pwv, pdf)
    return start


def rerun_t9(env, target, aps=(3, 4, 5)):
    """What ZLP_pipeline.sh --from_T9 --no_T12 writes: each frame's WCS written into it (the copy restage made), its
    per-frame files, the night's photometry and proper-motion catalogue, then T10's and T11's."""
    start = start_run(env, target, [9, 10, 11])
    d = env.v3('output', NIGHT, target)
    for f in sorted(glob.glob(os.path.join(d, '1', 'proc*.fits'))):
        with fits.open(f, mode='update') as h:                       # in place, as wcsfit does
            h[0].header['CD1_1'] = 9.7e-5
            h[0].header['WCSCOMPL'] = True
        for name in (f + '.phot', f + '.ell', f[:-len('.fits')] + '.cat', f[:-len('.fits')] + '_psf.png'):
            with open(name, 'w') as out:
                out.write('rerun')
    with open(os.path.join(d, '1', '{}_processed.dat_phot'.format(FILTER)), 'w') as out:
        out.write('rerun')
    gid = gaia_id(target)
    fits.PrimaryHDU(header=fits.Header([('RERUN', 9)])).writeto(
        os.path.join(d, '{}_{}_{}_output.fits'.format(gid, FILTER, NIGHT)), overwrite=True)
    fits.PrimaryHDU().writeto(os.path.join(d, '{}_stack_catalogue_pm_{}.fits'.format(target.upper(), FILTER)),
                              overwrite=True)
    write_t10(env, target, aps, True, True)
    return start


def write_t10(env, target, aps, pwv, pdf):
    d = env.v3('output', NIGHT, target)
    gid = gaia_id(target)
    os.makedirs(os.path.join(d, 'lightcurves'), exist_ok=True)
    for ap in aps:
        path = os.path.join(d, '{}_{}_{}_{}_diff.fits'.format(gid, FILTER, NIGHT, ap))
        write_diff(path, [gid], ap, seed=ap)
        if pwv:
            with_pwv(path, ap)
        with open(os.path.join(d, 'lightcurves', '{}_{}_{}_{}_MCMC'.format(gid, FILTER, NIGHT, ap)), 'w') as f:
            f.write('rerun {}\n'.format(ap))
    with open(os.path.join(d, 'lightcurves', '{}_{}_{}_5.png'.format(gid, FILTER, NIGHT)), 'w') as f:
        f.write('rerun png')
    with open(os.path.join(d, '{}_{}_{}_pwv_processed.fits'.format(gid, FILTER, NIGHT)), 'w') as f:
        f.write('pwv')
    if pdf:
        p = env.toy.pdf('v3', TEL, NIGHT)
        with open(p, 'a') as f:
            f.write('rerun with PWV\n')


def verdicts(env):
    j = Journal(env.db, readonly=True)
    try:
        return j.verdicts()
    finally:
        j.close()


# ---------------------------------------------------------------------------- restaging


def test_a_dry_run_changes_nothing(env, capsys):
    promoted(env)
    before = snapshot(env.base)
    assert stage(env) == 0
    assert snapshot(env.base) == before
    assert 'would be restaged from v2 for a rerun from T10' in capsys.readouterr().out


def test_restage_links_exactly_the_inputs_of_t10_and_t11_from_v2(env):
    promoted(env)
    v2_before = snapshot(env.v2())
    assert stage(env, '--apply') == 0
    d = env.v3('output', NIGHT, A)
    got = sorted(restage.walk_files(d))
    out = '{}_{}_{}_output.fits'.format(gaia_id(A), FILTER, NIGHT)
    assert got == sorted([out] + ['1/proc{}-{:04d}-{}.fits'.format(A, k, FILTER) for k in range(6)])
    for rel in got:                                   # hard links: v2's own files, no copies
        assert os.path.samefile(os.path.join(d, rel), env.v2('output', NIGHT, A, rel))
    assert not os.path.exists(os.path.join(d, 'lightcurves'))      # T10's own outputs are not inputs
    rec = restage.load_record(env.base, TEL, NIGHT, A)
    assert rec['state'] == 'staged' and rec['source'] == 'v2' and rec['from_stage'] == 10
    assert rec['stages'] == [10, 11] and len(rec['files']) == len(got)
    assert {e['rel']: e['how'] for e in rec['reduction']}['1_image_{}.list'.format(A)] == 'present'
    assert snapshot(env.v2()) == v2_before            # v2 is only read
    # again: recognised, checked, left alone
    assert stage(env, '--apply') == 0
    assert restage.load_record(env.base, TEL, NIGHT, A)['made_at'] == rec['made_at']


def test_restage_fills_a_missing_v3_reduction_from_v2(env):
    promoted(env)
    import shutil
    shutil.rmtree(env.v3('output', NIGHT, 'reduction'))
    assert stage(env, '--apply') == 0
    red = env.v3('output', NIGHT, 'reduction')
    assert os.path.samefile(os.path.join(red, '1_image_{}.list'.format(A)),
                            env.v2('output', NIGHT, 'reduction', '1_image_{}.list'.format(A)))
    assert os.path.exists(os.path.join(red, '1_MasterBias.fits'))


def test_restage_refuses_what_it_cannot_do(env, capsys):
    promoted(env)
    assert stage(env, '--apply', first='8') == 2                     # not restaged from T8 back
    assert 'cannot be restaged' in capsys.readouterr().out
    os.makedirs(env.v3('output', NIGHT, A, '1'))                      # a v3 directory without the inputs
    assert stage(env, '--apply', '--target', A) == 1
    assert 'not this tool' in capsys.readouterr().out
    assert stage(env, '--apply', '--target', B) == 1
    assert restage.load_record(env.base, TEL, NIGHT, B) is None


def test_a_night_still_in_v3_is_restaged_in_place_with_its_old_outputs_set_aside(env):
    start = now() - dt.timedelta(hours=2)
    env.toy.night(TEL, NIGHT, [A], start)
    d = env.v3('output', NIGHT, A)
    open(os.path.join(d, 'notes.txt'), 'w').close()
    old = sorted(r for r in restage.walk_files(d) if r.endswith('_diff.fits') or r.startswith('lightcurves'))
    assert stage(env, '--apply') == 0
    rec = restage.load_record(env.base, TEL, NIGHT, A)
    assert rec['source'] == 'v3' and rec['set_aside']['files'] == old
    left = restage.walk_files(d)
    assert not [r for r in left if r.endswith('_diff.fits') or r.startswith('lightcurves')]
    aside = rec['set_aside']['dir']
    assert aside.startswith(env.sup(NIGHT, A)) and sorted(
        r for r in restage.walk_files(aside) if r != 'superseded.json') == old
    # the inputs it keeps are the earlier stages' products; a file no stage writes is not one
    names = {e['rel'] for e in rec['files']}
    assert '1/proc{}-0000-{}.fits'.format(A, FILTER) in names and 'notes.txt' not in names


# ---------------------------------------------------------------------------- the overlay


def test_a_rerun_from_t10_is_overlaid_on_v2_and_nothing_else_in_v2_changes(env):
    promoted(env, aps=(3, 4, 5, 8))
    v2dir = env.v2('output', NIGHT, A)
    ino = os.stat(v2dir).st_ino
    before = snapshot(v2dir)
    assert stage(env, '--apply') == 0
    rerun_t10(env, A)
    d = env.v3('output', NIGHT, A)
    written = {r: restage.file_sha1(os.path.join(d, r)) for r in restage.walk_files(d)
               if r.endswith('_diff.fits') or r.startswith('lightcurves') or 'pwv' in r}
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'PROMOTE' and v['state'] == 'promoted', v['reasons']
    assert 'overlay' in v['notes']
    assert os.stat(v2dir).st_ino == ino                               # v2 keeps its directory
    after = snapshot(v2dir)
    for rel, sha in written.items():                                  # the rerun's files, as written
        assert restage.file_sha1(os.path.join(v2dir, rel)) == sha
    assert promote.pwv_state(os.path.join(v2dir, '{}_{}_{}_5_diff.fits'.format(gaia_id(A), FILTER, NIGHT))) == \
        'applied'
    files = {r for r in before if not os.path.isdir(os.path.join(v2dir, r))}
    changed = {r for r in files if r in after and before[r] != after[r]}
    replaced = {r for r in written if r in before}
    assert changed == replaced                                        # nothing else in v2 changed
    gone = set(before) - set(after)
    ap8 = '{}_{}_{}_8_diff.fits'.format(gaia_id(A), FILTER, NIGHT)
    assert gone == {ap8}                                              # the old ap 8 the rerun did not remake
    # what it replaced, and the ap 8 curve, are kept under <stamp>-overlay
    (stamp_dir,) = [x for x in os.listdir(env.sup(NIGHT, A)) if x.endswith(promote.OVERLAY)]
    kept = set(restage.walk_files(env.sup(NIGHT, A, stamp_dir))) - {'superseded.json'}
    assert kept == replaced | {ap8}
    # the restaged inputs left v3 with their record; the links to v2's files were dropped
    assert not os.path.exists(env.v3('output', NIGHT, A))
    assert not os.path.exists(restage.record_path(env.base, TEL, NIGHT, A))
    (left,) = [x for x in os.listdir(env.sup(NIGHT, A)) if x.endswith(promote.RESTAGED)]
    assert os.path.exists(env.sup(NIGHT, A, left, 'restage.json'))
    frame = os.path.join(v2dir, '1', 'proc{}-0000-{}.fits'.format(A, FILTER))
    assert os.stat(frame).st_nlink == 1
    # the run's PDF covers the night's only target: it takes the plain name; the old one is kept
    pdf = env.v2('reports', '{}_{}.pdf'.format(TEL, NIGHT))
    assert 'rerun with PWV' in open(pdf).read()
    (old_pdf,) = [os.path.join(dp, n) for dp, _, ns in os.walk(env.sup(NIGHT, 'reports')) for n in ns]
    assert 'rerun with PWV' not in open(old_pdf).read()
    assert any(n.startswith(NIGHT + '_1_v3_' + A + '.run') for n in os.listdir(env.v2('logs')))
    j = Journal(env.db, readonly=True)
    (p,) = j.promotions()
    assert [s['step'] for s in j.steps(p['id'])] == list(promote.OVERLAY_STEPS) and p['state'] == 'done'
    j.close()


def test_on_a_night_of_two_targets_the_other_and_the_pdf_stay(env):
    promoted(env, targets=(A, B))
    b_before = snapshot(env.v2('output', NIGHT, B))
    pdf = env.v2('reports', '{}_{}.pdf'.format(TEL, NIGHT))
    pdf_text = open(pdf).read()
    assert stage(env, '--apply', '--target', A) == 0
    assert not os.path.exists(env.v3('output', NIGHT, B))
    rerun_t10(env, A)
    assert promote_night(env) == 0
    assert snapshot(env.v2('output', NIGHT, B)) == b_before
    assert open(pdf).read() == pdf_text                               # the rerun's PDF lacks B: a .run copy only
    assert any('.run' in n for n in os.listdir(env.v2('reports')))


def test_a_linked_input_changed_by_the_rerun_holds_for_a_person(env):
    promoted(env)
    assert stage(env, '--apply') == 0
    rerun_t10(env, A)
    out = env.v3('output', NIGHT, A, '{}_{}_{}_output.fits'.format(gaia_id(A), FILTER, NIGHT))
    with open(out, 'ab') as f:                                        # an in-place write through the hard link
        f.write(b'\0' * 2880)
    before = snapshot(env.base)
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'O1' in v['reasons'] and "v2's changed too" in v['reasons']
    j = Journal(env.db, readonly=True)
    assert j.latest_policy_decision(v['id'])['outcome'] == 'person'
    j.close()
    moved = {k for k in set(before) ^ set(snapshot(env.base)) if 'PipelineOutput/v' in k}
    assert not moved                                                  # nothing in v2 or v3 moved


def test_v2_changed_after_the_restage_is_left_for_a_rerun(env):
    promoted(env)
    assert stage(env, '--apply') == 0
    out = env.v2('output', NIGHT, A, '{}_{}_{}_output.fits'.format(gaia_id(A), FILTER, NIGHT))
    data = open(out, 'rb').read()
    os.unlink(out)                                                    # v2's file replaced (a new inode)
    with open(out, 'wb') as f:
        f.write(data)
    rerun_t10(env, A)
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'O2 v2 changed after the restage' in v['reasons']
    j = Journal(env.db, readonly=True)
    d = j.latest_policy_decision(v['id'])
    j.close()
    assert d['outcome'] == 'reprocess' and d['rerun']
    assert os.path.isdir(env.v3('output', NIGHT, A))


def test_a_restaged_night_whose_rerun_has_not_run_is_not_judged(env, capsys):
    promoted(env)
    assert stage(env, '--apply', first='11') == 0                    # T11's inputs include T10's curves
    assert promote_night(env) == 0
    assert 'which has not run yet' in capsys.readouterr().out
    assert not os.path.exists(env.db) or not verdicts(env)


def test_a_failed_rerun_is_reported_and_can_be_cleared(env, capsys):
    promoted(env)
    assert stage(env, '--apply') == 0
    assert promote_night(env) == 0
    assert 'v3 holds no light curves of it' in capsys.readouterr().out
    v2_before = snapshot(env.v2())
    assert restage.main(['clear', '--base', env.base, '--queue-db', env.qdb, '--telescope', TEL, '--night', NIGHT,
                         '--reason', 'rerun failed', '--apply']) == 0
    assert not os.path.exists(env.v3('output', NIGHT, A))
    assert restage.load_record(env.base, TEL, NIGHT, A) is None
    assert snapshot(env.v2()) == v2_before
    (d,) = os.listdir(env.sup(NIGHT, A))
    assert restage.walk_files(env.sup(NIGHT, A, d)) == ['restage.json', 'superseded.json']


def test_a_night_restaged_in_place_is_published_whole_without_s1_holds(env):
    start = now() - dt.timedelta(hours=2)
    env.toy.night(TEL, NIGHT, [A], start, pdf=False)
    assert stage(env, '--apply') == 0
    rerun_t10(env, A)
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'PROMOTE' and v['state'] == 'promoted', v['reasons']
    assert os.path.isdir(env.v2('output', NIGHT, A)) and not os.path.exists(env.v3('output', NIGHT, A))
    assert restage.load_record(env.base, TEL, NIGHT, A) is None       # the record went with the promotion
    assert any(os.path.exists(os.path.join(env.sup(NIGHT, A), d, 'restage.json'))
               for d in os.listdir(env.sup(NIGHT, A)))


# ---------------------------------------------------------------------------- reruns from T9


def frames(d):
    return sorted(r for r in restage.walk_files(d) if r.startswith('1/proc') and r.endswith('.fits'))


def test_a_restage_from_t9_copies_the_frames_and_points_the_frame_lists_at_them(env):
    promoted(env)
    v2d = env.v2('output', NIGHT, A)
    pm = os.path.join(v2d, '{}_stack_catalogue_pm_{}.fits'.format(A.upper(), FILTER))
    open(pm, 'w').close()                                               # an old T9 product: never an input
    v2_before = snapshot(env.v2())
    assert stage(env, '--apply', first='9') == 0
    d = env.v3('output', NIGHT, A)
    got = set(restage.walk_files(d))
    assert set(frames(d)) == set(frames(v2d)) and frames(d)
    for rel in frames(d):                                               # copies: T9 writes WCS into them
        assert not os.path.samefile(os.path.join(d, rel), os.path.join(v2d, rel))
    cat = '{}_stack_catalogue_{}.fits'.format(A.upper(), FILTER)
    assert os.path.samefile(os.path.join(d, cat), os.path.join(v2d, cat))   # read only: a link
    assert not got & {os.path.basename(pm), '{}_{}_{}_output.fits'.format(gaia_id(A), FILTER, NIGHT)}
    assert not [r for r in got if r.endswith(('.phot', '.ell', '.cat', '.png', '_diff.fits'))]
    lines = open(os.path.join(d, '1', '{}_processed.dat'.format(FILTER))).read().split()
    assert lines and all(x.startswith(os.path.join(d, '1') + os.sep) and os.path.exists(x) for x in lines)
    # the toy's catalogue is target-named, as before 2026-10: T9 is told which one the night's products used
    assert open(os.path.join(d, '.active_stack_catalogue_' + FILTER)).read().strip() == cat
    rec = restage.load_record(env.base, TEL, NIGHT, A)
    how = {e['rel']: (e['how'], e['edited']) for e in rec['files']}
    assert how['.active_stack_catalogue_' + FILTER] == ('made', False)
    assert all(how[r] == ('copy', True) for r in frames(d))
    assert how['1/{}_processed.dat'.format(FILTER)] == ('rewritten', False)
    assert snapshot(env.v2()) == v2_before


def test_a_restage_from_t9_will_not_guess_between_stack_catalogues(env, capsys):
    promoted(env)
    other = env.v2('output', NIGHT, A, 'SP0451+0305_OLD_stack_catalogue_{}.fits'.format(FILTER))
    open(other, 'w').close()
    assert stage(env, '--apply', first='9') == 1
    assert 'more than one' in capsys.readouterr().out and not os.path.exists(env.v3('output', NIGHT, A))
    # a Gaia-named one the pipeline's own rule picks needs no record
    open(env.v2('output', NIGHT, A, '4123_Europa_andor_{}_stack_catalogue_{}.fits'.format(NIGHT, FILTER)), 'w').close()
    assert stage(env, '--apply', first='9') == 0
    assert not os.path.exists(env.v3('output', NIGHT, A, '.active_stack_catalogue_' + FILTER))


def test_a_rerun_from_t9_replaces_v2s_frames_and_their_files_and_nothing_else(env):
    promoted(env)
    v2d = env.v2('output', NIGHT, A)
    before = snapshot(v2d)
    old_frames = {r: os.stat(os.path.join(v2d, r)).st_ino for r in frames(v2d)}
    assert stage(env, '--apply', first='9') == 0
    copies = {r: os.stat(env.v3('output', NIGHT, A, r)).st_ino for r in frames(v2d)}
    rerun_t9(env, A)
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'PROMOTE' and v['state'] == 'promoted', v['reasons']
    after = snapshot(v2d)
    for rel, ino in copies.items():                                    # the rerun's copies, with the new WCS
        assert os.stat(os.path.join(v2d, rel)).st_ino == ino != old_frames[rel]
        assert fits.getheader(os.path.join(v2d, rel))['CD1_1'] == 9.7e-5
    rerun = {r for r in after if not os.path.isdir(os.path.join(v2d, r))} - set(
        r for r in before if before[r] == after.get(r))
    for rel in rerun:
        assert rel.startswith(('1/proc', 'lightcurves', '1/' + FILTER + '_processed.dat_phot')) or \
            rel.endswith(('_diff.fits', '_output.fits', '_pwv_processed.fits', '_pm_' + FILTER + '.fits')), rel
    same = {r for r in before if before[r] == after.get(r) and not os.path.isdir(os.path.join(v2d, r))}
    assert {'1/{}_processed.dat'.format(FILTER), '1/{}_stacked.dat'.format(FILTER),
            '{}_stack_catalogue_{}.fits'.format(A.upper(), FILTER),
            '{}_outstack_{}.fits'.format(A.upper(), FILTER)} <= same
    # the old frames and per-frame files, set aside under <stamp>-overlay/1/, were thinned like any superseded copy
    (od,) = [x for x in os.listdir(env.sup(NIGHT, A)) if x.endswith(promote.OVERLAY)]
    left = restage.walk_files(env.sup(NIGHT, A, od, '1'))
    assert left == [] and sum(1 for _ in open(env.sup(NIGHT, A, od, 'thinned.tsv'))) == 6 * 5
    assert not os.path.exists(env.v3('output', NIGHT, A))


def test_a_read_only_input_changed_by_a_rerun_from_t9_holds(env):
    promoted(env)
    assert stage(env, '--apply', first='9') == 0
    rerun_t9(env, A)
    cat = env.v3('output', NIGHT, A, '{}_stack_catalogue_{}.fits'.format(A.upper(), FILTER))
    with open(cat, 'ab') as f:                                          # written through the hard link
        f.write(b'\0' * 2880)
    assert promote_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'O1 1 of the' in v['reasons'], v['reasons']


# ---------------------------------------------------------------------------- S1 after a run from a later stage


def test_s1_does_not_hold_the_inputs_of_a_run_that_started_at_t10(env):
    # a --no_T1 ... --no_T9 rerun in place, without a restage (jobs 2064-2077 were queued so)
    start = now() - dt.timedelta(hours=2)
    env.toy.night(TEL, NIGHT, [A], start, pdf=False)
    d = env.v3('output', NIGHT, A)
    for rel in restage.walk_files(d):                                 # the T10 run rewrites its own outputs
        if rel.endswith('_diff.fits') or rel.startswith('lightcurves'):
            os.unlink(os.path.join(d, rel))
    rerun_t10(env, A, specified=False)
    assert promote_night(env, '--no-policy') == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'PROMOTE', v['reasons']


def test_s1_still_holds_a_t10_output_the_run_from_t10_did_not_rewrite(env):
    start = now() - dt.timedelta(hours=2)
    env.toy.night(TEL, NIGHT, [A], start, pdf=False)
    d = env.v3('output', NIGHT, A)
    env.toy.target_dir('v3', TEL, NIGHT, A, aps=(8,))                 # an old aperture-8 curve
    for rel in restage.walk_files(d):
        if (rel.endswith('_diff.fits') and not rel.endswith('_8_diff.fits')) or rel.startswith('lightcurves'):
            os.unlink(os.path.join(d, rel))
    rerun_t10(env, A, specified=False)
    assert promote_night(env, '--no-policy') == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'S1 1 of its' in v['reasons'] and '_8_diff.fits' in v['reasons']


def test_a_logs_stages_line_says_whether_the_run_made_light_curves(tmp_path):
    toy = Toy(str(tmp_path))
    start = now()
    p = toy.log('v3', TEL, NIGHT, start, [A], specified=True, stages=[11])
    run = promote.parse_log(p)
    assert run['stages'] == [11] and run['science'] is False
    p = toy.log('v3', TEL, '20010102', start, [A], stages=[10, 11])
    assert promote.parse_log(p)['science'] is True
    assert promote.science_argv(['./main/ZLP_pipeline.sh', '--from_T10', '--no_T12', '1', '/d', NIGHT])
    assert not promote.science_argv(['./main/ZLP_pipeline.sh', '--from_T11', '--no_T12', '1', '/d', NIGHT])
