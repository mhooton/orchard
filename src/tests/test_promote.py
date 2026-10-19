"""
Tests for the promotion step (utils/promote.py, utils/promote_journal.py) on toy data roots
(tests/promote_toys.py). Nights are in 2001, which no real pipeline run can hold, and /proc is pointed at an
empty directory, so runs in the container never see the queue's own jobs.
"""
import datetime as dt
import json
import os
import sqlite3
import sys
import time
import types

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from astropy.io import fits  # noqa: E402

from promote_toys import Toy, gaia_id, snapshot  # noqa: E402
from utils import promote  # noqa: E402
from utils.promote_journal import Journal  # noqa: E402

TEL, NIGHT = 'Europa', '20010101'
A, B = 'Sp0451+0305', 'Sp2209-3707'


def now():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'queue'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    monkeypatch.setenv('GAIADATABASEPATH', str(tmp_path / 'no-gaia.db'))
    monkeypatch.delenv('TARGET_LIST', raising=False)
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


def run(env, *args):
    # every path explicit: nothing may fall back to the production queue root
    return promote.main(list(args) + ['--base', env.base, '--db', env.db, '--queue-db', env.qdb])


def apply_night(env, night=NIGHT, *extra):
    return run(env, 'run', '--apply', '--telescope', TEL, '--night', night, *extra)


def journal(env):
    return Journal(env.db, readonly=True)


def verdicts(env, **kw):
    j = journal(env)
    try:
        return j.verdicts(**kw)
    finally:
        j.close()


def master_tag(path):
    return fits.getheader(os.path.join(path, '1_MasterBias.fits'))['TAG']


# ---------------------------------------------------------------------------- the four harness cases


def test_dry_run_changes_nothing(env, tmp_path):
    env.toy.night(TEL, NIGHT, [A, B], now() - dt.timedelta(minutes=30))
    env.toy.target_dir('v2', TEL, NIGHT, A)
    before = snapshot(env.base)
    assert run(env, 'run', '--report', str(tmp_path / 'dry.md'), '--csv', str(tmp_path / 'dry.csv')) == 0
    assert snapshot(env.base) == before
    assert not os.path.exists(env.db)
    text = open(str(tmp_path / 'dry.md')).read()
    assert '2 PROMOTE' in text and 'nothing moved' in text


def test_new_target_night_is_published_with_reduction_and_records(env):
    t0 = time.time()
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    ino = os.stat(env.v3('output', NIGHT, A)).st_ino
    assert apply_night(env) == 0
    assert not os.path.exists(env.v3('output', NIGHT, A))
    assert os.stat(env.v2('output', NIGHT, A)).st_ino == ino          # a rename, not a copy
    diff5 = [n for n in os.listdir(env.v2('output', NIGHT, A)) if n.endswith('_5_diff.fits')]
    assert diff5 and all(os.path.getmtime(env.v2('output', NIGHT, A, n)) >= t0 - 1 for n in diff5)
    assert master_tag(env.v2('output', NIGHT, 'reduction')) == 'new'
    assert os.path.isdir(env.v3('output', NIGHT, 'reduction'))         # copied, v3 keeps its own
    logs = sorted(os.listdir(env.v2('logs')))
    assert logs[0] == NIGHT + '_1_v3.log' and logs[1].startswith(NIGHT + '_1_v3.run')
    assert os.path.samefile(env.v2('logs', logs[0]), env.v2('logs', logs[1]))
    pdfs = sorted(os.listdir(env.v2('reports')))
    assert pdfs[0] == '{}_{}.pdf'.format(TEL, NIGHT) and pdfs[1].startswith('{}_{}.run'.format(TEL, NIGHT))
    assert os.path.exists(env.v3('logs', NIGHT + '_1_v3.log'))         # copied, never moved
    (v,) = verdicts(env)
    assert v['state'] == 'promoted' and v['verdict'] == 'PROMOTE'


def test_repeated_supersession_of_one_night_never_nests(env):
    stamps = []
    for k in range(4):
        start = now() - dt.timedelta(minutes=30)
        env.toy.night(TEL, NIGHT, [A, B], start)
        assert apply_night(env) == 0
        assert os.listdir(env.v3('output', NIGHT)) == ['reduction']
        stamps.append(sorted(os.listdir(env.sup(NIGHT, A))) if os.path.isdir(env.sup(NIGHT, A)) else [])
        time.sleep(1.1)       # each rerun is a new run: a later start, and new stamps
    assert [len(s) for s in stamps] == [0, 1, 2, 3]
    for stamp in stamps[-1]:
        d = env.sup(NIGHT, A, stamp)
        names = os.listdir(d)
        assert A not in names                                           # nothing nested inside
        assert any(n.endswith('_5_diff.fits') for n in names) and 'superseded.json' in names
        assert sorted(os.listdir(os.path.join(d, '1'))) == ['I+z_processed.dat', 'I+z_stacked.dat']
        assert sum(1 for _ in open(os.path.join(d, 'thinned.tsv'))) == 6 * 5      # 6 frames, 4 sidecars each
        assert os.listdir(os.path.join(d, 'lightcurves'))
    # the night's reduction: installed by the first run, replaced by each later full run
    assert len(os.listdir(env.sup(NIGHT, 'reduction'))) == 3
    assert len([n for n in os.listdir(env.v2('logs')) if '.run' in n]) == 4
    j = journal(env)
    assert [p['state'] for p in j.promotions()] == ['done'] * 8
    j.close()


def test_single_target_run_on_a_multi_target_night_keeps_the_night_reduction(env):
    old = now() - dt.timedelta(days=3)
    for t in (A, B):
        env.toy.raw(TEL, NIGHT, t, 121)
        env.toy.plan(TEL, NIGHT, t, gaia_id(t))
        env.toy.target_dir('v2', TEL, NIGHT, t)
    env.toy.reduction('v2', TEL, NIGHT, [A, B], tag='old')
    env.toy.log('v2', TEL, NIGHT, old, [A, B])
    a_ino = os.stat(env.v2('output', NIGHT, A)).st_ino
    # a rerun of B alone: T1 still writes the whole night's reduction into v3
    start = now() - dt.timedelta(minutes=20)
    env.toy.log('v3', TEL, NIGHT, start, [B], specified=True)
    env.toy.target_dir('v3', TEL, NIGHT, B)
    env.toy.reduction('v3', TEL, NIGHT, [A, B], tag='new')
    env.toy.pdf('v3', TEL, NIGHT)
    assert apply_night(env) == 0
    assert os.stat(env.v2('output', NIGHT, A)).st_ino == a_ino          # the other target untouched
    assert os.path.isdir(env.v2('output', NIGHT, B)) and len(os.listdir(env.sup(NIGHT, B))) == 1
    assert master_tag(env.v2('output', NIGHT, 'reduction')) == 'old'   # not swapped under A's products
    assert not os.path.exists(env.sup(NIGHT, 'reduction'))
    assert any(n.startswith(NIGHT + '_1_v3_' + B + '.run') for n in os.listdir(env.v2('logs')))
    j = journal(env)
    (p,) = j.promotions()
    red = [s for s in j.steps(p['id']) if s['step'] == 'reduction'][0]
    j.close()
    assert red['state'] == 'skipped' and 'covered 1 of the 2 targets' in red['detail']['why']


def test_stale_files_in_a_reused_v3_directory_are_held_until_accepted(env):
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    d = env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4, 5, 8))     # an earlier run's directory
    time.sleep(2.2)
    start = now()
    time.sleep(1.1)
    env.toy.log('v3', TEL, NIGHT, start, [A])
    env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4, 5))            # the rerun rewrites all but ap 8
    env.toy.reduction('v3', TEL, NIGHT, [A])
    before = snapshot(env.base)
    assert apply_night(env, NIGHT, '--no-policy') == 0                          # holds wait for a person
    assert snapshot(env.base) == before                                         # a HOLD moves nothing
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and v['state'] == 'pending'
    assert 'S1 1 of its' in v['reasons'] and '_8_diff.fits' in v['reasons']
    # unchanged on the next run: the same verdict, not a new one
    assert apply_night(env, NIGHT, '--no-policy') == 0
    assert [x['id'] for x in verdicts(env)] == [v['id']]
    # a person accepts it
    assert run(env, 'accept', str(v['id']), '--reason', 'ap 8 is not used', '--telescope', TEL, '--night', NIGHT) == 0
    assert os.path.isdir(d)                                                     # accept without --apply: dry run
    assert run(env, 'accept', str(v['id']), '--reason', 'ap 8 is not used', '--apply', '--by', 'mh',
               '--telescope', TEL, '--night', NIGHT) == 0
    assert not os.path.exists(d) and os.path.isdir(env.v2('output', NIGHT, A))
    assert not [n for n in os.listdir(env.v2('output', NIGHT, A)) if n.endswith('_8_diff.fits')]   # set aside
    (v,) = verdicts(env)
    assert (v['state'], v['decided_by'], v['decision_reason']) == ('accepted', 'mh', 'ap 8 is not used')


def test_a_directory_the_run_did_not_write_is_held(env):
    env.toy.night(TEL, NIGHT, [A], now() + dt.timedelta(minutes=5))      # its run started after every file
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'wrote none of its' in v['reasons']
    assert os.path.isdir(env.v3('output', NIGHT, A))


def test_incumbent_with_aperture_5_against_a_candidate_without_is_rejected(env):
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    env.toy.target_dir('v2', TEL, NIGHT, A, aps=(3, 4, 5))
    env.toy.log('v3', TEL, NIGHT, now() - dt.timedelta(minutes=20), [A])
    env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4))
    before = snapshot(env.base)
    assert apply_night(env, NIGHT, '--no-policy') == 0
    assert snapshot(env.base) == before
    (v,) = verdicts(env)
    assert v['verdict'] == 'REJECT' and v['state'] == 'pending'
    assert 'A1 no aperture-5 light curve' in v['reasons']
    assert run(env, 'reject', str(v['id']), '--reason', 'v2 keeps its ap5 curve', '--by', 'mh') == 0
    assert apply_night(env) == 0                                  # rejected: left alone, not judged anew
    (v,) = verdicts(env)
    assert v['state'] == 'rejected' and snapshot(env.base) == before


# ---------------------------------------------------------------------------- crash and resume


def test_a_crash_after_publishing_is_finished_by_resume(env, monkeypatch):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    env.toy.target_dir('v2', TEL, NIGHT, A)

    def killed(self, plan, prom):
        raise KeyboardInterrupt('killed')

    original = promote.Promoter.step_touch
    monkeypatch.setattr(promote.Promoter, 'step_touch', killed)
    with pytest.raises(KeyboardInterrupt):
        apply_night(env)
    j = journal(env)
    (p,) = j.promotions()
    steps = {s['step']: s['state'] for s in j.steps(p['id'])}
    j.close()
    assert p['state'] == 'in_progress'
    assert steps == {'check': 'done', 'aside': 'skipped', 'supersede': 'done', 'publish': 'done', 'touch': 'intent'}
    assert os.path.isdir(env.v2('output', NIGHT, A)) and not os.path.exists(env.v2('output', NIGHT, 'reduction'))
    monkeypatch.setattr(promote.Promoter, 'step_touch', original)
    assert run(env, 'resume', '--telescope', TEL, '--night', NIGHT) == 0       # dry run
    assert run(env, 'resume', '--apply', '--telescope', TEL, '--night', NIGHT) == 0
    j = journal(env)
    (p,) = j.promotions()
    j.close()
    assert p['state'] == 'done'
    assert master_tag(env.v2('output', NIGHT, 'reduction')) == 'new'
    (stamp,) = os.listdir(env.sup(NIGHT, A))
    assert sorted(os.listdir(env.sup(NIGHT, A, stamp, '1'))) == ['I+z_processed.dat', 'I+z_stacked.dat']


def test_a_crash_between_the_two_renames_is_finished_by_the_next_run(env, monkeypatch):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    env.toy.target_dir('v2', TEL, NIGHT, A)

    def killed(self, plan, prom):
        raise KeyboardInterrupt('killed')

    original = promote.Promoter.step_publish
    monkeypatch.setattr(promote.Promoter, 'step_publish', killed)
    with pytest.raises(KeyboardInterrupt):
        apply_night(env)
    assert not os.path.exists(env.v2('output', NIGHT, A))                       # the moment v2 lacks the target
    monkeypatch.setattr(promote.Promoter, 'step_publish', original)
    assert apply_night(env) == 0                                                # run --apply resumes it first
    assert os.path.isdir(env.v2('output', NIGHT, A)) and not os.path.exists(env.v3('output', NIGHT, A))
    j = journal(env)
    assert [p['state'] for p in j.promotions()] == ['done']
    j.close()


def test_accept_refuses_a_candidate_that_changed_after_its_verdict(env):
    env.toy.night(TEL, NIGHT, [A], now() + dt.timedelta(minutes=5))
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    with open(env.v3('output', NIGHT, A, 'lightcurves', 'new.png'), 'w') as f:
        f.write('x')
    before = snapshot(env.base)
    assert run(env, 'accept', str(v['id']), '--reason', 'x', '--apply', '--telescope', TEL, '--night', NIGHT) == 2
    assert snapshot(env.base) == before


def test_a_failed_publish_puts_the_incumbent_back(env, monkeypatch):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    inc = env.toy.target_dir('v2', TEL, NIGHT, A)
    ino = os.stat(inc).st_ino

    def broken(self, plan, prom):
        raise OSError('no permission')

    monkeypatch.setattr(promote.Promoter, 'step_publish', broken)
    assert apply_night(env) == 1
    assert os.stat(inc).st_ino == ino and os.path.isdir(env.v3('output', NIGHT, A))
    j = journal(env)
    (p,) = j.promotions()
    j.close()
    assert p['state'] == 'failed' and 'v2 restored' in p['error']
    assert apply_night(env) == 0                         # the night waits for a person while that is unresolved
    assert os.path.isdir(env.v3('output', NIGHT, A))
    assert run(env, 'abandon', str(p['id']), '--reason', 'checked by hand', '--by', 'mh') == 0
    (v,) = verdicts(env)
    assert v['state'] == 'outdated'


# ---------------------------------------------------------------------------- queue, runs, busy nights


def queue_db(env, jobs):
    from jobqueue.store import Store
    path = env.qdb
    s = Store(path, check_local=False)
    ids = []
    for j in jobs:
        jid, _ = s.add_job(j.get('cls', 'P2'), j.get('kind', 'pipeline'), j.get('argv', ['./main/ZLP_pipeline.sh']),
                           telescope=TEL,
                           night=j.get('night', NIGHT), targets=j.get('targets', ''),
                           locks=['night:{}:{}'.format(TEL, j.get('night', NIGHT))])
        s.conn.execute('UPDATE jobs SET state = ?, started_at = ?, finished_at = ?, exit_code = ? WHERE id = ?',
                       (j['state'], j['start'].strftime('%Y-%m-%dT%H:%M:%SZ'),
                        j['end'].strftime('%Y-%m-%dT%H:%M:%SZ') if j.get('end') else None, j.get('exit'), jid))
        ids.append(jid)
    s.close()
    return path, ids


def test_queue_mode_judges_the_jobs_targets_and_a_busy_night_waits(env):
    start = now() - dt.timedelta(minutes=30)
    env.toy.night(TEL, NIGHT, [A, B], start)
    path, (done_job, running_job) = queue_db(env, [
        {'state': 'failed', 'start': start, 'end': now(), 'targets': B, 'exit': 1},
        {'state': 'running', 'start': now(), 'night': '20010102'}])
    assert run(env, 'run', '--apply', '--job', str(done_job), '--telescope', TEL, '--night', NIGHT) == 0
    assert os.path.isdir(env.v2('output', NIGHT, B)) and os.path.isdir(env.v3('output', NIGHT, A))
    (v,) = verdicts(env)
    assert v['target'] == B and v['run']['job'] == done_job
    # the running job holds 20010102: that night is left alone, and a job for it is told to retry
    env.toy.night(TEL, '20010102', [A], start)
    assert apply_night(env, '20010102') == 0
    assert os.path.isdir(env.v3('output', '20010102', A)) and len(verdicts(env)) == 1
    assert run(env, 'run', '--apply', '--job', str(running_job), '--telescope', TEL, '--night', '20010102') == \
        promote.EX_TEMPFAIL


def test_runs_come_from_logs_and_queue_jobs(env):
    start = now() - dt.timedelta(hours=2)
    log = env.toy.log('v3', TEL, NIGHT, start, [A, B])
    env.toy.log('v2', TEL, NIGHT, start - dt.timedelta(days=5), [B], specified=True)
    path, (jid,) = queue_db(env, [{'state': 'done', 'start': start - dt.timedelta(seconds=5), 'end': now(),
                                   'exit': 0}])
    runs = promote.discover_runs(promote.Tree(env.base), promote.open_queue(path), TEL, NIGHT)
    assert len(runs) == 2
    job = [r for r in runs if r['source'] == 'job'][0]
    assert job['log'] == log and job['targets'] == [A, B] and job['complete']
    assert promote.latest_run(runs, B)['job'] == jid
    assert promote.latest_run(runs, 'Sp0000+0000') is None
    single = [r for r in runs if r['source'] == 'log'][0]
    assert single['specified'] and single['targets'] == [B]


def test_a_crash_before_the_promotion_started_does_not_block_the_target(env):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    os.makedirs(os.path.dirname(env.db))
    j = Journal(env.db)
    j.add_verdict(None, TEL, NIGHT, A, env.v3('output', NIGHT, A), '', 'PROMOTE', '', '', {}, {}, 'x', None,
                  'promoting')
    j.close()
    assert apply_night(env) == 0
    states = sorted(v['state'] for v in verdicts(env))
    assert states == ['outdated', 'promoted']


# ---------------------------------------------------------------------------- pieces


def test_raw_frame_counts_are_cached_and_invalidated_when_the_night_changes(env):
    env.toy.raw(TEL, NIGHT, A, 10)
    cache = str(env.tmp / 'raw.sqlite')
    r = promote.CachedRawFrames(env.base, cache)
    assert r.count(TEL, NIGHT, A, 'I+z') == 10 and r.stats['counted'] == 1
    r = promote.CachedRawFrames(env.base, cache)
    assert r.count(TEL, NIGHT, A, 'I+z') == 10 and r.stats['cached'] == 1
    time.sleep(0.05)
    env.toy.raw(TEL, NIGHT, A, 12)                       # the backfill adds frames
    r = promote.CachedRawFrames(env.base, cache)
    assert r.count(TEL, NIGHT, A, 'I+z') == 12 and r.stats['counted'] == 1
    # a seed in promote_check's JSON format counts as made when that file was written
    seed = str(env.tmp / 'raw_cache.json')
    with open(seed, 'w') as f:
        json.dump({'{}/{}'.format(TEL, NIGHT): {'fits': {'{}\ti+z'.format(A.upper()): 3}}}, f)
    old = os.stat(os.path.join(env.base, 'Observations', TEL, 'images', NIGHT)).st_mtime - 100
    os.utime(seed, (old, old))
    r = promote.CachedRawFrames(env.base, None, seed=seed)
    assert r.count(TEL, NIGHT, A, 'I+z') == 12 and r.stats['seed_invalidated'] == 1
    os.utime(seed, None)
    r = promote.CachedRawFrames(env.base, None, seed=seed)
    assert r.count(TEL, NIGHT, A, 'I+z') == 3 and r.stats['seed'] == 1


def test_identity_from_the_target_table_holds_where_the_plan_disagrees(env):
    true, wrong = gaia_id('true star'), gaia_id('background star')
    with open(os.path.join(env.base, 'ml_40pc_v2.csv'), 'w') as f:
        f.write('Sp_ID,Gaia_DR2_ID,Gaia_DR3_ID,RA,DEC,T_eff\n{},,{},15.0,10.0,2600\n{},,{},16.0,11.0,2700\n'.format(
            A, true, B, gaia_id(B)))
    start = now() - dt.timedelta(minutes=30)
    env.toy.night(TEL, NIGHT, [A, B], start, gids={A: true})
    env.toy.plan(TEL, NIGHT, A, wrong)                   # the plan names the wrong star for A
    assert run(env, 'run', '--telescope', TEL, '--night', NIGHT, '--csv', str(env.tmp / 'plan.csv')) == 0
    assert run(env, 'run', '--identity', 'table', '--telescope', TEL, '--night', NIGHT, '--csv',
               str(env.tmp / 'table.csv')) == 0
    import csv
    plan = {r['target']: r for r in csv.DictReader(open(str(env.tmp / 'plan.csv')))}
    table = {r['target']: r for r in csv.DictReader(open(str(env.tmp / 'table.csv')))}
    assert plan[A]['verdict'] == 'HOLD' and 'A3 light curve is' in plan[A]['reasons']
    assert table[A]['verdict'] == 'HOLD' and 'A3 the plan file names {}'.format(wrong) in table[A]['reasons']
    assert table[A]['expected_id'] == true
    assert plan[B]['verdict'] == table[B]['verdict'] == 'PROMOTE'


def test_parse_log_reads_start_end_and_targets(env):
    start = dt.datetime(2026, 10, 17, 18, 33, 6)
    path = env.toy.log('v3', TEL, NIGHT, start, [A, B], end=start + dt.timedelta(minutes=50))
    r = promote.parse_log(path)
    assert (r['start'], r['end'], r['targets'], r['specified'], r['complete']) == \
        (start, start + dt.timedelta(minutes=50), [A, B], False, True)
    with open(path, 'w') as f:
        f.write('Sat Oct 17 18:33:06 BST 2026\n')                     # only UTC is trusted
    assert promote.parse_log(path) is None


def test_thin_keeps_products_and_frame_lists(env):
    d = env.toy.target_dir('v2', TEL, NIGHT, A)
    os.makedirs(os.path.join(d, '1', 'sub'))
    out = promote.thin(d)
    assert out['deleted'] == 30 and out['kept_lists'] == 2
    assert sorted(os.listdir(os.path.join(d, '1'))) == ['I+z_processed.dat', 'I+z_stacked.dat', 'sub']
    assert any(n.endswith('_5_diff.fits') for n in os.listdir(d))
    assert len(open(out['manifest']).read().splitlines()) == 30


def test_copy_no_overwrite_never_replaces(env):
    a, b = str(env.tmp / 'a.log'), str(env.tmp / 'out' / 'x.log')
    with open(a, 'w') as f:
        f.write('one')
    assert promote.copy_no_overwrite(a, b) == (b, True)
    assert promote.copy_no_overwrite(a, b) == (b, False)
    with open(a, 'w') as f:
        f.write('two')
    assert promote.copy_no_overwrite(a, b) == (str(env.tmp / 'out' / 'x.2.log'), True)
    assert open(b).read() == 'one'


def test_run_apply_over_several_nights_runs_one_child_per_night(env):
    start = now() - dt.timedelta(minutes=30)
    env.toy.night(TEL, '20010103', [A], start)
    env.toy.night(TEL, '20010104', [B], start)
    assert run(env, 'run', '--apply') == 0
    assert os.path.isdir(env.v2('output', '20010103', A)) and os.path.isdir(env.v2('output', '20010104', B))
    j = journal(env)
    argvs = [json.loads(r['argv']) for r in j.conn.execute('SELECT argv FROM invocations ORDER BY id')]
    j.close()
    assert [a[a.index('--night') + 1] for a in argvs] == ['20010103', '20010104']


def test_a_directory_that_cannot_be_moved_aborts_before_anything_moves(env):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    inc = env.toy.target_dir('v2', TEL, NIGHT, A)
    os.chmod(inc, 0o555)                   # moving a directory to a new parent rewrites its '..'
    try:
        before = snapshot(env.base)
        assert run(env, 'run', '--telescope', TEL, '--night', NIGHT, '--csv', str(env.tmp / 'dry.csv')) == 0
        import csv
        (row,) = list(csv.DictReader(open(str(env.tmp / 'dry.csv'))))
        assert 'BUT cannot write in ' + inc in row['action']
        assert apply_night(env) == 0
        assert snapshot(env.base) == before
        j = journal(env)
        (p,) = j.promotions()
        j.close()
        assert p['state'] == 'aborted' and 'cannot write in' in p['error']
        assert verdicts(env)[0]['state'] == 'outdated'
    finally:
        os.chmod(inc, 0o755)


def test_resuming_waits_while_the_night_is_busy(env, monkeypatch):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))

    def killed(self, plan, prom):
        raise KeyboardInterrupt('killed')

    original = promote.Promoter.step_records
    monkeypatch.setattr(promote.Promoter, 'step_records', killed)
    with pytest.raises(KeyboardInterrupt):
        apply_night(env)
    monkeypatch.setattr(promote.Promoter, 'step_records', original)
    path, (jid,) = queue_db(env, [{'state': 'running', 'start': now()}])        # a rerun of the night starts
    assert run(env, 'resume', '--apply', '--telescope', TEL, '--night', NIGHT) == promote.EX_TEMPFAIL
    assert not os.listdir(env.v2('logs')) if os.path.isdir(env.v2('logs')) else True
    con = sqlite3.connect(path)
    con.execute("UPDATE jobs SET state = 'done', finished_at = started_at WHERE id = ?", (jid,))
    con.commit()
    con.close()
    assert run(env, 'resume', '--apply', '--telescope', TEL, '--night', NIGHT) == 0
    assert any('.run' in n for n in os.listdir(env.v2('logs')))


def test_night_busy_names_queue_jobs_and_outside_processes_apart(env, tmp_path):
    from jobqueue.tests.test_misc import fake_proc
    path, (jid,) = queue_db(env, [{'state': 'running', 'start': now()}])
    con = sqlite3.connect(path)
    con.execute('UPDATE jobs SET pid = 700 WHERE id = ?', (jid,))
    con.commit()
    con.close()
    runner = ['python', '-m', 'jobqueue.runner', '--log', 'x', '--', './main/ZLP_pipeline.sh', '1', '/data', NIGHT, '8',
              '2', TEL]
    manual = ['bash', './main/ZLP_pipeline.sh', '1', '/data', '20010105', '8', '2', TEL]
    (tmp_path / 'fake').mkdir()                     # the env fixture's empty proc/ is already in tmp_path
    proc = fake_proc(tmp_path / 'fake', [(700, 700, runner), (800, 800, manual)])
    conn = promote.open_queue(path)
    assert promote.night_busy(conn, TEL, NIGHT, proc=proc) == 'queue job {} (pipeline) is running on it'.format(jid)
    assert promote.night_busy(conn, TEL, NIGHT, own_job=jid, proc=proc) == ''
    assert promote.night_busy(conn, TEL, '20010105', proc=proc).startswith('a process outside the queue holds it: pid 800')


def test_holds_can_be_listed_by_reason_and_decided_in_bulk(env, capsys):
    for night in ('20010106', '20010107', '20010108'):
        env.toy.night(TEL, night, [A], now() + dt.timedelta(minutes=5))         # S1: the run wrote none of it
        assert apply_night(env, night) == 0
    ids = [v['id'] for v in verdicts(env)]
    capsys.readouterr()
    assert run(env, 'list', '--grep', 's1 the run') == 0
    assert '3 verdict(s)' in capsys.readouterr().out
    assert run(env, 'reject', str(ids[0]), '--reason', 'old run', '--by', 'mh') == 0
    assert run(env, 'accept', '{}-{}'.format(ids[1], ids[2]), '--reason', 'checked', '--apply', '--by', 'mh') == 0
    states = {v['night']: v['state'] for v in verdicts(env)}
    assert states == {'20010106': 'rejected', '20010107': 'accepted', '20010108': 'accepted'}
    assert os.path.isdir(env.v2('output', '20010107', A)) and os.path.isdir(env.v2('output', '20010108', A))
    assert os.path.isdir(env.v3('output', '20010106', A))


def test_a_binary_with_two_table_rows_agrees_with_a_plan_naming_either(env):
    """Resolved binaries share an Sp name and have a row each: the plan naming the second is no disagreement."""
    bright, faint = gaia_id('bright component'), gaia_id('faint component')
    with open(os.path.join(env.base, 'ml_40pc_v2.csv'), 'w') as f:
        f.write('Sp_ID,Gaia_DR2_ID,Gaia_DR3_ID,RA,DEC,T_eff\n{0},{1},{1},15.0,10.0,3000\n{0},{2},{2},15.0,10.0,2600\n'
                .format(A, bright, faint))
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30), gids={A: faint})
    d = env.v3('output', NIGHT, A)
    env.toy.target_dir('v3', TEL, NIGHT, A, gid=faint, others=(bright,))       # both components have curves
    assert run(env, 'run', '--identity', 'table', '--telescope', TEL, '--night', NIGHT, '--csv',
               str(env.tmp / 'table.csv')) == 0
    import csv
    (row,) = list(csv.DictReader(open(str(env.tmp / 'table.csv'))))
    assert row['verdict'] == 'PROMOTE', row['reasons']
    assert row['expected_id'] == faint and row['expected_src'] == 'target table, plan agrees'
    assert os.path.isdir(d)


# ---------------------------------------------------------------------------- the approval policy


def policy_rows(env, vid):
    j = journal(env)
    try:
        return j.policy_decisions(vid)
    finally:
        j.close()


def stale_rerun(env):
    """An earlier run's directory, then a rerun that rewrites all but its aperture-8 curve: S1 holds that file."""
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    d = env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4, 5, 8))
    time.sleep(2.2)
    start = now()
    time.sleep(1.1)
    env.toy.log('v3', TEL, NIGHT, start, [A])
    env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4, 5))
    env.toy.reduction('v3', TEL, NIGHT, [A])
    return d


def write_solves(log, target, solved, failed, sources, error):
    """T7's plate-solve summaries for a target, as astrometry.py prints them, appended to a toy log."""
    total = solved + failed
    lines = ['TARGET: {}'.format(target)]
    lines += ['[{0}/{1}] proc{2}-ok{0:04d}.fits SUCCESS - Sources: 200→16 Gaia: 300→32 Matches: 12 '
              'CRPIX=(1.00,1.00) CRVAL=(15.000000,10.000000)'.format(k + 1, total, target) for k in range(solved)]
    lines += ['[{0}/{1}] proc{2}-bad{0:04d}.fits FAILED - Sources: {3}→16 Gaia: 300→32 Matches: 0 '
              'Error: {4}'.format(solved + k + 1, total, target, sources, error) for k in range(failed)]
    with open(log, 'a', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


def partial_night(env, sources, error):
    """121 raw frames, a light curve of 60: the other 61 failed to plate solve, for the given reason."""
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    log = env.toy.log('v3', TEL, NIGHT, now() - dt.timedelta(minutes=30), [A])
    write_solves(log, A, 60, 61, sources, error)
    env.toy.target_dir('v3', TEL, NIGHT, A, frames=60)
    env.toy.reduction('v3', TEL, NIGHT, [A])


def test_policy_promotes_a_candidate_whose_stale_files_are_clutter_and_sets_them_aside(env):
    stale_rerun(env)
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'S1 1 of its' in v['reasons']
    assert (v['state'], v['decided_by']) == ('accepted', 'policy:v1:clutter-only')
    names = os.listdir(env.v2('output', NIGHT, A))
    assert any(n.endswith('_5_diff.fits') for n in names) and not any(n.endswith('_8_diff.fits') for n in names)
    (aside,) = [s for s in os.listdir(env.sup(NIGHT, A)) if s.endswith(promote.STALE)]
    kept = os.listdir(env.sup(NIGHT, A, aside))
    assert 'superseded.json' in kept and any(n.endswith('_8_diff.fits') for n in kept)
    assert [(d['outcome'], d['rule']) for d in policy_rows(env, v['id'])] == [('promote', 'clutter-only')]


def test_policy_rejects_a_candidate_without_aperture_5_and_v2_keeps_its_own(env):
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    env.toy.target_dir('v2', TEL, NIGHT, A, aps=(3, 4, 5))
    env.toy.log('v3', TEL, NIGHT, now() - dt.timedelta(minutes=20), [A])
    env.toy.target_dir('v3', TEL, NIGHT, A, aps=(3, 4))
    before = snapshot(env.base)
    assert apply_night(env) == 0
    assert snapshot(env.base) == before
    (v,) = verdicts(env)
    assert v['verdict'] == 'REJECT' and v['state'] == 'rejected'
    assert v['decided_by'] == 'policy:v1:no-aperture-5' and 'v2 keeps its own' in v['decision_reason']
    assert apply_night(env) == 0 and len(policy_rows(env, v['id'])) == 1     # left alone from then on


def test_policy_publishes_a_night_the_pipeline_lost_frames_of_and_marks_it_for_a_rerun(env, capsys):
    partial_night(env, 250, 'Plate scale sanity check failed: solved=25.674 arcsec/px, expected=0.348 arcsec/px, '
                            'deviation=7277.6%')
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'A4' in v['reasons']
    assert (v['state'], v['decided_by']) == ('accepted', 'policy:v1:partial-pipeline')    # nothing in v2: publish
    (d,) = policy_rows(env, v['id'])
    assert (d['outcome'], d['rule'], d['rerun']) == ('promote', 'partial-pipeline', 1)
    assert d['evidence']['losses']['processing'] == 61 and d['evidence']['losses']['solved'] == 60
    assert os.path.isdir(env.v2('output', NIGHT, A)) and not os.path.exists(env.v3('output', NIGHT, A))
    capsys.readouterr()
    assert run(env, 'list', '--needs', 'rerun') == 0
    assert '1 verdict(s)' in capsys.readouterr().out
    assert run(env, 'list', '--needs', 'person') == 0
    assert '0 verdict(s)' in capsys.readouterr().out


def test_policy_promotes_a_partial_night_the_sky_explains(env):
    partial_night(env, 3, 'Not enough stars detected for plate solve: 3 (more than 4 are needed to match an 4-star '
                          'asterism)')
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and 'A4' in v['reasons']
    assert (v['state'], v['decided_by']) == ('accepted', 'policy:v1:partial-night-data')
    assert os.path.isdir(env.v2('output', NIGHT, A))


def test_a_t11_only_run_after_the_science_run_is_not_the_run_s1_judges(env):
    science = now() - dt.timedelta(minutes=40)
    env.toy.night(TEL, NIGHT, [A], science, pdf=False)          # the science run: it died before T11
    time.sleep(1.1)
    t11 = now()
    time.sleep(1.1)
    env.toy.pdf('v3', TEL, NIGHT)                               # an --only_T11 run writes the night's PDF
    path, (sci_job, t11_job) = queue_db(env, [
        {'state': 'failed', 'start': science - dt.timedelta(seconds=5), 'end': t11 - dt.timedelta(seconds=1),
         'exit': 1},
        {'state': 'done', 'start': t11, 'end': now() + dt.timedelta(seconds=5), 'exit': 0,
         'argv': ['./main/ZLP_pipeline.sh', '--no_T12', '--only_T11', '1', '/data', NIGHT, '8', '2', TEL]}])
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['verdict'] == 'PROMOTE' and v['state'] == 'promoted', v['reasons']
    assert v['run']['job'] == sci_job and [r['job'] for r in v['run']['later']] == [t11_job]
    assert [n for n in os.listdir(env.v2('reports')) if '-j{}'.format(t11_job) in n]


def test_dry_run_says_what_the_policy_would_do(env, capsys):
    stale_rerun(env)
    before = snapshot(env.base)
    assert run(env, 'run', '--telescope', TEL, '--night', NIGHT) == 0
    assert 'policy would promote (clutter-only)' in capsys.readouterr().out
    assert snapshot(env.base) == before and not os.path.exists(env.db)


def test_holds_recorded_before_the_policy_are_decided_by_the_next_run(env):
    stale_rerun(env)
    assert apply_night(env, NIGHT, '--no-policy') == 0
    (v,) = verdicts(env)
    assert v['state'] == 'pending'
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert (v['state'], v['decided_by']) == ('accepted', 'policy:v1:clutter-only')


def test_science_argv_tells_runs_that_make_light_curves_from_those_that_do_not():
    assert promote.science_argv(['./main/ZLP_pipeline.sh', '--no_T12', '1', '/d', NIGHT, '8', '2', TEL])
    assert not promote.science_argv(['./main/ZLP_pipeline.sh', '--no_T12', '--only_T11', '1', '/d', NIGHT])
    assert not promote.science_argv(['./main/ZLP_pipeline.sh', '--only_T12', '1', '/d', NIGHT])
    assert promote.science_argv(['./main/ZLP_pipeline.sh', '--only_T10', '1', '/d', NIGHT])
    assert not promote.science_argv(['./main/ZLP_pipeline.sh', '--no_T10', '1', '/d', NIGHT])


def add_pwv(d, target, applied, ap=5):
    """Give the target's aperture-ap curve in d T10's PWV HDUs: with data, or empty as when the LHATPRO fetch fails."""
    import numpy as np
    (path,) = [os.path.join(d, n) for n in os.listdir(d)
               if n.startswith(gaia_id(target)) and n.endswith('_{}_diff.fits'.format(ap))]
    with fits.open(path, mode='update') as h:
        shape = h['LIGHTCURVE_{}'.format(ap)].data.shape
        table = np.recarray(3 if applied else 0, dtype=[('JD', np.float64), ('PWV', np.float64)])
        if applied:
            table['JD'], table['PWV'] = [1.0, 2.0, 3.0], [2.0, 2.1, 2.2]
        h.append(fits.BinTableHDU(table, name='PWV'))
        h.append(fits.ImageHDU(np.full(shape, 1.0 if applied else np.nan), name='PWV_LIGHTCURVE_{}'.format(ap)))
    return path


def test_pwv_state_reads_the_pwv_table_and_curve(env):
    d = env.toy.target_dir('v3', TEL, NIGHT, A)
    plain = [os.path.join(d, n) for n in os.listdir(d) if n.endswith('_3_diff.fits')][0]
    assert promote.pwv_state(plain) == 'missing'                   # no PWV HDUs at all
    assert promote.pwv_state(add_pwv(d, A, applied=True)) == 'applied'
    assert promote.pwv_state(add_pwv(d, A, applied=False, ap=4)) == 'missing'
    assert promote.pwv_state(os.path.join(d, 'absent_diff.fits')) is None


def test_p1_holds_a_replacement_that_lost_the_pwv_correction_v2_has(env):
    env.toy.raw(TEL, NIGHT, A, 121)
    env.toy.plan(TEL, NIGHT, A, gaia_id(A))
    add_pwv(env.toy.target_dir('v2', TEL, NIGHT, A), A, applied=True)
    env.toy.log('v3', TEL, NIGHT, now() - dt.timedelta(minutes=20), [A])
    add_pwv(env.toy.target_dir('v3', TEL, NIGHT, A), A, applied=False)
    env.toy.reduction('v3', TEL, NIGHT, [A])
    before = snapshot(env.base)
    assert apply_night(env) == 0
    assert snapshot(env.base) == before                             # v2 keeps its PWV-corrected curve
    (v,) = verdicts(env)
    assert v['verdict'] == 'HOLD' and v['state'] == 'pending'
    assert v['reasons'].startswith('P1 the candidate has no PWV correction')
    assert [(d['outcome'], d['rule']) for d in policy_rows(env, v['id'])] == [('reprocess', 'pwv-lost')]


def test_a_gap_fill_without_pwv_is_promoted_with_a_note(env):
    env.toy.night(TEL, NIGHT, [A], now() - dt.timedelta(minutes=30))
    add_pwv(env.v3('output', NIGHT, A), A, applied=False)
    assert apply_night(env) == 0
    (v,) = verdicts(env)
    assert v['state'] == 'promoted' and 'no PWV correction in the candidate' in v['notes']
