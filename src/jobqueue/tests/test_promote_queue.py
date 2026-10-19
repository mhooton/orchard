"""The queue's side of the promotion step: the promote job, chaining it after pipeline jobs, requeueing, the lock a
promotion run outside the queue holds, `q add --promote`, the daily report's section and the nightly switch."""
import datetime as dt
import json
import os

import pytest

from jobqueue import cli, jobspec
from jobqueue.config import queue_paths
from jobqueue.resources import external_locks, parse_promote_argv
from jobqueue.scheduler import eligible
from jobqueue.shadow_report import build_report
from jobqueue.tests.fakes import PROGS, FakeEso
from jobqueue.tests.test_misc import fake_proc
from jobqueue.watcher import Watcher
from utils.promote_journal import Journal

NIGHT = '20250101'


def pipeline(cfg, store, promote=True, targets=(), cls='P2'):
    return jobspec.add_pipeline(store, cfg, jobspec.pipeline_job(cfg, cls, 'Europa', NIGHT, run_targets=targets,
                                                                 lock_targets=['Sp0001+0001'], no_t12=promote),
                                promote=promote)


def test_the_promote_job_follows_its_pipeline_job_whatever_its_outcome(cfg, store):
    pl, created, pr = pipeline(cfg, store, targets=['Sp0001+0001'])
    assert pipeline(cfg, store, targets=['Sp0001+0001']) == (pl, False, pr)     # queued twice: one of each
    p, j = store.get_job(pl), store.get_job(pr)
    assert '--no_T12' in p['argv'] and created
    assert j['kind'] == 'promote' and j['depends_on'] == pl and not j['dep_requires_success']
    assert j['argv'][1:] == ['utils/promote.py', 'run', '--apply', '--telescope', 'Europa', '--night', NIGHT,
                             '--job', str(pl)]
    assert j['locks'] == ['night:Europa:20250101'] and not j['disk_heavy'] and j['cores'] == 2
    assert j['targets'] == 'Sp0001+0001' and j['class'] == 'P2' and j['priority'] >= 1e9
    assert j['env']['ORCHARD_QUEUE_ROOT'] == cfg['queue_root'] and j['cwd'] == cfg['src_dir']
    assert j['env']['OPENBLAS_NUM_THREADS'] == '1'
    now = dt.datetime(2026, 10, 18, 12)
    assert eligible([j], now, set(), {pl: 'queued'}) == []
    store.mark_running(pl)
    store.finish(pl, 1, signature='ValueError: boom')
    assert store.get_job(pr)['state'] == 'queued'                       # a failed run is still judged
    assert eligible([store.get_job(pr)], now, set(), store.states_of([pl])) != []


def test_a_cancelled_pipeline_job_cancels_its_promote_job(cfg, store):
    pl, _, pr = pipeline(cfg, store)
    store.cancel(pl)
    assert store.get_job(pr)['state'] == 'cancelled'


def test_promotion_needs_the_pipeline_job_without_T12(cfg, store):
    with pytest.raises(ValueError):
        jobspec.add_pipeline(store, cfg, jobspec.pipeline_job(cfg, 'P2', 'Europa', NIGHT), promote=True)
    assert store.jobs() == []


def test_requeueing_a_failed_pipeline_job_requeues_its_promote_job(cfg, store):
    pl, _, pr = pipeline(cfg, store)
    store.mark_running(pl)
    store.finish(pl, 1, signature="KeyError: 'gaia_dr3_id'")
    store.mark_running(pr)
    store.finish(pr, 0)
    out = store.requeue(signature='gaia_dr3_id', dry_run=True)
    assert out == [(pl, 'would requeue'), (pr, 'would requeue (promotes job {})'.format(pl))]
    assert store.get_job(pr)['state'] == 'done'
    out = store.requeue(signature='gaia_dr3_id')
    assert [x[0] for x in out] == [pl, pr]
    assert store.get_job(pl)['state'] == store.get_job(pr)['state'] == 'queued'
    assert eligible([store.get_job(pr)], dt.datetime(2026, 10, 18), set(), store.states_of([pl])) == []


def test_a_promotion_run_outside_the_queue_holds_its_night(tmp_path):
    apply = ['python', '/opt/orchard/src/utils/promote.py', 'accept', '12', '--reason', 'ok', '--apply',
             '--telescope', 'Europa', '--night', NIGHT]
    dry = ['python', 'utils/promote.py', 'run', '--telescope', 'Io', '--night', NIGHT]
    wide = ['python', 'utils/promote.py', 'run', '--apply']                  # the parent: its children hold the nights
    check = ['python', 'utils/promote_check.py', 'scan', '--apply', '--telescope', 'Io', '--night', NIGHT]
    proc = fake_proc(tmp_path, [(20, 20, apply), (21, 21, dry), (22, 22, wide), (23, 23, check)])
    assert external_locks(own_sessions=[], proc=proc) == {'night:Europa:20250101': 'pid 20 ' + ' '.join(apply)}
    assert external_locks(own_sessions=[20], proc=proc) == {}               # the queue's own promote job
    assert parse_promote_argv(dry) is None and parse_promote_argv(wide) is None


def test_q_add_promote_queues_a_pipeline_run_without_T12_and_a_promote_job(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'q'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    assert cli.main(['add', '--class', 'P2', '--telescope', 'Europa', '--night', NIGHT, '--lock-targets',
                     'Sp0001+0001', '--promote', '--campaign', 'ESO backfill']) == 0
    assert 'added job 1: P2 pipeline Europa 20250101 with --no_T12, then promote job 2' in capsys.readouterr().out
    assert cli.main(['add', '--class', 'P2', '--kind', 'promote', '--after', '1,7']) == 0
    out = capsys.readouterr().out
    assert 'already queued as promote job 2 after pipeline job 1' in out and 'job 7: not a pipeline job' in out
    assert cli.main(['add', '--class', 'P2', '--kind', 'promote', '--telescope', 'Io', '--night', '20250102']) == 0
    assert 'added promote job 3: P2 Io 20250102' in capsys.readouterr().out
    assert cli.main(['add', '--class', 'P2', '--kind', 'sweep', '--promote', '--', 'true']) == 2
    from jobqueue.store import Store
    s = Store(os.path.join(str(tmp_path / 'q'), 'queue.sqlite'))
    p, j, k = s.get_job(1), s.get_job(2), s.get_job(3)
    s.close()
    assert '--no_T12' in p['argv'] and j['depends_on'] == 1 and j['meta']['campaign'] == 'ESO backfill'
    assert k['depends_on'] is None and '--job' not in k['argv'] and k['locks'] == ['night:Io:20250102']
    # q promote is utils/promote.py
    assert cli.main(['promote', 'list']) == 0
    assert 'no journal at' in capsys.readouterr().out


def test_backlog_files_queue_promote_jobs_too(cfg, store):
    req = {'telescope': 'Europa', 'night': NIGHT, 'download': False, 'process': True, 'replace': 'never',
           'run_targets': ['Sp0001+0001']}
    lines = []
    cli.queue_night(cfg, store, None, req, 'P2', promote=True, log=lines.append)
    assert 'queued pipeline job 1, then promote job 2' in lines[0]
    assert '--no_T12' in store.get_job(1)['argv'] and store.get_job(2)['depends_on'] == 1
    cli.queue_night(cfg, store, None, dict(req, night='20250102'), 'P2', promote=True, dry_run=True, log=lines.append)
    assert lines[1].endswith('would run the pipeline for Sp0001+0001 with --no_T12, then a promote job')


def test_the_daily_report_lists_holds_rejects_and_unfinished_promotions(cfg, store):
    clock = lambda: dt.datetime(2026, 10, 18, 6, 0)                                        # noqa: E731
    j = Journal(queue_paths(cfg)['promote_db'], now_fn=clock, check_local=False)
    held = j.add_verdict(None, 'Europa', '20250922', 'Sp0352+0210', '/v3/x', '', 'HOLD',
                         'S1 the run (queue job 9, started 2026-10-17T10:00:00Z) wrote none of its 12 files', '', {},
                         {}, 'fp', None, 'pending')
    rejected = j.add_verdict(None, 'Ganymede', '20250410', 'TOI-1227', '/v3/y', '/v2/y', 'REJECT',
                             'A1 no aperture-5 light curve (apertures 3)', '', {}, {}, 'fp', 'fp', 'pending')
    done = j.add_verdict(None, 'Io', '20251023', 'Sp0335-0829', '/v3/z', '', 'PROMOTE', '', '', {}, {}, 'fp', None,
                         'promoting')
    pid = j.start_promotion(done, None, 'Io', '20251023', 'Sp0335-0829', 'auto', lambda pid: {'stamp': 'x'})
    j.finish_promotion(pid, 'failed', error='publish step: OSError: no permission; v2 restored from /sup/x')
    j.close()
    text = build_report(dict(cfg, mode='live'), store, days=3, now=dt.datetime(2026, 10, 18, 6, 30))
    assert '## Promotion: holds and rejects' in text
    assert '2 candidate(s) wait for a decision' in text and '3 promotion item(s) wait for a person' in text
    assert '| {} | HOLD | Europa | 20250922 | Sp0352+0210 | 2026-10-18 | S1 the run'.format(held) in text
    assert '| {} | REJECT | Ganymede | 20250410 | TOI-1227 |'.format(rejected) in text
    assert '**Promotion {} did not finish** (Io 20251023 Sp0335-0829, verdict {}, failed)'.format(pid, done) in text
    assert 'q promote resume --promotion {} --apply'.format(pid) in text


def test_the_daily_report_lists_reprocessing_apart_from_decisions(cfg, store):
    clock = lambda: dt.datetime(2026, 10, 18, 6, 0)                                        # noqa: E731
    j = Journal(queue_paths(cfg)['promote_db'], now_fn=clock, check_local=False)
    lost = j.add_verdict(None, 'Callisto', '20241117', 'Sp0409-0605', '/v3/x', '/v2/x', 'HOLD',
                         'A4 light curve has 80 of 243 raw science frames (33%)', '', {}, {}, 'fp', 'fp', 'pending')
    j.add_policy_decision(lost, None, 1, 'reprocess', 'frames-lost-to-pipeline', '162 of 243 raw frames lost to '
                          'plate-solve failures on frames with stars', {})
    mixed = j.add_verdict(None, 'Callisto', '20250202', 'TOI-6714', '/v3/y', '/v2/y', 'HOLD',
                          'R1 444 finite points against 483 in v2', '', {}, {}, 'fp', 'fp', 'pending')
    j.add_policy_decision(mixed, None, 1, 'person', 'mixed-against-v2', 'mixed against v2: x0.92 the points, x1.16 '
                          'the scatter', {})
    gone = j.add_verdict(None, 'Europa', '20260523', 'Sp1746-3214', '/v3/z', '', 'REJECT',
                         'A1 no aperture-5 light curve', '', {}, {}, 'fp', None, 'pending')
    j.add_policy_decision(gone, None, 1, 'reject', 'no-aperture-5', 'no aperture-5 light curve', {})
    j.set_verdict_state(gone, 'rejected', decided_by='policy:v1:no-aperture-5', decided_at='2026-10-18T06:00:00Z',
                        decision_reason='no aperture-5 light curve')
    j.close()
    text = build_report(dict(cfg, mode='live'), store, days=3, now=dt.datetime(2026, 10, 18, 6, 30))
    assert '1 candidate(s) wait for a decision' in text
    assert '| {} | HOLD | Callisto | 20250202 | TOI-6714 | 2026-10-18 | mixed against v2'.format(mixed) in text
    assert '1 candidate(s) wait for reprocessing' in text
    assert '| {} | Callisto | 20241117 | Sp0409-0605 | frames-lost-to-pipeline | 162 of 243'.format(lost) in text
    assert 'The approval policy decided 3 verdict(s) in the last 3 days: 1 person, 1 reject, 1 reprocess.' in text
    assert '0 by the policy, 0 accepted by a person, 1 rejected' in text
    assert 'Sp1746-3214 (verdict' not in text            # policy rejects are counted, not listed one by one


def test_no_journal_no_section(cfg, store):
    text = build_report(dict(cfg, mode='live'), store, days=1, now=dt.datetime(2026, 10, 18, 6, 30))
    assert 'Promotion' not in text


@pytest.mark.parametrize('nightly', [False, True])
def test_the_nightly_switch_moves_the_watchers_runs_to_promotion(cfg, store, nightly):
    cfg = dict(cfg, telescopes={'Europa': cfg['telescopes']['Europa']}, promote=dict(cfg['promote'], nightly=nightly))
    eso = FakeEso()
    eso.objects[(PROGS['Europa'], '20261008')] = ['Sp0246+1625']
    for t in ((14, 0, 100, '2026-10-09T13:55:00.000Z'), (14, 30, 300, '2026-10-09T14:25:00.000Z'),
              (15, 0, 300, '2026-10-09T14:25:00.000Z')):
        eso.set('20261008', 'Europa', t[2], t[3])
        Watcher(cfg, store, eso, now=dt.datetime(2026, 10, 9, t[0], t[1]), log=lambda m: None).run()
    kinds = [j['kind'] for j in store.jobs()]
    pl = [j for j in store.jobs() if j['kind'] == 'pipeline'][0]
    if nightly:
        assert kinds == ['download', 'pipeline', 'promote'] and '--no_T12' in pl['argv']
        assert store.jobs(kind='promote')[0]['depends_on'] == pl['id']
    else:
        assert kinds == ['download', 'pipeline'] and '--no_T12' not in pl['argv']
