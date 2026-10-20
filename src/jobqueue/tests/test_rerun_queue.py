"""The queue's side of reruns from a late stage (docs/rerun-from-stage.md): `q add --from T10` queues a restage job,
the pipeline job with --from_T10 --no_T12 after it, and a promote job; a restage outside the queue holds its night."""
import datetime as dt
import os

import pytest

from jobqueue import cli, jobspec
from jobqueue.resources import parse_promote_argv
from jobqueue.scheduler import eligible

NIGHT = '20010101'
TARGET = 'Sp0451+0305'


def promoted_night(cfg, tel='Europa', targets=(TARGET,)):
    """A night promotion has moved into v2: each target's *_output.fits, and the reduction's image lists."""
    base = cfg['basedir']
    for t in targets:
        d = os.path.join(base, 'PipelineOutput', 'v2', tel, 'output', NIGHT, t)
        os.makedirs(d)
        open(os.path.join(d, '4123_I+z_{}_output.fits'.format(NIGHT)), 'w').close()
    red = os.path.join(base, 'PipelineOutput', 'v2', tel, 'output', NIGHT, 'reduction')
    os.makedirs(red)
    for t in targets:
        open(os.path.join(red, '1_image_{}.list'.format(t)), 'w').close()


def add(cfg, *args):
    a = cli.build_parser().parse_args(['add', '--class', 'P3'] + list(args))
    from jobqueue.store import Store
    s = Store(os.path.join(cfg['queue_root'], 'queue.sqlite'), check_local=False)
    try:
        return cli.cmd_add(cfg, s, a), s.jobs()
    finally:
        s.close()


def test_a_rerun_from_t10_is_a_restage_then_the_pipeline_then_the_promotion(cfg, store):
    promoted_night(cfg)
    spec = jobspec.pipeline_job(cfg, 'P3', 'Europa', NIGHT, run_targets=[TARGET], no_t12=True,
                                extra_flags=['--from_T10'])
    rid, pl, created, pr = jobspec.add_rerun(store, cfg, spec, 10, promote=True, campaign='PWV outage reruns')
    r, p, q = store.get_job(rid), store.get_job(pl), store.get_job(pr)
    assert r['kind'] == 'restage' and r['argv'][1:] == ['utils/restage.py', 'stage', '--apply', '--telescope',
                                                         'Europa', '--night', NIGHT, '--from', '10', '--target', TARGET]
    assert r['locks'] == ['night:Europa:{}'.format(NIGHT)] and r['cores'] == 1 and not r['disk_heavy']
    assert p['depends_on'] == rid and p['dep_requires_success']
    assert '--from_T10' in p['argv'] and '--no_T12' in p['argv'] and p['argv'][-1] == TARGET
    assert q['kind'] == 'promote' and q['depends_on'] == pl
    assert (r['meta'] or {}).get('campaign') == (p['meta'] or {}).get('campaign') == 'PWV outage reruns'
    # a normal run of the night is another job: the rerun's dedupe key is its own
    other, created2, _ = jobspec.add_pipeline(store, cfg, jobspec.pipeline_job(cfg, 'P3', 'Europa', NIGHT,
                                                                               run_targets=[TARGET]))
    assert created2 and other != pl
    # the pipeline waits for the restage
    now = dt.datetime(2026, 10, 20, 12)
    assert eligible([p], now, set(), {rid: 'queued'}) == []
    assert eligible([p], now, set(), {rid: 'done'}) != []


def test_a_rerun_must_not_run_t12(cfg, store):
    with pytest.raises(ValueError):
        jobspec.add_rerun(store, cfg, jobspec.pipeline_job(cfg, 'P3', 'Europa', NIGHT, extra_flags=['--from_T10']), 10)
    with pytest.raises(ValueError):
        jobspec.add_rerun(store, cfg, jobspec.pipeline_job(cfg, 'P3', 'Europa', NIGHT, no_t12=True), 10)


def test_q_add_from_t10_queues_the_restaged_targets(cfg, capsys):
    promoted_night(cfg, targets=(TARGET, 'Sp2209-3707'))
    rc, jobs = add(cfg, '--telescope', 'Europa', '--night', NIGHT, '--from', 'pwv', '--promote',
                   '--campaign', 'PWV outage reruns')
    assert rc == 0
    kinds = [j['kind'] for j in jobs]
    assert kinds == ['restage', 'pipeline', 'promote']
    pl = jobs[1]
    assert pl['targets'] == '{} Sp2209-3707'.format(TARGET) and '--from_T10' in pl['argv']
    assert pl['est_minutes'] >= cfg['min_estimate_minutes']
    assert pl['cores'] == min(jobspec.RERUN_CORES[10], cfg['pipeline_cores'])
    assert not jobs[0]['disk_heavy']
    assert 'restage job' in capsys.readouterr().out


def test_a_rerun_from_t9_copies_frames_so_its_restage_is_disk_heavy(cfg, store):
    spec = jobspec.pipeline_job(cfg, 'P3', 'Europa', NIGHT, run_targets=[TARGET], no_t12=True,
                                extra_flags=['--from_T9'])
    rid, pl, _, _ = jobspec.add_rerun(store, cfg, spec, 9, promote=True)
    r = store.get_job(rid)
    assert r['disk_heavy'] and r['argv'][-3:] == ['9', '--target', TARGET]
    assert '--from_T9' in store.get_job(pl)['argv']


def test_q_add_from_refuses_what_it_cannot_do(cfg, capsys):
    promoted_night(cfg)
    for args, why in ((['--from', 'T8', '--promote'], 'not supported'),
                      (['--from', 'T10'], 'needs --promote'),
                      (['--from', 'T10', '--promote', '--download'], 'does not go with --download')):
        rc, jobs = add(cfg, '--telescope', 'Europa', '--night', NIGHT, *args)
        assert rc == 2 and jobs == [] and why in capsys.readouterr().out
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(['add', '--class', 'P3', '--from', 'gremlins'])
    rc, jobs = add(cfg, '--telescope', 'Europa', '--night', '20010102', '--from', 'T10', '--promote')
    assert rc == 0 and jobs == [] and 'no target of the night' in capsys.readouterr().out


def test_q_add_from_dry_run_queues_nothing(cfg, capsys):
    promoted_night(cfg)
    rc, jobs = add(cfg, '--telescope', 'Europa', '--night', NIGHT, '--from', 'T10', '--promote', '--dry-run')
    assert rc == 0 and jobs == []
    out = capsys.readouterr().out
    assert 'would restage' in out and 'would be restaged from v2' in out


def test_q_add_from_a_file_of_nights(cfg, tmp_path):
    promoted_night(cfg)
    f = tmp_path / 'nights.txt'
    f.write_text('Europa,{},0,0,1,{}\n'.format(NIGHT, TARGET))
    rc, jobs = add(cfg, '--from-file', str(f), '--from', 'T10', '--promote')
    assert rc == 0 and [j['kind'] for j in jobs] == ['restage', 'pipeline', 'promote']
    f.write_text('Europa,{},1,0,1\n'.format(NIGHT))
    rc, jobs = add(cfg, '--from-file', str(f), '--from', 'T10', '--promote')
    assert rc == 2


def test_a_restage_outside_the_queue_holds_its_night():
    assert parse_promote_argv(['python', 'utils/restage.py', 'stage', '--apply', '--telescope', 'Io', '--night',
                               NIGHT, '--from', '10']) == {'night:Io:{}'.format(NIGHT)}
    assert parse_promote_argv(['python', 'utils/restage.py', 'stage', '--telescope', 'Io', '--night', NIGHT]) is None
