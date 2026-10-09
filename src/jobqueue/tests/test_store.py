"""Job state transitions in the SQLite store."""
import pytest

from jobqueue.util import parse_iso


def add(store, **kw):
    kw.setdefault('telescope', 'Europa')
    kw.setdefault('night', '20261008')
    jid, created = store.add_job(kw.pop('cls', 'P2'), kw.pop('kind', 'pipeline'), kw.pop('argv', ['true']), **kw)
    return jid


def run(store, jid):
    assert store.mark_running(jid, log_path='/tmp/x.log', status_path='/tmp/x.status', code_sha='abc')
    return store.get_job(jid)


def test_add_is_queued_and_dedupes_while_active(store):
    a, created = store.add_job('P0', 'pipeline', ['x'], dedupe_key='pipeline:Europa:20261008:*')
    assert created and store.get_job(a)['state'] == 'queued'
    b, created = store.add_job('P0', 'pipeline', ['x'], dedupe_key='pipeline:Europa:20261008:*')
    assert (b, created) == (a, False)
    run(store, a)
    store.finish(a, 0)
    c, created = store.add_job('P0', 'pipeline', ['x'], dedupe_key='pipeline:Europa:20261008:*')
    assert created and c != a


def test_bad_class_rejected(store):
    with pytest.raises(ValueError):
        store.add_job('P9', 'pipeline', ['x'])


def test_start_records_attempt_and_only_from_queued(store):
    jid = add(store)
    j = run(store, jid)
    assert j['state'] == 'running' and j['attempts'] == 1 and j['code_sha'] == 'abc'
    assert not store.mark_running(jid)


def test_success(store):
    jid = add(store)
    run(store, jid)
    j = store.finish(jid, 0)
    assert j['state'] == 'done' and j['exit_code'] == 0 and j['failure_signature'] is None


def test_transient_failure_retried_twice_after_an_hour_then_failed(store, clock):
    jid = add(store, max_retries=2)
    for attempt in (1, 2):
        run(store, jid)
        j = store.finish(jid, 1, signature='OSError: [Errno #] No space left on device', transient=True,
                         retry_delay_minutes=60)
        assert j['state'] == 'queued'
        assert parse_iso(j['not_before']) == clock() + __import__('datetime').timedelta(minutes=60)
        clock.advance(minutes=61)
    run(store, jid)
    j = store.finish(jid, 1, signature='OSError: no space', transient=True)
    assert j['state'] == 'failed' and j['attempts'] == 3 and j['failure_kind'] == 'transient'


def test_deterministic_failure_stops_with_signature(store):
    jid = add(store)
    run(store, jid)
    j = store.finish(jid, 1, signature="KeyError: 'gaia_dr3_id'")
    assert j['state'] == 'failed' and j['failure_kind'] == 'deterministic'
    assert j['failure_signature'] == "KeyError: 'gaia_dr3_id'"


def test_vanished_job_requeued_once(store):
    jid = add(store)
    run(store, jid)
    assert store.vanish(jid)['state'] == 'queued'
    run(store, jid)
    j = store.vanish(jid)
    assert j['state'] == 'failed' and j['failure_kind'] == 'vanished'


def test_vanished_restart_does_not_use_up_a_transient_retry(store):
    jid = add(store, max_retries=1)
    run(store, jid)
    store.vanish(jid)
    run(store, jid)
    assert store.finish(jid, 1, signature='x', transient=True)['state'] == 'queued'


def test_cancel(store):
    a = add(store, kind='download')
    b = add(store, depends_on=a)
    assert store.cancel(a) == 'cancelled'
    assert store.get_job(b)['state'] == 'cancelled'       # dependents go with it
    c = add(store)
    run(store, c)
    with pytest.raises(ValueError):
        store.cancel(c)
    assert store.cancel(c, kill=True) == 'kill requested'
    assert store.get_job(c)['cancel_requested'] == 1
    assert store.finish(c, -15, failure_kind='cancelled')['state'] == 'cancelled'


def test_failed_dependency_cancels_dependents_that_need_success(store):
    a = add(store, kind='fetch')
    b = add(store, depends_on=a, dep_requires_success=True)
    c = add(store, depends_on=a)
    run(store, a)
    store.finish(a, 1, signature='boom')
    assert store.get_job(b)['state'] == 'cancelled'
    assert store.get_job(c)['state'] == 'queued'          # e.g. the pipeline after a failed download still runs


def test_requeue_by_signature(store):
    ids = []
    for sig in ("KeyError: 'gaia_dr3_id'", "KeyError: 'gaia_dr3_id'", 'ValueError: broadcast'):
        jid = add(store, night=str(20261000 + len(ids)))
        run(store, jid)
        store.finish(jid, 1, signature=sig)
        ids.append(jid)
    assert [i for i, _ in store.requeue(signature='gaia_dr3_id', dry_run=True)] == ids[:2]
    assert store.get_job(ids[0])['state'] == 'failed'
    out = store.requeue(signature='gaia_dr3_id')
    assert [o for _, o in out] == ['requeued', 'requeued']
    j = store.get_job(ids[0])
    assert j['state'] == 'queued' and j['attempts'] == 0 and j['failure_signature'] is None
    assert store.get_job(ids[2])['state'] == 'failed'
    assert any(e['kind'] == 'requeued' for e in store.events(job_id=ids[0]))


def test_requeue_skips_a_job_already_active_again(store):
    a, _ = store.add_job('P1', 'pipeline', ['x'], dedupe_key='k')
    run(store, a)
    store.finish(a, 1, signature='boom')
    store.add_job('P1', 'pipeline', ['x'], dedupe_key='k')
    assert store.requeue(job_ids=[a])[0][1].startswith('skipped')


def test_pause_resume_drain(store):
    store.pause('P2', 'testing')
    assert store.paused() == {'P2'}
    store.pause('all')
    assert store.paused() == {'P0', 'P1', 'P2', 'P3'}
    store.resume('all')
    assert store.paused() == set()
    store.set_drain(True)
    assert store.draining()
    store.set_drain(False)
    assert not store.draining()


def test_every_transition_is_an_event(store):
    jid = add(store)
    run(store, jid)
    store.finish(jid, 0)
    assert [e['kind'] for e in store.events(job_id=jid)] == ['added', 'started', 'finished']
