"""The dispatcher with real processes: launch, exit codes, retries, timeouts, restarts, shadow mode."""
import json
import os
import signal
import sys
import time

import pytest

from jobqueue.dispatcher import (AlreadyRunning, Dispatcher, acquire_single_instance, deployed_sha, lock_is_held,
                                 watchdog)
from jobqueue.tests.conftest import FixedResources

PY = sys.executable


def make(cfg, store, clock, **kw):
    return Dispatcher(cfg, store, resources=FixedResources(**kw), now_fn=clock, log=lambda m: None)


def add(store, argv, cls='P2', est=10.0, **kw):
    kw.setdefault('disk_heavy', False)
    jid, _ = store.add_job(cls, kw.pop('kind', 'command'), argv, cores=kw.pop('cores', 2), est_minutes=est, **kw)
    return jid


def wait_for(fn, timeout=15.0):
    end = time.time() + timeout
    while time.time() < end:
        if fn():
            return True
        time.sleep(0.1)
    return False


def poll_until(d, store, jid, states, timeout=15.0):
    def check():
        d.poll()
        return store.get_job(jid)['state'] in states
    assert wait_for(check, timeout), store.get_job(jid)
    return store.get_job(jid)


def test_runs_a_job_in_its_own_session_and_records_everything(cfg, store, clock):
    probe = ('import os, sys; print("sid-is-runner", os.getsid(0) == os.getppid()); '
             'print("cores", os.environ.get("N_CORES")); sys.exit(0)')
    jid = add(store, [PY, '-c', probe], kind='pipeline', cores=7)
    d = make(cfg, store, clock)
    d.poll()
    j = store.get_job(jid)
    assert j['state'] == 'running' and j['pid'] and j['attempts'] == 1 and j['code_sha']
    j = poll_until(d, store, jid, ('done',))
    assert j['exit_code'] == 0 and j['finished_at']
    log = open(j['log_path']).read()
    assert 'sid-is-runner True' in log          # the runner leads the job's session
    assert 'cores 7' in log                     # pipeline jobs get N_CORES = declared cores
    assert json.load(open(j['status_path']))['exit_code'] == 0
    kinds = [e['kind'] for e in store.events(job_id=jid)]
    assert kinds[:2] == ['added', 'started'] and kinds[-1] == 'finished'


def test_transient_failure_is_retried_after_the_delay(cfg, store, clock):
    jid = add(store, [PY, '-c', 'import sys; print("OSError: [Errno 28] No space left on device: /x/y"); sys.exit(1)'])
    d = make(cfg, store, clock)
    d.poll()
    j = poll_until(d, store, jid, ('queued',))
    assert j['failure_kind'] == 'transient' and 'No space left' in j['failure_signature']
    d.poll()
    assert store.get_job(jid)['state'] == 'queued'           # not before the hour is up
    clock.advance(minutes=61)
    d.poll()
    assert store.get_job(jid)['attempts'] == 2


def test_deterministic_failure_stops_and_keeps_the_signature(cfg, store, clock):
    script = ('import sys\nprint("Traceback (most recent call last):")\n'
              'print("  File \\"/opt/orchard/src/photometry/x.py\\", line 3")\n'
              'print("KeyError: \'gaia_dr3_id\'")\nsys.exit(1)\n')
    jid = add(store, [PY, '-c', script])
    d = make(cfg, store, clock)
    d.poll()
    j = poll_until(d, store, jid, ('failed',))
    assert j['failure_kind'] == 'deterministic'
    assert j['failure_signature'] == "KeyError: 'gaia_dr3_id'"


def test_job_over_three_times_its_estimate_is_killed(cfg, store, clock):
    cfg = dict(cfg, min_timeout_minutes=0, kill_grace_seconds=0)
    jid = add(store, [PY, '-c', 'import time; time.sleep(60)'], est=10.0)
    d = make(cfg, store, clock)
    d.poll()
    assert store.get_job(jid)['state'] == 'running'
    clock.advance(minutes=25)
    d.poll()
    assert store.get_job(jid)['state'] == 'running'          # 25 min < 3 x 10 min
    clock.advance(minutes=10)
    d.poll()
    j = store.get_job(jid)
    assert j['kill_reason'] == 'timeout' and j['kill_sent_at']
    j = poll_until(d, store, jid, ('failed',))
    assert j['failure_kind'] == 'timeout' and j['failure_signature'].startswith('timeout')


def test_running_job_is_never_killed_to_make_room(cfg, store, clock):
    low = add(store, [PY, '-c', 'import time; time.sleep(3)'], cls='P3', cores=70)
    d = make(cfg, store, clock)
    d.poll()
    high = add(store, [PY, '-c', 'pass'], cls='P0', cores=50)    # 70 + 50 > 112
    d.poll()
    assert store.get_job(low)['state'] == 'running' and store.get_job(low)['kill_sent_at'] is None
    assert store.get_job(high)['state'] == 'queued'
    poll_until(d, store, low, ('done',))
    poll_until(d, store, high, ('done',))


def test_cancel_with_kill_stops_a_running_job(cfg, store, clock):
    jid = add(store, [PY, '-c', 'import time; time.sleep(60)'])
    d = make(cfg, store, clock)
    d.poll()
    store.cancel(jid, kill=True)
    d.poll()
    j = poll_until(d, store, jid, ('cancelled',))
    assert j['failure_kind'] == 'cancelled'


def test_restarted_dispatcher_adopts_a_job_that_is_still_running(cfg, store, clock):
    jid = add(store, [PY, '-c', 'import time; time.sleep(1.5)'])
    make(cfg, store, clock).poll()
    assert store.get_job(jid)['state'] == 'running'
    d2 = make(cfg, store, clock)                             # the old dispatcher is gone
    d2.poll()
    assert store.get_job(jid)['state'] == 'running'
    j = poll_until(d2, store, jid, ('done',))
    assert j['exit_code'] == 0 and j['vanished'] == 0


def test_job_whose_process_vanished_is_requeued_once(cfg, store, clock):
    jid = add(store, [PY, '-c', 'import time; time.sleep(60)'])
    for expected in ('queued', 'failed'):
        make(cfg, store, clock).poll()
        j = store.get_job(jid)
        assert j['state'] == 'running'
        os.killpg(j['pid'], signal.SIGKILL)                  # like a container restart: no exit status
        d2 = make(cfg, store, clock)
        clock.advance(seconds=31)
        assert wait_for(lambda: (d2.check_running(store.get_job(jid), clock()) or True)
                        and store.get_job(jid)['state'] != 'running')
        assert store.get_job(jid)['state'] == expected
    assert store.get_job(jid)['failure_kind'] == 'vanished'


def test_dependent_job_waits_for_its_dependency(cfg, store, clock):
    marker = os.path.join(cfg['queue_root'], 'marker')
    a = add(store, [PY, '-c', 'import time; time.sleep(0.5); open({!r}, "w").write("x")'.format(marker)], cls='P0')
    b = add(store, [PY, '-c', 'import os, sys; sys.exit(0 if os.path.exists({!r}) else 3)'.format(marker)], cls='P0',
            depends_on=a)
    d = make(cfg, store, clock)
    d.poll()
    assert store.get_job(b)['state'] == 'queued'
    poll_until(d, store, a, ('done',))
    assert poll_until(d, store, b, ('done',))['exit_code'] == 0


def test_drain_starts_nothing(cfg, store, clock):
    jid = add(store, [PY, '-c', 'pass'])
    store.set_drain(True)
    d = make(cfg, store, clock)
    d.poll()
    assert store.get_job(jid)['state'] == 'queued'
    store.set_drain(False)
    d.poll()
    assert store.get_job(jid)['state'] in ('running', 'done')


def test_shadow_mode_launches_nothing_and_simulates_the_run(cfg, store, clock, tmp_path):
    marker = tmp_path / 'ran'
    cfg = dict(cfg, mode='shadow')
    jid = add(store, [PY, '-c', 'open({!r}, "w")'.format(str(marker))], est=30.0)
    other = add(store, [PY, '-c', 'pass'], est=30.0, locks=['tel:Europa'])
    blocked = add(store, [PY, '-c', 'pass'], est=30.0, locks=['tel:Europa'])
    d = make(cfg, store, clock)
    d.poll()
    j = store.get_job(jid)
    assert j['state'] == 'running' and j['shadow'] == 1 and j['pid'] is None
    assert store.get_job(blocked)['state'] == 'queued'        # shadow jobs hold their locks
    clock.advance(minutes=31)
    d.poll()
    assert store.get_job(jid)['state'] == 'done' and store.get_job(other)['state'] == 'done'
    assert not marker.exists()
    kinds = {e['kind'] for e in store.events()}
    assert {'shadow-start', 'shadow-finished', 'waiting'} <= kinds


def test_single_instance_lock_and_watchdog(cfg, store, clock):
    lock = os.path.join(cfg['queue_root'], 'dispatcher.lock')
    os.makedirs(cfg['queue_root'], exist_ok=True)
    assert watchdog(cfg, store, now=clock(), log=lambda m: None) == 3      # nobody holds the lock
    fd = acquire_single_instance(lock)
    try:
        assert lock_is_held(lock)
        with pytest.raises(AlreadyRunning):
            acquire_single_instance(lock)
        Dispatcher(cfg, store, resources=FixedResources(), now_fn=clock, log=lambda m: None).heartbeat(clock())
        assert watchdog(cfg, store, now=clock(), log=lambda m: None) == 0
        killed = []
        clock.advance(minutes=11)                                          # heartbeat now stale
        rc = watchdog(cfg, store, now=clock(), kill=lambda pid, sig: killed.append((pid, sig)), log=lambda m: None)
        assert rc == 4 and killed == [(os.getpid(), signal.SIGTERM)]
    finally:
        os.close(fd)
    assert not lock_is_held(lock)


def test_deployed_sha_prefers_version_file(tmp_path):
    (tmp_path / 'main').mkdir()
    (tmp_path / 'main' / 'ZLP_pipeline.sh').write_text('echo hi\n')
    a = deployed_sha(str(tmp_path))
    assert a.startswith('tree-')
    (tmp_path / 'main' / 'ZLP_pipeline.sh').write_text('echo bye\n')
    assert deployed_sha(str(tmp_path)) != a
    (tmp_path / 'VERSION').write_text('96eb235 2026-10-09\n')
    assert deployed_sha(str(tmp_path)) == '96eb235'


def test_command_outliving_its_runner_is_not_requeued_alongside_itself(cfg, store, clock, tmp_path):
    # the runner is SIGKILLed alone (as the OOM killer might); its command keeps running in the same group
    jid = add(store, [PY, '-c', 'import time; time.sleep(3)'])
    d = make(cfg, store, clock)
    d.poll()
    j = store.get_job(jid)
    assert wait_for(lambda: os.path.exists(j['status_path'] + '.started'))
    time.sleep(0.5)                                          # let the runner start its command
    os.kill(j['pid'], signal.SIGKILL)
    clock.advance(seconds=31)
    assert wait_for(lambda: (d.poll() or True) and bool(store.events(kinds=['runner-gone'])))
    assert store.get_job(jid)['state'] == 'running'
    # once the whole group has gone, with no exit status: vanished, requeued once (the same poll restarts it)
    assert wait_for(lambda: (d.poll() or True) and bool(store.events(kinds=['vanished-requeued'], job_id=jid)), 20)
    assert store.get_job(jid)['vanished'] == 1
