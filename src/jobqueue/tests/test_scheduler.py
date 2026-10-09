"""Admission, strict priority and locks: scheduler.plan is pure, so these build job dicts directly."""
import datetime as dt

import pytest

from jobqueue.config import DEFAULTS
from jobqueue.resources import Snapshot
from jobqueue.scheduler import plan

NOON = dt.datetime(2026, 10, 9, 10, 30)      # outside the P0 reserve window (11:00-01:00 UTC)
AFTERNOON = dt.datetime(2026, 10, 9, 14, 0)  # inside it


def snap(load5=1.0, mem=400.0, sda=5.0, sdb=5.0):
    return Snapshot(load5, mem, {'sda': sda, 'sdb': sdb}, 128)


def job(i, cls='P2', cores=20, disk_heavy=True, locks=(), est=30, telescope=None, created=None, priority=0.0,
        depends_on=None, not_before=None, state='queued'):
    return {'id': i, 'class': cls, 'cores': cores, 'disk_heavy': disk_heavy, 'locks': sorted(locks),
            'est_minutes': est, 'telescope': telescope, 'night': '20261008', 'priority': priority,
            'created_at': created or '2026-10-09T09:00:{:02d}Z'.format(i % 60), 'depends_on': depends_on,
            'not_before': not_before, 'state': state}


def started(decisions):
    return [d.job['id'] for d in decisions if d.start]


def waiting(decisions, i):
    return next(d.waiting for d in decisions if d.job['id'] == i)


@pytest.fixture
def cfg():
    c = dict(DEFAULTS)
    c['p0_guard'] = {'enabled': False}
    return c


# ---------------------------------------------------------------- admission


def test_declared_cores_capped_at_112(cfg):
    queued = [job(i) for i in range(1, 7)]
    d = plan(queued, [], snap(), dict(cfg, max_disk_heavy_starts_per_poll=10), NOON)
    assert started(d) == [1, 2, 3, 4, 5]          # 100 cores; a sixth would make 120
    assert 'cores' in waiting(d, 6)               # (and the 100 cores just started count as load)


def test_running_jobs_count_against_the_cap(cfg):
    running = [job(90 + i, state='running') for i in range(5)]
    d = plan([job(1, cores=12, disk_heavy=False), job(2, cores=13, disk_heavy=False)], running, snap(), cfg, NOON)
    assert started(d) == [1]
    assert waiting(d, 2) == ['cores']


def test_p0_reserve_in_the_nightly_window(cfg):
    queued = [job(i, disk_heavy=False) for i in range(1, 6)]
    d = plan(queued, [], snap(), cfg, AFTERNOON)
    assert started(d) == [1, 2, 3]                 # 112 - 40 reserved = 72 cores for P1-P3
    d = plan(queued, [], snap(), cfg, NOON)
    assert started(d) == [1, 2, 3, 4, 5]
    p0 = [job(10 + i, cls='P0') for i in range(5)]
    d = plan(p0, [], snap(), cfg, AFTERNOON)
    assert started(d) == [10, 11, 12, 13, 14]      # P0 may use the reserve


def test_load_average_blocks_everything_including_p0(cfg):
    d = plan([job(1, cls='P0'), job(2)], [], snap(load5=100.0), cfg, NOON)
    assert started(d) == []
    assert 'load' in waiting(d, 1)


def test_load_counts_cores_started_in_the_same_pass(cfg):
    d = plan([job(1, cls='P0'), job(2, cls='P0')], [], snap(load5=85.0), cfg, NOON)
    assert started(d) == [1]
    assert waiting(d, 2) == ['load']


def test_memory_floor(cfg):
    d = plan([job(1, cls='P0'), job(2)], [], snap(mem=50.0), cfg, NOON)
    assert started(d) == []
    assert 'mem' in waiting(d, 1)


def test_disk_util_checked_except_for_p0(cfg):
    q = [job(1, cls='P0'), job(2, cls='P2'), job(3, cls='P2', disk_heavy=False)]
    d = plan(q, [], snap(sdb=80.0), cfg, NOON)
    assert started(d) == [1, 3]
    assert waiting(d, 2) == ['disk']


def test_unknown_disk_util_counts_as_busy(cfg):
    d = plan([job(1)], [], Snapshot(1.0, 400.0, {'sda': None, 'sdb': 5.0}, 128), cfg, NOON)
    assert waiting(d, 1) == ['disk']


def test_one_disk_heavy_start_per_poll(cfg):
    d = plan([job(1), job(2), job(3, disk_heavy=False, cores=1)], [], snap(), cfg, NOON)
    assert started(d) == [1, 3]
    assert waiting(d, 2) == ['disk']


# ---------------------------------------------------------------- strict priority


def test_lower_class_never_takes_what_a_higher_class_waits_for(cfg):
    running = [job(90 + i, state='running', disk_heavy=False) for i in range(5)]   # 100 cores busy
    q = [job(1, cls='P0', cores=20), job(2, cls='P1', cores=10, disk_heavy=False),
         job(3, cls='P0', cores=2, disk_heavy=False)]
    d = plan(q, running, snap(), cfg, NOON)
    assert waiting(d, 1) == ['cores']
    assert 'cores' in waiting(d, 2)            # 10 cores would fit, but P0 is waiting for cores
    assert started(d) == [3]                   # same class may use what is left (a download)


def test_lower_class_may_use_a_different_resource(cfg):
    running = [job(90, state='running', telescope='Callisto', locks=['tel:Callisto'])]
    q = [job(1, cls='P0', telescope='Callisto', locks=['tel:Callisto', 'night:Callisto:20261008']),
         job(2, cls='P1', telescope='Europa', locks=['tel:Europa', 'night:Europa:20261001'])]
    d = plan(q, running, snap(), cfg, NOON)
    assert waiting(d, 1) == ['tel:Callisto']
    assert started(d) == [2]


def test_higher_class_waiting_for_disk_blocks_lower_disk_jobs_only(cfg):
    q = [job(1, cls='P1'), job(2, cls='P2'), job(3, cls='P3', disk_heavy=False, cores=2)]
    d = plan(q, [], snap(sda=75.0), cfg, NOON)
    assert waiting(d, 1) == ['disk']
    assert waiting(d, 2) == ['disk']
    assert started(d) == [3]


def test_classes_go_in_order_then_priority_then_age(cfg):
    q = [job(1, cls='P3', cores=60, disk_heavy=False), job(2, cls='P2', cores=60, disk_heavy=False),
         job(3, cls='P2', cores=60, disk_heavy=False, priority=5.0)]
    d = plan(q, [], snap(), cfg, NOON)
    assert [x.job['id'] for x in d] == [3, 2, 1]
    assert started(d) == [3]


def test_nothing_running_is_ever_preempted(cfg):
    running = [job(90 + i, cls='P3', state='running', disk_heavy=False) for i in range(5)]
    d = plan([job(1, cls='P0')], running, snap(), cfg, NOON)
    assert started(d) == []                    # the P0 job waits; plan has no way to stop anything
    assert all(not hasattr(x, 'kill') for x in d)


# ---------------------------------------------------------------- locks


def test_one_pipeline_job_per_telescope(cfg):
    q = [job(1, telescope='Europa', locks=['tel:Europa', 'night:Europa:20261001']),
         job(2, telescope='Europa', locks=['tel:Europa', 'night:Europa:20261002'], disk_heavy=False)]
    d = plan(q, [], snap(), cfg, NOON)
    assert started(d) == [1]
    assert waiting(d, 2) == ['tel:Europa']


def test_telescope_night_is_exclusive(cfg):
    q = [job(1, cls='P0', telescope='Europa', locks=['night:Europa:20261008', 'sem:eso'], disk_heavy=False, cores=2),
         job(2, cls='P1', telescope='Europa', locks=['night:Europa:20261008'], disk_heavy=False)]
    d = plan(q, [], snap(), cfg, NOON)
    assert started(d) == [1]
    assert waiting(d, 2) == ['night:Europa:20261008']


def test_target_is_exclusive_across_telescopes(cfg):
    q = [job(1, cls='P0', telescope='Europa', locks=['tel:Europa', 'target:SP0001+0001']),
         job(2, cls='P0', telescope='Ganymede', locks=['tel:Ganymede', 'target:SP0001+0001']),
         job(3, cls='P0', telescope='Callisto', locks=['tel:Callisto', 'target:SP0002+0002'])]
    d = plan(q, [], snap(), cfg, NOON)
    assert started(d) == [1, 3]
    assert waiting(d, 2) == ['target:SP0001+0001']


def test_processes_the_queue_did_not_start_hold_locks(cfg):
    q = [job(1, cls='P0', telescope='Europa', locks=['tel:Europa', 'night:Europa:20261008'])]
    d = plan(q, [], snap(), cfg, NOON, external={'tel:Europa': 'pid 4 ZLP_pipeline.sh ... Europa'})
    assert waiting(d, 1) == ['tel:Europa']


def test_semaphore_counts_holders(cfg):
    q = [job(i, cls='P0', cores=2, disk_heavy=False, locks=['sem:eso', 'night:T{}:1'.format(i)]) for i in range(1, 4)]
    d = plan(q, [], snap(), cfg, NOON)
    assert started(d) == [1, 2]                # semaphores.eso = 2
    assert waiting(d, 3) == ['sem:eso']


# ---------------------------------------------------------------- eligibility


def test_dependency_must_be_settled(cfg):
    q = [job(2, cls='P0', depends_on=1)]
    assert plan(q, [], snap(), cfg, NOON, dep_states={1: 'running'}) == []
    assert started(plan(q, [], snap(), cfg, NOON, dep_states={1: 'done'})) == [2]
    assert started(plan(q, [], snap(), cfg, NOON, dep_states={1: 'failed'})) == [2]


def test_retry_backoff_pause_and_drain(cfg):
    q = [job(1, not_before='2026-10-09T11:00:00Z', disk_heavy=False), job(2, cls='P3', disk_heavy=False)]
    assert plan(q, [], snap(), cfg, NOON, paused={'P3'}) == []
    assert started(plan(q, [], snap(), cfg, NOON + dt.timedelta(hours=1), paused={'P3'})) == [1]
    assert plan(q, [], snap(), cfg, NOON + dt.timedelta(hours=1), draining=True) == []


def test_p0_guard_keeps_backlog_off_a_telescope_due_tonight(cfg):
    cfg = dict(cfg, p0_guard={'enabled': True, 'windows_utc': {'Callisto': ['13:00', '17:00']}})
    j = job(1, telescope='Callisto', locks=['tel:Callisto'], est=60, disk_heavy=False)
    t = dt.datetime(2026, 10, 9, 12, 30)
    assert waiting(plan([j], [], snap(), cfg, t), 1) == ['p0-guard:Callisto']
    assert started(plan([j], [], snap(), cfg, t, p0_settled={'Callisto'})) == [1]
    short = dict(j, est_minutes=20)
    assert started(plan([short], [], snap(), cfg, t)) == [1]
    assert started(plan([j], [], snap(), cfg, dt.datetime(2026, 10, 9, 17, 30))) == [1]
