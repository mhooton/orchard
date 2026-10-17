"""Code releases: building, pinning jobs at launch, activating, and bringing the live tree up to date."""
import json
import os
import shutil
import stat
import sys
import time

import pytest

from jobqueue import cli
from jobqueue.config import load_config
from jobqueue.dispatcher import Dispatcher
from jobqueue.releases import (SYNC_KEY, ReleaseError, activate, build, code_dir, current, current_id, live_diff,
                               live_tree_users, pin, read_release, releases_dir, sync_live)
from jobqueue.tests.conftest import Clock, FixedResources

PY = sys.executable
JOBQUEUE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def put_jobqueue(src):
    """The real package: a job's runner is `python -m jobqueue.runner` from the tree it is pinned to."""
    shutil.copytree(JOBQUEUE, os.path.join(src, 'jobqueue'), ignore=shutil.ignore_patterns('tests', '__pycache__'),
                    dirs_exist_ok=True)


def write(path, text, mode=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(text)
    if mode is not None:
        os.chmod(path, mode)


def read(path):
    with open(path) as f:
        return f.read()


@pytest.fixture
def live(cfg, monkeypatch):
    """A small production tree: code, a deployed-only secret and data file, and the clutter that lives there.
    As in orchard-server, ORCHARD_PATH and PYTHONPATH name it."""
    src = cfg['src_dir']
    monkeypatch.setenv('ORCHARD_PATH', src)
    monkeypatch.setenv('PYTHONPATH', os.pathsep.join([src, os.environ.get('PYTHONPATH', '')]))
    write(os.path.join(src, 'main', 'ZLP_pipeline.sh'), '#!/bin/bash\necho live\n', 0o755)
    put_jobqueue(src)
    write(os.path.join(src, 'photometry', 'x.py'), 'VERSION = 1\n')
    write(os.path.join(src, 'marker.txt'), 'old\n')
    write(os.path.join(src, 'reporting', '.env'), 'ESO_PASSWORD=x\n', 0o600)
    write(os.path.join(src, 'utils', 'ml_40pc_v2.csv'), 'a,b\n')
    write(os.path.join(src, 'catcache', 'cch_00000001'), 'cache')
    write(os.path.join(src, 'photometry', '__pycache__', 'x.cpython-39.pyc'), 'pyc')
    write(os.path.join(src, 'core.113496'), 'core dump')
    write(os.path.join(src, 'observation_summary.log'), 'log')
    write(os.path.join(src, '1', 'PipelineOutput', 'v3', 'junk'), 'junk')
    write(os.path.join(src, 'utils', '1', 'keep.py'), 'kept: only the top-level 1/ is clutter\n')
    return src


def stage_tree(tmp_path, name, files):
    root = tmp_path / name / 'src'
    base = {'main/ZLP_pipeline.sh': '#!/bin/bash\necho live\n'}
    base.update(files)
    for rel, text in base.items():
        write(str(root / rel), text, 0o755 if rel.endswith('.sh') else None)
    put_jobqueue(str(root))
    return str(root)


def test_build_copies_the_live_tree_lays_the_commit_over_it_and_freezes_it(cfg, live, tmp_path):
    stage = stage_tree(tmp_path, 'stage1', {'photometry/x.py': 'VERSION = 2\n', 'photometry/new.py': 'NEW = 1\n',
                                             'utils/ml_40pc_v2.csv': 'a,b\n'})
    os.chmod(os.path.join(stage, 'utils', 'ml_40pc_v2.csv'), 0o755)       # git's mode, not the live tree's
    rel = build(cfg, stage=stage, sha='70f3708abcdef', log=lambda m: None)
    src = rel['src']
    assert rel['id'].endswith('-70f3708') and rel['sha'] == '70f3708abcdef'
    assert read(os.path.join(src, 'photometry', 'x.py')) == 'VERSION = 2\n'        # the commit wins
    assert read(os.path.join(src, 'photometry', 'new.py')) == 'NEW = 1\n'
    assert read(os.path.join(src, 'reporting', '.env')) == 'ESO_PASSWORD=x\n'      # deployed-only files come along
    assert read(os.path.join(src, 'utils', 'ml_40pc_v2.csv')) == 'a,b\n'
    assert read(os.path.join(src, 'utils', '1', 'keep.py')).startswith('kept')
    for clutter in ('catcache', 'core.113496', 'observation_summary.log', '1',
                    os.path.join('photometry', '__pycache__')):
        assert not os.path.exists(os.path.join(src, clutter)), clutter
    assert read(os.path.join(src, 'VERSION')).split()[0] == rel['id']
    assert sorted(rel['differs_from_live']) == ['photometry/new.py', 'photometry/x.py']   # a mode alone is no change
    for rel_path in ('photometry/x.py', 'main/ZLP_pipeline.sh', 'reporting/.env'):
        assert not os.stat(os.path.join(src, rel_path)).st_mode & 0o222, rel_path   # read-only
    assert os.stat(os.path.join(src, 'main', 'ZLP_pipeline.sh')).st_mode & 0o100     # still executable
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 1\n'          # building changes nothing live
    assert current_id(cfg) is None                                                    # nor what jobs run


def test_build_without_a_stage_snapshots_the_live_tree(cfg, live):
    rel = build(cfg, log=lambda m: None)
    assert rel['id'].endswith('-live') and rel['differs_from_live'] == []
    assert read(os.path.join(rel['src'], 'photometry', 'x.py')) == 'VERSION = 1\n'


def test_two_builds_in_the_same_second_get_different_ids(cfg, live, clock):
    a = build(cfg, now=clock(), log=lambda m: None)
    b = build(cfg, now=clock(), log=lambda m: None)
    assert a['id'] != b['id'] and b['id'] == a['id'] + '.2'


def test_build_refuses_a_stage_that_is_not_a_src_tree(cfg, live, tmp_path):
    (tmp_path / 'notsrc').mkdir()
    with pytest.raises(ReleaseError):
        build(cfg, stage=str(tmp_path / 'notsrc'), sha='abc', log=lambda m: None)
    stage = stage_tree(tmp_path, 'stage', {})
    with pytest.raises(ReleaseError):
        build(cfg, stage=stage, sha=None, log=lambda m: None)
    assert not os.path.exists(releases_dir(cfg)) or not os.listdir(releases_dir(cfg))


def test_next_release_links_unchanged_files_and_leaves_the_old_one_untouched(cfg, live, tmp_path, monkeypatch):
    r1 = build(cfg, stage=stage_tree(tmp_path, 's1', {'photometry/x.py': 'VERSION = 2\n'}), sha='aaaaaaa',
               log=lambda m: None)
    activate(cfg, r1['id'], log=lambda m: None)
    r2 = build(cfg, stage=stage_tree(tmp_path, 's2', {'photometry/x.py': 'VERSION = 3\n'}), sha='bbbbbbb',
               log=lambda m: None)
    assert r1['id'] != r2['id'] and r2['previous'] == r1['id']
    same = os.path.join('utils', 'ml_40pc_v2.csv')
    assert os.stat(os.path.join(r1['src'], same)).st_ino == os.stat(os.path.join(r2['src'], same)).st_ino
    x1, x2 = (os.path.join(r['src'], 'photometry', 'x.py') for r in (r1, r2))
    assert os.stat(x1).st_ino != os.stat(x2).st_ino
    assert read(x1) == 'VERSION = 2\n' and read(x2) == 'VERSION = 3\n'
    assert r2['changed_since_previous'] == ['photometry/x.py']


def test_activate_repoints_current_atomically_and_records_it(cfg, live, store):
    r1 = build(cfg, log=lambda m: None)
    assert activate(cfg, r1['id'], store=store, log=lambda m: None) is None
    link = os.path.join(releases_dir(cfg), 'current')
    assert os.readlink(link) == r1['id']                       # relative: resolves on the host and in the container
    assert current(cfg)['id'] == r1['id']
    with pytest.raises(ReleaseError):
        activate(cfg, 'no-such-release', log=lambda m: None)
    assert current_id(cfg) == r1['id']
    ev = store.events(kinds=['release'])
    assert ev and ev[-1]['message'].endswith(r1['id'])
    with open(os.path.join(releases_dir(cfg), 'history.jsonl')) as f:
        assert json.loads(f.readlines()[-1])['to'] == r1['id']


def test_pin_points_a_job_at_the_release_and_leaves_other_paths_alone(cfg, live, tmp_path):
    env = {'ORCHARD_PATH': cfg['src_dir'], 'PYTHONPATH': os.pathsep.join(['/scratch/pylib', cfg['src_dir']])}
    assert pin(cfg, cfg['src_dir'], env) == (cfg['src_dir'], env, None)        # no release: as before
    r1 = build(cfg, log=lambda m: None)
    activate(cfg, r1['id'], log=lambda m: None)
    cwd, out, rid = pin(cfg, cfg['src_dir'], env)
    assert rid == r1['id'] and cwd == r1['src']
    assert out['ORCHARD_PATH'] == r1['src'] and out['JOBQUEUE_SRC'] == r1['src']
    assert out['PYTHONPATH'].split(os.pathsep) == ['/scratch/pylib', r1['src']]
    assert pin(cfg, os.path.join(cfg['src_dir'], 'photometry'), env)[0] == os.path.join(r1['src'], 'photometry')
    scratch = str(tmp_path)
    assert pin(cfg, scratch, {'ORCHARD_PATH': '/scratch/src'})[0] == scratch
    assert pin(cfg, scratch, {'ORCHARD_PATH': '/scratch/src'})[1]['ORCHARD_PATH'] == '/scratch/src'
    r2 = build(cfg, log=lambda m: None)
    activate(cfg, r2['id'], log=lambda m: None)
    assert pin(cfg, r1['src'], {})[0] == r2['src']             # a path into an older release moves to the current one


def test_code_dir_is_the_jobs_release(cfg, monkeypatch):
    monkeypatch.delenv('JOBQUEUE_SRC', raising=False)
    assert code_dir(cfg) == cfg['src_dir']
    monkeypatch.setenv('JOBQUEUE_SRC', '/r/src')
    assert code_dir(cfg) == '/r/src'


def wait_for(fn, timeout=15.0):
    end = time.time() + timeout
    while time.time() < end:
        if fn():
            return True
        time.sleep(0.1)
    return False


def test_running_job_keeps_its_release_while_the_next_job_gets_the_new_one(cfg, live, store, tmp_path):
    """The point of releases: deploy now; a job already running finishes on the code it started with."""
    go = str(tmp_path / 'go')
    probe = ('import os, time\n'
             'while not os.path.exists({!r}): time.sleep(0.05)\n'
             'print("marker", open("marker.txt").read().strip())\n'
             'print("orchard", os.environ["ORCHARD_PATH"])\n').format(go)
    r1 = build(cfg, log=lambda m: None)
    activate(cfg, r1['id'], log=lambda m: None)
    clock = Clock(store.now())
    d = Dispatcher(cfg, store, resources=FixedResources(), now_fn=clock, log=lambda m: None)
    a, _ = store.add_job('P2', 'command', [PY, '-c', probe], cwd=cfg['src_dir'], cores=1, est_minutes=5,
                         disk_heavy=False)
    d.poll()
    assert store.get_job(a)['state'] == 'running' and store.get_job(a)['code_sha'] == r1['id']

    r2 = build(cfg, stage=stage_tree(tmp_path, 's2', {'marker.txt': 'new\n'}), sha='cafe123', log=lambda m: None)
    activate(cfg, r2['id'], log=lambda m: None)
    assert sync_live(cfg, r2['id'], force=True, log=lambda m: None)['written'] == ['VERSION', 'marker.txt']
    b, _ = store.add_job('P2', 'command', [PY, '-c', probe], cwd=cfg['src_dir'], cores=1, est_minutes=5,
                         disk_heavy=False)
    d.poll()
    assert store.get_job(b)['code_sha'] == r2['id']
    open(go, 'w').close()

    def done():
        d.poll()
        return all(store.get_job(j)['state'] == 'done' for j in (a, b))
    assert wait_for(done)
    log_a, log_b = (read(store.get_job(j)['log_path']) for j in (a, b))
    assert 'marker old' in log_a and 'orchard ' + r1['src'] in log_a
    assert 'marker new' in log_b and 'orchard ' + r2['src'] in log_b
    assert 'ORCHARD_PATH ' + r1['src'] in log_a                  # the runner's header says where it ran from


def test_without_a_release_jobs_run_from_the_live_tree(cfg, live, store):
    clock = Clock(store.now())
    d = Dispatcher(cfg, store, resources=FixedResources(), now_fn=clock, log=lambda m: None)
    jid, _ = store.add_job('P2', 'command', [PY, '-c', 'print(open("marker.txt").read())'], cwd=cfg['src_dir'],
                           cores=1, est_minutes=5, disk_heavy=False)
    d.poll()
    assert wait_for(lambda: d.poll() is not None and store.get_job(jid)['state'] == 'done')
    j = store.get_job(jid)
    assert 'old' in read(j['log_path']) and not j['code_sha'].startswith('20')


def fake_proc(tmp_path, procs):
    root = tmp_path / 'proc'
    for pid, (argv, env) in procs.items():
        d = root / str(pid)
        d.mkdir(parents=True)
        (d / 'cmdline').write_bytes(b'\0'.join(a.encode() for a in argv) + b'\0')
        if env is not None:
            (d / 'environ').write_bytes(b'\0'.join('{}={}'.format(k, v).encode() for k, v in env.items()) + b'\0')
    return str(root)


def test_live_tree_users_are_unpinned_pipeline_download_and_job_processes(cfg, live, tmp_path):
    rel = os.path.join(releases_dir(cfg), 'r1', 'src')
    zlp = ['bash', './main/ZLP_pipeline.sh', '1', '/data', '20240408', '8', '2', 'Callisto']
    proc = fake_proc(tmp_path, {
        101: (zlp, {'ORCHARD_PATH': cfg['src_dir']}),                                   # the old way: live tree
        102: (zlp, {'ORCHARD_PATH': rel}),                                              # pinned
        103: (['python', 'download/SSO_download.py', '--telescope', 'Io'], None),      # unreadable: counted
        104: (['python', '-m', 'jobqueue', 'jobs', 'fetch'], {'ORCHARD_PATH': cfg['src_dir'] + '/'}),
        105: (['python', '-m', 'jobqueue', 'dispatcher'], {'ORCHARD_PATH': cfg['src_dir']}),  # not pipeline work
        106: (['python', '-m', 'jobqueue.runner', '--log', 'x'], {'ORCHARD_PATH': '/scratch/src'}),  # someone's copy
    })
    assert sorted(live_tree_users(cfg, proc)) == [101, 103, 104]


def test_sync_live_writes_what_differs_backs_it_up_and_deletes_nothing(cfg, live, tmp_path):
    os.chmod(os.path.join(live, 'photometry', 'x.py'), 0o664)
    stage = stage_tree(tmp_path, 's1', {'photometry/x.py': 'VERSION = 2\n', 'tools/run.sh': '#!/bin/bash\n'})
    rel = build(cfg, stage=stage, sha='abcdef1', log=lambda m: None)
    empty = fake_proc(tmp_path / 'p0', {})
    dry = sync_live(cfg, rel['id'], dry_run=True, proc=empty, log=lambda m: None)
    assert dry['written'] == ['VERSION', 'photometry/x.py', 'tools/run.sh'] and dry['dry_run']
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 1\n'
    res = sync_live(cfg, rel['id'], proc=empty, log=lambda m: None)
    assert res['written'] == ['VERSION', 'photometry/x.py', 'tools/run.sh']
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 2\n'
    assert stat.S_IMODE(os.stat(os.path.join(live, 'photometry', 'x.py')).st_mode) == 0o664   # keeps its mode
    assert stat.S_IMODE(os.stat(os.path.join(live, 'tools', 'run.sh')).st_mode) == 0o755     # new: writable
    assert read(os.path.join(res['backup'], 'photometry', 'x.py')) == 'VERSION = 1\n'
    assert os.path.exists(os.path.join(live, 'catcache', 'cch_00000001'))                     # nothing deleted
    assert os.path.exists(os.path.join(live, 'core.113496'))
    assert live_diff(cfg, rel['id']) == []
    assert not os.stat(os.path.join(rel['src'], 'photometry', 'x.py')).st_mode & 0o222       # release untouched
    assert sync_live(cfg, rel['id'], proc=empty, log=lambda m: None)['written'] == []


def test_sync_live_waits_for_processes_running_from_the_live_tree(cfg, live, tmp_path):
    rel = build(cfg, stage=stage_tree(tmp_path, 's1', {'photometry/x.py': 'VERSION = 2\n'}), sha='abcdef1',
                log=lambda m: None)
    zlp = ['bash', './main/ZLP_pipeline.sh', '1', '/data', '20240408', '8', '2', 'Callisto']
    busy = fake_proc(tmp_path / 'p1', {101: (zlp, {'ORCHARD_PATH': cfg['src_dir']})})
    res = sync_live(cfg, rel['id'], proc=busy, log=lambda m: None)
    assert list(res['blocked']) == [101]
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 1\n'


def test_dispatcher_finishes_the_pending_sync_once_the_live_tree_is_free(cfg, live, store, tmp_path):
    stage = stage_tree(tmp_path, 's1', {'photometry/x.py': 'VERSION = 2\n', 'jobqueue/extra.py': 'X = 1\n'})
    rel = build(cfg, stage=stage, sha='abcdef1', log=lambda m: None)
    activate(cfg, rel['id'], store=store, log=lambda m: None)
    assert store.kv_get(SYNC_KEY) == rel['id']
    zlp = ['bash', './main/ZLP_pipeline.sh', '1', '/data', '20240408', '8', '2', 'Callisto']
    proc = fake_proc(tmp_path / 'p1', {101: (zlp, {'ORCHARD_PATH': cfg['src_dir']})})
    logs = []
    d = Dispatcher(cfg, store, resources=FixedResources(), now_fn=Clock(store.now()), log=logs.append, proc=proc)
    d.poll()
    d.poll()
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 1\n'
    assert store.kv_get(SYNC_KEY) and sum('waits for 1 process' in m for m in logs) == 1   # said once
    assert not d.stop
    os.remove(os.path.join(proc, '101', 'cmdline'))                                    # the old run ends
    d.poll()
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 2\n'
    assert store.kv_get(SYNC_KEY) is None and store.events(kinds=['release-sync'])
    assert d.stop                      # jobqueue/extra.py is new: the dispatcher exits for the watchdog to restart

    activate(cfg, rel['id'], store=store, log=lambda m: None)                          # nothing left to write
    d2 = Dispatcher(cfg, store, resources=FixedResources(), now_fn=Clock(store.now()), log=logs.append, proc=proc)
    d2.poll()
    assert store.kv_get(SYNC_KEY) is None and not d2.stop


def test_cli_deploy_builds_activates_and_syncs(cfg, live, store, tmp_path, capsys, monkeypatch):
    conf = tmp_path / 'config.json'
    conf.write_text(json.dumps({k: cfg[k] for k in ('queue_root', 'basedir', 'src_dir', 'python')}))
    monkeypatch.setattr('jobqueue.releases.live_tree_users', lambda cfg, proc='/proc', exclude=(): {})
    stage = stage_tree(tmp_path, 's1', {'photometry/x.py': 'VERSION = 2\n'})
    assert cli.main(['--config', str(conf), 'release', 'deploy', '--from', stage, '--sha', 'feed123',
                     '--skip-tests']) == 0
    out = capsys.readouterr().out
    rid = current_id(load_config(str(conf)))
    assert rid and rid.endswith('-feed123') and 'wrote photometry/x.py' in out
    assert read(os.path.join(live, 'photometry', 'x.py')) == 'VERSION = 2\n'
    assert cli.main(['--config', str(conf), 'release']) == 0
    out = capsys.readouterr().out
    assert '* ' + rid in out and 'matches ' + rid in out
    assert cli.main(['--config', str(conf), 'status']) == 0
    assert 'release ' + rid in capsys.readouterr().out
    assert read_release(load_config(str(conf)), rid)['sha'] == 'feed123'
