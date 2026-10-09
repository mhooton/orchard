"""Signatures, /proc readings (against a fake /proc), ESO credential handling, the CLI and the shadow report."""
import datetime as dt
import os

import pytest

from jobqueue import cli
from jobqueue.config import DEFAULTS, load_config
from jobqueue.eso import EsoError, read_eso_credentials, scrub
from jobqueue.resources import ProcResources, external_locks, parse_sso_download_argv, parse_zlp_argv
from jobqueue.shadow_report import build_report
from jobqueue.signatures import extract, is_transient
from jobqueue.util import ensure_local_fs

# ---------------------------------------------------------------- signatures

TRACE = """START T8
[TIMEOUT] /data/SPECULOOSPipeline/Observations/Io/images/20240426/SPECU1.2024-04-26T03:00:00.000.fits: timed out
Traceback (most recent call last):
  File "/opt/orchard/src/photometry/catalogue_fov.py", line 812, in <module>
    main()
KeyError: "Column 'gaia_dr3_id' not found in /data/SPECULOOSPipeline/PipelineOutput/v2/StackImages/Sp1056+0700_I+z.fits"
"""


def test_signature_is_the_last_exception_normalised():
    sig = extract(TRACE, 1)
    assert sig == 'KeyError: "Column \'gaia_dr3_id\' not found in <path>"'
    assert not is_transient(sig, TRACE, 1, DEFAULTS)       # the per-frame [TIMEOUT] line is not the cause


def test_signature_falls_back_to_error_lines_and_exit_code():
    assert extract('all good\nERROR: wcsfit failed for 12 frames\ndone\n', 2) == 'ERROR: wcsfit failed for # frames'
    assert extract('nothing useful\n', 3) == 'exit 3'
    assert extract('', -9) == 'killed by signal 9'


def test_transient_patterns_and_exit_codes():
    t = 'writing...\nOSError: [Errno 28] No space left on device\n'
    assert is_transient(extract(t, 1), t, 1, DEFAULTS)
    assert is_transient('exit 75', '', 75, DEFAULTS)
    assert is_transient('killed by signal 9', '', 137, DEFAULTS)
    n = "requests.exceptions.ConnectionError: HTTPSConnectionPool(host='archive.eso.org'): Max retries exceeded"
    assert is_transient(extract(n, 1), n, 1, DEFAULTS)


# ---------------------------------------------------------------- /proc


def fake_proc(tmp_path, procs=(), loadavg='12.0 34.5 20.0 3/2890 38119', diskstats=None):
    proc = tmp_path / 'proc'
    proc.mkdir()
    (proc / 'loadavg').write_text(loadavg + '\n')
    (proc / 'meminfo').write_text('MemTotal:       527921052 kB\nMemAvailable:   472838808 kB\n')
    (proc / 'diskstats').write_text(diskstats or (
        '   8      16 sdb 130716918 1360196 11684728994 805362464 127999739 1795935 46286150941 144516072 0 261368273 949723701\n'
        '   8       0 sda 157733935 4335553 6689317591 795294270 3507632665 3041729 532999216703 976278357 0 956591916 1758555317\n'))
    for pid, sid, argv in procs:
        d = proc / str(pid)
        d.mkdir()
        (d / 'cmdline').write_bytes(b'\0'.join(a.encode() for a in argv) + b'\0')
        fields = ['S', '1', str(sid), str(sid)] + ['0'] * 15 + ['123456']
        (d / 'stat').write_text('{} (bash) {}\n'.format(pid, ' '.join(fields)))
    return str(proc)


def test_proc_readings_and_disk_util(tmp_path):
    proc = fake_proc(tmp_path)
    t = [100.0]
    r = ProcResources(['sda', 'sdb'], proc=proc, clock=lambda: t[0], sleep=lambda s: t.__setitem__(0, t[0] + s))
    assert r.load5() == 34.5
    assert round(r.mem_available_gb()) == 451
    u = r.disk_util()                        # first call samples twice; no change in between -> 0 %
    assert u == {'sda': 0.0, 'sdb': 0.0}
    with open(os.path.join(proc, 'diskstats'), 'w') as f:
        f.write('   8      16 sdb 0 0 0 0 0 0 0 0 0 261398273 0\n   8       0 sda 0 0 0 0 0 0 0 0 0 956591916 0\n')
    t[0] += 50.0                             # 30 s of io_ticks over 50 s on sdb
    assert r.disk_util() == {'sda': 0.0, 'sdb': 60.0}


def test_external_pipeline_and_download_processes_hold_locks(tmp_path):
    zlp = ['bash', './main/ZLP_pipeline.sh', '--force-platesolve', '1', '/data/SPECULOOSPipeline', '20261008', '8', '2',
           'Europa']
    sso = ['python', 'download/SSO_download.py', '--dir', '/data/SPECULOOSPipeline', '--telescope', 'Io', '--sdate',
           '20261007', '--edate', '20261008', '--max-retries', '25']
    ours = ['bash', './main/ZLP_pipeline.sh', '1', '/data/SPECULOOSPipeline', '20261008', '8', '2', 'Callisto']
    proc = fake_proc(tmp_path, [(10, 10, zlp), (11, 11, sso), (12, 500, ours), (13, 13, ['sleep', '5'])])
    held = external_locks(own_sessions=[500], proc=proc)
    assert set(held) == {'tel:Europa', 'night:Europa:20261008', 'night:Io:20261007'}
    assert parse_zlp_argv(['./main/ZLP_pipeline.sh', '--cores', '8', '1', 'B', '20261001 20261002', '8', '2', 'Io',
                           'Sp0001 Sp0002']) == {'tel:Io', 'night:Io:20261001', 'night:Io:20261002',
                                                 'target:SP0001', 'target:SP0002'}
    assert parse_sso_download_argv(['SSO_download.py', '--telescope', 'all', '--sdate', '20261007']) == {
        'night:Io:20261007', 'night:Europa:20261007', 'night:Ganymede:20261007', 'night:Callisto:20261007'}


def test_queue_database_refuses_nfs(tmp_path):
    mounts = tmp_path / 'mounts'
    mounts.write_text('appct:/export/data /appct/data nfs4 rw 0 0\n/dev/sdb1 /data/SPECULOOSPipeline xfs rw 0 0\n')
    with pytest.raises(RuntimeError, match='NFS|nfs'):
        ensure_local_fs('/appct/data/SPECULOOSPipeline/queue', str(mounts))
    assert ensure_local_fs('/data/SPECULOOSPipeline/queue', str(mounts)) == 'xfs'


# ---------------------------------------------------------------- ESO credentials


def test_credentials_read_without_touching_the_environment(tmp_path, monkeypatch):
    monkeypatch.delenv('ESO_USERNAME', raising=False)
    monkeypatch.delenv('EMAIL_PASSWORD', raising=False)
    env = tmp_path / '.env'
    env.write_text('EMAIL_PASSWORD=mailsecret\nESO_USERNAME=speculoos\nESO_PASSWORD="p@ss w0rd&x"\n')
    assert read_eso_credentials(str(env)) == ('speculoos', 'p@ss w0rd&x')
    assert 'EMAIL_PASSWORD' not in os.environ and 'ESO_USERNAME' not in os.environ
    with pytest.raises(EsoError):
        read_eso_credentials(str(tmp_path / 'missing.env'))


def test_scrub_removes_raw_and_url_encoded_secrets():
    msg = ("Error getting OAuth2.0 token: HTTPSConnectionPool(host='www.eso.org', port=443): Max retries exceeded "
           "with url: /sso/oidc/token?username=speculoos&password=p%40ss+w0rd%26x (Caused by ...)")
    out = scrub(msg, ['p@ss w0rd&x', 'speculoos'])
    assert 'w0rd' not in out and 'speculoos' not in out and 'password=***' in out


# ---------------------------------------------------------------- CLI


def test_cli_end_to_end(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'q'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    assert cli.main(['add', '--class', 'P2', '--telescope', 'Europa', '--night', '20250101',
                     '--lock-targets', 'Sp0001+0001', '--estimate', '30']) == 0
    assert 'added job 1' in capsys.readouterr().out
    assert cli.main(['add', '--class', 'P3', '--kind', 'sweep', '--cores', '1', '--no-disk-heavy', '--',
                     'echo', 'hi']) == 0
    assert cli.main(['pause', 'P2']) == 0
    assert cli.main(['drain']) == 0
    assert cli.main(['status']) == 0
    out = capsys.readouterr().out
    assert 'drain      ON' in out and 'paused: P2' in out and 'dispatcher NOT RUNNING' in out
    assert cli.main(['undrain']) == 0 and cli.main(['resume', 'all']) == 0
    assert cli.main(['cancel', '2']) == 0
    assert 'job 2: cancelled' in capsys.readouterr().out
    assert cli.main(['requeue', '--signature', 'nothing-like-this']) == 0
    assert cli.main(['show', '1']) == 0
    out = capsys.readouterr().out
    assert "'tel:Europa'" in out and 'target:SP0001+0001' in out and '--force-platesolve' in out


def test_config_file_and_env(tmp_path, monkeypatch):
    p = tmp_path / 'c.json'
    p.write_text('{"mode": "shadow", "watcher": {"stable_minutes": 40}}')
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'root'))
    cfg = load_config(str(p))
    assert cfg['mode'] == 'shadow' and cfg['watcher']['stable_minutes'] == 40
    assert cfg['watcher']['watch_days'] == 2 and cfg['queue_root'] == str(tmp_path / 'root')
    p.write_text('{"mode": "sideways"}')
    with pytest.raises(ValueError):
        load_config(str(p))


# ---------------------------------------------------------------- shadow report


def test_shadow_report_compares_with_the_cron(cfg, store, clock):
    base = cfg['basedir']
    os.makedirs(os.path.join(base, 'ESO_logs'))
    with open(os.path.join(base, 'ESO_logs', '20261008Europa.log'), 'w') as f:
        f.write("['/opt/orchard/src']\n2026-10-09 19:00:02.123456\nINITIAL ATTEMPT\nRETRY ATTEMPT 2/25\n")
    logdir = os.path.join(base, 'PipelineOutput', 'v2', 'Europa', 'logs')
    os.makedirs(logdir)
    with open(os.path.join(logdir, '20261008_1_v3.log'), 'w') as f:
        f.write('PIPELINE VERSION v3\nFri Oct  9 19:10:00 UTC 2026\n...\nFri Oct  9 19:52:00 UTC 2026\nPIPELINE COMPLETE\n')
    os.makedirs(os.path.join(base, 'Observations', 'Europa'))
    with open(os.path.join(base, 'Observations', 'Europa', 'download_log.csv'), 'w') as f:
        f.write('Night,Telescope,ESO_Archive,Transferred,Downloaded\n20261008,Europa,874,874,874\n')

    dl, _ = store.add_job('P0', 'download', ['x'], telescope='Europa', night='20261008', est_minutes=9)
    pl, _ = store.add_job('P0', 'pipeline', ['x'], telescope='Europa', night='20261008', est_minutes=60,
                          depends_on=dl)
    clock.t = dt.datetime(2026, 10, 9, 15, 0)
    store.mark_running(dl, shadow=True)
    clock.advance(minutes=10)
    store.finish(dl, 0, shadow=True)
    store.mark_running(pl, shadow=True)
    clock.advance(minutes=60)
    store.finish(pl, 0, shadow=True)
    store.watch_put({'telescope': 'Europa', 'night': '20261008', 'state': 'enqueued', 'ready_at': '2026-10-09T15:00:00Z',
                     'ready_reason': 'ESO has 874 frames, transfer log 874', 'download_job': dl, 'pipeline_job': pl,
                     'eso_rows': 874, 'transfer_count': 874, 'noted': 0, 'polls_unchanged': 0})
    store.watch_put({'telescope': 'Io', 'night': '20261008', 'state': 'no_data', 'noted': 1, 'polls_unchanged': 3,
                     'note': 'Io 20261008: no frames at ESO and no transfer-log entry; nothing to download or process'})
    text = build_report(dict(cfg, mode='shadow'), store, days=1, now=dt.datetime(2026, 10, 10, 7, 0))
    assert '| Europa | ready 09 15:00' in text
    assert '+4.0 h' in text                                    # would start 15:10, cron started 19:10
    assert 'cron download made 2 attempts' in text and 'cron counts ESO 874 / transfer 874 / downloaded 874' in text
    assert 'cron pipeline completed (log in v2)' in text
    assert 'Io 20261008: no frames at ESO' in text
    assert 'median +4.0 h earlier' in text


def test_live_report_shows_jobs_results_and_delay(cfg, store, clock):
    base = cfg['basedir']
    lc = os.path.join(base, 'PipelineOutput', 'v2', 'Ganymede', 'output', '20261008', 'Sp0055-3052')
    os.makedirs(lc)
    open(os.path.join(lc, 'Sp0055-3052_I+z_5_diff.fits'), 'w').close()
    clock.t = dt.datetime(2026, 10, 9, 15, 0)
    dl, _ = store.add_job('P0', 'download', ['x'], telescope='Ganymede', night='20261008', est_minutes=10)
    pl, _ = store.add_job('P0', 'pipeline', ['x'], telescope='Ganymede', night='20261008', est_minutes=30,
                          depends_on=dl)
    store.mark_running(dl)
    clock.advance(minutes=8)
    store.finish(dl, 0)
    store.mark_running(pl)
    clock.advance(minutes=22)
    store.finish(pl, 0)
    store.watch_put({'telescope': 'Ganymede', 'night': '20261008', 'state': 'enqueued', 'noted': 0,
                     'polls_unchanged': 1, 'first_seen': '2026-10-09T14:00:00Z', 'ready_at': '2026-10-09T15:00:00Z',
                     'download_job': dl, 'pipeline_job': pl})
    rerun, _ = store.add_job('P1', 'pipeline', ['x'], telescope='Callisto', night='20260918', note='no light curves')
    store.mark_running(rerun)
    store.finish(rerun, 1, signature="KeyError: 'gaia_dr3_id'")
    text = build_report(dict(cfg, mode='live'), store, days=1, now=dt.datetime(2026, 10, 10, 6, 30))
    assert text.startswith('# Job queue daily report')
    assert '| Ganymede | 09 15:00 | 09 15:00–09 15:08 | 09 15:08–09 15:30 | exit 0 | v2 | 0.5 h |' in text
    assert 'first seen at ESO 09 14:00' in text and 'jobs 1 + 2' in text
    assert 'Ready to done: median 0.5 h over 1 telescope-nights' in text and '1 failed job(s)' in text
    assert "job 3 P1 pipeline Callisto 20260918: [deterministic] KeyError: 'gaia_dr3_id'" in text
    assert '- job 3 pipeline Callisto 20260918: failed (no light curves)' in text
