"""The download and fetch job bodies, with the ESO side faked."""
import datetime as dt
import os

import pytest

from jobqueue import jobs
from jobqueue.tests.fakes import write_fits
from jobqueue.util import utcnow, write_json_atomic


class Archive:
    def __init__(self, rows):
        self.rows = rows

    def frames(self, prog, first, last):
        return [dict(r) for r in self.rows]


def recent(days=1):
    return (utcnow() - dt.timedelta(days=days)).strftime('%Y%m%d')


@pytest.fixture
def no_sso_download(monkeypatch):
    calls = []
    monkeypatch.setattr(jobs.subprocess, 'call', lambda cmd, **kw: calls.append(cmd) or 0)
    return calls


def test_download_refuses_old_nights(cfg, monkeypatch, no_sso_download):
    monkeypatch.setattr(jobs, '_archive', lambda cfg: pytest.fail('must not even talk to ESO'))
    assert jobs.download(cfg, 'Callisto', '20190415') == 2
    assert no_sso_download == []


def test_download_runs_sso_download_for_a_recent_night(cfg, monkeypatch, no_sso_download):
    night = recent()
    rows = [{'dp_id': 'SPECU2.{}T01:00:00.000'.format(dt.datetime.strptime(night, '%Y%m%d').strftime('%Y-%m-%d')),
             'dp_type': 'OBJECT', 'origfile': ''}]
    monkeypatch.setattr(jobs, '_archive', lambda cfg: Archive(rows))
    os.makedirs(os.path.join(cfg['basedir'], 'ESO_logs'))
    d = os.path.join(cfg['basedir'], 'Observations', 'Europa', 'images', night)
    write_fits(os.path.join(d, rows[0]['dp_id'] + '.fits'), OBSERVER='Astra', ASTRAROT=90)
    assert jobs.download(cfg, 'Europa', night) == 0
    (cmd,) = no_sso_download
    assert cmd[1:] == ['download/SSO_download.py', '--dir', cfg['basedir'], '--telescope', 'Europa', '--sdate', night,
                       '--edate', (dt.datetime.strptime(night, '%Y%m%d') + dt.timedelta(days=1)).strftime('%Y%m%d'),
                       '--max-retries', '1']
    assert os.path.exists(os.path.join(cfg['basedir'], 'ESO_logs', '{}Europa.log'.format(night)))


def test_download_retries_later_when_nothing_landed(cfg, monkeypatch, no_sso_download):
    monkeypatch.setattr(jobs, '_archive', lambda cfg: Archive([{'dp_id': 'SPECU2.X', 'dp_type': 'OBJECT',
                                                                 'origfile': ''}]))
    os.makedirs(os.path.join(cfg['basedir'], 'ESO_logs'))
    assert jobs.download(cfg, 'Europa', recent()) == jobs.EX_TEMPFAIL


def test_download_never_lets_transformation_check_delete_frames(cfg, monkeypatch, no_sso_download):
    night = recent()
    d = os.path.join(cfg['basedir'], 'Observations', 'Callisto', 'images', night)
    write_fits(os.path.join(d, 'SPECU4.2026-01-01T01:00:00.000.fits'), OBSERVER='speculoos')   # no ASTRAMIR
    fetched = []
    monkeypatch.setattr(jobs, '_archive', lambda cfg: Archive([]))
    monkeypatch.setattr(jobs, 'fetch_rows', lambda *a, **k: fetched.append(a) or {})
    assert jobs.download(cfg, 'Callisto', night) == 0
    assert no_sso_download == [] and len(fetched) == 1          # add-only fetch instead
    assert os.path.exists(os.path.join(d, 'SPECU4.2026-01-01T01:00:00.000.fits'))


def test_fetch_queues_a_rerun_when_science_was_added(cfg, store, monkeypatch):
    rows_file = os.path.join(cfg['queue_root'], 'rows.json')
    write_json_atomic(rows_file, [{'dp_id': 'SPECU4.A', 'dp_type': 'OBJECT', 'object': 'SP1945-2557'}])
    result = {'requested': 1, 'downloaded': 1, 'added': 120, 'added_science': 120, 'settled': ['SPECU4.Z']}
    monkeypatch.setattr(jobs, '_archive', lambda cfg: None)
    monkeypatch.setattr(jobs, 'fetch_rows', lambda *a, **k: dict(result))
    assert jobs.fetch(cfg, store, 'Callisto', '20260911', rows_file) == 0
    (rerun,) = store.jobs(kind='pipeline')
    assert rerun['class'] == 'P1' and rerun['night'] == '20260911' and 'target:SP1945-2557' in rerun['locks']
    lb = store.lookback_get('Callisto', '20260911')
    assert lb['reruns'] == 1 and lb['data']['settled'] == ['SPECU4.Z']
    assert os.path.exists(rows_file + '.result.json')


def test_fetch_trial_into_scratch_records_nothing(cfg, store, monkeypatch, tmp_path):
    rows_file = os.path.join(cfg['queue_root'], 'rows.json')
    write_json_atomic(rows_file, [{'dp_id': 'SPECU4.A', 'dp_type': 'OBJECT', 'object': 'X'}])
    seen = {}
    monkeypatch.setattr(jobs, '_archive', lambda cfg: None)
    monkeypatch.setattr(jobs, 'fetch_rows', lambda *a, **k: seen.update(k) or {'added_science': 5})
    assert jobs.fetch(cfg, store, 'Callisto', '20260911', rows_file, dest=str(tmp_path / 'trial')) == 0
    assert seen['dest'] == str(tmp_path / 'trial')
    assert store.jobs() == [] and store.lookback_rows() == []
