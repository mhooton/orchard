"""The notifier against a fake production tree: digest lines, incident rules, threading, acks, delivery."""
import datetime as dt
import email
import email.policy
import glob
import json
import os

import pytest

from jobqueue import digest as D
from jobqueue import esoreach, pwvlog
from jobqueue.notify import Notifier, read_smtp_settings, NotifyError
from jobqueue.util import iso, parse_night

NOW = dt.datetime(2026, 10, 18, 7, 45)        # last_night = 20261016
NIGHT = '20261016'
BIG = (1, 50 * 2 ** 40, 100 * 2 ** 40)        # (st_dev, free bytes, total bytes)


def usage_ok(path):
    return BIG


def night_before(night, k=1):
    return (parse_night(night) - dt.timedelta(days=k)).strftime('%Y%m%d')


def frames(cfg, tel, night, n, start=None):
    d = os.path.join(cfg['basedir'], 'Observations', tel, 'images', night)
    os.makedirs(d, exist_ok=True)
    t0 = start or dt.datetime.combine(parse_night(night), dt.time(23, 0))
    for k in range(n):
        stamp = (t0 + dt.timedelta(seconds=60 * k)).strftime('%Y-%m-%dT%H:%M:%S.000')
        open(os.path.join(d, 'SPECU1.{}.fits'.format(stamp)), 'w').close()


def lightcurve(cfg, tel, night, target, version='v2'):
    d = os.path.join(cfg['basedir'], 'PipelineOutput', version, tel, 'output', night, target)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, '{}_I+z_{}_5_diff.fits'.format(target, night)), 'w').close()


def blocks(cfg, tel, night, rows):
    d = os.path.join(cfg['basedir'], 'Observations', tel, 'schedule', 'Archive_night_blocks')
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'night_blocks_{}_{}.txt'.format(tel, parse_night(night).isoformat())), 'w') as f:
        f.write('target "start time (UTC)" "end time (UTC)" "duration (minutes)" "ra (h)" configuration\n')
        for target, minutes in rows:
            f.write('{} "2026-10-16 23:30:00.000" "2026-10-17 09:47:00.000" {} 14.0 "{{\'texp\': 60}}"\n'.format(
                target, minutes))


def watch(store, tel, night, state, rows=0, transfer=None, targets=(), last_poll=None):
    store.watch_put({'telescope': tel, 'night': night, 'state': state, 'eso_rows': rows, 'transfer_count': transfer,
                     'polls_unchanged': 0, 'noted': 0, 'last_poll': last_poll,
                     'data': {'targets': list(targets), 'by_type': {'OBJECT': rows} if rows else {}}})


def job(store, tel, night, kind='pipeline', cls='P0', state='done', signature=None, note=None, campaign=None,
        est=60.0, not_before=None):
    jid, _ = store.add_job(cls, kind, ['true'], telescope=tel, night=night, note=note, est_minutes=est,
                           meta={'campaign': campaign} if campaign else None, not_before=not_before)
    if state in ('running', 'done', 'failed'):
        store.mark_running(jid)
    if state == 'done':
        store.finish(jid, 0)
    elif state == 'failed':
        store.finish(jid, 1, failure_kind='deterministic', signature=signature)
    return jid


def heartbeat(store, at=NOW):
    store.kv_set('dispatcher', json.dumps({'heartbeat': iso(at)}))


def observed_night(cfg, store, tel, night, n=400):
    watch(store, tel, night, 'enqueued', rows=n)


def notifier(cfg, store, now=NOW, **kw):
    kw.setdefault('usage', usage_ok)
    kw.setdefault('smtp_factory', NoSmtp)
    kw.setdefault('log', lambda m: None)
    return Notifier(cfg, store, now=now, **kw)


def outbox(cfg):
    return sorted(glob.glob(os.path.join(cfg['queue_root'], 'notify', 'outbox', '*.eml')))


def read(path):
    with open(path, 'rb') as f:
        return email.message_from_binary_file(f, policy=email.policy.default)


class NoSmtp:
    def __init__(self, *a, **k):
        raise AssertionError('SMTP must not be used')


def artemis_observing(store):
    """Artemis has no sister telescope, so any run of dark nights would be an incident: keep it observing."""
    for k in range(10):
        watch(store, 'Artemis', night_before('20261018', k), 'enqueued', rows=300)


@pytest.fixture
def clean(cfg, store, clock):
    """A healthy queue at NOW: a fresh dispatcher heartbeat, no jobs, Artemis observing."""
    clock.t = NOW
    heartbeat(store)
    artemis_observing(store)
    return cfg


# ---------------------------------------------------------------------------------------------- the digest

def test_digest_lines_for_each_kind_of_night(clean, store):
    cfg = clean
    # Europa: done, light curves for both targets
    frames(cfg, 'Europa', NIGHT, 300)
    watch(store, 'Europa', NIGHT, 'enqueued', rows=300, targets=['SP0246+1625', 'SP0020+3305'])
    job(store, 'Europa', NIGHT, 'download')
    job(store, 'Europa', NIGHT)
    lightcurve(cfg, 'Europa', NIGHT, '111')
    lightcurve(cfg, 'Europa', NIGHT, '222', version='v3')
    # Ganymede: the pipeline failed
    frames(cfg, 'Ganymede', NIGHT, 200)
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=200)
    job(store, 'Ganymede', NIGHT, 'download')
    job(store, 'Ganymede', NIGHT, state='failed', signature="KeyError: 'gaia_dr3_id'")
    # Callisto: a Chilean night, nothing for us
    watch(store, 'Callisto', NIGHT, 'no_data')
    blocks(cfg, 'Callisto', NIGHT, [('ch_em01_223138_020_ML00002', 617.0)])
    # Io: dark while two others observed (one night only: no incident)
    watch(store, 'Io', NIGHT, 'no_data')
    # Artemis: still arriving
    watch(store, 'Artemis', NIGHT, 'arriving', rows=12)

    n = notifier(cfg, store, smtp_factory=NoSmtp)
    n.digest()
    files = outbox(cfg)
    digests = [read(p) for p in files if '_digest-' in p]
    assert len(digests) == 1
    msg = digests[0]
    body = msg.get_content()
    assert msg['Subject'] == '[orchard] Night of 16 Oct 2026: 1 needs attention'
    assert 'Europa    300 frames, 5.0 h · light curves for 2 of 2 targets · done' in body
    assert "✗ Ganymede  200 frames, 3.3 h · pipeline failed (job" in body and "KeyError: 'gaia_dr3_id'" in body
    assert 'Callisto  Chilean night (10.3 h planned for Chile; not ours to download)' in body
    assert '  Io        no data; Europa and Ganymede observed' in body
    assert 'Artemis   waiting for the night to arrive' in body
    assert 'Open incidents (1)' in body and 'Ganymede, night of 16 Oct 2026: pipeline job failed' in body
    # the failure also got its own alert
    alerts = [read(p) for p in files if '_alert-' in p]
    assert [a['Subject'] for a in alerts] == ['[orchard] Ganymede, night of 16 Oct 2026: pipeline job failed']


def test_all_well_digest_still_goes_out(clean, store):
    cfg = clean
    for tel in cfg['telescopes']:
        frames(cfg, tel, NIGHT, 200)
        watch(store, tel, NIGHT, 'enqueued', rows=200)
        job(store, tel, NIGHT)
        lightcurve(cfg, tel, NIGHT, '42')
    n = notifier(cfg, store)
    out = n.digest()
    assert 'written via outbox' in out
    [path] = outbox(cfg)
    assert read(path)['Subject'] == '[orchard] Night of 16 Oct 2026: all well'
    assert store.kv_get('notify_last_digest') == iso(NOW)


def test_print_changes_and_delivers_nothing(clean, store):
    cfg = clean
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=10)
    job(store, 'Ganymede', NIGHT, state='failed', signature='boom')
    text = notifier(cfg, store).digest(dry_run=True)
    assert text.startswith('[orchard] Night of 16 Oct 2026:') and 'pipeline job failed' in text
    assert outbox(cfg) == []
    assert notifier(cfg, store).incidents() == []


# ---------------------------------------------------------------------------------------------- incidents

def dark_io(cfg, store, nights, chilean=()):
    """Io dark on each of `nights`; Europa, Ganymede and Callisto observed; Io observed the night before."""
    for night in nights:
        watch(store, 'Io', night, 'no_data')
        if night in chilean:
            blocks(cfg, 'Io', night, [('ch_Em01_121132_020_ML00007', 600.0)])
        for tel in ('Europa', 'Ganymede', 'Callisto'):
            observed_night(cfg, store, tel, night)
    before = night_before(nights[-1])
    observed_night(cfg, store, 'Io', before)


def last_four(night=NIGHT):
    return [night_before(night, k) for k in range(4)]


def test_dark_telescope_opens_one_incident_and_alerts_once(clean, store):
    cfg = clean
    dark_io(cfg, store, last_four())
    n = notifier(cfg, store)
    n.check()
    [inc] = n.incidents()
    assert inc['key'] == 'dark:Io' and inc['title'] == 'Io: no data since 13 Oct 2026'
    assert 'On 4 of them at least two of the other SSO telescopes observed' in inc['detail']
    assert len(outbox(cfg)) == 1
    notifier(cfg, store, now=NOW + dt.timedelta(minutes=30)).check()
    assert len(outbox(cfg)) == 1                    # still open: no second alert


def test_dark_incident_resolves_in_the_same_thread(clean, store):
    cfg = clean
    dark_io(cfg, store, last_four())
    notifier(cfg, store).check()
    [opening] = [read(p) for p in outbox(cfg)]
    # the next night Io observes again
    later = NOW + dt.timedelta(days=1)
    observed_night(cfg, store, 'Io', '20261017')
    heartbeat(store, later)
    n = notifier(cfg, store, now=later)
    n.check()
    resolved = [read(p) for p in outbox(cfg) if '_resolved-' in p]
    assert len(resolved) == 1
    assert resolved[0]['Subject'] == 'Re: ' + opening['Subject']
    assert resolved[0]['In-Reply-To'] == opening['Message-ID']
    assert n.incidents() == []


def test_chilean_nights_do_not_count_towards_dark(clean, store):
    cfg = clean
    nights = last_four()
    dark_io(cfg, store, nights, chilean=nights[:2])
    n = notifier(cfg, store)
    assert n.check(dry_run=True) == ['no conditions']


def test_site_wide_dark_is_the_weather_not_an_incident(clean, store):
    cfg = clean
    for night in last_four():
        for tel in ('Io', 'Europa', 'Ganymede', 'Callisto'):
            watch(store, tel, night, 'no_data')
    assert notifier(cfg, store).check(dry_run=True) == ['no conditions']


def test_failed_nightly_job_resolves_when_a_rerun_succeeds(clean, store):
    cfg = clean
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=10)
    job(store, 'Ganymede', NIGHT, state='failed', signature='boom')
    n = notifier(cfg, store)
    n.check()
    assert [i['key'] for i in n.incidents()] == ['failed:Ganymede:20261016:pipeline']
    job(store, 'Ganymede', NIGHT, cls='P1')
    later = NOW + dt.timedelta(hours=2)
    heartbeat(store, later)
    n = notifier(cfg, store, now=later)
    n.check()
    assert n.incidents() == []
    assert len([p for p in outbox(cfg) if '_resolved-' in p]) == 1


def test_backlog_failures_alert_only_in_bulk(clean, store):
    cfg = clean
    for k in range(9):
        job(store, 'Io', '2024010{}'.format(k % 9), cls='P2', state='failed', signature='ValueError: x')
    assert notifier(cfg, store).check(dry_run=True) == ['no conditions']
    job(store, 'Io', '20240110', cls='P2', state='failed', signature='ValueError: x')
    [line] = notifier(cfg, store).check(dry_run=True)
    assert line.startswith('signature:') and 'Repeated job failure: ValueError: x' in line


def test_frames_missing_at_eso_after_the_deadline(clean, store):
    cfg = clean
    night = '20261017'
    watch(store, 'Europa', night, 'arriving', rows=300, transfer=812, last_poll='2026-10-18T20:30:00Z')
    early = dt.datetime(2026, 10, 18, 20, 0)
    heartbeat(store, early)
    assert notifier(cfg, store, now=early).check(dry_run=True) == ['no conditions']
    late = dt.datetime(2026, 10, 18, 21, 30)
    heartbeat(store, late)
    [line] = notifier(cfg, store, now=late).check(dry_run=True)
    assert line == 'late:Europa:20261017  Europa, night of 17 Oct 2026: frames missing at ESO'


def test_stalled_queue_and_dead_dispatcher(cfg, store, clock):
    clock.t = NOW
    artemis_observing(store)
    jid, _ = store.add_job('P2', 'pipeline', ['true'], telescope='Io', night='20240101')
    store.kv_set('dispatcher', json.dumps({'heartbeat': iso(NOW - dt.timedelta(hours=7))}))
    keys = [l.split()[0] for l in notifier(cfg, store).check(dry_run=True)]
    assert keys == ['dispatcher', 'queue-stalled']
    store.pause('P2')
    heartbeat(store)
    assert notifier(cfg, store).check(dry_run=True) == ['no conditions']


def test_disk_nearly_full(clean, store):
    cfg = clean
    n = notifier(cfg, store, usage=lambda p: (1, int(1.5 * 2 ** 40), 100 * 2 ** 40))
    [line] = n.check(dry_run=True)
    assert line == 'disk:raw data  raw data disk nearly full'


def test_ack_holds_the_alert_but_the_digest_lists_it(clean, store):
    cfg = clean
    dark_io(cfg, store, last_four())
    n = notifier(cfg, store)
    n.sync(n.evaluate(statuses=False)[-1])          # open it without alerting
    [inc] = n.incidents()
    n.ack(inc['id'], '2026-11-01', 'Seba knows')
    n.digest()
    assert [p for p in outbox(cfg) if '_alert-' in p] == []
    [digest] = [read(p).get_content() for p in outbox(cfg) if '_digest-' in p]
    assert 'acknowledged until 2026-11-01: Seba knows' in digest
    n.unack(inc['id'])
    notifier(cfg, store, now=NOW + dt.timedelta(minutes=30)).check()
    assert len([p for p in outbox(cfg) if '_alert-' in p]) == 1


# ---------------------------------------------------------------------------------------------- delivery

class FakeSmtp:
    sent = []

    def __init__(self, host, port, timeout=None):
        self.host, self.port, self.calls = host, port, []

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def starttls(self):
        self.calls.append('starttls')

    def login(self, user, password):
        self.calls.append(('login', user))

    def send_message(self, msg, from_addr=None, to_addrs=None):
        FakeSmtp.sent.append((msg, from_addr, list(to_addrs), self.host, self.port, list(self.calls)))


def smtp_cfg(cfg, tmp_path, env=True):
    path = tmp_path / 'notify.env'
    if env:
        path.write_text('EMAIL_ADDRESS=speculoos.pipeline@gmail.com\nEMAIL_PASSWORD="app pass"\n'
                        'SMTP_SERVER=smtp.example.org:2525\n')
    return dict(cfg, notify=dict(cfg['notify'], transport='smtp', env_file=str(path), to=['someone@example.org']))


def test_smtp_sends_utf8_and_keeps_an_outbox_copy(clean, store, tmp_path):
    FakeSmtp.sent = []
    cfg = smtp_cfg(clean, tmp_path)
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=10)
    job(store, 'Ganymede', NIGHT, state='failed', signature='boom')
    n = notifier(cfg, store, smtp_factory=FakeSmtp)
    out = n.digest()
    assert 'sent via smtp' in out
    subjects = [m['Subject'] for m, *_ in FakeSmtp.sent]
    assert subjects == ['[orchard] Ganymede, night of 16 Oct 2026: pipeline job failed',
                        '[orchard] Night of 16 Oct 2026: 1 needs attention']
    msg, frm, to, host, port, calls = FakeSmtp.sent[1]
    assert (frm, to, host, port) == ('speculoos.pipeline@gmail.com', ['someone@example.org'], 'smtp.example.org', 2525)
    assert calls == ['starttls', ('login', 'speculoos.pipeline@gmail.com')]
    assert msg.as_bytes().decode('ascii')                       # quoted-printable: the wire stays ASCII
    assert '✗ Ganymede' in msg.get_content()
    assert len(outbox(cfg)) == 2
    assert {m['status'] for m in n.messages()} == {'sent'}


def test_smtp_without_settings_still_writes_the_outbox(clean, store, tmp_path):
    cfg = smtp_cfg(clean, tmp_path, env=False)
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=10)
    job(store, 'Ganymede', NIGHT, state='failed', signature='boom')
    n = notifier(cfg, store, smtp_factory=NoSmtp)
    n.digest()
    assert len(outbox(cfg)) == 2
    assert {m['status'] for m in n.messages()} == {'failed'}
    [inc] = n.incidents()
    assert inc['alerted_at'] is None                             # tried again on the next run


def test_outbox_flag_overrides_smtp(clean, store, tmp_path):
    cfg = smtp_cfg(clean, tmp_path)
    notifier(cfg, store, outbox_only=True, smtp_factory=NoSmtp).digest()
    assert len(outbox(cfg)) == 1


def test_read_smtp_settings(tmp_path):
    p = tmp_path / '.env'
    p.write_text('# comment\nexport EMAIL_ADDRESS=a@b.c\nEMAIL_PASSWORD=\'x y\'\nESO_PASSWORD=zzz\n')
    assert read_smtp_settings(str(p)) == {'user': 'a@b.c', 'password': 'x y', 'host': 'smtp.gmail.com', 'port': 587}
    assert 'ESO_PASSWORD' not in os.environ
    with pytest.raises(NotifyError):
        read_smtp_settings(str(tmp_path / 'missing'))


def test_heartbeat_only_when_configured(clean, store):
    calls = []

    class Resp:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, n):
            return b'OK'

    def opener(url, timeout=None):
        calls.append(url)
        return Resp()

    notifier(clean, store, opener=opener).digest()
    assert calls == []
    cfg = dict(clean, notify=dict(clean['notify'], heartbeat_url='https://hc.example/ping/abc'))
    out = notifier(cfg, store, opener=opener, now=NOW + dt.timedelta(days=1)).digest()
    assert calls == ['https://hc.example/ping/abc'] and 'heartbeat ok' in out


# ---------------------------------------------------------------------------------------------- helpers

def test_planned_minutes_split_chilean_from_survey(cfg):
    blocks(cfg, 'Io', NIGHT, [('Sp1045-0149', 200.0), ('ch_em01_225147_020_ML00015', 400.0), ('"ch_x"', 10.0)])
    assert D.planned_minutes(cfg['basedir'], 'Io', NIGHT) == {'survey': 200.0, 'chilean': 410.0}
    assert D.planned_minutes(cfg['basedir'], 'Io', '20261015') is None


def test_calibrations_alone_are_not_observing(clean, store):
    cfg = clean
    frames(cfg, 'Io', '20260901', 80)                # no watch row: an old night on disk
    frames(cfg, 'Io', '20260902', 400)
    ctx = D.Context(cfg, store, NOW)
    assert ctx.observed('Io', '20260901') is False
    assert ctx.observed('Io', '20260902') is True
    lightcurve(cfg, 'Io', '20260901', '7')
    assert D.Context(cfg, store, NOW).observed('Io', '20260901') is True


# ---------------------------------------------------------------------------------------------- the queue's work

def test_campaign_labels():
    def label(cls, kind='pipeline', campaign=None, note=None):
        return D.campaign_label({'class': cls, 'kind': kind, 'meta': {'campaign': campaign} if campaign else {},
                                 'note': note})
    assert label('P3', campaign='WCS audit reruns', note='anything') == 'WCS audit reruns'
    assert label('P0') == 'nightly' and label('P1', 'fetch', note='late frames') == 'look-back'
    assert label('P3', 'framepack', note='16-bit archive sweep: v2 510 frames in 2 folders') == '16-bit archive sweep'
    assert label('P3', note='never reduced; raw frames lack WCS (WCS audit 2026-10-16)') == \
        'never reduced; raw frames lack WCS'
    assert label('P2') == 'reruns' and label('P2', 'framepack') == 'frame packing'
    assert len(label('P3', note='x' * 80)) == 40


def test_digest_summarises_processed_and_remaining_work(clean, store):
    cfg = clean
    for night in ('20240118', '20240119', '20240120'):
        job(store, 'Callisto', night, cls='P2', campaign='WCS audit reruns')
    lightcurve(cfg, 'Callisto', '20240118', '9')
    lightcurve(cfg, 'Callisto', '20240119', '9')
    job(store, 'Ganymede', '20240121', cls='P2', state='failed', signature='boom')
    job(store, 'Callisto', '20240122', cls='P2', state='running', campaign='WCS audit reruns')
    for night in ('20240123', '20240124'):
        job(store, 'Callisto', night, cls='P2', state='queued', campaign='WCS audit reruns', est=120.0)
    job(store, 'Artemis', '20190506', kind='framepack', cls='P3', state='queued', est=6.0,
        note='16-bit archive sweep: v2 510 frames in 2 folders', not_before=NOW + dt.timedelta(hours=1))
    store.pause('P3')
    text = notifier(cfg, store).digest(dry_run=True)
    lines = text.splitlines()
    i = lines.index('Processed in the last 24 h: 4 jobs')
    processed = lines[i + 1:i + 3]                   # by class, then the bigger group first
    assert processed[0].startswith('  P2  WCS audit reruns') and \
        processed[0].endswith('3 done · Callisto · Jan 2024 · light curves on 2 of 3 nights')
    assert processed[1].startswith('  P2  reruns') and processed[1].endswith('0 done, 1 failed · Ganymede · 21 Jan 2024')
    j = lines.index("Remaining: 4 jobs, ≈5 h of estimated work, ≈1.7 days at the last 24 h's pace; paused P3")
    remaining = lines[j + 1:j + 3]
    assert remaining[0].startswith('  P2  WCS audit reruns') and \
        remaining[0].endswith('1 running, 2 queued · ≈5 h · Callisto · Jan 2024')
    assert remaining[1].startswith('  P3  16-bit archive sweep') and \
        remaining[1].endswith('1 queued · ≈0 h · Artemis 6 May 2019 · paused · 1 waiting to retry')


def test_empty_queue_says_so(clean, store):
    text = notifier(clean, store).digest(dry_run=True)
    assert 'Processed in the last 24 h: nothing' in text and 'Remaining: nothing queued' in text


def test_q_add_campaign_and_q_campaign(tmp_path, monkeypatch, capsys):
    from jobqueue import cli
    from jobqueue.store import Store
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'q'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    assert cli.main(['add', '--class', 'P3', '--kind', 'sweep', '--campaign', '16-bit sweep', '--', 'true']) == 0
    assert cli.main(['add', '--class', 'P3', '--kind', 'sweep', '--note',
                     'never reduced; raw frames lack WCS (WCS audit 2026-10-16)', '--', 'true']) == 0
    assert cli.main(['campaign', 'WCS audit reruns', '--note', 'wcs audit', '--dry-run']) == 0
    assert 'would update 1 job(s): campaign "WCS audit reruns" (2-2)' in capsys.readouterr().out
    assert cli.main(['campaign', 'WCS audit reruns', '--note', 'wcs audit']) == 0
    assert cli.main(['campaign', '--clear', '--ids', '1']) == 0
    assert cli.main(['campaign', 'x']) == 2
    s = Store(str(tmp_path / 'q' / 'queue.sqlite'))
    try:
        assert 'campaign' not in s.get_job(1)['meta']
        assert s.get_job(2)['meta']['campaign'] == 'WCS audit reruns'
    finally:
        s.close()


# ---------------------------------------------------------------------------------------------- recipients

def test_recipient_list_lives_in_the_queue_config(tmp_path, monkeypatch, capsys):
    from jobqueue import cli
    root = tmp_path / 'q'
    root.mkdir()
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(root))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    conf = root / 'config.json'
    conf.write_text('{"notify": {"transport": "outbox", "to": ["mh2143@cam.ac.uk"]}, "max_cores": 100}')
    assert cli.main(['notify', 'recipients', '--add', 'cj467@cam.ac.uk']) == 0
    out = capsys.readouterr().out
    assert 'mh2143@cam.ac.uk, cj467@cam.ac.uk' in out and 'nothing is sent until it is smtp' in out
    data = json.loads(conf.read_text())
    assert data['notify']['to'] == ['mh2143@cam.ac.uk', 'cj467@cam.ac.uk'] and data['max_cores'] == 100
    assert data['notify']['transport'] == 'outbox'
    assert cli.main(['notify', 'recipients', '--add', 'CJ467@cam.ac.uk']) == 0          # already there
    assert cli.main(['notify', 'recipients', '--remove', 'mh2143@cam.ac.uk']) == 0
    assert cli.main(['notify', 'recipients', '--add', 'not-an-address']) == 2
    assert cli.main(['notify', 'recipients', '--remove', 'nobody@example.org']) == 2
    assert json.loads(conf.read_text())['notify']['to'] == ['cj467@cam.ac.uk']
    capsys.readouterr()
    assert cli.main(['notify', 'recipients']) == 0
    assert 'recipients ({}): cj467@cam.ac.uk'.format(conf) in capsys.readouterr().out


def test_recipients_without_a_config_file(tmp_path, monkeypatch):
    from jobqueue import cli
    monkeypatch.setenv('ORCHARD_QUEUE_ROOT', str(tmp_path / 'q'))
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    assert cli.main(['notify', 'recipients', '--add', 'mh2143@cam.ac.uk']) == 0
    assert json.loads((tmp_path / 'q' / 'config.json').read_text()) == {'notify': {'to': ['mh2143@cam.ac.uk']}}


# ---------------------------------------------------------------------------------------------- never real email

def test_tests_read_no_production_config(cfg):
    """The cfg fixture must not merge /data/SPECULOOSPipeline/queue/config.json (on 2026-10-18 it did, inside the
    container, and a release's test run emailed fake digests through the real account)."""
    assert '_config_path' not in cfg and cfg['notify']['transport'] == 'outbox'


def test_real_smtp_is_refused_inside_a_test(clean, store, tmp_path):
    cfg = smtp_cfg(clean, tmp_path)
    watch(store, 'Ganymede', NIGHT, 'enqueued', rows=10)
    job(store, 'Ganymede', NIGHT, state='failed', signature='boom')
    n = Notifier(cfg, store, now=NOW, usage=usage_ok, log=lambda m: None)     # the real smtplib.SMTP
    n.digest()
    assert {(m['status'], m['error']) for m in n.messages()} == {('failed', 'refused: a test tried to send real email')}
    assert len(outbox(cfg)) == 2


# ---------------------------------------------------------------------------------------------- outages (2026-10-18)

DNS_ERROR = ('ESO authentication failed:  with url: /sso/oidc/token?response_type=id_token+token&grant_type=password'
             '&client_id=clientid&username=***&password=*** (Caused by NameResolutionError("<urllib3.connection.'
             "HTTPSConnection object at 0x7f893eb1a400>: Failed to resolve 'www.eso.org' ([Errno -3] Temporary "
             'failure in name resolution)"))')


def eso_polls(store, start, outcomes, step=30, source='watcher'):
    """Record watcher polls `step` minutes apart from `start`: True answered, False every query failed."""
    t = start
    for ok in outcomes:
        esoreach.record(store, t, ok=ok, error=None if ok else DNS_ERROR, source=source)
        t += dt.timedelta(minutes=step)


def check_at(cfg, store, t, **kw):
    heartbeat(store, t)
    n = notifier(cfg, store, now=t, **kw)
    n.check()
    return n


def keys(n):
    return [i['key'] for i in n.incidents()]


def test_eso_unreachable_opens_after_four_failed_polls_and_resolves(clean, store):
    cfg = clean
    t = dt.datetime(2026, 10, 18, 10, 30)
    eso_polls(store, t, [True, False, False, False])            # answered at 10:30; failed 11:00, 11:30, 12:00
    assert notifier(cfg, store).check(dry_run=True) == ['no conditions']
    eso_polls(store, dt.datetime(2026, 10, 18, 12, 30), [False])
    n = check_at(cfg, store, dt.datetime(2026, 10, 18, 12, 40))
    [inc] = n.incidents()
    assert (inc['key'], inc['title']) == ('eso-unreachable', 'ESO queries failing since 18 Oct 11:00')
    assert inc['detail'].startswith('Every ESO query has failed on the last 4 polls, since 18 Oct 11:00 UTC: 2.0 h '
                                    'of watcher time, the latest by the watcher at 18 Oct 12:30 UTC. ESO last answered '
                                    '18 Oct 10:30 UTC.\nLast error: cannot resolve www.eso.org (Temporary failure in '
                                    'name resolution)\n')
    assert 'getent hosts www.eso.org' in inc['detail'] and 'password' not in inc['detail']
    [alert] = [read(p) for p in outbox(cfg) if '_alert-' in p]
    assert alert['Subject'] == '[orchard] ESO queries failing since 18 Oct 11:00'
    # ESO answers again: one resolved email in the same thread
    eso_polls(store, dt.datetime(2026, 10, 18, 13, 0), [True])
    n = check_at(cfg, store, dt.datetime(2026, 10, 18, 13, 10))
    assert n.incidents() == []
    [resolved] = [read(p) for p in outbox(cfg) if '_resolved-' in p]
    assert resolved['In-Reply-To'] == alert['Message-ID']
    assert [e['message'] for e in store.events(kinds=['eso-recovered'])] == [
        'ESO answered the watcher at 2026-10-18T13:00:00Z after 4 failed polls since 2026-10-18T11:00:00Z']


def test_eso_failures_count_watcher_time_not_the_overnight_break(clean, store):
    cfg = clean
    # hourly polls: two failures are two hours of failure
    eso_polls(store, dt.datetime(2026, 10, 18, 10, 0), [True, False, False], step=60)
    assert keys(notifier(cfg, store)) == [] and \
        [l.split()[0] for l in notifier(cfg, store).check(dry_run=True)] == ['eso-unreachable']
    # answered at 23:00, failed at 23:30 and at 11:00 the next morning: 30 + 60 min, not 12 hours
    eso_polls(store, dt.datetime(2026, 10, 18, 23, 0), [True, False])
    t = dt.datetime(2026, 10, 19, 11, 0)
    eso_polls(store, t, [False])
    assert esoreach.get(store)['failing_minutes'] == 90.0
    heartbeat(store, t)
    assert notifier(cfg, store, now=t).check(dry_run=True) == ['no conditions']


def test_a_look_back_failure_counts_as_a_poll(clean, store):
    eso_polls(store, dt.datetime(2026, 10, 18, 20, 0), [False, False])           # watcher 20:00, 20:30
    eso_polls(store, dt.datetime(2026, 10, 18, 20, 45), [False], source='look-back')
    eso_polls(store, dt.datetime(2026, 10, 18, 21, 0), [False])
    s = esoreach.get(store)
    assert (s['failed_polls'], s['failing_minutes'], s['last_source']) == (4, 30.0 + 15.0 + 15.0, 'watcher')
    # the look-back finishes its run after the watcher's next poll: time never goes backwards
    esoreach.record(store, dt.datetime(2026, 10, 18, 20, 50), ok=False, error='x', source='look-back')
    assert esoreach.get(store)['last_poll'] == '2026-10-18T21:00:00Z'


LC_LINE = ('python /data/SPECULOOSPipeline/queue/releases/20261018T102613Z-f854bc4/src/diff_photometry/SPlightcurve.py '
           '/data/SPECULOOSPipeline/PipelineOutput/v3/{tel}/output/{night}/Sp0439+1615/3310083401849536384_zYJ_{night}_'
           'output.fits -o /data/SPECULOOSPipeline/PipelineOutput/v3/{tel}/output -f zYJ -b 5 -d {night} -a {ap}\n')
LC_OUTCOME = {
    'ok': 'Performing PWV Correction\nPWV correction successful!\n',
    # releases before the LHATPRO fix: no cause, and numpy's message after it
    'nodata': 'Performing PWV Correction\nNo PWV data for: 2025-12-24T12:30:16.819 - 2025-12-25T16:41:39.941\n'
              'PWV Failed: All-NaN slice encountered\n',
    # since the fix: the cause after the window
    'dns': ('Performing PWV Correction\nNo PWV data for: 2025-12-24T12:30:16.819 - 2025-12-25T16:41:39.941 '
            "(ConnectionError: HTTPSConnectionPool(host='archive.eso.org', port=443): Max retries exceeded with url: "
            '/wdb/wdb/asm/lhatpro_paranal/query?wdbo=csv (Caused by NameResolutionError("<urllib3.connection.'
            "HTTPSConnection object at 0x7f1c>: Failed to resolve 'archive.eso.org' ([Errno -3] Temporary failure in "
            'name resolution)")))\nPWV Failed: no LHATPRO PWV data (ConnectionError: ...)\n'),
    'off': 'WARNING: Effective Temperature is too low for reliable PWV correction: 1900K\n'
           'PWV Failed: No T_eff for target\n',
    'other': 'Performing PWV Correction\nPWV Failed: output array is read-only\n',
    'none': 'Traceback (most recent call last):\nValueError: no stars\n',
}


def pwv_log(cfg, name, tel, night, outcomes):
    path = os.path.join(cfg['queue_root'], 'logs', 'jobs', '2026-10', '{}.log'.format(name))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('=== jobqueue job | started 2026-10-18T15:04:55Z\nRunning T1\nPWV correction successful! (not T10)\n')
        for ap, outcome in enumerate(outcomes, 3):
            f.write(LC_LINE.format(tel=tel, night=night, ap=ap))
            f.write('FILES: /data/x/3310083401849536384_zYJ_{}_output.fits\n'.format(night))
            f.write(LC_OUTCOME[outcome])
        f.write('Running T11\n')
    return path


def pwv_job(cfg, store, clock, tel, night, outcomes, minutes=20, state='done'):
    """A pipeline job finishing `minutes` after the clock, whose log holds one light curve per outcome."""
    clock.advance(minutes=minutes)
    jid, _ = store.add_job('P2', 'pipeline', ['./main/ZLP_pipeline.sh'], telescope=tel, night=night)
    store.mark_running(jid, log_path=pwv_log(cfg, 'job{}'.format(jid), tel, night, outcomes))
    if state == 'done':
        store.finish(jid, 0)
    else:
        store.finish(jid, 1, failure_kind='deterministic', signature='boom')
    return jid


def test_pwv_outcome_of_each_light_curve(cfg):
    def scan(*outcomes):
        c = pwvlog.scan_log(pwv_log(cfg, 'x', 'Callisto', '20251224', outcomes))
        return c, pwvlog.verdict(c)
    c, v = scan('ok', 'nodata', 'dns', 'off', 'other', 'none')
    assert c == {'n_lc': 6, 'n_ok': 1, 'n_nodata': 2, 'n_off': 1, 'n_other': 1, 'n_none': 1,
                 'reason': 'cannot resolve archive.eso.org (Temporary failure in name resolution)'}
    assert v == 'ok'                                  # 2 of the 4 that tried had no data: not most
    assert scan('ok', 'nodata', 'nodata', 'off', 'off')[1] == 'nodata'
    assert scan('nodata')[0]['reason'] is None        # no cause logged, and numpy's message is not one
    assert scan('off', 'off')[1] == 'neutral' and scan('none')[1] == 'neutral'
    c, v = scan('other', 'other')
    assert (c['reason'], v) == ('output array is read-only', 'neutral')
    assert scan()[0]['n_lc'] == 0


@pytest.fixture
def outage(clean, store, clock):
    """At NOW: a job that corrected PWV, then two SSO jobs without LHATPRO data, an Artemis job and a job whose
    targets are too cool for PWV; one SSO job short of an incident."""
    cfg = clean
    clock.t = NOW - dt.timedelta(hours=6)
    ok = pwv_job(cfg, store, clock, 'Callisto', '20251223', ['ok'] * 6)
    a = pwv_job(cfg, store, clock, 'Callisto', '20251224', ['nodata'] * 6)
    artemis = pwv_job(cfg, store, clock, 'Artemis', '20261016', ['nodata'] * 6)     # not SSO: not looked at
    b = pwv_job(cfg, store, clock, 'Europa', '20240101', ['ok'] + ['dns'] * 5)
    cool = pwv_job(cfg, store, clock, 'Io', '20240102', ['off'] * 4)                # neither success nor failure
    return {'ok': ok, 'a': a, 'artemis': artemis, 'b': b, 'cool': cool}


def test_pwv_failing_opens_one_incident_naming_the_jobs(outage, cfg, store, clock):
    jobs = outage
    assert notifier(cfg, store).check(dry_run=True) == ['no conditions']
    c = pwv_job(cfg, store, clock, 'Ganymede', '20240103', ['dns'] * 6, state='failed')
    n = notifier(cfg, store)
    n.check()
    [inc] = [i for i in n.incidents() if i['key'] == 'pwv-nodata']
    assert inc['title'] == 'PWV correction failing: no LHATPRO data'
    detail = inc['detail']
    assert detail.startswith("3 SSO pipeline jobs have finished since 18 Oct 02:25 UTC with no LHATPRO data for most "
                             "of their light curves ('No PWV data for:' in the job log): 17 of their 18 light curves "
                             'were written without the PWV correction. The last job that corrected PWV was job {} '
                             '(Callisto 20251223, finished 18 Oct 02:05 UTC).'.format(jobs['ok']))
    assert 'Latest: job {} (Ganymede 20240103, finished 18 Oct 03:45 UTC), 6 of 6 light curves.'.format(c) in detail
    assert 'Cause logged: cannot resolve archive.eso.org (Temporary failure in name resolution)' in detail
    assert 'Jobs: {}, {}, {}\n'.format(jobs['a'], jobs['b'], c) in detail
    [alert] = [read(p) for p in outbox(cfg) if '_alert-' in p]
    assert alert['Subject'] == '[orchard] PWV correction failing: no LHATPRO data'
    # more jobs without data: the same incident, its detail updated, no second email
    d = pwv_job(cfg, store, clock, 'Io', '20240104', ['nodata'] * 2)
    n = check_at(cfg, store, NOW + dt.timedelta(minutes=30))
    [inc] = [i for i in n.incidents() if i['key'] == 'pwv-nodata']
    assert inc['detail'].startswith('4 SSO pipeline jobs') and ', {}\n'.format(d) in inc['detail']
    assert len([p for p in outbox(cfg) if '_alert-' in p]) == 1


def test_pwv_incident_resolves_only_when_a_later_job_corrects_pwv(outage, cfg, store, clock):
    pwv_job(cfg, store, clock, 'Ganymede', '20240103', ['dns'] * 6)
    n = check_at(cfg, store, NOW)
    assert 'pwv-nodata' in keys(n)
    # no SSO job for four days (the queue paused, say): still open, not resolved by the jobs ageing out
    later = NOW + dt.timedelta(days=4)
    n = check_at(cfg, store, later)
    assert 'pwv-nodata' in keys(n)
    # a job mostly without data, though one light curve was corrected: still failing
    clock.t = later
    pwv_job(cfg, store, clock, 'Europa', '20240105', ['ok', 'dns', 'dns'])
    assert 'pwv-nodata' in keys(check_at(cfg, store, later + dt.timedelta(hours=1)))
    pwv_job(cfg, store, clock, 'Europa', '20240106', ['ok'] * 6)
    n = check_at(cfg, store, later + dt.timedelta(hours=2))
    assert 'pwv-nodata' not in keys(n)
    [resolved] = [read(p) for p in outbox(cfg) if '_resolved-' in p and 'PWV' in read(p)['Subject']]
    assert resolved['Subject'] == 'Re: [orchard] PWV correction failing: no LHATPRO data'


def test_pwv_logs_are_read_once_and_again_after_a_rerun(outage, cfg, store, clock, monkeypatch):
    read_logs = []
    scan = pwvlog.scan_log
    monkeypatch.setattr(pwvlog, 'scan_log', lambda path: read_logs.append(path) or scan(path))
    notifier(cfg, store).check(dry_run=True)
    assert len(read_logs) == 4                                   # the four SSO jobs; Artemis is not read
    assert store.conn.execute('SELECT COUNT(*) FROM job_pwv').fetchone()[0] == 0     # a dry run keeps nothing
    check_at(cfg, store, NOW)
    check_at(cfg, store, NOW + dt.timedelta(minutes=30))
    assert len(read_logs) == 8                                   # read once more for real, then never again
    # job a fails, is requeued and runs again, this time with PWV: its new log is read
    a = outage['a']
    clock.t = NOW + dt.timedelta(minutes=40)
    store.conn.execute("UPDATE jobs SET state = 'failed' WHERE id = ?", (a,))
    store.requeue([a])
    store.mark_running(a, log_path=pwv_log(cfg, 'rerun', 'Callisto', '20251224', ['ok'] * 6))
    store.finish(a, 0)
    check_at(cfg, store, NOW + dt.timedelta(hours=1))
    assert read_logs[8:] == [store.get_job(a)['log_path']]
    row = store.conn.execute('SELECT verdict, n_ok FROM job_pwv WHERE job_id = ?', (a,)).fetchone()
    assert tuple(row) == ('ok', 6)


def test_pwv_reads_at_most_pwv_max_scans_logs_a_run_newest_first(outage, cfg, store):
    jobs = outage
    cfg = dict(cfg, notify=dict(cfg['notify'], pwv_max_scans=2))

    def read():
        return {r[0] for r in store.conn.execute('SELECT job_id FROM job_pwv')}
    check_at(cfg, store, NOW)
    assert read() == {jobs['b'], jobs['cool']}
    check_at(cfg, store, NOW + dt.timedelta(minutes=30))
    assert read() == {jobs['ok'], jobs['a'], jobs['b'], jobs['cool']}


def test_digest_counts_pwv_outcomes(outage, cfg, store, clock):
    pwv_job(cfg, store, clock, 'Ganymede', '20240103', ['dns'] * 6)
    text = notifier(cfg, store).digest(dry_run=True)
    assert ('PWV, 28 light curves from 5 SSO pipeline jobs in the last 24 h: 7 corrected, 17 without LHATPRO data '
            '(3 jobs), 4 off by design (Teff or telescope)') in text
    assert 'PWV correction failing: no LHATPRO data' in text


def test_digest_says_nothing_of_pwv_without_sso_jobs(clean, store):
    assert 'PWV,' not in notifier(clean, store).digest(dry_run=True)


class DnsDownSmtp:
    """smtplib.SMTP while appct cannot resolve names: smtp.gmail.com fails before anything is sent."""

    def __init__(self, host, port, timeout=None):
        import socket
        raise socket.gaierror(-3, 'Temporary failure in name resolution')


def test_mail_that_cannot_be_sent_stays_in_the_outbox_and_the_next_digest_says_so(clean, store, tmp_path):
    cfg = smtp_cfg(clean, tmp_path)
    eso_polls(store, NOW - dt.timedelta(hours=2), [False] * 4)
    # the alert: written to the outbox, not sent, tried again on the next run
    n = check_at(cfg, store, NOW + dt.timedelta(minutes=10), smtp_factory=DnsDownSmtp)
    [inc] = n.incidents()
    assert inc['key'] == 'eso-unreachable' and inc['alerted_at'] is None
    [msg] = n.messages()
    assert (msg['status'], msg['error']) == ('failed',
                                             'smtp.example.org:2525: [Errno -3] Temporary failure in name resolution')
    assert [read(p)['Subject'] for p in outbox(cfg)] == ['[orchard] ESO queries failing since 18 Oct 05:45']
    # the next morning's digest cannot be sent either; its outbox copy lists the alert
    day2 = NOW + dt.timedelta(days=1)
    heartbeat(store, day2)
    out = notifier(cfg, store, now=day2, smtp_factory=DnsDownSmtp).digest()
    assert 'failed via smtp' in out
    [digest2] = [read(p).get_content() for p in outbox(cfg) if '_digest-' in p]
    assert 'Sending failed for 1 message since 12 Oct 07:45 UTC (2 attempts); each is in {}:'.format(
        os.path.join(cfg['queue_root'], 'notify', 'outbox')) in digest2
    assert '  • [orchard] ESO queries failing since 18 Oct 05:45 — not sent (2 attempts)' in digest2
    assert '  Last error: smtp.example.org:2525: [Errno -3] Temporary failure in name resolution' in digest2
    # mail works again: the alert goes out, and the digest says what was missed
    FakeSmtp.sent = []
    day3 = NOW + dt.timedelta(days=2)
    heartbeat(store, day3)
    notifier(cfg, store, now=day3, smtp_factory=FakeSmtp).digest()
    assert [m['Subject'] for m, *_ in FakeSmtp.sent] == ['[orchard] ESO queries failing since 18 Oct 05:45',
                                                          '[orchard] Night of 18 Oct 2026: 1 needs attention']
    body = FakeSmtp.sent[-1][0].get_content()
    assert 'Sending failed for 2 messages since' in body
    assert '  • [orchard] ESO queries failing since 18 Oct 05:45 — sent on a later attempt' in body
    assert '  • [orchard] Night of 17 Oct 2026: 1 needs attention — not sent' in body
    # once a digest is through, the next one starts afresh
    day4 = NOW + dt.timedelta(days=3)
    heartbeat(store, day4)
    assert 'Sending failed' not in notifier(cfg, store, now=day4, smtp_factory=FakeSmtp).digest(dry_run=True)
