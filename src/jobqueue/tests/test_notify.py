"""The notifier against a fake production tree: digest lines, incident rules, threading, acks, delivery."""
import datetime as dt
import email
import email.policy
import glob
import json
import os

import pytest

from jobqueue import digest as D
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
