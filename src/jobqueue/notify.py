"""The notifier: the pipeline's only sender of email.

Each run (`q notify check` from cron every 30 minutes, `q notify digest` every morning) works out the conditions
that need a person (digest.conditions) and keeps one incident per condition in the queue database. An incident
gets one email when it opens and one when it resolves, in the same thread; in between it is listed in every
digest. `q notify ack` silences one until a date. The digest goes out every morning even when all is well, so
its absence is a signal too, and `notify.heartbeat_url` lets an outside service raise the alarm if it stops.

Delivery: every message is written to <queue_root>/notify/outbox as an .eml file. Only with notify.transport =
'smtp' is it also sent, to the recipient list (notify.to, kept in the queue's config.json by
`q notify recipients`), with EMAIL_ADDRESS / EMAIL_PASSWORD / SMTP_SERVER read from notify.env_file. If the SMTP
settings are missing or sending fails, the message still lands in the outbox and the failure is recorded; an alert
that failed is tried again on the next run.
"""
import datetime as dt
import email.utils
import os
import re
import smtplib
import urllib.request
from email.message import EmailMessage

from . import digest as D
from .config import queue_paths
from .util import iso, last_night, night_str, parse_iso, read_json, utcnow, write_json_atomic

SCHEMA = """
CREATE TABLE IF NOT EXISTS notify_incidents (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    key TEXT NOT NULL,
    kind TEXT NOT NULL,
    telescope TEXT,
    night TEXT,
    title TEXT NOT NULL,
    detail TEXT,
    state TEXT NOT NULL,
    opened_at TEXT NOT NULL,
    last_seen TEXT NOT NULL,
    resolved_at TEXT,
    alerted_at TEXT,
    resolve_alerted_at TEXT,
    thread_id TEXT,
    ack_until TEXT,
    ack_note TEXT
);
CREATE UNIQUE INDEX IF NOT EXISTS notify_incidents_open ON notify_incidents(key) WHERE state = 'open';
CREATE TABLE IF NOT EXISTS notify_messages (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    ts TEXT NOT NULL,
    kind TEXT NOT NULL,
    incident_id INTEGER,
    subject TEXT NOT NULL,
    message_id TEXT NOT NULL,
    transport TEXT NOT NULL,
    path TEXT,
    status TEXT NOT NULL,
    error TEXT
);
"""

MSGID_DOMAIN = 'orchard.appct'
DELIVERED = ('written', 'sent')


class NotifyError(RuntimeError):
    pass


ADDRESS = re.compile(r'^[^@\s,;<>"]+@[^@\s,;<>"]+\.[^@\s,;<>"]+$')


def config_path(cfg):
    """The JSON file the queue's settings come from: --config / $ORCHARD_QUEUE_CONFIG, else <queue_root>/config.json."""
    return cfg.get('_config_path') or os.path.join(cfg['queue_root'], 'config.json')


def edit_recipients(cfg, add=(), remove=()):
    """Add or remove addresses on the recipient list (notify.to in the queue's config file, created if missing;
    every other setting in it is kept). Returns (the new list, the file)."""
    path = config_path(cfg)
    data = read_json(path, default=None) if os.path.exists(path) else {}
    if not isinstance(data, dict):
        raise NotifyError('{} is not a JSON object; fix it by hand first'.format(path))
    section = data.setdefault('notify', {})
    current = list(section.get('to', cfg['notify'].get('to') or []))
    for address in add:
        if not ADDRESS.match(address):
            raise NotifyError('not an email address: {!r}'.format(address))
        if address.lower() not in [a.lower() for a in current]:
            current.append(address)
    for address in remove:
        kept = [a for a in current if a.lower() != address.lower()]
        if len(kept) == len(current):
            raise NotifyError('{} is not on the recipient list'.format(address))
        current = kept
    section['to'] = current
    write_json_atomic(path, data)
    return current, path


def read_smtp_settings(env_path):
    """EMAIL_ADDRESS, EMAIL_PASSWORD and SMTP_SERVER from a dotenv file, without touching the environment."""
    if not env_path:
        raise NotifyError('notify.env_file is not set')
    found = {}
    try:
        with open(env_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                key = key.strip()
                if key.startswith('export '):
                    key = key[7:].strip()
                value = value.strip()
                if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
                    value = value[1:-1]
                if key in ('EMAIL_ADDRESS', 'EMAIL_PASSWORD', 'SMTP_SERVER'):
                    found[key] = value
    except OSError as e:
        raise NotifyError('cannot read {}: {}'.format(env_path, e.strerror))
    if not found.get('EMAIL_ADDRESS') or not found.get('EMAIL_PASSWORD'):
        raise NotifyError('EMAIL_ADDRESS / EMAIL_PASSWORD not set in {}'.format(env_path))
    host, _, port = (found.get('SMTP_SERVER') or 'smtp.gmail.com:587').partition(':')
    return {'user': found['EMAIL_ADDRESS'], 'password': found['EMAIL_PASSWORD'], 'host': host,
            'port': int(port or 587)}


def build_message(subject, body, kind, now, sender, to, thread_id=None):
    m = EmailMessage()
    m['Subject'] = subject
    m['From'] = 'orchard <{}>'.format(sender)
    m['To'] = ', '.join(to) if to else 'outbox@localhost'
    m['Date'] = email.utils.format_datetime(now.replace(tzinfo=dt.timezone.utc))
    m['Message-ID'] = email.utils.make_msgid(domain=MSGID_DOMAIN)
    if thread_id:
        m['In-Reply-To'] = thread_id
        m['References'] = thread_id
    m['X-Orchard-Kind'] = kind
    # quoted-printable keeps the wire ASCII whatever the text holds (the old download email died on a ✓)
    m.set_content(body, cte='quoted-printable')
    return m


class Delivery:
    """Writes every message to the outbox; with notify.transport = 'smtp' also sends it."""

    def __init__(self, cfg, outbox_only=False, smtp_factory=smtplib.SMTP):
        n = cfg['notify']
        if n['transport'] not in ('outbox', 'smtp'):
            raise NotifyError("notify.transport must be 'outbox' or 'smtp', not {!r}".format(n['transport']))
        self.outbox = queue_paths(cfg)['outbox']
        self.smtp = n['transport'] == 'smtp' and not outbox_only
        self.factory = smtp_factory
        self.to = [a for a in (n.get('to') or []) if a]
        self.settings, self.problem = None, None
        if self.smtp:
            try:
                if not self.to:
                    raise NotifyError('notify.to is empty')
                self.settings = read_smtp_settings(n.get('env_file'))
            except NotifyError as e:
                self.problem = 'not sent, SMTP is not configured: {}'.format(e)

    @property
    def name(self):
        return 'smtp' if self.smtp else 'outbox'

    @property
    def sender(self):
        return self.settings['user'] if self.settings else 'orchard-notifier@localhost'

    def write(self, msg, slug, now):
        os.makedirs(self.outbox, exist_ok=True)
        base = os.path.join(self.outbox, '{}_{}'.format(now.strftime('%Y%m%dT%H%M%SZ'), slug))
        path, k = base + '.eml', 1
        while os.path.exists(path):
            k += 1
            path = '{}_{}.eml'.format(base, k)
        tmp = path + '.tmp'
        with open(tmp, 'wb') as f:
            f.write(msg.as_bytes())
        os.replace(tmp, path)
        return path

    def deliver(self, msg, slug, now):
        """(outbox path, status, error): status 'written' (outbox only), 'sent', or 'failed'."""
        path = self.write(msg, slug, now)
        if not self.smtp:
            return path, 'written', None
        if self.problem:
            return path, 'failed', self.problem
        if self.factory is smtplib.SMTP and os.environ.get('PYTEST_CURRENT_TEST'):
            return path, 'failed', 'refused: a test tried to send real email'
        s = self.settings
        try:
            with self.factory(s['host'], s['port'], timeout=60) as conn:
                conn.starttls()
                conn.login(s['user'], s['password'])
                conn.send_message(msg, from_addr=s['user'], to_addrs=self.to)
        except (smtplib.SMTPException, OSError) as e:
            return path, 'failed', str(e).replace(s['password'], '***')[:500]
        return path, 'sent', None


def _slug(text, limit=48):
    out = ''.join(c if c.isalnum() else '-' for c in text.lower())
    while '--' in out:
        out = out.replace('--', '-')
    return out.strip('-')[:limit] or 'message'


class Notifier:
    def __init__(self, cfg, store, now=None, outbox_only=False, smtp_factory=smtplib.SMTP, opener=None,
                 usage=D.disk_usage, log=print):
        self.cfg, self.store, self.log = cfg, store, log
        self.now = now or utcnow()
        self.n = cfg['notify']
        self.usage = usage
        self.opener = opener or urllib.request.urlopen
        store.conn.executescript(SCHEMA)
        self.delivery = Delivery(cfg, outbox_only=outbox_only, smtp_factory=smtp_factory)

    # ------------------------------------------------------------------ evaluating

    def evaluate(self, night=None, statuses=True):
        ctx = D.Context(self.cfg, self.store, self.now)
        night = night or night_str(last_night(self.now))
        st = [D.night_status(ctx, t, night) for t in self.cfg['telescopes']] if statuses else []
        q = D.queue_summary(ctx)
        if statuses:
            q['work'] = D.work_summary(ctx, q)
        disk_list = D.disks(ctx, usage=self.usage)
        return night, st, q, disk_list, D.conditions(ctx, q, disk_list)

    def incidents(self, state='open', limit=None):
        q = 'SELECT * FROM notify_incidents WHERE state = ? ORDER BY id' + (' DESC LIMIT {:d}'.format(limit)
                                                                             if limit else '')
        return [dict(r) for r in self.store.conn.execute(q, (state,))]

    def sync(self, conds):
        """Open an incident for each new condition, refresh the ones still true, resolve the ones that cleared.
        One transaction, reading the open incidents inside it, so a check and a digest that overlap cannot both
        open the same one."""
        now = iso(self.now)
        current = {c['key']: c for c in conds}
        opened, resolved = [], []
        with self.store.tx() as c:
            open_rows = {r['key']: dict(r) for r in c.execute("SELECT * FROM notify_incidents WHERE state = 'open'")}
            for key, cond in current.items():
                if key in open_rows:
                    c.execute('UPDATE notify_incidents SET last_seen = ?, detail = ? WHERE id = ?',
                              (now, cond['detail'], open_rows[key]['id']))
                else:
                    c.execute('INSERT INTO notify_incidents (key, kind, telescope, night, title, detail, state, '
                              "opened_at, last_seen) VALUES (?, ?, ?, ?, ?, ?, 'open', ?, ?)",
                              (key, cond['kind'], cond['telescope'], cond['night'], cond['title'], cond['detail'],
                               now, now))
                    opened.append(cond)
            for key, row in open_rows.items():
                if key not in current:
                    c.execute("UPDATE notify_incidents SET state = 'resolved', resolved_at = ? WHERE id = ?",
                              (now, row['id']))
                    resolved.append(row)
        for cond in opened:
            self.store.event('incident-opened', source='notify', telescope=cond['telescope'], night=cond['night'],
                             message=cond['title'])
        for row in resolved:
            self.store.event('incident-resolved', source='notify', telescope=row['telescope'], night=row['night'],
                             message=row['title'])
        return opened, resolved

    # ------------------------------------------------------------------ delivering

    def _deliver(self, kind, subject, body, incident_id=None, thread_id=None):
        msg = build_message(subject, body, kind, self.now, self.delivery.sender, self.delivery.to, thread_id)
        path, status, error = self.delivery.deliver(msg, '{}-{}'.format(kind, _slug(subject)), self.now)
        self.store.conn.execute(
            'INSERT INTO notify_messages (ts, kind, incident_id, subject, message_id, transport, path, status, error) '
            'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
            (iso(self.now), kind, incident_id, subject, msg['Message-ID'], self.delivery.name, path, status, error))
        if error:
            self.log('{} "{}" {}: {}'.format(kind, subject, status, error))
        return status, msg['Message-ID']

    def _claim(self, incident_id, column):
        """Mark an alert as being sent, atomically, so two overlapping runs never both send it."""
        cur = self.store.conn.execute('UPDATE notify_incidents SET {0} = ? WHERE id = ? AND {0} IS NULL'.format(column),
                                      (iso(self.now), incident_id))
        return cur.rowcount == 1

    def _release(self, incident_id, column):
        self.store.conn.execute('UPDATE notify_incidents SET {} = NULL WHERE id = ?'.format(column), (incident_id,))

    def send_alerts(self):
        """One email when an incident opens (unless acknowledged), one when it resolves, in the same thread."""
        now = iso(self.now)
        prefix = self.n['subject_prefix']
        out = []
        for r in self.store.conn.execute(
                "SELECT * FROM notify_incidents WHERE state = 'open' AND alerted_at IS NULL ORDER BY id").fetchall():
            if (r['ack_until'] and r['ack_until'] > now) or not self._claim(r['id'], 'alerted_at'):
                continue
            body = ('{}\n\nOpened {} UTC. This thread gets one more email when it resolves; until then it is listed '
                    'in every morning digest.\nTo silence it: q notify ack {} --until <YYYY-MM-DD> --note <text>\n'
                    .format(r['detail'], D.fmt_time(r['opened_at']), r['id']))
            status, msgid = self._deliver('alert', '{} {}'.format(prefix, r['title']), body, incident_id=r['id'])
            if status in DELIVERED:
                self.store.conn.execute('UPDATE notify_incidents SET thread_id = ? WHERE id = ?', (msgid, r['id']))
            else:
                self._release(r['id'], 'alerted_at')          # try again on the next run
            out.append('alert {}: {} ({})'.format(r['id'], r['title'], status))
        for r in self.store.conn.execute(
                "SELECT * FROM notify_incidents WHERE state = 'resolved' AND alerted_at IS NOT NULL "
                'AND thread_id IS NOT NULL AND resolve_alerted_at IS NULL ORDER BY id').fetchall():
            if not self._claim(r['id'], 'resolve_alerted_at'):
                continue
            opened, resolved = parse_iso(r['opened_at']), parse_iso(r['resolved_at'])
            body = ('Resolved {} UTC, {} after it opened.\n\nLast seen before it cleared:\n{}\n'.format(
                D.fmt_time(r['resolved_at']), _duration(resolved - opened), r['detail']))
            status, _ = self._deliver('resolved', 'Re: {} {}'.format(prefix, r['title']), body, incident_id=r['id'],
                                      thread_id=r['thread_id'])
            if status not in DELIVERED:
                self._release(r['id'], 'resolve_alerted_at')
            out.append('resolved {}: {} ({})'.format(r['id'], r['title'], status))
        return out

    def heartbeat(self):
        url = self.n.get('heartbeat_url')
        if not url:
            return None
        try:
            with self.opener(url, timeout=15) as r:
                r.read(256)
        except Exception as e:                               # noqa: BLE001 - any failure is just reported
            self.log('heartbeat failed: {}'.format(e))
            return 'failed'
        self.store.kv_set('notify_heartbeat', iso(self.now))
        return 'ok'

    # ------------------------------------------------------------------ the two cron entry points

    def check(self, dry_run=False):
        """Bring incidents up to date and deliver their alerts. dry_run: list the conditions, change nothing."""
        _, _, _, _, conds = self.evaluate(statuses=False)
        if dry_run:
            return ['{}  {}'.format(c['key'], c['title']) for c in conds] or ['no conditions']
        opened, resolved = self.sync(conds)
        lines = ['{} notify check: {} condition(s), {} opened, {} resolved'.format(
            iso(self.now), len(conds), len(opened), len(resolved))]
        return lines + self.send_alerts()

    def digest(self, night=None, dry_run=False):
        """The morning digest of `night` (default last night), after bringing incidents up to date. dry_run: return
        the text only; no incident changes, nothing delivered."""
        night, statuses, q, disk_list, conds = self.evaluate(night)
        lines = []
        if dry_run:
            known = {i['key']: i for i in self.incidents()}
            open_now = [known.get(c['key']) or dict(c, opened_at=iso(self.now), id=None) for c in conds]
            resolved = []
        else:
            self.sync(conds)
            lines = self.send_alerts()
            open_now = self.incidents()
            since = self.store.kv_get('notify_last_digest') or ''
            resolved = [dict(r) for r in self.store.conn.execute(
                "SELECT * FROM notify_incidents WHERE state = 'resolved' AND resolved_at > ? ORDER BY id", (since,))]
        report = os.path.join(queue_paths(self.cfg)['reports'], 'report_{}.md'.format(self.now.strftime('%Y%m%d')))
        shutter = D.shutter_notes(D.Context(self.cfg, self.store, self.now))
        subject, body = D.render_digest(self.cfg, night, statuses, q, disk_list, open_now, resolved, self.now,
                                        report_path=report if os.path.exists(report) else None, shutter=shutter)
        if dry_run:
            return subject + '\n\n' + body
        status, _ = self._deliver('digest', subject, body)
        if status in DELIVERED:
            self.store.kv_set('notify_last_digest', iso(self.now))
        beat = self.heartbeat()
        lines.insert(0, '{} digest "{}" {} via {}{}'.format(iso(self.now), subject, status, self.delivery.name,
                                                           '; heartbeat {}'.format(beat) if beat else ''))
        return '\n'.join(lines)

    # ------------------------------------------------------------------ people

    def _find(self, ident):
        r = self.store.conn.execute('SELECT * FROM notify_incidents WHERE id = ?', (ident,)).fetchone() \
            if str(ident).isdigit() else \
            self.store.conn.execute("SELECT * FROM notify_incidents WHERE key = ? AND state = 'open'",
                                    (ident,)).fetchone()
        if r is None:
            raise NotifyError('no incident {}'.format(ident))
        return dict(r)

    def ack(self, ident, until, note=''):
        r = self._find(ident)
        day = dt.datetime.strptime(until, '%Y-%m-%d')
        self.store.conn.execute('UPDATE notify_incidents SET ack_until = ?, ack_note = ? WHERE id = ?',
                                (iso(day), note or None, r['id']))
        self.store.event('incident-acked', source='notify', telescope=r['telescope'], night=r['night'],
                         message='{} until {}: {}'.format(r['title'], until, note or '-'))
        return r

    def unack(self, ident):
        r = self._find(ident)
        self.store.conn.execute('UPDATE notify_incidents SET ack_until = NULL, ack_note = NULL WHERE id = ?',
                                (r['id'],))
        return r

    def messages(self, last=20):
        return [dict(r) for r in self.store.conn.execute(
            'SELECT * FROM notify_messages ORDER BY id DESC LIMIT {:d}'.format(int(last)))]


def _duration(td):
    hours = td.total_seconds() / 3600.0
    return '{:.0f} min'.format(hours * 60) if hours < 2 else '{:.0f} h'.format(hours) if hours < 48 else \
        '{:.0f} days'.format(hours / 24)


def format_incidents(rows, now):
    if not rows:
        return 'no incidents'
    out = []
    for r in rows:
        flags = []
        if r['state'] == 'resolved':
            flags.append('resolved {}'.format(D.fmt_time(r['resolved_at'])))
        if r.get('alerted_at'):
            flags.append('alerted')
        if r.get('ack_until') and r['ack_until'] > iso(now):
            flags.append('acknowledged until {}{}'.format(r['ack_until'][:10], ': ' + r['ack_note']
                                                          if r.get('ack_note') else ''))
        out.append('{:>4}  opened {}  {}{}'.format(r['id'], D.fmt_time(r['opened_at']), r['title'],
                                                  '  [{}]'.format('; '.join(flags)) if flags else ''))
        out.append('      {}'.format(r['key']))
    return '\n'.join(out)
