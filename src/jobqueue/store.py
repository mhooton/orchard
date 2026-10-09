"""SQLite queue store: jobs, an append-only event log, watcher and look-back state.

The database must live on appct's local disk and be opened only from inside
orchard-server on appct (Store refuses a network filesystem). WAL mode lets the
dispatcher, the CLI, the watcher and running fetch jobs share it.

Job states: queued -> running -> done | failed | cancelled. A transient failure
goes back to queued with a not_before time; so does a job whose process vanished
(once). Every transition is recorded in the events table.
"""
import contextlib
import datetime as dt
import json
import os
import sqlite3

from . import ACTIVE_STATES, CLASSES
from .util import ensure_local_fs, iso, parse_iso, utcnow

SCHEMA_VERSION = 1

SCHEMA = """
CREATE TABLE IF NOT EXISTS jobs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    class TEXT NOT NULL,
    kind TEXT NOT NULL,
    telescope TEXT,
    night TEXT,
    targets TEXT NOT NULL DEFAULT '',
    argv TEXT NOT NULL,
    cwd TEXT,
    env TEXT NOT NULL DEFAULT '{}',
    cores INTEGER NOT NULL,
    disk_heavy INTEGER NOT NULL DEFAULT 1,
    est_minutes REAL NOT NULL,
    locks TEXT NOT NULL DEFAULT '[]',
    depends_on INTEGER,
    dep_requires_success INTEGER NOT NULL DEFAULT 0,
    priority REAL NOT NULL DEFAULT 0,
    dedupe_key TEXT,
    state TEXT NOT NULL DEFAULT 'queued',
    attempts INTEGER NOT NULL DEFAULT 0,
    max_retries INTEGER NOT NULL DEFAULT 2,
    vanished INTEGER NOT NULL DEFAULT 0,
    not_before TEXT,
    created_at TEXT NOT NULL,
    started_at TEXT,
    finished_at TEXT,
    exit_code INTEGER,
    log_path TEXT,
    status_path TEXT,
    code_sha TEXT,
    failure_kind TEXT,
    failure_signature TEXT,
    pid INTEGER,
    proc_start TEXT,
    kill_sent_at TEXT,
    kill_reason TEXT,
    cancel_requested INTEGER NOT NULL DEFAULT 0,
    shadow INTEGER NOT NULL DEFAULT 0,
    source TEXT NOT NULL DEFAULT 'manual',
    meta TEXT NOT NULL DEFAULT '{}',
    note TEXT
);
CREATE UNIQUE INDEX IF NOT EXISTS jobs_active_dedupe ON jobs(dedupe_key)
    WHERE dedupe_key IS NOT NULL AND state IN ('queued', 'running');
CREATE INDEX IF NOT EXISTS jobs_state ON jobs(state, class);
CREATE INDEX IF NOT EXISTS jobs_tel_night ON jobs(telescope, night);
CREATE INDEX IF NOT EXISTS jobs_signature ON jobs(failure_signature);

CREATE TABLE IF NOT EXISTS events (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    ts TEXT NOT NULL,
    job_id INTEGER,
    source TEXT NOT NULL,
    kind TEXT NOT NULL,
    telescope TEXT,
    night TEXT,
    message TEXT,
    data TEXT
);
CREATE INDEX IF NOT EXISTS events_ts ON events(ts);
CREATE INDEX IF NOT EXISTS events_job ON events(job_id);

CREATE TABLE IF NOT EXISTS kv (
    key TEXT PRIMARY KEY,
    value TEXT,
    updated_at TEXT
);

CREATE TABLE IF NOT EXISTS watch (
    telescope TEXT NOT NULL,
    night TEXT NOT NULL,
    state TEXT NOT NULL,
    eso_rows INTEGER,
    eso_last_mod TEXT,
    transfer_count INTEGER,
    first_seen TEXT,
    last_change TEXT,
    last_poll TEXT,
    polls_unchanged INTEGER NOT NULL DEFAULT 0,
    ready_at TEXT,
    ready_reason TEXT,
    enqueued_at TEXT,
    download_job INTEGER,
    pipeline_job INTEGER,
    noted INTEGER NOT NULL DEFAULT 0,
    note TEXT,
    data TEXT,
    PRIMARY KEY (telescope, night)
);

CREATE TABLE IF NOT EXISTS lookback (
    telescope TEXT NOT NULL,
    night TEXT NOT NULL,
    checked_at TEXT,
    eso_rows INTEGER,
    disk_files INTEGER,
    missing_rows INTEGER,
    missing_science INTEGER,
    products TEXT,
    fetches INTEGER NOT NULL DEFAULT 0,
    reruns INTEGER NOT NULL DEFAULT 0,
    last_action TEXT,
    data TEXT,
    PRIMARY KEY (telescope, night)
);
"""

JOB_JSON = {'argv': list, 'env': dict, 'locks': list, 'meta': dict}
WATCH_COLS = ('telescope', 'night', 'state', 'eso_rows', 'eso_last_mod', 'transfer_count', 'first_seen',
              'last_change', 'last_poll', 'polls_unchanged', 'ready_at', 'ready_reason', 'enqueued_at',
              'download_job', 'pipeline_job', 'noted', 'note', 'data')
LOOKBACK_COLS = ('telescope', 'night', 'checked_at', 'eso_rows', 'disk_files', 'missing_rows',
                 'missing_science', 'products', 'fetches', 'reruns', 'last_action', 'data')


def _job(row):
    if row is None:
        return None
    d = dict(row)
    for k, typ in JOB_JSON.items():
        d[k] = json.loads(d[k]) if d.get(k) else typ()
    return d


def _plain(row, json_cols=('data',)):
    if row is None:
        return None
    d = dict(row)
    for k in json_cols:
        if k in d:
            d[k] = json.loads(d[k]) if d[k] else {}
    return d


class Store:
    def __init__(self, path, now_fn=utcnow, check_local=True):
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)
        if check_local:
            ensure_local_fs(directory)
        self.path = path
        self.now = now_fn
        self.conn = sqlite3.connect(path, timeout=60, isolation_level=None)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute('PRAGMA busy_timeout=60000')
        self.conn.execute('PRAGMA journal_mode=WAL')
        self.conn.execute('PRAGMA synchronous=NORMAL')
        self.conn.executescript(SCHEMA)
        if self.kv_get('schema_version') is None:
            self.kv_set('schema_version', str(SCHEMA_VERSION))

    def close(self):
        self.conn.close()

    @contextlib.contextmanager
    def tx(self):
        self.conn.execute('BEGIN IMMEDIATE')
        try:
            yield self.conn
        except BaseException:
            self.conn.execute('ROLLBACK')
            raise
        self.conn.execute('COMMIT')

    # ---------------------------------------------------------------- events and key/value

    def _event(self, c, kind, job_id=None, source='queue', telescope=None, night=None, message=None, data=None):
        c.execute('INSERT INTO events (ts, job_id, source, kind, telescope, night, message, data) '
                  'VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                  (iso(self.now()), job_id, source, kind, telescope, night, message,
                   json.dumps(data, sort_keys=True, default=str) if data is not None else None))

    def event(self, kind, **kw):
        with self.tx() as c:
            self._event(c, kind, **kw)

    def events(self, since=None, until=None, kinds=None, job_id=None, source=None, limit=None):
        q, args = 'SELECT * FROM events WHERE 1=1', []
        if since:
            q += ' AND ts >= ?'
            args.append(iso(since) if isinstance(since, dt.datetime) else since)
        if until:
            q += ' AND ts < ?'
            args.append(iso(until) if isinstance(until, dt.datetime) else until)
        if kinds:
            q += ' AND kind IN ({})'.format(','.join('?' * len(kinds)))
            args.extend(kinds)
        if job_id is not None:
            q += ' AND job_id = ?'
            args.append(job_id)
        if source:
            q += ' AND source = ?'
            args.append(source)
        q += ' ORDER BY id'
        if limit:
            q += ' LIMIT {:d}'.format(limit)
        return [_plain(r) for r in self.conn.execute(q, args)]

    def kv_get(self, key, default=None):
        r = self.conn.execute('SELECT value FROM kv WHERE key = ?', (key,)).fetchone()
        return r['value'] if r else default

    def kv_set(self, key, value):
        self.conn.execute('INSERT OR REPLACE INTO kv (key, value, updated_at) VALUES (?, ?, ?)',
                          (key, value, iso(self.now())))

    def kv_delete(self, key):
        self.conn.execute('DELETE FROM kv WHERE key = ?', (key,))

    def pause(self, cls, reason='', source='cli'):
        classes = CLASSES if cls == 'all' else (cls,)
        with self.tx() as c:
            for k in classes:
                if k not in CLASSES:
                    raise ValueError('unknown class {}'.format(k))
                c.execute('INSERT OR REPLACE INTO kv (key, value, updated_at) VALUES (?, ?, ?)',
                          ('paused:' + k, reason or 'paused', iso(self.now())))
                self._event(c, 'paused', source=source, message='{} paused: {}'.format(k, reason or '-'))

    def resume(self, cls, source='cli'):
        classes = CLASSES if cls == 'all' else (cls,)
        with self.tx() as c:
            for k in classes:
                c.execute('DELETE FROM kv WHERE key = ?', ('paused:' + k,))
                self._event(c, 'resumed', source=source, message='{} resumed'.format(k))

    def paused(self):
        rows = self.conn.execute("SELECT key FROM kv WHERE key LIKE 'paused:%'").fetchall()
        return {r['key'].split(':', 1)[1] for r in rows}

    def set_drain(self, on, source='cli'):
        with self.tx() as c:
            if on:
                c.execute("INSERT OR REPLACE INTO kv (key, value, updated_at) VALUES ('drain', '1', ?)",
                          (iso(self.now()),))
            else:
                c.execute("DELETE FROM kv WHERE key = 'drain'")
            self._event(c, 'drain' if on else 'undrain', source=source)

    def draining(self):
        return self.kv_get('drain') == '1'

    # ---------------------------------------------------------------- jobs

    def add_job(self, cls, kind, argv, telescope=None, night=None, targets='', cwd=None, env=None, cores=1,
                disk_heavy=True, est_minutes=60.0, locks=None, depends_on=None, dep_requires_success=False,
                priority=0.0, dedupe_key=None, max_retries=2, source='manual', meta=None, note=None,
                not_before=None):
        """Queue a job. Returns (job_id, created); an active job with the same dedupe_key is returned
        instead of a duplicate."""
        if cls not in CLASSES:
            raise ValueError('class must be one of {}'.format(', '.join(CLASSES)))
        if not argv:
            raise ValueError('a job needs a command')
        now = iso(self.now())
        with self.tx() as c:
            if dedupe_key:
                r = c.execute("SELECT id FROM jobs WHERE dedupe_key = ? AND state IN ('queued', 'running')",
                              (dedupe_key,)).fetchone()
                if r:
                    return r['id'], False
            cur = c.execute(
                'INSERT INTO jobs (class, kind, telescope, night, targets, argv, cwd, env, cores, disk_heavy, '
                'est_minutes, locks, depends_on, dep_requires_success, priority, dedupe_key, state, max_retries, '
                'created_at, source, meta, note, not_before) '
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'queued', ?, ?, ?, ?, ?, ?)",
                (cls, kind, telescope, night, targets or '', json.dumps(list(argv)), cwd, json.dumps(env or {}),
                 int(cores), int(bool(disk_heavy)), float(est_minutes), json.dumps(sorted(set(locks or []))),
                 depends_on, int(bool(dep_requires_success)), float(priority), dedupe_key, int(max_retries), now,
                 source, json.dumps(meta or {}, sort_keys=True, default=str), note,
                 iso(not_before) if isinstance(not_before, dt.datetime) else not_before))
            jid = cur.lastrowid
            self._event(c, 'added', job_id=jid, source=source, telescope=telescope, night=night,
                        message='{} {} {} {} cores={} est={:.0f}m'.format(cls, kind, telescope or '-', night or '-',
                                                                       cores, est_minutes))
        return jid, True

    def get_job(self, job_id):
        return _job(self.conn.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone())

    def jobs(self, states=None, classes=None, telescope=None, night=None, kind=None, limit=None, newest_first=False):
        q, args = 'SELECT * FROM jobs WHERE 1=1', []
        for col, vals in (('state', states), ('class', classes)):
            if vals:
                q += ' AND {} IN ({})'.format(col, ','.join('?' * len(vals)))
                args.extend(vals)
        for col, val in (('telescope', telescope), ('night', night), ('kind', kind)):
            if val is not None:
                q += ' AND {} = ?'.format(col)
                args.append(val)
        q += ' ORDER BY id DESC' if newest_first else ' ORDER BY id'
        if limit:
            q += ' LIMIT {:d}'.format(limit)
        return [_job(r) for r in self.conn.execute(q, args)]

    def active_for(self, telescope, night):
        return self.jobs(states=ACTIVE_STATES, telescope=telescope, night=night)

    def states_of(self, ids):
        ids = [i for i in ids if i is not None]
        if not ids:
            return {}
        rows = self.conn.execute('SELECT id, state FROM jobs WHERE id IN ({})'.format(','.join('?' * len(ids))), ids)
        return {r['id']: r['state'] for r in rows}

    def update_job(self, job_id, **fields):
        if not fields:
            return
        for k in JOB_JSON:
            if k in fields:
                fields[k] = json.dumps(fields[k], sort_keys=True, default=str)
        cols = ', '.join('{} = ?'.format(k) for k in fields)
        self.conn.execute('UPDATE jobs SET {} WHERE id = ?'.format(cols), list(fields.values()) + [job_id])

    def mark_running(self, job_id, log_path=None, status_path=None, code_sha=None, est_minutes=None,
                     shadow=False, source='dispatcher', message=None):
        """queued -> running. Returns False if the job was no longer queued."""
        with self.tx() as c:
            job = _job(c.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone())
            if job is None or job['state'] != 'queued':
                return False
            c.execute("UPDATE jobs SET state = 'running', attempts = attempts + 1, started_at = ?, "
                      'finished_at = NULL, exit_code = NULL, pid = NULL, proc_start = NULL, kill_sent_at = NULL, '
                      'kill_reason = NULL, log_path = ?, status_path = ?, code_sha = ?, est_minutes = ?, shadow = ? '
                      'WHERE id = ?',
                      (iso(self.now()), log_path, status_path, code_sha,
                       float(est_minutes if est_minutes is not None else job['est_minutes']), int(shadow), job_id))
            self._event(c, 'shadow-start' if shadow else 'started', job_id=job_id, source=source,
                        telescope=job['telescope'], night=job['night'],
                        message=message or 'attempt {}'.format(job['attempts'] + 1))
        return True

    def set_process(self, job_id, pid, proc_start):
        self.conn.execute('UPDATE jobs SET pid = ?, proc_start = ? WHERE id = ?', (pid, proc_start, job_id))

    def mark_kill_sent(self, job_id, reason, signame, source='dispatcher'):
        with self.tx() as c:
            c.execute('UPDATE jobs SET kill_sent_at = COALESCE(kill_sent_at, ?), kill_reason = ? WHERE id = ?',
                      (iso(self.now()), reason, job_id))
            self._event(c, 'kill', job_id=job_id, source=source, message='{} sent ({})'.format(signame, reason))

    def finish(self, job_id, exit_code, failure_kind=None, signature=None, transient=False,
               retry_delay_minutes=60, source='dispatcher', message=None, shadow=False):
        """running -> done, failed, cancelled, or back to queued for a transient failure with retries left."""
        with self.tx() as c:
            job = _job(c.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone())
            if job is None or job['state'] != 'running':
                return job
            now = self.now()
            not_before = None
            if exit_code == 0 and not failure_kind:
                state = 'done'
            elif failure_kind == 'cancelled' or job['cancel_requested']:
                state, failure_kind = 'cancelled', 'cancelled'
            else:
                retries_used = job['attempts'] - 1 - job['vanished']
                if transient and retries_used < job['max_retries']:
                    state = 'queued'
                    not_before = iso(now + dt.timedelta(minutes=retry_delay_minutes))
                else:
                    state = 'failed'
                failure_kind = failure_kind or ('transient' if transient else 'deterministic')
            c.execute('UPDATE jobs SET state = ?, finished_at = ?, exit_code = ?, failure_kind = ?, '
                      'failure_signature = ?, not_before = ?, pid = NULL, proc_start = NULL WHERE id = ?',
                      (state, iso(now), exit_code, failure_kind if state != 'done' else None,
                       signature if state != 'done' else None, not_before, job_id))
            kind = {'done': 'finished', 'queued': 'retry', 'failed': 'failed', 'cancelled': 'cancelled'}[state]
            if shadow:
                kind = 'shadow-' + kind
            self._event(c, kind, job_id=job_id, source=source, telescope=job['telescope'], night=job['night'],
                        message=message or 'exit {} {}'.format(exit_code, signature or ''),
                        data={'exit_code': exit_code, 'failure_kind': failure_kind, 'signature': signature,
                              'not_before': not_before})
            if state in ('failed', 'cancelled'):
                self._cascade(c, job_id, state, source)
        return self.get_job(job_id)

    def vanish(self, job_id, max_requeues=1, source='dispatcher'):
        """A running job whose process is gone and left no status file (dispatcher or container restart)."""
        with self.tx() as c:
            job = _job(c.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone())
            if job is None or job['state'] != 'running':
                return job
            if job['vanished'] < max_requeues and not job['cancel_requested']:
                c.execute("UPDATE jobs SET state = 'queued', vanished = vanished + 1, not_before = NULL, pid = NULL, "
                          'proc_start = NULL, kill_sent_at = NULL, kill_reason = NULL WHERE id = ?', (job_id,))
                self._event(c, 'vanished-requeued', job_id=job_id, source=source, telescope=job['telescope'],
                            night=job['night'], message='process gone with no exit status; requeued once')
            else:
                sig = 'process vanished (dispatcher or container restart)'
                c.execute("UPDATE jobs SET state = 'failed', finished_at = ?, failure_kind = 'vanished', "
                          'failure_signature = ?, pid = NULL, proc_start = NULL WHERE id = ?',
                          (iso(self.now()), sig, job_id))
                self._event(c, 'failed', job_id=job_id, source=source, telescope=job['telescope'],
                            night=job['night'], message=sig)
                self._cascade(c, job_id, 'failed', source)
        return self.get_job(job_id)

    def _cascade(self, c, job_id, state, source):
        """Dependents of a cancelled job are cancelled; dependents that need success are cancelled when it fails."""
        rows = c.execute("SELECT id, dep_requires_success FROM jobs WHERE depends_on = ? AND state = 'queued'",
                         (job_id,)).fetchall()
        for r in rows:
            if state == 'cancelled' or r['dep_requires_success']:
                c.execute("UPDATE jobs SET state = 'cancelled', finished_at = ?, failure_kind = 'cancelled', "
                          'failure_signature = ? WHERE id = ?',
                          (iso(self.now()), 'dependency {} {}'.format(job_id, state), r['id']))
                self._event(c, 'cancelled', job_id=r['id'], source=source,
                            message='dependency {} {}'.format(job_id, state))
                self._cascade(c, r['id'], 'cancelled', source)

    def cancel(self, job_id, kill=False, source='cli'):
        with self.tx() as c:
            job = _job(c.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone())
            if job is None:
                raise KeyError('no job {}'.format(job_id))
            if job['state'] == 'queued':
                c.execute("UPDATE jobs SET state = 'cancelled', finished_at = ?, failure_kind = 'cancelled' "
                          'WHERE id = ?', (iso(self.now()), job_id))
                self._event(c, 'cancelled', job_id=job_id, source=source, telescope=job['telescope'],
                            night=job['night'])
                self._cascade(c, job_id, 'cancelled', source)
                return 'cancelled'
            if job['state'] == 'running':
                if not kill:
                    raise ValueError('job {} is running; cancelling it kills it, so pass --kill'.format(job_id))
                c.execute('UPDATE jobs SET cancel_requested = 1 WHERE id = ?', (job_id,))
                self._event(c, 'cancel-requested', job_id=job_id, source=source)
                return 'kill requested'
            return 'already ' + job['state']

    def requeue(self, job_ids=None, signature=None, exact=False, dry_run=False, source='cli'):
        """Put failed jobs back in the queue with fresh retries. Returns [(id, outcome)]."""
        q, args = "SELECT * FROM jobs WHERE state = 'failed'", []
        if job_ids:
            q += ' AND id IN ({})'.format(','.join('?' * len(job_ids)))
            args.extend(job_ids)
        if signature is not None:
            if exact:
                q += ' AND failure_signature = ?'
                args.append(signature)
            else:
                q += ' AND failure_signature LIKE ?'
                args.append('%{}%'.format(signature))
        if not job_ids and signature is None:
            raise ValueError('requeue needs job ids or a signature')
        out = []
        with self.tx() as c:
            for job in [_job(r) for r in c.execute(q + ' ORDER BY id', args).fetchall()]:
                if job['dedupe_key']:
                    dup = c.execute("SELECT id FROM jobs WHERE dedupe_key = ? AND state IN ('queued', 'running')",
                                    (job['dedupe_key'],)).fetchone()
                    if dup:
                        out.append((job['id'], 'skipped: job {} already active'.format(dup['id'])))
                        continue
                if not dry_run:
                    c.execute("UPDATE jobs SET state = 'queued', attempts = 0, vanished = 0, not_before = NULL, "
                              'exit_code = NULL, failure_kind = NULL, failure_signature = NULL, finished_at = NULL, '
                              'cancel_requested = 0 WHERE id = ?', (job['id'],))
                    self._event(c, 'requeued', job_id=job['id'], source=source, telescope=job['telescope'],
                                night=job['night'], message='was: {}'.format(job['failure_signature']))
                out.append((job['id'], 'would requeue' if dry_run else 'requeued'))
        return out

    # ---------------------------------------------------------------- watcher and look-back state

    def watch_get(self, telescope, night):
        return _plain(self.conn.execute('SELECT * FROM watch WHERE telescope = ? AND night = ?',
                                        (telescope, night)).fetchone())

    def watch_put(self, row):
        row = dict(row)
        row['data'] = json.dumps(row.get('data') or {}, sort_keys=True, default=str)
        vals = [row.get(k) for k in WATCH_COLS]
        self.conn.execute('INSERT OR REPLACE INTO watch ({}) VALUES ({})'.format(
            ', '.join(WATCH_COLS), ', '.join('?' * len(WATCH_COLS))), vals)

    def watch_rows(self, since_night=None):
        q, args = 'SELECT * FROM watch', []
        if since_night:
            q += ' WHERE night >= ?'
            args.append(since_night)
        return [_plain(r) for r in self.conn.execute(q + ' ORDER BY night, telescope', args)]

    def lookback_get(self, telescope, night):
        return _plain(self.conn.execute('SELECT * FROM lookback WHERE telescope = ? AND night = ?',
                                        (telescope, night)).fetchone())

    def lookback_put(self, row):
        row = dict(row)
        row['data'] = json.dumps(row.get('data') or {}, sort_keys=True, default=str)
        vals = [row.get(k) for k in LOOKBACK_COLS]
        self.conn.execute('INSERT OR REPLACE INTO lookback ({}) VALUES ({})'.format(
            ', '.join(LOOKBACK_COLS), ', '.join('?' * len(LOOKBACK_COLS))), vals)

    def lookback_rows(self, since_night=None):
        q, args = 'SELECT * FROM lookback', []
        if since_night:
            q += ' WHERE night >= ?'
            args.append(since_night)
        return [_plain(r) for r in self.conn.execute(q + ' ORDER BY night, telescope', args)]


def job_is_due(job, now):
    nb = parse_iso(job.get('not_before'))
    return nb is None or nb <= now
