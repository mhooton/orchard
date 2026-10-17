"""
promote_journal.py
==================

The promotion journal: one SQLite file (by default <queue_root>/promote.sqlite) that records every
verdict promote.py reached with --apply, every promotion it carried out, and each step of each
promotion before and after it ran. A crash in the middle of a promotion therefore leaves a row in
state 'in_progress' with its plan and the steps already done, and `promote.py resume` finishes it.

Standard library only, so the queue's daily report can read it without numpy or astropy. Like the
queue database it must live on appct's local disk and is only opened inside orchard-server: SQLite
locking cannot be trusted over NFS, and appcs sees that disk over NFS.

Tables
    invocations  one row per promote.py command that wrote to the journal (run, accept, reject, resume)
    verdicts     one row per judged candidate: a v3 target directory with a given fingerprint
    promotions   one row per promotion (automatic or by accept), with every path decided up front
    steps        the intent and the outcome of each step of a promotion

Verdict states
    pending    HOLD or REJECT, waiting for a person: `promote.py accept` or `promote.py reject`
    promoting  a promotion of it started and has not finished, or failed (see promotions)
    promoted   promoted automatically (verdict PROMOTE)
    accepted   promoted by `accept`
    rejected   closed by `reject`; nothing was moved
    outdated   the candidate changed or left v3 before it was promoted or decided; a newer verdict, if any,
               replaces it

Promotion states
    in_progress  started; `promote.py resume` finishes it after a crash
    done         every step done or skipped
    aborted      the pre-flight check failed, so nothing was moved
    failed       a step failed after something was moved; it needs a person (`resume --promotion`, `abandon`)
    abandoned    closed by `abandon` after being sorted out by hand
"""
import contextlib
import datetime as dt
import json
import os
import sqlite3

ISO = '%Y-%m-%dT%H:%M:%SZ'

SCHEMA = """
CREATE TABLE IF NOT EXISTS invocations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    started_at TEXT NOT NULL,
    finished_at TEXT,
    command TEXT NOT NULL,
    argv TEXT,
    host TEXT,
    pid INTEGER,
    code TEXT,
    job_id INTEGER,
    pipeline_job INTEGER,
    summary TEXT
);

CREATE TABLE IF NOT EXISTS verdicts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    invocation_id INTEGER,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    telescope TEXT NOT NULL,
    night TEXT NOT NULL,
    target TEXT NOT NULL,
    candidate TEXT NOT NULL,
    incumbent TEXT,
    verdict TEXT NOT NULL,
    reasons TEXT NOT NULL DEFAULT '',
    notes TEXT NOT NULL DEFAULT '',
    metrics TEXT NOT NULL DEFAULT '{}',
    run TEXT NOT NULL DEFAULT '{}',
    cand_fp TEXT,
    inc_fp TEXT,
    state TEXT NOT NULL,
    decided_by TEXT,
    decided_at TEXT,
    decision_reason TEXT,
    promotion_id INTEGER
);
CREATE INDEX IF NOT EXISTS verdicts_target ON verdicts(telescope, night, target);
CREATE INDEX IF NOT EXISTS verdicts_state ON verdicts(state);

CREATE TABLE IF NOT EXISTS promotions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    verdict_id INTEGER NOT NULL,
    invocation_id INTEGER,
    telescope TEXT NOT NULL,
    night TEXT NOT NULL,
    target TEXT NOT NULL,
    how TEXT NOT NULL,
    started_at TEXT NOT NULL,
    finished_at TEXT,
    state TEXT NOT NULL,
    plan TEXT NOT NULL,
    error TEXT
);
CREATE INDEX IF NOT EXISTS promotions_state ON promotions(state);

CREATE TABLE IF NOT EXISTS steps (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    promotion_id INTEGER NOT NULL,
    step TEXT NOT NULL,
    state TEXT NOT NULL,
    intent_at TEXT NOT NULL,
    done_at TEXT,
    detail TEXT NOT NULL DEFAULT '{}'
);
CREATE INDEX IF NOT EXISTS steps_promotion ON steps(promotion_id);
"""

OPEN_STATES = ('pending', 'promoting')
PROMOTED_STATES = ('promoted', 'accepted')


def utcnow():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


def iso(t):
    return t.strftime(ISO) if t else None


def parse_iso(s):
    if not s:
        return None
    return dt.datetime.strptime(s.rstrip('Z').split('.')[0], '%Y-%m-%dT%H:%M:%S')


def _dumps(obj):
    return json.dumps(obj, sort_keys=True, default=str)


def _row(r, json_cols=('metrics', 'run', 'plan', 'detail', 'summary')):
    if r is None:
        return None
    d = dict(r)
    for k in json_cols:
        if k in d:
            d[k] = json.loads(d[k]) if d[k] else {}
    return d


class Journal(object):
    """The promotion journal. readonly=True opens an existing file without creating or changing it."""

    def __init__(self, path, now_fn=utcnow, readonly=False, check_local=True):
        self.path = path
        self.now = now_fn
        self.readonly = readonly
        if readonly:
            self.conn = sqlite3.connect('file:{}?mode=ro'.format(path), uri=True, timeout=60, isolation_level=None)
        else:
            directory = os.path.dirname(os.path.abspath(path))
            os.makedirs(directory, exist_ok=True)
            if check_local:
                from jobqueue.util import ensure_local_fs
                ensure_local_fs(directory)
            self.conn = sqlite3.connect(path, timeout=60, isolation_level=None)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute('PRAGMA busy_timeout=60000')
        if not readonly:
            self.conn.execute('PRAGMA journal_mode=WAL')
            self.conn.execute('PRAGMA synchronous=FULL')
            self.conn.executescript(SCHEMA)

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

    # ------------------------------------------------------------------ invocations

    def invocation_start(self, command, argv=(), code=None, job_id=None, pipeline_job=None):
        import socket
        with self.tx() as c:
            cur = c.execute('INSERT INTO invocations (started_at, command, argv, host, pid, code, job_id, '
                            'pipeline_job) VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
                            (iso(self.now()), command, _dumps(list(argv)), socket.gethostname(), os.getpid(), code,
                             job_id, pipeline_job))
            return cur.lastrowid

    def invocation_end(self, inv_id, summary):
        with self.tx() as c:
            c.execute('UPDATE invocations SET finished_at = ?, summary = ? WHERE id = ?',
                      (iso(self.now()), _dumps(summary), inv_id))

    # ------------------------------------------------------------------ verdicts

    def verdict(self, vid):
        return _row(self.conn.execute('SELECT * FROM verdicts WHERE id = ?', (vid,)).fetchone())

    def latest_verdict(self, telescope, night, target):
        return _row(self.conn.execute('SELECT * FROM verdicts WHERE telescope = ? AND night = ? AND target = ? '
                                      'ORDER BY id DESC LIMIT 1', (telescope, night, target)).fetchone())

    def verdicts(self, states=None, telescope=None, night=None, since=None, limit=None):
        q, args = 'SELECT * FROM verdicts WHERE 1=1', []
        if states:
            q += ' AND state IN ({})'.format(','.join('?' * len(states)))
            args.extend(states)
        for col, val in (('telescope', telescope), ('night', night)):
            if val:
                q += ' AND {} = ?'.format(col)
                args.append(val)
        if since:
            q += ' AND updated_at >= ?'
            args.append(since)
        q += ' ORDER BY id'
        if limit:
            q += ' LIMIT {:d}'.format(limit)
        return [_row(r) for r in self.conn.execute(q, args)]

    def add_verdict(self, invocation_id, telescope, night, target, candidate, incumbent, verdict, reasons, notes,
                    metrics, run, cand_fp, inc_fp, state):
        now = iso(self.now())
        with self.tx() as c:
            cur = c.execute(
                'INSERT INTO verdicts (invocation_id, created_at, updated_at, telescope, night, target, candidate, '
                'incumbent, verdict, reasons, notes, metrics, run, cand_fp, inc_fp, state) '
                'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                (invocation_id, now, now, telescope, night, target, candidate, incumbent, verdict, reasons, notes,
                 _dumps(metrics), _dumps(run), cand_fp, inc_fp, state))
            return cur.lastrowid

    def touch_verdict(self, vid):
        self.conn.execute('UPDATE verdicts SET updated_at = ? WHERE id = ?', (iso(self.now()), vid))

    def set_verdict_state(self, vid, state, expect=None, **fields):
        """Change a verdict's state (and other columns); with expect, only from one of those states.
        Returns True if the row changed."""
        fields = dict(fields, state=state, updated_at=iso(self.now()))
        cols = ', '.join('{} = ?'.format(k) for k in fields)
        q, args = 'UPDATE verdicts SET {} WHERE id = ?'.format(cols), list(fields.values()) + [vid]
        if expect:
            q += ' AND state IN ({})'.format(','.join('?' * len(expect)))
            args.extend(expect)
        with self.tx() as c:
            return c.execute(q, args).rowcount == 1

    # ------------------------------------------------------------------ promotions and steps

    def start_promotion(self, verdict_id, invocation_id, telescope, night, target, how, make_plan):
        """Claim an open verdict for promotion. make_plan(promotion_id) returns the plan, every path decided up
        front; it is stored in the same transaction, so a promotion never exists without one. Returns the
        promotion id, or None if the verdict was no longer open."""
        now = iso(self.now())
        with self.tx() as c:
            row = c.execute('SELECT state FROM verdicts WHERE id = ?', (verdict_id,)).fetchone()
            if row is None or row['state'] not in ('pending', 'promoting'):
                return None
            busy = c.execute("SELECT id FROM promotions WHERE verdict_id = ? AND state IN ('in_progress', 'failed')",
                             (verdict_id,)).fetchone()
            if busy:
                return None
            cur = c.execute('INSERT INTO promotions (verdict_id, invocation_id, telescope, night, target, how, '
                            "started_at, state, plan) VALUES (?, ?, ?, ?, ?, ?, ?, 'in_progress', '{}')",
                            (verdict_id, invocation_id, telescope, night, target, how, now))
            pid = cur.lastrowid
            c.execute('UPDATE promotions SET plan = ? WHERE id = ?', (_dumps(make_plan(pid)), pid))
            c.execute("UPDATE verdicts SET state = 'promoting', promotion_id = ?, updated_at = ? WHERE id = ?",
                      (pid, now, verdict_id))
            return pid

    def promotion(self, pid):
        return _row(self.conn.execute('SELECT * FROM promotions WHERE id = ?', (pid,)).fetchone())

    def promotions(self, states=None, since=None, telescope=None, night=None):
        q, args = 'SELECT * FROM promotions WHERE 1=1', []
        if states:
            q += ' AND state IN ({})'.format(','.join('?' * len(states)))
            args.extend(states)
        if since:
            q += ' AND started_at >= ?'
            args.append(since)
        for col, val in (('telescope', telescope), ('night', night)):
            if val:
                q += ' AND {} = ?'.format(col)
                args.append(val)
        return [_row(r) for r in self.conn.execute(q + ' ORDER BY id', args)]

    def update_plan(self, pid, plan):
        with self.tx() as c:
            c.execute('UPDATE promotions SET plan = ? WHERE id = ?', (_dumps(plan), pid))

    def finish_promotion(self, pid, state, verdict_state=None, error=None):
        """Close a promotion; verdict_state moves its verdict too (back to pending clears an accept's decision)."""
        now = iso(self.now())
        with self.tx() as c:
            c.execute('UPDATE promotions SET state = ?, finished_at = ?, error = ? WHERE id = ?',
                      (state, now, error, pid))
            if verdict_state == 'pending':
                c.execute("UPDATE verdicts SET state = 'pending', promotion_id = NULL, decided_by = NULL, "
                          'decided_at = NULL, decision_reason = NULL, updated_at = ? WHERE promotion_id = ?',
                          (now, pid))
            elif verdict_state:
                c.execute('UPDATE verdicts SET state = ?, updated_at = ? WHERE promotion_id = ?',
                          (verdict_state, now, pid))

    def step_intent(self, pid, step, detail):
        with self.tx() as c:
            cur = c.execute("INSERT INTO steps (promotion_id, step, state, intent_at, detail) VALUES "
                            "(?, ?, 'intent', ?, ?)", (pid, step, iso(self.now()), _dumps(detail)))
            return cur.lastrowid

    def step_end(self, sid, state, detail):
        with self.tx() as c:
            c.execute('UPDATE steps SET state = ?, done_at = ?, detail = ? WHERE id = ?',
                      (state, iso(self.now()), _dumps(detail), sid))

    def steps(self, pid):
        return [_row(r) for r in self.conn.execute('SELECT * FROM steps WHERE promotion_id = ? ORDER BY id', (pid,))]


def review_summary(path, days=3, now=None):
    """What the daily report says about promotion: verdicts waiting for a person, promotions that did not
    finish, and what was promoted or decided in the last `days` days. None if there is no journal yet."""
    if not path or not os.path.exists(path):
        return None
    now = now or utcnow()
    since = iso(now - dt.timedelta(days=days))
    j = Journal(path, readonly=True)
    try:
        out = {
            'pending': j.verdicts(states=['pending']),
            'stuck': j.promotions(states=['in_progress', 'failed']),
            'promoted': j.verdicts(states=list(PROMOTED_STATES), since=since),
            'rejected': j.verdicts(states=['rejected'], since=since),
        }
    finally:
        j.close()
    return out
