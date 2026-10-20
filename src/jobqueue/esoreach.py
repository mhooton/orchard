"""Whether ESO answers: the outcome of each watcher poll and look-back run, kept for the notifier.

Every watcher poll and look-back run that asks ESO anything records one outcome in the queue database (kv
'eso_reach'): answered if any of its queries came back, failed if every one raised. A run of failures is counted
in polls and in watcher time: the minutes between consecutive polls, at most MAX_GAP_MINUTES each, so the watcher's
overnight break (it polls 11:00-23:30 local) is not counted as hours of outage. digest.eso_conditions opens an
incident from it, and the first answered query clears it.

On 2026-10-18 appct's first nameserver began answering SERVFAIL, no name resolved inside the containers, and the
watcher logged a failed query every 30 minutes for a day and a half without anyone being told.
"""
import json

from .util import brief_error, iso, parse_iso

KEY = 'eso_reach'
MAX_GAP_MINUTES = 60.0


def get(store):
    """The latest state: {} before anything was recorded, else state ('ok' or 'failing'), last_ok, last_failure,
    failing_since, failed_polls, failing_minutes, last_poll, last_source and error (brief, credentials scrubbed)."""
    try:
        return json.loads(store.kv_get(KEY) or '{}')
    except ValueError:
        return {}


def record(store, now, ok, error=None, source='watcher'):
    """Record one poll's outcome; returns the new state. One transaction, so the watcher and a look-back that
    overlap do not lose each other's outcome."""
    with store.tx() as c:
        row = c.execute('SELECT value FROM kv WHERE key = ?', (KEY,)).fetchone()
        try:
            s = json.loads(row['value']) if row and row['value'] else {}
        except ValueError:
            s = {}
        prev = parse_iso(s.get('last_poll'))
        gap = min((now - prev).total_seconds() / 60.0, MAX_GAP_MINUTES) if prev and now > prev else 0.0
        was_failing = s.get('state') == 'failing'
        since, polls = s.get('failing_since'), s.get('failed_polls')
        if ok:
            s.update(state='ok', last_ok=iso(now), failing_since=None, failed_polls=0, failing_minutes=0.0)
        else:
            if not was_failing:
                s.update(failing_since=iso(now), failed_polls=0, failing_minutes=0.0)
            s.update(state='failing', last_failure=iso(now), failed_polls=int(s.get('failed_polls') or 0) + 1,
                     failing_minutes=float(s.get('failing_minutes') or 0.0) + gap, error=brief_error(error, 300))
        s.update(last_poll=iso(max(now, prev) if prev else now), last_source=source)
        c.execute('INSERT OR REPLACE INTO kv (key, value, updated_at) VALUES (?, ?, ?)',
                  (KEY, json.dumps(s, sort_keys=True), iso(store.now())))
    if ok and was_failing:
        store.event('eso-recovered', source=source, message='ESO answered the {} at {} after {} failed poll{} since '
                    '{}'.format(source, iso(now), polls, '' if polls == 1 else 's', since))
    return s
