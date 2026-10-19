"""What became of the PWV correction in each SSO pipeline job, read from the job's log, for the notifier.

ZLP_pipeline.sh echoes one `python .../SPlightcurve.py ...` line per light curve (T10), and SPlightcurve.py then
prints the outcome of that light curve's PWV correction:

  ok      'PWV correction successful!'
  nodata  'No PWV data for: <window>': the LHATPRO fetch raised, twice. From 15:50 UTC on 2026-10-18 that was
          name resolution inside the container, and every light curve went out without PWV for a day.
  off     'PWV Failed: No T_eff for target': never corrected, by design (Artemis, SAINT-EX, Teff < 2000 K)
  other   any other 'PWV Failed: ...' (an empty LHATPRO window, a bug)
  none    the light curve never reached the PWV step

This is the classification of mh_scratch/pwv_outage/scan_pwv.py. A job is 'nodata' when most of its light curves
that tried a correction got no LHATPRO data, 'ok' when any were corrected otherwise, else 'neutral'.

Job logs reach 20 MB, so each finished job is read once and its counts kept in the queue database (job_pwv),
keyed by job and finish time: a requeued job that runs again is read again.
"""
import datetime as dt
import re

from .config import eso_telescopes
from .util import brief_error, iso

SCHEMA = """
CREATE TABLE IF NOT EXISTS job_pwv (
    job_id INTEGER PRIMARY KEY,
    telescope TEXT,
    night TEXT,
    finished_at TEXT NOT NULL,
    log_path TEXT,
    n_lc INTEGER NOT NULL,
    n_ok INTEGER NOT NULL,
    n_nodata INTEGER NOT NULL,
    n_off INTEGER NOT NULL,
    n_other INTEGER NOT NULL,
    n_none INTEGER NOT NULL,
    verdict TEXT NOT NULL,
    reason TEXT,
    scanned_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS job_pwv_finished ON job_pwv(finished_at);
"""
COLS = ('job_id', 'telescope', 'night', 'finished_at', 'log_path', 'n_lc', 'n_ok', 'n_nodata', 'n_off', 'n_other',
        'n_none', 'verdict', 'reason', 'scanned_at')

SPLC = re.compile(r'^python \S*SPlightcurve\.py ')
NODATA = 'No PWV data for:'
# since the LHATPRO fix the cause follows the window: 'No PWV data for: <start> - <end> (ConnectionError: ...)'
NODATA_CAUSE = re.compile(r'^No PWV data for: \S+ - \S+ \((.*)\)\s*$')


def scan_log(path):
    """{'n_lc', 'n_ok', 'n_nodata', 'n_off', 'n_other', 'n_none', 'reason'} for one log; reason is the first cause
    given for missing LHATPRO data, or with none missing the first other PWV failure. Raises OSError if the log
    cannot be read."""
    counts = dict(n_lc=0, n_ok=0, n_nodata=0, n_off=0, n_other=0, n_none=0)
    causes, others = [], []
    cur = None

    def close(lc):
        if lc is None:
            return
        counts['n_lc'] += 1
        outcome = 'ok' if lc['ok'] else 'nodata' if lc['nodata'] else 'off' if lc['off'] else \
            'other' if lc['other'] else 'none'
        counts['n_' + outcome] += 1

    with open(path, errors='replace') as f:
        for line in f:
            if line.startswith('python ') and SPLC.match(line):
                close(cur)
                cur = {'ok': False, 'nodata': False, 'off': False, 'other': False}
            elif cur is None or not line.startswith(('PWV ', NODATA)):
                continue
            elif line.startswith('PWV correction successful'):
                cur['ok'] = True
            elif line.startswith(NODATA):
                cur['nodata'] = True
                m = NODATA_CAUSE.match(line)
                if m and not causes:
                    causes.append(m.group(1))
            elif line.startswith('PWV Failed: No T_eff'):
                cur['off'] = True
            elif line.startswith('PWV Failed:'):
                cur['other'] = True
                if not others:
                    others.append(line[len('PWV Failed:'):].strip())
    close(cur)
    # a job without LHATPRO data that logged no cause gets none: what follows ('All-NaN slice encountered') is not one
    counts['reason'] = brief_error(causes[0], 300) if causes else \
        others[0][:300] if others and not counts['n_nodata'] else None
    return counts


def verdict(c):
    tried = c['n_ok'] + c['n_nodata'] + c['n_other']
    if c['n_nodata'] and 2 * c['n_nodata'] > tried:
        return 'nodata'
    return 'ok' if c['n_ok'] else 'neutral'


class PwvLog:
    """The job_pwv table, brought up to date with SSO pipeline jobs that finished in the last scan_hours.

    record=False (the notifier's dry runs) reads the new logs but keeps their counts in memory only."""

    def __init__(self, cfg, store, now, record=True):
        self.cfg, self.store, self.now, self.record = cfg, store, now, record
        self.n = cfg['notify']
        self.fresh = {}                               # job_id -> row scanned this run and not saved
        store.conn.executescript(SCHEMA)

    def refresh(self):
        """Read the logs of finished SSO pipeline jobs not yet read (newest first, at most pwv_max_scans)."""
        since = iso(self.now - dt.timedelta(hours=float(self.n['pwv_scan_hours'])))
        tels = eso_telescopes(self.cfg)
        if not tels:
            return 0
        known = {r['job_id']: r['finished_at'] for r in self.store.conn.execute(
            'SELECT job_id, finished_at FROM job_pwv WHERE finished_at >= ?', (since,))}
        todo = [r for r in self.store.conn.execute(
            "SELECT id, telescope, night, finished_at, log_path FROM jobs WHERE kind = 'pipeline' AND "
            "state IN ('done', 'failed') AND shadow = 0 AND log_path IS NOT NULL AND finished_at >= ? AND "
            'telescope IN ({}) ORDER BY finished_at DESC, id DESC'.format(','.join('?' * len(tels))),
            [since] + tels) if known.get(r['id']) != r['finished_at']]
        scanned = 0
        for r in todo[:int(self.n['pwv_max_scans'])]:
            try:
                c = scan_log(r['log_path'])
            except OSError as e:                     # kept as a job with no light curves, so it is not retried
                c = dict(n_lc=0, n_ok=0, n_nodata=0, n_off=0, n_other=0, n_none=0,
                         reason='cannot read the log: {}'.format(e.strerror or e))
            row = dict(c, job_id=r['id'], telescope=r['telescope'], night=r['night'], finished_at=r['finished_at'],
                       log_path=r['log_path'], verdict=verdict(c), scanned_at=iso(self.now))
            if self.record:
                self.store.conn.execute('INSERT OR REPLACE INTO job_pwv ({}) VALUES ({})'.format(
                    ', '.join(COLS), ', '.join('?' * len(COLS))), [row[k] for k in COLS])
            else:
                self.fresh[r['id']] = row
            scanned += 1
        return scanned

    def _rows(self, where, args):
        return {r['job_id']: dict(r) for r in self.store.conn.execute(
            'SELECT * FROM job_pwv WHERE {} ORDER BY finished_at, job_id'.format(where), args)}

    def finished_since(self, since):
        """Rows of the jobs that finished at or after `since`, oldest first."""
        rows = self._rows('finished_at >= ?', (iso(since),))
        rows.update({k: v for k, v in self.fresh.items() if v['finished_at'] >= iso(since)})
        return sorted(rows.values(), key=lambda r: (r['finished_at'], r['job_id']))

    def last_ok(self):
        """The newest job that corrected PWV, or None."""
        r = self.store.conn.execute("SELECT * FROM job_pwv WHERE verdict = 'ok' ORDER BY finished_at DESC, job_id DESC "
                                    'LIMIT 1').fetchone()
        best = dict(r) if r else None
        for row in self.fresh.values():
            if row['verdict'] == 'ok' and (best is None or (row['finished_at'], row['job_id']) >
                                           (best['finished_at'], best['job_id'])):
                best = row
        return best

    def failing_run(self):
        """(jobs without LHATPRO data since the last job that corrected PWV, oldest first; that job or None)."""
        ok = self.last_ok()
        after = ok['finished_at'] if ok else ''
        rows = self._rows("verdict = 'nodata' AND finished_at > ?", (after,))
        for k, v in self.fresh.items():
            rows.pop(k, None)
            if v['verdict'] == 'nodata' and v['finished_at'] > after:
                rows[k] = v
        return sorted(rows.values(), key=lambda r: (r['finished_at'], r['job_id'])), ok


def summary(rows):
    """Light-curve totals over job rows, and how many jobs had no LHATPRO data."""
    out = {k: sum(int(r[k]) for r in rows) for k in ('n_lc', 'n_ok', 'n_nodata', 'n_off', 'n_other', 'n_none')}
    out['jobs'] = len(rows)
    out['jobs_nodata'] = sum(1 for r in rows if r['verdict'] == 'nodata')
    return out

