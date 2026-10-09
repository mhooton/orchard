"""Daily report. Live mode: per telescope-night, when the night was ready at ESO, when its download and
pipeline jobs ran, how they ended, whether light curves reached v2, and the delay from ready to done; then failed
jobs, look-back actions and watcher notes. Shadow mode: what the queue would have done, next to what the cron
actually did.

Shadow side (the shadow database): when the watcher found each night ready, when the shadow dispatcher
would have started its download and pipeline jobs and when they would have finished (simulated from
their estimates), notes, waiting reasons, and the look-back's actions.

Actual side (read-only, from the production tree): the cron's download log ESO_logs/<NIGHT><TEL>.log
(start, end, attempts), Observations/<TEL>/download_log.csv (counts), and the pipeline log
PipelineOutput/v3|v2/<TEL>/logs/<NIGHT>_1_v3.log (first and last `date` lines, in UTC inside the container).
"""
import csv
import datetime as dt
import os
import re
import statistics

from .nights import products_state
from .util import iso, last_night, night_str, parse_iso, utc_from_ts, utcnow

PY_TS = re.compile(r'^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)')
DATE_LINE = re.compile(r'^(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) (\w{3}) +(\d{1,2}) (\d\d:\d\d:\d\d) (\w+) (\d{4})\s*$')


def _mtime(path):
    return utc_from_ts(os.stat(path).st_mtime)


def actual_download(basedir, tel, night):
    path = os.path.join(basedir, 'ESO_logs', '{}{}.log'.format(night, tel))
    if not os.path.exists(path):
        return None
    start = None
    attempts = 0
    with open(path, errors='replace') as f:
        for line in f:
            if start is None:
                m = PY_TS.match(line)
                if m:
                    start = dt.datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S')
            if 'INITIAL ATTEMPT' in line or 'RETRY ATTEMPT' in line:
                attempts += 1
    return {'start': start, 'end': _mtime(path), 'attempts': attempts}


def actual_pipeline(basedir, tel, night):
    for version in ('v3', 'v2'):
        path = os.path.join(basedir, 'PipelineOutput', version, tel, 'logs', '{}_1_v3.log'.format(night))
        if not os.path.exists(path):
            continue
        stamps, complete = [], False
        with open(path, errors='replace') as f:
            for line in f:
                m = DATE_LINE.match(line.strip())
                if m:
                    try:
                        stamps.append(dt.datetime.strptime(' '.join(m.group(1, 2, 3, 5)), '%b %d %H:%M:%S %Y'))
                    except ValueError:
                        pass
                if 'PIPELINE COMPLETE' in line:
                    complete = True
        return {'start': stamps[0] if stamps else None, 'end': stamps[-1] if len(stamps) > 1 else _mtime(path),
                'complete': complete, 'where': version}
    return None


def download_counts(basedir, tel, night):
    path = os.path.join(basedir, 'Observations', tel, 'download_log.csv')
    try:
        with open(path) as f:
            for row in csv.DictReader(f):
                if row.get('Night') == night:
                    return row
    except OSError:
        pass
    return None


def _hm(t):
    return t.strftime('%d %H:%M') if t else '-'


def _hours(a, b):
    return (a - b).total_seconds() / 3600.0 if a and b else None


def _span(job):
    if not job:
        return '-'
    if job['state'] == 'queued':
        return 'queued' + (' (retry)' if job['attempts'] else '')
    start, end = parse_iso(job['started_at']), parse_iso(job['finished_at'])
    if job['state'] == 'running':
        return 'running since {}'.format(_hm(start))
    return '{}–{}{}'.format(_hm(start), _hm(end), '' if job['state'] == 'done' else ' ' + job['state'])


def _live_nights(cfg, store, nights, out):
    basedir = cfg['basedir']
    delays = []
    for night in nights:
        out.append('## Night {}'.format(night))
        out.append('')
        out.append('| Telescope | Ready at ESO | Download | Pipeline | Result | Light curves | Ready to done | Notes |')
        out.append('|---|---|---|---|---|---|---|---|')
        for tel in cfg['telescopes']:
            w = store.watch_get(tel, night) or {}
            pl = store.get_job(w['pipeline_job']) if w.get('pipeline_job') else None
            dl = store.get_job(w['download_job']) if w.get('download_job') else None
            ready = parse_iso(w.get('ready_at'))
            if w.get('state') == 'enqueued':
                ready_txt = _hm(ready)
            elif w.get('state') == 'processed':
                ready_txt = 'run outside the queue'
            elif w:
                ready_txt = '{} (ESO {}{})'.format(w['state'], w.get('eso_rows'),
                                                   ' of {}'.format(w['transfer_count']) if w.get('transfer_count') else '')
            else:
                ready_txt = 'not seen'
            result = '-'
            if pl and pl['state'] in ('done', 'failed', 'cancelled'):
                result = 'exit {}'.format(pl['exit_code']) if pl['state'] == 'done' else '{}: {}'.format(
                    pl['state'], pl['failure_signature'] or '')
            lcs = '-'
            if pl and pl['state'] in ('done', 'failed'):
                lcs = {'v2': 'v2', 'v3_only': 'v3 only', 'none': 'none'}[products_state(basedir, tel, night)]
            delay = _hours(parse_iso(pl['finished_at']), ready) if pl and pl['state'] == 'done' and ready else None
            if delay is not None:
                delays.append(delay)
            notes = []
            if w.get('first_seen') and ready:
                notes.append('first seen at ESO {}'.format(_hm(parse_iso(w['first_seen']))))
            if w.get('note'):
                notes.append(w['note'])
            counts = download_counts(basedir, tel, night)
            if counts and dl:
                notes.append('download log: ESO {} / transfer {} / on disk {}'.format(
                    counts.get('ESO_Archive'), counts.get('Transferred'), counts.get('Downloaded')))
            if pl:
                notes.append('jobs {}{}'.format('{} + '.format(dl['id']) if dl else '', pl['id']))
            out.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
                tel, ready_txt, _span(dl) if dl else ('-' if not pl else 'none'), _span(pl), result, lcs,
                '{:.1f} h'.format(delay) if delay is not None else '-', '; '.join(notes)))
        out.append('')
    return delays


def build_report(cfg, store, days=3, now=None):
    now = now or utcnow()
    basedir = cfg['basedir']
    ln = last_night(now)
    nights = [night_str(ln - dt.timedelta(days=k)) for k in range(days)]
    since = iso(now - dt.timedelta(days=days))
    if cfg['mode'] == 'live':
        out = ['# Job queue daily report, {}'.format(iso(now)), '',
               'Times are UTC (day hh:mm). "Ready to done" runs from the watcher finding the night complete at '
               'ESO (or Artemis on disk) to the end of its pipeline job.', '']
        delays = _live_nights(cfg, store, nights, out)
        failed = store.conn.execute("SELECT * FROM jobs WHERE state = 'failed' AND finished_at >= ? ORDER BY id",
                                    (since,)).fetchall()
        if failed:
            out.append('## Failed jobs')
            out.append('')
            for j in failed:
                out.append('- job {} {} {} {} {}: [{}] {} (log {})'.format(
                    j['id'], j['class'], j['kind'], j['telescope'] or '', j['night'] or '', j['failure_kind'],
                    j['failure_signature'] or '', j['log_path'] or '-'))
            out.append('')
        p1 = store.conn.execute("SELECT * FROM jobs WHERE class = 'P1' AND created_at >= ? ORDER BY id",
                                (since,)).fetchall()
        if p1:
            out.append('## P1 jobs (look-back)')
            out.append('')
            for j in p1:
                out.append('- job {} {} {} {}: {}{}'.format(
                    j['id'], j['kind'], j['telescope'] or '', j['night'] or '', j['state'],
                    ' ({})'.format(j['note']) if j['note'] else ''))
            out.append('')
        _trailer(store, since, out, waits_title='Jobs that had to wait')
        summary = '{} failed job(s) in the last {} days.'.format(len(failed), days)
        if delays:
            summary = ('Ready to done: median {:.1f} h over {} telescope-nights (range {:.1f} to {:.1f} h). '
                       .format(statistics.median(delays), len(delays), min(delays), max(delays)) + summary)
        out.insert(4, summary)
        out.insert(5, '')
        return '\n'.join(out) + '\n'

    out = ['# Job queue shadow report, {}'.format(iso(now)), '',
           'Times are UTC (day hh:mm). "Would" times come from the shadow dispatcher, which starts nothing and '
           'treats a job as finished after its estimate. "Cron" times come from the production logs.', '']
    gains = []
    for night in nights:
        out.append('## Night {}'.format(night))
        out.append('')
        out.append('| Telescope | Watcher | Would start pipeline | Would finish | Cron download | '
                   'Cron pipeline | Gain | Notes |')
        out.append('|---|---|---|---|---|---|---|---|')
        for tel in cfg['telescopes']:
            w = store.watch_get(tel, night) or {}
            pl = store.get_job(w['pipeline_job']) if w.get('pipeline_job') else None
            dl = store.get_job(w['download_job']) if w.get('download_job') else None
            ready = parse_iso(w.get('ready_at'))
            would_start = parse_iso(pl['started_at']) if pl else None
            would_end = parse_iso(pl['finished_at']) if pl and pl['state'] == 'done' else None
            if would_start and not would_end and pl:
                would_end = would_start + dt.timedelta(minutes=float(pl['est_minutes']))
            ad = actual_download(basedir, tel, night)
            ap = actual_pipeline(basedir, tel, night)
            counts = download_counts(basedir, tel, night)
            gain = _hours(ap['start'], would_start) if ap and would_start else None
            if gain is not None:
                gains.append(gain)
            if w.get('state') == 'enqueued':
                watcher = 'ready {} ({})'.format(_hm(ready), w.get('ready_reason') or '')
            elif w:
                watcher = '{} (ESO {}, transfer log {})'.format(w.get('state'), w.get('eso_rows'), w.get('transfer_count'))
            else:
                watcher = 'not seen'
            notes = []
            if w.get('note'):
                notes.append(w['note'])
            if dl and dl['state'] not in ('done', 'running'):
                notes.append('download job {}'.format(dl['state']))
            if counts:
                notes.append('cron counts ESO {} / transfer {} / downloaded {}'.format(
                    counts.get('ESO_Archive'), counts.get('Transferred'), counts.get('Downloaded')))
            if ad and ad['attempts'] > 1:
                notes.append('cron download made {} attempts over {:.1f} h'.format(
                    ad['attempts'], _hours(ad['end'], ad['start']) or 0))
            if ap:
                notes.append('cron pipeline {} (log in {})'.format('completed' if ap['complete'] else 'did not complete',
                                                                   ap['where']))
            out.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
                tel, watcher, _hm(would_start), _hm(would_end),
                '{}–{}'.format(_hm(ad['start']), _hm(ad['end'])) if ad else '-',
                '{}–{}'.format(_hm(ap['start']), _hm(ap['end'])) if ap else '-',
                '{:+.1f} h'.format(gain) if gain is not None else '-', '; '.join(notes) or ''))
        out.append('')

    _trailer(store, since, out, waits_title='Shadow jobs that had to wait')
    if gains:
        out.insert(4, 'Pipeline start, shadow versus cron: median {:+.1f} h earlier over {} telescope-nights '
                      '(range {:+.1f} to {:+.1f} h).'.format(statistics.median(gains), len(gains), min(gains),
                                                             max(gains)))
        out.insert(5, '')
    return '\n'.join(out) + '\n'


def _trailer(store, since, out, waits_title):
    waits = store.events(since=since, kinds=['waiting'])
    if waits:
        out.append('## ' + waits_title)
        out.append('')
        for e in waits[-40:]:
            out.append('- {} job {} {} {}: {}'.format(e['ts'], e['job_id'], e['telescope'] or '', e['night'] or '',
                                                     e['message']))
        out.append('')
    actions = [r for r in store.lookback_rows()
               if (r.get('checked_at') or '') >= since and r.get('last_action') not in (None, 'ok')]
    if actions:
        out.append('## Look-back actions')
        out.append('')
        for r in actions:
            out.append('- {} {}: {} (ESO {} rows, disk {} files, products {})'.format(
                r['telescope'], r['night'], r['last_action'], r.get('eso_rows'), r.get('disk_files'),
                r.get('products')))
        out.append('')
    notes = store.events(since=since, kinds=['note', 'eso-error'])
    if notes:
        out.append('## Watcher notes and errors')
        out.append('')
        for e in notes:
            out.append('- {} {}'.format(e['ts'], e['message']))
        out.append('')
