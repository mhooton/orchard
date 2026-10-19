"""What the notifier says: each telescope's night, the queue, the disks, and the conditions that open incidents.

Read-only. Everything comes from the queue database (the watcher's per-night rows, jobs and their failure
signatures, look-back rows) and the production tree (raw night directories, SPOCK's night blocks, light curves).
Nothing here sends anything; notify.py delivers.

The digest describes last_night(now). At the morning digest that is the night before last: its frames reached
ESO yesterday afternoon and went through the pipeline then.

A night with no survey frames is attributed, not explained:
  - SPOCK's Archive_night_blocks name Chilean-time blocks ch_...; a night planned mostly for Chile has no survey
    frames by design (the downloader only asks ESO for 60.A-9009), so it is never a gap.
  - The four SSO telescopes share the weather. One that is dark while two others observe is down itself; when
    none observes it is the weather (or the site), and nobody here can act on it.
"""
import datetime as dt
import glob
import hashlib
import json
import os
import re
from collections import Counter

from .config import eso_telescopes
from .nights import frame_time, list_frames, night_dir
from .store import job_is_due
from .util import hhmm, iso, last_night, night_str, parse_iso, parse_night

BLOCK_TOKENS = re.compile(r'"[^"]*"|\S+')
TIB = float(2 ** 40)


def fmt_night(night):
    d = parse_night(night)
    return '{} {}'.format(d.day, d.strftime('%b %Y'))


def fmt_time(ts):
    t = parse_iso(ts) if isinstance(ts, str) else ts
    return '{} {}'.format(t.day, t.strftime('%b %H:%M')) if t else '-'


def _and(names):
    names = list(names)
    return names[0] if len(names) == 1 else '{} and {}'.format(', '.join(names[:-1]), names[-1])


def planned_minutes(basedir, tel, night):
    """SPOCK's plan for the night as {'survey': minutes, 'chilean': minutes}, or None without a night-blocks file."""
    path = os.path.join(basedir, 'Observations', tel, 'schedule', 'Archive_night_blocks',
                        'night_blocks_{}_{}.txt'.format(tel, parse_night(night).isoformat()))
    out = {'survey': 0.0, 'chilean': 0.0}
    try:
        with open(path, errors='replace') as f:
            for line in f:
                tok = BLOCK_TOKENS.findall(line)
                if len(tok) < 4 or tok[0] == 'target':
                    continue
                try:
                    minutes = float(tok[3])
                except ValueError:
                    continue
                out['chilean' if tok[0].strip('"').lower().startswith('ch_') else 'survey'] += minutes
    except OSError:
        return None
    return out


def lc_targets(basedir, tel, night):
    """Targets with a differential light curve for the night, in v3 or (after T12) v2."""
    found = set()
    for version in ('v3', 'v2'):
        pattern = os.path.join(basedir, 'PipelineOutput', version, tel, 'output', night, '*', '*_diff.fits')
        found.update(os.path.basename(os.path.dirname(p)) for p in glob.glob(pattern))
    return found


def span_hours(names):
    """First to last exposure from SPECU<n>.<timestamp> names (SSO); None where names carry no time (Artemis)."""
    times = [t for t in (frame_time(n) for n in names) if t]
    return (max(times) - min(times)).total_seconds() / 3600.0 if len(times) > 1 else None


def disk_usage(path):
    st = os.statvfs(path)
    return os.stat(path).st_dev, st.f_bavail * st.f_frsize, st.f_blocks * st.f_frsize


class Context:
    """Cached reads of the queue and the production tree for one notifier run."""

    def __init__(self, cfg, store, now):
        self.cfg, self.store, self.now = cfg, store, now
        self.n = cfg['notify']
        self.basedir = cfg['basedir']
        self.sso = eso_telescopes(cfg)
        self._cache = {}

    def _get(self, key, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    def watch(self, tel, night):
        return self._get(('watch', tel, night), lambda: self.store.watch_get(tel, night))

    def frames(self, tel, night):
        return self._get(('frames', tel, night), lambda: list_frames(night_dir(self.basedir, tel, night)))

    def plan(self, tel, night):
        return self._get(('plan', tel, night), lambda: planned_minutes(self.basedir, tel, night))

    def lcs(self, tel, night):
        return self._get(('lcs', tel, night), lambda: lc_targets(self.basedir, tel, night))

    def chilean(self, tel, night):
        p = self.plan(tel, night) or {}
        total = p.get('survey', 0.0) + p.get('chilean', 0.0)
        return total > 0 and p['chilean'] / total >= float(self.n['chilean_min_fraction'])

    def observed(self, tel, night):
        """True if the telescope took survey frames that night, False if not, None while the watcher still
        waits for them at ESO. The watcher's ESO counts first; older nights from light curves and frames on
        disk (a few frames and no light curve are calibrations, not observing)."""
        w = self.watch(tel, night)
        if w:
            by_type = (w.get('data') or {}).get('by_type') or {}
            if (by_type.get('OBJECT', 0) if by_type else w.get('eso_rows')):
                return True
            if w.get('state') in ('waiting', 'arriving'):
                return None
            if w.get('state') == 'no_data':
                return False
        if self.lcs(tel, night):
            return True
        return len(self.frames(tel, night)) >= int(self.n['dark_min_frames'])

    def siblings_observed(self, tel, night):
        return [t for t in self.sso if t != tel and self.observed(t, night)]


# ---------------------------------------------------------------------------------------------- one night

def _frames_text(s):
    if not s['frames']:
        return ''
    text = '{:,} frames'.format(s['frames'])
    if s['span_h']:
        text += ', {:.1f} h'.format(s['span_h'])
    return text


def night_status(ctx, tel, night):
    """One telescope's night: what arrived, what ran, what came out, and the digest line describing it."""
    w = ctx.watch(tel, night) or {}
    jobs = ctx.store.jobs(telescope=tel, night=night)
    latest = {}
    for j in jobs:                                   # ascending id: the last of each kind wins
        latest[j['kind']] = j
    active = [j for j in jobs if j['state'] in ('queued', 'running')]
    names = ctx.frames(tel, night)
    lcs = ctx.lcs(tel, night)
    targets = list((w.get('data') or {}).get('targets') or [])
    s = {'telescope': tel, 'night': night, 'watch': w.get('state'), 'eso_rows': w.get('eso_rows'),
         'transfer': w.get('transfer_count'), 'frames': len(names), 'span_h': span_hours(names), 'targets': targets,
         'lcs': sorted(lcs), 'jobs': {k: j['id'] for k, j in latest.items()}, 'attention': False}
    pl, dl = latest.get('pipeline'), latest.get('download')
    sso = ctx.cfg['telescopes'][tel].get('source') == 'eso'
    parts = [p for p in [_frames_text(s)] if p]

    if active:
        state = 'processing'
        parts.append('{} job{} {} ({})'.format(len(active), '' if len(active) == 1 else 's',
                                                  'running' if any(j['state'] == 'running' for j in active)
                                                  else 'queued', ', '.join(str(j['id']) for j in active)))
    elif dl and dl['state'] == 'failed' and (not pl or pl['state'] != 'done'):
        state, s['attention'] = 'failed', True
        parts.append('download failed (job {}): {}'.format(dl['id'], dl['failure_signature'] or dl['failure_kind']))
    elif pl and pl['state'] == 'failed':
        state, s['attention'] = 'failed', True
        parts.append('pipeline failed (job {}): {}'.format(pl['id'], pl['failure_signature'] or pl['failure_kind']))
    elif pl and pl['state'] == 'cancelled':
        state, s['attention'] = 'cancelled', True
        parts.append('pipeline job {} cancelled'.format(pl['id']))
    elif (pl and pl['state'] == 'done') or w.get('state') == 'processed' or (lcs and not pl):
        when = ' · done {} (job {})'.format(fmt_time(pl['finished_at']), pl['id']) if pl else ''
        if lcs:
            state = 'done'
            # ESO object names and Gaia-ID output directories are different namespaces: compare counts only
            n_t = max(len(targets), len(lcs))
            parts.append('light curves for {} of {} target{}{}'.format(len(lcs), n_t, '' if n_t == 1 else 's', when))
        else:
            state, s['attention'] = 'no-lc', True
            parts.append('pipeline finished, no light curves' + when)
    elif w.get('state') in ('waiting', 'arriving'):
        state = 'waiting'
        late = ctx.now >= dt.datetime.combine(parse_night(night) + dt.timedelta(days=1),
                                              hhmm(ctx.n['late_after_utc']))
        s['attention'] = bool(late and w.get('transfer_count'))
        if w.get('transfer_count'):
            parts.append('{} of the {} frames in the transfer log at ESO so far'.format(
                w.get('eso_rows') or 0, w['transfer_count']))
        else:
            parts.append('waiting for frames at ESO' if sso else 'waiting for the night to arrive')
    elif not ctx.observed(tel, night) and ctx.chilean(tel, night):
        state = 'chilean'
        parts.append('Chilean night ({:.1f} h planned for Chile; not ours to download)'.format(
            ctx.plan(tel, night)['chilean'] / 60.0))
    elif not ctx.observed(tel, night):
        state = 'no-data'
        if sso:
            sib = ctx.siblings_observed(tel, night)
            parts.append('no data; {} observed'.format(_and(sib)) if sib else
                         'no data; no SSO telescope observed (weather or the site)')
        else:
            parts.append('no data')
    elif names and not pl:
        state, s['attention'] = 'not-run', True
        parts.append('no pipeline run')
    else:
        state = 'unknown'
        parts.append('state unknown (watcher: {})'.format(w.get('state') or 'not seen'))
    s['state'], s['text'] = state, ' · '.join(parts)
    return s


# ---------------------------------------------------------------------------------------------- queue, disks

def queue_summary(ctx, hours=24):
    store, now = ctx.store, ctx.now
    since = iso(now - dt.timedelta(hours=hours))
    finished = []
    for r in store.conn.execute("SELECT * FROM jobs WHERE finished_at >= ? AND state IN ('done', 'failed', "
                                "'cancelled') ORDER BY id", (since,)):
        d = dict(r)
        d['meta'] = json.loads(d['meta'] or '{}')
        finished.append(d)
    failures = Counter((r['failure_signature'] or r['failure_kind'] or 'unknown') for r in finished
                       if r['state'] == 'failed')
    new = set()
    for sig in failures:
        seen = store.conn.execute('SELECT 1 FROM jobs WHERE failure_signature = ? AND finished_at < ? LIMIT 1',
                                  (sig, since)).fetchone()
        if not seen:
            new.add(sig)
    info = json.loads(store.kv_get('dispatcher') or '{}')
    last_start = store.conn.execute('SELECT MAX(started_at) AS t FROM jobs').fetchone()['t']
    actions = [r for r in store.lookback_rows()
               if (r.get('checked_at') or '') >= since and r.get('last_action') not in (None, 'ok')]
    return {
        'since': since, 'hours': hours, 'finished_jobs': finished,
        'finished': Counter(r['state'] for r in finished),
        'by_class': Counter((r['class'], r['state']) for r in finished),
        'failures': failures, 'new_signatures': new,
        'queued': store.jobs(states=['queued']), 'running': store.jobs(states=['running']),
        'paused': store.paused(), 'draining': store.draining(),
        'heartbeat': parse_iso(info.get('heartbeat')), 'last_start': parse_iso(last_start),
        'lookback_actions': actions,
        'eso_errors': store.events(since=since, kinds=['eso-error']),
    }


DEFAULT_LABELS = {'pipeline': 'reruns', 'download': 'downloads', 'fetch': 'downloads', 'framepack': 'frame packing'}


def tidy_note(note, limit=40):
    """A note as a group label: without parenthesised detail or what follows a colon ('16-bit archive sweep: v2 510
    frames in 2 folders' -> '16-bit archive sweep')."""
    s = re.sub(r'\([^)]*\)', '', note or '')
    s = ' '.join(s.split(':')[0].split()).strip(' ;,.-')
    return s if len(s) <= limit else s[:limit - 1].rstrip() + '…'


def campaign_label(job):
    """What the digest groups a job under: its campaign (`q add --campaign`, `q campaign`), else nightly (P0) or
    look-back (P1), else its tidied note, else its kind."""
    label = str((job.get('meta') or {}).get('campaign') or '').strip()
    if label:
        return label
    if job['class'] == 'P0':
        return 'nightly'
    if job['class'] == 'P1':
        return 'look-back'
    return tidy_note(job.get('note')) or DEFAULT_LABELS.get(job['kind'], job['kind'])


def night_span(nights):
    nights = sorted({n for n in nights if n})
    if not nights:
        return ''
    if nights[0] == nights[-1]:
        return fmt_night(nights[0])
    a, b = (parse_night(n).strftime('%b %Y') for n in (nights[0], nights[-1]))
    return a if a == b else '{} – {}'.format(a, b)


def work_summary(ctx, q):
    """The queue's work over the window and what is left, grouped by (class, campaign_label). Light curves are
    counted per telescope-night of the finished pipeline jobs."""
    now = ctx.now
    done, left = {}, {}
    for j in q['finished_jobs']:
        g = done.setdefault((j['class'], campaign_label(j)), {
            'done': 0, 'failed': 0, 'cancelled': 0, 'telescopes': set(), 'nights': [], 'runs': set()})
        g[j['state']] += 1
        if j['telescope']:
            g['telescopes'].add(j['telescope'])
        g['nights'].append(j['night'])
        if j['kind'] == 'pipeline' and j['state'] == 'done' and j['telescope'] and j['night']:
            g['runs'].add((j['telescope'], j['night']))
    for g in done.values():
        g['lc_runs'] = sum(1 for tel, night in g['runs'] if ctx.lcs(tel, night))
    for j in q['running'] + q['queued']:
        g = left.setdefault((j['class'], campaign_label(j)), {
            'running': 0, 'queued': 0, 'minutes': 0.0, 'telescopes': set(), 'nights': [], 'retry': 0, 'jobs': []})
        g[j['state']] += 1
        g['minutes'] += float(j['est_minutes'] or 0)
        if j['telescope']:
            g['telescopes'].add(j['telescope'])
        g['nights'].append(j['night'])
        g['jobs'].append(j)
        if j['state'] == 'queued' and not job_is_due(j, now):
            g['retry'] += 1
    done_h = sum(float(j['est_minutes'] or 0) for j in q['finished_jobs'] if j['state'] == 'done') / 60.0
    left_h = sum(g['minutes'] for g in left.values()) / 60.0
    return {'done': done, 'left': left, 'done_h': done_h, 'left_h': left_h,
            'days': left_h / done_h if done_h > 0 else None}


def disks(ctx, usage=disk_usage):
    out, seen = [], set()
    for label, path in (('raw data', os.path.join(ctx.basedir, 'Observations')),
                        ('PipelineOutput', os.path.join(ctx.basedir, 'PipelineOutput')),
                        ('queue', ctx.cfg['queue_root'])):
        try:
            dev, free, total = usage(path)
        except OSError:
            continue
        if dev in seen or not total:
            continue
        seen.add(dev)
        out.append({'label': label, 'path': path, 'free_tib': free / TIB, 'total_tib': total / TIB,
                    'free_pct': 100.0 * free / total})
    return out


# ---------------------------------------------------------------------------------------------- conditions

# ---------------------------------------------------------------------------------------------- shutter state

SHUTTER_KINDS = {
    'CAMSHUTT': lambda d: 'but CAMSHUTT = {} on {} of {} frames'.format(d['header'], d['n_frames'], d['of']),
    'override': lambda d: 'but the override window {} to {} says {}'.format(
        d['window'][0][:16], d['window'][1][:16], 'shutter used' if d['override'] == 'T' else
        'shutterless' + (' at {} s'.format(d['override_t_read']) if d.get('override_t_read') else '')),
    'READOUTM': lambda d: 'but READOUTM names {}'.format(', '.join(d['header'])),
    'no-header': lambda d: 'with no CAMSHUTT and no override window, so not de-smeared',
}


def _epoch(t):
    return (t - dt.datetime(1970, 1, 1)).total_seconds()


def shutter_notes(ctx, hours=24, limit=12):
    """Digest lines for nights reduced in the last `hours` whose measured shutter state (reduction/shutter.json,
    written at T2 by calibration/shutter_detect.py) disagrees with CAMSHUTT, the override table or READOUTM, or
    is shutterless where no header says so. One line per disagreement, at most `limit`, then a count."""
    cutoff = _epoch(ctx.now - dt.timedelta(hours=hours))
    lines = []
    for version in ('v3', 'v2'):
        pattern = os.path.join(ctx.basedir, 'PipelineOutput', version, '*', 'output', '*', 'reduction', 'shutter.json')
        for path in sorted(glob.glob(pattern)):
            try:
                if os.path.getmtime(path) < cutoff:
                    continue
                with open(path) as f:
                    rec = json.load(f)
            except (OSError, ValueError):
                continue
            tel = path.split(os.sep)[-5]
            for d in rec.get('disagreements') or []:
                say = SHUTTER_KINDS.get(d.get('kind'))
                if say is None:
                    continue
                measured = d.get('measured') or '?'
                if measured == 'shutterless' and d.get('measured_t_read'):
                    measured += ' ({} s)'.format(d['measured_t_read'])
                try:
                    detail = say(d)
                except (KeyError, IndexError, TypeError):
                    detail = 'disagrees with {}'.format(d.get('kind'))
                lines.append('  • {}, night of {}: frames measure {}{} {}'.format(
                    tel, fmt_night(rec.get('night') or path.split(os.sep)[-3]), measured,
                    ' [{}]'.format(d['confidence']) if d.get('confidence') else '', detail))
    if len(lines) > limit:
        lines = lines[:limit] + ['  … and {} more'.format(len(lines) - limit)]
    return lines


def condition(key, kind, title, detail, telescope=None, night=None):
    return {'key': key, 'kind': kind, 'title': title, 'detail': detail, 'telescope': telescope, 'night': night}


def failed_conditions(ctx):
    """The newest job of each (telescope, night, kind) failed, and the night had a P0 job: nightly work that needs
    a person. A later success (a look-back rerun, a requeue) resolves it."""
    first = night_str(last_night(ctx.now) - dt.timedelta(days=int(ctx.n['failed_lookback_nights'])))
    latest, had_p0 = {}, set()
    for r in ctx.store.conn.execute('SELECT * FROM jobs WHERE night >= ? ORDER BY id', (first,)):
        key = (r['telescope'], r['night'], r['kind'])
        latest[key] = r
        if r['class'] == 'P0':
            had_p0.add(key)
    out = []
    for (tel, night, kind), r in sorted(latest.items()):
        if r['state'] != 'failed' or (tel, night, kind) not in had_p0:
            continue
        detail = ('Job {} ({}) failed {}: [{}] {}\nLog: {}\nOnce fixed: q requeue {}'.format(
            r['id'], r['class'], fmt_time(r['finished_at']), r['failure_kind'] or '-',
            r['failure_signature'] or 'exit {}'.format(r['exit_code']), r['log_path'] or '-', r['id']))
        out.append(condition('failed:{}:{}:{}'.format(tel, night, kind), 'failed',
                             '{}, night of {}: {} job failed'.format(tel, fmt_night(night), kind), detail, tel, night))
    return out


def late_conditions(ctx):
    """The observatory's transfer log lists frames that have still not all reached ESO by late_after_utc the next
    day, and the watcher is still waiting."""
    out = []
    ln = last_night(ctx.now)
    for night in (night_str(ln), night_str(ln - dt.timedelta(days=1))):
        deadline = dt.datetime.combine(parse_night(night) + dt.timedelta(days=1), hhmm(ctx.n['late_after_utc']))
        if ctx.now < deadline:
            continue
        for tel in ctx.sso:
            w = ctx.watch(tel, night) or {}
            transfer, rows = w.get('transfer_count'), w.get('eso_rows') or 0
            if w.get('state') in ('waiting', 'arriving') and transfer and rows < transfer:
                out.append(condition(
                    'late:{}:{}'.format(tel, night), 'late',
                    '{}, night of {}: frames missing at ESO'.format(tel, fmt_night(night)),
                    "The observatory's transfer log lists {} frames; ESO had {} at the last poll ({} UTC). The watcher "
                    'queues what is there once the count stops changing, and the look-back fetches the rest.'.format(
                        transfer, rows, fmt_time(w.get('last_poll'))), tel, night))
    return out


def dark_conditions(ctx):
    """A telescope with no survey data for dark_min_nights nights on which at least two other SSO telescopes
    observed (any number of nights for Artemis, which has no sister telescope). Chilean nights do not count
    and do not break the run."""
    out = []
    ln = last_night(ctx.now)
    for tel, tcfg in ctx.cfg['telescopes'].items():
        run, day, exhausted = [], ln, True
        for _ in range(int(ctx.n['dark_lookback_nights'])):
            night = night_str(day)
            day -= dt.timedelta(days=1)
            obs = ctx.observed(tel, night)
            if obs is None:
                continue
            if obs:
                exhausted = False
                break
            run.append(night)
        if not run:
            continue
        sso = tcfg.get('source') == 'eso'
        chile = [n for n in run if ctx.chilean(tel, n)]
        specific = [n for n in run if n not in chile and (not sso or len(ctx.siblings_observed(tel, n)) >= 2)]
        if len(specific) < int(ctx.n['dark_min_nights']):
            continue
        first = run[-1]
        title = '{}: no data since {}{}'.format(tel, 'at least ' if exhausted else '', fmt_night(first))
        if sso:
            detail = ('{} nights without survey data, up to the night of {}. On {} of them at least two of the other '
                      'SSO telescopes observed, so it is not the weather.'.format(len(run), fmt_night(run[0]),
                                                                                    len(specific)))
        else:
            detail = ('{} nights without data, up to the night of {}. {} has no sister telescope, so the weather is '
                      'not ruled out.'.format(len(run), fmt_night(run[0]), tel))
        if chile:
            detail += ' {} of the nights were Chilean time.'.format(len(chile))
        out.append(condition('dark:{}'.format(tel), 'dark', title, detail, tel))
    return out


def queue_conditions(ctx, q):
    out = []
    now = ctx.now
    beat = q['heartbeat']
    if beat is None or now - beat > dt.timedelta(minutes=float(ctx.n['dispatcher_stale_minutes'])):
        out.append(condition('dispatcher', 'queue', 'Job queue dispatcher not running',
                             'Last heartbeat: {}. The watchdog (cron, every 5 minutes) should restart it; see '
                             'queue/logs/watchdog.log.'.format(fmt_time(beat) if beat else 'never')))
    due = [j for j in q['queued'] if j['class'] not in q['paused'] and job_is_due(j, now)]
    stall = dt.timedelta(hours=float(ctx.n['stall_hours']))
    if due and not q['running'] and not q['draining'] and (q['last_start'] is None or now - q['last_start'] >= stall):
        out.append(condition('queue-stalled', 'queue', 'Job queue stalled',
                             '{} jobs are due, none is running, and none has started since {}. `q status` shows '
                             'what each is waiting for.'.format(len(due), fmt_time(q['last_start'])
                                                                if q['last_start'] else 'the queue began')))
    return out


def signature_conditions(ctx, q):
    out = []
    for sig, count in q['failures'].items():
        if count >= int(ctx.n['signature_alert_count']):
            out.append(condition('signature:' + hashlib.sha1(sig.encode()).hexdigest()[:10], 'signature',
                                 'Repeated job failure: {}'.format(sig[:80]),
                                 '{} jobs failed with this signature in the last {} hours. `q signatures` groups '
                                 'them; `q requeue --signature` reruns them once fixed.'.format(count, q['hours'])))
    return out


def disk_conditions(ctx, disk_list):
    out = []
    for d in disk_list:
        if d['free_tib'] < float(ctx.n['disk_min_free_tb']) or d['free_pct'] < float(ctx.n['disk_min_free_pct']):
            out.append(condition('disk:{}'.format(d['label']), 'disk', '{} disk nearly full'.format(d['label']),
                                 '{:.1f} TiB free of {:.0f} TiB ({:.1f}%) at {}.'.format(
                                     d['free_tib'], d['total_tib'], d['free_pct'], d['path'])))
    return out


def conditions(ctx, q, disk_list):
    return (failed_conditions(ctx) + late_conditions(ctx) + dark_conditions(ctx) + queue_conditions(ctx, q) +
            signature_conditions(ctx, q) + disk_conditions(ctx, disk_list))


# ---------------------------------------------------------------------------------------------- the digest

def _incident_line(i, now):
    line = '  • {} (open since {})'.format(i['title'], fmt_time(i['opened_at']))
    if i.get('ack_until') and i['ack_until'] > iso(now):
        line += ' — acknowledged until {}{}'.format(i['ack_until'][:10], ': ' + i['ack_note'] if i.get('ack_note')
                                                    else '')
    return line


def _fit(label, width):
    if len(label) <= width:
        return label
    cut = label[:width - 1]
    return (cut[:cut.rfind(' ')] if ' ' in cut else cut).rstrip(' ;,.-') + '…'


def _top(groups, size, top):
    """The `top` biggest groups, in class order then size, and a line for the rest."""
    keys = sorted(groups, key=lambda k: -size(groups[k]))
    shown = sorted(keys[:top], key=lambda k: (k[0], -size(groups[k])))
    rest = keys[top:]
    more = '  … and {} smaller group{} ({} jobs)'.format(len(rest), '' if len(rest) == 1 else 's',
                                                        sum(size(groups[k]) for k in rest)) if rest else None
    return shown, more


def render_work(q, now, top=8):
    """'Processed in the last 24 h' and 'Remaining': the queue's work grouped by class and campaign."""
    w = q['work']
    out = []
    labels = [k[1] for k in list(w['done']) + list(w['left'])]
    width = min(max([len(l) for l in labels] + [8]), 32)
    row = '  {:<3} {:<' + str(width) + '}  {}'

    n_done = len(q['finished_jobs'])
    out.append('Processed in the last {} h: {}'.format(q['hours'], '{} job{}'.format(n_done, '' if n_done == 1 else 's')
                                                       if n_done else 'nothing'))
    shown, more = _top(w['done'], lambda g: g['done'] + g['failed'] + g['cancelled'], top)
    for cls, label in shown:
        g = w['done'][(cls, label)]
        parts = ['{} done'.format(g['done'])]
        parts += ['{} {}'.format(g[s], s) for s in ('failed', 'cancelled') if g[s]]
        text = ', '.join(parts)
        if g['telescopes']:
            text += ' · ' + ', '.join(sorted(g['telescopes']))
        if night_span(g['nights']):
            text += ' · ' + night_span(g['nights'])
        if g['runs']:
            text += ' · light curves on {} of {} night{}'.format(g['lc_runs'], len(g['runs']),
                                                                 '' if len(g['runs']) == 1 else 's')
        out.append(row.format(cls, _fit(label, width), text))
    if more:
        out.append(more)

    n_left = len(q['running']) + len(q['queued'])
    if not n_left:
        out.append('Remaining: nothing queued')
        return out
    head = 'Remaining: {} job{}, ≈{:.0f} h of estimated work'.format(n_left, '' if n_left == 1 else 's', w['left_h'])
    head += (', ≈{:.1f} days at the last {} h\'s pace'.format(w['days'], q['hours']) if w['days'] is not None
             else ', and nothing finished in the last {} h'.format(q['hours']))
    if q['paused']:
        head += '; paused {}'.format(', '.join(sorted(q['paused'])))
    if q['draining']:
        head += '; draining'
    out.append(head)
    shown, more = _top(w['left'], lambda g: g['running'] + g['queued'], top)
    for cls, label in shown:
        g = w['left'][(cls, label)]
        text = ', '.join(['{} running'.format(g['running'])] * bool(g['running']) +
                         ['{} queued'.format(g['queued'])] * bool(g['queued']))
        text += ' · ≈{:.0f} h'.format(g['minutes'] / 60.0)
        if len(g['jobs']) <= 2:
            text += ' · ' + ', '.join('{} {}'.format(j['telescope'] or '-', fmt_night(j['night']) if j['night']
                                                     else '') .strip() for j in g['jobs'])
        else:
            text += ' · {}{}'.format(', '.join(sorted(g['telescopes'])) + ' · ' if g['telescopes'] else '',
                                     night_span(g['nights']))
        if cls in q['paused']:
            text += ' · paused'
        if g['retry']:
            text += ' · {} waiting to retry'.format(g['retry'])
        out.append(row.format(cls, _fit(label, width), text.rstrip(' ·')))
    if more:
        out.append(more)
    return out


def render_digest(cfg, night, statuses, q, disk_list, open_incidents, resolved, now, report_path=None,
                  shutter=None):
    """Subject and plain-text body of the morning digest. shutter: lines from shutter_notes()."""
    incident_tels = {i['telescope'] for i in open_incidents if i.get('telescope')}
    lines, flagged = [], 0
    for s in statuses:
        mark = s['attention'] or s['telescope'] in incident_tels
        flagged += mark
        lines.append('{} {:<9} {}'.format('✗' if mark else ' ', s['telescope'], s['text']))
    others = [i for i in open_incidents if not i.get('telescope')]
    n = flagged + len(others)
    verdict = 'all well' if n == 0 else '{} need{} attention'.format(n, 's' if n == 1 else '')
    subject = '{} Night of {}: {}'.format(cfg['notify']['subject_prefix'], fmt_night(night), verdict)

    start = dt.datetime.combine(parse_night(night), dt.time(15, 0))
    body = ['Night of {} · {}'.format(fmt_night(night), verdict),
            '(the night from {} to {} UTC; written {} UTC)'.format(fmt_time(start), fmt_time(start + dt.timedelta(1)),
                                                                  fmt_time(now)),
            ''] + lines + ['']
    if open_incidents:
        body.append('Open incidents ({})'.format(len(open_incidents)))
        body += [_incident_line(i, now) for i in open_incidents]
        body.append('')
    if resolved:
        body.append('Resolved since the last digest')
        body += ['  • {} (resolved {})'.format(i['title'], fmt_time(i['resolved_at'])) for i in resolved]
        body.append('')
    if shutter:
        body.append('Shutter state measured against the headers (nights reduced in the last 24 h)')
        body += shutter
        body.append('')

    body += render_work(q, now)
    if q['failures']:
        body.append('Failures by signature:')
        for sig, count in q['failures'].most_common(8):
            body.append('  {} × {}{}'.format(count, sig[:110], '  [new]' if sig in q['new_signatures'] else ''))
    if q['lookback_actions']:
        body.append('Look-back: ' + '; '.join('{} {} {}'.format(r['telescope'], r['night'], r['last_action'])
                                              for r in q['lookback_actions'][:10]))
    if q['eso_errors']:
        body.append('Watcher: {} ESO query error{} in the last {} h'.format(
            len(q['eso_errors']), '' if len(q['eso_errors']) == 1 else 's', q['hours']))
    body.append('')
    if disk_list:
        body.append('Disks: ' + ' · '.join('{} {:.1f} TiB free ({:.0f}%)'.format(d['label'], d['free_tib'],
                                                                               d['free_pct']) for d in disk_list))
        body.append('')
    if report_path:
        body.append('Daily report: {}'.format(report_path))
    body.append('q notify incidents · q notify ack <id> --until <YYYY-MM-DD> --note <text> · q status')
    return subject, '\n'.join(body) + '\n'
