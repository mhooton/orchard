"""The ESO watcher: queue each night as soon as its frames are complete, instead of at a fixed cron slot.

Run from cron every 30 minutes through the afternoon. For the last `watch_days` nights and each SSO
telescope it asks ESO (authenticated, read-only) how many survey frames the night has. A night is ready
when ESO's count matches Observations/<TEL>/transfer_log.txt, or, with no transfer-log entry, when the count
and ESO's last-modified time have not moved for two polls. Then it queues two P0 jobs: an SSO_download.py
run for that one night, and the pipeline run after it.

  - A night with no frames at ESO and no transfer-log entry (Io since 2026-08-29) gets one quiet note,
    not 24 download attempts.
  - Below the transfer-log count it keeps waiting; after partial_stable_minutes without change it queues
    what is there, and the daily look-back fetches the rest.
  - Artemis does not come through ESO: its night directory is ready when Data_Download.txt (written at the
    end of each transfer, ~10:45-11:45 UTC) lists only files that are present, or when the directory has
    stopped changing; by artemis_deadline_utc it is queued with whatever is there.
  - A night that already has a pipeline log (the cron or someone else ran it) is left to the look-back,
    so switching over never reprocesses the night before last.
  - Shadow mode needs nothing special here: the shadow configuration points at a separate database, and
    the shadow dispatcher only pretends to run what is queued.
"""
import datetime as dt
import os

from . import jobspec
from .config import camera, eso_telescopes
from .nights import artemis_targets, list_frames, night_dir, pipeline_logs, read_transfer_log
from .util import hhmm, iso, last_night, night_str, parse_iso, parse_night, utc_from_ts, utcnow


FINAL = ('enqueued', 'processed')


def _after(night, hhmm_utc, now, days=1):
    """True once now is past hh:mm UTC on the day `days` after the night."""
    day = parse_night(night) + dt.timedelta(days=days)
    return now >= dt.datetime.combine(day, hhmm(hhmm_utc))


def _new_row(telescope, night):
    return {'telescope': telescope, 'night': night, 'state': 'waiting', 'polls_unchanged': 0, 'noted': 0,
            'data': {}}


def _minutes(a, b):
    return (a - b).total_seconds() / 60.0


def _track(w, count, marker, now):
    """Update change tracking; returns minutes since the count or marker last changed."""
    if count != w.get('eso_rows') or marker != w.get('eso_last_mod'):
        w['last_change'] = iso(now)
        w['polls_unchanged'] = 0
        if count and not w.get('first_seen'):
            w['first_seen'] = iso(now)
    else:
        w['polls_unchanged'] = int(w.get('polls_unchanged') or 0) + 1
    w['eso_rows'], w['eso_last_mod'], w['last_poll'] = count, marker, iso(now)
    return _minutes(now, parse_iso(w['last_change']) or now)


def _note(store, w, message, log):
    if w.get('noted'):
        return
    w['noted'], w['note'] = 1, message
    store.event('note', source='watcher', telescope=w['telescope'], night=w['night'], message=message)
    log(message)


class Watcher:
    def __init__(self, cfg, store, eso, now=None, log=print):
        self.cfg = cfg
        self.w = cfg['watcher']
        self.store = store
        self.eso = eso
        self.now = now or utcnow()
        self.log = log

    def nights(self):
        ln = last_night(self.now)
        return [night_str(ln - dt.timedelta(days=k)) for k in range(int(self.w['watch_days']))]

    def run(self):
        out = []
        tels = eso_telescopes(self.cfg)
        for night in self.nights():
            pending = [t for t in tels if (self.store.watch_get(t, night) or {}).get('state') not in FINAL]
            if pending:
                try:
                    summary = self.eso.night_summary(night)
                except Exception as e:
                    msg = 'ESO query for {} failed ({}); nothing changed, the next poll retries'.format(night, e)
                    self.store.event('eso-error', source='watcher', night=night, message=str(e)[:500])
                    self.log(msg)
                    out.append(msg)
                    summary = None
                if summary is not None:
                    for tel in pending:
                        prog = self.cfg['telescopes'][tel]['prog_id']
                        out.append(self.poll_eso(tel, night, summary.get(prog) or {'rows': 0, 'by_type': {},
                                                                                     'last_mod': None}))
            for tel, tcfg in self.cfg['telescopes'].items():
                if tcfg.get('source') == 'raw_dir':
                    out.append(self.poll_raw_dir(tel, night))
        return out

    # ------------------------------------------------------------------ SSO telescopes via ESO

    def poll_eso(self, tel, night, s):
        now, w = self.now, self.store.watch_get(tel, night) or _new_row(tel, night)
        if w['state'] in FINAL:
            return '{} {}: already {}'.format(tel, night, w['state'])
        rows, last_mod = int(s.get('rows') or 0), s.get('last_mod')
        transfer = read_transfer_log(self.cfg['basedir'], tel, night)
        since_change = _track(w, rows, last_mod, now)
        w['transfer_count'] = transfer
        w['data'] = dict(w.get('data') or {}, by_type=s.get('by_type') or {})

        if rows == 0:
            if transfer:
                w['state'] = 'waiting'
                line = '{} {}: transfer log says {} frames, none at ESO yet'.format(tel, night, transfer)
            else:
                late = _after(night, self.w['no_data_note_after_utc'], now)
                w['state'] = 'no_data' if late else 'waiting'
                line = '{} {}: no frames at ESO{}'.format(tel, night, '' if transfer is None else ' (transfer log 0)')
                if late:
                    _note(self.store, w, '{} {}: no frames at ESO and {}; nothing to download or process'.format(
                        tel, night, 'no transfer-log entry' if transfer is None else 'the transfer log says 0'),
                        self.log)
            self.store.watch_put(w)
            return line

        last_mod_age = _minutes(now, parse_iso(last_mod)) if last_mod else None
        quiet = last_mod_age is None or last_mod_age >= float(self.w['quiet_minutes'])
        stable = int(w['polls_unchanged']) >= 1 and since_change >= float(self.w['stable_minutes']) and quiet
        ready, reason = False, None
        if transfer and rows >= transfer:
            ready, reason = True, 'ESO has {} frames, transfer log {}'.format(rows, transfer)
        elif transfer and rows < transfer:
            if stable and since_change >= float(self.w['partial_stable_minutes']):
                ready, reason = True, ('ESO has {} of the {} frames in the transfer log, unchanged for {:.0f} min; '
                                       'the look-back fetches the rest'.format(rows, transfer, since_change))
        elif stable:
            ready, reason = True, 'ESO count {} unchanged for {:.0f} min (no transfer-log entry)'.format(
                rows, since_change)

        if not ready:
            w['state'] = 'arriving'
            self.store.watch_put(w)
            return '{} {}: {} frames at ESO{}, waiting ({} unchanged polls, last ESO change {} min ago)'.format(
                tel, night, rows, '' if not transfer else ' of {}'.format(transfer), w['polls_unchanged'],
                'n/a' if last_mod_age is None else '{:.0f}'.format(last_mod_age))

        try:
            objects = self.eso.night_objects(self.cfg['telescopes'][tel]['prog_id'], night)
        except Exception as e:
            self.store.watch_put(w)
            msg = '{} {}: ready but the target query failed ({}); retrying next poll'.format(tel, night, e)
            self.log(msg)
            return msg
        return self._enqueue(w, tel, night, rows, objects, reason, download=True)

    # ------------------------------------------------------------------ Artemis via its raw directory

    def poll_raw_dir(self, tel, night):
        now, w = self.now, self.store.watch_get(tel, night) or _new_row(tel, night)
        if w['state'] in FINAL:
            return '{} {}: already {}'.format(tel, night, w['state'])
        path = night_dir(self.cfg['basedir'], tel, night)
        deadline = _after(night, self.w['artemis_deadline_utc'], now)
        if not os.path.isdir(path):
            w['state'] = 'no_data' if deadline else 'waiting'
            if deadline:
                _note(self.store, w, '{} {}: no night directory by {} UTC; nothing to process'.format(
                    tel, night, self.w['artemis_deadline_utc']), self.log)
            w['last_poll'] = iso(now)
            self.store.watch_put(w)
            return '{} {}: no night directory yet'.format(tel, night)

        names = list_frames(path)
        dir_mtime = utc_from_ts(os.stat(path).st_mtime)
        since_change = _track(w, len(names), iso(dir_mtime), now)
        quiet = _minutes(now, dir_mtime) >= float(self.w['quiet_minutes'])
        ready, reason = False, None

        manifest = os.path.join(path, self.w['artemis_manifest'])
        if os.path.exists(manifest):
            listed = []
            with open(manifest, errors='replace') as f:
                listed = [l.strip() for l in f if l.strip().endswith(('.fts', '.fits', '.fit'))]
            present = sum(1 for l in listed if os.path.exists(os.path.join(path, l)))
            age = _minutes(now, utc_from_ts(os.stat(manifest).st_mtime))
            w['data'] = dict(w.get('data') or {}, manifest_listed=len(listed), manifest_present=present)
            if age >= float(self.w['quiet_minutes']) and listed and present >= len(listed):
                ready, reason = True, '{} lists {} files, all present'.format(self.w['artemis_manifest'], len(listed))
        if not ready and names and quiet and int(w['polls_unchanged']) >= 1 and \
                since_change >= float(self.w['stable_minutes']):
            ready, reason = True, 'directory unchanged for {:.0f} min ({} frames)'.format(since_change, len(names))
        if not ready and names and deadline:
            ready, reason = True, 'deadline {} UTC reached with {} frames'.format(self.w['artemis_deadline_utc'],
                                                                                len(names))
        if not ready:
            w['state'] = 'arriving' if names else 'waiting'
            self.store.watch_put(w)
            return '{} {}: {} frames on disk, waiting'.format(tel, night, len(names))
        return self._enqueue(w, tel, night, len(names), artemis_targets(names), reason, download=False)

    # ------------------------------------------------------------------ queueing

    def _enqueue(self, w, tel, night, count, targets, reason, download):
        logs = pipeline_logs(self.cfg['basedir'], tel, night)
        if logs:
            w.update(state='processed', ready_at=iso(self.now), ready_reason=reason)
            w['data'] = dict(w.get('data') or {}, pipeline_log=logs[-1])
            self.store.watch_put(w)
            msg = '{} {}: ready ({}), but already run outside the queue ({}); left to the look-back'.format(
                tel, night, reason, os.path.basename(logs[-1]))
            self.store.event('night-processed-elsewhere', source='watcher', telescope=tel, night=night, message=msg)
            self.log(msg)
            return msg
        dl_id = None
        if download:
            dl_id, _ = jobspec.add(self.store, jobspec.download_job(self.cfg, tel, night, count))
        # ESO rows are datacubes on SPIRIT nights, so they say nothing about frames; the dispatcher re-estimates
        # from the frames on disk when the job starts
        frames = None if download and camera(self.cfg, tel, night) == 'SPIRIT' else count
        pl_id, _ = jobspec.add(self.store, jobspec.pipeline_job(
            self.cfg, 'P0', tel, night, lock_targets=targets, frames=frames, source='watcher', depends_on=dl_id,
            note=reason))
        w.update(state='enqueued', ready_at=iso(self.now), ready_reason=reason, enqueued_at=iso(self.now),
                 download_job=dl_id, pipeline_job=pl_id)
        w['data'] = dict(w.get('data') or {}, targets=targets)
        self.store.watch_put(w)
        msg = '{} {}: ready ({}); queued {}pipeline job {} for targets {}'.format(
            tel, night, reason, 'download job {} then '.format(dl_id) if dl_id else '', pl_id,
            ' '.join(targets) or '(none listed)')
        self.store.event('night-ready', source='watcher', telescope=tel, night=night, message=msg,
                         data={'reason': reason, 'download_job': dl_id, 'pipeline_job': pl_id, 'count': count})
        self.log(msg)
        return msg
