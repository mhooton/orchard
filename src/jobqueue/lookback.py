"""Daily look-back over the last 30 nights: frames that reached ESO late, and nights with no products.

For each SSO telescope one ESO query lists every frame of the window; each night is compared with its
directory on disk (cube-aware, see nights.diff_eso_disk). Then, per night:

  - frames missing on disk  -> a P1 fetch job (add-only, fetch.py), which queues a P1 rerun of the night
                               itself if it added science frames;
  - science frames on disk but no light curve in v2 or v3 -> a P1 rerun, at most max_reruns_per_night;
  - light curves only in v3 -> reported (a T12 promotion question, never automatic).

Nights with a queued or running job are left alone, and the newest night of all belongs to the watcher.
Artemis has no ESO side, so only its products are checked.
"""
import datetime as dt
import json
import os

from . import jobspec
from .config import eso_telescopes, queue_paths
from .eso import group_by_night
from .nights import (artemis_targets, diff_eso_disk, list_frames, night_dir, normalise_target, products_state,
                     science_count)
from .util import iso, last_night, night_str, utcnow, write_json_atomic


class Lookback:
    def __init__(self, cfg, store, eso, now=None, log=print, dry_run=False):
        self.cfg = cfg
        self.lb = cfg['lookback']
        self.store = store
        self.eso = eso
        self.now = now or utcnow()
        self.log = log
        self.dry_run = dry_run
        self.budget = int(self.lb['max_fetch_rows_per_run'])

    def nights(self):
        ln = last_night(self.now)
        return [night_str(ln - dt.timedelta(days=k)) for k in range(1, int(self.lb['nights']) + 1)]

    def run(self):
        nights = self.nights()
        out = []
        for tel in eso_telescopes(self.cfg):
            prog = self.cfg['telescopes'][tel]['prog_id']
            try:
                rows = self.eso.frames(prog, nights[-1], nights[0])
            except Exception as e:
                msg = '{}: ESO query failed ({}); skipped today'.format(tel, e)
                self.log(msg)
                out.append({'telescope': tel, 'night': None, 'action': msg})
                continue
            by_night = group_by_night(rows)
            for night in nights:
                out.append(self.check_eso_night(tel, night, by_night.get(night, [])))
        for tel, tcfg in self.cfg['telescopes'].items():
            if tcfg.get('source') == 'raw_dir':
                for night in nights:
                    out.append(self.check_raw_night(tel, night))
        return out

    # ------------------------------------------------------------------

    def _row(self, tel, night):
        return self.store.lookback_get(tel, night) or {'telescope': tel, 'night': night, 'fetches': 0, 'reruns': 0,
                                                        'data': {}}

    def _save(self, lb):
        if not self.dry_run:
            self.store.lookback_put(lb)

    def _busy(self, tel, night):
        return bool(self.store.active_for(tel, night))

    def check_eso_night(self, tel, night, rows):
        lb = self._row(tel, night)
        result = {'telescope': tel, 'night': night, 'eso_rows': len(rows)}
        if self._busy(tel, night):
            result['action'] = 'skipped: a job for this night is queued or running'
            return result
        path = night_dir(self.cfg['basedir'], tel, night)
        names = list_frames(path)
        settled = set((lb.get('data') or {}).get('settled', []))
        present, missing, _ = diff_eso_disk(rows, names)
        missing = [r for r in missing if r['dp_id'] not in settled]
        products = products_state(self.cfg['basedir'], tel, night) if names else 'none'
        sci_present = science_count(present)
        result.update(disk_files=len(names), missing_rows=len(missing), missing_science=science_count(missing),
                      products=products)
        lb.update(checked_at=iso(self.now), eso_rows=len(rows), disk_files=len(names), missing_rows=len(missing),
                  missing_science=science_count(missing), products=products)

        action = 'ok'
        if missing:
            if int(lb['fetches']) >= int(self.lb['max_fetches_per_night']):
                action = 'missing {} rows, but {} fetches already tried'.format(len(missing), lb['fetches'])
            elif self.budget <= 0:
                action = 'missing {} rows; run budget used up, tomorrow'.format(len(missing))
            else:
                take = missing[:min(int(self.lb['max_fetch_rows_per_night']), self.budget)]
                self.budget -= len(take)
                action = self._queue_fetch(lb, tel, night, take)
        elif sci_present and products == 'none':
            if int(lb['reruns']) >= int(self.lb['max_reruns_per_night']):
                action = 'no light curves; {} automatic reruns already tried'.format(lb['reruns'])
            else:
                objects = sorted({normalise_target(r['object']) for r in present
                                  if (r.get('dp_type') or '').upper() == 'OBJECT' and (r.get('object') or '').strip()})
                action = self._queue_rerun(lb, tel, night, objects, len(names), 'no light curves in v2 or v3')
        elif products == 'v3_only':
            action = 'light curves only in v3: needs a T12 promotion (not automatic)'
        lb['last_action'] = action
        self._save(lb)
        result['action'] = action
        return result

    def check_raw_night(self, tel, night):
        lb = self._row(tel, night)
        result = {'telescope': tel, 'night': night}
        path = night_dir(self.cfg['basedir'], tel, night)
        if not os.path.isdir(path):
            result['action'] = 'no night directory'
            return result
        if self._busy(tel, night):
            result['action'] = 'skipped: a job for this night is queued or running'
            return result
        names = list_frames(path)
        targets = artemis_targets(names)
        products = products_state(self.cfg['basedir'], tel, night)
        result.update(disk_files=len(names), products=products)
        lb.update(checked_at=iso(self.now), disk_files=len(names), products=products)
        action = 'ok'
        if targets and products == 'none':
            if int(lb['reruns']) >= int(self.lb['max_reruns_per_night']):
                action = 'no light curves; {} automatic reruns already tried'.format(lb['reruns'])
            else:
                action = self._queue_rerun(lb, tel, night, targets, len(names), 'no light curves in v2 or v3')
        elif products == 'v3_only':
            action = 'light curves only in v3: needs a T12 promotion (not automatic)'
        lb['last_action'] = action
        self._save(lb)
        result['action'] = action
        return result

    # ------------------------------------------------------------------

    def _queue_fetch(self, lb, tel, night, rows):
        n_sci = science_count(rows)
        if self.dry_run:
            return 'would fetch {} rows ({} science), add-only'.format(len(rows), n_sci)
        fetch_dir = queue_paths(self.cfg)['fetch']
        os.makedirs(fetch_dir, exist_ok=True)
        path = os.path.join(fetch_dir, '{}_{}_{}.json'.format(tel, night, self.now.strftime('%Y%m%dT%H%M%S')))
        write_json_atomic(path, [{k: r.get(k) for k in ('dp_id', 'dp_type', 'object', 'origfile', 'mjd_obs')}
                                 for r in rows])
        jid, created = jobspec.add(self.store, jobspec.fetch_job(self.cfg, 'P1', tel, night, path, len(rows), n_sci))
        if created:
            lb['fetches'] = int(lb['fetches']) + 1
        return 'fetch {} rows ({} science): job {}{}'.format(len(rows), n_sci, jid, '' if created else ' (existing)')

    def _queue_rerun(self, lb, tel, night, targets, frames, why):
        if self.dry_run:
            return 'would rerun ({}), targets {}'.format(why, ' '.join(targets) or '-')
        jid, created = jobspec.add(self.store, jobspec.pipeline_job(
            self.cfg, 'P1', tel, night, lock_targets=targets, frames=frames, source='lookback',
            priority=float(night), note=why))
        if created:
            lb['reruns'] = int(lb['reruns']) + 1
        return 'rerun ({}): job {}{}'.format(why, jid, '' if created else ' (existing)')


QUIET = ('ok', 'no night directory')


def format_report(results):
    lines = []
    for r in results:
        if r.get('action', 'ok') in QUIET:
            continue
        lines.append('{:9} {:8} eso={:<5} disk={:<5} missing={:<4} sci_missing={:<4} products={:8} {}'.format(
            r['telescope'], r.get('night') or '-', r.get('eso_rows', '-'), r.get('disk_files', '-'),
            r.get('missing_rows', '-'), r.get('missing_science', '-'), r.get('products', '-'), r['action']))
    ok = sum(1 for r in results if r.get('action') in QUIET)
    lines.append('{} telescope-nights checked, {} need nothing'.format(len(results), ok))
    return '\n'.join(lines)


def results_json(results):
    return json.dumps(results, indent=1, default=str)
