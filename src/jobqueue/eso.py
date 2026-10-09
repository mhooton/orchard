"""Authenticated, read-only ESO archive queries for the watcher, the look-back and fetch jobs.

Authentication goes through download.request_eso.ESODownloader with the pipeline's credentials from
src/reporting/.env (anonymous TAP only sees public calibrations, so an unauthenticated answer would look
like a nearly empty night). Credentials are never printed: ESODownloader prints request errors, which can
contain the token URL with the password in its query string, so its output is captured and scrubbed.

ESO facts this relies on (checked 2026-10-09):
  - survey frames are prog_id 60.A-9009(A/B/C/D) = Io/Europa/Ganymede/Callisto, dp_id LIKE 'SPECU%';
    dp_type is OBJECT/BIAS/DARK/FLAT (OBJECT rows have come back as dp_cat SCIENCE and as CALIB, so dp_cat
    is never filtered on);
  - night D is [D 15:00, D+1 15:00) UTC in mjd_obs;
  - Callisto SPIRIT nights arrive as datacubes, one row per cube; det_ndit is empty, and origfile
    (SPECU4.20261009T015815_S_<target>_<filter>_<exp>s.fits) carries the cube's start time;
  - the ADQL service rejects FLOOR() and aliases in GROUP BY, so grouping is by plain columns.
"""
import contextlib
import csv
import io
import os
import time

from .util import mjd_to_night, night_bounds, night_str, to_mjd

TAP_TIMEOUT = 600


class EsoError(RuntimeError):
    pass


def read_eso_credentials(env_path):
    """ESO_USERNAME and ESO_PASSWORD from a dotenv file, without loading anything else into the environment."""
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
                if key in ('ESO_USERNAME', 'ESO_PASSWORD'):
                    value = value.strip()
                    if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
                        value = value[1:-1]
                    found[key] = value
    except OSError as e:
        raise EsoError('cannot read ESO credentials file {}: {}'.format(env_path, e.strerror))
    user = found.get('ESO_USERNAME') or os.getenv('ESO_USERNAME')
    password = found.get('ESO_PASSWORD') or os.getenv('ESO_PASSWORD')
    if not user or not password:
        raise EsoError('ESO_USERNAME / ESO_PASSWORD not set in {}'.format(env_path))
    return user, password


def scrub(text, secrets):
    from urllib.parse import quote, quote_plus
    for s in secrets:
        if not s:
            continue
        for form in {s, quote(s, safe=''), quote_plus(s)}:
            text = text.replace(form, '***')
    return text


def adql_str(s):
    return "'" + str(s).replace("'", "''") + "'"


class EsoArchive:
    """Thin read-only client. Construct inside the container; tests use a fake with the same methods."""

    def __init__(self, env_path, log=print, prog_ids=None):
        user, password = read_eso_credentials(env_path)
        self._secrets = [password, user]
        self.log = log
        from download.request_eso import ESODownloader  # the pipeline's own client
        import urllib3
        urllib3.disable_warnings()
        self._d = ESODownloader(user, password)
        self.prog_ids = prog_ids or ['60.A-9009(A)', '60.A-9009(B)', '60.A-9009(C)', '60.A-9009(D)']

    def _ensure_token(self):
        if self._d.token and time.time() + 1800 < self._d.token_expires:
            return
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ok = self._d.get_new_token()
        out = scrub(buf.getvalue().strip(), self._secrets)
        if not ok or not self._d.token:
            raise EsoError('ESO authentication failed: {}'.format(out[-300:]))

    def tap(self, query, retries=2):
        """Run a synchronous ADQL query; return rows as dicts. Raises EsoError rather than returning nothing."""
        import requests
        last = None
        for attempt in range(retries + 1):
            try:
                self._ensure_token()
                r = requests.post(self._d.tap_url + '/sync',
                                  data={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'MAXREC': '200000',
                                        'QUERY': query},
                                  headers={'Authorization': 'Bearer {}'.format(self._d.token)},
                                  verify=False, timeout=TAP_TIMEOUT)
                if r.status_code != 200:
                    raise EsoError('TAP HTTP {}: {}'.format(r.status_code, scrub(r.text[:300], self._secrets)))
                text = r.text
                if text.lstrip().startswith('<'):
                    raise EsoError('TAP error: {}'.format(scrub(text[:400], self._secrets)))
                return list(csv.DictReader(io.StringIO(text)))
            except EsoError as e:
                last = e
            except Exception as e:  # network errors from requests
                last = EsoError(scrub('{}: {}'.format(type(e).__name__, e), self._secrets))
            if attempt < retries:
                time.sleep(10 * (attempt + 1))
        raise last

    def night_summary(self, night):
        """{prog_id: {'rows': n, 'by_type': {dp_type: n}, 'last_mod': iso}} for one night, every programme."""
        a, b = (to_mjd(t) for t in night_bounds(night))
        q = ('SELECT prog_id, dp_type, COUNT(*) AS n, MAX(last_mod_date) AS last_mod FROM dbo.raw '
             'WHERE prog_id IN ({}) AND dp_id LIKE {} AND mjd_obs >= {:.6f} AND mjd_obs < {:.6f} '
             'GROUP BY prog_id, dp_type').format(','.join(adql_str(p) for p in self.prog_ids), adql_str('SPECU%'), a, b)
        out = {p: {'rows': 0, 'by_type': {}, 'last_mod': None} for p in self.prog_ids}
        for row in self.tap(q):
            s = out.setdefault(row['prog_id'], {'rows': 0, 'by_type': {}, 'last_mod': None})
            n = int(row['n'])
            s['rows'] += n
            s['by_type'][row['dp_type']] = s['by_type'].get(row['dp_type'], 0) + n
            if row['last_mod'] and (s['last_mod'] is None or row['last_mod'] > s['last_mod']):
                s['last_mod'] = row['last_mod']
        return out

    def night_objects(self, prog_id, night):
        a, b = (to_mjd(t) for t in night_bounds(night))
        q = ("SELECT DISTINCT object FROM dbo.raw WHERE prog_id = {} AND dp_id LIKE {} AND dp_type = 'OBJECT' "
             'AND mjd_obs >= {:.6f} AND mjd_obs < {:.6f}').format(adql_str(prog_id), adql_str('SPECU%'), a, b)
        return sorted({r['object'].strip() for r in self.tap(q) if r.get('object', '').strip()})

    def frames(self, prog_id, first_night, last_night):
        """Every frame row for the nights first..last inclusive, each with a 'night' key added."""
        a = to_mjd(night_bounds(first_night)[0])
        b = to_mjd(night_bounds(last_night)[1])
        q = ('SELECT dp_id, dp_type, object, origfile, mjd_obs, access_estsize, last_mod_date FROM dbo.raw '
             'WHERE prog_id = {} AND dp_id LIKE {} AND mjd_obs >= {:.6f} AND mjd_obs < {:.6f} ORDER BY mjd_obs'
             ).format(adql_str(prog_id), adql_str('SPECU%'), a, b)
        rows = self.tap(q)
        for r in rows:
            r['night'] = night_str(mjd_to_night(float(r['mjd_obs'])))
        return rows


def group_by_night(rows):
    out = {}
    for r in rows:
        out.setdefault(r['night'], []).append(r)
    return out
