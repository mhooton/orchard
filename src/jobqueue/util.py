"""Shared helpers: UTC timestamps, night arithmetic, atomic writes, filesystem checks.

All times in the queue are naive UTC datetimes, stored as ISO strings ending in Z.
A night D runs from D 15:00 UTC to D+1 15:00 UTC, the convention of SSO_download.py
and of the ESO archive queries.
"""
import datetime as dt
import json
import os

ISO = '%Y-%m-%dT%H:%M:%SZ'
MJD0 = dt.datetime(1858, 11, 17)
NIGHT_START = dt.timedelta(hours=15)

# Filesystems on which SQLite locking cannot be trusted.
NETWORK_FS = {'nfs', 'nfs4', 'cifs', 'smb3', 'smbfs', 'fuse.sshfs', '9p', 'afs', 'lustre', 'gpfs'}


def utcnow():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


def utc_from_ts(ts):
    return dt.datetime.fromtimestamp(ts, dt.timezone.utc).replace(tzinfo=None, microsecond=0)


def iso(t):
    return t.strftime(ISO) if t else None


def parse_iso(s):
    if not s:
        return None
    s = s.rstrip('Z').replace(' ', 'T')
    if '.' in s:
        s = s.split('.')[0]
    return dt.datetime.strptime(s, '%Y-%m-%dT%H:%M:%S')


def parse_night(s):
    return dt.datetime.strptime(str(s), '%Y%m%d').date()


def night_str(d):
    return d.strftime('%Y%m%d')


def last_night(now):
    """The most recent night whose frames can be complete.

    From 11:00 UTC on day X (after dawn in Chile and before ESO receives the
    frames) the last night is X-1; before 11:00 it is still X-2.
    """
    return (now - dt.timedelta(hours=11)).date() - dt.timedelta(days=1)


def night_bounds(night):
    """[start, end) of a night in UTC datetimes."""
    d = night if isinstance(night, dt.date) else parse_night(night)
    start = dt.datetime(d.year, d.month, d.day) + NIGHT_START
    return start, start + dt.timedelta(days=1)


def to_mjd(t):
    return (t - MJD0).total_seconds() / 86400.0


def from_mjd(mjd):
    return MJD0 + dt.timedelta(days=float(mjd))


def mjd_to_night(mjd):
    """The night (date) an exposure at this MJD belongs to."""
    return (from_mjd(mjd) - NIGHT_START).date()


def hhmm(s):
    h, m = s.split(':')
    return dt.time(int(h), int(m))


def in_daily_window(now, start, end):
    """True if now's time of day lies in [start, end); windows may wrap midnight."""
    t = now.time()
    a, b = hhmm(start), hhmm(end)
    if a <= b:
        return a <= t < b
    return t >= a or t < b


def write_json_atomic(path, obj):
    tmp = '{}.tmp{}'.format(path, os.getpid())
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=True, default=str)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def read_json(path, default=None):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return default


def mount_fstype(path, mounts_file='/proc/self/mounts'):
    """Filesystem type of the mount holding path, or None where /proc is absent (macOS)."""
    try:
        with open(mounts_file) as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    real = os.path.realpath(path)
    best, fstype = '', None
    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        mnt = parts[1].replace('\\040', ' ')
        if (real == mnt or real.startswith(mnt.rstrip('/') + '/') or mnt == '/') and len(mnt) >= len(best):
            best, fstype = mnt, parts[2]
    return fstype


def ensure_local_fs(path, mounts_file='/proc/self/mounts'):
    """Refuse to use a SQLite file on a network filesystem (appcs sees appct's disks over NFS)."""
    fstype = mount_fstype(path, mounts_file)
    if fstype in NETWORK_FS:
        raise RuntimeError('{} is on a {} filesystem; the queue database must only be opened on appct '
                           'itself, inside orchard-server (SQLite locking is unreliable over NFS)'
                           .format(path, fstype))
    return fstype
