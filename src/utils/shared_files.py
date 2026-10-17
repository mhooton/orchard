"""Files that every pipeline run on a telescope reads, changes and writes back.

reports/Pipeline_Report and reports/QC (one table per telescope, updated by every stage),
output/StackImages/already_run.txt (the telescope's run targets, which rename later ones) and
output/<gaia_id>_<filter>_output.fits (one per star, gathering every night). Two runs on the same
telescope at once (two nights) must not interleave a read-change-write of one of these, or one run's
change is lost, and a reader must never see a half-written table.
"""
import contextlib
import fcntl
import os

from astropy.io import ascii


def acquire(fname):
    """Block until this process holds the exclusive lock on <fname>.lock; os.close() the result to release it."""
    fd = os.open(fname + '.lock', os.O_RDWR | os.O_CREAT, 0o664)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
    except BaseException:
        os.close(fd)
        raise
    return fd


@contextlib.contextmanager
def locked(fname):
    """Hold the lock on <fname>.lock for one read-change-write of fname."""
    fd = acquire(fname)
    try:
        yield
    finally:
        os.close(fd)


def write_table(table, fname):
    """Write a fixed-width table through a temporary file, so a reader sees the old table or the new one."""
    tmp = '{}.tmp{}'.format(fname, os.getpid())
    try:
        ascii.write(table, tmp, format='fixed_width', overwrite=True)
        os.replace(tmp, fname)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
