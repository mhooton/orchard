"""
Priority job queue for the orchard pipeline on appct.

A SQLite queue, a dispatcher daemon that admits jobs by class, cores, load,
memory, disk utilisation and locks, an ESO watcher that queues each night as
soon as its frames are at ESO, and a daily look-back over recent nights.

Everything runs inside the orchard-server container on appct, with the
standard library only (the ESO side reuses download.request_eso). Entry point:

    python -m jobqueue <command>        # see `python -m jobqueue --help`

Design and operations: docs/job-queue.md
"""

CLASSES = ('P0', 'P1', 'P2', 'P3')
CLASS_RANK = {c: i for i, c in enumerate(CLASSES)}

ACTIVE_STATES = ('queued', 'running')
TERMINAL_STATES = ('done', 'failed', 'cancelled')
