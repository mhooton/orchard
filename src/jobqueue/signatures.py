"""Failure signatures from job logs, and the transient-or-deterministic call.

A signature is the last exception line in the log (else the last line that reports an error),
with paths, numbers and target names normalised so that the same fault on different nights groups
together. `q requeue --signature` matches on it after a fix.
"""
import re

TAIL_BYTES = 256 * 1024

EXC_LINE = re.compile(r'^\s*(?:[A-Za-z_][\w.]*\.)?[A-Z]\w*(?:Error|Exception|Exit|Interrupt|Warning|Failure)\b.*')
ERROR_LINE = re.compile(r'(?i)\b(error|failed|failure|fatal|killed|traceback|segmentation fault|aborted)\b')
PATH = re.compile(r"(/[\w.+@=:~-]+)+/?")
SPID = re.compile(r'\b(?:Sp|SP|sp)\d{4}[+-]\d{2,4}\w*')
DPID = re.compile(r'\bSPECU\d+\.\d{4}-\d\d-\d\dT[\d:.]+')
HEXNUM = re.compile(r'\b0x[0-9a-fA-F]+\b')
NUM = re.compile(r'(?<![A-Za-z_\d])\d+(\.\d+)?(?![A-Za-z_\d])')  # numbers, not digits inside names
SPACE = re.compile(r'\s+')


def read_tail(path, nbytes=TAIL_BYTES):
    try:
        with open(path, 'rb') as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - nbytes))
            return f.read().decode('utf-8', 'replace')
    except OSError:
        return ''


def normalise(line, limit=200):
    s = DPID.sub('<frame>', line.strip())
    s = PATH.sub('<path>', s)
    s = SPID.sub('<target>', s)
    s = HEXNUM.sub('#', s)
    s = NUM.sub('#', s)
    s = SPACE.sub(' ', s)
    return s[:limit]


def extract(text, exit_code=None):
    """Normalised signature of the failure recorded in a log tail."""
    lines = [l for l in text.splitlines() if l.strip()]
    for line in reversed(lines):
        if EXC_LINE.match(line) and not line.strip().startswith(('Warning', 'UserWarning', 'DeprecationWarning',
                                                                  'RuntimeWarning', 'FutureWarning')):
            return normalise(line)
    for line in reversed(lines):
        if ERROR_LINE.search(line):
            return normalise(line)
    if exit_code is not None and exit_code < 0:
        return 'killed by signal {}'.format(-exit_code)
    return 'exit {}'.format(exit_code) if exit_code is not None else 'no output'


def is_transient(signature, tail, exit_code, cfg):
    if exit_code in (cfg.get('transient_exit_codes') or []):
        return True
    patterns = cfg.get('transient_patterns') or []
    # Only the end of the log: earlier lines hold per-frame warnings (e.g. the plate solver's
    # [TIMEOUT] lines) that say nothing about why the job stopped.
    window = '\n'.join(tail.splitlines()[-30:])
    for p in patterns:
        if re.search(p, signature or '') or re.search(p, window):
            return True
    return False
