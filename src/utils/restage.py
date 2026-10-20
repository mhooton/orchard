#!/usr/bin/env python3
"""
restage.py
==========

Put back into v3 exactly what a rerun from a late stage reads, for a night whose products the promotion step has
already moved into v2, so that fixing T10 (the PWV correction, say) costs a rerun of T10 and T11 rather than one
from raw. docs/rerun-from-stage.md has the design; utils/stages.py says what each stage reads and writes.

For each target of the night, `stage`

  - from v2 (the night was promoted): makes a fresh v3 target directory and fills it with the inputs of the rerun,
    as hard links to v2's own files where no stage of the rerun edits them in place (every input from T10 on:
    the night's *_output.fits, the calibrated frames for T11's stack panel) and as copies otherwise. The night's
    v3 reduction gets the target's image list and the masters if it lacks them (it keeps its own if it has them:
    promotion leaves v3's reduction where it is). v2 is only read.
  - in v3 (never promoted, or held): sets aside what the stages of the rerun wrote there before, into
    PipelineOutput/superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-prerun/, so that nothing an earlier run left (a
    best-aperture file the new T10 would read, a curve the new one does not replace) can mix with the rerun's.
  - writes a record, v3/<TEL>/output/<NIGHT>/.restage/<TARGET>.json, of every input: where it came from, how it
    was placed, and its inode, size and mtime. The promotion step (promote.py) reads it: the files of the
    directory that are not unchanged inputs are what the rerun wrote (its manifest), S1 is judged on those only,
    and a night restaged from v2 is promoted as an overlay that replaces only those files in v2.

The record is written before anything is placed (state 'staging') and marked 'staged' at the end, so a restage cut
short is finished by running it again; one already done is checked and left alone.

Usage (inside orchard-server)
    python utils/restage.py stage --telescope T --night N --from 10 [--target X ...] [--apply]
    python utils/restage.py show  --telescope T --night N [--target X] [--manifest]
    python utils/restage.py clear --telescope T --night N [--target X] --reason TEXT [--apply]

Without --apply nothing changes. A process with --apply --telescope T --night N holds that night against the
queue (jobqueue.resources.parse_promote_argv), as a promotion does.
"""
import argparse
import datetime as dt
import errno
import fnmatch
import getpass
import hashlib
import json
import os
import shutil
import stat
import sys

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import stages  # noqa: E402

RUNNAME = '1'
RECORD_DIR = '.restage'          # v3/<TEL>/output/<NIGHT>/.restage/<TARGET>.json: not a target directory
VERSION = 1
SUPERSEDED = 'superseded'
EX_TEMPFAIL = 75
FITS = ('.fits', '.fit', '.fts')
CARDS = ('ORCHVER', 'ORCHREL', 'ORCHMOD')
FRAME_LISTS = ('1/*_processed.dat',)   # lists of calibrated frames by absolute path: copied, pointed at the copies


class RestageError(Exception):
    pass


def utcnow():
    return dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)


def iso(t):
    return t.strftime('%Y-%m-%dT%H:%M:%SZ')


def parse_iso(s):
    return dt.datetime.strptime(s, '%Y-%m-%dT%H:%M:%SZ')


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def output_dir(base, version, tel, night, *parts):
    return os.path.join(base, 'PipelineOutput', version, tel, 'output', night, *parts)


def record_path(base, tel, night, target):
    return output_dir(base, 'v3', tel, night, RECORD_DIR, target + '.json')


def superseded_dir(base, tel, night, name, stamp):
    return os.path.join(base, 'PipelineOutput', SUPERSEDED, tel, night, name, stamp)


def ident(path):
    """(inode, size, mtime_ns) of a path itself, or None."""
    try:
        st = os.lstat(path)
    except OSError:
        return None
    return [st.st_ino, st.st_size, st.st_mtime_ns]


def file_sha1(path):
    h = hashlib.sha1()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def walk_files(root):
    """Relative paths of the regular files (and links) under root, sorted."""
    out = []
    for dirpath, dirs, names in os.walk(root):
        dirs.sort()
        for name in sorted(names):
            out.append(os.path.relpath(os.path.join(dirpath, name), root))
    return out


def fits_cards(path):
    """ORCHVER, ORCHREL and ORCHMOD of a FITS file's primary header (the provenance orchard stamps on what it makes,
    utils/orchard_version.py), as far as they are there."""
    try:
        from astropy.io import fits
        h = fits.getheader(path, 0)
    except Exception:
        return {}
    return {k.lower(): str(h[k]) for k in CARDS if k in h}


# ---------------------------------------------------------------------------
# What a target needs, and where it is
# ---------------------------------------------------------------------------

def target_inputs(first, root, night, target):
    """[(rel, link)] of the files in a target directory that a rerun from `first` reads (made by earlier stages). A
    file a stage of the rerun writes is never one, whatever read pattern it matches: T9 makes the
    *_stack_catalogue_pm_* that T9's read of T8's *_stack_catalogue_* would otherwise take in, and writes it anew."""
    pats = [(p, ln) for scope, p, _, ln in stages.restage_inputs(first) if scope == stages.TARGET]
    out = []
    for rel in walk_files(root) if os.path.isdir(root) else []:
        if stages.is_output(rel, first, night=night, target=target):
            continue
        for p, ln in pats:
            if stages.matches(rel, p, night=night, target=target):
                out.append((rel, ln))
                break
    return out


def reduction_inputs(first, root, night, target):
    pats = [(p, ln) for scope, p, _, ln in stages.restage_inputs(first) if scope == stages.REDUCTION]
    out = []
    for name in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        if not os.path.isfile(os.path.join(root, name)):
            continue
        for p, ln in pats:
            if stages.matches(name, p, night=night, target=target):
                out.append((name, ln))
                break
    return out


def has_inputs(first, root, night, target):
    """Whether a target directory holds what the first stage of the rerun needs to run at all: for T8 a non-empty
    frame list, for T9 that and T8's stack catalogue; from T10 on the night's *_output.fits (T11 also needs T10's
    curves)."""
    names = os.listdir(root) if os.path.isdir(root) else []
    if first <= 9:
        lists = [n for n in (os.listdir(os.path.join(root, RUNNAME)) if os.path.isdir(os.path.join(root, RUNNAME))
                             else []) if n.endswith('_processed.dat')]
        if not any(os.path.getsize(os.path.join(root, RUNNAME, n)) for n in lists):
            return False
        return first <= 8 or \
            any('_stack_catalogue_' in n and '_stack_catalogue_pm_' not in n and n.endswith('.fits') for n in names)
    if not any(n.endswith('_{}_output.fits'.format(night)) for n in names):
        return False
    return first < 11 or any(n.endswith('_diff.fits') for n in names)


def night_targets(base, tel, night, first=10):
    """Targets of the night that a rerun from `first` can start from: v3 or v2 target directories with its inputs."""
    out = set()
    for version in ('v3', 'v2'):
        d = output_dir(base, version, tel, night)
        for name in os.listdir(d) if os.path.isdir(d) else []:
            if name == 'reduction' or name.startswith('.'):
                continue
            if has_inputs(first, os.path.join(d, name), night, name):
                out.add(name)
    return sorted(out)


def reduction_targets(base, tel, night, version='v2'):
    """The night's targets according to a reduction: one 1_image_<target>.list each (T1 writes them)."""
    d = output_dir(base, version, tel, night, 'reduction')
    pre, suf = RUNNAME + '_image_', '.list'
    return sorted(n[len(pre):-len(suf)] for n in (os.listdir(d) if os.path.isdir(d) else [])
                  if n.startswith(pre) and n.endswith(suf))


# ---------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------

def load_record(base, tel, night, target):
    path = record_path(base, tel, night, target)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        rec = json.load(f)
    rec['_path'] = path
    return rec


def save_record(rec):
    path = rec['_path']
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = '{}.tmp-{}'.format(path, os.getpid())
    with open(tmp, 'w') as f:
        json.dump({k: v for k, v in rec.items() if not k.startswith('_')}, f, indent=1, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def check(rec):
    """How the record's inputs are now. unchanged: as placed. edited: a copy the stage map says the rerun edits (a
    frame's WCS), so part of what it wrote. changed / missing: an input the rerun should only have read (a linked
    input changed in place means v2's file changed too). src_changed: the v2 file it came from is no longer the one
    restaged (v2 changed since). reduction_changed: the night's v3 reduction files it relies on."""
    out = {'unchanged': [], 'changed': [], 'edited': [], 'missing': [], 'src_changed': [], 'reduction_changed': []}
    for e in rec.get('files', []):
        now = ident(os.path.join(rec['dest'], e['rel']))
        if now is None:
            out['missing'].append(e['rel'])
        elif now == [e['ino'], e['size'], e['mtime_ns']]:
            out['unchanged'].append(e['rel'])
        elif e.get('edited'):
            out['edited'].append(e['rel'])          # a copy the rerun edits by design (a frame's WCS): its output
        else:
            out['changed'].append(e['rel'])
        was = [e['src_ino'], e.get('src_size', e['size']), e.get('src_mtime_ns', e['mtime_ns'])]
        if e['how'] in ('link', 'copy', 'rewritten') and ident(e['src']) != was:      # not 'made' or 'in-place'
            out['src_changed'].append(e['rel'])
    red = output_dir(rec['base'], 'v3', rec['telescope'], rec['night'], 'reduction')
    for e in rec.get('reduction', []):
        if ident(os.path.join(red, e['rel'])) != [e['ino'], e['size'], e['mtime_ns']]:
            out['reduction_changed'].append(e['rel'])
    return out


def manifest(rec, run_start=None, sha=True):
    """The rerun's manifest: every file of the restaged target directory that is not an unchanged input, i.e. what
    the rerun wrote, with size, mtime, ctime, sha1 and, for FITS files, the ORCHVER/ORCHREL/ORCHMOD orchard stamped
    on them. 'early': those of them whose inode last changed before run_start (UTC datetime), which the run cannot
    have written."""
    st = check(rec)
    inputs = set(st['unchanged'])
    written, early = [], []
    limit = None if run_start is None else (run_start - dt.datetime(1970, 1, 1)).total_seconds() - 1.0
    for rel in walk_files(rec['dest']) if os.path.isdir(rec['dest']) else []:
        if rel in inputs:
            continue
        p = os.path.join(rec['dest'], rel)
        s = os.lstat(p)
        if not stat.S_ISREG(s.st_mode):
            continue
        e = {'rel': rel, 'size': s.st_size, 'mtime_ns': s.st_mtime_ns, 'ctime': s.st_ctime, 'ino': s.st_ino,
             'nlink': s.st_nlink}
        if sha:
            e['sha1'] = file_sha1(p)
        if rel.lower().endswith(FITS):
            e.update(fits_cards(p))
        written.append(e)
        if limit is not None and s.st_ctime < limit:
            early.append(rel)
    return dict(st, inputs=sorted(inputs), written=written, early=early)


# ---------------------------------------------------------------------------
# Restaging
# ---------------------------------------------------------------------------

def place(src, dst, link):
    """Hard-link (or, across filesystems or when the rerun edits it, copy) src to dst. Returns 'link' or 'copy'.
    Recognises its own earlier result, so a restage cut short can be run again."""
    if os.path.lexists(dst):
        if link and os.path.samefile(src, dst):
            return 'link'
        s, d = os.lstat(src), os.lstat(dst)
        if not link and (s.st_size, s.st_mtime_ns) == (d.st_size, d.st_mtime_ns) and s.st_ino != d.st_ino:
            return 'copy'
        raise RestageError('{} exists and is not the restaged {}'.format(dst, src))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if link:
        try:
            os.link(src, dst)
            return 'link'
        except OSError as e:
            if e.errno != errno.EXDEV:
                raise
    tmp = os.path.join(os.path.dirname(dst), '.{}.restage-{}'.format(os.path.basename(dst), os.getpid()))
    shutil.copy2(src, tmp)
    os.rename(tmp, dst)
    return 'copy'


def place_list(src, dst, src_dir, dst_dir):
    """Copy a list of frames (absolute paths), pointing each line at the restaged copy of that frame when there is one.
    In production the paths are already the v3 ones restaging recreates, so the copy is the same; a scratch tree, or
    a night whose lists name another place, needs them rewritten. Returns 'rewritten' or 'copy'."""
    with open(src) as f:
        lines = f.read().splitlines()
    out = []
    for line in lines:
        name = os.path.basename(line.strip())
        if line.strip() and os.path.lexists(os.path.join(dst_dir, name)):
            out.append(os.path.join(dst_dir, name))
        else:
            out.append(line)
    text = '\n'.join(out) + ('\n' if lines else '')
    if os.path.lexists(dst):
        with open(dst) as f:
            if f.read() == text:
                return 'rewritten' if out != lines else 'copy'
        raise RestageError('{} exists and is not the restaged {}'.format(dst, src))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = os.path.join(os.path.dirname(dst), '.{}.restage-{}'.format(os.path.basename(dst), os.getpid()))
    with open(tmp, 'w') as f:
        f.write(text)
    shutil.copystat(src, tmp)
    os.rename(tmp, dst)
    return 'rewritten' if out != lines else 'copy'


def zlp_catalogue(root, filt):
    """The stack catalogue ZLP_pipeline.sh's choose_stack_catalogue would give T9 for a filter in a target directory, or
    None: the one .active_stack_catalogue_<filter> names, else (with .expected_stack_catalogue) a Gaia-named one that
    is not the expected one, else the expected one, else (without either record) any Gaia-named one."""
    active = os.path.join(root, '.active_stack_catalogue_' + filt)
    if os.path.isfile(active):
        name = open(active).read().strip()
        if name and os.path.isfile(os.path.join(root, name)):
            return name
    expected = None
    if os.path.isfile(os.path.join(root, '.expected_stack_catalogue')):
        expected = open(os.path.join(root, '.expected_stack_catalogue')).read().strip() or None
    for n in sorted(os.listdir(root)):
        if n[:1].isdigit() and fnmatch.fnmatchcase(n, '[0-9]*_stack_catalogue_[!p]*.fits') and n != expected:
            return n
    return expected


def filters_of(root):
    run = os.path.join(root, RUNNAME)
    return sorted(n[:-len('_processed.dat')] for n in (os.listdir(run) if os.path.isdir(run) else [])
                  if n.endswith('_processed.dat'))


def catalogue_records(root, night, target):
    """For a night reduced before T8 recorded the catalogue it validated (target-named catalogues, before 2026-10):
    {filter: name} of the .active_stack_catalogue_<filter> records T9 needs, naming the only stack catalogue of that
    filter in the directory, the one the night's products were made with. Raises when there is more than one."""
    out = {}
    for filt in filters_of(root):
        if zlp_catalogue(root, filt):
            continue
        cands = sorted(n for n in os.listdir(root) if n.endswith('_stack_catalogue_{}.fits'.format(filt)) and
                       '_stack_catalogue_pm_' not in n)
        if len(cands) != 1:
            raise RestageError('T9 would find no stack catalogue for {} {} and there is {} to choose from: {}'.format(
                target, filt, 'none' if not cands else 'more than one', ', '.join(cands) or '-'))
        out[filt] = cands[0]
    return out


def has_target_roles(path):
    """Whether a stack catalogue has the TARGET_ROLE column condense reads (T8 writes it since 2026-03-31)."""
    try:
        from astropy.io import fits
        with fits.open(path) as h:
            for hdu in h[1:]:
                if hdu.name.upper() == 'GAIA_CROSSMATCH':
                    return 'TARGET_ROLE' in [c.upper() for c in hdu.columns.names]
    except Exception:
        return False
    return False


def entry(rel, how, src, dst, edited=False):
    d, s = os.lstat(dst), os.lstat(src)
    return {'rel': rel, 'how': how, 'src': src, 'size': d.st_size, 'mtime_ns': d.st_mtime_ns, 'ino': d.st_ino,
            'src_ino': s.st_ino, 'src_size': s.st_size, 'src_mtime_ns': s.st_mtime_ns, 'edited': bool(edited)}


class Restager(object):
    def __init__(self, base, apply=False, out=print, now_fn=utcnow, by=None, job=None):
        self.base = base
        self.apply = apply
        self.out = out
        self.now = now_fn
        self.by = by
        self.job = job

    def plan(self, tel, night, target, first):
        """What restaging one target would do: {'action': 'restage'|'in-place'|'done'|'resume', ...}, or raise
        RestageError when it cannot be done."""
        v3 = output_dir(self.base, 'v3', tel, night, target)
        v2 = output_dir(self.base, 'v2', tel, night, target)
        rec = load_record(self.base, tel, night, target)
        if rec is not None:
            if rec['from_stage'] != first:
                raise RestageError('{} was restaged for a rerun from T{} ({}); clear it first'.format(
                    target, rec['from_stage'], rec['_path']))
            if rec.get('state') == 'staged':
                st = check(rec)
                bad = st['changed'] + st['missing'] + st['reduction_changed']
                if bad:
                    raise RestageError('{} is restaged, but {} of its inputs changed since ({}); clear it and restage '
                                       'again'.format(target, len(bad), ', '.join(bad[:3])))
                return {'action': 'done', 'record': rec}
            return {'action': 'resume', 'record': rec, 'source': rec['source']}
        if os.path.isdir(v3):
            if not has_inputs(first, v3, night, target):
                raise RestageError('v3 already has a directory for {} without the inputs of a rerun from T{} ({}); '
                                   'it is not this tool\'s to move'.format(target, first, v3))
            return {'action': 'in-place', 'source': 'v3', 'src_dir': v3}
        if not has_inputs(first, v2, night, target):
            raise RestageError('neither v3 nor v2 has the inputs of a rerun from T{} for {} {} {}'.format(
                first, tel, night, target))
        return {'action': 'restage', 'source': 'v2', 'src_dir': v2}

    def target(self, tel, night, target, first):
        """Restage one target. Returns a one-line summary."""
        p = self.plan(tel, night, target, first)
        if p['action'] == 'done':
            return 'already restaged for a rerun from T{} at {} ({} inputs unchanged)'.format(
                first, p['record']['made_at'], len(p['record']['files']))
        source = p['source']
        v3 = output_dir(self.base, 'v3', tel, night, target)
        src_dir = p.get('src_dir') or p['record']['source_dir']
        files = target_inputs(first, src_dir, night, target)
        red_v3 = output_dir(self.base, 'v3', tel, night, 'reduction')
        red_v2 = output_dir(self.base, 'v2', tel, night, 'reduction')
        red = []
        for name, ln in reduction_inputs(first, red_v2, night, target):
            red.append((name, ln, 'present' if os.path.isfile(os.path.join(red_v3, name)) else 'place'))
        have = {name for name, _, _ in red}
        red += [(name, ln, 'present') for name, ln in reduction_inputs(first, red_v3, night, target)
                if name not in have]
        if first <= 10 and not any(n == '{}_image_{}.list'.format(RUNNAME, target) for n, _, _ in red):
            raise RestageError('no {}_image_{}.list in the v3 or v2 reduction of the night: T10 needs it'.format(
                RUNNAME, target))
        aside = []
        if source == 'v3':
            aside = [rel for rel in walk_files(v3) if stages.is_output(rel, first, night=night, target=target)]
        records = catalogue_records(src_dir, night, target) if first == 9 else {}     # T8 makes its own
        if first == 9:
            for filt in sorted(set(filters_of(src_dir)) - set(records)) + sorted(records):
                name = records.get(filt) or zlp_catalogue(src_dir, filt)
                if name and not has_target_roles(os.path.join(src_dir, name)):
                    raise RestageError('{} predates the TARGET_ROLE column condense reads (T8 writes it since '
                                       '2026-03-31): rerun this night from T8'.format(name))
        links = sum(1 for _, ln in files if ln)
        copied = sum(os.path.getsize(os.path.join(src_dir, rel)) for rel, ln in files if not ln) if source == 'v2' else 0
        summary = '{} for a rerun from T{}: {} inputs ({} as hard links{}), reduction {} ({} placed){}'.format(
            'restaged from v2' if source == 'v2' else 'restaged in place in v3', first, len(files), links,
            ', {} copied: {:.2f} GB'.format(len(files) - links, copied / 1e9) if copied else '',
            ', '.join(sorted(n for n, _, _ in red)), sum(1 for _, _, how in red if how == 'place'),
            '; {} earlier outputs set aside'.format(len(aside)) if aside else '') + ''.join(
            '; T9 told to use {} for {}'.format(name, filt) for filt, name in sorted(records.items()))
        if not self.apply:
            return 'would be ' + summary
        stamp = self.now().strftime('%Y%m%dT%H%M%SZ')
        rec = p.get('record') or {
            '_path': record_path(self.base, tel, night, target), 'version': VERSION, 'stage_map': stages.VERSION,
            'base': self.base, 'telescope': tel, 'night': night, 'target': target, 'from_stage': first,
            'stages': stages.rerun_stages(first), 'source': source, 'source_dir': src_dir, 'dest': v3,
            'state': 'staging', 'made_at': iso(self.now()), 'by': self.by, 'job': self.job, 'code': code_line(),
            'files': [], 'reduction': [], 'set_aside': None}
        save_record(rec)                 # first: a restage cut short is recognised and finished, never mistaken
        if source == 'v3':
            rec['set_aside'] = rec.get('set_aside') or {
                'dir': superseded_dir(self.base, tel, night, target, stamp + '-prerun'), 'files': aside}
            save_record(rec)
            set_aside(v3, rec['set_aside'])
            # the inputs are the products of the earlier stages; anything no stage writes (clutter) is not one, and
            # the promotion's S1 sees it as it would in any run
            rec['files'] = [entry(rel, 'in-place', os.path.join(v3, rel), os.path.join(v3, rel),
                                  edited=stages.edited_by_rerun(rel, first, night=night, target=target))
                            for rel in walk_files(v3) if stages.made_before(rel, first, night=night, target=target)]
        else:
            done = {e['rel'] for e in rec['files']}
            lists = [(rel, ln) for rel, ln in files if any(stages.matches(rel, p) for p in FRAME_LISTS)]
            for rel, ln in [x for x in files if x not in lists] + lists:     # frames first: the lists point at them
                if rel in done:
                    continue
                src, dst = os.path.join(src_dir, rel), os.path.join(v3, rel)
                if (rel, ln) in lists:
                    how = place_list(src, dst, os.path.dirname(src), os.path.dirname(dst))
                else:
                    how = place(src, dst, ln)
                rec['files'].append(entry(rel, how, src, dst, edited=not ln))
                if len(rec['files']) % 500 == 0:
                    save_record(rec)                                          # a long copy can be resumed
        for filt, name in sorted(records.items()):
            rel = '.active_stack_catalogue_' + filt
            if rel in {e['rel'] for e in rec['files']}:
                continue
            with open(os.path.join(v3, rel), 'w') as f:                      # what T8 writes since 2026-10-17
                f.write(name + '\n')
            path = os.path.join(v3, rel)
            rec['files'].append(entry(rel, 'made', path, path))
        rec['reduction'] = []
        for name, ln, how in red:
            dst = os.path.join(red_v3, name)
            if how == 'place':
                src = os.path.join(red_v2, name)
                how = place(src, dst, ln)
            else:
                src = dst
            rec['reduction'].append(entry(name, how, src, dst))
        rec['state'] = 'staged'
        rec['staged_at'] = iso(self.now())
        save_record(rec)
        return summary

    def night(self, tel, night, first, targets=None):
        """Restage the given targets (default: every target of the night with inputs). Returns (rc, lines)."""
        if first not in stages.RESTAGE:
            raise RestageError('a rerun from T{} cannot be restaged: only from {} (utils/stages.py restage_inputs '
                               'says why)'.format(first, ', '.join('T{}'.format(n) for n in stages.RESTAGE)))
        explicit = bool(targets)
        targets = list(targets or night_targets(self.base, tel, night, first))
        if not targets:
            return 1, ['{} {}: no target of the night has the inputs of a rerun from T{} in v3 or v2'.format(
                tel, night, first)]
        rc, lines = 0, []
        for t in targets:
            try:
                lines.append('{} {} {}: {}'.format(tel, night, t, self.target(tel, night, t, first)))
            except (RestageError, OSError) as e:
                rc = 1
                lines.append('{} {} {}: {}{}'.format(tel, night, t, 'NOT restaged: ' if explicit else 'skipped: ', e))
        return rc, lines


def set_aside(root, aside):
    """Move a v3 target directory's earlier outputs of the rerun's stages into the superseded tree, keeping their
    relative paths (renames on one filesystem); recognises files moved before an interruption."""
    for rel in aside['files']:
        src, dst = os.path.join(root, rel), os.path.join(aside['dir'], rel)
        if not os.path.lexists(src) and os.path.lexists(dst):
            continue
        if os.path.lexists(dst):
            raise RestageError('{} is set aside already and {} is still there'.format(dst, src))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.rename(src, dst)
    with open(os.path.join(aside['dir'], 'superseded.json'), 'w') as f:
        json.dump({'set_aside_at': iso(utcnow()), 'from': root, 'files': len(aside['files']),
                   'why': 'outputs of an earlier run, set aside before a rerun of the stages that write them'},
                  f, indent=1, sort_keys=True)


def code_line():
    try:
        from utils import orchard_version
        return orchard_version.line()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Undoing a restage that will not be promoted
# ---------------------------------------------------------------------------

def clear(base, tel, night, target, reason, apply=False, now_fn=utcnow):
    """Retire a restage: a v2 one's directory (inputs and whatever the rerun wrote) goes to
    superseded/<TEL>/<NIGHT>/<TARGET>/<stamp>-restage/ with the record, its links to v2's files dropped; a v3 one
    keeps its directory and gets its record moved beside the outputs it set aside. Returns a summary."""
    from utils.promote import claim_and_rename, thin, write_note
    rec = load_record(base, tel, night, target)
    if rec is None:
        raise RestageError('{} {} {} is not restaged'.format(tel, night, target))
    stamp = now_fn().strftime('%Y%m%dT%H%M%SZ')
    if rec['source'] == 'v3':
        dst = (rec.get('set_aside') or {}).get('dir') or superseded_dir(base, tel, night, target, stamp + '-prerun')
        if not apply:
            return 'would move the record to {} (the directory stays in v3; what was set aside stays there)'.format(dst)
        os.makedirs(dst, exist_ok=True)
        os.rename(rec['_path'], os.path.join(dst, 'restage.json'))
        try:
            os.rmdir(os.path.dirname(rec['_path']))
        except OSError:
            pass
        return 'record moved to {}'.format(dst)
    dst = superseded_dir(base, tel, night, target, stamp + '-restage')
    if not apply:
        return 'would move {} to {} and drop its links to v2'.format(rec['dest'], dst)
    thinned = None
    if os.path.isdir(rec['dest']):
        claim_and_rename(rec['dest'], dst)
        dropped = drop_links(dst, rec)
        thinned = thin(dst)['deleted']          # copies of v2's frames (a restage from T9 or T8), and what a rerun made
    else:
        os.makedirs(dst, exist_ok=True)
        dropped = 0
    os.rename(rec['_path'], os.path.join(dst, 'restage.json'))
    try:
        os.rmdir(os.path.dirname(rec['_path']))         # .restage/, if no other target of the night is restaged
    except OSError:
        pass
    write_note(dst, {'cleared_at': iso(now_fn()), 'reason': reason, 'was': rec['dest'], 'links_dropped': dropped,
                     'frames_thinned': thinned})
    return 'moved to {} ({} links to v2 dropped, {} files thinned)'.format(dst, dropped, thinned or 0)


def drop_links(root, rec):
    """Unlink the restaged hard links under root that are still v2's own file (nothing is lost: v2 keeps it)."""
    n = 0
    for e in rec.get('files', []):
        p = os.path.join(root, e['rel'])
        if e['how'] == 'link' and os.path.lexists(p) and os.path.lexists(e['src']) and os.path.samefile(p, e['src']):
            os.unlink(p)
            n += 1
    return n


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _defaults(a):
    from jobqueue.config import load_config, queue_paths
    cfg = load_config(a.config)
    a.base = a.base or cfg['basedir']
    a.queue_db = a.queue_db or queue_paths(cfg)['db']
    return cfg


def _busy(a):
    from utils.promote import night_busy, open_queue
    own = os.environ.get('JOBQUEUE_JOB_ID')
    conn = open_queue(a.queue_db)
    try:
        return night_busy(conn, a.telescope, a.night, int(own) if own and own.isdigit() else None)
    finally:
        if conn is not None:
            conn.close()


def _first(x):
    return int(str(x).upper().lstrip('T'))


def cmd_stage(a):
    if a.apply:
        busy = _busy(a)
        if busy:
            print('{} {}: not now: {}'.format(a.telescope, a.night, busy))
            return EX_TEMPFAIL
    job = os.environ.get('JOBQUEUE_JOB_ID')
    r = Restager(a.base, apply=a.apply, by=a.by or _user(), job=int(job) if job and job.isdigit() else None)
    try:
        rc, lines = r.night(a.telescope, a.night, a.first, a.target)
    except RestageError as e:
        print(e)
        return 2
    print('\n'.join(lines))
    if not a.apply:
        print('dry run: nothing changed; add --apply to restage')
    return rc


def restaged_targets(base, tel, night):
    d = output_dir(base, 'v3', tel, night, RECORD_DIR)
    return sorted(n[:-len('.json')] for n in os.listdir(d) if n.endswith('.json')) if os.path.isdir(d) else []


def cmd_show(a):
    targets = a.target or restaged_targets(a.base, a.telescope, a.night)
    if not targets:
        print('{} {}: nothing restaged'.format(a.telescope, a.night))
        return 0
    for t in targets:
        rec = load_record(a.base, a.telescope, a.night, t)
        if rec is None:
            print('{}: not restaged'.format(t))
            continue
        st = check(rec)
        print('{} {} {}: {} for a rerun from T{} ({}), from {} at {}{}'.format(
            a.telescope, a.night, t, rec['state'], rec['from_stage'], ' '.join('T{}'.format(n) for n in rec['stages']),
            rec['source'], rec['made_at'], ', job {}'.format(rec['job']) if rec.get('job') else ''))
        print('  inputs: {} unchanged, {} changed, {} missing; v2 sources changed since: {}; reduction changed: '
              '{}'.format(len(st['unchanged']), len(st['changed']), len(st['missing']), len(st['src_changed']),
                          len(st['reduction_changed'])))
        if a.manifest:
            m = manifest(rec, sha=False)
            print('  written by the rerun: {} files'.format(len(m['written'])))
            for e in m['written']:
                print('    {}  {} bytes{}'.format(e['rel'], e['size'], ''.join(
                    ' {}={}'.format(k.upper(), e[k]) for k in ('orchver', 'orchrel') if e.get(k))))
    return 0


def cmd_clear(a):
    if a.apply:
        busy = _busy(a)
        if busy:
            print('{} {}: not now: {}'.format(a.telescope, a.night, busy))
            return EX_TEMPFAIL
    targets = a.target or restaged_targets(a.base, a.telescope, a.night)
    rc = 0
    for t in targets:
        try:
            print('{} {} {}: {}'.format(a.telescope, a.night, t, clear(a.base, a.telescope, a.night, t, a.reason,
                                                                     apply=a.apply)))
        except (RestageError, OSError) as e:
            print('{} {} {}: {}'.format(a.telescope, a.night, t, e))
            rc = 1
    return rc


def _user():
    try:
        return getpass.getuser()
    except Exception:
        return 'uid{}'.format(os.getuid())


def build_parser():
    p = argparse.ArgumentParser(prog='restage.py', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    def common(s):
        s.add_argument('--config', help='queue config (default $ORCHARD_QUEUE_CONFIG, then <queue_root>/config.json)')
        s.add_argument('--base', help='data root (default: the queue config basedir)')
        s.add_argument('--queue-db', help='the queue database, read-only (default <queue_root>/queue.sqlite)')
        s.add_argument('--telescope', required=True)
        s.add_argument('--night', required=True)
        s.add_argument('--target', action='append', help='only this target (repeatable; default: every target of the '
                                                         'night with the inputs)')

    s = sub.add_parser('stage', help='put back into v3 what a rerun from a stage reads')
    common(s)
    s.add_argument('--from', dest='first', type=_first, required=True, help='the first stage of the rerun: T10 or T11')
    s.add_argument('--apply', action='store_true')
    s.add_argument('--by', help='who restages (default: the user)')
    s.set_defaults(fn=cmd_stage)

    s = sub.add_parser('show', help='a night\'s restage records and how their inputs are now')
    common(s)
    s.add_argument('--manifest', action='store_true', help='also list what the rerun wrote')
    s.set_defaults(fn=cmd_show)

    s = sub.add_parser('clear', help='retire a restage that will not be promoted')
    common(s)
    s.add_argument('--reason', required=True)
    s.add_argument('--apply', action='store_true')
    s.set_defaults(fn=cmd_clear)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    _defaults(a)
    return a.fn(a)


if __name__ == '__main__':
    sys.exit(main())
