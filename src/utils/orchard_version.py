"""
Which orchard made a data product: the version the pipeline stamps into what it writes.

Versions are git tags on main, vMAJOR.MINOR.PATCH (docs/versioning.md says when each part goes up). The
version of a piece of code is `git describe` of the commit it came from, without the leading v:

    4.1.0             the tagged commit itself
    4.1.0-3-gab12cd3  three commits after v4.1.0, at commit ab12cd3
    ab12cd3           no version tag reachable from the commit
    ...-dirty         a checkout with uncommitted changes (never production)

Production code is not a git checkout. Deploys stage `git archive <sha> src`, and git archive fills in the
$Format placeholders in ARCHIVE below (this file is export-subst in .gitattributes), so the version travels
with the code into every release and every scratch tree made the same way. In a checkout (development, the
tests) the placeholders are left alone and the version comes from `git describe` instead.

A code release (src/VERSION, written by `q release deploy`) adds its id, which names the exact tree a queue
job ran from, including the deployed-only files no commit holds.

Every FITS file a stage creates is stamped with ORCHVER (and ORCHREL from a release): stamp() or stamp_file().
Every later in-place edit of a file (by another stage, a later run or a tool) sets ORCHMOD, the version of the
latest edit, and adds one HISTORY card 'orchard <version> <script>: <what>': record_edit(), edit_file() or
edit_fitsio_hdu(). ORCHVER is never changed after creation.

    python <src>/utils/orchard_version.py            prints the version
    python <src>/utils/orchard_version.py --line     prints the one-line description the pipeline logs
"""
import functools
import os
import re
import subprocess
import sys

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Filled in by git archive. Outside an archive each value stays as written here, starting with a dollar sign,
# which no describe, sha or date does. (A git older than 2.35 leaves the describe's '%(' in place.)
ARCHIVE = {
    'describe': '$Format:%(describe:tags=true,match=v[0-9]*)$',
    'commit': '$Format:%H$',
    'date': '$Format:%cI$',
}

RELEASE_FILE = 'VERSION'           # src/VERSION in a release: "<release id> <sha>"

# set when a file is created: keyword, info() field, comment
KEYWORDS = (
    ('ORCHVER', 'version', 'orchard version that made this file'),
    ('ORCHREL', 'release', 'orchard code release that made this file'),
)
EDIT_KEYWORD = 'ORCHMOD'
EDIT_COMMENT = 'orchard version of the last in-place edit'
HISTORY_LEN = 72            # one HISTORY card; a longer text is split over several, which readers would not expect
# Text that code reads back out of HISTORY (repair_orientation, repair_astra_wcs, framepack): an edit card must not
# contain it, or a frame would look transformed or packed once more than it was.
RESERVED_HISTORY = ('horizontal mirroring', 'anticlockwise rotation', 'orient repair', 'framepack:')

_TAGGED = re.compile(r'^v(\d+\.\d+\.\d+\S*)$')
_SHA = re.compile(r'^[0-9a-f]{7,40}$')


def _filled(value):
    return bool(value) and not value.startswith('$') and '%(' not in value


def _strip_v(describe):
    """'v1.2.0-3-gab12cd3' -> '1.2.0-3-gab12cd3'; anything else (a bare sha from --always) as it is."""
    m = _TAGGED.match(describe)
    return m.group(1) if m else describe


def _git(src, *args):
    try:
        out = subprocess.run(('git',) + args, cwd=src, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else None


def _release(src):
    """(release id, sha) from src/VERSION, or (None, None) outside a release."""
    try:
        with open(os.path.join(src, RELEASE_FILE)) as f:
            first = f.readline().split()
    except OSError:
        return None, None
    if not first:
        return None, None
    return first[0], (first[1] if len(first) > 1 and _SHA.match(first[1]) else None)


def _checkout(src):
    """True when `src` is the src directory of an orchard checkout (and not, say, a tree inside some other repo)."""
    top = _git(src, 'rev-parse', '--show-toplevel')
    return bool(top) and os.path.realpath(os.path.join(top, 'src')) == os.path.realpath(src)


def compute(src=SRC, archive=None):
    """{'version', 'commit', 'date', 'release', 'source'} for the code in `src`. `source` says where the version
    came from: 'archive' (a git archive, so every deployed tree), 'git' (a checkout) or 'none'."""
    archive = ARCHIVE if archive is None else archive
    release, release_sha = _release(src)
    if _filled(archive.get('commit')):
        commit = archive['commit']
        describe = archive.get('describe')
        return {'version': _strip_v(describe) if _filled(describe) else commit[:7], 'commit': commit,
                'date': archive['date'] if _filled(archive.get('date')) else None,
                'release': release, 'source': 'archive'}
    if _checkout(src):
        describe = _git(src, 'describe', '--tags', '--match', 'v[0-9]*', '--dirty', '--always')
        commit = _git(src, 'rev-parse', 'HEAD')
        if describe and commit:
            return {'version': _strip_v(describe), 'commit': commit,
                    'date': _git(src, 'show', '-s', '--format=%cI', 'HEAD'), 'release': release, 'source': 'git'}
    return {'version': release_sha[:7] if release_sha else 'unknown', 'commit': release_sha, 'date': None,
            'release': release, 'source': 'none'}


@functools.lru_cache(maxsize=None)
def info():
    """compute() for this tree, worked out once per process."""
    return compute()


def version():
    """The version string, e.g. '1.2.0' or '1.2.0-3-gab12cd3'."""
    return info()['version']


def header_cards():
    """[(keyword, value, comment)] for a new file's header: ORCHVER, and ORCHREL when running from a release."""
    i = info()
    return [(key, i[field], comment) for key, field, comment in KEYWORDS if i[field]]


def stamp(header):
    """Mark an astropy Header as made by this orchard, and return it: set ORCHVER (and ORCHREL from a release) and
    drop what a header copied from another file brings along (its ORCHMOD, an ORCHREL this run did not set)."""
    cards = header_cards()
    for key, value, comment in cards:
        header[key] = (value, comment)
    written = {c[0] for c in cards}
    for key in [k for k, _, _ in KEYWORDS if k not in written] + [EDIT_KEYWORD]:
        header.remove(key, ignore_missing=True, remove_all=True)
    return header


def fitsio_records():
    """The creation keywords as fitsio header records, for FITS.write(..., header=...) or hdu.write_keys()."""
    return [{'name': key, 'value': value, 'comment': comment} for key, value, comment in header_cards()]


def edit_card(script, what=''):
    """The HISTORY card of an in-place edit: 'orchard <version> <script>: <what>', cut to one card."""
    text = 'orchard {} {}{}'.format(version(), script, ': ' + what if what else '')
    for reserved in RESERVED_HISTORY:
        if reserved in text:
            raise ValueError('an edit card must not contain {!r}: {!r}'.format(reserved, text))
    return text[:HISTORY_LEN]


def record_edit(header, script, what='', history=True):
    """Record an in-place edit of an existing file in its astropy Header, and return it: ORCHMOD, and the
    edit's HISTORY card unless the header already holds that card (a stage that edits a file twice in one run,
    as T8 and T9 both refit a frame's WCS, leaves one). `history=False` for a file whose HISTORY means something
    else (output.fits: the list of nights)."""
    header[EDIT_KEYWORD] = (version(), EDIT_COMMENT)
    if history:
        card = edit_card(script, what)
        if card not in [str(v) for v in header.get('HISTORY', [])]:
            header.add_history(card)
    return header


# The file-level helpers below never raise: a provenance stamp that fails is reported and the stage goes on, so it
# can never turn a good frame into a failed one (casu_solve marks a frame unsolved on any exception).

def _warn(action, target, err):
    print('orchard_version: could not {} {}: {}'.format(action, target, err), flush=True)


def stamp_file(path, hdu=0):
    """stamp() a file a tool has just created (CASU's catalogues and stacks copy their input frame's header), in
    place. The pixels are never read, so a scaled integer image stays as it is. Returns True if stamped."""
    from astropy.io import fits
    try:
        with fits.open(path, mode='update', do_not_scale_image_data=True) as hdul:
            stamp(hdul[hdu].header)
        return True
    except Exception as err:
        _warn('stamp', path, err)
        return False


def edit_file(path, script, what='', hdu=0, history=True):
    """record_edit() an existing file in place, header only (see stamp_file). Returns True if recorded."""
    from astropy.io import fits
    try:
        with fits.open(path, mode='update', do_not_scale_image_data=True) as hdul:
            record_edit(hdul[hdu].header, script, what, history=history)
        return True
    except Exception as err:
        _warn('record an edit of', path, err)
        return False


def edit_fitsio_hdu(hdu, script, what='', history=True):
    """record_edit() through an HDU of a file already open in fitsio 'rw' mode. Returns True if recorded."""
    try:
        hdu.write_key(EDIT_KEYWORD, version(), comment=EDIT_COMMENT)
        if history:
            card = edit_card(script, what)
            if card not in [str(r.get('value', '')).strip() for r in hdu.read_header().records()
                            if r.get('name') == 'HISTORY']:
                hdu.write_history(card)
        return True
    except Exception as err:
        _warn('record an edit of', getattr(hdu, '_filename', 'an HDU'), err)
        return False


def line():
    """One line for a log or a text product: 'orchard 1.2.0 (commit ab12cd3ef456, release 20261019T...)'."""
    i = info()
    extra = (['commit ' + i['commit'][:12]] if i['commit'] else []) + \
            (['release ' + i['release']] if i['release'] else [])
    return 'orchard {}{}'.format(i['version'], ' ({})'.format(', '.join(extra)) if extra else '')


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    print(line() if '--line' in argv else version())
    return 0


if __name__ == '__main__':
    sys.exit(main())
