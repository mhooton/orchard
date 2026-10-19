"""
Which orchard made a data product: the version the pipeline stamps into what it writes.

Versions are git tags on main, vMAJOR.MINOR.PATCH (docs/versioning.md says when each part goes up). The
version of a piece of code is `git describe` of the commit it came from, without the leading v:

    1.2.0             the tagged commit itself
    1.2.0-3-gab12cd3  three commits after v1.2.0, at commit ab12cd3
    ab12cd3           no version tag reachable from the commit
    ...-dirty         a checkout with uncommitted changes (never production)

Production code is not a git checkout. Deploys stage `git archive <sha> src`, and git archive fills in the
$Format placeholders in ARCHIVE below (this file is export-subst in .gitattributes), so the version travels
with the code into every release and every scratch tree made the same way. In a checkout (development, the
tests) the placeholders are left alone and the version comes from `git describe` instead.

A code release (src/VERSION, written by `q release deploy`) adds its id, which names the exact tree a queue
job ran from, including the deployed-only files no commit holds.

FITS products get ORCHVER (and ORCHREL from a release); text products a '# orchard ...' comment line.

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

# keyword, info() field, comment, comment in a file many runs add to (where the last run's version is the one kept)
KEYWORDS = (
    ('ORCHVER', 'version', 'orchard version (git describe)', 'orchard version of the last update'),
    ('ORCHREL', 'release', 'orchard code release that ran', 'orchard code release of the last update'),
)

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


def header_cards(updated=False):
    """[(keyword, value, comment)] for a FITS header: ORCHVER, and ORCHREL when running from a release.
    `updated`: for a file that gathers several runs' work, so the comments say these are the last run's."""
    i = info()
    return [(key, i[field], last if updated else comment) for key, field, comment, last in KEYWORDS if i[field]]


def stamp(header, updated=False):
    """Set the version keywords in an astropy Header and return it. Without a release ORCHREL is removed, so a
    header copied from another product never keeps a release that did not write this one."""
    cards = header_cards(updated)
    for key, value, comment in cards:
        header[key] = (value, comment)
    for key, _, _, _ in KEYWORDS:
        if key not in [c[0] for c in cards]:
            header.remove(key, ignore_missing=True, remove_all=True)
    return header


def fitsio_records():
    """The version keywords as fitsio header records, for FITS.write(..., header=...) or hdu.write_keys()."""
    return [{'name': key, 'value': value, 'comment': comment} for key, value, comment in header_cards()]


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
