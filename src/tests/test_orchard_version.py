"""
Tests for the orchard version (utils/orchard_version.py): where it comes from in a git archive, a checkout and
a release, and the keywords and comment line it puts into data products.

Run:  pytest src/tests/test_orchard_version.py -q
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys

import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import orchard_version as ov            # noqa: E402

MODULE = os.path.join(SRC, 'utils', 'orchard_version.py')
COMMIT = 'ab12cd3ef4567890ab12cd3ef4567890ab12cd3e'
# the module's ARCHIVE as committed (this file is not export-subst, so git archive leaves these alone)
PLACEHOLDERS = {'describe': '$Format:%(describe:tags=true,match=v[0-9]*)$', 'commit': '$Format:%H$',
                'date': '$Format:%cI$'}
ATTRIBUTES = 'src/utils/orchard_version.py export-subst\n'


def module_text(**values):
    """utils/orchard_version.py with its ARCHIVE set to `values`, the placeholders by default: the file as
    committed, whichever this tree holds (the suite also runs from deployed trees, where they are filled in)."""
    with open(MODULE) as f:
        text = f.read()
    for key, placeholder in PLACEHOLDERS.items():
        text, n = re.subn(r"(\n    '{}': ')[^'\n]*(',\n)".format(key),
                          lambda m: m.group(1) + values.get(key, placeholder) + m.group(2), text)
        assert n == 1, key
    return text


def git_ok():
    try:
        out = subprocess.run(['git', 'version'], stdout=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return False
    parts = out.stdout.split()[2].split('.')
    return (int(parts[0]), int(parts[1])) >= (2, 35)     # %(describe:tags=...) in git archive


needs_git = pytest.mark.skipif(not git_ok(), reason='needs git 2.35 or later')


def git(cwd, *args):
    env = dict(os.environ, GIT_AUTHOR_NAME='t', GIT_AUTHOR_EMAIL='t@t', GIT_COMMITTER_NAME='t',
               GIT_COMMITTER_EMAIL='t@t', GIT_CONFIG_NOSYSTEM='1', HOME=str(cwd))
    return subprocess.run(('git',) + args, cwd=str(cwd), env=env, check=True, stdout=subprocess.PIPE,
                          universal_newlines=True).stdout.strip()


def load(path, name='_orchard_version_under_test'):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def repo_with_module(tmp_path):
    """A small orchard repo: the real version module and .gitattributes, one commit."""
    repo = tmp_path / 'repo'
    (repo / 'src' / 'utils').mkdir(parents=True)
    (repo / 'src' / 'utils' / 'orchard_version.py').write_text(module_text())
    (repo / '.gitattributes').write_text(ATTRIBUTES)
    (repo / 'src' / 'x.py').write_text('X = 1\n')
    git(repo, 'init', '-q')
    git(repo, 'add', '-A')
    git(repo, 'commit', '-qm', 'first')
    return repo


# --------------------------------------------------------------------------- a git archive (every deployed tree)

@pytest.mark.parametrize('describe, expected', [
    ('v1.2.0', '1.2.0'),
    ('v1.2.0-3-gab12cd3', '1.2.0-3-gab12cd3'),
    ('', 'ab12cd3'),                                             # no version tag reachable from the commit
    ('%(describe:tags=true,match=v[0-9]*)', 'ab12cd3'),          # a git too old for %(describe) leaves it
])
def test_an_archived_tree_takes_its_version_from_the_placeholders(tmp_path, describe, expected):
    i = ov.compute(src=str(tmp_path), archive={'describe': describe, 'commit': COMMIT, 'date': '2026-10-19T12:00:00Z'})
    assert i['version'] == expected and i['commit'] == COMMIT and i['source'] == 'archive'
    assert i['date'] == '2026-10-19T12:00:00Z' and i['release'] is None


def test_a_release_adds_its_id(tmp_path):
    (tmp_path / 'VERSION').write_text('20261019T112705Z-ab12cd3 {}\n'.format(COMMIT))
    i = ov.compute(src=str(tmp_path), archive={'describe': 'v1.0.0', 'commit': COMMIT, 'date': ''})
    assert i['version'] == '1.0.0' and i['release'] == '20261019T112705Z-ab12cd3' and i['date'] is None


def test_a_tree_neither_archived_nor_a_checkout_says_unknown_or_its_release(tmp_path):
    unfilled = dict(PLACEHOLDERS)
    assert ov.compute(src=str(tmp_path), archive=unfilled)['version'] == 'unknown'
    (tmp_path / 'VERSION').write_text('20261017T085429Z-live live-snapshot\n')
    i = ov.compute(src=str(tmp_path), archive=unfilled)
    assert i['version'] == 'unknown' and i['release'] == '20261017T085429Z-live' and i['source'] == 'none'
    (tmp_path / 'VERSION').write_text('20261019T112705Z-ab12cd3 {}\n'.format(COMMIT))
    assert ov.compute(src=str(tmp_path), archive=unfilled)['version'] == 'ab12cd3'


def test_the_repository_holds_the_placeholders_and_git_archive_fills_them():
    """A committed filled-in copy would make every tree claim that commit; without the attribute git archive
    would leave every deployed tree 'unknown'."""
    if ov.info()['source'] == 'archive':
        pytest.skip('a deployed tree: git archive filled them in')
    assert ov.ARCHIVE == PLACEHOLDERS
    if ov.info()['source'] == 'git':
        attr = subprocess.run(['git', 'check-attr', 'export-subst', 'src/utils/orchard_version.py'],
                              cwd=os.path.dirname(SRC), stdout=subprocess.PIPE, universal_newlines=True).stdout
        assert attr.strip().endswith('export-subst: set')


@needs_git
def test_git_archive_fills_in_the_version(tmp_path):
    """The deploy route: `git archive <sha> src` of a tagged commit and of a later one."""
    repo = repo_with_module(tmp_path)
    git(repo, 'tag', '-a', 'v1.0.0', '-m', 'orchard 1.0.0')
    for name, expected in (('a', '1.0.0'), ('b', None)):
        if expected is None:
            (repo / 'src' / 'x.py').write_text('X = 2\n')
            git(repo, 'commit', '-qam', 'second')
            expected = '1.0.0-1-g' + git(repo, 'rev-parse', '--short=7', 'HEAD')
        out = tmp_path / name
        out.mkdir()
        tar = subprocess.run(['git', 'archive', 'HEAD', 'src'], cwd=str(repo), stdout=subprocess.PIPE, check=True)
        subprocess.run(['tar', '-x', '-C', str(out)], input=tar.stdout, check=True)
        module = load(out / 'src' / 'utils' / 'orchard_version.py', '_ov_' + name)
        i = module.compute()
        assert i['source'] == 'archive' and i['version'] == expected
        assert i['commit'] == git(repo, 'rev-parse', 'HEAD')


@needs_git
def test_a_lightweight_tag_counts_too(tmp_path):
    """GitHub's release page makes lightweight tags."""
    repo = repo_with_module(tmp_path)
    git(repo, 'tag', 'v2.0.0')
    tar = subprocess.run(['git', 'archive', 'HEAD', 'src'], cwd=str(repo), stdout=subprocess.PIPE, check=True)
    subprocess.run(['tar', '-x', '-C', str(tmp_path)], input=tar.stdout, check=True)
    assert load(tmp_path / 'src' / 'utils' / 'orchard_version.py').compute()['version'] == '2.0.0'


# --------------------------------------------------------------------------- a checkout (development, tests)

@needs_git
def test_a_checkout_uses_git_describe_and_marks_uncommitted_changes(tmp_path):
    repo = repo_with_module(tmp_path)
    module = load(repo / 'src' / 'utils' / 'orchard_version.py')
    sha = git(repo, 'rev-parse', 'HEAD')
    assert module.compute()['version'] == sha[:7] and module.compute()['source'] == 'git'
    git(repo, 'tag', '-a', 'v1.1.0', '-m', 'orchard 1.1.0')
    assert module.compute()['version'] == '1.1.0' and module.compute()['commit'] == sha
    (repo / 'src' / 'x.py').write_text('X = 3\n')
    assert module.compute()['version'] == '1.1.0-dirty'


@needs_git
def test_a_tree_inside_some_other_repository_is_not_a_checkout(tmp_path):
    outer = tmp_path / 'outer'
    (outer / 'releases' / 'r1' / 'src' / 'utils').mkdir(parents=True)
    (outer / 'releases' / 'r1' / 'src' / 'utils' / 'orchard_version.py').write_text(module_text())
    git(outer, 'init', '-q')
    git(outer, 'commit', '-q', '--allow-empty', '-m', 'unrelated')
    module = load(outer / 'releases' / 'r1' / 'src' / 'utils' / 'orchard_version.py')
    assert module.compute()['version'] == 'unknown'


# --------------------------------------------------------------------------- what products carry

@pytest.fixture
def released(monkeypatch):
    i = {'version': '1.2.0-3-gab12cd3', 'commit': COMMIT, 'date': None, 'release': '20261019T112705Z-ab12cd3',
         'source': 'archive'}
    monkeypatch.setattr(ov, 'info', lambda: i)
    return i


def test_stamp_sets_the_keywords_with_comments(released):
    h = ov.stamp(fits.Header())
    assert h['ORCHVER'] == '1.2.0-3-gab12cd3' and h['ORCHREL'] == '20261019T112705Z-ab12cd3'
    assert 'orchard' in h.comments['ORCHVER'] and 'release' in h.comments['ORCHREL']


def test_outside_a_release_there_is_no_orchrel(monkeypatch):
    monkeypatch.setattr(ov, 'info', lambda: {'version': '1.0.0', 'commit': COMMIT, 'date': None, 'release': None,
                                             'source': 'git'})
    h = ov.stamp(fits.Header())
    assert h['ORCHVER'] == '1.0.0' and 'ORCHREL' not in h


def test_fitsio_records_round_trip(released, tmp_path):
    fitsio = pytest.importorskip('fitsio')
    import numpy as np
    path = str(tmp_path / 'x.fits')
    with fitsio.FITS(path, 'rw', clobber=True) as f:
        f.write(np.zeros((2, 2)), header=ov.fitsio_records())
    h = fits.getheader(path)
    assert h['ORCHVER'] == '1.2.0-3-gab12cd3' and h['ORCHREL'] == '20261019T112705Z-ab12cd3'


def test_line(released, monkeypatch):
    assert ov.line() == 'orchard 1.2.0-3-gab12cd3 (commit ab12cd3ef456, release 20261019T112705Z-ab12cd3)'
    monkeypatch.setattr(ov, 'info', lambda: {'version': 'unknown', 'commit': None, 'date': None, 'release': None,
                                             'source': 'none'})
    assert ov.line() == 'orchard unknown'


def test_command_line(tmp_path):
    out = subprocess.run([sys.executable, MODULE, '--line'], stdout=subprocess.PIPE, universal_newlines=True,
                         check=True).stdout
    assert out.startswith('orchard ')
