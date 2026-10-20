"""
Unit tests for reporting/completeness_check.py, the 05:00 cron job that writes <data root>/completeness_check.

The script ran under Pipeline_v2's Python 2.7 until October 2026.  Its Python 3 port had to keep two
behaviours that astropy 5 changes: values longer than a column are cut to its width (numpy did that
under Python 2; astropy 5 widens the column, changing every line), and each day's run replaces the
previous day's report (astropy 5's ascii.write refuses to overwrite by default).

Run:  pytest src/tests/test_completeness_check.py -q
"""
import datetime
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from reporting import completeness_check as cc      # noqa: E402

TELS = ['Europa', 'Io', 'Callisto', 'Ganymede', 'Artemis']


class _Date(datetime.date):
    @classmethod
    def today(cls):
        return cls(2017, 9, 1)


@pytest.fixture(autouse=True)
def short_survey(monkeypatch):
    """End the report on 2017-09-01: running it up to today takes ~20 s per call."""
    monkeypatch.setattr(datetime, 'date', _Date)


def _touch(path, lines=()):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.writelines(line + '\n' for line in lines)


def _data_root(tmp_path):
    root = str(tmp_path / 'SPECULOOSPipeline')
    for tel in TELS:
        os.makedirs(os.path.join(root, 'Observations', tel, 'images'))
        os.makedirs(os.path.join(root, 'PipelineOutput', 'v2', tel, 'output'))
    return root


def _night(root, tel, night, targets):
    """A fully processed night: two raw frames per target plus a bias, lists, output, light curves, PDF."""
    images = os.path.join(root, 'Observations', tel, 'images', night)
    out = os.path.join(root, 'PipelineOutput', 'v2', tel, 'output', night)
    _touch(os.path.join(images, 'bias', 'bias_%s.fits' % night))
    _touch(os.path.join(out, 'reduction', '1_bias.list'), ['/raw/%s/bias_%s.fits' % (night, night)])
    for target in targets:
        frames = ['%s_%d.fits' % (target, k) for k in (1, 2)]
        for frame in frames:
            _touch(os.path.join(images, target, frame))
        _touch(os.path.join(out, 'reduction', '1_image_%s.list' % target), frames)
        _touch(os.path.join(out, target, '%s_1_output.fits' % target))
        _touch(os.path.join(out, target, '%s_5_diff.fits' % target))
    _touch(os.path.join(root, 'PipelineOutput', 'v2', tel, 'reports', '%s_%s.pdf' % (tel, night)))


def _rows(path):
    with open(path) as f:
        lines = f.read().splitlines()
    return lines, [[cell.strip() for cell in line.split('|')[1:-1]] for line in lines]


def test_processed_night_row(tmp_path):
    root = _data_root(tmp_path)
    _night(root, 'Europa', '20170807', ['Sp0000-0000'])
    cc.create_file('report', root, cc.completeness_check(root, 'v2'), 'v2')

    lines, rows = _rows(os.path.join(root, 'report'))
    assert rows[0][4] == 'Pipeline (v2) Started'
    # the dashed test row, cut to the column widths as under Python 2 (16 and 11 dashes become 11 and 7)
    assert [len(cell) for cell in rows[1]] == [8, 9, 50, 12, 11, 7, 7, 10, 16, 11, 7, 30]
    assert set(''.join(rows[1])) == {'-'}
    row = [r for r in rows if r[:2] == ['20170807', 'Europa']]
    assert row == [['20170807', 'Europa', 'Sp0000-0000', '1', '1', '1', '1', '1', '1', '1', '1', '']]
    # nights with no images after a telescope's first night get a blank row
    assert ['20170808', 'Europa'] + [''] * 10 in rows
    assert len({len(line) for line in lines}) == 1


def test_long_values_are_cut_to_the_column_width(tmp_path):
    root = _data_root(tmp_path)
    targets = ['SPEC151047-28181%d' % k for k in range(4)]     # 4 x 17 characters + 3 spaces = 71
    _night(root, 'Europa', '20170807', targets)
    cc.create_file('report', root, cc.completeness_check(root, 'v2'), 'v2')

    lines, rows = _rows(os.path.join(root, 'report'))
    row = [r for r in rows if r[:2] == ['20170807', 'Europa']][0]
    assert len(row[2]) == 50 and row[2].split()[0] in targets     # in directory order, as in production
    assert row[8] == '1 1 1 1'
    assert len({len(line) for line in lines}) == 1
    assert len(lines[0].split('|')[3]) == 50 + 2


def test_second_run_replaces_the_report(tmp_path):
    root = _data_root(tmp_path)
    _night(root, 'Io', '20170807', ['Sp0000-0000'])
    cc.create_file('report', root, cc.completeness_check(root, 'v2'), 'v2')
    with open(os.path.join(root, 'report')) as f:
        first = f.read()

    _night(root, 'Io', '20170808', ['Sp1111-1111'])
    cc.create_file('report', root, cc.completeness_check(root, 'v2'), 'v2')
    with open(os.path.join(root, 'report')) as f:
        second = f.read()
    assert second != first
    assert '| 20170808 |        Io |' in second and 'Sp1111-1111' in second
