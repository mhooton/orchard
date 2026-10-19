"""
Tests for the stack catalogue T9 and condense use (ZLP_pipeline.sh choose_stack_catalogue, and the
.active_stack_catalogue_<filter> record catalogue_fov.py writes).

On 2026-10-17 a rerun's T8 validated tonight's stack, but T9 and condense took a backup catalogue an earlier,
failed run had restored into the same directory (the old rule prefers any catalogue that is not tonight's):
condense crashed on 2024 catalogues without gaia_dr3_id (jobs 1502, 1625, 1627, 1649) and one night's photometry
ran on another night's catalogue (job 1696).

Run:  pytest src/tests/test_active_stack_catalogue.py -q
"""
import os
import re
import subprocess
import sys

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

TONIGHT = '3357958852543597056_Callisto_andor_20240126_stack_catalogue_I+z.fits'
STALE = '3357958852543597056_Callisto_andor_20240127_stack_catalogue_I+z.fits'
TONIGHT_ZYJ = '3357958852543597056_Callisto_spirit_20240126_stack_catalogue_zYJ.fits'


def shell_function():
    text = open(os.path.join(SRC, 'main', 'ZLP_pipeline.sh')).read()
    m = re.search(r'^choose_stack_catalogue\(\) \{\n.*?^\}\n', text, re.S | re.M)
    assert m, 'choose_stack_catalogue not found in ZLP_pipeline.sh'
    return m.group(0)


def choose(directory, filt='I+z'):
    script = 'set -e\nOUTEXT=fits\n' + shell_function() + 'CATFILE=$(choose_stack_catalogue "%s" "%s")\necho "[$CATFILE]"\n' % (
        directory, filt)
    out = subprocess.run(['bash', '-c', script], capture_output=True, text=True, check=True).stdout.strip()
    return os.path.basename(out.strip('[]')) if out != '[]' else ''


def make(tmp_path, files, expected=TONIGHT, active=None):
    for f in files:
        (tmp_path / f).write_text('x')
    if expected:
        (tmp_path / '.expected_stack_catalogue').write_text(expected + '\n')
    for filt, name in (active or {}).items():
        (tmp_path / ('.active_stack_catalogue_' + filt)).write_text(name + '\n')
    return str(tmp_path)


def test_stale_backup_is_not_taken_when_t8_validated_tonight(tmp_path):
    d = make(tmp_path, [TONIGHT, STALE], active={'I+z': TONIGHT})
    assert choose(d) == TONIGHT


def test_a_backup_restored_tonight_is_used(tmp_path):
    d = make(tmp_path, [TONIGHT, STALE], active={'I+z': STALE})
    assert choose(d) == STALE


def test_directories_without_the_record_keep_the_old_rule(tmp_path):
    assert choose(make(tmp_path, [TONIGHT, STALE])) == STALE
    (tmp_path / STALE).unlink()
    assert choose(str(tmp_path)) == TONIGHT


def test_a_record_naming_a_missing_file_falls_back(tmp_path):
    d = make(tmp_path, [TONIGHT], active={'I+z': STALE})
    assert choose(d) == TONIGHT


def test_records_are_per_filter(tmp_path):
    d = make(tmp_path, [TONIGHT, STALE, TONIGHT_ZYJ], active={'zYJ': TONIGHT_ZYJ, 'I+z': TONIGHT})
    assert choose(d, 'I+z') == TONIGHT
    assert choose(d, 'zYJ') == TONIGHT_ZYJ


def test_nothing_to_choose_prints_nothing_and_succeeds(tmp_path):
    assert choose(make(tmp_path, [], expected=None)) == ''


def test_catalogue_fov_records_the_basename(tmp_path):
    import pytest
    try:
        from photometry import catalogue_fov as cf
    except ImportError as e:      # heavy imports unavailable outside the container
        pytest.skip('catalogue_fov not importable here: %s' % e)
    sidecar = tmp_path / '.active_stack_catalogue_I+z'
    cf.write_active_catalogue(str(sidecar), str(tmp_path / TONIGHT))
    assert sidecar.read_text() == TONIGHT + '\n'
    src = open(cf.__file__).read()
    assert src.index("os.remove(active_sidecar)") < src.index('Step 2: image selection'), \
        'an earlier run\'s record must be removed before this run can fail'
