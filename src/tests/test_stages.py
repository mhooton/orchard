"""
Tests for the stage map (utils/stages.py) and the stage flags of main/ZLP_pipeline.sh (--from_TN, the 'Stages run:'
line): which stages a run runs, what a rerun from a late stage reads and rewrites, and that docs/rerun-from-stage.md
shows the map the code holds.
"""
import os
import re
import subprocess
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import stages  # noqa: E402

ZLP = os.path.join(SRC, 'main', 'ZLP_pipeline.sh')
DOC = os.path.join(os.path.dirname(SRC), 'docs', 'rerun-from-stage.md')
NIGHT = '20250317'


def zlp_stages(*flags):
    """The stages ZLP_pipeline.sh's own argument handling runs for these flags: its parsing block, cut from the
    script, run in bash, then its 'Stages run:' line."""
    text = open(ZLP, encoding='utf-8').read()
    start = text.index('T1="1"')
    end = text.index('\n', text.index('readonly T1 T2'))
    func = re.search(r'^stages_run\(\) \{\n.*?^\}\n', text, re.S | re.M).group(0)
    script = text[start:end] + '\n' + func + 'stages_run\n'
    p = subprocess.run(['bash', '-c', 'set -o nounset\n' + script, 'zlp'] + list(flags) +
                       ['1', '/data', NIGHT, '8', '2', 'Io'], capture_output=True, text=True)
    if p.returncode:
        return p.returncode, p.stderr.strip()
    return 0, stages.parse_stages_line(p.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize('flags', [
    [], ['--no_T12'], ['--from_T10', '--no_T12'], ['--from_T11'], ['--from_T10', '--no_T11', '--no_T12'],
    ['--only_T11'], ['--only_T12'], ['--no_T1', '--no_T2', '--no_T3', '--no_T4', '--no_T5', '--no_T6', '--no_T7',
                                     '--no_T8', '--no_T9'], ['--from_T7', '--force-platesolve'], ['--from_T1']])
def test_the_pipeline_and_the_stage_map_agree_on_the_stages_a_run_runs(flags):
    rc, got = zlp_stages(*flags)
    assert rc == 0, got
    assert got == stages.stages_from_argv(['./main/ZLP_pipeline.sh'] + flags + ['1', '/data', NIGHT])


def test_from_and_only_cannot_be_combined():
    rc, err = zlp_stages('--from_T10', '--only_T11')
    assert rc == 1 and 'cannot be combined' in err


def test_stages_from_argv():
    assert stages.stages_from_argv(['--from_T10', '--no_T12']) == [10, 11]
    assert stages.stages_from_argv(['--only_T10']) == [10]
    assert stages.stages_from_argv(['--cores', '8', '--frame-format', 'rice']) == list(range(1, 13))
    assert stages.parse_stages_line('Stages run: T10 T11') == [10, 11]
    assert stages.parse_stages_line('Using 20 cores') is None


def test_patterns_match_path_components():
    assert stages.matches('a_I+z_{}_5_diff.fits'.format(NIGHT), '*_diff.fits')
    assert not stages.matches('lightcurves/a_5_diff.fits', '*_diff.fits')         # '*' never crosses a '/'
    assert stages.matches('lightcurves/test_background_fix/x.png', 'lightcurves/**')
    assert stages.matches('1_image_Sp0451+0305.list', '1_image_{target}.list', target='Sp0451+0305')
    assert not stages.matches('1_image_Sp0451+0305.list', '1_image_{target}.list', target='Sp2209-3707')
    assert stages.matches('a_zYJ_{}_output.fits'.format(NIGHT), '*_{night}_output.fits', night=NIGHT)
    assert not stages.matches('a_zYJ_20250318_output.fits', '*_{night}_output.fits', night=NIGHT)


@pytest.mark.parametrize('rel,stage', [
    ('1/procX.fits', 6), ('1/zYJ_processed.dat', 6), ('1/procX.fits.phot', 9), ('1/procX.cat', 9),
    ('1/zYJ_processed.dat_phot', 9), ('1/zYJ_stacked.dat', 8), ('599_Callisto_spirit_{}_outstack_zYJ.fits', 8),
    ('599_Callisto_spirit_{}_stack_catalogue_zYJ.fits', 8), ('599_Callisto_spirit_{}_stack_catalogue_pm_zYJ.fits', 9),
    ('.active_stack_catalogue_zYJ', 8), ('599_zYJ_{}_output.fits', 9), ('599_zYJ_{}_5_diff.fits', 10),
    ('599_zYJ_{}_pwv_processed.fits', 10), ('lightcurves/599_zYJ_{}_5_MCMC', 10), ('1/procX_psf.png', None)])
def test_the_last_stage_to_write_a_file_is_its_maker(rel, stage):
    assert stages.writer(rel.format(NIGHT), night=NIGHT) == stage


def test_a_rerun_from_t10_reads_the_photometry_and_shares_every_input_with_v2():
    got = {(s, p): (m, ln) for s, p, m, ln in stages.restage_inputs(10)}
    assert got[(stages.TARGET, '*_{night}_output.fits')] == (9, True)
    assert got[(stages.REDUCTION, '1_image_{target}.list')] == (1, True)
    assert got[(stages.TARGET, stages.FRAMES)] == (6, True)              # T11's stack panel only reads them
    assert (stages.TARGET, '*_diff.fits') not in got                    # T10 makes those itself
    assert all(ln for ln in (v[1] for v in got.values()))
    assert stages.is_output('599_zYJ_{}_8_diff.fits'.format(NIGHT), 10, night=NIGHT)
    assert stages.is_output('lightcurves/599_zYJ_{}_bestap.txt'.format(NIGHT), 10, night=NIGHT)
    assert not stages.is_output('599_zYJ_{}_output.fits'.format(NIGHT), 10, night=NIGHT)
    assert not stages.is_output('1/procX.fits', 10)


def test_a_rerun_from_t11_also_reads_the_light_curves():
    got = {(s, p) for s, p, _, _ in stages.restage_inputs(11)}
    assert {(stages.TARGET, '*_diff.fits'), (stages.TARGET, 'lightcurves/**')} <= got
    assert stages.rerun_outputs(11) == []                              # the PDF is not in the target directory


def test_from_t9_back_frames_would_be_edited_in_place_so_they_cannot_be_linked():
    got = {(s, p): ln for s, p, _, ln in stages.restage_inputs(9)}
    assert got[(stages.TARGET, stages.FRAMES)] is False
    assert 9 not in stages.RESTAGE


def test_causes_resolve_to_their_first_stage():
    assert stages.cause_stage('PWV') == 10
    plan = stages.rerun_plan(stages.cause_stage('pwv'))
    assert plan['stages'] == [10, 11] and plan['flags'] == ['--from_T10', '--no_T12'] and plan['restage']
    assert stages.cause_stage('pdf') == 11
    with pytest.raises(KeyError):
        stages.cause_stage('gremlins')


def test_the_design_note_shows_the_map_the_code_holds():
    if not os.path.exists(DOC):
        pytest.skip('no docs/ beside src/ (a release or scratch tree is `git archive ... src`)')
    doc = open(DOC, encoding='utf-8').read()
    assert stages.markdown() in doc, 'docs/rerun-from-stage.md: paste `python utils/stages.py --markdown` again'
