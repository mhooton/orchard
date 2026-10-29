"""
utils/already_run_targs.py: the per-telescope list output/StackImages/already_run.txt of the targets T8 has
stacked, one row per target (Gaia ID or N, SPECULOOS name), read and rewritten by every run on the telescope.
"""
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)


def test_a_numeric_target_in_the_list_does_not_break_later_runs(tmp_path, capsys):
    # Artemis 20221028 (asteroid 03200) left the list with the single row 'N 03200'. Read with type inference
    # the SPECULOOS_ID column was integer, and `i in targ` raised TypeError in every later Artemis run
    # (136 jobs, 2026-10-26..29). Both columns are names, whatever they look like.
    pytest.importorskip('fitsio')
    from utils.already_run_targs import already_run_targs, get_runtargetlist
    outdir = tmp_path / 'output'
    runlist = outdir / 'StackImages' / 'already_run.txt'
    runlist.parent.mkdir(parents=True)
    runlist.write_text('Gaia_ID SPECULOOS_ID\nN 03200\n')
    data = get_runtargetlist(str(runlist))
    assert list(data['SPECULOOS_ID']) == ['03200'] and list(data['Gaia_ID']) == ['N']

    already_run_targs('unused', str(outdir), '03200', 'N', add_targ=False)
    assert capsys.readouterr().out.strip() == '03200'          # found, not re-added

    already_run_targs('unused', str(outdir), 'Sp0219+5912', '465423095041562752')
    assert capsys.readouterr().out.strip() == 'SP0219+5912'    # a new target comes back upper case
    assert '465423095041562752 SP0219+5912' in runlist.read_text()
