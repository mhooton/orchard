"""
Unit tests for the 40 pc target table: the format-agnostic reader, the
generator that builds ml_40pc_v2.csv, and the pipeline lookups that were
repointed at it.

Covers in particular:
  - reading the legacy comma-header/whitespace-data hybrid and the new CSV
    through one code path
  - the '--' that astropy yields for an empty CSV integer cell, which used
    to be passed on as though it were a Gaia identifier
  - matching a field catalogue on DR2 *or* DR3 identifiers without
    observing one star twice
  - the generator copying every astrophysical value across verbatim

Run:  pytest src/tests/test_target_table.py -q
"""
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import target_list as tl                # noqa: E402
from utils import target_management as tm          # noqa: E402
from utils import build_target_table as btt        # noqa: E402

WOLF3 = "3864972938605115520"       # Wolf 359: DR3 only, no DR2 entry
DUAL2 = "2523085654797153280"       # Sp0049-0635 DR2 ...
DUAL3 = "2523085654797795072"       # ... and its different DR3 identifier

LEGACY_HEADER = (
    "Sp_ID, 2MASS_ID, Gaia_ID, RA, DEC, G, I, J, H, K, Dis, e_Dis, M, e_M, "
    "R, e_R,  T_eff, e_Teff, SpT, e_Spt, SNR_TESS_temp, SNR_Spec_temp, "
    "SNR_TESS_HZ, SNR_Spec_HZ, SNR_JWST_HZ_tr,  SNR_JWST_HZ_occ, "
    "SNR_JWST_temp_occ, Program\n")
LEGACY_TAIL = ("50.00 16.39 13.34 50.00 12.10 21.7400  4.7300 0.08 0.01 0.10 "
               "0.03 %s  145.  9.0  1.0   1.64 10.06  1.16  7.11  4.39  1.51  "
               "4.58 3")


def write_legacy(path, rows):
    """rows: (sp_id, gaia_id, ra, dec, teff)"""
    with open(path, "w") as f:
        f.write(LEGACY_HEADER)
        for sp, gid, ra, dec, teff in rows:
            f.write("%s  00000000-0000000 %s %.7f %.7f %s\n"
                    % (sp, gid, ra, dec, LEGACY_TAIL % teff))


def write_v2(path, rows):
    """rows: (sp_id, dr2, dr3, ra, dec, teff) — dr2/dr3 may be ''."""
    cols = list(btt.OUTPUT_COLUMNS)
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for sp, dr2, dr3, ra, dec, teff in rows:
            vals = {"Sp_ID": sp, "2MASS_ID": "00000000-0000000",
                    "Gaia_DR2_ID": dr2, "Gaia_DR3_ID": dr3,
                    "RA": "%.7f" % ra, "DEC": "%.7f" % dec, "T_eff": teff,
                    "J": "13.34", "Program": "3"}
            f.write(",".join(str(vals.get(c, "")) for c in cols) + "\n")


# --------------------------------------------------------------------------
# clean_id

@pytest.mark.parametrize("value", [
    "", " ", "nan", "NaN", "none", "None", "null", "n/a",
    "--",                                  # astropy's masked cell
    "0", "0000000000000000000",            # the legacy sentinel
    None,
])
def test_clean_id_absent_spellings(value):
    assert tl.clean_id(value) == ""


def test_clean_id_keeps_real_identifiers():
    assert tl.clean_id(" %s " % WOLF3) == WOLF3
    assert tl.clean_id(101) == "101"


# --------------------------------------------------------------------------
# reader: both formats through one path

def test_reads_legacy_format(tmp_path):
    p = tmp_path / "ml_40pc.txt"
    write_legacy(str(p), [("Sp1056+0700", WOLF3, 164.1207917, 7.0144444, "2831.")])
    t = tl.read_target_list(str(p))
    assert len(t) == 1
    assert t.sp_id(0) == "Sp1056+0700"
    assert t.ids(0) == [WOLF3]
    assert t.teff(0) == 2831
    ra, dec = t.coords(0)
    assert abs(ra - 164.1207917) < 1e-9 and abs(dec - 7.0144444) < 1e-9
    # the trailing comma the legacy header leaves is invisible to callers
    assert t.column("Program")[0] == 3


def test_reads_v2_csv(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp0049-0635", DUAL2, DUAL3, 12.3615417, -6.5963083, "2395.")])
    t = tl.read_target_list(str(p))
    assert t.sp_id(0) == "Sp0049-0635"
    assert t.teff(0) == 2395
    # DR2 first, so legacy output filenames are unchanged where DR2 exists
    assert t.ids(0) == [DUAL2, DUAL3]


def test_empty_dr2_cell_is_not_an_identifier(tmp_path):
    """
    An empty CSV field in an integer column comes back from astropy as a
    masked cell whose str() is '--'.  Before this was handled, '--' was
    passed on as though it were a Gaia ID.
    """
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp1056+0700", "", WOLF3, 164.1207917, 7.0144444, "2831.")])
    t = tl.read_target_list(str(p))
    assert t.ids(0) == [WOLF3]
    assert "--" not in t.ids(0)


def test_identical_dr2_and_dr3_listed_once(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("SpA", WOLF3, WOLF3, 1.0, 2.0, "2500.")])
    assert tl.read_target_list(str(p)).ids(0) == [WOLF3]


def test_id_rows_and_placeholder(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("SpA", DUAL2, DUAL3, 1.0, 2.0, "2500."),
                      ("SpB", "", "", 3.0, 4.0, "2400.")])
    t = tl.read_target_list(str(p))
    assert t.id_rows() == [(DUAL2, 0), (DUAL3, 0)]
    # every row present, so callers may still walk the list positionally
    assert t.id_rows(placeholder="0") == [(DUAL2, 0), (DUAL3, 0), ("0", 1)]
    assert t.row_for_id(DUAL3) == 0 and t.row_for_id("nope") is None
    assert t.teff_for_id(DUAL2) == 2500


def test_missing_id_column_is_an_error(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("Sp_ID,RA,DEC\nSpA,1.0,2.0\n")
    with pytest.raises(KeyError):
        tl.read_target_list(str(p))
    assert tl.load_target_list(str(p)) is None      # logs instead of raising


# --------------------------------------------------------------------------
# pipeline lookups against the new table

def test_crossmatch_matches_dr2_or_dr3(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp0049-0635", DUAL2, DUAL3, 12.36, -6.59, "2395."),
                      ("Sp1056+0700", "", WOLF3, 164.12, 7.01, "2831.")])
    # catalogue holds only the DR3 identifier for both stars
    got = dict((g, t) for g, t in
               tm.get_targets_from_target_list([DUAL3, WOLF3], str(p)))
    assert got == {DUAL3: 2395, WOLF3: 2831}


def test_crossmatch_does_not_observe_one_star_twice(tmp_path):
    """
    The FOV catalogue carries a DR2 *and* a DR3 column, so both of a row's
    identifiers can be in it at once.  That must still be one target.
    """
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp0049-0635", DUAL2, DUAL3, 12.36, -6.59, "2395.")])
    got = tm.get_targets_from_target_list([DUAL2, DUAL3], str(p))
    assert got == [(DUAL2, 2395)]       # DR2 preferred, reported once


def test_name_lookup_on_v2_returns_dr3_when_no_dr2(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp1056+0700", "", WOLF3, 164.1207917, 7.0144444, "2831.")])
    got = tm.get_target_from_target_list_by_name("Sp1056+0700", str(p))
    assert got[0] == WOLF3 and got[3] == 2831


def test_name_lookup_prefers_catalogue_identifier(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("Sp0049-0635", DUAL2, DUAL3, 12.36, -6.59, "2395.")])
    # whichever identifier the field catalogue holds is the one carried on
    assert tm.get_target_from_target_list_by_name(
        "Sp0049-0635", str(p), catalogue_ids={DUAL3})[0] == DUAL3
    assert tm.get_target_from_target_list_by_name(
        "Sp0049-0635", str(p))[0] == DUAL2


def test_name_lookup_skips_rows_with_no_identifier(tmp_path):
    p = tmp_path / "ml_40pc_v2.csv"
    write_v2(str(p), [("SpA", "", "", 1.0, 2.0, "2500.")])
    assert tm.get_target_from_target_list_by_name("SpA", str(p)) is None


# --------------------------------------------------------------------------
# generator

def _build(tmp_path, legacy_rows, resolved=None, id_map=None, **kw):
    src = tmp_path / "ml_40pc.txt"
    write_legacy(str(src), legacy_rows)
    warnings = []
    header, rows, stats = btt.build(
        str(src), resolved or {}, id_map or {},
        expect_rows=len(legacy_rows), warn=warnings.append, **kw)
    return header, rows, stats, warnings


def test_generator_splits_identifiers_and_cleans_names(tmp_path):
    header, rows, _, _ = _build(
        tmp_path, [("Sp0049-0635", "0", 12.3615417, -6.5963083, "2395.")],
        resolved={("Sp0049-0635", 12.36154, -6.59631): {
            "Sp_ID": "Sp0049-0635", "gaia_dr2_id": DUAL2,
            "gaia_dr3_id": DUAL3, "match_strength": "strong",
            "dr3_match": "differs", "parallax": "", "parallax_error": ""}})
    assert "Gaia_ID" not in header
    assert header[header.index("Gaia_DR2_ID")] == "Gaia_DR2_ID"
    # no column name carries the trailing comma of the legacy header
    assert not any(c.endswith(",") for c in header)
    row = dict(zip(header, rows[0]))
    assert row["Gaia_DR2_ID"] == DUAL2 and row["Gaia_DR3_ID"] == DUAL3
    assert row["dr3_match"] == "differs"


def test_generator_copies_parameters_verbatim(tmp_path):
    """Values must be byte-identical to the source: no float round-trip."""
    header, rows, _, _ = _build(
        tmp_path, [("SpA", "123", 1.2345678, -9.8765432, "2395.")],
        id_map={"123": {"dr2": "123", "dr3": "123", "kind": "identical",
                        "parallax": "", "parallax_error": ""}})
    row = dict(zip(header, rows[0]))
    assert row["T_eff"] == "2395."          # not '2395.0'
    assert row["Dis"] == "21.7400"          # not '21.74'
    assert row["RA"] == "1.2345678"
    assert row["G"] == "50.00"


def test_generator_leaves_dr2_empty_when_star_has_no_dr2_entry(tmp_path):
    header, rows, _, _ = _build(
        tmp_path, [("Sp1056+0700", "0", 164.1207917, 7.0144444, "2831.")],
        resolved={("Sp1056+0700", 164.12079, 7.01444): {
            "Sp_ID": "Sp1056+0700", "gaia_dr2_id": "", "gaia_dr3_id": WOLF3,
            "match_strength": "strong", "dr3_match": "no_dr2",
            "parallax": "415.179", "parallax_error": "0.068"}})
    row = dict(zip(header, rows[0]))
    assert row["Gaia_DR2_ID"] == "" and row["Gaia_DR3_ID"] == WOLF3
    assert row["dr3_match"] == "no_dr2"
    assert float(row["DR3_dist_pc"]) == pytest.approx(2.409, abs=1e-3)


def test_generator_excludes_weak_matches_unless_asked(tmp_path):
    resolved = {("SpW", 1.0, 2.0): {
        "Sp_ID": "SpW", "gaia_dr2_id": "", "gaia_dr3_id": "999",
        "match_strength": "weak", "dr3_match": "no_dr2",
        "parallax": "", "parallax_error": ""}}
    legacy = [("SpW", "0", 1.0, 2.0, "2500.")]

    header, rows, _, warnings = _build(tmp_path, legacy, resolved=resolved)
    row = dict(zip(header, rows[0]))
    assert row["Gaia_DR3_ID"] == "" and row["dr3_match"] == "weak_excluded"
    assert any("weak" in w for w in warnings)

    header, rows, _, _ = _build(tmp_path, legacy, resolved=resolved,
                                allow_weak=True)
    assert dict(zip(header, rows[0]))["Gaia_DR3_ID"] == "999"


def test_generator_does_not_adopt_an_ambiguous_identifier(tmp_path):
    header, rows, _, _ = _build(
        tmp_path, [("SpA", "123", 1.0, 2.0, "2500.")],
        id_map={"123": {"dr2": "123", "dr3": "456", "kind": "ambiguous",
                        "parallax": "1.0", "parallax_error": "0.1"}})
    row = dict(zip(header, rows[0]))
    assert row["Gaia_DR2_ID"] == "123"
    assert row["Gaia_DR3_ID"] == ""     # would be adopted at face value
    assert row["parallax"] == ""


def test_generator_warns_when_the_sample_changes_size(tmp_path):
    src = tmp_path / "ml_40pc.txt"
    write_legacy(str(src), [("SpA", "123", 1.0, 2.0, "2500.")])
    warnings = []
    btt.build(str(src), {}, {}, expect_rows=14168, warn=warnings.append)
    assert any("row count changed" in w for w in warnings)


def test_generator_keeps_every_row_including_distant_ones(tmp_path):
    """
    A star beyond the 40 pc selection boundary is flagged, never dropped:
    the table is a tool for finding M and L dwarfs, not a membership list.
    """
    header, rows, _, _ = _build(
        tmp_path, [("SpNear", "1", 1.0, 2.0, "2500."),
                   ("SpFar", "2", 3.0, 4.0, "2400.")],
        id_map={"1": {"dr2": "1", "dr3": "1", "kind": "identical",
                      "parallax": "100.0", "parallax_error": "0.1"},
                "2": {"dr2": "2", "dr3": "2", "kind": "identical",
                      "parallax": "10.0", "parallax_error": "0.1"}})
    assert len(rows) == 2
    by_name = {r[0]: dict(zip(header, r)) for r in rows}
    assert by_name["SpNear"]["flagged"] == "False"
    assert by_name["SpFar"]["dist_gt_40pc"] == "True"      # 100 pc
    assert by_name["SpFar"]["flagged"] == "True"


@pytest.mark.parametrize("plx,err,low,gt40,flagged", [
    ("100.0", "0.1", "False", "False", "False"),   # 10 pc, crisp
    ("10.0", "0.1", "False", "True", "True"),      # 100 pc, > 40 pc
    ("10.0", "5.0", "True", "False", "True"),      # 100 pc but < 3 sigma
])
def test_distance_flag_criteria(plx, err, low, gt40, flagged):
    _, _, l, g, f = btt.distance_fields(plx, err)
    assert (l, g, f) == (low, gt40, flagged)


@pytest.mark.parametrize("plx,err", [("-1.0", "0.1"), ("0", "0.1"),
                                     ("", ""), ("nan", "nan")])
def test_no_distance_without_a_usable_parallax(plx, err):
    """Unknown is left empty, not False: it is not the same as unflagged."""
    assert btt.distance_fields(plx, err) == ("", "", "", "", "")
