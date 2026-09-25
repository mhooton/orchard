"""
Unit tests for Gaia-ID-independent target identification.

Covers the paths added for stars that are missing from Gaia DR2, from the
local database, or from Gaia altogether:
  - plan-file coordinate parsing
  - ID alias table
  - DR3-only ID matching
  - coordinate fallback with proper-motion propagation and ID injection
  - refusal of stale backup catalogues for targets without proper motion
  - supplementary (hand-added) sources

Run:  pytest src/tests/test_target_identification.py -q
"""
import os
import sys
import sqlite3

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import target_management as tm          # noqa: E402
from utils import gaia_id_from_schedule as gs      # noqa: E402
from utils import supplementary_sources as ss      # noqa: E402

WOLF = "3864972938605115520"        # Wolf 359, Gaia DR3 only
NEIGH = "3864972728150826368"       # G=16 neighbour SPOCK wrote into 2021-22 plans
WOLF_RA, WOLF_DEC = 164.10319030755974, 7.002726940984864   # J2016.0
WOLF_PMRA, WOLF_PMDEC = -3866.338, -2699.215


def make_plan(obsdir, date, targname, gaia_id, ra="10 56 28.99", dec="+07 00 52.00"):
    iso = "%s-%s-%s" % (date[:4], date[4:6], date[6:])
    d = os.path.join(obsdir, "schedule", "Plans_by_date", iso)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "Obj_%s.txt" % targname)
    with open(p, "w") as f:
        f.write(";\n; %s\n;\n; %s\n;\n#waituntil 1, 23:16\n#filter r'\n%s\t%s\t%s\n;\n#quitat 04:01\n"
                % (targname, gaia_id, targname, ra, dec))
    return p


def make_catalogue(path, ra, dec, flux, dr2, dr3, pmra=None, pmdec=None):
    n = len(ra)
    apm = fits.BinTableHDU.from_columns([
        fits.Column(name="Sequence_number", format="J", array=np.arange(n)),
        fits.Column(name="RA", format="D", array=np.radians(ra)),
        fits.Column(name="DEC", format="D", array=np.radians(dec)),
        fits.Column(name="Isophotal_flux", format="E", array=np.asarray(flux, float)),
    ], name="APM-BINARYTABLE")
    nanv = np.full(n, np.nan)
    pmra = nanv if pmra is None else np.asarray(pmra, float)
    pmdec = nanv if pmdec is None else np.asarray(pmdec, float)
    xm = fits.BinTableHDU.from_columns([
        fits.Column(name="GAIA_DR2_ID", format="26A", array=np.array(dr2)),
        fits.Column(name="GAIA_DR3_ID", format="26A", array=np.array(dr3)),
        fits.Column(name="PARALLAX", format="D", array=nanv),
        fits.Column(name="GMAG", format="D", array=nanv),
        fits.Column(name="G_RP", format="D", array=nanv),
        fits.Column(name="BP_RP", format="D", array=nanv),
        fits.Column(name="TEFF", format="D", array=nanv),
        fits.Column(name="PMRA", format="D", array=pmra),
        fits.Column(name="PMDEC", format="D", array=pmdec),
    ], name="GAIA_CROSSMATCH")
    fits.HDUList([fits.PrimaryHDU(), apm, xm]).writeto(path, overwrite=True)
    return path


def make_db(path, rows):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE '7_8' (ra REAL, dec REAL, pmra REAL, pmdec REAL, "
                 "phot_g_mean_mag REAL, g_rp REAL, bp_rp REAL, parallax REAL, "
                 "teff_gspphot REAL, source_id TEXT, dr2_source_id TEXT, j_m REAL)")
    conn.executemany("INSERT INTO '7_8' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    conn.commit()
    conn.close()
    return path


WOLF_DB_ROW = (WOLF_RA, WOLF_DEC, WOLF_PMRA, WOLF_PMDEC, 11.04, 1.45, 4.18, 415.18, None, WOLF, None, 7.085)


# --------------------------------------------------------------------------
def test_read_coords(tmp_path):
    p = make_plan(str(tmp_path), "20220104", "Sp1056+0700", NEIGH)
    ra, dec = gs.read_coords(p)
    assert abs(ra - 164.12079) < 1e-4 and abs(dec - 7.01444) < 1e-4
    assert gs.read_file(p, silent=True) == NEIGH


def test_alias_table_from_repo_and_matching():
    aliases = tm.load_id_aliases()
    assert any(a["schedule_id"] == NEIGH and a["use_id"] == WOLF for a in aliases)
    assert tm.apply_id_alias(NEIGH, "Sp1056+0700", "20220104", aliases)[0] == WOLF
    assert tm.apply_id_alias(NEIGH, "Sp1056--0700".replace("--", "+"), "20210401", aliases)[0] == WOLF
    # outside the date range, or a different target: untouched
    assert tm.apply_id_alias(NEIGH, "Sp1056+0700", "20230114", aliases)[0] == NEIGH
    assert tm.apply_id_alias(NEIGH, "Sp0000+0000", "20220104", aliases)[0] == NEIGH
    assert tm.apply_id_alias(None, "x", "20220104", aliases) == (None, None)


def test_schedule_id_matches_dr3_only_row(tmp_path):
    make_plan(str(tmp_path), "20240422", "Sp1056+0700", WOLF)
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), "20240422", "Sp1056+0700",
        catalogue_dr2_ids=["111", "nan", "222"],
        target_list_path=str(tmp_path / "missing.txt"),
        toi_table_path=str(tmp_path / "missing.csv"),
        catalogue_dr3_ids=["111", WOLF, "222"],
        info=info)
    assert targets == [(WOLF, "primary", None)]
    assert info["status"] == "OK" and info["primary_method"] == "schedule"
    assert info["primary_row"] == 1


def test_alias_applied_in_identification(tmp_path):
    make_plan(str(tmp_path), "20220104", "Sp1056+0700", NEIGH)
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), "20220104", "Sp1056+0700",
        catalogue_dr2_ids=[NEIGH, "nan"],
        target_list_path="/nonexistent", toi_table_path="/nonexistent",
        catalogue_dr3_ids=[NEIGH, WOLF], info=info)
    assert [t for t in targets if t[1] == "primary"] == [(WOLF, "primary", None)]
    assert info["primary_row"] == 1


def test_coordinate_fallback_with_db_pm(tmp_path):
    date = "20240422"
    make_plan(str(tmp_path), date, "Sp1056+0700", WOLF)
    db = make_db(str(tmp_path / "gaia.db"), [WOLF_DB_ROW])
    ra_now, dec_now = tm._propagate(WOLF_RA, WOLF_DEC, WOLF_PMRA, WOLF_PMDEC, 2016.0,
                                    tm._date_to_epoch(date))
    coords = np.array([[164.0900, 7.0000], [ra_now + 0.4 / 3600, dec_now - 0.2 / 3600],
                       [164.1500, 7.0500]])
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), date, "Sp1056+0700",
        catalogue_dr2_ids=["123", "nan", "456"],
        target_list_path="/nonexistent", toi_table_path="/nonexistent",
        catalogue_dr3_ids=["123", "nan", "456"],
        catalogue_coords=coords, catalogue_fluxes=[100.0, 9e5, 200.0],
        db_path=db, info=info)
    assert targets == [(WOLF, "primary", None)]
    assert info["status"] == "OK_COORDINATES" and info["primary_method"] == "coordinates"
    assert info["inject"]["row"] == 1 and info["inject"]["gaia_id"] == WOLF
    assert info["inject"]["separation_arcsec"] < 1.0
    assert info["db_row"]["source_id"] == WOLF


def test_coordinate_fallback_not_in_db_fast_star_is_missed(tmp_path):
    """Without a proper motion the J2000 plan position is 2' off a 4.7''/yr star."""
    date = "20240422"
    make_plan(str(tmp_path), date, "Sp1056+0700", WOLF)
    db = make_db(str(tmp_path / "gaia.db"), [])
    ra_now, dec_now = tm._propagate(WOLF_RA, WOLF_DEC, WOLF_PMRA, WOLF_PMDEC, 2016.0,
                                    tm._date_to_epoch(date))
    coords = np.array([[ra_now, dec_now]])
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), date, "Sp1056+0700",
        catalogue_dr2_ids=["nan"], target_list_path="/nonexistent",
        toi_table_path="/nonexistent", catalogue_dr3_ids=["nan"],
        catalogue_coords=coords, db_path=db, info=info)
    assert targets == []
    assert info["status"].startswith("PRIMARY_NOT_IN_DB") and "COORD_NO_MATCH" in info["status"]


def test_coordinate_fallback_slow_star_from_plan_coords(tmp_path):
    date = "20240422"
    make_plan(str(tmp_path), date, "Sp0000+0000", "5555555555555555555",
              ra="00 40 00.00", dec="+10 00 00.0")
    db = make_db(str(tmp_path / "gaia.db"), [])
    coords = np.array([[10.0 + 0.3 / 3600, 10.0], [10.2, 10.1]])
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), date, "Sp0000+0000",
        catalogue_dr2_ids=["nan", "9"], target_list_path="/nonexistent",
        toi_table_path="/nonexistent", catalogue_dr3_ids=["nan", "9"],
        catalogue_coords=coords, db_path=db, info=info)
    assert targets == [("5555555555555555555", "primary", None)]
    assert info["primary_method"] == "coordinates" and info["inject"]["row"] == 0


def test_coordinate_conflict_refuses_override(tmp_path):
    date = "20240422"
    make_plan(str(tmp_path), date, "Sp1056+0700", WOLF)
    db = make_db(str(tmp_path / "gaia.db"), [WOLF_DB_ROW])
    ra_now, dec_now = tm._propagate(WOLF_RA, WOLF_DEC, WOLF_PMRA, WOLF_PMDEC, 2016.0,
                                    tm._date_to_epoch(date))
    info = {}
    targets = tm.identify_targets(
        str(tmp_path), date, "Sp1056+0700",
        catalogue_dr2_ids=["777"], target_list_path="/nonexistent",
        toi_table_path="/nonexistent", catalogue_dr3_ids=["777"],
        catalogue_coords=np.array([[ra_now, dec_now]]), db_path=db, info=info)
    assert targets == [] and "COORD_CONFLICT" in info["status"]


def test_catalogue_is_valid_pm_requirement(tmp_path):
    cat = make_catalogue(str(tmp_path / "c.fits"), [164.1, 164.2], [7.0, 7.1], [1e5, 1e3],
                         ["nan", "111"], [WOLF, "111"])
    assert tm.catalogue_is_valid(cat, WOLF) is True
    assert tm.catalogue_is_valid(cat, WOLF, require_pm=True) is False
    cat2 = make_catalogue(str(tmp_path / "c2.fits"), [164.1], [7.0], [1e5], ["nan"], [WOLF],
                          pmra=[WOLF_PMRA], pmdec=[WOLF_PMDEC])
    assert tm.catalogue_is_valid(cat2, WOLF, require_pm=True) is True
    assert tm.catalogue_is_valid(cat2, "999") is False


def test_find_backup_prefers_valid_and_rejects_stale_without_pm(tmp_path):
    bdir = tmp_path / "StackImages"
    bdir.mkdir()

    def add(date, with_pm):
        make_catalogue(str(bdir / ("%s_Io_andor_%s_stack_catalogue_r.fits" % (WOLF, date))),
                       [164.1], [7.0], [1e5], ["nan"], [WOLF],
                       pmra=[WOLF_PMRA] if with_pm else None,
                       pmdec=[WOLF_PMDEC] if with_pm else None)
        (bdir / ("%s_Io_andor_%s_outstack_r.fits" % (WOLF, date))).write_bytes(b"")

    add("20210401", False)                       # 3 years old, no PM: unusable
    stack, cat = tm.find_backup(str(bdir), WOLF, "andor", "20240422")
    assert stack is None and cat is None
    add("20240301", True)                        # 52 days, has PM: fine
    stack, cat = tm.find_backup(str(bdir), WOLF, "andor", "20240422")
    assert cat.endswith("20240301_stack_catalogue_r.fits")
    add("20240410", False)                       # 12 days, no PM: still fine
    stack, cat = tm.find_backup(str(bdir), WOLF, "andor", "20240422")
    assert cat.endswith("20240410_stack_catalogue_r.fits")


def test_supplementary_sources_roundtrip(tmp_path):
    csvp = tmp_path / "sup.csv"
    csvp.write_text(
        "# comment\nname,source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,g_rp,bp_rp,teff_gspphot,j_m,provenance\n"
        "Fake-Y1,9000000000000000001,10.5,-20.25,100,-50,80,21.5,1.6,4.5,500,15.9,test\n")
    rows = ss.read_csv(str(csvp))
    assert len(rows) == 1 and ss.is_synthetic_id(rows[0]["source_id"])
    bad = tmp_path / "bad.csv"
    bad.write_text("name,source_id,ra,dec\nX,123,1,2\n")
    with pytest.raises(ValueError):
        ss.read_csv(str(bad))
    db = make_db(str(tmp_path / "gaia.db"), [])
    assert ss.load(db, str(csvp)) == 1
    row = tm.lookup_gaia_id_in_db("9000000000000000001", dec_hint=None, db_path=db)
    assert row is not None and row["table"] == "custom_sources" and row["ra"] == 10.5
    assert tm.lookup_gaia_id_in_db("9000000000000000002", db_path=db) is None
    # repo CSV must be loadable even when empty
    assert ss.read_csv(ss.DEFAULT_CSV) == []


def test_db_query_includes_custom_sources(tmp_path, monkeypatch):
    fitsio = pytest.importorskip("fitsio")  # noqa: F841
    from photometry import gaia_dr2_test as gx
    db = make_db(str(tmp_path / "gaia.db"), [WOLF_DB_ROW])
    csvp = tmp_path / "sup.csv"
    csvp.write_text("name,source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,g_rp,bp_rp,teff_gspphot,j_m,provenance\n"
                    "Fake,9000000000000000001,164.2,7.05,0,0,10,18,1,2,900,15,test\n")
    ss.load(db, str(csvp))
    monkeypatch.setenv("GAIADATABASEPATH", db)
    rows = gx._db_query(6.9, 7.2, 164.0, 164.3)
    ids = {r["source_id"] for r in rows}
    assert WOLF in ids and "9000000000000000001" in ids
