"""
Tests for utils/gaia_db_manual_rows.py: hand-made rows for real Gaia stars the database build left out.

Covers: CSV validation (band, method codes, real IDs, duplicates); check is read-only; apply inserts a
missing row and records it in manual_rows, records (without inserting) a row already present with the
same values, leaves alone a row whose ID is present with other values or that has another star within
0.5"; a second apply changes nothing; the pipeline's lookup finds an inserted row; the repo CSV is valid;
the ID an alias row maps a wrong plan-file ID to has its repo row (Sp1444-2019).

Run:  pytest src/tests/test_gaia_db_manual_rows.py -q
"""
import hashlib
import os
import sqlite3
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import gaia_db_manual_rows as mr        # noqa: E402
from utils import target_management as tm         # noqa: E402

SCHEMA = ("(ra REAL, dec REAL, pmra REAL, pmdec REAL, phot_g_mean_mag REAL, g_rp REAL, bp_rp REAL, "
          "parallax REAL, teff_gspphot REAL, source_id TEXT, dr2_source_id TEXT, j_m REAL)")
HEADER = ",".join(mr.CSV_FIELDS)

# (sp_id, band, values in COLUMNS order, methods)
NEW = ("Sp0643+1631", "16_17", (100.80829905243309, 16.528876182618006, 109.15062237550545, 65.14202763413617,
                                15.712870597839355, 1.4469108581542969, 4.202862739562988, 49.968943470333315,
                                None, "3357958852543597056", "3357958852543597056", 11.855), "DR2_TMASS_MISS")
PRESENT = ("Sp0049-0635", "-7_-6", (12.36106489018951, -6.596845239036842, -108.16708878164654, -102.62934870424121,
                                    18.256824493408203, 1.9091949462890625, 4.95197868347168, 29.844896202594548,
                                    None, "2523085654797795072", "2523085654797153280", 13.342),
           "SIX_PARAM;DR2_ID_FROM_PLANS")
CONFLICT = ("Sp1042-2416", "-25_-24", (160.67254182930927, -24.267106365293188, 40.78695328389358, 156.58201148037386,
                                       14.069744110107422, 1.4312849044799805, 4.0899248123168945, 96.14931043355692,
                                       None, "5471345889049823744", "5471345889049823744", 10.277), "PM_FROM_DR2")
NEAR = ("Sp9999-2500", "-25_-24", (161.0, -24.5, 10.0, 10.0, 15.0, 1.4, 4.0, 50.0, None,
                                    "5471345889049800000", "5471345889049800000", 11.0), "DR2_TMASS_MISS")
OTHER_STAR = (161.0 + 0.2 / 3600, -24.5, 0.0, 0.0, 17.0, 0.6, 1.0, 1.0, 5000.0, "5471345889049811111",
              "5471345889049811111", 15.0)


def csv_line(sp, band, values, methods):
    cells = [sp, band] + ["" if v is None else repr(v) if isinstance(v, float) else str(v) for v in values]
    return ",".join(cells + [methods, "test sources", "test note", "2026-10-17"])


def write_csv(path, rows, comment=True):
    text = ("# comment line\n" if comment else "") + HEADER + "\n" + "\n".join(csv_line(*r) for r in rows) + "\n"
    path.write_text(text)
    return str(path)


def make_db(path):
    conn = sqlite3.connect(path)
    for band in ("16_17", "-7_-6", "-25_-24"):
        conn.execute("CREATE TABLE '%s' %s" % (band, SCHEMA))
    conn.execute("INSERT INTO '-7_-6' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", PRESENT[2])
    changed = list(CONFLICT[2]); changed[2] = None; changed[3] = None    # same ID, other values
    conn.execute("INSERT INTO '-25_-24' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", changed)
    conn.execute("INSERT INTO '-25_-24' VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", OTHER_STAR)
    conn.commit()
    conn.close()
    return path


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def statuses(result):
    return {sp: (status, rid) for sp, band, status, rid in result}


def test_check_is_read_only_and_classifies(tmp_path):
    db = make_db(str(tmp_path / "gaia.db"))
    rows = mr.read_csv(write_csv(tmp_path / "rows.csv", [NEW, PRESENT, CONFLICT, NEAR]))
    before = md5(db)
    got = statuses(mr.run(db, rows, write=False, log=lambda *a: None))
    assert md5(db) == before
    assert got["Sp0643+1631"][0] == "missing"
    assert got["Sp0049-0635"][0] == "present, not yet recorded"
    assert got["Sp1042-2416"][0] == "conflict"
    assert got["Sp9999-2500"][0] == "near"


def test_apply_inserts_records_and_is_idempotent(tmp_path):
    db = make_db(str(tmp_path / "gaia.db"))
    rows = mr.read_csv(write_csv(tmp_path / "rows.csv", [NEW, PRESENT, CONFLICT, NEAR]))
    got = statuses(mr.run(db, rows, write=True, log=lambda *a: None))
    assert got["Sp0643+1631"][0] == "inserted"
    assert got["Sp0049-0635"] == ("recorded", 1)
    assert got["Sp1042-2416"][0] == "conflict" and got["Sp9999-2500"][0] == "near"
    conn = sqlite3.connect(db)
    assert conn.execute("SELECT count(*) FROM '16_17'").fetchone()[0] == 1
    assert conn.execute("SELECT count(*) FROM '-25_-24'").fetchone()[0] == 2     # nothing added beside them
    prov = conn.execute("SELECT band, row_id, sp_id, methods FROM manual_rows ORDER BY sp_id").fetchall()
    assert prov == [("-7_-6", 1, "Sp0049-0635", "SIX_PARAM;DR2_ID_FROM_PLANS"),
                    ("16_17", 1, "Sp0643+1631", "DR2_TMASS_MISS")]
    conn.close()
    assert mr.list_rows(db)[0][2] == "Sp0049-0635"
    # the pipeline's own lookup now finds the inserted star
    found = tm.lookup_gaia_id_in_db("3357958852543597056", dec_hint=16.53, db_path=db)
    assert found is not None and found["table"] == "16_17" and found["pmra"] == NEW[2][2]
    # a second apply changes nothing
    before = md5(db)
    again = statuses(mr.run(db, rows, write=True, log=lambda *a: None))
    assert md5(db) == before
    assert again["Sp0643+1631"][0] == "present" and again["Sp0049-0635"][0] == "present"


@pytest.mark.parametrize("change, message", [
    (lambda r: (r[0], "15_16", r[2], r[3]), "belongs in 16_17"),
    (lambda r: (r[0], r[1], r[2], "MADE_UP"), "unknown"),
    (lambda r: (r[0], r[1], r[2][:9] + ("9000000000000000001",) + r[2][10:], r[3]), "not a real Gaia"),
])
def test_read_csv_refuses_bad_rows(tmp_path, change, message):
    with pytest.raises(ValueError, match=message):
        mr.read_csv(write_csv(tmp_path / "bad.csv", [change(NEW)]))


def test_read_csv_refuses_duplicates(tmp_path):
    with pytest.raises(ValueError, match="duplicate"):
        mr.read_csv(write_csv(tmp_path / "dup.csv", [NEW, NEW]))


def test_repo_csv_is_valid():
    rows = mr.read_csv(mr.DEFAULT_CSV)
    assert rows, "the repo CSV holds the rows added so far"
    assert all(r["added"] and r["sources"] for r in rows)


def test_sp1444_2019_alias_target_has_its_row(tmp_path):
    # Sp1444-2019's plan files name a G=20.7 background star, and the alias table maps that ID to SSSPM
    # J1444-2019, which the build left out. Applied, the repo row is what T8's lookup finds in the band of
    # the plan file's J2000 dec.
    use_id, _ = tm.apply_id_alias("6281432379554902656", "Sp1444-2019", "20250401", tm.load_id_aliases())
    rows = [r for r in mr.read_csv(mr.DEFAULT_CSV) if use_id in r["values"][9:11]]
    assert [r["sp_id"] for r in rows] == ["Sp1444-2019"]
    db = str(tmp_path / "gaia.db")
    conn = sqlite3.connect(db)
    conn.execute("CREATE TABLE '%s' %s" % (rows[0]["band"], SCHEMA))
    conn.commit()
    conn.close()
    assert statuses(mr.run(db, rows, write=True, log=lambda *a: None))["Sp1444-2019"][0] == "inserted"
    found = tm.lookup_gaia_id_in_db(use_id, dec_hint=-20.3227833, db_path=db)
    assert found is not None and (found["source_id"], found["dr2_source_id"]) == (
        "6281432246412503424", "6281432246412144384")
