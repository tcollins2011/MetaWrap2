"""A database key whose value names a file, not a directory holding files.

Every check of a configured database asked "is this a directory". GUNC takes ``--db_file`` and nothing
else, so install-db correctly points GUNC_DB at the ``.dmnd`` itself - which meant the one database
whose correct value cannot be a directory was reported as missing for as long as it existed, in
`doctor` and anywhere else that described a database.
"""

from __future__ import annotations

from metawrap2 import config, dbcheck


def test_gunc_is_declared_as_a_file_valued_key():
    assert "GUNC_DB" in config.FILE_DB_KEYS
    # And it is the only one, so this is a narrow exception rather than a new convention.
    assert config.FILE_DB_KEYS == {"GUNC_DB"}
    # A directory-valued key must not be in it, or the check would invert for that one.
    assert "CHECKM_DB" not in config.FILE_DB_KEYS


def test_a_file_valued_database_pointing_at_a_file_is_usable(tmp_path):
    db = tmp_path / "gunc_db_progenomes2.1.dmnd"
    db.write_bytes(b"not really diamond, but a file")
    described = dbcheck.describe(str(db), "GUNC_DB")
    assert described["status"] == dbcheck.OK
    assert described["files"] == 1
    assert described["key_file_bytes"] == db.stat().st_size


def test_a_file_valued_database_pointing_at_a_directory_is_not(tmp_path):
    """The mistake worth catching is the opposite one: GUNC cannot be given a folder."""
    described = dbcheck.describe(str(tmp_path), "GUNC_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a file"


def test_a_missing_file_valued_database_is_missing(tmp_path):
    described = dbcheck.describe(str(tmp_path / "absent.dmnd"), "GUNC_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a file"


def test_a_directory_valued_database_still_wants_a_directory(tmp_path):
    """The change must not have loosened the rule for every other key."""
    stray = tmp_path / "checkm_is_a_file_somehow"
    stray.write_text("x")
    described = dbcheck.describe(str(stray), "CHECKM_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a directory"
