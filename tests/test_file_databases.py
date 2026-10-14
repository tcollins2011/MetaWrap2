"""A database key whose value names a file, not a directory holding files.

Every check of a configured database asked "is this a directory". GUNC takes ``--db_file`` and nothing
else, so install-db correctly points GUNC_DB at the ``.dmnd`` itself - which meant the one database
whose correct value cannot be a directory was reported as missing for as long as it existed, in
`doctor` and anywhere else that described a database.
"""

from __future__ import annotations

from metawrap2 import dbcheck


def test_a_file_valued_database_pointing_at_a_file_is_usable(tmp_path):
    db = tmp_path / "gunc_db_progenomes2.1.dmnd"
    db.write_bytes(b"not really diamond, but a file")
    described = dbcheck.describe(str(db), "GUNC_DB")
    assert described["status"] == dbcheck.OK
    assert described["files"] == 1
    assert described["key_file_bytes"] == db.stat().st_size

    # The mistake worth catching is the opposite one: GUNC cannot be given a folder.
    described = dbcheck.describe(str(tmp_path), "GUNC_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a file"

    # The change must not have loosened the rule for every other key.
    described = dbcheck.describe(str(db), "CHECKM_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a directory"
