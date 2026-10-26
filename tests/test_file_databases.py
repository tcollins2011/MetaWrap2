"""A database key whose value names a file, not a directory holding files.

Database checks only accepted directories. GUNC takes ``--db_file``, so install-db points GUNC_DB at
the ``.dmnd`` file, and `doctor` and other database reports showed it as missing.
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

    # GUNC_DB pointing at a folder is reported.
    described = dbcheck.describe(str(tmp_path), "GUNC_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a file"

    # Other keys still require a directory.
    described = dbcheck.describe(str(db), "CHECKM_DB")
    assert described["status"] == dbcheck.MISSING
    assert described["problem"] == "not a directory"
