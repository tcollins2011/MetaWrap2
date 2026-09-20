"""Locating and validating a BLAST nucleotide database.

Both ``blobology`` and ``classify_bins`` MEGABLAST contigs against a nucleotide database.
They used to hardcode two assumptions: that the database is named exactly ``nt``, and that
``<BLASTDB>/nt.00.nhd`` exists. Neither holds in general:

* BLAST v5 databases number volumes with three digits (``nt.000.nhd``), and a database
  small enough to fit in one volume has no volume suffix at all (``nt.nhd``);
* ``.nhd`` only exists when the database was built with hashed sequence ids, so a perfectly
  good ``makeblastdb`` output can lack it entirely;
* users legitimately want to search something other than ``nt`` - core_nt, a RefSeq subset,
  or their own curated database.

So the database *name* is configurable (``BLASTDB_NAME``, default ``nt``) and presence is
established by looking for any of the files BLAST actually requires, whatever the volume
numbering.
"""

from __future__ import annotations

import glob
import os
from typing import Optional, Tuple

#: Extensions that every blast nucleotide DB volume has, regardless of build options.
_REQUIRED_EXTS = (".nin", ".nsq")
#: An alias file stands in for a set of volumes and has no .nin/.nsq of its own.
_ALIAS_EXTS = (".nal", ".nvl")

DEFAULT_DB_NAME = "nt"


def db_name(settings) -> str:
    """The BLAST database name to search inside the BLASTDB directory."""
    return settings.db("BLASTDB_NAME") or DEFAULT_DB_NAME


def db_path(settings) -> str:
    """The ``-db`` argument to pass to blastn: ``<BLASTDB>/<name>``."""
    return os.path.join(settings.db("BLASTDB"), db_name(settings))


def find(settings) -> Tuple[bool, Optional[str]]:
    """Return ``(present, problem)`` for the configured BLAST database.

    ``problem`` is a ready-to-print explanation when the database cannot be used, else None.
    """
    directory = settings.db("BLASTDB")
    name = db_name(settings)
    if not directory:
        return False, (
            "BLASTDB is not set. Point it at a directory holding a BLAST "
            "nucleotide database in [databases] of metawrap2.toml."
        )
    if not os.path.isdir(directory):
        return False, "BLASTDB directory %s does not exist." % directory

    base = os.path.join(directory, name)
    for ext in _ALIAS_EXTS:
        if os.path.isfile(base + ext):
            return True, None
    for ext in _REQUIRED_EXTS:
        # single-volume (nt.nsq) or multi-volume, 2- or 3-digit (nt.00.nsq / nt.000.nsq)
        if os.path.isfile(base + ext) or glob.glob(base + ".*" + ext):
            return True, None

    present = sorted(os.path.basename(p) for p in glob.glob(os.path.join(directory, "*")))[:8]
    return False, (
        "No BLAST nucleotide database named '%s' was found in %s (looked for %s%s and "
        "volume-numbered variants). Set BLASTDB to the database directory and, if your "
        "database is not called 'nt', set BLASTDB_NAME. Directory currently contains: %s"
        % (
            name,
            directory,
            name,
            "/".join(_REQUIRED_EXTS + _ALIAS_EXTS),
            ", ".join(present) or "<nothing>",
        )
    )
