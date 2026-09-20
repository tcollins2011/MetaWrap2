"""Identifying a database, not just pointing at it.

Provenance used to record a database as a *path*: ``KRAKEN2_DB = /db/kraken2_standard``. A path
is not an identity. Replace the contents of that directory with a newer build - which is exactly
what happens when a database is updated in place - and every previous run's record now names
something that no longer exists, while claiming to describe what produced the data. "Which of my
analyses used the nt I just replaced?" has no answer.

This module gives a database a cheap, stable identity: which of its key files are present, how
big they are, and when they were last modified, hashed into one short string. That fingerprint
changes when the database is rebuilt or replaced and not when it is merely re-read, which is what
makes it worth recording. Only the handful of files that define each database are looked at -
fingerprinting a 300 GB BLAST database by content would cost more than the analysis.

The same definitions answer "is this database usable?" for :mod:`metawrap2.commands.doctor`, so
there is one place that knows what a complete CheckM or Kraken2 database looks like rather than
two that can disagree.
"""

from __future__ import annotations

import hashlib
import os
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .config import NON_PATH_DB_KEYS, Settings

OK, MISSING, INCOMPLETE = "OK", "MISSING", "INCOMPLETE"

#: What each database must contain to be usable, and how to recognise it.
#:
#: ``required`` files must all be present. ``any_of`` is satisfied by one match, for databases
#: whose filenames vary by build (the BLAST database's name is configurable; bmtagger's index is
#: named after the host genome). A key absent from this table can still be fingerprinted - it is
#: simply checked for existence only, since nothing here knows what it should contain.
DB_CONTENTS: Dict[str, Dict[str, Sequence[str]]] = {
    "KRAKEN2_DB": {
        "required": ("hash.k2d", "opts.k2d", "taxo.k2d"),
        "any_of": (),
    },
    "TAXDUMP": {
        "required": ("names.dmp", "nodes.dmp"),
        "any_of": (),
    },
    "CHECKM_DB": {
        "required": ("taxon_marker_sets.tsv",),
        # CheckM's own layout varies between the 2015 data release and later repackagings.
        "any_of": ("hmms", "genome_tree", "pfam"),
    },
    "BMTAGGER_DB": {
        "required": (),
        "any_of": ("*.bitmask",),
    },
    "BLASTDB": {
        "required": (),
        # A BLAST database is a set of numbered volumes; any one .nin/.nal proves it was built.
        "any_of": ("*.nin", "*.nal", "*.nsq"),
    },
    "GTDBTK_DATA_PATH": {
        "required": (),
        "any_of": ("taxonomy", "markers", "metadata"),
    },
    "BAKTA_DB": {
        "required": (),
        "any_of": ("bakta.db", "version.json"),
    },
}


def _matches(directory: str, pattern: str) -> List[str]:
    """Names in *directory* matching *pattern* (a literal name, or one leading ``*.ext``)."""
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return []
    if pattern.startswith("*."):
        suffix = pattern[1:]
        return [n for n in names if n.endswith(suffix)]
    return [n for n in names if n == pattern]


def key_files(path: str, key: str) -> List[str]:
    """The files/directories that define this database, as absolute paths, sorted.

    Empty when *path* is not a directory or holds none of them. Used both for the completeness
    check and as the input to the fingerprint, so the two can never describe different files.
    """
    if not path or not os.path.isdir(path):
        return []
    spec = DB_CONTENTS.get(key)
    if spec is None:
        # Unknown database: its identity is its top-level listing, which is still better than a
        # bare path and costs one readdir.
        try:
            return sorted(os.path.join(path, n) for n in os.listdir(path) if not n.startswith("."))[
                :50
            ]
        except OSError:
            return []
    found: List[str] = []
    for name in spec["required"]:
        found.extend(os.path.join(path, n) for n in _matches(path, name))
    for pattern in spec["any_of"]:
        found.extend(os.path.join(path, n) for n in _matches(path, pattern))
    return sorted(set(found))


def missing_files(path: str, key: str) -> List[str]:
    """Which of this database's defining files are absent. Empty means it looks complete."""
    spec = DB_CONTENTS.get(key)
    if spec is None or not os.path.isdir(path):
        return []
    absent = [name for name in spec["required"] if not _matches(path, name)]
    if spec["any_of"] and not any(_matches(path, pattern) for pattern in spec["any_of"]):
        absent.append("one of: " + ", ".join(spec["any_of"]))
    return absent


def fingerprint(path: str, key: str) -> Optional[str]:
    """A short, stable identity for the database at *path*, or None if there is nothing there.

    Built from each defining file's name, size and mtime - not its content. A database rebuilt or
    re-downloaded changes all three; one that is merely read changes none. That is the exact
    distinction worth recording, at the cost of a few stat calls rather than hours of hashing.
    """
    files = key_files(path, key)
    if not files:
        return None
    digest = hashlib.sha256()
    for full in files:
        try:
            stat = os.stat(full)
        except OSError:
            continue
        digest.update(
            ("%s:%d:%d\n" % (os.path.basename(full), stat.st_size, int(stat.st_mtime))).encode()
        )
    return digest.hexdigest()[:16]


def describe(path: str, key: str) -> Dict[str, Any]:
    """Everything worth recording about one configured database."""
    entry: Dict[str, Any] = {"path": path}
    if key in NON_PATH_DB_KEYS:
        # e.g. BLASTDB_NAME is a name, not a location; checking it for existence is meaningless.
        entry["kind"] = "setting"
        entry["status"] = OK if path else MISSING
        return entry
    entry["kind"] = "database"
    if not path:
        entry["status"] = MISSING
        return entry
    if not os.path.isdir(path):
        entry["status"] = MISSING
        entry["problem"] = "not a directory"
        return entry
    absent = missing_files(path, key)
    entry["status"] = INCOMPLETE if absent else OK
    if absent:
        entry["missing"] = absent
    entry["files"] = len(key_files(path, key))
    entry["fingerprint"] = fingerprint(path, key)
    total = 0
    for full in key_files(path, key):
        try:
            total += os.path.getsize(full)
        except OSError:
            continue
    entry["key_file_bytes"] = total
    return entry


def describe_all(settings: Settings) -> Dict[str, Dict[str, Any]]:
    """Every configured database, described. Unconfigured keys are left out entirely."""
    return {key: describe(value, key) for key, value in settings.databases.items() if value}


def check(settings: Settings, keys: Optional[Iterable[str]] = None) -> List[Tuple[str, str, str]]:
    """``[(key, status, detail), ...]`` for the requested database keys.

    Unlike :func:`describe_all` this includes keys that are *not* configured, reporting them as
    MISSING - because "you never set KRAKEN2_DB" is the answer a user needs when kraken2 fails,
    and an absent key cannot report itself.
    """
    results: List[Tuple[str, str, str]] = []
    for key in keys if keys is not None else sorted(settings.databases):
        value = settings.db(key)
        if not value:
            results.append((key, MISSING, "not set in the config"))
            continue
        entry = describe(value, key)
        if entry["status"] == OK:
            detail = value
            if entry.get("fingerprint"):
                detail = "%s (%s)" % (value, entry["fingerprint"])
            results.append((key, OK, detail))
        elif entry["status"] == INCOMPLETE:
            results.append(
                (key, INCOMPLETE, "%s is missing %s" % (value, ", ".join(entry["missing"])))
            )
        else:
            results.append(
                (key, MISSING, "%s %s" % (value, entry.get("problem", "does not exist")))
            )
    return results
