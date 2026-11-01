"""Identifying a database, not just pointing at it.

Provenance used to record a database only as a path (``KRAKEN2_DB = /db/kraken2_standard``). If
the directory is updated in place, earlier run records point at contents that no longer exist,
and there is no way to tell which runs used the old build.

This module fingerprints a database cheaply: the presence, size and mtime of its key files,
hashed into one short string. The fingerprint changes when the database is rebuilt or replaced,
not when it is read. Only the defining files are looked at; hashing a 300 GB BLAST database by
content would cost more than the analysis.

The same definitions are used by :mod:`metawrap2.commands.doctor` to check whether a database is
complete.

It also decides which databases a module invocation needs (from
:data:`metawrap2.config.MODULE_DB_NEEDS` and the flags given), for `check`, `doctor` and the
sample-sheet preflight in `run`/`cohort`.
"""

from __future__ import annotations

import difflib
import hashlib
import os
import shlex
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

from .config import FILE_DB_KEYS, MODULE_DB_NEEDS, NON_PATH_DB_KEYS, DbNeed, Settings

OK, MISSING, INCOMPLETE = "OK", "MISSING", "INCOMPLETE"

#: What each database must contain to be usable, and how to recognise it.
#:
#: ``required`` files must all be present. ``any_of`` is satisfied by one match, for databases
#: whose filenames vary by build (the BLAST database's name is configurable; bmtagger's index is
#: named after the host genome). A key absent from this table is only checked for existence.
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
    check and as the fingerprint input, so both describe the same files.
    """
    if not path or not os.path.isdir(path):
        return []
    spec = DB_CONTENTS.get(key)
    if spec is None:
        # Unknown database: use its top-level listing (one readdir).
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

    Built from each defining file's name, size and mtime, not its content. A rebuild or
    re-download changes them; reading does not. This costs a few stat calls instead of hashing.
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
    """Provenance record for one configured database."""
    entry: Dict[str, Any] = {"path": path}
    if key in NON_PATH_DB_KEYS:
        # e.g. BLASTDB_NAME is a name, not a location, so it is not checked for existence.
        entry["kind"] = "setting"
        entry["status"] = OK if path else MISSING
        return entry
    entry["kind"] = "database"
    if not path:
        entry["status"] = MISSING
        return entry
    if key in FILE_DB_KEYS:
        # This key names the database file itself (GUNC's --db_file), not a directory.
        if not os.path.isfile(path):
            entry["status"] = MISSING
            entry["problem"] = "not a file"
            return entry
        entry["status"] = OK
        entry["files"] = 1
        try:
            entry["key_file_bytes"] = os.path.getsize(path)
        except OSError:
            entry["key_file_bytes"] = 0
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

    Unlike :func:`describe_all` this includes keys that are not configured, reported as MISSING,
    so a failing module can be traced to an unset key.
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


# -- which databases a module needs ---------------------------------------------------------------


def install_name(key: str) -> Optional[str]:
    """The `metawrap2 install-db` name that installs *key*, or None if it has no installer."""
    from .commands.install_db import DATABASES

    for name, database in DATABASES.items():
        if database.key == key:
            return name
    return None


def install_hint(key: str) -> str:
    """How to get *key*: the install-db command, or where to set it."""
    name = install_name(key)
    if name:
        return "install with 'metawrap2 install-db %s'" % name
    return "set it under [databases] in metawrap2.toml"


def is_optional(need: DbNeed) -> bool:
    """True when the module can run without this database."""
    return bool(need.when or need.skip or need.warn_only)


def need_label(need: DbNeed) -> str:
    """``""`` for a required database, else why it is optional, naming the flag."""
    parts = []
    if need.when:
        parts.append("only with %s" % need.when)
    if need.skip:
        parts.append("skip with %s" % " or ".join(need.skip))
    if need.warn_only and need.note and not parts:
        parts.append("without it %s" % need.note)
    return "(optional - %s)" % "; ".join(parts) if parts else ""


def label_for(key: str, modules: Iterable[str]) -> str:
    """The label for *key* across *modules*: ``""`` if any of them requires it.

    Otherwise the label of a need that applies by default (``--skip-gunc``) is preferred over one
    that a flag turns on (``--checkm1``).
    """
    needs = [n for m in modules for n in module_needs(m) if n.key == key]
    if not needs or any(not is_optional(n) for n in needs):
        return ""
    default_on = [n for n in needs if not n.when]
    return need_label((default_on or needs)[0])


def module_needs(module: str) -> Tuple[DbNeed, ...]:
    """Every database *module* can read, needed or not."""
    return MODULE_DB_NEEDS.get(module, ())


def keys_for(modules: Iterable[str]) -> List[str]:
    """The database keys *modules* can read, in first-seen order."""
    keys: List[str] = []
    for module in modules:
        for need in module_needs(module):
            if need.key not in keys:
                keys.append(need.key)
    return keys


def _flags(argv: Sequence[str]) -> Dict[str, Optional[str]]:
    """``--flag`` -> its value (the next token, or after ``=``), or None for a switch."""
    flags: Dict[str, Optional[str]] = {}
    for i, token in enumerate(argv):
        if not token.startswith("--"):
            continue
        name, eq, value = token.partition("=")
        if eq:
            flags[name] = value
        elif i + 1 < len(argv) and not argv[i + 1].startswith("-"):
            flags[name] = argv[i + 1]
        else:
            flags[name] = None
    return flags


def active_needs(module: str, argv: Sequence[str]) -> List[DbNeed]:
    """The databases this invocation of *module* reads, given its arguments."""
    flags = _flags(argv)
    return [
        need
        for need in module_needs(module)
        if (not need.when or need.when in flags) and not any(f in flags for f in need.skip)
    ]


def problem(settings: Settings, key: str, path: Optional[str] = None) -> Optional[str]:
    """Why database *key* cannot be used, or None if it can. *path* overrides the config."""
    value = path if path is not None else settings.db(key)
    if not value:
        return "%s is not set (%s)" % (key, install_hint(key))
    entry = describe(value, key)
    if entry["status"] == OK:
        return None
    if entry["status"] == INCOMPLETE:
        return "%s at %s is missing %s (%s)" % (
            key,
            value,
            ", ".join(entry["missing"]),
            install_hint(key),
        )
    return "%s points at %s, which does not exist (%s)" % (key, value, install_hint(key))


class Resolved(NamedTuple):
    """One database an invocation reads, where it was looked for, and whether it is usable."""

    key: str
    module: str
    #: The path checked: from the module's path flag, the config, or the fallback key. ``""`` when
    #: none is set.
    path: str
    #: The path flag the path came from (``--kofam-db``), or ``""`` when it came from the config.
    flag: str
    #: Why the database cannot be used, or None when it can.
    problem: Optional[str]
    #: True when the module only warns about this database being unset, and runs without it.
    warn_only: bool


def resolve(settings: Settings, invocations: Iterable[Tuple[str, Sequence[str]]]) -> List[Resolved]:
    """Each database a list of ``(module, argv)`` invocations reads, once per key and path.

    A database passed with its path flag is checked at that path. Only stat calls: nothing is run.
    """
    found: List[Resolved] = []
    seen = set()
    for module, argv in invocations:
        flags = _flags(argv)
        for need in active_needs(module, argv):
            path = flags.get(need.path_flag) if need.path_flag else None
            from_flag = path is not None
            if path is None and not settings.db(need.key) and need.fallback:
                path = settings.db(need.fallback) or None
            ident = (need.key, path)
            if ident in seen:
                continue
            seen.add(ident)
            reason = problem(settings, need.key, path)
            if reason is not None and from_flag:
                reason = reason.replace(need.key, need.path_flag, 1)
            unset = path is None and not settings.db(need.key)
            found.append(
                Resolved(
                    key=need.key,
                    module=module,
                    path=path if path is not None else settings.db(need.key),
                    flag=need.path_flag if from_flag else "",
                    problem=reason,
                    warn_only=need.warn_only and unset,
                )
            )
    return found


def preflight(
    settings: Settings, invocations: Iterable[Tuple[str, Sequence[str]]]
) -> Tuple[List[str], List[str]]:
    """``(errors, warnings)`` for the databases a list of ``(module, argv)`` invocations need.

    Each database is reported once, however many invocations need it. A database passed with its
    path flag is checked at that path. An unset database the module only warns about is a warning.
    Only stat calls: nothing is run.
    """
    errors: List[str] = []
    warnings: List[str] = []
    for entry in resolve(settings, invocations):
        if entry.problem is None:
            continue
        message = "%s (%s)" % (entry.problem, entry.module)
        if entry.warn_only:
            warnings.append(message)
        else:
            errors.append(message)
    return errors, warnings


# -- reading a typed command line ------------------------------------------------------------------

#: Global flags without a value, accepted before the module name.
_GLOBAL_SWITCHES = {
    "--dry-run",
    "--force",
    "--resume",
    "--verbose-logs",
    "--skip-space-check",
    "--skip-validation",
    "--gpu-strict",
}


def split_command(text: str) -> Tuple[str, List[str]]:
    """``(module, argv)`` from a typed module command such as ``"annotate --kofam -o OUT"``.

    A leading ``metawrap2`` and the global flags before the module name are dropped. Raises
    ValueError when the text cannot be split or names no module.
    """
    from .cli import _LIMIT_FLAGS, _SETTING_FLAGS, MODULES

    try:
        tokens = shlex.split(text)
    except ValueError as exc:
        raise ValueError("cannot split %r: %s" % (text, exc)) from None
    if tokens and os.path.basename(tokens[0]) == "metawrap2":
        tokens = tokens[1:]
    while tokens:
        head = tokens[0]
        if (
            head in _GLOBAL_SWITCHES
            or head.split("=", 1)[0] in _LIMIT_FLAGS.keys() | _SETTING_FLAGS
        ):
            tokens = tokens[1 if "=" in head or head in _GLOBAL_SWITCHES else 2 :]
        else:
            break
    if not tokens:
        raise ValueError("%r names no module" % text)
    module, argv = tokens[0], tokens[1:]
    if module not in MODULES:
        close = difflib.get_close_matches(module, list(MODULES), n=1)
        hint = " Did you mean '%s'?" % close[0] if close else ""
        raise ValueError(
            "'%s' is not a module.%s Modules: %s" % (module, hint, ", ".join(sorted(MODULES)))
        )
    return module, argv
