"""Checking that inputs and intermediates are actually usable, before a tool is handed them.

A metagenomics pipeline is a long chain, and a file that is present but empty or truncated
propagates quietly: an empty FASTQ produces an empty assembly, which produces zero bins,
which fails in CheckM with a message about file extensions. The cause is then hours and
several modules behind the symptom.

These checks are cheap - a stat, or reading the last few kilobytes - and are applied both to
what the user supplies and to what one module hands the next:

* :func:`require_file` - exists, is a file, and is not zero bytes
* :func:`require_nonempty_dir` - exists and contains something
* :func:`check_fastq` / :func:`check_fasta` - the format is what the extension claims, and
  the file does not end mid-record
* :func:`require_inputs` - the module-level entry point, reporting *all* the problems at once
  rather than one per re-run

Truncation is detected by reading only the tail of the file, so the cost does not scale with a
100 GB library. Compressed files are read through :mod:`metawrap2.io.seqio`, which also means
a truncated gzip stream surfaces here as a clear message rather than a CRC error from a tool.
"""

from __future__ import annotations

import gzip
import os
import zlib
from typing import List, Optional, Sequence

from .io.seqio import detect_compression, smart_open

#: How much of the end of a file to read when checking for truncation.
TAIL_BYTES = 128 * 1024

#: Everything a damaged or truncated (possibly compressed) file can raise on read. zlib.error
#: is the one a corrupt deflate stream produces, and it is not an OSError subclass.
_READ_ERRORS = (OSError, EOFError, zlib.error, gzip.BadGzipFile)


class InputProblem(Exception):
    """One or more inputs cannot be used. The message lists every problem found."""


def _describe(path: str) -> str:
    return path if len(path) < 100 else "..." + path[-97:]


def require_file(path: str, what: str = "file") -> Optional[str]:
    """Return a problem description if *path* is missing, not a file, or empty; else None."""
    if not path:
        return "no %s was given" % what
    if not os.path.exists(path):
        return "%s does not exist (%s)" % (_describe(path), what)
    if os.path.isdir(path):
        return "%s is a directory, but a %s was expected" % (_describe(path), what)
    if not os.path.isfile(path):
        return "%s is not a regular file (%s)" % (_describe(path), what)
    if os.path.getsize(path) == 0:
        return "%s is empty (0 bytes) - %s" % (_describe(path), what)
    return None


def require_nonempty_dir(
    path: str, what: str = "directory", extensions: Optional[Sequence[str]] = None
) -> Optional[str]:
    """Return a problem description if *path* is not a directory with matching files in it."""
    if not path:
        return "no %s was given" % what
    if not os.path.isdir(path):
        return "%s is not a directory (%s)" % (_describe(path), what)
    entries = os.listdir(path)
    if not entries:
        return "%s is empty (%s)" % (_describe(path), what)
    if extensions is not None:
        matching = [e for e in entries if e.endswith(tuple(extensions))]
        if not matching:
            return "%s contains no %s files (%s); it holds %d other entr%s" % (
                _describe(path),
                "/".join(extensions),
                what,
                len(entries),
                "y" if len(entries) == 1 else "ies",
            )
    return None


def _tail_text(path: str) -> str:
    """The last :data:`TAIL_BYTES` of *path* as text, decompressing if needed.

    For a compressed file the whole stream has to be read to reach the end, but only the tail
    is retained - and a truncated compressed stream raises here, which is itself the answer.
    """
    if detect_compression(path) == "none":
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            fh.seek(max(0, size - TAIL_BYTES))
            data = fh.read()
        return data.decode("utf-8", errors="replace")
    chunk = ""
    with smart_open(path, "rt") as fh:
        while True:
            block = fh.read(TAIL_BYTES)
            if not block:
                break
            chunk = (chunk + block)[-TAIL_BYTES:]
    return chunk


def check_fastq(path: str, what: str = "FASTQ file") -> Optional[str]:
    """Validate that *path* looks like a complete FASTQ file."""
    problem = require_file(path, what)
    if problem:
        return problem
    try:
        with smart_open(path, "rt") as fh:
            first = fh.readline()
    except _READ_ERRORS as exc:
        return "%s could not be read (%s): %s" % (_describe(path), what, exc)
    if not first:
        return "%s has no content (%s)" % (_describe(path), what)
    if not first.startswith("@"):
        return "%s does not look like FASTQ - the first line starts with %r, not '@' (%s)" % (
            _describe(path),
            first[:1],
            what,
        )
    try:
        tail = _tail_text(path)
    except _READ_ERRORS as exc:
        return "%s is truncated or corrupt - reading the end of it failed: %s (%s)" % (
            _describe(path),
            exc,
            what,
        )
    # A complete FASTQ ends on a record boundary: a multiple of four lines.
    lines = tail.split("\n")
    if tail and not tail.endswith("\n"):
        return "%s does not end with a newline, so its last record is truncated (%s)" % (
            _describe(path),
            what,
        )
    complete = [ln for ln in lines if ln != ""] if len(lines) < 5 else lines[:-1]
    if len(complete) % 4 != 0:
        return (
            "%s ends mid-record: the tail holds %d lines, which is not a multiple of 4 "
            "(%s). The file was probably truncated by a full disk or an interrupted "
            "download." % (_describe(path), len(complete), what)
        )
    return None


def check_fasta(path: str, what: str = "FASTA file") -> Optional[str]:
    """Validate that *path* looks like a complete FASTA file with at least one sequence."""
    problem = require_file(path, what)
    if problem:
        return problem
    try:
        with smart_open(path, "rt") as fh:
            first = fh.readline()
    except _READ_ERRORS as exc:
        return "%s could not be read (%s): %s" % (_describe(path), what, exc)
    if not first.startswith(">"):
        return "%s does not look like FASTA - the first line starts with %r, not '>' (%s)" % (
            _describe(path),
            first[:1],
            what,
        )
    try:
        tail = _tail_text(path)
    except _READ_ERRORS as exc:
        return "%s is truncated or corrupt - reading the end of it failed: %s (%s)" % (
            _describe(path),
            exc,
            what,
        )
    # A header with no sequence after it means the file stops immediately after a '>' line.
    # rpartition covers the single-line case too: with no newline it returns ("", "", tail), so
    # a file that is nothing but ">c1" is caught by this one condition. (It used to be followed
    # by an `or ... and ...` clause for that case, which Python read as `A or (B and C)` and
    # which was therefore unreachable.)
    stripped = tail.rstrip("\n")
    if stripped.rpartition("\n")[2].startswith(">"):
        return "%s ends with a header line and no sequence, so it is truncated (%s)" % (
            _describe(path),
            what,
        )
    return None


def require_inputs(checks: Sequence[tuple], context: str = "") -> None:
    """Run several checks and raise :class:`InputProblem` listing every failure.

    *checks* is a sequence of ``(callable, path, description)``. Reporting all the problems at
    once matters: a user with five samples and three bad files should learn that in one run,
    not across three successive failures.
    """
    problems: List[str] = []
    for check, path, description in checks:
        problem = check(path, description) if path is not None else "no %s was given" % description
        if problem:
            problems.append(problem)
    if problems:
        message = ["Cannot start: %d input problem(s) found." % len(problems)]
        if context:
            message[0] = "Cannot start %s: %d input problem(s) found." % (context, len(problems))
        message += ["  - " + p for p in problems]
        raise InputProblem("\n".join(message))
