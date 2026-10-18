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
import shutil
import subprocess
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


#: How much is read at a time when counting a FASTQ's lines. Big enough that the newline count is
#: one memchr per chunk rather than a Python loop per line, small enough to stay in cache.
SCAN_CHUNK = 8 * 1024 * 1024


def _scan_lines(path: str):
    """Stream *path* once and return ``(first_line, line_count, ends_with_newline)``.

    Counting newlines over binary chunks rather than iterating lines in Python, and through `pigz`
    when it is installed. The check itself is worth keeping - it catches a truncated gzip member,
    which is what a concatenated download leaves behind, and a user's run died on exactly that -
    but it was costing ~17 minutes of one core on a 21 GB + 23 GB pair before MEGAHIT even started:
    the decompression is C, the per-line loop was not.

    `pigz -dc` parallelises the decompression itself; without it this falls back to the same
    streaming read, which is still several times faster than a line loop. Either way the whole
    stream is read, so a partial member still raises.
    """
    first = b""
    count = 0
    last_byte = b""
    use_pigz = detect_compression(path) == "gzip" and shutil.which("pigz")
    process = None
    if use_pigz:
        process = subprocess.Popen(
            ["pigz", "-dc", path], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        assert process.stdout is not None
        stream = process.stdout
    else:
        stream = smart_open(path, "rb")
    try:
        while True:
            chunk = stream.read(SCAN_CHUNK)
            if not chunk:
                break
            if not first:
                head, _nl, _rest = chunk.partition(b"\n")
                first = head
            count += chunk.count(b"\n")
            last_byte = chunk[-1:]
    finally:
        stream.close()
        if process is not None:
            stderr = process.stderr.read() if process.stderr is not None else b""
            if process.stderr is not None:
                process.stderr.close()
            if process.wait() != 0:
                # Raised as a read error so the caller reports it the way it reports a bad gzip:
                # pigz exits non-zero on a truncated member, which is the case this check is for.
                raise OSError(
                    "pigz could not read %s: %s" % (path, stderr.decode("utf-8", "replace").strip())
                )
    ends_with_newline = last_byte == b"\n"
    if last_byte and not ends_with_newline:
        # The unterminated tail is a line too, which is what makes the multiple-of-four check see
        # a truncated record rather than a round number.
        count += 1
    return first.decode("utf-8", "replace"), count, ends_with_newline


def check_fastq(path: str, what: str = "FASTQ file") -> Optional[str]:
    """Validate that *path* looks like a complete FASTQ file."""
    problem = require_file(path, what)
    if problem:
        return problem
    # One streaming pass: validate the first line, count the total lines (a complete FASTQ is a
    # multiple of four), and remember the last line. Reading the whole stream also surfaces a
    # truncated gzip (a partial member raises here). We must count *all* lines, not a tail window:
    # the number of lines in an arbitrary byte-tail is not a multiple of four even for a valid file.
    try:
        first, n_lines, ends_with_newline = _scan_lines(path)
    except _READ_ERRORS as exc:
        return "%s is truncated or corrupt - reading it failed: %s (%s)" % (
            _describe(path),
            exc,
            what,
        )
    if not n_lines and not first:
        return "%s has no content (%s)" % (_describe(path), what)
    if not first.startswith("@"):
        return "%s does not look like FASTQ - the first line starts with %r, not '@' (%s)" % (
            _describe(path),
            first[:1],
            what,
        )
    if not ends_with_newline:
        return "%s does not end with a newline, so its last record is truncated (%s)" % (
            _describe(path),
            what,
        )
    if n_lines % 4 != 0:
        return (
            "%s ends mid-record: it holds %d lines, which is not a multiple of 4 (%s). The file "
            "was probably truncated by a full disk or an interrupted download."
            % (_describe(path), n_lines, what)
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
