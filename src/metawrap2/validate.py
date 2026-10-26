"""Checking that inputs and intermediates are usable before a tool is given them.

A file that is present but empty or truncated propagates down the pipeline: an empty FASTQ gives an
empty assembly, then zero bins, then a CheckM error about file extensions, several modules later.

The checks are cheap (a stat, or reading the last few kilobytes) and are applied to user inputs and
to what one module passes to the next:

* :func:`require_file` - exists, is a file, and is not zero bytes
* :func:`require_nonempty_dir` - exists and contains something
* :func:`check_fastq` / :func:`check_fasta` - the format matches the extension, and the file does
  not end mid-record
* :func:`require_inputs` - the module-level entry point, reporting all problems at once

FASTA truncation is detected from the tail of the file. Compressed files are read through
:mod:`metawrap2.io.seqio`, so a truncated gzip stream is reported here with a clear message.
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
#: (corrupt deflate stream) is not an OSError subclass.
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

    For a compressed file the whole stream is read to reach the end, keeping only the tail. A
    truncated compressed stream raises here.
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


#: Read size when counting a FASTQ's lines. Large enough for one memchr per chunk, small enough to
#: stay in cache.
SCAN_CHUNK = 8 * 1024 * 1024


def _scan_lines(path: str):
    """Stream *path* once and return ``(first_line, line_count, ends_with_newline)``.

    Counts newlines over binary chunks, decompressing through `pigz -dc` when it is installed. The
    check catches truncated gzip members (left by interrupted or concatenated downloads). A
    per-line Python loop took ~17 minutes of one core on a 21 GB + 23 GB pair before MEGAHIT
    started.

    Without pigz it falls back to a streaming read, still several times faster than a line loop.
    The whole stream is read in both cases, so a partial member raises.
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
                # Raised as a read error so it is reported like a bad gzip: pigz exits non-zero on
                # a truncated member.
                raise OSError(
                    "pigz could not read %s: %s" % (path, stderr.decode("utf-8", "replace").strip())
                )
    ends_with_newline = last_byte == b"\n"
    if last_byte and not ends_with_newline:
        # Count the unterminated tail as a line, so a truncated record fails the multiple-of-four
        # check.
        count += 1
    return first.decode("utf-8", "replace"), count, ends_with_newline


def check_fastq(path: str, what: str = "FASTQ file") -> Optional[str]:
    """Validate that *path* looks like a complete FASTQ file."""
    problem = require_file(path, what)
    if problem:
        return problem
    # One streaming pass: validate the first line, count all lines (a complete FASTQ has a multiple
    # of four), and keep the last line. Reading the whole stream also detects a truncated gzip. A
    # tail window cannot be used, because its line count is arbitrary even for a valid file.
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
    # A header with no sequence after it means the file ends right after a '>' line. rpartition
    # also covers the single-line case: with no newline it returns ("", "", tail), so a file that
    # is only ">c1" is caught too.
    stripped = tail.rstrip("\n")
    if stripped.rpartition("\n")[2].startswith(">"):
        return "%s ends with a header line and no sequence, so it is truncated (%s)" % (
            _describe(path),
            what,
        )
    return None


def require_inputs(checks: Sequence[tuple], context: str = "") -> None:
    """Run several checks and raise :class:`InputProblem` listing every failure.

    *checks* is a sequence of ``(callable, path, description)``. All problems are reported in one
    run.
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
