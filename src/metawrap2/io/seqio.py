"""Transparent handling of compressed and uncompressed sequence files.

Addresses the most-requested metaWRAP feature (issue #159): accept ``.gz`` / ``.bz2``
input anywhere a plain FASTA/FASTQ was accepted, without forcing users to decompress
terabytes of data first. Compression is detected by magic bytes, not just extension,
so mislabeled files still work.
"""

from __future__ import annotations

import bz2
import gzip
import io
from typing import IO

from ..constants import BZIP2_MAGIC, GZIP_MAGIC


def detect_compression(path: str) -> str:
    """Return ``"gzip"``, ``"bzip2"`` or ``"none"`` by inspecting the file header."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(3)
    except OSError:
        # Fall back to extension if we cannot read (e.g. writing a new file).
        lower = path.lower()
        if lower.endswith(".gz"):
            return "gzip"
        if lower.endswith(".bz2"):
            return "bzip2"
        return "none"
    if head[:2] == GZIP_MAGIC:
        return "gzip"
    if head[:3] == BZIP2_MAGIC:
        return "bzip2"
    return "none"


def _compression_for_write(path: str) -> str:
    lower = path.lower()
    if lower.endswith(".gz"):
        return "gzip"
    if lower.endswith(".bz2"):
        return "bzip2"
    return "none"


def smart_open(path: str, mode: str = "rt", **kwargs) -> IO:
    """Open *path* for reading or writing, transparently (de)compressing.

    For reads, compression is detected from the file contents. For writes, it is chosen
    from the extension (``.gz`` -> gzip, ``.bz2`` -> bzip2, else plain). Text mode is the
    default so callers deal in ``str``.
    """
    reading = "r" in mode
    comp = detect_compression(path) if reading else _compression_for_write(path)

    if comp == "gzip":
        return gzip.open(path, mode, **kwargs)
    if comp == "bzip2":
        return bz2.open(path, mode, **kwargs)
    return io.open(path, mode, **kwargs)


def iter_fasta(path: str):
    """Yield ``(header, sequence)`` tuples from a (optionally compressed) FASTA file.

    ``header`` excludes the leading ``>`` and trailing newline but is otherwise the exact
    identifier line. Contig identity is treated as the whole string; no field is parsed
    out of it (the legacy code assumed SPAdes-style ``NODE_x_length_...`` headers, which
    silently corrupted results on other assemblers).
    """
    header = None
    seq_parts: list[str] = []
    with smart_open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)
