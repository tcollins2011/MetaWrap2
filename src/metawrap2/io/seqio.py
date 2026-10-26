"""Transparent handling of compressed and uncompressed sequence files.

Accepts ``.gz`` / ``.bz2`` input anywhere a plain FASTA/FASTQ is accepted, so input does not
have to be decompressed first. Compression is detected by magic bytes, so mislabeled files work.
"""

from __future__ import annotations

import bz2
import gzip
import os
from typing import IO, Sequence, cast

from ..constants import BZIP2_MAGIC, FASTA_EXTENSIONS_WITH_GZ, GZIP_MAGIC, is_unbinned_name

#: Bases that carry composition signal for a binner. A record with none of them (only N, IUPAC
#: ambiguity codes or gaps) has an all-zero tetranucleotide vector and no read aligns to it.
#: Shared by `binning` and `coverage` so both filter the same records.
INFORMATIVE_BASES = frozenset("ACGTUacgtu")


def is_binnable(seq: str) -> bool:
    """Whether *seq* carries any signal a binner could use. Empty and all-N sequences do not."""
    return bool(INFORMATIVE_BASES.intersection(seq))


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
        # gzip/bz2/io.open all return file-like objects; the union of their concrete types
        # is not IO[Any] to mypy, but every caller only uses the IO interface.
        return cast(IO, gzip.open(path, mode, **kwargs))
    if comp == "bzip2":
        return cast(IO, bz2.open(path, mode, **kwargs))
    return open(path, mode, **kwargs)


def contig_id(header: str) -> str:
    """The identity of a contig: its header up to the first whitespace.

    All code that compares contigs between files uses this, because tools append their own
    annotations to the description field. Current metaBAT2, for example, writes

        >NODE_2_length_158684_cov_2.764955 total_depth=41.98 sample_depths=18.3,23.5,0.2

    while the assembly and other binners write only the first token. Keying on the whole header
    in some places made metaBAT2 bins match nothing: bin_refinement found 0% overlap with the
    same bin from another binner, and quant_bins aborted with "None of the contigs in the
    assembly were present in the bin files".

    The first token is also the SAM/BAM reference name and the BLAST query id.
    """
    return header.split()[0] if header else header


def fasta_files(folder: str, extensions: Sequence[str] = FASTA_EXTENSIONS_WITH_GZ) -> list[str]:
    """Paths of the FASTA files in *folder*, sorted, compressed ones included.

    A bin folder can also hold ``.stats.tsv`` tables and tool subdirectories (``checkm/``,
    ``gunc/``) once another module has run on it; only FASTA files are returned.

    The unbinned catch-all (``concoct_bins/unbinned.fasta``, the contigs no CONCOCT cluster
    claimed) is excluded too, since callers (CheckM, prokka, taxator, SPAdes, CoverM, sylph, the
    AMBER scorer) want candidate genomes. See :func:`metawrap2.constants.is_unbinned_name`.

    Sorted, because several callers number or order their output by read order.
    """
    return sorted(
        os.path.join(folder, name)
        for name in os.listdir(folder)
        if name.endswith(tuple(extensions)) and not is_unbinned_name(name)
    )


def iter_fasta(path: str):
    """Yield ``(header, sequence)`` tuples from a (optionally compressed) FASTA file.

    ``header`` excludes the leading ``>`` and trailing newline but is otherwise the exact
    identifier line. No field is parsed out of it (the legacy code assumed SPAdes-style
    ``NODE_x_length_...`` headers, which broke results from other assemblers).
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
