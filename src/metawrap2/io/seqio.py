"""Transparent handling of compressed and uncompressed sequence files.

Accepts ``.gz`` / ``.bz2``
input anywhere a plain FASTA/FASTQ was accepted, without forcing users to decompress
terabytes of data first. Compression is detected by magic bytes, not just extension,
so mislabeled files still work.
"""

from __future__ import annotations

import bz2
import gzip
import os
from typing import IO, Sequence, cast

from ..constants import BZIP2_MAGIC, FASTA_EXTENSIONS_WITH_GZ, GZIP_MAGIC, is_unbinned_name

#: Bases a binner can learn a composition from. Anything else in a record - N, an IUPAC ambiguity
#: code, a gap character - carries no k-mer signal, so a record made of nothing but these is not a
#: contig for binning purposes however long it is: its tetranucleotide vector is all zeros and no
#: read aligns to it. Shared by `binning` and `coverage` so the two filter the same thing.
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

    Everything that compares contigs between files must agree on this, because tools append
    their own annotations to the description field. Current metaBAT2, for example, writes

        >NODE_2_length_158684_cov_2.764955 total_depth=41.98 sample_depths=18.3,23.5,0.2

    while the assembly and every other binner write just the first token. MetaWrap2 used the
    *whole* header as the key in some places and the first token in others, so metaBAT2 bins
    matched nothing: bin_refinement measured 0%% overlap between a metaBAT2 bin and the very
    same bin from another binner (silently skipping consolidation), and quant_bins aborted
    with "None of the contigs in the assembly were present in the bin files".

    This is also what SAM/BAM reference names are, and what BLAST uses as the query id, so
    the first token is the only identity that is consistent end to end.
    """
    return header.split()[0] if header else header


def fasta_files(folder: str, extensions: Sequence[str] = FASTA_EXTENSIONS_WITH_GZ) -> list[str]:
    """Paths of the FASTA files in *folder*, sorted, compressed ones included.

    A bin folder is not guaranteed to hold only bins. Once a bin set has been through another
    module it also holds ``.stats.tsv`` tables and whole tool subdirectories (``checkm/``,
    ``gunc/``), and code that walked the folder with a bare ``os.listdir`` handed those to a FASTA
    reader: a stray table read as an empty bin, and a subdirectory raised ``IsADirectoryError``
    with no indication of which folder or which entry was at fault.

    The unbinned catch-all is a FASTA and still not a bin, so it is excluded too. `binning` writes the
    contigs no CONCOCT cluster claimed to ``concoct_bins/unbinned.fasta`` on purpose, so a user can see
    what was left over; every caller here is asking for candidate genomes, and each of them - CheckM,
    prokka, taxator, SPAdes, CoverM, sylph, the AMBER scorer - will answer about the leftovers as
    readily as about a genome. See :func:`metawrap2.constants.is_unbinned_name`.

    Sorted, because several callers number or order their output by the order bins are read in.
    """
    return sorted(
        os.path.join(folder, name)
        for name in os.listdir(folder)
        if name.endswith(tuple(extensions)) and not is_unbinned_name(name)
    )


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
