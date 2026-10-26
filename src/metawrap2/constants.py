"""Shared constants for MetaWrap2.

Numbers and literals used by more than one module. Import them from here instead of
hardcoding a value in a module.
"""

from __future__ import annotations

import os

# --- Bin quality scoring -------------------------------------------------------------
# A bin is ranked by completeness penalised by contamination.
CONTAMINATION_WEIGHT = 5.0
# Tiny tiebreaker added during dereplication so equal-scoring bins order deterministically.
DEREPLICATION_TIEBREAK = 1e-10

# Default quality thresholds (percent). Mirror the bin_refinement CLI defaults.
DEFAULT_MIN_COMPLETION = 70.0
DEFAULT_MAX_CONTAMINATION = 10.0

# --- Bin refinement ------------------------------------------------------------------
# Minimum length-overlap (percent) for two bins to be considered the same bin.
OVERLAP_THRESHOLD = 80.0
# bin_refinement ignores bins outside this size range (bytes) to save time.
MIN_BIN_SIZE = 50_000  # 50 kb
MAX_BIN_SIZE = 20_000_000  # 20 Mb

# --- Output file naming --------------------------------------------------------------
# Each output's final extension names its format. Defined once so all modules use the same suffix.
STATS_SUFFIX = ".stats.tsv"  # per-bin completeness/contamination table (was ".stats")

#: Columns :mod:`metawrap2.mimag` appends to a ``.stats.tsv``, after the fixed header and in this
#: order. Shared by the writer and the readers.
MIMAG_STATS_COLUMNS = ("mimag_tag", "mimag_rrna", "mimag_trnas")
CONTIGS_SUFFIX = ".contigs.tsv"  # contig -> bin membership table (was ".contigs")
TSV_SUFFIX = ".tsv"  # anything else tab-separated (was ".tab"/".tax")

# --- Bin file naming -----------------------------------------------------------------
# Bins from every binner are renamed ``bin_001.fasta``: zero-padded, counted from 1. The binners
# name their output differently (metaBAT2 counts from 1, MaxBin2 writes ``.001.fasta``, CONCOCT
# uses cluster ids). Bin sets are iterated in sorted filename order, so padding keeps that order
# the same as numeric order (unpadded, ``bin.10`` sorts before ``bin.2``).
BIN_PREFIX = "bin_"
BIN_NUMBER_WIDTH = 3  # bin_001 .. bin_999; widens automatically past 999 (see bin_filename)
BIN_EXTENSION = ".fasta"  # was ".fa"
UNBINNED_NAME = "unbinned" + BIN_EXTENSION


def bin_stem(number: int) -> str:
    """Return the bin's name without extension (``bin_001``), as used in .stats.tsv rows."""
    return "%s%0*d" % (BIN_PREFIX, BIN_NUMBER_WIDTH, number)


def bin_filename(number: int) -> str:
    """Return the filename for bin *number*, counting from 1: ``bin_001.fasta``.

    Numbers wider than the padding are not truncated (``bin_1000.fasta``), so names stay unique
    but sort order breaks past 999 bins.
    """
    return bin_stem(number) + BIN_EXTENSION


# --- File conventions ----------------------------------------------------------------
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fas")
FASTQ_EXTENSIONS = (".fastq", ".fq")

#: FASTA_EXTENSIONS plus their gzipped forms.
#:
#: Kept separate because gzip support depends on the reading tool. MetaWrap2's readers
#: (:func:`metawrap2.io.seqio.smart_open`) and CoverM accept it; prokka and CheckM1 do not. Modules
#: use this set only where every tool they call reads gzipped input.
FASTA_EXTENSIONS_WITH_GZ = FASTA_EXTENSIONS + tuple(e + ".gz" for e in FASTA_EXTENSIONS)


def bin_label(filename: str) -> str:
    """Return the name of a bin in a table: ``bin_001.fasta.gz`` -> ``bin_001``.

    ``quant_bins`` writes its abundance table from salmon or CoverM, and the two name gzipped bins
    differently; this gives both the same key. An unrecognised extension is kept, so
    ``E.coli.fasta`` becomes ``E.coli``, not ``E``.
    """
    base = os.path.basename(filename)
    stem = base.removesuffix(".gz")
    root, ext = os.path.splitext(stem)
    return root if ext.lower() in FASTA_EXTENSIONS else stem


def is_unbinned_name(filename: str) -> bool:
    """Return True if *filename* is the unbinned catch-all file, not a bin.

    `binning` writes contigs no CONCOCT cluster claimed to ``concoct_bins/unbinned.fasta``. Code
    that treats a directory of FASTAs as genomes must skip it. Compared on :func:`bin_label`, so
    renamed (``unbinned.fa``) and gzipped (``unbinned.fasta.gz``) forms also match.
    """
    return bin_label(filename) == bin_label(UNBINNED_NAME)


# Separator Binning_refiner uses to fold folder/bin/contig into one contig id.
REFINER_SEPARATOR = "__"


# --- The consolidated bin set --------------------------------------------------------
# bin_refinement writes one directory named for its thresholds. Callers build and recognise the
# name through the two functions below. (The benchmark once looked for a literal
# "metawrap2_bins", which never exists, and treated refined runs as unconsolidated.)
REFINED_BINS_PREFIX = "metawrap_"
REFINED_BINS_SUFFIX = "_bins"


def refined_bins_dirname(completeness: float, contamination: float) -> str:
    """Return the directory name for consolidated bins: ``metawrap_70_10_bins``."""
    return "%s%d_%d%s" % (
        REFINED_BINS_PREFIX,
        int(completeness),
        int(contamination),
        REFINED_BINS_SUFFIX,
    )


def is_refined_bins_dirname(name: str) -> bool:
    """Return True if *name* is a consolidated bin set directory, at any thresholds.

    Matched on the name pattern, since users can set any thresholds (for example -c 50 -x 5).
    """
    base = os.path.basename(name.rstrip("/"))
    if not (base.startswith(REFINED_BINS_PREFIX) and base.endswith(REFINED_BINS_SUFFIX)):
        return False
    middle = base[len(REFINED_BINS_PREFIX) : -len(REFINED_BINS_SUFFIX)]
    parts = middle.split("_")
    return len(parts) == 2 and all(p.isdigit() for p in parts)


# --- Compression ---------------------------------------------------------------------
GZIP_MAGIC = b"\x1f\x8b"
BZIP2_MAGIC = b"BZh"

# --- Console banner rendering (print_comment style) ----------------------------------
BANNER_WIDTH = 120
BANNER_MAX_LINE = 90
