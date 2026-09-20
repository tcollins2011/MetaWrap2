"""Shared constants for MetaWrap2.

Single source of truth for the magic numbers and literals that were previously
duplicated (and occasionally inconsistent) across the individual scripts. Import from
here rather than re-hardcoding a value in a module.
"""

from __future__ import annotations

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
# A file's final extension should say what format it is, so that editors, spreadsheets and
# other tools recognise it and a reader can tell at a glance. These are defined once so every
# module and helper agrees, rather than each hardcoding a suffix of its own.
STATS_SUFFIX = ".stats.tsv"  # per-bin completeness/contamination table (was ".stats")
CONTIGS_SUFFIX = ".contigs.tsv"  # contig -> bin membership table (was ".contigs")
TSV_SUFFIX = ".tsv"  # anything else tab-separated (was ".tab"/".tax")

# --- Bin file naming -----------------------------------------------------------------
# Bins are named ``bin_001.fasta``: one convention for every binner, zero-padded, counted
# from 1.
#
# Every binner names its output differently - metaBAT2 counts from 1, MaxBin2 writes
# zero-padded ``.001.fasta``, CONCOCT emits its own cluster ids - so without normalising, the
# same number means a different thing in each folder and ``bin.0`` exists in some and not others.
#
# Unpadded numbers also sort wrongly: ``bin.1``, ``bin.10``, ``bin.2``. That is not cosmetic,
# because bin sets are iterated in sorted filename order and refined bins are numbered in that
# order - so with ten or more bins the numbering depended on a lexicographic accident. Padding
# to a fixed width makes filename order and numeric order the same thing.
BIN_PREFIX = "bin_"
BIN_NUMBER_WIDTH = 3  # bin_001 .. bin_999; widens automatically past 999 (see bin_filename)
BIN_EXTENSION = ".fasta"  # was ".fa": an extension should name the format in full
UNBINNED_NAME = "unbinned" + BIN_EXTENSION


def bin_stem(number: int) -> str:
    """The bin's name without its extension: ``bin_001``. This is what .stats.tsv rows hold."""
    return "%s%0*d" % (BIN_PREFIX, BIN_NUMBER_WIDTH, number)


def bin_filename(number: int) -> str:
    """The filename for bin *number*, counting from 1: ``bin_001.fasta``.

    Numbers past the padding width are not truncated or wrapped - they simply get longer
    (``bin_1000.fasta``), which keeps names unique at the cost of sort order in the rare
    thousand-bin case.
    """
    return bin_stem(number) + BIN_EXTENSION


# --- File conventions ----------------------------------------------------------------
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fas")
FASTQ_EXTENSIONS = (".fastq", ".fq")
# Separator Binning_refiner uses to fold folder/bin/contig into one contig id.
REFINER_SEPARATOR = "__"

# --- Compression ---------------------------------------------------------------------
GZIP_MAGIC = b"\x1f\x8b"
BZIP2_MAGIC = b"BZh"

# --- Console banner rendering (print_comment style) ----------------------------------
BANNER_WIDTH = 120
BANNER_MAX_LINE = 90
