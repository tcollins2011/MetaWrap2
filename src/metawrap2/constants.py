"""Shared constants for MetaWrap2.

Single source of truth for the magic numbers and literals that were previously
duplicated (and occasionally inconsistent) across the individual scripts. Import from
here rather than re-hardcoding a value in a module.
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

#: The same extensions plus their gzipped forms - a superset, so anything matching FASTA_EXTENSIONS
#: matches this too.
#:
#: Deliberately not folded into FASTA_EXTENSIONS, because whether a gzipped bin is acceptable is a
#: question about the tool that will read it, not about the file. MetaWrap2's own readers
#: (:func:`metawrap2.io.seqio.smart_open` and everything built on it) do not care about compression,
#: and neither does CoverM - but prokka and CheckM1 do, so for the modules that hand a bin straight
#: to one of those, accepting ``.fasta.gz`` at validation would only move the failure from the first
#: second of the run into the middle of it. Modules opt in where compression is known to work.
FASTA_EXTENSIONS_WITH_GZ = FASTA_EXTENSIONS + tuple(e + ".gz" for e in FASTA_EXTENSIONS)


def bin_label(filename: str) -> str:
    """The name a bin goes by in a table: ``bin_001.fasta.gz`` -> ``bin_001``.

    One function because two tables have to agree. ``quant_bins`` can write its abundance table
    from salmon or from CoverM, and the two tools name a genome differently once compression is
    involved: stripping only the last dotted component leaves ``bin_001.fasta`` for a gzipped bin,
    while CoverM (given ``-x fasta.gz``) calls the same bin ``bin_001``. Tables that key the same
    bin two ways cannot be joined, compared, or read by the same script.

    An unrecognised extension is left alone rather than guessed at: a bin called ``E.coli.fasta``
    must not become ``E``.
    """
    base = os.path.basename(filename)
    stem = base.removesuffix(".gz")
    root, ext = os.path.splitext(stem)
    return root if ext.lower() in FASTA_EXTENSIONS else stem


def is_unbinned_name(filename: str) -> bool:
    """True if *filename* is the unbinned catch-all rather than a bin.

    A bin set may legitimately contain one: `binning` writes the contigs no CONCOCT cluster claimed to
    ``concoct_bins/unbinned.fasta`` on purpose, so a user can see what was left over. Nothing that
    treats a directory of FASTAs as a set of candidate genomes should pick it up - it is not a genome,
    and every tool that scores, refines, annotates, classifies or reassembles one will do it anyway.

    Compared on :func:`bin_label` rather than the filename, so it still holds once an extension has
    been normalised (``unbinned.fa`` -> ``unbinned.fasta``) or a bin set has been gzipped
    (``unbinned.fasta.gz``) - each of which is a point at which a check against the full name would
    quietly start letting the file through.
    """
    return bin_label(filename) == bin_label(UNBINNED_NAME)


# Separator Binning_refiner uses to fold folder/bin/contig into one contig id.
REFINER_SEPARATOR = "__"


# --- The consolidated bin set --------------------------------------------------------
# bin_refinement's answer is one directory, named for the thresholds it was cut at, so a result
# carries the thresholds that produced it rather than relying on a note beside it.
#
# Written and recognised through the two functions below rather than by each caller spelling the
# name out. The benchmark looked for a literal "metawrap2_bins", which bin_refinement has never
# produced - so a refined run looked like three raw overlapping bin sets with nothing consolidated,
# and the harness declined to rank the one answer that was sitting right there.
REFINED_BINS_PREFIX = "metawrap_"
REFINED_BINS_SUFFIX = "_bins"


def refined_bins_dirname(completeness: float, contamination: float) -> str:
    """The directory bin_refinement leaves its consolidated bins in: ``metawrap_70_10_bins``."""
    return "%s%d_%d%s" % (
        REFINED_BINS_PREFIX,
        int(completeness),
        int(contamination),
        REFINED_BINS_SUFFIX,
    )


def is_refined_bins_dirname(name: str) -> bool:
    """Is *name* a consolidated bin set, at whatever thresholds it was cut at?

    Matched on shape rather than against a list, because the thresholds are the user's: -c 50 -x 5
    is as valid as the defaults, and a rule that only knew the common ones would silently treat an
    unusual run as unconsolidated.
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
