"""Shared constants for metaWRAP.

Single source of truth for the magic numbers and literals that were previously
duplicated (and occasionally inconsistent) across the individual scripts. Import from
here rather than re-hardcoding a value in a module.
"""

from __future__ import annotations

# --- Bin quality scoring -------------------------------------------------------------
# metaWRAP ranks a bin by completeness penalised by contamination.
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
MIN_BIN_SIZE = 50_000        # 50 kb
MAX_BIN_SIZE = 20_000_000    # 20 Mb

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
