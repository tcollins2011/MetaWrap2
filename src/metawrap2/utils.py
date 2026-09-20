"""Small utility functions shared across MetaWrap2 modules.

Kept deliberately dependency-light so any module or helper script can import from here
without pulling in heavy libraries.
"""

from __future__ import annotations

import os

from .constants import CONTAMINATION_WEIGHT, DEREPLICATION_TIEBREAK


def quality_score(completeness: float, contamination: float, tiebreak: float = 0.0) -> float:
    """Bin quality score: completeness - 5*contamination (+ tiny tiebreaker).

    The tiebreaker (typically genome size) only separates otherwise-equal bins and does
    not meaningfully shift ranking.
    """
    return completeness - CONTAMINATION_WEIGHT * contamination + DEREPLICATION_TIEBREAK * tiebreak


def bin_name_from_filename(filename: str) -> str:
    """Return a bin's name by stripping the directory and final extension.

    e.g. ``/path/bin.3.fa`` -> ``bin.3`` and ``some/dir/maxbin.007.fasta`` -> ``maxbin.007``.
    """
    base = os.path.basename(filename)
    stem, _ext = os.path.splitext(base)
    return stem


def strip_trailing_slash(path: str) -> str:
    """Remove a single trailing slash, matching the legacy scripts' path handling."""
    return path.removesuffix("/")
