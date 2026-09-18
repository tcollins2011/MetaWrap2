"""Detection and validation of read-file layouts.

metaWRAP modules historically assumed exactly two uncompressed FASTQ files named
``*_1.fastq`` / ``*_2.fastq`` and crashed downstream (often deep inside a binner or
assembler) when given interleaved or single-end data. This module validates layout up
front and returns a clear structure, so failures happen early with an actionable message.

Covers issues #115 (interleaved not recognized) and #32 (single-end binning).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from .seqio import smart_open


@dataclass
class ReadSet:
    """A validated set of sequencing reads."""

    layout: str  # "paired" | "interleaved" | "single"
    files: List[str]

    @property
    def r1(self) -> Optional[str]:
        return self.files[0] if self.layout in ("paired", "interleaved", "single") else None

    @property
    def r2(self) -> Optional[str]:
        return self.files[1] if self.layout == "paired" and len(self.files) > 1 else None


def _first_read_name(path: str) -> Optional[str]:
    """Return the base name of the first read (without /1, /2 or ' 1:'/' 2:' suffixes)."""
    with smart_open(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                name = line[1:].strip()
                # Illumina old style: "name/1"; new style: "name 1:N:0:..."
                name = name.split(" ")[0]
                if name.endswith("/1") or name.endswith("/2"):
                    name = name[:-2]
                return name
    return None


def classify(
    files: List[str],
    interleaved: bool = False,
    single: bool = False,
) -> ReadSet:
    """Classify and validate a list of read files into a :class:`ReadSet`.

    Raises ``ValueError`` with a human-readable message on any inconsistency, rather than
    letting a downstream tool fail cryptically.
    """
    files = [f for f in files if f]
    if not files:
        raise ValueError("No read files were provided.")

    if single:
        if len(files) != 1:
            raise ValueError(
                f"--single expects exactly one read file, got {len(files)}: {files}"
            )
        return ReadSet(layout="single", files=files)

    if interleaved:
        if len(files) != 1:
            raise ValueError(
                f"--interleaved expects exactly one read file, got {len(files)}: {files}"
            )
        return ReadSet(layout="interleaved", files=files)

    if len(files) == 1:
        raise ValueError(
            "Only one read file was given. For single-end data pass --single; "
            "for interleaved paired data pass --interleaved."
        )
    if len(files) != 2:
        raise ValueError(
            f"Paired mode expects exactly two read files (R1 and R2), got {len(files)}: {files}"
        )

    # Sanity-check that the two files look like a mate pair (same first read name).
    n1 = _first_read_name(files[0])
    n2 = _first_read_name(files[1])
    if n1 is not None and n2 is not None and n1 != n2:
        raise ValueError(
            "Forward and reverse read files do not appear to be a matching pair "
            f"(first read names differ: '{n1}' vs '{n2}'). "
            "Check the file order (R1 then R2) or that they belong to the same sample."
        )

    return ReadSet(layout="paired", files=files)
