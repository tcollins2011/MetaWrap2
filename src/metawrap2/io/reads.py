"""Recognising read-file layouts and pairing up mates by filename.

One place decides what a read filename means, so every module agrees. Two jobs:

**Pairing.** Only ``name_1.fastq`` / ``name_2.fastq`` used to be recognised. Real data off
a sequencer is very often ``name_R1_001.fastq.gz`` / ``name_R2_001.fastq.gz``, so users had to
rename everything before they could start. :data:`MATE_PATTERNS` lists the conventions that are
understood; add a pattern here and every module accepts it.

**Layout validation.** Paired, interleaved, and single-end data are checked up front by
:func:`classify`, so a mistake is reported with an actionable message instead of failing deep
inside a binner or assembler hours later.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from .seqio import smart_open

#: FASTQ extensions, longest first so ``.fastq.gz`` is matched before ``.fastq``.
FASTQ_SUFFIXES = (".fastq.gz", ".fastq.bz2", ".fq.gz", ".fq.bz2", ".fastq", ".fq")

#: Mate-marker conventions, as (regex matching the R1 marker, template for the R2 marker).
#: Ordered most specific first. ``\1`` in the replacement keeps whatever the pattern captured
#: (the lane/chunk suffix in Illumina's ``_R1_001`` form).
MATE_PATTERNS: Sequence[Tuple[str, str]] = (
    (r"_R1(_\d+)$", r"_R2\1"),  # Illumina bcl2fastq: sample_R1_001 / sample_R2_001
    (r"_R1$", r"_R2"),  # sample_R1 / sample_R2
    (r"_r1$", r"_r2"),  # sample_r1 / sample_r2
    (r"_1$", r"_2"),  # sample_1 / sample_2   (the long-standing MetaWrap2 form)
    (r"\.1$", r".2"),  # sample.1 / sample.2
    (r"_R1\.", r"_R2."),  # defensive: marker not at the end
)


def split_fastq_suffix(path: str) -> Tuple[str, str]:
    """Split a read filename into (stem, suffix), e.g. ``a_R1_001.fastq.gz`` -> the two parts.

    The suffix keeps its compression extension, so callers can reconstruct a mate's name.
    Returns ``(path, "")`` when it is not a FASTQ name at all.
    """
    base = os.path.basename(path)
    for suffix in FASTQ_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)], suffix
    return base, ""


def is_fastq(path: str) -> bool:
    """True if *path* looks like a FASTQ file (optionally compressed)."""
    return bool(split_fastq_suffix(path)[1])


def mate_of(path: str) -> Optional[str]:
    """The R2 filename matching R1 *path*, or None if *path* is not a recognised R1 name.

    Only the *filename* is rewritten; the directory is preserved.
    """
    stem, suffix = split_fastq_suffix(path)
    if not suffix:
        return None
    for pattern, replacement in MATE_PATTERNS:
        if re.search(pattern, stem):
            mate_stem = re.sub(pattern, replacement, stem)
            return os.path.join(os.path.dirname(path), mate_stem + suffix)
    return None


def sample_name(path: str) -> str:
    """Sample name for a read file: its stem with the mate marker removed.

    ``ERR011347_1.fastq.gz`` and ``ERR011347_R1_001.fastq.gz`` both give ``ERR011347``, so a
    sample keeps one name whichever convention its files use.
    """
    stem, suffix = split_fastq_suffix(path)
    if not suffix:
        stem = os.path.splitext(os.path.basename(path))[0]
    for pattern, _replacement in MATE_PATTERNS:
        if re.search(pattern, stem):
            return re.sub(pattern, "", stem)
    return stem


def describe_conventions() -> str:
    """Human-readable list of the accepted pair-naming conventions, for error messages."""
    examples = [
        "sample_1/sample_2",
        "sample_R1/sample_R2",
        "sample_R1_001/sample_R2_001",
        "sample.1/sample.2",
    ]
    return ", ".join(examples)


@dataclass
class ReadSet:
    """A validated set of sequencing reads."""

    layout: str  # "paired" | "interleaved" | "single"
    files: List[str]

    @property
    def r1(self) -> Optional[str]:
        return self.files[0] if self.files else None

    @property
    def r2(self) -> Optional[str]:
        return self.files[1] if self.layout == "paired" and len(self.files) > 1 else None

    @property
    def sample(self) -> str:
        return sample_name(self.files[0]) if self.files else ""


def _first_read_name(path: str) -> Optional[str]:
    """Base name of the first read, without /1, /2 or ' 1:'/' 2:' mate suffixes."""
    with smart_open(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                name = line[1:].strip().split(" ")[0]
                if name.endswith(("/1", "/2")):
                    name = name[:-2]
                return name
    return None


def classify(files: List[str], interleaved: bool = False, single: bool = False) -> ReadSet:
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
                "--single-end expects exactly one read file, got %d: %s" % (len(files), files)
            )
        return ReadSet(layout="single", files=files)

    if interleaved:
        if len(files) != 1:
            raise ValueError(
                "--interleaved expects exactly one read file, got %d: %s" % (len(files), files)
            )
        return ReadSet(layout="interleaved", files=files)

    if len(files) == 1:
        raise ValueError(
            "Only one read file was given. For single-end data pass --single-end; "
            "for interleaved paired data pass --interleaved."
        )
    if len(files) != 2:
        raise ValueError(
            "Paired mode expects exactly two read files (R1 and R2), got %d: %s"
            % (len(files), files)
        )

    # Sanity-check that the two files look like a mate pair (same first read name).
    n1 = _first_read_name(files[0])
    n2 = _first_read_name(files[1])
    if n1 is not None and n2 is not None and n1 != n2:
        raise ValueError(
            "Forward and reverse read files do not appear to be a matching pair "
            "(first read names differ: '%s' vs '%s'). Check the file order (R1 then R2) or "
            "that they belong to the same sample." % (n1, n2)
        )

    return ReadSet(layout="paired", files=files)


def collect_pairs(reads: List[str]) -> List[Tuple[str, str, str]]:
    """Pair up a positional list of read files into ``(sample, r1, r2)`` triples.

    Accepts any convention in :data:`MATE_PATTERNS`, and accepts being handed the R2 files
    too (they are recognised as mates of an R1 already seen, not treated as new samples).
    """
    pairs: List[Tuple[str, str, str]] = []
    seen_r2 = set()
    for path in reads:
        if not is_fastq(path):
            continue
        mate = mate_of(path)
        if mate is None or os.path.abspath(path) in seen_r2:
            continue
        if not os.path.isfile(path):
            raise ValueError("%s does not exist" % path)
        if not os.path.isfile(mate):
            raise ValueError("expected mate %s for %s but it does not exist" % (mate, path))
        seen_r2.add(os.path.abspath(mate))
        pairs.append((sample_name(path), path, mate))

    if not pairs:
        raise ValueError(
            "No paired read files found among: %s\n"
            "Recognised pair-naming conventions are: %s (.gz/.bz2 allowed).\n"
            "For single-end data use --single-end; for interleaved data use --interleaved."
            % (", ".join(os.path.basename(r) for r in reads) or "<nothing>", describe_conventions())
        )
    return pairs


def collect_singles(reads: List[str]) -> List[Tuple[str, str]]:
    """Return ``(sample, path)`` for each FASTQ in *reads* (single-end or interleaved)."""
    out: List[Tuple[str, str]] = []
    for path in reads:
        if not is_fastq(path):
            continue
        if not os.path.isfile(path):
            raise ValueError("%s does not exist" % path)
        out.append((sample_name(path), path))
    if not out:
        raise ValueError(
            "No read files found among: %s (expected %s)"
            % (
                ", ".join(os.path.basename(r) for r in reads) or "<nothing>",
                "/".join(FASTQ_SUFFIXES),
            )
        )
    return out


def deinterleave(path: str, out_r1: str, out_r2: str) -> int:
    """Split an interleaved FASTQ into two mate files. Returns the number of pairs written.

    Interleaved input is handled by splitting it up front and then running the ordinary paired
    path, rather than by teaching every tool about a third layout: Trim Galore has no
    interleaved mode at all, and treating interleaved data as single-end silently destroys the
    pairing that the assembler and binners depend on.

    Records are taken strictly two at a time, which is what "interleaved" means. A file with
    an odd number of records, or whose adjacent records are not mates, is a real problem worth
    reporting rather than quietly mis-pairing.
    """

    def _records(handle):
        block = []
        for line in handle:
            block.append(line if line.endswith("\n") else line + "\n")
            if len(block) == 4:
                yield block
                block = []
        if block:
            raise ValueError("%s ends in a truncated FASTQ record" % path)

    def _base(header: str) -> str:
        name = header[1:].strip().split(" ")[0]
        return name[:-2] if name.endswith(("/1", "/2")) else name

    pairs = 0
    with smart_open(path, "rt") as fh, open(out_r1, "w") as o1, open(out_r2, "w") as o2:
        records = _records(fh)
        for first in records:
            try:
                second = next(records)
            except StopIteration:
                raise ValueError(
                    "%s has an odd number of records, so it is not interleaved (the last read "
                    "has no mate). Did you mean --single-end?" % path
                ) from None
            if _base(first[0]) != _base(second[0]):
                raise ValueError(
                    "%s does not look interleaved: adjacent reads '%s' and '%s' are not mates. "
                    "Did you mean --single-end, or are these two concatenated files?"
                    % (path, _base(first[0]), _base(second[0]))
                )
            o1.writelines(first)
            o2.writelines(second)
            pairs += 1
    if pairs == 0:
        raise ValueError("%s contained no FASTQ records" % path)
    return pairs
