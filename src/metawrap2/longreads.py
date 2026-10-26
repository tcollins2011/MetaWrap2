"""Long reads: the flag each platform needs per tool, and a check of the declared platform.

Each long-read tool names the platform differently: Flye ``--nano-hq``, minimap2 ``-x map-ont``,
SemiBin2 ``--sequencing-type=long_read``. The error model differs between an ONT R9 read at 92%
accuracy and a PacBio HiFi read at 99.9%. A wrong flag does not fail; it gives a worse assembly,
depth table or bin set. :data:`READ_TYPES` holds the mapping.

The declared type is checked against the reads. The likeliest mistakes are declaring HiFi for ONT
reads (Flye then expects 99.9% accuracy from 92% data) and passing short reads by mistake.
:func:`scan_reads` measures N50, median and the fraction of reads over 10 kb, and
:func:`check_type` reports when these do not fit the declaration. It only warns, since some real
data (e.g. 5 kb HiFi reads) will not fit the typical distribution.

No polishing is done here. Flye polishes with the long reads (``--iterations``). Short-read
polishing of a metagenome assembly can rewrite one strain's consensus into a close relative's, so
it would need its own module and benchmark.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

#: Reads below this are not long reads. A common mistake is passing the short-read FASTQ, which
#: every tool here accepts and assembles badly.
SHORT_READ_CEILING = 1_000

#: An N50 below this with a long-read platform declared triggers a warning (a very old run or a
#: mis-declared file). Amplicon and targeted long-read data can be this short.
LOW_N50 = 3_000


@dataclass(frozen=True)
class ReadType:
    """One long-read platform, and the flag each tool wants for it."""

    name: str
    #: Flye's input flag, which encodes both platform and expected accuracy.
    flye_flag: str
    #: minimap2's preset, for mapping reads back to get coverage.
    minimap2_preset: str
    #: SemiBin2's sequencing type; it only distinguishes short from long.
    semibin2_type: str
    #: Approximate per-base accuracy, for the sanity check and provenance. Not used to filter.
    typical_accuracy: float
    description: str


#: The four long-read input types, named by platform so a user picks by what they sequenced.
READ_TYPES: Dict[str, ReadType] = {
    "ont": ReadType(
        name="ont",
        flye_flag="--nano-raw",
        minimap2_preset="map-ont",
        semibin2_type="long_read",
        typical_accuracy=0.92,
        description="Oxford Nanopore, pre-R10 or fast basecalling (~92% accuracy)",
    ),
    "ont-hq": ReadType(
        name="ont-hq",
        flye_flag="--nano-hq",
        minimap2_preset="map-ont",
        semibin2_type="long_read",
        typical_accuracy=0.99,
        description="Oxford Nanopore R10.4+ or super-accurate basecalling (>=99%)",
    ),
    "pacbio": ReadType(
        name="pacbio",
        flye_flag="--pacbio-raw",
        minimap2_preset="map-pb",
        semibin2_type="long_read",
        typical_accuracy=0.88,
        description="PacBio CLR, continuous long reads",
    ),
    "pacbio-hifi": ReadType(
        name="pacbio-hifi",
        flye_flag="--pacbio-hifi",
        minimap2_preset="map-hifi",
        semibin2_type="long_read",
        typical_accuracy=0.999,
        description="PacBio HiFi / CCS (>=99.9% accuracy)",
    ),
}

#: Most new long-read metagenomes are R10.4 nanopore. Flye handles lower-accuracy reads under
#: `--nano-hq` better than high-accuracy reads under the raw model, where it keeps errors it should
#: correct.
DEFAULT_READ_TYPE = "ont-hq"


@dataclass(frozen=True)
class ReadStats:
    """What a scan of a long-read file found."""

    reads: int
    bases: int
    n50: int
    median: int
    longest: int
    #: Fraction of reads at or above 10 kb; separates long-read from short-read runs better than
    #: the mean.
    fraction_over_10kb: float
    #: Reads below SHORT_READ_CEILING.
    short_reads: int

    @property
    def mean_length(self) -> Optional[float]:
        return self.bases / self.reads if self.reads else None

    @property
    def looks_like_short_reads(self) -> bool:
        """Whether this file is a short-read set that was passed by mistake."""
        return bool(self.reads) and self.n50 < SHORT_READ_CEILING


def _iter_lengths(path: str, limit: int = 0):
    """Read lengths from a FASTA or FASTQ file, gzipped or not."""
    opener: object = open
    if path.endswith(".gz"):
        import gzip

        opener = gzip.open
    seen = 0
    with opener(path, "rt") as handle:  # type: ignore[operator]
        first = handle.readline()
        if not first:
            return
        fastq = first.startswith("@")
        if fastq:
            while True:
                sequence = handle.readline().rstrip("\n")
                handle.readline()
                handle.readline()
                if not sequence:
                    return
                yield len(sequence)
                seen += 1
                if limit and seen >= limit:
                    return
                if not handle.readline():
                    return
        else:
            length = 0
            for line in handle:
                if line.startswith(">"):
                    if length:
                        yield length
                        seen += 1
                        if limit and seen >= limit:
                            return
                    length = 0
                else:
                    length += len(line.strip())
            if length:
                yield length


def scan_reads(path: str, limit: int = 100_000) -> ReadStats:
    """Measure a long-read file: N50, median, fraction of reads over 10 kb.

    Only the first *limit* reads are scanned, which is enough for the platform check and avoids a
    full pass over e.g. a 40 GB nanopore run.
    """
    lengths = sorted(_iter_lengths(path, limit), reverse=True)
    if not lengths:
        return ReadStats(0, 0, 0, 0, 0, 0.0, 0)
    total = sum(lengths)
    half = total / 2
    running = 0
    n50 = lengths[-1]
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return ReadStats(
        reads=len(lengths),
        bases=total,
        n50=n50,
        median=lengths[len(lengths) // 2],
        longest=lengths[0],
        fraction_over_10kb=sum(1 for length in lengths if length >= 10_000) / len(lengths),
        short_reads=sum(1 for length in lengths if length < SHORT_READ_CEILING),
    )


def check_type(stats: ReadStats, read_type: str) -> List[str]:
    """Warnings about mismatches between the measured reads and the declared platform.

    Only warnings, since real datasets can have any of these distributions. They are shown before
    assembly starts.
    """
    lines: List[str] = []
    if not stats.reads:
        return ["the read file appears to be empty"]
    if stats.looks_like_short_reads:
        lines.append(
            "N50 is %d bp, below %d: this looks like a short-read file. metaFlye accepts short "
            "reads without an error and produces a fragmented assembly"
            % (stats.n50, SHORT_READ_CEILING)
        )
        return lines
    if stats.n50 < LOW_N50:
        lines.append(
            "N50 is %d bp, which is short for a long-read run. Check that this is the right "
            "file before assembling" % stats.n50
        )
    if stats.short_reads:
        lines.append(
            "%.1f%% of the reads are under %d bp. Flye discards short reads, so the usable "
            "depth is lower than the file size suggests"
            % (100.0 * stats.short_reads / stats.reads, SHORT_READ_CEILING)
        )
    if read_type == "pacbio-hifi" and stats.fraction_over_10kb < 0.2:
        lines.append(
            "only %.0f%% of reads reach 10 kb, which is unusual for HiFi. If these are nanopore "
            "reads, --read-type pacbio-hifi makes Flye expect 99.9%% accuracy, and it will leave "
            "errors uncorrected" % (100.0 * stats.fraction_over_10kb)
        )
    return lines


def resolve(read_type: str) -> ReadType:
    """The flags for a declared platform, or a clear error listing what is known."""
    if read_type not in READ_TYPES:
        raise KeyError(
            "unknown read type %r. Known: %s"
            % (
                read_type,
                ", ".join(
                    "%s (%s)" % (name, spec.description)
                    for name, spec in sorted(READ_TYPES.items())
                ),
            )
        )
    return READ_TYPES[read_type]


#: Columns of the per-sample read report.
READ_REPORT_COLUMNS: Tuple[str, ...] = (
    "sample",
    "read_type",
    "reads_scanned",
    "bases_scanned",
    "n50",
    "median_length",
    "longest",
    "percent_over_10kb",
    "percent_under_1kb",
    "flye_flag",
    "minimap2_preset",
    "notes",
)


def read_report_row(
    sample: str, read_type: str, stats: ReadStats, notes: Sequence[str] = ()
) -> Dict[str, str]:
    spec = READ_TYPES.get(read_type)
    return {
        "sample": sample,
        "read_type": read_type,
        "reads_scanned": str(stats.reads),
        "bases_scanned": str(stats.bases),
        "n50": str(stats.n50),
        "median_length": str(stats.median),
        "longest": str(stats.longest),
        "percent_over_10kb": "%.1f" % (100.0 * stats.fraction_over_10kb),
        "percent_under_1kb": (
            "%.1f" % (100.0 * stats.short_reads / stats.reads) if stats.reads else ""
        ),
        # Record the flags used: assemblies made under different Flye models are not comparable.
        "flye_flag": spec.flye_flag if spec else "",
        "minimap2_preset": spec.minimap2_preset if spec else "",
        "notes": "; ".join(notes),
    }


#: Printed with every long-read assembly.
LONG_READ_CAVEAT = (
    "a long-read metagenome assembly is more contiguous but less accurate per base than a "
    "short-read one. Flye polishes with the long reads, which does not fix systematic "
    "basecalling errors. These errors carry into gene calls: a frameshift from an indel "
    "truncates a protein, and completeness estimates do not detect it"
)

#: Printed when both read kinds are given.
HYBRID_CAVEAT = (
    "with both read kinds, the long reads build the assembly and the short reads are used only "
    "where an option asks for them. MetaWrap2 does not polish a metagenome assembly with short "
    "reads by default: given reads from a close relative, a polisher can rewrite one strain's "
    "consensus into another's"
)


def coassembly_inputs(paths: Sequence[str]) -> List[str]:
    """Existing long-read files to pass to Flye as one co-assembly.

    Flye treats several read files as one dataset, so unlike the short-read path nothing is
    concatenated on disk.
    """
    return [path for path in paths if path and os.path.isfile(path)]
