"""Long reads: which flag each platform needs, and what the reads say about the flag you gave.

Every long-read tool in this pipeline takes a different name for the same four things. Flye wants
``--nano-hq``, minimap2 wants ``-x map-ont``, SemiBin2 wants ``--sequencing-type=long_read``, and
the error model behind all three differs between an ONT R9 read at 92% accuracy and a PacBio HiFi
read at 99.9%. Getting one of them wrong does not fail: it produces an assembly, a depth table or a
bin set that is quietly worse. :data:`READ_TYPES` is the one place that mapping lives.

**The declared type is checked against the reads.** A read set's length distribution and its
accuracy are observable, and the two things most likely to be wrong are declaring HiFi for ONT reads
(which tells Flye to expect 99.9% accuracy from 92% data) and pointing a long-read pipeline at short
reads by mistake. :func:`scan_reads` measures N50, median and the fraction of reads over 10 kb, and
:func:`check_type` says when the numbers do not fit the declaration. It warns rather than refusing -
a user with 5 kb HiFi reads knows something the distribution does not.

**What this module does not do: polishing.** Flye already polishes with the long reads themselves
(``--iterations``). Short-read polishing of a *metagenome* assembly is a separate decision with a
real failure mode - a polisher given reads from a close relative will happily rewrite one strain's
consensus into another's - and it belongs in its own module with its own measurement rather than
folded in here as a default.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

#: Reads below this are not long reads. Checked because the commonest long-read mistake is a path
#: that points at the short-read FASTQ, which every tool here will accept and assemble badly.
SHORT_READ_CEILING = 1_000

#: An N50 below this, with a long-read platform declared, is worth a word: it is either a very old
#: run or a mis-declared file. Not a refusal - amplicon and targeted long-read data really is short.
LOW_N50 = 3_000


@dataclass(frozen=True)
class ReadType:
    """One long-read platform, and the flag each tool wants for it."""

    name: str
    #: Flye's input flag, which encodes both platform and expected accuracy.
    flye_flag: str
    #: minimap2's preset, for mapping reads back to get coverage.
    minimap2_preset: str
    #: What SemiBin2 calls this; it distinguishes short from long and nothing finer.
    semibin2_type: str
    #: Roughly what per-base accuracy the platform delivers, for the sanity check and for the
    #: provenance record. Not used to filter reads.
    typical_accuracy: float
    description: str


#: The four long-read inputs worth distinguishing. Named as the platforms are, not as the tools'
#: flags are, so a user picks by what they sequenced.
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

#: The default, and the reason: R10.4 nanopore is what most new long-read metagenomes are, and
#: `--nano-hq` on genuinely raw reads is the less damaging of the two possible mistakes - Flye
#: handles a lower-accuracy read set under the high-accuracy model better than the reverse, where it
#: trusts errors it should have corrected.
DEFAULT_READ_TYPE = "ont-hq"


@dataclass(frozen=True)
class ReadStats:
    """What a scan of a long-read file found."""

    reads: int
    bases: int
    n50: int
    median: int
    longest: int
    #: Fraction of reads at or above 10 kb, which is what distinguishes a long-read run from a
    #: short-read one far more sharply than the mean does.
    fraction_over_10kb: float
    #: Reads below SHORT_READ_CEILING, which are not long reads whatever the file is called.
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
    """Measure a long-read file: N50, median, how much of it is genuinely long.

    A scan of the first *limit* reads rather than all of them. The decisions it feeds - is this the
    platform you said, is this a long-read file at all - do not change on read 100,001, and a full
    pass over a 40 GB nanopore run to answer them would cost more than the check is worth.
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
    """What the measured reads say about the declared platform. Lines to print, or none.

    Warnings rather than refusals throughout: every one of these distributions is something a real
    dataset can legitimately look like, and the point is that the user sees the mismatch before
    spending a day of assembly on it.
    """
    lines: List[str] = []
    if not stats.reads:
        return ["the read file appears to be empty"]
    if stats.looks_like_short_reads:
        lines.append(
            "N50 is %d bp, below %d: this looks like a short-read file. Every tool here will "
            "accept it and assemble it badly - metaFlye on short reads produces a fragmented "
            "assembly rather than an error" % (stats.n50, SHORT_READ_CEILING)
        )
        return lines
    if stats.n50 < LOW_N50:
        lines.append(
            "N50 is %d bp, which is short for a long-read run. Worth confirming this is the right "
            "file before assembling" % stats.n50
        )
    if stats.short_reads:
        lines.append(
            "%.1f%% of the reads are under %d bp. Flye discards short reads itself, so they cost "
            "nothing but they do mean the usable depth is lower than the file size suggests"
            % (100.0 * stats.short_reads / stats.reads, SHORT_READ_CEILING)
        )
    if read_type == "pacbio-hifi" and stats.fraction_over_10kb < 0.2:
        lines.append(
            "only %.0f%% of reads reach 10 kb, which is unusual for HiFi. If these are nanopore "
            "reads, --read-type pacbio-hifi tells Flye to expect 99.9%% accuracy from data that "
            "does not have it, and it will not correct errors it should"
            % (100.0 * stats.fraction_over_10kb)
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
        # The flags are in the row because they are what the run actually did, and two assemblies
        # made under different Flye models are not comparable.
        "flye_flag": spec.flye_flag if spec else "",
        "minimap2_preset": spec.minimap2_preset if spec else "",
        "notes": "; ".join(notes),
    }


#: Printed with every long-read assembly.
LONG_READ_CAVEAT = (
    "a long-read metagenome assembly is more contiguous and less accurate per base than a "
    "short-read one: Flye polishes with the long reads themselves, which cannot fix a systematic "
    "basecalling error. Gene calling and anything downstream of it inherits that - a frameshift "
    "from an indel truncates a protein, and no completeness estimate will notice"
)

#: Printed when both read kinds are given.
HYBRID_CAVEAT = (
    "with both read kinds, the long reads build the assembly and the short reads are used only "
    "where they are asked for. MetaWrap2 does not short-read-polish a metagenome assembly by "
    "default: a polisher given reads from a close relative rewrites one strain's consensus into "
    "another's, and that is a decision with a measurement behind it rather than a default"
)


def coassembly_inputs(paths: Sequence[str]) -> List[str]:
    """Long-read files to hand Flye as one co-assembly, after dropping what is not there.

    Flye takes several read files and treats them as one dataset, which is what makes a co-assembly
    a flag rather than a concatenation step - so unlike the short-read path, nothing has to be
    written to disk first.
    """
    return [path for path in paths if path and os.path.isfile(path)]
