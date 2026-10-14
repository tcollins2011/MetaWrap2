"""How much of an assembly reached a bin, and what a bin set does not say about itself.

A bin set reports what it found. The two things it does not report are what it missed and how much
of what it found is one organism, and both have been measured by other people on large numbers of
samples:

* **Most assembled sequence never reaches a MAG.** Across 92,187 metagenomes in SPIRE, only about
  10% of assembled base pairs end up in a medium- or high-quality MAG, and the unbinned remainder is
  not a random sample of the data - it holds an estimated 145 predicted undescribed phyla. A bin set
  presented without that fraction invites the reading that it describes the sample.
* **"High-quality" MAGs are chimeric more often than their scores suggest.** Treichel et al. 2026
  found only 54.5-81.8% of high-quality MAGs represent a single strain, with a mean of 2.3 reference
  genomes per MAG. CheckM2's completeness and contamination do not measure this; GUNC measures part
  of it, at the taxonomic level, and this pipeline already runs GUNC.

So this computes the first and states the second. Neither changes a bin - they change what a reader
concludes from one, which is the point. Both are cheap: the binned fraction is a pass over sequence
already on disk, and the chimerism caveat is a sentence attached to a number GUNC already produced.

The unit is **base pairs, not contigs**. A bin set can hold 2% of the contigs and 60% of the
sequence, because binners discard short contigs first and short contigs are numerous. Reporting
contigs would make every bin set look worse than it is and would not answer the question anyone has,
which is how much of the sample this describes.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .io.seqio import contig_id, iter_fasta

#: Printed with every binned-fraction figure. The number invites a conclusion it cannot support
#: without this, so the two travel together.
UNBINNED_CAVEAT = (
    "the unbinned fraction is not a random sample of the assembly: binners discard what they "
    "cannot place, which is biased toward short contigs, low coverage and novel lineages. Across "
    "92,187 metagenomes (SPIRE) only ~10% of assembled bp reach a medium- or high-quality MAG"
)

#: Printed with every quality-tier count. GUNC catches some of this and CheckM2 none of it.
CHIMERISM_CAVEAT = (
    "a 70/10 or 90/5 score says a bin has the markers of one genome, not that it is one: of MAGs "
    "called high-quality, 54.5-81.8% represent a single strain (mean 2.3 reference genomes per MAG; "
    "Treichel et al. 2026). GUNC detects the taxonomically incoherent subset, which is why it runs "
    "by default"
)


@dataclass(frozen=True)
class Recovery:
    """What share of an assembly a bin set accounts for."""

    assembly_contigs: int
    assembly_bp: int
    binned_contigs: int
    binned_bp: int
    #: Only the bins that pass a quality threshold, which is the number worth quoting: sequence in a
    #: bin nobody would keep has not really been recovered. None when no such set was given, which
    #: is distinct from zero - "not measured" and "none of it is clean" are different statements and
    #: reporting 0.0% for the first reads as the second.
    clean_contigs: Optional[int] = None
    clean_bp: Optional[int] = None

    @property
    def binned_fraction(self) -> float:
        return self.binned_bp / float(self.assembly_bp) if self.assembly_bp else 0.0

    @property
    def clean_fraction(self) -> Optional[float]:
        if self.clean_bp is None or not self.assembly_bp:
            return None
        return self.clean_bp / float(self.assembly_bp)

    @property
    def unbinned_bp(self) -> int:
        return max(self.assembly_bp - self.binned_bp, 0)

    def describe(self) -> str:
        """One line for the end of a run."""
        clean = self.clean_fraction
        quality = (
            ""
            if clean is None
            else " (%.1f%% in a bin that passes the quality thresholds)" % (100 * clean)
        )
        return "%.1f%% of assembled bp are in a bin%s; %s bp unbinned" % (
            100 * self.binned_fraction,
            quality,
            f"{self.unbinned_bp:,}",
        )


def contig_lengths(path: str) -> Dict[str, int]:
    """contig -> length, for one fasta."""
    return {contig_id(header): len(seq) for header, seq in iter_fasta(path)}


def binned_contigs(
    bins_dir: str, extensions: Sequence[str] = (".fa", ".fasta", ".fna")
) -> Set[str]:
    """Every contig name appearing in any bin in *bins_dir*.

    A set, because a contig in two bins is still one contig of the assembly. That happens: bin sets
    from different binners overlap, and a consolidated set can retain a contig that a size filter
    split out of a restored input.
    """
    found: Set[str] = set()
    if not os.path.isdir(bins_dir):
        return found
    for name in sorted(os.listdir(bins_dir)):
        if not name.endswith(tuple(extensions)):
            continue
        for header, _seq in iter_fasta(os.path.join(bins_dir, name)):
            found.add(contig_id(header))
    return found


def measure(
    assembly: str,
    bins_dir: str,
    clean_bins_dir: Optional[str] = None,
) -> Recovery:
    """How much of *assembly* is in *bins_dir*, and in *clean_bins_dir* if given.

    Lengths come from the assembly rather than from the bins, so a contig a binner trimmed or
    renamed does not change the denominator - and a contig in a bin that is *not* in the assembly is
    ignored rather than counted, because it cannot be a share of something it is not part of.
    """
    lengths = contig_lengths(assembly)
    binned = binned_contigs(bins_dir) & set(lengths)
    clean: Set[str] = (binned_contigs(clean_bins_dir) & set(lengths)) if clean_bins_dir else set()
    return Recovery(
        assembly_contigs=len(lengths),
        assembly_bp=sum(lengths.values()),
        binned_contigs=len(binned),
        binned_bp=sum(lengths[name] for name in binned),
        clean_contigs=len(clean) if clean_bins_dir else None,
        clean_bp=sum(lengths[name] for name in clean) if clean_bins_dir else None,
    )


#: Columns of the recovery table.
RECOVERY_COLUMNS: Tuple[str, ...] = (
    "assembly_contigs",
    "assembly_bp",
    "binned_contigs",
    "binned_bp",
    "binned_fraction_bp",
    "clean_binned_contigs",
    "clean_binned_bp",
    "clean_binned_fraction_bp",
    "unbinned_bp",
    "caveat",
)


def recovery_row(recovery: Recovery) -> Dict[str, object]:
    return {
        "assembly_contigs": recovery.assembly_contigs,
        "assembly_bp": recovery.assembly_bp,
        "binned_contigs": recovery.binned_contigs,
        "binned_bp": recovery.binned_bp,
        "binned_fraction_bp": "%.4f" % recovery.binned_fraction,
        "clean_binned_contigs": "" if recovery.clean_contigs is None else recovery.clean_contigs,
        "clean_binned_bp": "" if recovery.clean_bp is None else recovery.clean_bp,
        "clean_binned_fraction_bp": (
            "" if recovery.clean_fraction is None else "%.4f" % recovery.clean_fraction
        ),
        "unbinned_bp": recovery.unbinned_bp,
        "caveat": UNBINNED_CAVEAT,
    }


# ── what fraction of the *community* the bins represent ──────────────────────────────────────

#: Printed with every appraisal. The two fractions answer different questions and are routinely
#: confused, and the per-cell one is almost always the smaller and the more sobering.
APPRAISAL_CAVEAT = (
    "this is the fraction of the community's *cells* the genomes represent, estimated from "
    "single-copy marker reads - not the fraction of assembled base pairs that were binned. A bin "
    "set can hold 60% of the assembly and represent 30% of the community, because what assembles "
    "well and what is abundant are different things"
)

#: SingleM's own summary rows, which are not samples and must not be averaged with them.
APPRAISAL_SUMMARY_ROWS = ("total", "average")


@dataclass(frozen=True)
class Appraisal:
    """How much of one sample's community a set of genomes accounts for.

    Coverages rather than counts: SingleM weights each marker sequence by its read coverage, so a
    genome representing an abundant organism counts for more than one representing a rare one.
    That is the point - the question is what share of the cells in the sample these genomes stand
    for, not what share of the lineages.
    """

    sample: str
    domain: str
    binned: float
    not_found: float
    #: Only present when an assembly was appraised alongside the genomes: sequence that assembled
    #: but reached no bin. None when no assembly was given, which is a different thing from zero.
    assembled: Optional[float] = None

    @property
    def total(self) -> float:
        return self.binned + self.not_found + (self.assembled or 0.0)

    @property
    def percent_binned(self) -> Optional[float]:
        """Share of the community the genomes represent, or None when nothing was detected."""
        return 100.0 * self.binned / self.total if self.total else None

    @property
    def percent_assembled_not_binned(self) -> Optional[float]:
        if self.assembled is None or not self.total:
            return None
        return 100.0 * self.assembled / self.total


def parse_appraisal(path: str) -> List[Appraisal]:
    """Parse ``singlem appraise``'s summary table.

    Two things about its shape. The columns vary - ``coverage_assembled`` appears only when an
    assembly was appraised - so columns are read by name and the absent one stays None rather than
    being read as zero. And the table carries ``total`` and ``average`` rows beside the per-sample
    ones; those are SingleM's summaries, not samples, and averaging them back in would double-count
    every sample.
    """
    rows: List[Appraisal] = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("sample") or "").strip()
            if not sample or sample in APPRAISAL_SUMMARY_ROWS:
                continue
            assembled = row.get("coverage_assembled")
            rows.append(
                Appraisal(
                    sample=sample,
                    domain=(row.get("domain") or "").strip(),
                    binned=_as_coverage(row.get("coverage_binned")),
                    not_found=_as_coverage(row.get("coverage_not_found")),
                    assembled=None if assembled is None else _as_coverage(assembled),
                )
            )
    return rows


def _as_coverage(value: Optional[str]) -> float:
    """SingleM writes `nan` in rows where a domain was not seen; that is zero coverage, not data."""
    try:
        number = float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0
    return 0.0 if math.isnan(number) else number


#: Columns of the appraisal table.
APPRAISAL_COLUMNS: Tuple[str, ...] = (
    "sample",
    "domain",
    "coverage_binned",
    "coverage_assembled_not_binned",
    "coverage_not_found",
    "percent_community_binned",
    "percent_community_assembled_not_binned",
    "caveat",
)


def appraisal_row(appraisal: Appraisal) -> Dict[str, object]:
    return {
        "sample": appraisal.sample,
        "domain": appraisal.domain,
        "coverage_binned": "%.4g" % appraisal.binned,
        "coverage_assembled_not_binned": (
            "" if appraisal.assembled is None else "%.4g" % appraisal.assembled
        ),
        "coverage_not_found": "%.4g" % appraisal.not_found,
        "percent_community_binned": (
            "" if appraisal.percent_binned is None else "%.2f" % appraisal.percent_binned
        ),
        "percent_community_assembled_not_binned": (
            ""
            if appraisal.percent_assembled_not_binned is None
            else "%.2f" % appraisal.percent_assembled_not_binned
        ),
        "caveat": APPRAISAL_CAVEAT,
    }


@dataclass(frozen=True)
class UnrepresentedLineage:
    """A lineage present in the reads with no genome representing it, and how abundant it is."""

    taxonomy: str
    coverage: float
    markers: int

    @property
    def rank_reached(self) -> str:
        """The finest rank SingleM could assign, which bounds what can be said about the lineage."""
        parts = [part.strip() for part in self.taxonomy.split(";") if part.strip()]
        return parts[-1] if parts else "unclassified"


def unrepresented_lineages(path: str, limit: Optional[int] = None) -> List[UnrepresentedLineage]:
    """Aggregate ``--output-unaccounted-for-otu-table`` by lineage, most abundant first.

    Aggregated rather than passed through: that table has one row per marker gene per sequence, so
    a single unrepresented organism appears up to as many times as there are markers, and a reader
    counting rows counts markers. Summing coverage per lineage and counting the markers behind it
    gives one row per lineage with the evidence beside it.
    """
    by_taxon: Dict[str, Tuple[float, int]] = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxonomy = (row.get("taxonomy") or "").strip() or "unclassified"
            coverage = _as_coverage(row.get("coverage"))
            total, markers = by_taxon.get(taxonomy, (0.0, 0))
            by_taxon[taxonomy] = (total + coverage, markers + 1)
    rows = [
        UnrepresentedLineage(taxonomy=taxonomy, coverage=coverage, markers=markers)
        for taxonomy, (coverage, markers) in by_taxon.items()
    ]
    rows.sort(key=lambda row: (-row.coverage, row.taxonomy))
    return rows[:limit] if limit else rows


#: Columns of the unrepresented-lineage table.
UNREPRESENTED_COLUMNS: Tuple[str, ...] = (
    "taxonomy",
    "rank_reached",
    "coverage",
    "percent_of_unrepresented",
    "marker_sequences",
)


def unrepresented_rows(
    lineages: Sequence[UnrepresentedLineage],
) -> List[Dict[str, object]]:
    total = sum(lineage.coverage for lineage in lineages)
    return [
        {
            "taxonomy": lineage.taxonomy,
            "rank_reached": lineage.rank_reached,
            "coverage": "%.4g" % lineage.coverage,
            "percent_of_unrepresented": (
                "%.2f" % (100.0 * lineage.coverage / total) if total else ""
            ),
            "marker_sequences": lineage.markers,
        }
        for lineage in lineages
    ]
