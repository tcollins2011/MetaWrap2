"""Binned fraction of an assembly, community share of the bins, and the caveats printed with them.

This module computes the share of assembled sequence that is in a bin and attaches two caveats:

* Across 92,187 metagenomes in SPIRE, only about 10% of assembled base pairs reach a medium- or
  high-quality MAG, and the unbinned remainder is not a random sample: it holds an estimated 145
  predicted undescribed phyla.
* Treichel et al. 2026 found that 54.5-81.8% of high-quality MAGs represent a single strain, with
  a mean of 2.3 reference genomes per MAG. CheckM2 does not measure this; GUNC detects part of it
  at the taxonomic level and runs by default.

The binned fraction is measured in base pairs. A bin set can hold 2% of the contigs and 60% of the
sequence, because binners discard short contigs, which are numerous.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .io.seqio import contig_id, iter_fasta

#: Printed with every binned-fraction figure.
UNBINNED_CAVEAT = (
    "the unbinned fraction is not a random sample of the assembly: binners discard contigs they "
    "cannot place, mostly short, low-coverage and novel ones. Across 92,187 metagenomes "
    "(SPIRE) ~10% of assembled bp are in a medium- or high-quality MAG"
)

#: Printed with every quality-tier count.
CHIMERISM_CAVEAT = (
    "a 70/10 or 90/5 score means a bin has the markers of one genome, which does not show it is "
    "one genome. Of MAGs called high-quality, 54.5-81.8% represent a single strain (mean 2.3 "
    "reference genomes per MAG; Treichel et al. 2026). GUNC, which runs by default, detects the "
    "taxonomically incoherent ones"
)


@dataclass(frozen=True)
class Recovery:
    """What share of an assembly a bin set accounts for."""

    assembly_contigs: int
    assembly_bp: int
    binned_contigs: int
    binned_bp: int
    #: Contigs and bp in bins that pass the quality thresholds. None when no clean set was given,
    #: so that "not measured" is not reported as 0.0%.
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

    Returned as a set because a contig can appear in more than one bin, for example in overlapping
    bin sets from different binners.
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

    Lengths come from the assembly, so trimmed or renamed contigs in a bin do not change the totals.
    Bin contigs that are not in the assembly are ignored.
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

#: Printed with every appraisal. The community fraction is usually smaller than the binned
#: fraction of the assembly, and the two are often confused.
APPRAISAL_CAVEAT = (
    "this is the fraction of the community's cells that the genomes represent, estimated from "
    "single-copy marker reads. It is not the fraction of assembled base pairs that were binned. "
    "A bin set can hold 60% of the assembly and represent 30% of the community, because the "
    "genomes that assemble well are not always the abundant ones"
)

#: SingleM's summary rows, which are skipped when reading samples.
APPRAISAL_SUMMARY_ROWS = ("total", "average")


@dataclass(frozen=True)
class Appraisal:
    """How much of one sample's community a set of genomes accounts for.

    Values are SingleM coverages, weighted by read coverage, so they estimate the share of cells,
    not of lineages.
    """

    sample: str
    domain: str
    binned: float
    not_found: float
    #: Coverage of sequence that assembled but reached no bin. None when no assembly was
    #: appraised.
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

    Columns are read by name: ``coverage_assembled`` is present only when an assembly was
    appraised, and stays None otherwise. SingleM's ``total`` and ``average`` rows are skipped.
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
    """Parse a coverage value; empty, invalid and `nan` (domain not seen) give 0.0."""
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
        """The finest rank SingleM could assign."""
        parts = [part.strip() for part in self.taxonomy.split(";") if part.strip()]
        return parts[-1] if parts else "unclassified"


def unrepresented_lineages(path: str, limit: Optional[int] = None) -> List[UnrepresentedLineage]:
    """Aggregate ``--output-unaccounted-for-otu-table`` by lineage, most abundant first.

    The table has one row per marker gene per sequence, so coverage is summed per lineage and the
    number of marker sequences is kept beside it.
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
