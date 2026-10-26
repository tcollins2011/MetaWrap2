"""Gene and pathway abundance per genome, extended to whole samples.

For gene family *g* in sample *s*, summing over genomes *k*:

    A(g,s) = sum_k  n(g,k) * abundance(k,s)

where ``n(g,k)`` is the copy number of *g* in genome *k* and ``abundance(k,s)`` is that genome's
abundance in that sample.

Coverage-based abundances (CoverM ``relative_abundance``, ``trimmed_mean``) are already proportional
to cell abundance and must not be divided by genome size. Only read counts (raw counts, salmon
``NumReads``) need that division. For MAGs, incompleteness cancels in
``aligned_bases / MAG_length``, so it biases the gene content term ``n(g,k)``, not the abundance
term.

Three units are reported:

=================================  ==========================================================
unit                               what it is comparable across
=================================  ==========================================================
``relative``                       composition within one sample; compositional, so it needs
                                   CoDA methods and must not go into a raw t-test
``per_genome_equivalent``          samples, which is what makes it the unit to report; needs
                                   SingleM's prokaryotic_fraction as the denominator
``length_normalised``              genes, since a short gene otherwise collects fewer reads
=================================  ==========================================================

The module also reports the fraction of a sample's functional signal attributable to MAGs. Across
92,187 metagenomes about 10% of assembled base pairs reach a medium- or high-quality MAG, and the
unbinned remainder is not a random sample, so a MAG-only profile covers part of the community.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


class Unit(Enum):
    """The denominator a value is expressed over. Written on every row."""

    RELATIVE = "relative"
    PER_GENOME_EQUIVALENT = "per_genome_equivalent"
    LENGTH_NORMALISED = "length_normalised"


#: What each unit may be compared across. Written into the output beside each value.
UNIT_COMPARABLE: Mapping[str, str] = {
    Unit.RELATIVE.value: "within one sample only; compositional, needs CoDA methods",
    Unit.PER_GENOME_EQUIVALENT.value: "across samples",
    Unit.LENGTH_NORMALISED.value: "across gene families within one sample",
}


class AbundanceKind(Enum):
    """Whether an abundance number is coverage-proportional already, or a read count.

    Determines whether genome size is divided out. Treating coverage as counts would divide by
    length twice and under-estimate large genomes.
    """

    #: CoverM relative_abundance / trimmed_mean, sylph, any depth-derived figure.
    COVERAGE = "coverage"
    #: Raw mapped-read counts, salmon NumReads.
    READ_COUNT = "read_count"


@dataclass(frozen=True)
class GenomeAbundance:
    """One genome's abundance in one sample, and what kind of number it is."""

    genome: str
    sample: str
    value: float
    kind: AbundanceKind
    #: Needed only to convert a read count into something cell-proportional. None for coverage.
    genome_length_bp: Optional[int] = None

    @property
    def cell_proportional(self) -> float:
        """The abundance as a per-cell quantity, whichever kind it started as.

        Coverage is returned as is. A read count is divided by genome length (0.0 if unknown).
        """
        if self.kind is AbundanceKind.COVERAGE:
            return self.value
        if not self.genome_length_bp:
            return 0.0
        return self.value / float(self.genome_length_bp)


def feature_abundance(
    per_genome: Mapping[str, Mapping[str, float]],
    abundances: Sequence[GenomeAbundance],
) -> Dict[str, Dict[str, float]]:
    """sample -> feature -> abundance, from per-genome feature counts and genome abundances.

    *per_genome* is genome -> feature -> copy number (a KO count, a module completeness, whatever
    the caller is propagating). A genome with no abundance entry for a sample contributes nothing to
    that sample.
    """
    out: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for entry in abundances:
        features = per_genome.get(entry.genome)
        if not features:
            continue
        scale = entry.cell_proportional
        if not scale:
            continue
        for feature, copies in features.items():
            out[entry.sample][feature] += copies * scale
    return {sample: dict(features) for sample, features in out.items()}


def to_relative(features: Mapping[str, float]) -> Dict[str, float]:
    """Shares of the sample's total, summing to 1. Compositional: see UNIT_COMPARABLE."""
    total = sum(features.values())
    if not total:
        return {feature: 0.0 for feature in features}
    return {feature: value / total for feature, value in features.items()}


def to_genome_equivalents(
    features: Mapping[str, float], prokaryotic_fraction: float
) -> Dict[str, float]:
    """Copies per genome equivalent: the cross-sample comparable unit.

    *prokaryotic_fraction* is the share of the sample's sequence that is bacterial or archaeal, from
    SingleM. Dividing by total reads instead would make profiles depend on host or eukaryotic DNA
    content. SingleM replaces MicrobeCensus, which has had no release since 2016.
    """
    if prokaryotic_fraction <= 0:
        raise ValueError(
            "prokaryotic_fraction must be positive; %r would make every value infinite"
            % (prokaryotic_fraction,)
        )
    return {feature: value / prokaryotic_fraction for feature, value in features.items()}


def to_length_normalised(
    features: Mapping[str, float], lengths: Mapping[str, float]
) -> Dict[str, float]:
    """Per kilobase of the feature, so two gene families can be compared with each other.

    Needed for between-gene comparisons: at equal cellular abundance a 3 kb gene collects ten times
    the reads of a 300 bp gene. A feature with no known length is passed through unchanged.
    """
    out: Dict[str, float] = {}
    for feature, value in features.items():
        length = lengths.get(feature)
        out[feature] = value / (length / 1000.0) if length else value
    return out


@dataclass(frozen=True)
class Assignable:
    """How much of a sample's functional signal a MAG-based profile accounts for."""

    sample: str
    #: Total abundance of all genomes observed in the sample, binned or not.
    total_abundance: float
    #: Abundance of the genomes that have functional annotation.
    annotated_abundance: float
    #: Share of the community SingleM says the bins represent, when `singlem appraise` has run.
    binned_fraction_of_community: Optional[float] = None

    @property
    def assignable_fraction(self) -> float:
        if not self.total_abundance:
            return 0.0
        return self.annotated_abundance / self.total_abundance


def assignable(
    abundances: Sequence[GenomeAbundance],
    annotated_genomes: Iterable[str],
    binned_fraction: Optional[Mapping[str, float]] = None,
) -> List[Assignable]:
    """Per-sample fraction of abundance covered by annotated genomes.

    ``assignable_fraction`` is relative to the genomes in the abundance table (were all quantified
    MAGs annotated). ``binned_fraction_of_community``, from SingleM, is relative to the whole
    community (what fraction of cells a MAG represents), and is usually much smaller.
    """
    annotated = set(annotated_genomes)
    totals: Dict[str, float] = defaultdict(float)
    covered: Dict[str, float] = defaultdict(float)
    for entry in abundances:
        value = entry.cell_proportional
        totals[entry.sample] += value
        if entry.genome in annotated:
            covered[entry.sample] += value
    return [
        Assignable(
            sample=sample,
            total_abundance=total,
            annotated_abundance=covered.get(sample, 0.0),
            binned_fraction_of_community=(binned_fraction or {}).get(sample),
        )
        for sample, total in sorted(totals.items())
    ]


#: Columns of the feature abundance table.
ABUNDANCE_COLUMNS: Tuple[str, ...] = (
    "sample",
    "feature",
    "feature_type",
    "value",
    "unit",
    "comparable_across",
    "abundance_kind",
    "n_genomes_contributing",
)

#: Columns of the assignable-fraction table.
ASSIGNABLE_COLUMNS: Tuple[str, ...] = (
    "sample",
    "total_abundance",
    "annotated_abundance",
    "assignable_fraction",
    "binned_fraction_of_community",
    "caveat",
)

#: Printed in every assignable row.
ASSIGNABLE_CAVEAT = (
    "a MAG-based functional profile describes only the binned fraction of the sample; "
    "across 92,187 metagenomes ~10% of assembled bp reach a medium- or high-quality MAG, "
    "and the unbinned remainder is a biased subset"
)


def abundance_rows(
    sample: str,
    features: Mapping[str, float],
    unit: Unit,
    feature_type: str,
    kind: AbundanceKind,
    contributors: Optional[Mapping[str, int]] = None,
) -> List[Dict[str, object]]:
    return [
        {
            "sample": sample,
            "feature": feature,
            "feature_type": feature_type,
            "value": "%.6g" % value,
            "unit": unit.value,
            "comparable_across": UNIT_COMPARABLE[unit.value],
            "abundance_kind": kind.value,
            "n_genomes_contributing": (contributors or {}).get(feature, 0),
        }
        for feature, value in sorted(features.items(), key=lambda kv: (-kv[1], kv[0]))
    ]


def assignable_row(entry: Assignable) -> Dict[str, object]:
    return {
        "sample": entry.sample,
        "total_abundance": "%.6g" % entry.total_abundance,
        "annotated_abundance": "%.6g" % entry.annotated_abundance,
        "assignable_fraction": "%.4f" % entry.assignable_fraction,
        "binned_fraction_of_community": (
            ""
            if entry.binned_fraction_of_community is None
            else "%.4f" % entry.binned_fraction_of_community
        ),
        "caveat": ASSIGNABLE_CAVEAT,
    }


def contributors_per_feature(
    per_genome: Mapping[str, Mapping[str, float]],
    abundances: Sequence[GenomeAbundance],
    sample: str,
) -> Dict[str, int]:
    """How many genomes contribute to each feature in one sample.

    Distinguishes a feature carried by one genome from one carried by many at the same abundance.
    """
    counts: Dict[str, int] = defaultdict(int)
    present = {
        entry.genome for entry in abundances if entry.sample == sample and entry.cell_proportional
    }
    for genome in present:
        for feature in per_genome.get(genome, {}):
            counts[feature] += 1
    return dict(counts)


#: Column names `singlem microbial_fraction` uses for the figure, current and older. Checked
#: against singlem 0.21.4's read_fraction.py.
SINGLEM_FRACTION_COLUMNS: Tuple[str, ...] = (
    "read_fraction",
    "prokaryotic_fraction",
    "bacterial_archaeal_fraction",
)


def parse_singlem_fraction(path: str) -> Dict[str, float]:
    """sample -> prokaryotic fraction (0-1), from ``singlem microbial_fraction``.

    SingleM writes a percentage (``read_fraction.py`` computes ``account / metagenome_size * 100``,
    capped at 100), so the value is always divided by 100. Inferring the unit from the magnitude
    would misread a host-dominated sample at 0.8% as 0.8 (80%).
    """
    out: Dict[str, float] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        column = next((name for name in SINGLEM_FRACTION_COLUMNS if name in fields), None)
        if column is None:
            raise KeyError(
                "%s has none of the expected fraction columns %r; found %r"
                % (path, SINGLEM_FRACTION_COLUMNS, fields)
            )
        sample_column = "sample" if "sample" in fields else fields[0]
        for row in reader:
            raw = (row.get(column) or "").strip().rstrip("%")
            try:
                value = float(raw)
            except ValueError:
                continue
            if value > 100.0:
                raise ValueError(
                    "%s reports %s=%g for sample %r, above 100%%. The column is not a "
                    "percentage; check that this is the right file."
                    % (path, column, value, row[sample_column])
                )
            out[row[sample_column]] = value / 100.0
    return out
