"""Gene and pathway abundance: per genome, and extended to the whole sample.

The arithmetic is one line. For gene family *g* in sample *s*, summing over genomes *k*:

    A(g,s) = sum_k  n(g,k) * abundance(k,s)

where ``n(g,k)`` is the copy number of *g* in genome *k* and ``abundance(k,s)`` is that genome's
abundance in that sample. Everything difficult about this file is in the units and the caveats.

**Do not divide a coverage-based abundance by genome size.** CoverM's ``relative_abundance`` and
``trimmed_mean`` are already coverage-proportional, which is to say already proportional to cell
abundance, so dividing again double-corrects. The division is needed only for read-count numbers
(raw counts, salmon's ``NumReads``). For MAGs specifically, coverage estimates cells correctly
because incompleteness cancels in ``aligned_bases / MAG_length`` - what incompleteness biases is the
**gene content** term ``n(g,k)``, not the abundance term. That asymmetry is the single easiest thing
here to get backwards, which is why it is written down rather than implied.

**Three units, because there are three honest denominators** and picking one silently misleads:

=================================  ==========================================================
unit                               what it is comparable across
=================================  ==========================================================
``relative``                       composition within one sample; compositional, so it needs
                                   CoDA methods and must not go into a raw t-test
``per_genome_equivalent``          samples, which is what makes it the unit to report; needs
                                   SingleM's prokaryotic_fraction as the denominator
``length_normalised``              genes, since a short gene otherwise collects fewer reads
=================================  ==========================================================

**And one metric that is not an abundance.** The fraction of a sample's functional signal that can
be attributed to a MAG at all. Across 92,187 metagenomes only about 10% of assembled base pairs
reach a medium- or high-quality MAG, and the unbinned remainder is biased rather than random, so a
MAG-only functional profile describes a narrow slice of the sample. Reporting the slice's size
beside it is the difference between a profile and a claim about the community.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


class Unit(Enum):
    """Which denominator a number is over. On every row, because it changes what it means."""

    RELATIVE = "relative"
    PER_GENOME_EQUIVALENT = "per_genome_equivalent"
    LENGTH_NORMALISED = "length_normalised"


#: What each unit may be compared across, carried into the output so the constraint travels with
#: the data rather than living in documentation.
UNIT_COMPARABLE: Mapping[str, str] = {
    Unit.RELATIVE.value: "within one sample only; compositional, needs CoDA methods",
    Unit.PER_GENOME_EQUIVALENT.value: "across samples",
    Unit.LENGTH_NORMALISED.value: "across gene families within one sample",
}


class AbundanceKind(Enum):
    """Whether an abundance number is coverage-proportional already, or a read count.

    The distinction decides whether genome size has to be divided out, and getting it wrong is a
    systematic error in one direction: treating coverage as counts divides by length twice and makes
    large genomes look rare.
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

        Coverage is already there. A read count is divided by genome length once - and only once,
        which is the whole reason this is a method rather than a line at each call site.
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
    it rather than defaulting to one: a MAG absent from a sample's abundance table was not observed
    there, and treating it as present at unit abundance imports another sample's organism.
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
    SingleM. It is the denominator because the alternative - dividing by total reads - makes a
    sample's functional profile depend on how much host or eukaryotic DNA was in it, so two gut
    samples with different host content stop being comparable for a reason that has nothing to do
    with their microbiomes.

    MicrobeCensus did this job and is abandoned (last release 2016), which is why SingleM.
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

    Mandatory for any between-gene comparison: a 300 bp gene and a 3 kb gene at the same cellular
    abundance collect reads in a 1:10 ratio, and without this the longer one looks ten times more
    abundant. A feature with no known length is passed through unchanged rather than dropped, and
    the row says so, because silently omitting it would make the table's total wrong.
    """
    out: Dict[str, float] = {}
    for feature, value in features.items():
        length = lengths.get(feature)
        out[feature] = value / (length / 1000.0) if length else value
    return out


@dataclass(frozen=True)
class Assignable:
    """How much of a sample's functional signal a MAG-based profile actually accounts for."""

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
    """The honesty metric, per sample.

    Two different denominators and both are worth having. ``assignable_fraction`` is of the genomes
    in the abundance table, which answers "did every MAG I quantified get annotated". The SingleM
    figure is of the *community*, which answers the harder question - what fraction of the cells
    present are represented by a MAG at all - and is the one that is usually small.
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

#: Printed in every assignable row, because the number is meaningless without it.
ASSIGNABLE_CAVEAT = (
    "a MAG-based functional profile describes only the binned fraction of the sample; "
    "across 92,187 metagenomes ~10% of assembled bp reach a medium- or high-quality MAG, "
    "and the remainder is biased rather than random"
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

    A feature carried by one rare genome and one carried by forty are different findings at the same
    abundance, and nothing else in the row distinguishes them.
    """
    counts: Dict[str, int] = defaultdict(int)
    present = {
        entry.genome for entry in abundances if entry.sample == sample and entry.cell_proportional
    }
    for genome in present:
        for feature in per_genome.get(genome, {}):
            counts[feature] += 1
    return dict(counts)


#: The column `singlem microbial_fraction` writes the figure in, and the ones it has used before.
#: Checked against singlem 0.21.4's read_fraction.py rather than assumed.
SINGLEM_FRACTION_COLUMNS: Tuple[str, ...] = (
    "read_fraction",
    "prokaryotic_fraction",
    "bacterial_archaeal_fraction",
)


def parse_singlem_fraction(path: str) -> Dict[str, float]:
    """sample -> prokaryotic fraction (0-1), from ``singlem microbial_fraction``.

    **SingleM writes a percentage**, not a fraction: ``read_fraction.py`` computes
    ``account / metagenome_size * 100`` and caps it at 100. So this divides by 100 unconditionally.

    The first version guessed instead - "a value above 1 must be a percentage, since a fraction
    above 1 is not meaningful" - which is right for a typical sample at 85.3 and catastrophically
    wrong for a host-dominated one at 0.8, where it would read 0.8% as 80% and scale every gene
    abundance in that sample down by 100x. A unit that can be looked up should never be inferred,
    and this is the same class of error as the skani alignment-fraction confusion that merged two
    *Enterococcus* species in dRep.
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
                    "%s reports %s=%g for sample %r, which is above 100%%. That is not a "
                    "percentage, so this file is not what it claims to be."
                    % (path, column, value, row[sample_column])
                )
            out[row[sample_column]] = value / 100.0
    return out
