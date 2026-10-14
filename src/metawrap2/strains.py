"""Within-population variation: is this the same strain in two samples, and how varied is it?

A MAG is a consensus of a population, not a genome. Two samples can both carry "the same species"
at 99.9% ANI between their MAGs and be carrying different strains, and the pipeline has had no way
to say so - which matters for every transmission, persistence or engraftment claim anyone makes from
a cohort.

inStrain answers it from the alignments rather than from the consensus: it reads every position's
allele frequencies, so a population with two variants at 30/70 is visibly that rather than a
consensus of whichever won. Two numbers out of it carry the weight:

* **popANI** between two samples for one genome. Alleles are compared as *populations* - two samples
  share a position when their allele sets overlap at all - so a minor variant present in both does
  not count as a difference. inStrain's own threshold for "same strain" is popANI >= 99.999%, which
  is roughly one difference per 100 kb.
* **nucleotide diversity** within one sample for one genome, which says how varied the population
  is at all. A clonal population and a mixture of two strains can have the same consensus genome.

The thing that makes popANI misleading, and the reason this module reports it alongside, is
**breadth**: popANI is computed only over positions both samples covered deeply enough. Two samples
sharing 2% of a genome at 100% identity produce a popANI of 100%, and that is not evidence of
anything. inStrain's guidance is a floor of 50% of the genome compared, and
:data:`MIN_PERCENT_COMPARED` applies it rather than leaving it to the reader.

Also worth knowing before quoting a distance: Shaw & Yu's evaluation of strain-distance tools found
most of them overestimate SNP distances, and TRACS the least. inStrain is here because it is
maintained, installs cleanly and reports the breadth its own number depends on; a distance from it
is a comparison between samples processed identically, not an absolute SNP count.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

#: popANI at or above which inStrain considers two samples to carry the same strain. Its own
#: default, stated here so a verdict in this pipeline means the same thing as one in its paper.
SAME_STRAIN_POPANI = 99.999

#: Fraction of a genome that must be compared for a popANI to mean anything, as a percentage. Below
#: this the comparison is between slivers of two genomes: two samples sharing 2% of a genome at
#: 100% identity give popANI 100%, which says nothing about the strain. inStrain's own guidance.
MIN_PERCENT_COMPARED = 50.0

#: Minimum breadth (fraction of the genome covered at the profiler's depth floor) for a within-sample
#: diversity figure to be worth reporting. A genome covered over 10% of its length has a nucleotide
#: diversity measured on 10% of its length.
MIN_BREADTH = 0.5

#: Coverage below which inStrain's per-position allele frequencies are too thin to call a minor
#: variant at all - its own `min_cov` default.
MIN_COVERAGE = 5.0


class StrainVerdict(Enum):
    """What a comparison between two samples supports for one genome."""

    SAME = "same_strain"
    DIFFERENT = "different_strain"
    #: Compared, but over too little of the genome for either answer.
    UNDETERMINED = "undetermined"
    #: Not compared: the genome is not present deeply enough in both samples.
    ABSENT = "not_in_both"


#: What each verdict licenses saying, emitted in the table so the constraint travels with the data.
VERDICT_LANGUAGE: Dict[str, str] = {
    StrainVerdict.SAME.value: (
        "the populations in these two samples are indistinguishable at inStrain's threshold "
        "(popANI >= %g%% over at least %g%% of the genome). That is consistent with transmission, "
        "persistence or a shared source; it does not establish which"
        % (
            SAME_STRAIN_POPANI,
            MIN_PERCENT_COMPARED,
        )
    ),
    StrainVerdict.DIFFERENT.value: (
        "the populations differ by more than inStrain's threshold over enough of the genome to "
        "tell. Different strains of the same species"
    ),
    StrainVerdict.UNDETERMINED.value: (
        "too little of the genome was compared for the identity to mean anything. Not evidence of "
        "either answer - the usual cause is low coverage in one of the two samples"
    ),
    StrainVerdict.ABSENT.value: (
        "the genome was not detected deeply enough in both samples, so there was nothing to "
        "compare. Absence here is about this sequencing depth, not about the community"
    ),
}


@dataclass(frozen=True)
class GenomeProfile:
    """One genome's population in one sample, as inStrain's `genome_info.tsv` reports it."""

    sample: str
    genome: str
    #: Mean coverage over the genome's covered positions.
    coverage: Optional[float]
    #: Fraction of the genome covered at all, 0-1.
    breadth: Optional[float]
    #: Fraction covered at or above the profiler's depth floor, 0-1 - the positions SNVs are
    #: actually called over. This is the one a diversity figure has to be gated on: a genome can be
    #: 91% covered and have a quarter of that deep enough to call a variant, and the diversity is a
    #: statement about the second number.
    breadth_min_cov: Optional[float]
    #: Breadth inStrain expects given the coverage, 0-1. The ratio of the two says whether the
    #: genome is really present or whether reads from a relative are mapping to part of it.
    expected_breadth: Optional[float]
    #: Mean per-position nucleotide diversity over covered positions.
    nucleotide_diversity: Optional[float]
    #: Positions where a second allele is present above inStrain's thresholds.
    snvs: Optional[int]
    #: Positions whose major allele differs from the reference.
    consensus_snvs: Optional[int]
    #: Divergent sites per kb of compared sequence, as inStrain reports it.
    divergent_site_density: Optional[float]

    @property
    def breadth_ratio(self) -> Optional[float]:
        """Observed breadth over expected breadth, or None if either is unknown.

        Well below 1 means the genome is covered unevenly - reads piling on conserved regions of a
        relative rather than spread over a genome that is present. inStrain's own diagnostic, and
        the reason breadth alone does not establish presence.
        """
        if not self.breadth or not self.expected_breadth:
            return None
        return self.breadth / self.expected_breadth

    @property
    def usable_breadth(self) -> Optional[float]:
        """The breadth a diversity figure is measured over: at the depth floor, not at any depth."""
        return self.breadth_min_cov if self.breadth_min_cov is not None else self.breadth

    @property
    def diversity_is_reportable(self) -> bool:
        """Whether this genome is covered well enough for its diversity to describe the genome."""
        breadth = self.usable_breadth
        return (
            breadth is not None and breadth >= MIN_BREADTH and (self.coverage or 0) >= MIN_COVERAGE
        )


@dataclass(frozen=True)
class GenomeComparison:
    """One genome compared between two samples, as `genomeWide_compare.tsv` reports it."""

    genome: str
    sample_a: str
    sample_b: str
    #: Population ANI: alleles compared as sets, so a minor variant in both is not a difference.
    pop_ani: Optional[float]
    #: Consensus ANI: major alleles only, which is the stricter and less meaningful of the two for
    #: a population, and is reported because the two diverging is informative.
    con_ani: Optional[float]
    #: Percentage of the genome compared - the denominator popANI is computed over.
    percent_compared: Optional[float]
    population_snvs: Optional[int]
    consensus_snvs: Optional[int]

    @property
    def verdict(self) -> StrainVerdict:
        if self.pop_ani is None or self.percent_compared is None:
            return StrainVerdict.ABSENT
        if self.percent_compared < MIN_PERCENT_COMPARED:
            return StrainVerdict.UNDETERMINED
        return StrainVerdict.SAME if self.pop_ani >= SAME_STRAIN_POPANI else StrainVerdict.DIFFERENT

    @property
    def population_snv_distance(self) -> Optional[float]:
        """Population SNVs per megabase of compared sequence, or None without the breadth.

        A rate rather than a count, because the count is a function of how much was compared: 10
        differences over 2% of a genome and 10 over all of it are not the same finding. Still a
        comparison between samples processed identically rather than an absolute SNP distance - see
        the module docstring on what strain-distance tools overestimate.
        """
        if self.population_snvs is None or not self.percent_compared:
            return None
        return self.population_snvs / max(self.percent_compared / 100.0, 1e-9)


def _as_float(value: object) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return None if math.isnan(number) else number


def _as_percent(value: object) -> Optional[float]:
    """A fraction from inStrain's compare table, as a percentage.

    Unconditional, not inferred from the magnitude: popANI is ~0.99999 and percent_compared can be
    anything from 0.001 to 1.0, so "greater than 1 must already be a percentage" would be a rule
    that works on one column and silently fails on the other. The same mistake, and the same fix,
    as `parse_singlem_fraction`.
    """
    number = _as_float(value)
    return None if number is None else number * 100.0


def _as_int(value: object) -> Optional[int]:
    number = _as_float(value)
    return None if number is None else int(number)


def _strip_extension(name: str) -> str:
    """A genome name without its FASTA extension.

    ``os.path.splitext`` is wrong here for the same reason it is wrong in the clean-MAG join:
    metaWRAP names bins ``bin.13``, so splitext turns every bin in a set into ``bin``.
    """
    for extension in (".fasta", ".fa", ".fna", ".fas"):
        if name.endswith(extension):
            return name[: -len(extension)]
    return name


#: inStrain's profile columns. Named so a schema change is a KeyError rather than empty output.
PROFILE_COLUMNS = {
    "genome": "genome",
    "coverage": "coverage",
    "breadth": "breadth",
    "breadth_min_cov": "breadth_minCov",
    "expected_breadth": "breadth_expected",
    "diversity": "nucl_diversity",
    "snvs": "SNV_count",
    "consensus_snvs": "SNS_count",
    "density": "divergent_site_count",
}

#: inStrain's compare table stores popANI, conANI and the compared fraction as **fractions**, not
#: percentages, while its documentation and its paper quote popANI as "99.999%". Read literally
#: against the published threshold, every comparison is a different strain and every breadth is
#: below the floor - which is what the first real run reported. Parsed values are multiplied by 100
#: here so everything downstream, including the thresholds above, is in one unit.
FRACTION_COLUMNS = ("pop_ani", "con_ani", "percent_compared")

#: inStrain's compare columns.
COMPARE_COLUMNS = {
    "genome": "genome",
    "a": "name1",
    "b": "name2",
    "pop_ani": "popANI",
    "con_ani": "conANI",
    "percent_compared": "percent_compared",
    "population_snvs": "population_SNPs",
    "consensus_snvs": "consensus_SNPs",
}


def parse_profile(path: str, sample: str) -> List[GenomeProfile]:
    """Read one sample's `genome_info.tsv` from `inStrain profile`."""
    out: List[GenomeProfile] = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and PROFILE_COLUMNS["genome"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like inStrain's genome_info.tsv: no %r column"
                % (path, PROFILE_COLUMNS["genome"])
            )
        for row in reader:
            genome = _strip_extension((row.get(PROFILE_COLUMNS["genome"]) or "").strip())
            if not genome:
                continue
            out.append(
                GenomeProfile(
                    sample=sample,
                    genome=genome,
                    coverage=_as_float(row.get(PROFILE_COLUMNS["coverage"])),
                    breadth=_as_float(row.get(PROFILE_COLUMNS["breadth"])),
                    breadth_min_cov=_as_float(row.get(PROFILE_COLUMNS["breadth_min_cov"])),
                    expected_breadth=_as_float(row.get(PROFILE_COLUMNS["expected_breadth"])),
                    nucleotide_diversity=_as_float(row.get(PROFILE_COLUMNS["diversity"])),
                    snvs=_as_int(row.get(PROFILE_COLUMNS["snvs"])),
                    consensus_snvs=_as_int(row.get(PROFILE_COLUMNS["consensus_snvs"])),
                    divergent_site_density=_as_float(row.get(PROFILE_COLUMNS["density"])),
                )
            )
    return out


def parse_comparison(path: str) -> List[GenomeComparison]:
    """Read `genomeWide_compare.tsv` from `inStrain compare`.

    Sample names come back as inStrain writes them - the profile directory names - and the caller
    relabels them, exactly as the appraisal path does for SingleM: a table whose sample names do not
    match the rest of the pipeline's joins to nothing.
    """
    out: List[GenomeComparison] = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and COMPARE_COLUMNS["pop_ani"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like inStrain's genomeWide_compare.tsv: no %r column"
                % (path, COMPARE_COLUMNS["pop_ani"])
            )
        for row in reader:
            genome = _strip_extension((row.get(COMPARE_COLUMNS["genome"]) or "").strip())
            if not genome:
                continue
            out.append(
                GenomeComparison(
                    genome=genome,
                    sample_a=(row.get(COMPARE_COLUMNS["a"]) or "").strip(),
                    sample_b=(row.get(COMPARE_COLUMNS["b"]) or "").strip(),
                    pop_ani=_as_percent(row.get(COMPARE_COLUMNS["pop_ani"])),
                    con_ani=_as_percent(row.get(COMPARE_COLUMNS["con_ani"])),
                    percent_compared=_as_percent(row.get(COMPARE_COLUMNS["percent_compared"])),
                    population_snvs=_as_int(row.get(COMPARE_COLUMNS["population_snvs"])),
                    consensus_snvs=_as_int(row.get(COMPARE_COLUMNS["consensus_snvs"])),
                )
            )
    return out


def _label_key(name: str) -> str:
    """inStrain's own name for a sample, reduced to something a caller can map.

    `compare` labels a row with the BAM's basename (``ERR011347.bam``) while `profile` is given a
    directory, so the two halves of this module see different names for one sample. Stripping the
    extension is what lets the caller's own sample names be put back.
    """
    base = os.path.basename(name.strip())
    for extension in (".bam", ".sam", ".cram"):
        if base.endswith(extension):
            return base[: -len(extension)]
    return base


def relabel(
    comparisons: Sequence[GenomeComparison], labels: Mapping[str, str]
) -> List[GenomeComparison]:
    """Replace inStrain's sample names with the pipeline's, leaving unknown ones alone."""
    return [
        GenomeComparison(
            genome=c.genome,
            sample_a=labels.get(_label_key(c.sample_a), _label_key(c.sample_a)),
            sample_b=labels.get(_label_key(c.sample_b), _label_key(c.sample_b)),
            pop_ani=c.pop_ani,
            con_ani=c.con_ani,
            percent_compared=c.percent_compared,
            population_snvs=c.population_snvs,
            consensus_snvs=c.consensus_snvs,
        )
        for c in comparisons
    ]


#: Columns of the per-sample population table.
PROFILE_TABLE_COLUMNS: Tuple[str, ...] = (
    "sample",
    "genome",
    "coverage",
    "breadth",
    "breadth_at_min_coverage",
    "expected_breadth",
    "breadth_ratio",
    "nucleotide_diversity",
    "snvs",
    "consensus_snvs",
    "diversity_reportable",
)


def profile_row(profile: GenomeProfile) -> Dict[str, str]:
    def num(value: Optional[float], places: int = 4) -> str:
        return "" if value is None else "%.*f" % (places, value)

    return {
        "sample": profile.sample,
        "genome": profile.genome,
        "coverage": num(profile.coverage, 2),
        "breadth": num(profile.breadth),
        "breadth_at_min_coverage": num(profile.breadth_min_cov),
        "expected_breadth": num(profile.expected_breadth),
        "breadth_ratio": num(profile.breadth_ratio),
        # Blank rather than 0 when the genome is too thin to measure: a diversity of zero is a
        # clonal population, which is a finding, and this is not it.
        "nucleotide_diversity": (
            num(profile.nucleotide_diversity, 6) if profile.diversity_is_reportable else ""
        ),
        "snvs": "" if profile.snvs is None else str(profile.snvs),
        "consensus_snvs": "" if profile.consensus_snvs is None else str(profile.consensus_snvs),
        "diversity_reportable": "yes" if profile.diversity_is_reportable else "no",
    }


#: Columns of the sample-pair comparison table.
COMPARE_TABLE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "sample_a",
    "sample_b",
    "pop_ani",
    "con_ani",
    "percent_genome_compared",
    "population_snvs",
    "consensus_snvs",
    "population_snvs_per_genome_fraction",
    "verdict",
    "permitted_language",
    "popani_threshold",
    "min_percent_compared",
)


def comparison_row(comparison: GenomeComparison) -> Dict[str, str]:
    verdict = comparison.verdict

    def num(value: Optional[float], places: int = 4) -> str:
        return "" if value is None else "%.*f" % (places, value)

    return {
        "genome": comparison.genome,
        "sample_a": comparison.sample_a,
        "sample_b": comparison.sample_b,
        # Six decimal places on purpose: the threshold is 99.999%, so two would make every
        # comparison read as 100.00 and the verdict column look arbitrary.
        "pop_ani": num(comparison.pop_ani, 6),
        "con_ani": num(comparison.con_ani, 6),
        "percent_genome_compared": num(comparison.percent_compared, 2),
        "population_snvs": (
            "" if comparison.population_snvs is None else str(comparison.population_snvs)
        ),
        "consensus_snvs": (
            "" if comparison.consensus_snvs is None else str(comparison.consensus_snvs)
        ),
        "population_snvs_per_genome_fraction": num(comparison.population_snv_distance, 1),
        "verdict": verdict.value,
        "permitted_language": VERDICT_LANGUAGE[verdict.value],
        "popani_threshold": "%g" % SAME_STRAIN_POPANI,
        "min_percent_compared": "%g" % MIN_PERCENT_COMPARED,
    }


def summarise(comparisons: Sequence[GenomeComparison]) -> Dict[str, int]:
    """Counts per verdict, for the end-of-run report."""
    out: Dict[str, int] = {verdict.value: 0 for verdict in StrainVerdict}
    for comparison in comparisons:
        out[comparison.verdict.value] += 1
    return out


def shared_strains(
    comparisons: Iterable[GenomeComparison],
) -> List[GenomeComparison]:
    """The comparisons that support a same-strain claim, breadth included."""
    return [c for c in comparisons if c.verdict is StrainVerdict.SAME]


#: Printed with every popANI figure.
BREADTH_CAVEAT = (
    "popANI is computed only over positions both samples covered deeply enough, so it has to be "
    "read with `percent_genome_compared`: two samples sharing 2%% of a genome at 100%% identity "
    "produce popANI 100%%, which is not evidence of a shared strain. Comparisons below %g%% of the "
    "genome are reported as undetermined rather than as either answer" % MIN_PERCENT_COMPARED
)

#: Printed with every nucleotide-diversity figure.
DIVERSITY_CAVEAT = (
    "nucleotide diversity is measured over the positions covered at or above %gx, so a genome "
    "present at low coverage has its diversity measured on a fraction of itself. Genomes below "
    "%g breadth are reported without a diversity rather than with a thin one"
    % (
        MIN_COVERAGE,
        MIN_BREADTH,
    )
)


def comparable_genomes(
    profiles: Sequence[GenomeProfile], samples: Sequence[str]
) -> Dict[Tuple[str, str], int]:
    """How many genomes each *pair* of samples could even be compared for.

    Per pair, not across every sample at once. inStrain compares pairs, and the intersection over
    all samples is a different and much smaller number: on three real gut samples it was 0 while
    42 pairwise comparisons existed, so reporting it beside them read as a contradiction.

    "Comparable" means covered deeply enough in both to be worth a popANI - which is the first
    thing that explains a table of undetermined verdicts.
    """
    by_sample: Dict[str, set] = {}
    for profile in profiles:
        if profile.diversity_is_reportable:
            by_sample.setdefault(profile.sample, set()).add(profile.genome)
    out: Dict[Tuple[str, str], int] = {}
    ordered = list(samples)
    for index, first in enumerate(ordered):
        for second in ordered[index + 1 :]:
            out[(first, second)] = len(by_sample.get(first, set()) & by_sample.get(second, set()))
    return out


# ── TRACS, the second backend ────────────────────────────────────────────────────────────────

#: TRACS' own default SNP threshold for calling a transmission pair. Its paper's argument is that
#: most tools *overestimate* SNP distances, so a threshold calibrated on another tool's distances
#: does not transfer - which is the reason to keep its number rather than reuse inStrain's.
TRACS_SNP_THRESHOLD = 10

#: Fraction of a reference genome TRACS requires to have read coverage before it will use it. Its
#: own floor, and the same idea as MIN_PERCENT_COMPARED here: a distance computed over a sliver of
#: a genome is not a distance. Recorded because it is why a shallow sample yields no rows at all.
TRACS_MIN_BREADTH = 0.25

#: TRACS' transmission-distance columns, as `transmission_distances.csv` carries them.
TRACS_COLUMNS = {
    "a": "sampleA",
    "b": "sampleB",
    "snp_distance": "SNP distance",
    "transmission_distance": "transmission distance",
    "expected_k": "expected K",
    "filtered_snp_distance": "filtered SNP distance",
    "sites": "sites considered",
    "msa": "MSA file",
}


@dataclass(frozen=True)
class TracsDistance:
    """One pair of samples compared for one reference genome, as TRACS reports it."""

    genome: str
    sample_a: str
    sample_b: str
    #: Raw SNP distance between the two samples' consensus sequences for this genome.
    snp_distance: Optional[float]
    #: TRACS' own estimate after accounting for within-host diversity and the clock - the number its
    #: paper argues is the one not inflated by the usual pipelines.
    transmission_distance: Optional[float]
    #: Expected number of intermediate hosts, which is what a transmission claim rests on.
    expected_k: Optional[float]
    filtered_snp_distance: Optional[float]
    #: Sites with enough coverage in both samples - this backend's version of breadth.
    sites_considered: Optional[int]

    @property
    def verdict(self) -> StrainVerdict:
        """Same or different strain at TRACS' own threshold.

        Deliberately not re-derived from popANI: the two backends measure different quantities on
        different scales, and converting one threshold into the other would invent a number neither
        tool published.
        """
        distance = (
            self.filtered_snp_distance
            if self.filtered_snp_distance is not None
            else self.snp_distance
        )
        if distance is None or not self.sites_considered:
            return StrainVerdict.ABSENT
        return StrainVerdict.SAME if distance <= TRACS_SNP_THRESHOLD else StrainVerdict.DIFFERENT


def parse_tracs(path: str, samples: Sequence[str] = ()) -> List[TracsDistance]:
    """Read TRACS' `transmission_distances.csv`.

    The genome is not a column of its own, and the obvious place to look for it does not work:
    TRACS' `MSA file` column truncates the reference name at the first dot, so a run over `bin.1`
    and `bin.4` writes `combinedbin` in *both* rows. Measured, not assumed - the first version of
    this parser read that column and returned the genome `bin` for every row.

    What does carry it is the row's sample labels, which TRACS writes as `<sample>_<genome>`. Given
    the pipeline's own sample names, the genome is what remains after the matching prefix - split by
    the known sample rather than at the last underscore, because a genome called `bin_001` would
    take that split with it.

    An empty table is a legitimate result and the commonest one on shallow data: TRACS skips any
    reference with less than 25% of its genome covered, so a sample that does not carry a genome
    deeply produces no row for it rather than a row saying "far apart".
    """
    out: List[TracsDistance] = []
    if not os.path.isfile(path):
        return out
    known = sorted(samples, key=len, reverse=True)
    with open(path) as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames and TRACS_COLUMNS["snp_distance"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like TRACS' transmission_distances.csv: no %r column"
                % (path, TRACS_COLUMNS["snp_distance"])
            )
        for row in reader:
            sample_a, genome = _split_label((row.get(TRACS_COLUMNS["a"]) or ""), known)
            sample_b, _ = _split_label((row.get(TRACS_COLUMNS["b"]) or ""), known)
            out.append(
                TracsDistance(
                    genome=genome,
                    sample_a=sample_a,
                    sample_b=sample_b,
                    snp_distance=_as_float(row.get(TRACS_COLUMNS["snp_distance"])),
                    transmission_distance=_as_float(
                        row.get(TRACS_COLUMNS["transmission_distance"])
                    ),
                    expected_k=_as_float(row.get(TRACS_COLUMNS["expected_k"])),
                    filtered_snp_distance=_as_float(
                        row.get(TRACS_COLUMNS["filtered_snp_distance"])
                    ),
                    sites_considered=_as_int(row.get(TRACS_COLUMNS["sites"])),
                )
            )
    return out


def _split_label(label: str, samples: Sequence[str]) -> Tuple[str, str]:
    """`<sample>_<genome>` -> (sample, genome), split by a known sample name.

    *samples* is tried longest first, so one sample name that is a prefix of another cannot claim
    the other's rows. With no match the whole label is the sample and the genome is unknown, which
    is reported as an empty string rather than guessed at.
    """
    label = label.strip()
    for sample in samples:
        if label == sample:
            return sample, ""
        if label.startswith(sample + "_"):
            return sample, label[len(sample) + 1 :]
    return label, ""


#: Columns of the TRACS distance table.
TRACS_TABLE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "sample_a",
    "sample_b",
    "snp_distance",
    "filtered_snp_distance",
    "transmission_distance",
    "expected_intermediate_hosts",
    "sites_considered",
    "verdict",
    "permitted_language",
    "snp_threshold",
    "backend",
)


def tracs_row(distance: TracsDistance) -> Dict[str, str]:
    verdict = distance.verdict

    def num(value: Optional[float], places: int = 2) -> str:
        return "" if value is None else "%.*f" % (places, value)

    return {
        "genome": distance.genome,
        "sample_a": distance.sample_a,
        "sample_b": distance.sample_b,
        "snp_distance": num(distance.snp_distance),
        "filtered_snp_distance": num(distance.filtered_snp_distance),
        "transmission_distance": num(distance.transmission_distance, 4),
        "expected_intermediate_hosts": num(distance.expected_k, 3),
        "sites_considered": (
            "" if distance.sites_considered is None else str(distance.sites_considered)
        ),
        "verdict": verdict.value,
        "permitted_language": VERDICT_LANGUAGE[verdict.value],
        "snp_threshold": str(TRACS_SNP_THRESHOLD),
        # In the row because the two backends are not interchangeable and a table of distances with
        # no tool named cannot be compared with anything.
        "backend": "tracs",
    }


#: Printed with TRACS results.
TRACS_CAVEAT = (
    "TRACS reports a SNP distance and a transmission distance, and they are different claims: the "
    "first is between two consensus sequences, the second is its estimate after accounting for "
    "within-host diversity and the clock. Its paper's argument is that most tools overestimate the "
    "first, which is why its own threshold of %d SNPs is used here rather than one converted from "
    "another tool's scale. It also skips any reference with less than %.0f%% of its genome "
    "covered, so a shallow sample produces no row rather than a distant one"
    % (TRACS_SNP_THRESHOLD, TRACS_MIN_BREADTH * 100)
)


def compare_backends(
    instrain: Sequence[GenomeComparison], tracs: Sequence[TracsDistance]
) -> Dict[str, object]:
    """Where the two backends agree on same-versus-different, per genome and sample pair.

    Reported rather than reconciled. They measure different quantities - a population ANI and a
    consensus SNP distance - so a disagreement is information about the genome rather than a fault
    in either tool, and averaging them would produce a number neither one published.
    """

    def key(genome: str, a: str, b: str):
        return (genome, tuple(sorted((a, b))))

    first = {
        key(c.genome, c.sample_a, c.sample_b): c.verdict
        for c in instrain
        if c.verdict in (StrainVerdict.SAME, StrainVerdict.DIFFERENT)
    }
    second = {
        key(d.genome, d.sample_a, d.sample_b): d.verdict
        for d in tracs
        if d.verdict in (StrainVerdict.SAME, StrainVerdict.DIFFERENT)
    }
    shared = sorted(set(first) & set(second))
    agree = [k for k in shared if first[k] == second[k]]
    return {
        "comparable": len(shared),
        "agree": len(agree),
        "disagree": sorted(
            "%s (%s vs %s): inStrain %s, TRACS %s"
            % (genome, pair[0], pair[1], first[k].value, second[k].value)
            for k in shared
            if first[k] != second[k]
            for genome, pair in [k]
        ),
    }
