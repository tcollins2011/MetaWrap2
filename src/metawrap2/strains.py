"""Within-population variation: strain sharing between samples, and diversity within a sample.

A MAG is a population consensus. Two samples can carry MAGs at 99.9% ANI and still carry different
strains, which matters for transmission, persistence or engraftment claims.

inStrain works from per-position allele frequencies in the alignments. Two outputs are used:

* popANI between two samples for one genome. Alleles are compared as populations (a position is
  shared when the allele sets overlap), so a minor variant present in both is not a difference.
  inStrain's same-strain threshold is popANI >= 99.999%, about one difference per 100 kb.
* Nucleotide diversity within one sample for one genome. A clonal population and a mix of two
  strains can have the same consensus.

popANI is computed only over positions both samples cover deeply enough, so it must be read with
breadth: two samples sharing 2% of a genome at 100% identity give popANI 100%.
:data:`MIN_PERCENT_COMPARED` applies inStrain's recommended floor of 50% of the genome compared.

Shaw & Yu's evaluation of strain-distance tools found most overestimate SNP distances, TRACS the
least. inStrain distances are best used to compare samples processed the same way, not as absolute
SNP counts.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

#: popANI at or above which two samples carry the same strain (inStrain's default).
SAME_STRAIN_POPANI = 99.999

#: Minimum percentage of the genome compared for a popANI verdict (inStrain's guidance).
MIN_PERCENT_COMPARED = 50.0

#: Minimum breadth (fraction of the genome covered at the profiler's depth floor) for a
#: within-sample diversity figure to be reported.
MIN_BREADTH = 0.5

#: Coverage below which minor variants are not called (inStrain's `min_cov` default).
MIN_COVERAGE = 5.0


class StrainVerdict(Enum):
    """What a comparison between two samples supports for one genome."""

    SAME = "same_strain"
    DIFFERENT = "different_strain"
    #: Compared, but over too little of the genome for either answer.
    UNDETERMINED = "undetermined"
    #: Not compared: the genome is not present deeply enough in both samples.
    ABSENT = "not_in_both"


#: Interpretation of each verdict, written into the output table.
VERDICT_LANGUAGE: Dict[str, str] = {
    StrainVerdict.SAME.value: (
        "the populations in these two samples are indistinguishable at inStrain's threshold "
        "(popANI >= %g%% over at least %g%% of the genome). This is consistent with transmission, "
        "persistence or a shared source, and does not distinguish between them"
        % (
            SAME_STRAIN_POPANI,
            MIN_PERCENT_COMPARED,
        )
    ),
    StrainVerdict.DIFFERENT.value: (
        "the populations differ by more than inStrain's threshold over enough of the genome to "
        "compare: different strains of the same species"
    ),
    StrainVerdict.UNDETERMINED.value: (
        "too little of the genome was compared to call same or different strains. The usual "
        "cause is low coverage in one of the two samples"
    ),
    StrainVerdict.ABSENT.value: (
        "the genome was not covered deeply enough in both samples to compare. This reflects "
        "sequencing depth; the genome may still be present in the community"
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
    #: Fraction covered at or above the profiler's depth floor, 0-1: the positions SNVs are called
    #: over. Diversity is gated on this, since it can be much lower than ``breadth``.
    breadth_min_cov: Optional[float]
    #: Breadth inStrain expects given the coverage, 0-1. See ``breadth_ratio``.
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

        Well below 1 means uneven coverage, usually reads from a relative mapping to conserved
        regions, so the genome may not be present (inStrain's diagnostic).
        """
        if not self.breadth or not self.expected_breadth:
            return None
        return self.breadth / self.expected_breadth

    @property
    def usable_breadth(self) -> Optional[float]:
        """Breadth at the depth floor if known, else overall breadth."""
        return self.breadth_min_cov if self.breadth_min_cov is not None else self.breadth

    @property
    def diversity_is_reportable(self) -> bool:
        """Whether breadth and coverage meet MIN_BREADTH and MIN_COVERAGE."""
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
    #: Consensus ANI: major alleles only. Stricter than popANI; reported because a difference
    #: between the two is informative.
    con_ani: Optional[float]
    #: Percentage of the genome compared: the positions popANI is computed over.
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

        A rate, so comparisons over different breadths can be compared. Not an absolute SNP
        distance; see the module docstring.
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

    Always multiplied by 100, not inferred from magnitude: percent_compared ranges from 0.001 to
    1.0, so a magnitude rule would fail on it (as in `parse_singlem_fraction`).
    """
    number = _as_float(value)
    return None if number is None else number * 100.0


def _as_int(value: object) -> Optional[int]:
    number = _as_float(value)
    return None if number is None else int(number)


def _strip_extension(name: str) -> str:
    """A genome name without its FASTA extension.

    Not ``os.path.splitext``, which would turn metaWRAP bin names like ``bin.13`` into ``bin``.
    """
    for extension in (".fasta", ".fa", ".fna", ".fas"):
        if name.endswith(extension):
            return name[: -len(extension)]
    return name


#: inStrain's profile column names.
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

#: inStrain's compare table stores these as fractions, while its documentation quotes popANI as a
#: percentage ("99.999%"). They are multiplied by 100 on parsing to match the thresholds above.
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

    Sample names are returned as inStrain writes them; the caller maps them to the pipeline's
    sample names with :func:`relabel`.
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

    `compare` labels rows with the BAM basename (``ERR011347.bam``); the extension is stripped so
    the caller's sample names can be mapped back.
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
        # Blank, not 0, when coverage is too low: zero diversity would mean a clonal population.
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
        # Six decimal places, since the threshold is 99.999%.
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
    "popANI is computed only over positions both samples covered deeply enough, so read it with "
    "`percent_genome_compared`: two samples sharing 2%% of a genome at 100%% identity have popANI "
    "100%%, which is not evidence of a shared strain. Comparisons over less than %g%% of the "
    "genome are reported as undetermined" % MIN_PERCENT_COMPARED
)

#: Printed with every nucleotide-diversity figure.
DIVERSITY_CAVEAT = (
    "nucleotide diversity is measured over the positions covered at or above %gx, so for a "
    "genome at low coverage it uses only part of the genome. Genomes below %g breadth are "
    "reported without a diversity value"
    % (
        MIN_COVERAGE,
        MIN_BREADTH,
    )
)


def comparable_genomes(
    profiles: Sequence[GenomeProfile], samples: Sequence[str]
) -> Dict[Tuple[str, str], int]:
    """How many genomes each *pair* of samples could even be compared for.

    Counted per pair, as inStrain compares pairs; the intersection over all samples can be much
    smaller (0 on three gut samples with 42 pairwise comparisons).

    A genome is comparable when it is reportable (see ``diversity_is_reportable``) in both samples.
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

#: TRACS' default SNP threshold for a transmission pair. Thresholds from other tools do not
#: transfer, since most tools overestimate SNP distances.
TRACS_SNP_THRESHOLD = 10

#: Fraction of a reference genome TRACS requires to be covered before using it (its own floor).
#: A shallow sample therefore produces no rows.
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
    #: TRACS' estimate after accounting for within-host diversity and the molecular clock.
    transmission_distance: Optional[float]
    #: Expected number of intermediate hosts.
    expected_k: Optional[float]
    filtered_snp_distance: Optional[float]
    #: Sites with enough coverage in both samples (this backend's measure of breadth).
    sites_considered: Optional[int]

    @property
    def verdict(self) -> StrainVerdict:
        """Same or different strain at TRACS' own threshold.

        Not derived from popANI: the backends measure different quantities on different scales.
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

    There is no genome column, and the `MSA file` column truncates reference names at the first
    dot (`bin.1` and `bin.4` both give `combinedbin`). The genome is taken from the sample labels,
    which TRACS writes as `<sample>_<genome>`, by removing a known sample name prefix. Splitting at
    the last underscore would break genome names such as `bin_001`.

    An empty table is a valid result, common on shallow data: TRACS skips references with less than
    25% of the genome covered.
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

    *samples* is tried longest first, so a sample name that is a prefix of another does not match
    the other's rows. Without a match, the whole label is the sample and the genome is empty.
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
        # The backends are not interchangeable, so each row names its tool.
        "backend": "tracs",
    }


#: Printed with TRACS results.
TRACS_CAVEAT = (
    "TRACS reports a SNP distance and a transmission distance. The SNP distance is between two "
    "consensus sequences; the transmission distance is TRACS' estimate after accounting for "
    "within-host diversity and the molecular clock. The TRACS paper reports that most tools "
    "overestimate the SNP distance, so its own threshold of %d SNPs is used here. TRACS skips "
    "any reference with less than %.0f%% of its genome covered, so a shallow sample produces no "
    "row" % (TRACS_SNP_THRESHOLD, TRACS_MIN_BREADTH * 100)
)


def compare_backends(
    instrain: Sequence[GenomeComparison], tracs: Sequence[TracsDistance]
) -> Dict[str, object]:
    """Where the two backends agree on same-versus-different, per genome and sample pair.

    Disagreements are listed, not reconciled: the backends measure different quantities
    (population ANI and consensus SNP distance).
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
