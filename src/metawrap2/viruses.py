"""Viral contigs: what they are, how much of a genome each one is, and which populations they form.

MetaWrap2 bins prokaryotes and discards the phage, in datasets where phage is usually the most
abundant thing present. This is the viral counterpart of ``binning`` plus ``bin_refinement``, and it
needs both halves for the output to mean anything: geNomad says "this contig is a virus", CheckV says
"and it is 40% of one".

Four things this file owns, in the order a result depends on them:

**Provirus trimming is not optional.** An integrated prophage sits inside a bacterial contig. Reported
untrimmed it is one chimeric sequence that is neither a virus nor a chromosome, and its completeness
is meaningless. CheckV trims the host flanks; the trimmed coordinates are carried so a reader can see
how much was cut and from where.

**Completeness is estimated against a reference, so novelty looks like fragmentation.** CheckV's
confidence in its own estimate is part of the answer, not a detail: a contig with no close relative
in the database is reported as a fragment whether or not it is one. The confidence level is on every
row and :attr:`ViralContig.completeness_is_reliable` exists so a caller cannot silently treat a
low-confidence estimate as a measurement.

**Quality tiers follow MIUViG, which is a real standard**, unlike the 70/10 threshold this project
reports prokaryotic bins at. High-quality is >=90% complete; a genome fragment is anything under 50%.
Naming them after the standard means a count here is comparable with a count in a paper.

**A vOTU is the unit, not a contig.** Two contigs of the same viral population are one finding.
Clustering at 95% ANI over 85% of the shorter sequence is the convention the viral catalogues use -
the same shape of work ``dereplicate`` does for genomes, with different thresholds.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from .io.seqio import iter_fasta

if TYPE_CHECKING:  # the runtime import is local, to keep a library from importing a module
    from .modules.dereplicate import Edge

#: vOTU clustering thresholds, the convention the viral genome catalogues are built on (MIUViG;
#: Roux et al.). Stated rather than inherited from a tool's defaults so a population count here
#: means the same thing as one in a paper.
VOTU_ANI = 95.0
VOTU_COVERAGE = 85.0

#: Minimum length geNomad is worth asking about. Below this a "virus" call is a fragment of one and
#: CheckV cannot place it; the number is the floor the viral literature generally applies rather than
#: a measurement of ours, which is why it is a parameter on the module.
MIN_VIRAL_LENGTH = 2_000


class Tier(Enum):
    """MIUViG quality tiers for an uncultivated viral genome.

    The thresholds are the standard's, and unlike the prokaryotic 70/10 tier this project reports,
    they are not ours to choose. ``complete`` is CheckV's own determination from direct terminal
    repeat evidence rather than a completeness cutoff - a circular genome is complete because its
    ends join, not because an estimate crossed 100%.
    """

    COMPLETE = "complete"
    HIGH_QUALITY = "high_quality"
    MEDIUM_QUALITY = "medium_quality"
    LOW_QUALITY = "low_quality"
    UNDETERMINED = "undetermined"


#: How much CheckV trusts its own completeness estimate, best first, as it writes them in
#: ``completeness_method``: ``AAI-based (high-confidence)``, ``AAI-based (medium-confidence)``,
#: ``HMM-based (lower-bound)``, or ``NA``.
#:
#: ``lower-bound`` is the one that matters and the reason this is parsed rather than ignored. An
#: HMM-based figure is not a weak estimate of completeness, it is a *floor*: the genome is at least
#: this complete and may be much more. On a real assembly 35 of 161 viral contigs land there, so
#: treating those numbers as estimates would systematically under-report completeness for exactly
#: the viruses with no close relative in the database - the novel ones.
RELIABLE_CONFIDENCE: Tuple[str, ...] = ("high-confidence", "medium-confidence")
LOWER_BOUND = "lower-bound"
NO_ESTIMATE = "none"

#: CheckV's own ``checkv_quality``, which is its MIUViG determination. Read rather than re-derived
#: from the completeness number: CheckV applies the standard's rules including the terminal-repeat
#: evidence that makes a genome "Complete", and a second implementation of those rules here could
#: only disagree with it.
CHECKV_TIERS: Mapping[str, Tier] = {}  # populated below, once Tier exists


@dataclass(frozen=True)
class ViralContig:
    """One viral contig, as geNomad and CheckV together describe it."""

    genome: str
    contig: str
    length: int
    #: geNomad's call: virus, plasmid or chromosome. Only the first reaches this table, but the
    #: field is kept because a plasmid call on the same run is what the `amr` module's mobility
    #: work consumes, and the two should not need two passes over the same contigs.
    topology: str
    taxonomy: str
    #: geNomad's virus score, 0-1.
    virus_score: Optional[float]
    n_genes: int
    n_viral_genes: int
    n_host_genes: int
    #: CheckV completeness, 0-100, or None where CheckV would not estimate one.
    completeness: Optional[float]
    completeness_confidence: str
    contamination: Optional[float]
    #: True when CheckV found terminal repeats or a circular topology, which is direct evidence of a
    #: complete genome rather than an estimate.
    is_complete_genome: bool
    #: Set when the contig was a provirus and CheckV trimmed host sequence off it.
    provirus: bool
    trimmed_length: Optional[int]

    #: CheckV's own MIUViG tier for this contig, when CheckV scored it.
    checkv_tier: Optional[Tier] = None

    @property
    def completeness_is_reliable(self) -> bool:
        """Whether the completeness estimate is an estimate at all, rather than a floor.

        Exists so that treating a lower bound as a measurement is a decision rather than an
        oversight. An HMM-based figure says "at least this complete"; a novel virus with no close
        relative gets one, and reporting it as a fragment is a statement about the database rather
        than about the virus.
        """
        return self.completeness_confidence in RELIABLE_CONFIDENCE

    @property
    def host_gene_fraction(self) -> Optional[float]:
        """Share of this contig's genes that look bacterial, as a chimerism signal.

        A high value on something called a virus usually means an untrimmed provirus or a
        misclassified chromosomal region, which is the viral equivalent of what GUNC catches in a
        prokaryotic bin.
        """
        if not self.n_genes:
            return None
        return self.n_host_genes / float(self.n_genes)


def tier_of(contig: ViralContig) -> Tier:
    """MIUViG tier for one contig: CheckV's determination where there is one, ours otherwise.

    CheckV already applies the standard, including the terminal-repeat evidence that makes a genome
    "Complete" rather than merely 100% estimated, so re-deriving the tier from the completeness
    number could only disagree with the tool that produced it. The fallback is for a contig CheckV
    never scored - which happens, and must not silently become low quality.
    """
    if contig.checkv_tier is not None:
        return contig.checkv_tier
    if contig.is_complete_genome:
        return Tier.COMPLETE
    if contig.completeness is None or not contig.completeness_is_reliable:
        return Tier.UNDETERMINED
    if contig.completeness >= 90.0:
        return Tier.HIGH_QUALITY
    if contig.completeness >= 50.0:
        return Tier.MEDIUM_QUALITY
    return Tier.LOW_QUALITY


CHECKV_TIERS = {
    "complete": Tier.COMPLETE,
    "high-quality": Tier.HIGH_QUALITY,
    "medium-quality": Tier.MEDIUM_QUALITY,
    "low-quality": Tier.LOW_QUALITY,
    "not-determined": Tier.UNDETERMINED,
}


def parse_confidence(completeness_method: str) -> str:
    """The confidence CheckV reports inside ``completeness_method``.

    ``AAI-based (high-confidence)`` -> ``high-confidence``; ``HMM-based (lower-bound)`` ->
    ``lower-bound``; ``NA`` or empty -> ``none``. Parsed from the parenthesis rather than matched
    against whole strings, so a new method name with a known confidence still works.
    """
    found = re.search(r"\(([^)]+)\)", completeness_method or "")
    if found:
        return found.group(1).strip().lower()
    return NO_ESTIMATE


def _as_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _as_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_genomad(path: str) -> Dict[str, Dict[str, str]]:
    """geNomad's ``*_virus_summary.tsv``, keyed by its ``seq_name``.

    Note that geNomad renames a provirus ``<contig>|provirus_<start>_<end>``, so the key here is not
    always a contig id from the assembly. :func:`merge` keeps both the geNomad name and the original
    contig, because the first is what CheckV will echo and the second is what the assembly knows.
    """
    rows: Dict[str, Dict[str, str]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and "seq_name" not in reader.fieldnames:
            raise KeyError("%s is not a geNomad virus summary (no seq_name column)" % path)
        for row in reader:
            rows[row["seq_name"]] = row
    return rows


def parse_checkv(path: str) -> Dict[str, Dict[str, str]]:
    """CheckV's ``quality_summary.tsv``, keyed by ``contig_id``."""
    rows: Dict[str, Dict[str, str]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and "contig_id" not in reader.fieldnames:
            raise KeyError("%s is not a CheckV quality summary (no contig_id column)" % path)
        for row in reader:
            rows[row["contig_id"]] = row
    return rows


def original_contig(seq_name: str) -> str:
    """The assembly's contig id, with geNomad's provirus suffix removed if there is one."""
    return seq_name.split("|", 1)[0]


def merge(
    genome: str,
    genomad_rows: Mapping[str, Mapping[str, str]],
    checkv_rows: Mapping[str, Mapping[str, str]],
) -> List[ViralContig]:
    """Join geNomad's calls to CheckV's quality, one row per viral sequence.

    A sequence geNomad called and CheckV did not score still appears, with no completeness and
    ``undetermined`` quality. Dropping it would silently shrink the viral set to whatever CheckV
    managed, and "CheckV had nothing to compare this to" is a finding about a novel virus rather
    than a reason to discard it.
    """
    out: List[ViralContig] = []
    for seq_name, grow in genomad_rows.items():
        crow = checkv_rows.get(seq_name, {})
        completeness = _as_float(crow.get("completeness", ""))
        checkv_quality = (crow.get("checkv_quality", "") or "").strip().lower()
        provirus = (crow.get("provirus", "") or "").strip().lower()
        out.append(
            ViralContig(
                genome=genome,
                contig=original_contig(seq_name),
                length=_as_int(grow.get("length", "")),
                topology=grow.get("topology", "") or "",
                taxonomy=grow.get("taxonomy", "") or "",
                virus_score=_as_float(grow.get("virus_score", "")),
                n_genes=_as_int(grow.get("n_genes", "")),
                n_viral_genes=_as_int(grow.get("n_hallmarks", "")),
                n_host_genes=_as_int(crow.get("host_genes", "")),
                completeness=completeness,
                completeness_confidence=parse_confidence(crow.get("completeness_method", "")),
                contamination=_as_float(crow.get("contamination", "")),
                is_complete_genome=checkv_quality == "complete",
                # CheckV writes a literal "Yes"/"No" here; geNomad's topology is a different thing
                # (terminal repeats) and conflating them marked every linear contig a provirus.
                provirus=provirus == "yes",
                trimmed_length=_as_int(crow.get("proviral_length", ""), 0) or None,
                checkv_tier=CHECKV_TIERS.get(checkv_quality),
            )
        )
    return out


#: Columns of the viral contig table.
VIRAL_COLUMNS: Tuple[str, ...] = (
    "genome",
    "contig",
    "length",
    "tier",
    "topology",
    "taxonomy",
    "virus_score",
    "n_genes",
    "n_hallmark_genes",
    "n_host_genes",
    "host_gene_fraction",
    "completeness",
    "completeness_confidence",
    "completeness_is_reliable",
    "contamination",
    "is_complete_genome",
    "provirus",
    "trimmed_length",
)


def viral_row(contig: ViralContig) -> Dict[str, object]:
    fraction = contig.host_gene_fraction
    return {
        "genome": contig.genome,
        "contig": contig.contig,
        "length": contig.length,
        "tier": tier_of(contig).value,
        "topology": contig.topology,
        "taxonomy": contig.taxonomy,
        "virus_score": "" if contig.virus_score is None else "%.4f" % contig.virus_score,
        "n_genes": contig.n_genes,
        "n_hallmark_genes": contig.n_viral_genes,
        "n_host_genes": contig.n_host_genes,
        "host_gene_fraction": "" if fraction is None else "%.4f" % fraction,
        "completeness": "" if contig.completeness is None else "%.2f" % contig.completeness,
        "completeness_confidence": contig.completeness_confidence,
        "completeness_is_reliable": int(contig.completeness_is_reliable),
        "contamination": "" if contig.contamination is None else "%.2f" % contig.contamination,
        "is_complete_genome": int(contig.is_complete_genome),
        "provirus": int(contig.provirus),
        "trimmed_length": "" if contig.trimmed_length is None else contig.trimmed_length,
    }


def summarise(contigs: Sequence[ViralContig]) -> Dict[str, int]:
    """Counts per tier, plus the two that qualify them."""
    out: Dict[str, int] = {tier.value: 0 for tier in Tier}
    out["proviruses"] = 0
    out["unreliable_completeness"] = 0
    for contig in contigs:
        out[tier_of(contig).value] += 1
        if contig.provirus:
            out["proviruses"] += 1
        if not contig.completeness_is_reliable:
            out["unreliable_completeness"] += 1
    return out


#: Minimum spacer length to use for host linking. A CRISPR spacer is 25-40 bp; anything much
#: shorter is either a mis-called array or short enough to occur by chance in a viral genome, and an
#: exact match of it is not evidence of anything.
MIN_SPACER_LENGTH = 20


@dataclass(frozen=True)
class CrisprArray:
    """One CRISPR array found in a host genome, as cctyper reports it."""

    host: str
    contig: str
    array: str
    start: int
    end: int
    repeats: int
    #: cctyper's own verdict on the array. False means it kept the array but does not trust it.
    trusted: bool
    #: The subtype its repeat classifier predicted, which is informative about the host, not the
    #: virus - carried through so a host link can be read beside the system that recorded it.
    subtype: str
    subtype_probability: Optional[float]


def read_crispr_arrays(outdir: str, host: str) -> List[CrisprArray]:
    """Parse cctyper's ``crisprs_all.tab``. Returns [] when it found nothing, which is common.

    ``crisprs_all.tab`` rather than the per-category files: it holds every array *and* a ``Trusted``
    column, so one parse gives both the full set and cctyper's own judgement of each. Reading
    ``crisprs_near_cas.tab`` and ``crisprs_orphan.tab`` instead would reconstruct that column from
    which file a row came from, with more code and one more thing to get wrong.
    """
    path = os.path.join(outdir, "crisprs_all.tab")
    if not os.path.isfile(path):
        return []
    arrays: List[CrisprArray] = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            probability = row.get("Subtype_probability", "")
            arrays.append(
                CrisprArray(
                    host=host,
                    contig=row.get("Contig", ""),
                    array=row.get("CRISPR", ""),
                    start=int(float(row.get("Start") or 0)),
                    end=int(float(row.get("End") or 0)),
                    repeats=int(float(row.get("N_repeats") or 0)),
                    trusted=str(row.get("Trusted", "")).strip().lower() == "true",
                    subtype=row.get("Subtype", ""),
                    subtype_probability=float(probability) if probability else None,
                )
            )
    return arrays


def read_spacers(outdir: str, trusted_only: bool = True) -> Dict[str, List[str]]:
    """array id -> its spacer sequences, from cctyper's ``spacers/<array>.fa``.

    Untrusted arrays are excluded by default. They are cctyper's own low-confidence calls, and a
    spurious array contributes spacers that are really just genomic sequence - which then match
    whatever they came from and produce host links for nothing. Reachable with
    ``trusted_only=False`` for anyone who would rather filter later.
    """
    spacer_dir = os.path.join(outdir, "spacers")
    if not os.path.isdir(spacer_dir):
        return {}
    allowed: Optional[set] = None
    if trusted_only:
        allowed = {a.array for a in read_crispr_arrays(outdir, host="") if a.trusted}
    out: Dict[str, List[str]] = {}
    for name in sorted(os.listdir(spacer_dir)):
        if not name.endswith((".fa", ".fasta", ".fna")):
            continue
        array = name.rsplit(".", 1)[0]
        if allowed is not None and array not in allowed:
            continue
        sequences = [sequence for _header, sequence in iter_fasta(os.path.join(spacer_dir, name))]
        if sequences:
            out[array] = sequences
    return out


@dataclass(frozen=True)
class HostLink:
    """One piece of evidence that a host was infected by a virus.

    A spacer is a record of a past infection, so this is direct evidence rather than a prediction
    from composition - which is the reason to prefer it over the alternatives even though it only
    ever covers a fraction of the viruses found.
    """

    host: str
    array: str
    subtype: str
    viral_contig: str
    spacer: str
    #: Whether the match was to the reverse complement. Recorded because it says nothing about the
    #: biology and everything about whether a match is a duplicate of another one.
    reverse: bool

    @property
    def spacer_length(self) -> int:
        return len(self.spacer)


def host_links(
    spacers_by_host: Mapping[str, Mapping[str, Sequence[str]]],
    viral_sequences: Mapping[str, str],
    subtypes: Optional[Mapping[str, str]] = None,
    min_spacer_length: int = MIN_SPACER_LENGTH,
) -> List[HostLink]:
    """Every exact spacer match between a host's CRISPR arrays and a viral sequence.

    Exact matching only, in both orientations. A spacer is 25-40 bp and approximate matching at
    that length needs a model of how many mismatches still count as the same protospacer - which is
    a measurement nobody here has made, and guessing at it would turn direct evidence into a
    prediction while keeping the language of evidence.
    """
    subtypes = subtypes or {}
    links: List[HostLink] = []
    for host, arrays in spacers_by_host.items():
        for array, spacers in arrays.items():
            for spacer in spacers:
                if len(spacer) < min_spacer_length:
                    continue
                reverse_spacer = _reverse_complement(spacer)
                for contig, sequence in viral_sequences.items():
                    forward = spacer in sequence
                    if not forward and reverse_spacer not in sequence:
                        continue
                    links.append(
                        HostLink(
                            host=host,
                            array=array,
                            subtype=subtypes.get(array, ""),
                            viral_contig=contig,
                            spacer=spacer,
                            reverse=not forward,
                        )
                    )
    return links


#: Columns of the host-link table.
HOST_LINK_COLUMNS: Tuple[str, ...] = (
    "host",
    "viral_contig",
    "crispr_array",
    "crispr_subtype",
    "spacer",
    "spacer_length",
    "orientation",
)


def host_link_row(link: HostLink) -> Dict[str, str]:
    return {
        "host": link.host,
        "viral_contig": link.viral_contig,
        "crispr_array": link.array,
        "crispr_subtype": link.subtype,
        "spacer": link.spacer,
        "spacer_length": str(link.spacer_length),
        "orientation": "reverse" if link.reverse else "forward",
    }


def crispr_host_links(
    spacers: Mapping[str, Sequence[str]],
    viral_sequences: Mapping[str, str],
) -> Iterator[Tuple[str, str, str]]:
    """(host genome, viral contig, spacer) for every exact spacer match in a viral sequence.

    Host prediction for free. The prokaryotic bins are already in hand and ``cctyper`` reads CRISPR
    arrays off them; a spacer is a record of a past infection, so a spacer from bin X found in vOTU
    Y is direct evidence that X is a host of Y. Exact matching only - a spacer is 25-40 bp and
    approximate matching at that length needs a model of what counts as a match, which is a
    measurement nobody here has made.
    """
    for host, host_spacers in spacers.items():
        for spacer in host_spacers:
            if len(spacer) < 20:
                continue
            for contig, sequence in viral_sequences.items():
                if spacer in sequence or _reverse_complement(spacer) in sequence:
                    yield host, contig, spacer


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


# ── vOTU clustering ──────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Votu:
    """One viral operational taxonomic unit: a population, not a contig."""

    name: str
    representative: str
    members: Tuple[str, ...]
    representative_length: int

    @property
    def size(self) -> int:
        return len(self.members)


def choose_representative(members: Sequence[str], lengths: Mapping[str, int]) -> str:
    """The longest member, which is the viral convention.

    Not the highest-quality one, as ``dereplicate`` picks for genomes: CheckV's completeness is an
    estimate against a reference database, so on a cluster of equally novel contigs it is nearly
    uninformative, while length is observed. Ties go to the lexicographically first name so two runs
    over the same contigs give the same catalogue.
    """
    return min(members, key=lambda name: (-lengths.get(name, 0), name))


def build_votus(
    names: Sequence[str],
    edges: Iterable[Edge],
    lengths: Mapping[str, int],
    *,
    ani: float = VOTU_ANI,
    min_aligned_fraction: float = VOTU_COVERAGE,
) -> List[Votu]:
    """Cluster viral contigs into vOTUs at *ani* identity over *min_aligned_fraction* of the shorter.

    The clustering itself is :func:`metawrap2.modules.dereplicate.cluster_genomes` - single linkage
    over the pairs that clear both thresholds - because it is the same operation with different
    numbers, and a second implementation would be a second set of bugs. What differs is the
    thresholds (95/85 rather than 95/10) and how the representative is chosen.
    """
    # Imported here rather than at module scope: this is a library and dereplicate is a module, so
    # a top-level import would have the layering backwards even though nothing is circular today.
    from .modules.dereplicate import cluster_genomes

    clusters = cluster_genomes(names, edges, ani=ani, min_aligned_fraction=min_aligned_fraction)
    votus: List[Votu] = []
    for index, members in enumerate(clusters, start=1):
        representative = choose_representative(members, lengths)
        votus.append(
            Votu(
                name="vOTU_%04d" % index,
                representative=representative,
                members=tuple(members),
                representative_length=lengths.get(representative, 0),
            )
        )
    return votus


#: Columns of the vOTU table.
VOTU_COLUMNS: Tuple[str, ...] = (
    "votu",
    "representative",
    "representative_length",
    "n_members",
    "members",
    "ani_threshold",
    "aligned_fraction_threshold",
)


def votu_row(
    votu: Votu, ani: float = VOTU_ANI, coverage: float = VOTU_COVERAGE
) -> Dict[str, object]:
    return {
        "votu": votu.name,
        "representative": votu.representative,
        "representative_length": votu.representative_length,
        "n_members": votu.size,
        "members": ",".join(votu.members),
        "ani_threshold": ani,
        "aligned_fraction_threshold": coverage,
    }
