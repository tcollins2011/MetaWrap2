"""Viral contigs: classification, completeness, and clustering into populations.

The viral counterpart of ``binning`` plus ``bin_refinement``: geNomad identifies viral contigs and
CheckV estimates how complete each one is.

**Provirus trimming.** An integrated prophage sits inside a bacterial contig, and untrimmed it is a
chimera whose completeness has no meaning. CheckV trims the host flanks; the trimmed coordinates
are carried through.

**Completeness is estimated against a reference**, so a contig with no close relative in the
database is reported as a fragment whether or not it is one. CheckV's confidence level is on every
row, and :attr:`ViralContig.completeness_is_reliable` separates estimates from lower bounds.

**Quality tiers follow MIUViG**: high-quality is >=90% complete; a genome fragment is under 50%.
This makes counts comparable with published ones.

**The unit is the vOTU.** Contigs are clustered at 95% ANI over 85% of the shorter sequence, the
convention of the viral catalogues; ``dereplicate`` does the same for genomes with different
thresholds.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from .io.seqio import iter_fasta

if TYPE_CHECKING:  # imported locally at runtime; see build_votus
    from .modules.dereplicate import Edge

#: vOTU clustering thresholds, the convention viral genome catalogues use (MIUViG; Roux et al.).
VOTU_ANI = 95.0
VOTU_COVERAGE = 85.0

#: Minimum contig length sent to geNomad. Shorter "virus" calls are fragments CheckV cannot place.
#: The value is the floor commonly used in the viral literature, not a local measurement.
MIN_VIRAL_LENGTH = 2_000


class Tier(Enum):
    """MIUViG quality tiers for an uncultivated viral genome.

    The thresholds are the standard's. ``complete`` is CheckV's determination from direct terminal
    repeat evidence, not a completeness cutoff.
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
#: ``lower-bound`` (HMM-based) is a floor: the genome is at least this complete and may be more.
#: On a real assembly 35 of 161 viral contigs had one. Treating these as estimates under-reports
#: completeness for viruses with no close database relative.
RELIABLE_CONFIDENCE: Tuple[str, ...] = ("high-confidence", "medium-confidence")
LOWER_BOUND = "lower-bound"
NO_ESTIMATE = "none"

#: CheckV's ``checkv_quality``, its MIUViG determination. Read directly because CheckV applies the
#: standard's rules, including the terminal-repeat evidence for "Complete".
CHECKV_TIERS: Mapping[str, Tier] = {}  # populated below, once Tier exists


@dataclass(frozen=True)
class ViralContig:
    """One viral contig, as geNomad and CheckV together describe it."""

    genome: str
    contig: str
    length: int
    #: geNomad's call: virus, plasmid or chromosome. Only viruses reach this table; the field is
    #: kept because the `amr` module's mobility work uses the plasmid calls from the same run.
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
    #: True when CheckV found terminal repeats or a circular topology, which is direct evidence of
    #: a complete genome.
    is_complete_genome: bool
    #: Set when the contig was a provirus and CheckV trimmed host sequence off it.
    provirus: bool
    trimmed_length: Optional[int]

    #: CheckV's own MIUViG tier for this contig, when CheckV scored it.
    checkv_tier: Optional[Tier] = None

    @property
    def completeness_is_reliable(self) -> bool:
        """Whether the completeness figure is an estimate rather than a lower bound.

        HMM-based figures are lower bounds ("at least this complete"), which novel viruses with no
        close relative typically get.
        """
        return self.completeness_confidence in RELIABLE_CONFIDENCE

    @property
    def host_gene_fraction(self) -> Optional[float]:
        """Share of this contig's genes that look bacterial, as a chimerism signal.

        A high value on a viral contig usually means an untrimmed provirus or a misclassified
        chromosomal region, as GUNC detects in a prokaryotic bin.
        """
        if not self.n_genes:
            return None
        return self.n_host_genes / float(self.n_genes)


def tier_of(contig: ViralContig) -> Tier:
    """MIUViG tier for one contig: CheckV's where it has one, derived here otherwise.

    CheckV's determination includes terminal-repeat evidence, so it is used as is. The fallback
    covers contigs CheckV did not score, which are not treated as low quality.
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
    ``lower-bound``; ``NA`` or empty -> ``none``. Parsed from the parentheses, so a new method name
    with a known confidence still works.
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

    geNomad names a provirus ``<contig>|provirus_<start>_<end>``, so the key is not always an
    assembly contig id. :func:`merge` keeps both the geNomad name (which CheckV echoes) and the
    original contig.
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

    A sequence geNomad called and CheckV did not score is kept, with no completeness and
    ``undetermined`` quality, since it may be a novel virus.
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
                # CheckV writes "Yes"/"No" here. geNomad's topology (terminal repeats) is a
                # different field; using it marked every linear contig a provirus.
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


#: Minimum spacer length for host linking. CRISPR spacers are 25-40 bp; much shorter ones are
#: mis-called arrays or short enough to match a viral genome by chance.
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
    #: The subtype predicted by its repeat classifier. Describes the host's system, and is carried
    #: so a host link can be read alongside it.
    subtype: str
    subtype_probability: Optional[float]


def read_crispr_arrays(outdir: str, host: str) -> List[CrisprArray]:
    """Parse cctyper's ``crisprs_all.tab``. Returns [] when it found nothing, which is common.

    ``crisprs_all.tab`` holds every array plus a ``Trusted`` column, so one parse gives the full
    set and cctyper's judgement of each, without combining ``crisprs_near_cas.tab`` and
    ``crisprs_orphan.tab``.
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

    Untrusted arrays (cctyper's low-confidence calls) are excluded by default: a spurious array's
    spacers are genomic sequence that matches its source and produces false host links. Pass
    ``trusted_only=False`` to keep them.
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

    A spacer records a past infection, so this is direct evidence rather than a compositional
    prediction, though it covers only a fraction of the viruses found.
    """

    host: str
    array: str
    subtype: str
    viral_contig: str
    spacer: str
    #: Whether the match was to the reverse complement. Used to identify duplicate matches.
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

    Exact matching only, in both orientations. Approximate matching of 25-40 bp spacers would need
    a validated mismatch tolerance, which has not been measured.
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

    ``cctyper`` reads CRISPR arrays off the prokaryotic bins. A spacer from bin X found in vOTU Y is
    direct evidence that X is a host of Y. Exact matching only; approximate matching of 25-40 bp
    spacers would need a mismatch tolerance that has not been measured.
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
    """One viral operational taxonomic unit: a population of contigs."""

    name: str
    representative: str
    members: Tuple[str, ...]
    representative_length: int

    @property
    def size(self) -> int:
        return len(self.members)


def choose_representative(members: Sequence[str], lengths: Mapping[str, int]) -> str:
    """The longest member, as is conventional for viruses.

    ``dereplicate`` picks the highest-quality genome, but CheckV completeness is a reference-based
    estimate and says little within a cluster of equally novel contigs, while length is observed.
    Ties go to the lexicographically first name, so the result is deterministic.
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

    The clustering is :func:`metawrap2.modules.dereplicate.cluster_genomes` (single linkage over
    pairs that pass both thresholds). The thresholds differ (95/85 instead of 95/10), as does the
    choice of representative.
    """
    # Imported here because this library should not depend on a module at import time.
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
