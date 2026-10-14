"""Biosynthetic gene clusters from GECCO, and how much of one a fragmented MAG can show.

A BGC is a run of adjacent genes that together make a secondary metabolite - an antibiotic, a
siderophore, a signalling molecule. It is the one functional unit in this pipeline whose definition
is *positional*: a KO is a property of a gene, a CAZy family is a property of a domain, but a gene
cluster exists only as a neighbourhood. That makes it the annotation most damaged by assembly
fragmentation, and the reason this module reports where each cluster sits on its contig rather than
only that it was found.

GECCO rather than antiSMASH as the default. It is a conditional random field over Pyrodigal genes
and PyHMMER domains - the same two libraries this package already calls in-process - needs no
multi-gigabyte database, and takes minutes per genome rather than tens of minutes. antiSMASH remains
the deeper tool when a cluster is the result rather than a count; GECCO is what makes "how many
BGCs are in these 200 MAGs" a question with an answer.

Two things the output says that a cluster count does not:

* **Where the cluster sits.** A cluster whose edge is at the edge of its contig is a cluster that
  continues into sequence nobody assembled, so its gene content is a lower bound and its predicted
  product class is based on part of the evidence. On one real CAMI2 bin, two of four clusters start
  within 400 bp of a contig end.
* **What type it is, and how sure.** GECCO assigns a product class with a probability per class and
  calls the rest ``Unknown``, which on MAGs is the commonest answer. The probabilities are kept, so
  "Unknown" is visibly a weak classification of a real cluster rather than a cluster of unknown
  reality.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: GECCO's product classes, in the order its table carries them. ``Unknown`` is not one of these -
#: it is what GECCO writes in `type` when no class passes its threshold.
BGC_TYPES: Tuple[str, ...] = (
    "alkaloid",
    "nrp",
    "polyketide",
    "ripp",
    "saccharide",
    "terpene",
)

#: What GECCO writes in `type` when none of its classes is confident enough. Kept as a constant
#: because it is a legitimate and common answer on MAGs, not a parse failure.
UNKNOWN_TYPE = "Unknown"

#: A cluster whose nearer edge is within this many base pairs of its contig's edge is reported as
#: truncated. One gene's worth of sequence, roughly: a cluster that close to the edge is one whose
#: neighbourhood was not observed, so its gene content is a floor rather than a measurement. The
#: same reasoning as the mobility window in :mod:`metawrap2.amr`, and the same honesty about what
#: missing sequence can and cannot support.
EDGE_WINDOW = 1_000

#: GECCO's columns in ``*.clusters.tsv``. Named so a schema change is a KeyError at parse time.
CLUSTER_COLUMNS = {
    "contig": "sequence_id",
    "cluster": "cluster_id",
    "start": "start",
    "end": "end",
    "average_p": "average_p",
    "max_p": "max_p",
    "type": "type",
    "proteins": "proteins",
    "domains": "domains",
}


@dataclass(frozen=True)
class GeneCluster:
    """One predicted biosynthetic gene cluster, and how much of it the assembly shows."""

    genome: str
    contig: str
    cluster_id: str
    start: int
    end: int
    #: GECCO's product class, or ``Unknown``.
    type: str
    #: Mean and maximum per-gene cluster probability from the CRF.
    average_p: Optional[float]
    max_p: Optional[float]
    #: Per-class probabilities, so an ``Unknown`` call can be read rather than just accepted.
    probabilities: Mapping[str, float]
    #: The genes in the cluster, by the same ids `annotate` uses - which is why GECCO is given this
    #: module's GFF rather than calling its own genes.
    proteins: Tuple[str, ...]
    domains: Tuple[str, ...]
    #: 0 when the contig length is unknown, which is reported as such rather than as "at the edge".
    contig_length: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)

    @property
    def genes(self) -> int:
        return len(self.proteins)

    @property
    def distance_to_contig_end(self) -> Optional[int]:
        """Base pairs between the cluster's nearer edge and its contig's nearer end.

        None when the contig length is not known. A number here and a None are different claims,
        and the second cannot support "this cluster is complete".
        """
        if not self.contig_length:
            return None
        return max(0, min(self.start - 1, self.contig_length - self.end))

    @property
    def is_truncated(self) -> bool:
        """Whether the cluster runs into unobserved sequence, as far as can be told."""
        distance = self.distance_to_contig_end
        return distance is not None and distance < EDGE_WINDOW

    @property
    def best_type(self) -> Tuple[str, float]:
        """The most probable product class and its probability, whatever GECCO called it.

        Useful precisely when `type` is ``Unknown``: the class that came closest is informative,
        and it is in the table rather than only in the probabilities a reader has to compare by eye.
        """
        if not self.probabilities:
            return (UNKNOWN_TYPE, 0.0)
        best = max(self.probabilities, key=lambda name: self.probabilities[name])
        return (best, self.probabilities[best])


def _as_int(value: Optional[str]) -> int:
    try:
        return int(float(value)) if value not in (None, "") else 0
    except ValueError:
        return 0


def _as_float(value: Optional[str]) -> Optional[float]:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def _split(value: Optional[str]) -> Tuple[str, ...]:
    """Split a ``;``-joined field, discarding the empty-string literal GECCO writes.

    With no clusters to report, GECCO writes `""` in these columns rather than leaving them blank -
    so a naive split yields a protein whose name is two quote characters.
    """
    if not value:
        return ()
    return tuple(
        part.strip().strip('"')
        for part in value.split(";")
        if part.strip() and part.strip().strip('"')
    )


def parse_clusters(
    path: str, genome: str, contig_lengths: Optional[Mapping[str, int]] = None
) -> List[GeneCluster]:
    """Read GECCO's ``*.clusters.tsv``.

    *contig_lengths* is what turns a coordinate into a statement about completeness, exactly as in
    :func:`metawrap2.amr.parse_amrfinder`: without it, a cluster starting at 392 bp could be 392 bp
    from the end of a short contig or 392 bp into a megabase, and only the first makes it a
    fragment. GECCO does not report contig lengths, so they come from the assembly the caller has.

    An empty file with only a header is a legitimate result - most genomes have few clusters or
    none - which is why `annotate` passes GECCO ``--force-tsv``: without it a genome with no
    clusters writes no table at all, and "no clusters" becomes indistinguishable from "the run
    failed".
    """
    lengths = contig_lengths or {}
    clusters: List[GeneCluster] = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and CLUSTER_COLUMNS["cluster"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like a GECCO cluster table: no %r column"
                % (path, CLUSTER_COLUMNS["cluster"])
            )
        for row in reader:
            contig = (row.get(CLUSTER_COLUMNS["contig"]) or "").strip()
            cluster_id = (row.get(CLUSTER_COLUMNS["cluster"]) or "").strip()
            # GECCO's --force-tsv writes a placeholder row for a genome with no clusters: every
            # field empty except `type`, which says Unknown, and `proteins`/`domains`, which say
            # "". Read literally that is one unclassified cluster of zero genes on no contig, and
            # it would be counted as a BGC - so the row is recognised and dropped. The flag is
            # still passed, because the alternative is no file at all and no way to tell a genome
            # with no clusters from a run that failed.
            if not cluster_id or not contig:
                continue
            probabilities = {}
            for name in BGC_TYPES:
                value = _as_float(row.get("%s_probability" % name))
                if value is not None:
                    probabilities[name] = value
            clusters.append(
                GeneCluster(
                    genome=genome,
                    contig=contig,
                    cluster_id=cluster_id,
                    start=_as_int(row.get(CLUSTER_COLUMNS["start"])),
                    end=_as_int(row.get(CLUSTER_COLUMNS["end"])),
                    type=(row.get(CLUSTER_COLUMNS["type"]) or UNKNOWN_TYPE).strip() or UNKNOWN_TYPE,
                    average_p=_as_float(row.get(CLUSTER_COLUMNS["average_p"])),
                    max_p=_as_float(row.get(CLUSTER_COLUMNS["max_p"])),
                    probabilities=probabilities,
                    proteins=_split(row.get(CLUSTER_COLUMNS["proteins"])),
                    domains=_split(row.get(CLUSTER_COLUMNS["domains"])),
                    contig_length=lengths.get(contig, 0),
                )
            )
    return clusters


#: Columns of the per-cluster table.
CLUSTER_TABLE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "contig",
    "cluster_id",
    "start",
    "end",
    "length_bp",
    "genes",
    "type",
    "most_probable_type",
    "most_probable_type_p",
    "average_gene_probability",
    "max_gene_probability",
    "distance_to_contig_end",
    "truncated",
    "domains",
    "proteins",
)


def cluster_row(cluster: GeneCluster) -> Dict[str, str]:
    best_type, best_p = cluster.best_type
    distance = cluster.distance_to_contig_end
    return {
        "genome": cluster.genome,
        "contig": cluster.contig,
        "cluster_id": cluster.cluster_id,
        "start": str(cluster.start),
        "end": str(cluster.end),
        "length_bp": str(cluster.length),
        "genes": str(cluster.genes),
        "type": cluster.type,
        # Carried even when `type` is a real class, so the two can be compared: a Saccharide call
        # at p=0.23 and one at p=0.95 are not the same finding.
        "most_probable_type": best_type,
        "most_probable_type_p": "%.4g" % best_p,
        "average_gene_probability": (
            "" if cluster.average_p is None else "%.4g" % cluster.average_p
        ),
        "max_gene_probability": "" if cluster.max_p is None else "%.4g" % cluster.max_p,
        # Blank rather than 0 when the contig length is unknown: "at the edge" is a claim, and an
        # unknown contig length cannot support it.
        "distance_to_contig_end": "" if distance is None else str(distance),
        "truncated": "yes" if cluster.is_truncated else ("" if distance is None else "no"),
        "domains": ";".join(cluster.domains),
        "proteins": ";".join(cluster.proteins),
    }


#: Printed with every BGC count. A cluster is a neighbourhood, so it is the annotation that
#: fragmentation damages most, and a count of clusters in MAGs is a count of what survived binning.
TRUNCATION_CAVEAT = (
    "a gene cluster is defined by adjacency, so it is the annotation an assembly break damages "
    "most: a cluster whose edge is within %d bp of its contig's edge continues into sequence that "
    "was never assembled, and its gene content and product class are based on the part that was. "
    "Those rows are flagged `truncated`" % EDGE_WINDOW
)


@dataclass(frozen=True)
class GenomeBgcSummary:
    """One genome's cluster counts, split the two ways that change their meaning."""

    genome: str
    clusters: int
    truncated: int
    unknown_type: int
    genes_in_clusters: int
    types: Mapping[str, int]

    @property
    def complete_clusters(self) -> int:
        return self.clusters - self.truncated


def summarise(clusters: Sequence[GeneCluster], genome: str) -> GenomeBgcSummary:
    types: Dict[str, int] = {}
    for cluster in clusters:
        types[cluster.type] = types.get(cluster.type, 0) + 1
    return GenomeBgcSummary(
        genome=genome,
        clusters=len(clusters),
        truncated=sum(1 for cluster in clusters if cluster.is_truncated),
        unknown_type=sum(1 for cluster in clusters if cluster.type == UNKNOWN_TYPE),
        genes_in_clusters=sum(cluster.genes for cluster in clusters),
        types=types,
    )


#: Columns of the per-genome summary.
SUMMARY_COLUMNS: Tuple[str, ...] = (
    "genome",
    "clusters",
    "clusters_not_at_a_contig_edge",
    "clusters_truncated",
    "clusters_unclassified",
    "genes_in_clusters",
    "types",
)


def summary_row(summary: GenomeBgcSummary) -> Dict[str, str]:
    return {
        "genome": summary.genome,
        "clusters": str(summary.clusters),
        "clusters_not_at_a_contig_edge": str(summary.complete_clusters),
        "clusters_truncated": str(summary.truncated),
        "clusters_unclassified": str(summary.unknown_type),
        "genes_in_clusters": str(summary.genes_in_clusters),
        "types": ";".join("%s=%d" % (name, count) for name, count in sorted(summary.types.items())),
    }
