"""Biosynthetic gene clusters (BGCs) from GECCO, with their position on the contig.

A BGC is a run of adjacent genes that together make a secondary metabolite (an antibiotic, a
siderophore, a signalling molecule). Unlike a KO or a CAZy family, a cluster is defined by gene
adjacency, so assembly fragmentation affects it most. This module therefore reports where each
cluster sits on its contig.

GECCO is the default instead of antiSMASH. It is a conditional random field over Pyrodigal genes
and PyHMMER domains (libraries this package already calls in-process), needs no multi-gigabyte
database, and takes minutes per genome instead of tens of minutes. antiSMASH is the more detailed
tool for studying individual clusters; GECCO is practical for counting BGCs across hundreds of MAGs.

The output adds two things to a cluster count:

* **Position.** A cluster at the edge of its contig continues into unassembled sequence, so its gene
  content is a lower bound and its product class is based on partial evidence. On one CAMI2 bin, two
  of four clusters start within 400 bp of a contig end.
* **Type and confidence.** GECCO gives a probability per product class and calls the cluster
  ``Unknown`` when none passes its threshold, which is the commonest result on MAGs. The
  probabilities are kept so an ``Unknown`` call can be inspected.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: GECCO's product classes, in table order. ``Unknown`` is not a class; GECCO writes it in `type`
#: when no class passes its threshold.
BGC_TYPES: Tuple[str, ...] = (
    "alkaloid",
    "nrp",
    "polyketide",
    "ripp",
    "saccharide",
    "terpene",
)

#: What GECCO writes in `type` when no class passes its threshold. A valid, common result on MAGs.
UNKNOWN_TYPE = "Unknown"

#: A cluster whose nearer edge is within this many base pairs of its contig's edge is reported as
#: truncated. Roughly one gene of sequence; past that the neighbourhood was not observed, so the
#: gene content is a lower bound. Same reasoning as the mobility window in :mod:`metawrap2.amr`.
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
    #: Per-class probabilities, so an ``Unknown`` call can be inspected.
    probabilities: Mapping[str, float]
    #: The genes in the cluster, by the ids `annotate` uses. GECCO is given this module's GFF
    #: instead of calling its own genes so the ids match.
    proteins: Tuple[str, ...]
    domains: Tuple[str, ...]
    #: 0 when the contig length is unknown; such clusters are not reported as truncated.
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

        None when the contig length is not known, in which case completeness cannot be judged.
        """
        if not self.contig_length:
            return None
        return max(0, min(self.start - 1, self.contig_length - self.end))

    @property
    def is_truncated(self) -> bool:
        """Whether the cluster is within EDGE_WINDOW of a contig end (False if length unknown)."""
        distance = self.distance_to_contig_end
        return distance is not None and distance < EDGE_WINDOW

    @property
    def best_type(self) -> Tuple[str, float]:
        """The most probable product class and its probability, regardless of GECCO's `type`.

        Mainly useful when `type` is ``Unknown``, to show which class came closest.
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

    With no clusters to report, GECCO writes `""` in these columns instead of leaving them blank,
    which a plain split would turn into a protein named by two quote characters.
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

    *contig_lengths* is needed to judge truncation, as in :func:`metawrap2.amr.parse_amrfinder`:
    a cluster starting at 392 bp may be near the end of a short contig or deep inside a long one.
    GECCO does not report contig lengths, so the caller supplies them from the assembly.

    A header-only file is a valid result, since most genomes have few clusters or none. `annotate`
    passes GECCO ``--force-tsv`` because otherwise a genome with no clusters writes no table, which
    looks the same as a failed run.
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
            # With --force-tsv, GECCO writes a placeholder row for a genome with no clusters: all
            # fields empty except `type` (Unknown) and `proteins`/`domains` (""). Skip it so it is
            # not counted as a BGC.
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
        # Written even when `type` is a real class, so the confidence of the call is visible
        # (a Saccharide call at p=0.23 differs from one at p=0.95).
        "most_probable_type": best_type,
        "most_probable_type_p": "%.4g" % best_p,
        "average_gene_probability": (
            "" if cluster.average_p is None else "%.4g" % cluster.average_p
        ),
        "max_gene_probability": "" if cluster.max_p is None else "%.4g" % cluster.max_p,
        # Blank instead of 0 when the contig length is unknown, so it is not read as "at the edge".
        "distance_to_contig_end": "" if distance is None else str(distance),
        "truncated": "yes" if cluster.is_truncated else ("" if distance is None else "no"),
        "domains": ";".join(cluster.domains),
        "proteins": ";".join(cluster.proteins),
    }


#: Printed with every BGC count. Clusters are defined by adjacency, so fragmentation affects them
#: most, and a cluster count in MAGs only covers what was assembled and binned.
TRUNCATION_CAVEAT = (
    "gene clusters are defined by adjacent genes, so assembly breaks affect them most. A cluster "
    "within %d bp of its contig's end may continue into unassembled sequence, and its gene "
    "content and product class describe only the assembled part. These rows are flagged "
    "`truncated`" % EDGE_WINDOW
)


@dataclass(frozen=True)
class GenomeBgcSummary:
    """One genome's cluster counts, split by truncation and by classification."""

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
