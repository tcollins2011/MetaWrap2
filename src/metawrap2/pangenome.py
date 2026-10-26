"""Per-species pangenomes over a cohort's MAGs.

A pangenome is defined for one species. A dereplicated catalogue holds one genome per species, so
the input here is a dereplication cluster: the same organism recovered from several samples.
This module reads `dereplicate`'s cluster table and builds one pangenome per cluster with enough
members (:func:`groups_from_clusters`); it never treats the whole catalogue as one group.

PPanGGOLiN partitions gene families into persistent (in nearly every genome), shell (in some) and
cloud (in one or a few). The shell holds the strain-level differences between samples, which a
single representative does not show.

The partitioning is a statistical model over presence/absence. PPanGGOLiN's documentation
recommends at least 15 genomes (:data:`RELIABLE_GENOMES`); smaller groups are built but reported
as provisional, and read as core versus accessory among those genomes.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

#: Minimum genomes for a reliable persistent/shell/cloud partitioning, from PPanGGOLiN's guidance.
#: Smaller groups are still built and flagged.
RELIABLE_GENOMES = 15

#: Minimum genomes for a pangenome to be built at all.
MIN_GENOMES = 3

#: PPanGGOLiN's partition names, as it writes them.
PERSISTENT = "persistent"
SHELL = "shell"
CLOUD = "cloud"
PARTITIONS: Tuple[str, ...] = (PERSISTENT, SHELL, CLOUD)


@dataclass(frozen=True)
class GenomeGroup:
    """A set of genomes of one species, from one dereplication cluster."""

    name: str
    #: genome label -> path to its FASTA.
    genomes: Mapping[str, str]
    #: The cluster's representative, the genome kept in the catalogue.
    representative: Optional[str] = None

    @property
    def size(self) -> int:
        return len(self.genomes)

    @property
    def is_reliable(self) -> bool:
        """Whether the group has enough genomes for a reliable partitioning."""
        return self.size >= RELIABLE_GENOMES

    @property
    def is_a_pangenome(self) -> bool:
        return self.size >= MIN_GENOMES


def read_clusters(path: str) -> Dict[str, List[Tuple[str, bool]]]:
    """Return cluster -> [(genome, is_representative)] from `dereplicate`'s clusters.tsv."""
    out: Dict[str, List[Tuple[str, bool]]] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and "cluster" not in reader.fieldnames:
            raise KeyError(
                "%s does not look like dereplicate's clusters.tsv: no 'cluster' column" % path
            )
        for row in reader:
            cluster = (row.get("cluster") or "").strip()
            genome = (row.get("genome") or "").strip()
            if not cluster or not genome:
                continue
            representative = (row.get("representative") or "").strip().lower() in (
                "yes",
                "true",
                "1",
            )
            out.setdefault(cluster, []).append((genome, representative))
    return out


def groups_from_clusters(
    clusters: Mapping[str, Sequence[Tuple[str, bool]]],
    genome_paths: Mapping[str, str],
    min_genomes: int = MIN_GENOMES,
) -> Tuple[List[GenomeGroup], Dict[str, int]]:
    """Turn dereplication clusters into pangenome groups, with counts of what was skipped.

    In a small cohort most clusters are singletons, so the counts explain an empty result.
    Genomes whose FASTA is not in *genome_paths* (for example, a bin set no longer on disk) are
    dropped from their cluster and counted.
    """
    groups: List[GenomeGroup] = []
    census = {
        "clusters": len(clusters),
        "singletons": 0,
        "too_small": 0,
        "pangenomes": 0,
        "provisional": 0,
        "genomes_without_a_file": 0,
    }
    for cluster, members in sorted(clusters.items()):
        found: Dict[str, str] = {}
        representative: Optional[str] = None
        for genome, is_representative in members:
            path = genome_paths.get(genome)
            if path is None:
                census["genomes_without_a_file"] += 1
                continue
            found[genome] = path
            if is_representative:
                representative = genome
        if len(found) <= 1:
            census["singletons"] += 1
            continue
        if len(found) < min_genomes:
            census["too_small"] += 1
            continue
        group = GenomeGroup(name=cluster, genomes=found, representative=representative)
        groups.append(group)
        census["pangenomes"] += 1
        if not group.is_reliable:
            census["provisional"] += 1
    return groups, census


def write_genome_list(group: GenomeGroup, path: str) -> str:
    """Write PPanGGOLiN's ``--fasta`` input: genome name and path, tab separated, one per line."""
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for genome in sorted(group.genomes):
            writer.writerow([genome, group.genomes[genome]])
    return path


#: Directory where PPanGGOLiN writes one file of gene-family names per partition; family counts
#: come from here. The per-genome statistics table counts families per genome, so summing it over
#: a group counts each shared family once per genome.
PARTITION_DIR = "partitions"
PARTITION_FILES = {PERSISTENT: "persistent.txt", SHELL: "shell.txt", CLOUD: "cloud.txt"}

#: Columns of PPanGGOLiN's per-genome statistics, used only for genes per genome.
STATS_GENOME = "Genome_name"
STATS_GENES = "Genes"


@dataclass(frozen=True)
class PangenomeSummary:
    """One group's gene-content partitioning."""

    group: str
    genomes: int
    #: Gene families per partition, over the group.
    families: Mapping[str, int]
    #: Mean genes per genome. Differs from families, since a family can occur twice in a genome.
    mean_genes: Optional[float]
    reliable: bool

    @property
    def total_families(self) -> int:
        return sum(self.families.values())

    @property
    def shares(self) -> Dict[str, float]:
        """Return each partition as a percentage of the families; empty when there are none."""
        total = self.total_families
        if not total:
            return {}
        return {name: 100.0 * count / total for name, count in self.families.items()}

    @property
    def openness(self) -> Optional[float]:
        """Share of families outside the persistent partition, as a percentage.

        A summary of how much gene content varies across the genomes. This is not the pangenome
        openness index, which is a fitted Heaps-law exponent.
        """
        shares = self.shares
        if not shares:
            return None
        return 100.0 - shares.get(PERSISTENT, 0.0)


def read_genome_statistics(path: str) -> Dict[str, float]:
    """Return genome -> gene count from PPanGGOLiN's per-genome statistics table.

    The file starts with two `#`-prefixed settings lines, which are skipped so that csv.DictReader
    finds the real header.
    """
    out: Dict[str, float] = {}
    if not os.path.isfile(path):
        return out
    with open(path) as handle:
        lines = [line for line in handle if not line.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        genome = (row.get(STATS_GENOME) or "").strip()
        if not genome:
            continue
        try:
            out[genome] = float(row.get(STATS_GENES) or 0)
        except ValueError:
            continue
    return out


def parse_partitions(outdir: str, group: str, reliable: bool) -> Optional[PangenomeSummary]:
    """Read one group's partitioning from PPanGGOLiN's output directory.

    Families are counted from `partitions/*.txt` (one name per line). Genes per genome come from
    the statistics table.
    """
    counts: Dict[str, int] = {}
    for partition, filename in PARTITION_FILES.items():
        path = os.path.join(outdir, PARTITION_DIR, filename)
        if not os.path.isfile(path):
            continue
        with open(path) as handle:
            counts[partition] = sum(1 for line in handle if line.strip())
    if not counts:
        return None
    genes = read_genome_statistics(os.path.join(outdir, "genomes_statistics.tsv"))
    return PangenomeSummary(
        group=group,
        genomes=len(genes),
        families=counts,
        mean_genes=(sum(genes.values()) / len(genes)) if genes else None,
        reliable=reliable,
    )


#: Columns of the per-group table.
SUMMARY_COLUMNS: Tuple[str, ...] = (
    "group",
    "genomes",
    "mean_genes_per_genome",
    "persistent_genes",
    "shell_genes",
    "cloud_genes",
    "percent_persistent",
    "percent_shell",
    "percent_cloud",
    "percent_outside_persistent",
    "partitioning_reliable",
    "minimum_genomes_for_reliable_partitioning",
)


def summary_row(summary: PangenomeSummary) -> Dict[str, str]:
    shares = summary.shares

    def num(value: Optional[float], places: int = 2) -> str:
        return "" if value is None else "%.*f" % (places, value)

    return {
        "group": summary.group,
        "genomes": str(summary.genomes),
        "mean_genes_per_genome": num(summary.mean_genes, 1),
        "persistent_genes": str(summary.families.get(PERSISTENT, 0)),
        "shell_genes": str(summary.families.get(SHELL, 0)),
        "cloud_genes": str(summary.families.get(CLOUD, 0)),
        "percent_persistent": num(shares.get(PERSISTENT)),
        "percent_shell": num(shares.get(SHELL)),
        "percent_cloud": num(shares.get(CLOUD)),
        "percent_outside_persistent": num(summary.openness),
        # Stored per row so the table alone shows which partitionings are provisional.
        "partitioning_reliable": "yes" if summary.reliable else "no",
        "minimum_genomes_for_reliable_partitioning": str(RELIABLE_GENOMES),
    }


#: Persistent share (percent) below which a cluster is flagged as probably not one species. In one
#: three-sample cohort, six of seven clusters were 68-96% persistent and one was 20.7% (7,361
#: families over a mean of 3,214 genes per genome).
SUSPECT_PERSISTENT = 50.0

#: Printed beside a group whose persistent share is below SUSPECT_PERSISTENT.
SUSPECT_CLUSTER_CAVEAT = (
    "a persistent share this low suggests the genomes in the group are not strains of one "
    "organism. Check the grouping, starting with the cluster's "
    "ANI and the contamination of its members"
)


def suspect_groups(summaries: Sequence[PangenomeSummary]) -> List[PangenomeSummary]:
    """Return groups whose persistent share suggests they are not one species."""
    return [
        summary
        for summary in summaries
        if (summary.shares.get(PERSISTENT, 100.0) < SUSPECT_PERSISTENT)
    ]


#: Printed with every pangenome built from fewer than RELIABLE_GENOMES genomes.
SMALL_GROUP_CAVEAT = (
    "the persistent/shell/cloud partitioning is a statistical model of presence and absence "
    "across genomes; PPanGGOLiN's guidance is that below %d genomes there is too little signal "
    "to fit it. Read these rows as 'core versus accessory among these genomes', not "
    "as a species pangenome; they are flagged `partitioning_reliable=no`" % RELIABLE_GENOMES
)

#: Printed when a cohort yields no group, which is usual for a few samples.
NO_GROUPS_CAVEAT = (
    "no cluster had enough genomes for a pangenome. This is expected with few "
    "samples: a pangenome needs the same species recovered from several of them, and a "
    "dereplicated catalogue keeps one genome per species. More samples are needed"
)


def census_lines(census: Mapping[str, int]) -> List[str]:
    """Return the cluster counts as lines to print."""
    lines = [
        "%d cluster(s): %d pangenome(s) built, %d singleton(s), %d below the group floor"
        % (
            census.get("clusters", 0),
            census.get("pangenomes", 0),
            census.get("singletons", 0),
            census.get("too_small", 0),
        )
    ]
    if census.get("provisional"):
        lines.append(
            "%d of them have fewer than %d genomes, so their partitioning is provisional"
            % (census["provisional"], RELIABLE_GENOMES)
        )
    if census.get("genomes_without_a_file"):
        lines.append(
            "%d genome(s) named in the cluster table have no FASTA in the given bin folder and "
            "were left out of their cluster" % census["genomes_without_a_file"]
        )
    return lines


def genome_paths(folders: Iterable[str]) -> Dict[str, str]:
    """Return genome label -> path, with labels as `dereplicate` writes them in its cluster table.

    Uses dereplicate's `collect_genomes`, since labels are ``<source-directory>__<bin>`` with a
    numeric suffix when bin directories share a basename (every sample has its own
    ``metawrap_70_10_bins``).
    """
    # Imported here: this is a library and dereplicate is a module, which should not be imported
    # at the top level.
    from .modules.dereplicate import collect_genomes

    return {genome.name: genome.path for genome in collect_genomes(list(folders))}
