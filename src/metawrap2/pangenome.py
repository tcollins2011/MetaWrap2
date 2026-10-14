"""Pangenomes over a cohort's MAGs: which genes are in every strain of a species and which are not.

The unit of a pangenome is **one species**, and that is the thing most easily got wrong here. A
dereplicated catalogue holds one genome per species, so a pangenome computed over a catalogue is a
pangenome of a community - a quantity with no definition. The input that makes sense is a
*dereplication cluster*: the same organism recovered from several samples, which is exactly what
`dereplicate` groups and then discards everything but a representative of.

So this module reads `dereplicate`'s cluster table and builds one pangenome per cluster that has
enough members, and refuses to treat the catalogue as one group. :func:`groups_from_clusters` is
that decision, and it is where the honest part lives.

What a pangenome gives you that a catalogue does not: PPanGGOLiN partitions genes into **persistent**
(in essentially every genome of the species), **shell** (in some) and **cloud** (in one or a few).
The shell is where the strain-level differences live - the gene content that distinguishes the
strain in one sample from the strain in another - and it is invisible to anything computed on a
single representative.

The floor is the part to state plainly. PPanGGOLiN's partitioning is a statistical model over
presence/absence across genomes, and with a handful of genomes there is not enough signal to fit it:
its own documentation recommends at least 15, and below that the partitions are better read as "core
versus accessory in these few genomes" than as a species pangenome. :data:`RELIABLE_GENOMES` is that
number, and groups below it are built if asked but reported as provisional.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

#: Genomes below which PPanGGOLiN's persistent/shell/cloud partitioning is not worth reading as a
#: species pangenome. Its own guidance. Not a refusal: a 5-genome group still answers "which genes
#: are in all five", which is a smaller but real question.
RELIABLE_GENOMES = 15

#: Genomes below which there is no pangenome at all - one genome has no gene-content variation and
#: two have a comparison rather than a partitioning.
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
    #: The cluster's representative, which is the genome the catalogue kept.
    representative: Optional[str] = None

    @property
    def size(self) -> int:
        return len(self.genomes)

    @property
    def is_reliable(self) -> bool:
        """Whether a persistent/shell/cloud partitioning over this group is worth reading."""
        return self.size >= RELIABLE_GENOMES

    @property
    def is_a_pangenome(self) -> bool:
        return self.size >= MIN_GENOMES


def read_clusters(path: str) -> Dict[str, List[Tuple[str, bool]]]:
    """cluster -> [(genome, is_representative)], from `dereplicate`'s clusters.tsv."""
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
    """Turn dereplication clusters into pangenome groups, with a census of what was skipped.

    The census is the output that matters as much as the groups: on a three-sample cohort almost
    every cluster is a singleton, so "0 pangenomes" needs to come with "23 clusters had one genome"
    rather than looking like a failure.

    Genomes whose FASTA is not in *genome_paths* are dropped from their cluster and counted, because
    a cluster table can name genomes from a bin set that is no longer on disk.
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
    """PPanGGOLiN's ``--fasta`` input: genome name then path, tab separated, one line each."""
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for genome in sorted(group.genomes):
            writer.writerow([genome, group.genomes[genome]])
    return path


#: Where PPanGGOLiN writes one file of gene-family names per partition. This is the group-level
#: answer and the per-genome statistics table is not: that table's `Persistent_families` column is
#: how many persistent families *this genome* carries, so summing it over a group multiplies every
#: shared family by the number of genomes that have it.
PARTITION_DIR = "partitions"
PARTITION_FILES = {PERSISTENT: "persistent.txt", SHELL: "shell.txt", CLOUD: "cloud.txt"}

#: PPanGGOLiN's per-genome statistics, used only for genes per genome. Its first two lines are
#: `#`-prefixed settings rather than a header, so a plain DictReader takes `#soft_core=0.95` as the
#: column names and every lookup returns None.
STATS_GENOME = "Genome_name"
STATS_GENES = "Genes"


@dataclass(frozen=True)
class PangenomeSummary:
    """One group's gene-content partitioning."""

    group: str
    genomes: int
    #: Gene families per partition, over the group.
    families: Mapping[str, int]
    #: Mean genes per genome, which is not the same as families: a family can be present twice.
    mean_genes: Optional[float]
    reliable: bool

    @property
    def total_families(self) -> int:
        return sum(self.families.values())

    @property
    def shares(self) -> Dict[str, float]:
        """Each partition as a percentage of the families, or an empty dict when there are none."""
        total = self.total_families
        if not total:
            return {}
        return {name: 100.0 * count / total for name, count in self.families.items()}

    @property
    def openness(self) -> Optional[float]:
        """Share of families outside the persistent partition, as a percentage.

        The one-number summary of how varied the species' gene content is across these genomes.
        Deliberately not called a "pangenome openness index": that term means a fitted Heaps-law
        exponent, which needs a rarefaction curve this does not compute.
        """
        shares = self.shares
        if not shares:
            return None
        return 100.0 - shares.get(PERSISTENT, 0.0)


def read_genome_statistics(path: str) -> Dict[str, float]:
    """genome -> gene count, from PPanGGOLiN's per-genome statistics table.

    Comment lines are skipped explicitly. The file opens with two `#`-prefixed settings lines, so
    csv.DictReader handed the file directly keys every row by `#soft_core=0.95` and returns None for
    every column - which reads as a pangenome with no genes rather than as a parse failure.
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

    Families are counted from `partitions/*.txt`, one name per line, which is the group-level
    quantity. Genes per genome come from the statistics table, which is the only thing it is used
    for here.
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
        # Stated per row rather than once in a log line: a provisional partitioning read out of a
        # table months later is indistinguishable from a sound one otherwise.
        "partitioning_reliable": "yes" if summary.reliable else "no",
        "minimum_genomes_for_reliable_partitioning": str(RELIABLE_GENOMES),
    }


#: A persistent share below this is evidence about the *cluster* rather than about the species.
#: Genomes of one species share most of their gene content; three genomes sharing a fifth of their
#: families are not three strains of one organism, they are a dereplication cluster that should be
#: looked at. Measured: of seven real clusters from a three-sample cohort, six sat at 68-96%
#: persistent and one at 20.7% with 7,361 families over a mean of 3,214 genes per genome.
SUSPECT_PERSISTENT = 50.0

#: Printed beside a group whose persistent share is below SUSPECT_PERSISTENT.
SUSPECT_CLUSTER_CAVEAT = (
    "a persistent share this low says the genomes in the group are not strains of one organism. "
    "The pangenome is not wrong - the grouping is the thing to check, starting with the cluster's "
    "ANI and the contamination of its members"
)


def suspect_groups(summaries: Sequence[PangenomeSummary]) -> List[PangenomeSummary]:
    """Groups whose gene content says they are not one species."""
    return [
        summary
        for summary in summaries
        if (summary.shares.get(PERSISTENT, 100.0) < SUSPECT_PERSISTENT)
    ]


#: Printed with every pangenome built from fewer than RELIABLE_GENOMES genomes.
SMALL_GROUP_CAVEAT = (
    "a persistent/shell/cloud partitioning is a statistical model over presence and absence across "
    "genomes, and below %d genomes there is not enough signal to fit it - PPanGGOLiN's own "
    "guidance. These rows are better read as 'core versus accessory among these few genomes' than "
    "as a species pangenome, and they are flagged `partitioning_reliable=no`" % RELIABLE_GENOMES
)

#: Printed when a cohort yields no group at all, which is the normal result for a few samples.
NO_GROUPS_CAVEAT = (
    "no cluster had enough genomes for a pangenome. That is the expected result for a handful of "
    "samples: a pangenome needs the same species recovered from several of them, and a "
    "dereplicated catalogue keeps one genome per species precisely because the rest were "
    "redundant. It is a statement about the cohort's size, not about the species"
)


def census_lines(census: Mapping[str, int]) -> List[str]:
    """The skipped-cluster census as lines to print, so "0 pangenomes" is never bare."""
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
            "%d genome(s) named in the cluster table have no FASTA under the bin folder given, and "
            "were left out of their cluster" % census["genomes_without_a_file"]
        )
    return lines


def genome_paths(folders: Iterable[str]) -> Dict[str, str]:
    """genome label -> path, labelled exactly as `dereplicate` labels its cluster rows.

    Through dereplicate's own `collect_genomes` rather than by stripping extensions here, because
    the label is not just the filename: it is ``<source-directory>__<bin>``, with a numeric suffix
    when two bin directories share a basename - which is the normal cohort layout, where every
    sample has its own ``metawrap_70_10_bins``. Reimplementing that would produce labels that look
    right and join to nothing, and the cluster table would appear to name genomes that are not
    there.
    """
    # Imported here rather than at module scope: this is a library and dereplicate is a module, so
    # a top-level import would have the layering backwards.
    from .modules.dereplicate import collect_genomes

    return {genome.name: genome.path for genome in collect_genomes(list(folders))}
