"""Distribute a genome's or sample's KOs over the KEGG BRITE hierarchy, as tables and Krona input.

BRITE ko00001 has four levels:

=========  =========================================  =======================================
level      example                                    what it is
=========  =========================================  =======================================
A          ``Metabolism``                             the broadest division, 8 of them
B          ``Carbohydrate metabolism``                the chemical process, 57 of them
C          ``Glycolysis / Gluconeogenesis``           the pathway, 570 of them
D          ``K00844`` hexokinase                      the gene, 28,545 of them
=========  =========================================  =======================================

Two weightings are supported:

* Per genome, every gene weighted equally: the share of a MAG's annotated genes in each category.
* Per sample, each genome weighted by its abundance: the share of the community's functional
  capacity. This is not the average of the per-genome numbers; rare genomes contribute little.

A KO can belong to several pathways (65k KO-pathway pairs over 28.5k KOs). :func:`distribute`
splits each KO's weight equally across its paths, so the output sums to the input total and
percentages add up to 100.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .io.seqio import smart_open

#: Where BRITE ko00001 is downloaded from by `metawrap2 install-db kegg`.
BRITE_URL = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir="

#: Files `install-db kegg` writes into the KEGG_BRITE directory: the projection this module reads
#: (one row per (KO, pathway) pair with its level A and B), and a record of the download.
HIERARCHY_FILE = "ko_hierarchy.tsv.gz"
SNAPSHOT_FILE = "SNAPSHOT.json"

#: The projection's columns, in order.
HIERARCHY_FILE_COLUMNS: Tuple[str, ...] = (
    "ko",
    "level_a_id",
    "level_a",
    "level_b_id",
    "level_b",
    "pathway_id",
    "pathway",
    "gene_symbol",
    "gene_name",
    "ec",
)

#: The levels, outermost first. Shared by the Krona input and the long and wide tables.
LEVELS: Tuple[str, ...] = ("level_a", "level_b", "pathway")

#: BRITE top-level branches that group proteins instead of describing pathways. ``09180 Brite
#: Hierarchies`` is the largest branch (22,768 of 64,956 rows), with entries such as "Protein
#: families: genetic information processing". ``09190 Not Included in Pathway or Brite`` holds the
#: rest.
#:
#: They are kept in every table, so percentages include genes KEGG does not place in a pathway, and
#: flagged so pathway-only views can exclude them.
NON_PATHWAY_BRANCHES: Tuple[str, ...] = ("09180", "09190")

#: KEGG's pathway branches, in BRITE order. ``Organismal Systems`` and ``Human Diseases`` are
#: annotated from host biology, so hits from a bacterial MAG usually reflect a KO shared with a host
#: pathway.
PATHWAY_BRANCHES: Tuple[str, ...] = ("09100", "09120", "09130", "09140", "09150", "09160")

#: Label for KOs not in the hierarchy. They are kept so percentages are not inflated.
UNCLASSIFIED = "Unclassified"


@dataclass(frozen=True)
class Hierarchy:
    """KO -> the (level_a, level_b, pathway) paths it belongs to, and the names of things."""

    paths: Mapping[str, Tuple[Tuple[str, str, str], ...]]
    #: KO -> gene symbol and name, for labelling the deepest level of a table.
    gene_names: Mapping[str, str]
    #: Level-A name -> BRITE id, so branches are tested by id instead of by display name.
    level_a_ids: Mapping[str, str]
    #: Download date and checksum, written into every output row, since the hierarchy is revised.
    snapshot: str = ""

    def is_pathway_branch(self, level_a: str) -> bool:
        """True when this top-level branch is one of PATHWAY_BRANCHES."""
        return self.level_a_ids.get(level_a, "") in PATHWAY_BRANCHES

    def __contains__(self, ko: str) -> bool:
        return ko in self.paths

    def paths_for(self, ko: str) -> Tuple[Tuple[str, str, str], ...]:
        """Every path this KO sits on, or one Unclassified path if KEGG does not place it."""
        return self.paths.get(ko, ((UNCLASSIFIED, UNCLASSIFIED, UNCLASSIFIED),))


def _require(hierarchy: Optional[Hierarchy]) -> Hierarchy:
    """Return *hierarchy*, or raise if None. No default is bundled; load one from KEGG_BRITE."""
    if hierarchy is None:
        raise ValueError(
            "no KEGG hierarchy given; load one with load_hierarchy(<KEGG_BRITE directory>)"
        )
    return hierarchy


def project_htext(htext: str, out_path: str) -> int:
    """Project BRITE's ko00001 htext into the per-(KO, pathway) table. Returns the rows written.

    The htext is indented by level: ``A09100 Metabolism``, ``B  09101 Carbohydrate metabolism``,
    ``C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]``, then one ``D`` line per KO,
    ``K00844  HK; hexokinase [EC:2.7.1.1]``. A KO appears once under every pathway it belongs to.
    """
    level_a = level_b = pathway = None
    rows = 0
    with open(htext, encoding="utf-8") as src, gzip.open(out_path, "wt", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(HIERARCHY_FILE_COLUMNS)
        for line in src:
            line = line.rstrip("\n")
            if not line or line[0] not in "ABCD" or not line[1:].strip():
                continue
            ident, _, name = line[1:].strip().partition(" ")
            name = name.strip()
            if line[0] == "A":
                level_a, level_b, pathway = (ident, name), None, None
            elif line[0] == "B":
                level_b, pathway = (ident, name), None
            elif line[0] == "C":
                pathway = (ident, re.sub(r"\s*\[[A-Z]+:[^\]]*\]$", "", name))
            elif level_a and level_b and pathway:
                ec = ""
                found = re.search(r"\s*\[EC:([^\]]*)\]$", name)
                if found:
                    ec, name = found.group(1), name[: found.start()]
                symbol, separator, gene_name = name.partition("; ")
                if not separator:
                    symbol, gene_name = "", name
                writer.writerow(
                    [ident, *level_a, *level_b, *pathway, symbol.strip(), gene_name.strip(), ec]
                )
                rows += 1
    return rows


def write_snapshot(directory: str, source: str, url: str, rows: int) -> None:
    """Record what was downloaded and when, next to the projection."""
    digest = hashlib.md5()
    with open(source, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    import datetime

    record = {
        "source": url,
        "downloaded": datetime.datetime.now(datetime.timezone.utc).date().isoformat(),
        "md5": digest.hexdigest(),
        "rows": rows,
    }
    with open(os.path.join(directory, SNAPSHOT_FILE), "w") as fh:
        json.dump(record, fh, indent=2)


def snapshot_label(directory: str) -> str:
    """``ko00001 2026-10-14 md5:1a2b3c4d``, or "unknown" when the download was not recorded."""
    try:
        with open(os.path.join(directory, SNAPSHOT_FILE)) as fh:
            record = json.load(fh)
        return "ko00001 %s md5:%s" % (record["downloaded"], record["md5"][:8])
    except (OSError, ValueError, KeyError):
        return "unknown"


def hierarchy_path(location: str) -> str:
    """The projection inside a KEGG_BRITE directory, or *location* itself if it is the file."""
    if os.path.isdir(location):
        return os.path.join(location, HIERARCHY_FILE)
    return location


def load_hierarchy(location: str) -> Hierarchy:
    """Read the projection `install-db kegg` wrote. *location* is the directory or the file."""
    path = hierarchy_path(location)
    paths: Dict[str, List[Tuple[str, str, str]]] = defaultdict(list)
    names: Dict[str, str] = {}
    level_a_ids: Dict[str, str] = {}
    with smart_open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ko = row["ko"]
            paths[ko].append((row["level_a"], row["level_b"], row["pathway"]))
            level_a_ids.setdefault(row["level_a"], row["level_a_id"])
            if ko not in names:
                symbol, description = row.get("gene_symbol", ""), row.get("gene_name", "")
                names[ko] = ("%s; %s" % (symbol, description)).strip("; ") or ko
    return Hierarchy(
        paths={ko: tuple(dict.fromkeys(ps)) for ko, ps in paths.items()},
        gene_names=names,
        level_a_ids=level_a_ids,
        snapshot=snapshot_label(os.path.dirname(os.path.abspath(path))),
    )


def distribute(
    weights: Mapping[str, float], hierarchy: Optional[Hierarchy] = None
) -> Dict[Tuple[str, str, str], float]:
    """Spread each KO's weight over the hierarchy paths it belongs to.

    *weights* is KO -> weight: a gene count for a per-genome distribution, or an abundance-weighted
    gene count for a per-sample one. The weight of a KO on *n* paths is split *n* ways, so the
    returned values sum to ``sum(weights.values())``.
    """
    tree = _require(hierarchy)
    out: Dict[Tuple[str, str, str], float] = defaultdict(float)
    for ko, weight in weights.items():
        if not weight:
            continue
        paths = tree.paths_for(ko)
        share = weight / float(len(paths))
        for path in paths:
            out[path] += share
    return dict(out)


def collapse(
    distribution: Mapping[Tuple[str, str, str], float], depth: int
) -> Dict[Tuple[str, ...], float]:
    """Sum a distribution up to *depth* levels (1 = level A, 2 = + level B, 3 = + pathway)."""
    if not 1 <= depth <= len(LEVELS):
        raise ValueError("depth must be 1..%d, not %r" % (len(LEVELS), depth))
    out: Dict[Tuple[str, ...], float] = defaultdict(float)
    for path, weight in distribution.items():
        out[tuple(path[:depth])] += weight
    return dict(out)


#: Columns of the long hierarchy table.
HIERARCHY_COLUMNS: Tuple[str, ...] = (
    "unit",
    "name",
    "level",
    "level_a",
    "level_b",
    "pathway",
    "weight",
    "percent",
    "percent_of_pathway_branches",
    "is_pathway_branch",
    # Which abundance the weights came from. MAG-abundance weighting describes the binned part of
    # the community; gene-abundance weighting (reads mapped to a gene catalogue) describes the whole
    # sample. Recorded so outputs from different runs can be told apart.
    "weighting",
    "hierarchy_snapshot",
)

#: Values of `weighting`. One genome's genes counted equally (the per-genome unit).
WEIGHTING_GENES = "gene_count"
#: MAG abundance x per-genome gene counts: the binned community's potential.
WEIGHTING_MAG_ABUNDANCE = "mag_abundance"
#: Per-gene abundance from mapping reads to a gene catalogue: the whole sample's potential.
WEIGHTING_GENE_ABUNDANCE = "gene_abundance"
#: Per-contig depth, each gene weighted by its contig's depth. Like gene abundance it covers the
#: whole sample, but it reuses the depth table `binning` wrote instead of mapping reads again. Less
#: precise per gene.
WEIGHTING_CONTIG_DEPTH = "contig_depth"


def hierarchy_rows(
    unit: str,
    name: str,
    distribution: Mapping[Tuple[str, str, str], float],
    hierarchy: Optional[Hierarchy] = None,
    weighting: str = WEIGHTING_GENES,
) -> List[Dict[str, object]]:
    """Long rows for one genome or one sample, at every depth.

    All three depths are emitted, each summing to 100%, so callers need not re-aggregate.

    Each row has two percentages:

    * ``percent``: of everything annotated, including NON_PATHWAY_BRANCHES. Sums to 100 per level.
    * ``percent_of_pathway_branches``: of PATHWAY_BRANCHES only. Empty, not zero, for rows outside
      them.
    """
    tree = _require(hierarchy)
    total = sum(distribution.values())
    pathway_total = sum(
        weight for path, weight in distribution.items() if tree.is_pathway_branch(path[0])
    )
    rows: List[Dict[str, object]] = []
    for depth in range(1, len(LEVELS) + 1):
        for path, weight in sorted(
            collapse(distribution, depth).items(), key=lambda kv: (-kv[1], kv[0])
        ):
            padded = list(path) + [""] * (len(LEVELS) - len(path))
            in_pathways = tree.is_pathway_branch(padded[0])
            rows.append(
                {
                    "unit": unit,
                    "name": name,
                    "level": LEVELS[depth - 1],
                    "level_a": padded[0],
                    "level_b": padded[1],
                    "pathway": padded[2],
                    "weight": "%.4f" % weight,
                    "percent": "%.4f" % (100.0 * weight / total) if total else "0.0000",
                    "percent_of_pathway_branches": (
                        "%.4f" % (100.0 * weight / pathway_total)
                        if in_pathways and pathway_total
                        else ""
                    ),
                    "is_pathway_branch": int(in_pathways),
                    "weighting": weighting,
                    "hierarchy_snapshot": tree.snapshot,
                }
            )
    return rows


def write_krona_input(
    distribution: Mapping[Tuple[str, str, str], float],
    path: str,
    min_weight: float = 0.0,
    hierarchy: Optional[Hierarchy] = None,
    pathways_only: bool = False,
) -> int:
    """Write the TSV ``ktImportText`` takes: a weight, then one column per level.

    Krona is already used by ``kraken2``. Each line is ``<weight>\\tlevel_a\\tlevel_b\\tpathway``;
    Krona computes the nesting and percentages.

    *pathways_only* drops NON_PATHWAY_BRANCHES. The full plot shows how much KEGG cannot place in a
    pathway, but that branch is the largest and dominates the chart.
    """
    tree = _require(hierarchy) if pathways_only else hierarchy
    written = 0
    with open(path, "w") as fh:
        for path_tuple, weight in sorted(distribution.items(), key=lambda kv: -kv[1]):
            if weight <= min_weight:
                continue
            if pathways_only and tree is not None and not tree.is_pathway_branch(path_tuple[0]):
                continue
            fh.write("%.6f\t%s\n" % (weight, "\t".join(path_tuple)))
            written += 1
    return written


def weights_from_gene_kos(gene_kos: Iterable[Sequence[str]]) -> Dict[str, float]:
    """KO -> gene count, from an iterable of each gene's KO assignments.

    A gene with several KOs contributes an equal fraction to each, so the total is the gene count.
    """
    out: Dict[str, float] = defaultdict(float)
    for kos in gene_kos:
        unique = list(dict.fromkeys(kos))
        if not unique:
            continue
        share = 1.0 / len(unique)
        for ko in unique:
            out[ko] += share
    return dict(out)


def weights_from_gene_abundance(
    gene_kos: Mapping[str, Sequence[str]],
    abundance: Mapping[str, float],
) -> Dict[str, float]:
    """KO -> weight for one sample, from *per-gene* abundance rather than per-genome abundance.

    ``A(ko) = sum over genes of abundance(gene) * 1/|KOs(gene)|``: the share-splitting of
    :func:`weights_from_gene_kos`, weighted by each gene's abundance.

    This is the whole-sample functional potential. A gene catalogue from the assembly includes the
    ~90% of assembled sequence that never reaches a MAG, so pathways on unbinned contigs appear here
    but not in the MAG-weighted profile.

    Genes with no abundance entry contribute nothing.
    """
    out: Dict[str, float] = defaultdict(float)
    for gene, kos in gene_kos.items():
        value = abundance.get(gene)
        if not value:
            continue
        unique = list(dict.fromkeys(kos))
        if not unique:
            continue
        share = value / len(unique)
        for ko in unique:
            out[ko] += share
    return dict(out)


def weights_from_contig_depth(
    gene_kos: Mapping[str, Sequence[str]],
    depth: Mapping[str, float],
    contig_of: Optional[Mapping[str, str]] = None,
) -> Dict[str, float]:
    """KO -> weight for one sample, from *per-contig* depth.

    Uses the per-contig depth table `binning` writes, so no extra read mapping is needed. Every gene
    on a contig gets that contig's depth, which assumes uniform depth along the contig.

    *contig_of* maps a gene id to its contig. Without it the contig is taken from the gene id
    (``<contig>_<n>``, as :mod:`metawrap2.genecall` writes it); see :func:`contig_from_gene`.

    Genes whose contig has no depth entry contribute nothing.
    """
    out: Dict[str, float] = defaultdict(float)
    for gene, kos in gene_kos.items():
        contig = (contig_of or {}).get(gene) or contig_from_gene(gene)
        value = depth.get(contig)
        if not value:
            continue
        unique = list(dict.fromkeys(kos))
        if not unique:
            continue
        share = value / len(unique)
        for ko in unique:
            out[ko] += share
    return dict(out)


def contig_from_gene(gene: str) -> str:
    """The contig a gene id belongs to, for ids of the form ``<contig>_<n>``.

    Splits at the last underscore, because contig names contain underscores:
    ``k141_73537_length_916682_cov_50.9988_12`` is gene 12 of contig
    ``k141_73537_length_916682_cov_50.9988``.

    A catalogue gene's source prefix (``bin_003~k141_..._12``) is kept, because depth tables built
    over the same catalogue carry it too.
    """
    head, _, tail = gene.rpartition("_")
    return head if head and tail.isdigit() else gene


def weighted_sample_profile(
    per_genome: Mapping[str, Mapping[str, float]],
    abundance: Mapping[str, float],
) -> Dict[str, float]:
    """Community KO weights: each genome's KO weights scaled by that genome's abundance.

    ``A(ko, sample) = sum over genomes of n(ko, genome) * abundance(genome, sample)``. *abundance*
    must be coverage-proportional (CoverM relative_abundance or trimmed_mean) and not divided by
    genome size.

    Genomes with no abundance entry contribute nothing.
    """
    out: Dict[str, float] = defaultdict(float)
    for genome, kos in per_genome.items():
        scale = abundance.get(genome, 0.0)
        if not scale:
            continue
        for ko, weight in kos.items():
            out[ko] += weight * scale
    return dict(out)
