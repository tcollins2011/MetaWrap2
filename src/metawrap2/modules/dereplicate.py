"""MetaWrap2 dereplicate module: collapse many samples' bins into one non-redundant catalogue.

Binning each sample separately produces the same genome once per sample it occurs in. Twenty
samples of one gut cohort can yield two thousand bins for a few hundred organisms. Near-identical
genomes then split reads in abundance estimation, annotation repeats per copy, and the MAG count
overstates the number of organisms.

The bins are clustered by whole-genome identity and one representative is kept per cluster.
Identity comes from skani (Shaw & Yu, *Nat Methods* 2023), which estimates ANI from sparse k-mer
sketches and reports the aligned fraction. Clustering and representative choice are done here.

The result depends on two flags, SPECIES_ANI and MIN_ALIGNED_FRACTION. Neither has been tuned on
this project's data yet.

This file is meant to be read and edited: the command skani runs lives in the COMMANDS block
below, and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List, NamedTuple, Sequence, Set

from ..cliformat import make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, STATS_SUFFIX, TSV_SUFFIX
from ..utils import bin_name_from_filename, quality_score
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    list_bins,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["dereplicate"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
#: ``-E`` writes one row per pair instead of an N^2 matrix. ``--min-af`` drops pairs that share only
#: a short region (e.g. a plasmid) before clustering. ``-s`` is skani's pre-filter; pairs below it
#: are not compared.
SKANI_TRIANGLE = (
    "skani triangle -E -l {genome_list} -o {output} -t {threads} "
    "-s {screen} --min-af {min_af}{preset}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: ANI at which two bins are called the same species. 95% is the convention used by the large MAG
#: catalogues, so species counts are comparable with published ones. It was not tuned here; 99
#: separates strains, as an inStrain-style analysis needs.
SPECIES_ANI = 95.0
#: Minimum percentage of a genome that must align for its ANI estimate to count. Permissive: the
#: ANI does the species call, and this only stops a short shared region (a prophage, a plasmid, a
#: conserved operon) from standing in for the genome. Raising it gives smaller clusters and a larger
#: catalogue.
MIN_ALIGNED_FRACTION = 10.0
#: skani's pre-filter, in approximate percent identity. Pairs screened out cannot cluster, so it
#: must sit well below SPECIES_ANI. skani's default of 80 is used.
SKANI_SCREEN = 80
#: skani presets trade compression factor for accuracy on fragmented assemblies. skani's help
#: recommends ``medium`` for an N50 below 10 kb, at 2x the cost. Which one MAGs need has not been
#: measured, so skani's default is used.
SKANI_PRESET = "default"
SKANI_PRESETS = ("default", "fast", "medium", "slow")
#: Where the kept genomes and the record of what was collapsed into them are written.
CATALOGUE_DIR = "dereplicated_genomes"
CLUSTERS_TABLE = "clusters" + TSV_SUFFIX
ANI_TABLE = "ani" + TSV_SUFFIX
#: The cluster table's columns. ``ani_to_representative`` shows over-merged clusters.
CLUSTER_COLUMNS = (
    "cluster",
    "genome",
    "source",
    "representative",
    "ani_to_representative",
    "quality_score",
    "size_bp",
)
# ───────────────────────────────────────────────────────────────────────────────────────────


class Edge(NamedTuple):
    """One skani comparison: two genomes, their ANI, and how much of each aligned."""

    ref: str
    query: str
    ani: float
    af_ref: float
    af_query: float

    @property
    def aligned_fraction(self) -> float:
        """The larger of the two aligned fractions.

        A 60% complete bin aligns almost entirely to a complete bin of the same organism, which
        aligns only 60% to it. Using the smaller fraction would keep such duplicates apart.
        """
        return max(self.af_ref, self.af_query)


class Genome(NamedTuple):
    """One input bin: where it came from, and what it is called in the catalogue."""

    name: str  # unique within the run, "<source>__<bin>"
    path: str
    source: str  # the bin directory it came from, as a label


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 dereplicate",
        usage="metawrap2 dereplicate [options] -o output_dir -b bins_1 [-b bins_2 ...]",
        description="Collapse several samples' bins into one non-redundant genome catalogue.",
    )
    p.add_argument(
        "-b",
        "--bins",
        action="append",
        default=[],
        metavar="DIR",
        required=True,
        help="folder of bins to include; repeat once per sample",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "--ani",
        type=percent(),
        default=SPECIES_ANI,
        help="ANI%% at which two bins count as the same genome (the species boundary); "
        "99 gives strain-level clusters",
    )
    p.add_argument(
        "--min-aligned-fraction",
        type=percent(),
        default=MIN_ALIGNED_FRACTION,
        help="minimum %% of either genome that must align for an ANI value to be used",
    )
    p.add_argument(
        "--preset",
        choices=SKANI_PRESETS,
        default=SKANI_PRESET,
        help="skani speed/accuracy preset (default %s)" % SKANI_PRESET,
    )
    p.add_argument(
        "--copy",
        action="store_true",
        help="copy the representative genomes into the catalogue instead of symlinking them",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


# --- the parts that are pure, and therefore testable without skani ------------------------


def collect_genomes(bin_dirs: Sequence[str]) -> List[Genome]:
    """Every bin in every input directory, under a name unique across directories.

    Bin names repeat across samples (``bin_001.fasta`` exists in each), so the source directory's
    name is used as a prefix. The catalogue's filenames then also show where each genome came from.
    """
    genomes: List[Genome] = []
    seen: Set[str] = set()
    for directory in bin_dirs:
        source = os.path.basename(os.path.normpath(directory))
        for filename in list_bins(directory, FASTA_EXTENSIONS):
            name = "%s__%s" % (source, bin_name_from_filename(filename))
            # Input directories can share a basename (one per sample, all called
            # "metawrap_70_10_bins"); add a numeric suffix instead of refusing.
            if name in seen:
                suffix = 2
                while "%s_%d" % (name, suffix) in seen:
                    suffix += 1
                name = "%s_%d" % (name, suffix)
            seen.add(name)
            genomes.append(Genome(name, os.path.join(directory, filename), source))
    return genomes


def parse_skani(path: str) -> List[Edge]:
    """Read skani's sparse table into edges, keyed by file path.

    Columns are looked up by header name because skani's ``--ci`` and ``--detailed`` flags insert
    extra columns.
    """
    edges: List[Edge] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            columns = {
                field: header.index(field)
                for field in (
                    "Ref_file",
                    "Query_file",
                    "ANI",
                    "Align_fraction_ref",
                    "Align_fraction_query",
                )
            }
        except ValueError as exc:
            raise ValueError("%s is not a skani table: %s" % (path, exc)) from exc
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if len(cut) <= max(columns.values()):
                continue
            edges.append(
                Edge(
                    ref=cut[columns["Ref_file"]],
                    query=cut[columns["Query_file"]],
                    ani=float(cut[columns["ANI"]]),
                    af_ref=float(cut[columns["Align_fraction_ref"]]),
                    af_query=float(cut[columns["Align_fraction_query"]]),
                )
            )
    return edges


def cluster_genomes(
    names: Sequence[str],
    edges: Iterable[Edge],
    *,
    ani: float = SPECIES_ANI,
    min_aligned_fraction: float = MIN_ALIGNED_FRACTION,
) -> List[List[str]]:
    """Group *names* so that every pair above *ani* ends up together. Single linkage.

    Identity is not transitive, so single linkage can chain genomes into one cluster whose ends are
    below the threshold. Requiring every pair to pass would need the full matrix and make the result
    depend on input order. Single linkage never splits a true duplicate, and chaining is visible in
    the cluster table as a member with low ANI to its representative.

    Genomes with no edge form their own cluster, so every name is accounted for. Clusters and their
    members are returned in input order, so the result is reproducible.
    """
    index = {name: i for i, name in enumerate(names)}
    parent = list(range(len(names)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for edge in edges:
        if edge.ref == edge.query:  # skani's --diagonal, if anyone turns it on
            continue
        if edge.ani < ani or edge.aligned_fraction < min_aligned_fraction:
            continue
        a, b = index.get(edge.ref), index.get(edge.query)
        if a is None or b is None:
            continue
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            # Lower index wins, so a cluster's root is its earliest member and the order below
            # follows the input, not the edge order.
            parent[max(root_a, root_b)] = min(root_a, root_b)

    groups: Dict[int, List[str]] = {}
    for i, name in enumerate(names):
        groups.setdefault(find(i), []).append(name)
    return [groups[root] for root in sorted(groups)]


def read_quality(bin_dirs: Sequence[str], genomes: Sequence[Genome]) -> Dict[str, float]:
    """Each genome's quality score, from the CheckM table beside its bin directory if there is one.

    CheckM runs in this project write ``<bins>.stats.tsv`` next to the bins, so bins from
    bin_refinement need no re-run. Genomes with no row are left out instead of scored zero, so the
    caller can tell "scored badly" from "not scored".
    """
    by_source: Dict[str, Dict[str, float]] = {}
    for directory in bin_dirs:
        stats = os.path.normpath(directory) + STATS_SUFFIX
        if not os.path.isfile(stats):
            continue
        scores: Dict[str, float] = {}
        with open(stats) as fh:
            for line in fh:
                cut = line.rstrip("\n").split("\t")
                if len(cut) < 3 or cut[1] == "completeness":
                    continue
                try:
                    size = float(cut[6]) if len(cut) > 6 else 0.0
                    scores[cut[0]] = quality_score(float(cut[1]), float(cut[2]), size)
                except ValueError:  # a header this project did not write, or a partial row
                    continue
        by_source[os.path.abspath(directory)] = scores

    out: Dict[str, float] = {}
    for genome in genomes:
        scores = by_source.get(os.path.abspath(os.path.dirname(genome.path)), {})
        bin_name = bin_name_from_filename(genome.path)
        if bin_name in scores:
            out[genome.name] = scores[bin_name]
    return out


def genome_size(path: str) -> int:
    """Total bases in a fasta, used to break ties and to rank unscored genomes."""
    total = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def choose_representative(
    members: Sequence[str],
    scores: Dict[str, float],
    sizes: Dict[str, int],
) -> str:
    """The genome that stands for its cluster: best CheckM score, else largest.

    A scored genome always beats an unscored one regardless of size; otherwise a 6 Mb chimera could
    represent a cluster that contains a clean 3 Mb genome. Ties go to input order.
    """
    return max(
        members,
        key=lambda name: (
            name in scores,
            scores.get(name, 0.0),
            sizes.get(name, 0),
            -list(members).index(name),
        ),
    )


def ani_to_representative(
    members: Sequence[str],
    representative: str,
    edges: Iterable[Edge],
) -> Dict[str, float]:
    """Each member's ANI to its cluster's representative, where there is an edge to it.

    Used to spot over-merged clusters: a member with ANI below the threshold, or no edge to the
    representative, was joined through another member by single linkage.
    """
    wanted = set(members)
    out: Dict[str, float] = {}
    for edge in edges:
        if edge.ref == representative and edge.query in wanted:
            out[edge.query] = max(out.get(edge.query, 0.0), edge.ani)
        elif edge.query == representative and edge.ref in wanted:
            out[edge.ref] = max(out.get(edge.ref, 0.0), edge.ani)
    out[representative] = 100.0
    return out


# --- the orchestration --------------------------------------------------------------------


def _write_genome_list(genomes: Sequence[Genome], path: str) -> None:
    with open(path, "w") as fh:
        fh.writelines(genome.path + "\n" for genome in genomes)


def _preset_flag(preset: str) -> str:
    return "" if preset == "default" else " --" + preset


def _link_or_copy(src: str, dest: str, copy: bool) -> None:
    if os.path.lexists(dest):
        os.remove(dest)
    if copy:
        import shutil

        shutil.copyfile(src, dest)
    else:
        os.symlink(os.path.abspath(src), dest)


def _write_clusters(
    path: str,
    clusters: Sequence[Sequence[str]],
    representatives: Dict[int, str],
    genomes: Dict[str, Genome],
    scores: Dict[str, float],
    sizes: Dict[str, int],
    ani: Dict[int, Dict[str, float]],
) -> None:
    """One row per input genome, so the table accounts for every bin that went in."""
    with open(path, "w") as fh:
        fh.write("\t".join(CLUSTER_COLUMNS) + "\n")
        for number, members in enumerate(clusters, start=1):
            rep = representatives[number]
            for member in members:
                identity = ani[number].get(member)
                fh.write(
                    "%s\t%s\t%s\t%s\t%s\t%s\t%d\n"
                    % (
                        "cluster_%04d" % number,
                        member,
                        genomes[member].source,
                        "yes" if member == rep else "no",
                        "%.2f" % identity if identity is not None else "NA",
                        "%.2f" % scores[member] if member in scores else "NA",
                        sizes.get(member, 0),
                    )
                )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    for directory in args.bins:
        if not os.path.isdir(directory):
            error("%s does not exist. Exiting..." % directory)

    # A sample whose binning or refinement found nothing leaves an empty folder, which is common in
    # a cohort. Empty folders are skipped with a warning and the other folders are dereplicated.
    # The run stops only if every folder is empty.
    kept = []
    for directory in args.bins:
        if any(name.endswith(tuple(FASTA_EXTENSIONS)) for name in os.listdir(directory)):
            kept.append(directory)
        else:
            warning("%s has no bins in it; dereplicating without it." % directory)
    if not kept:
        error(
            "Every bin folder given is empty (%s), so there is nothing to dereplicate."
            % ", ".join(args.bins)
        )
    args.bins = kept

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "dereplicate",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                path,
                "bin folder (-b)",
            )
            for path in args.bins
        ],
    )
    env = env_for("dereplicate", settings)
    rec = start_run("dereplicate", args, env, settings, inputs=list(args.bins))

    ckpt = make_checkpoint(args.output)
    try:
        announcement("GATHERING BINS FROM %d FOLDER(S)" % len(args.bins))
        ensure_dir(args.output)
        genomes = collect_genomes(args.bins)
        if not genomes:
            error("no bins found in %s" % ", ".join(args.bins))
        comm("%d bins from %d folder(s)" % (len(genomes), len(args.bins)))

        genome_list = os.path.join(args.output, "genomes.txt")
        ani_table = os.path.join(args.output, ANI_TABLE)
        if not dry_run():
            _write_genome_list(genomes, genome_list)

        if ckpt.todo("ani"):
            announcement("COMPARING EVERY PAIR OF BINS WITH SKANI")
            comm(
                "clustering at %g%% ANI with at least %g%% of one genome aligned"
                % (args.ani, args.min_aligned_fraction)
            )
            run(
                SKANI_TRIANGLE.format(
                    genome_list=genome_list,
                    output=ani_table,
                    threads=args.threads,
                    screen=SKANI_SCREEN,
                    min_af=args.min_aligned_fraction,
                    preset=_preset_flag(args.preset),
                ),
                env=env,
                tool="skani",
                hint="Failed to compute ANI between bins with skani.",
            )
            ckpt.done("ani")
        else:
            comm("skipping skani (already done; --resume)")

        if dry_run():
            comm(
                "(dry run) would cluster the ANI table and write the catalogue to %s"
                % os.path.join(args.output, CATALOGUE_DIR)
            )
            announcement("DEREPLICATION PIPELINE FINISHED SUCCESSFULLY!!!")
            return 0

        announcement("CLUSTERING AND CHOOSING A REPRESENTATIVE PER CLUSTER")
        # skani names genomes by the path it was given, which is the path in genomes.txt.
        by_path = {os.path.abspath(g.path): g.name for g in genomes}
        edges = [
            edge._replace(
                ref=by_path.get(os.path.abspath(edge.ref), edge.ref),
                query=by_path.get(os.path.abspath(edge.query), edge.query),
            )
            for edge in parse_skani(ani_table)
        ]
        names = [g.name for g in genomes]
        clusters = cluster_genomes(
            names, edges, ani=args.ani, min_aligned_fraction=args.min_aligned_fraction
        )

        by_name = {g.name: g for g in genomes}
        scores = read_quality(args.bins, genomes)
        if not scores:
            warning(
                "no CheckM table beside any input folder, so representatives are chosen by\n"
                "genome size alone, and a large chimeric bin can be picked. Run bin_refinement\n"
                "(or CheckM) over these bins first; this module then uses the scores it writes\n"
                "to <bins>%s." % STATS_SUFFIX
            )
        sizes = {g.name: genome_size(g.path) for g in genomes}

        representatives: Dict[int, str] = {}
        identities: Dict[int, Dict[str, float]] = {}
        for number, members in enumerate(clusters, start=1):
            rep = choose_representative(members, scores, sizes)
            representatives[number] = rep
            identities[number] = ani_to_representative(members, rep, edges)

        catalogue = os.path.join(args.output, CATALOGUE_DIR)
        ensure_dir(catalogue)
        for number, members in enumerate(clusters, start=1):
            rep = representatives[number]
            _link_or_copy(
                by_name[rep].path,
                os.path.join(catalogue, rep + BIN_EXTENSION),
                args.copy,
            )

        clusters_table = os.path.join(args.output, CLUSTERS_TABLE)
        _write_clusters(
            clusters_table, clusters, representatives, by_name, scores, sizes, identities
        )

        collapsed = len(genomes) - len(clusters)
        comm(
            "%d bins collapsed into %d genomes (%d duplicates removed); catalogue in %s"
            % (len(genomes), len(clusters), collapsed, catalogue)
        )
        comm("what went into each cluster, and how close it was, is in %s" % clusters_table)
        announcement("DEREPLICATION PIPELINE FINISHED SUCCESSFULLY!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
