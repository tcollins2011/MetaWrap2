"""MetaWrap2 dereplicate module: collapse many samples' bins into one non-redundant catalogue.

Binning one cohort per sample produces the same genome once per sample it occurs in. Twenty
samples of one gut cohort can yield two thousand bins describing a few hundred organisms, and
every downstream step pays for the duplication: abundance estimation has near-identical genomes
stealing each other's reads, annotation runs on the same genome twenty times, and a MAG count
reported off that set is not a count of organisms.

So the bins are clustered by whole-genome identity and one representative is kept per cluster.
Identity comes from skani (Shaw & Yu, *Nat Methods* 2023), which estimates ANI from sparse k-mer
sketches rather than alignment and reports the aligned fraction alongside it; the clustering and
the choice of representative are done here, where they can be read and tested.

Two numbers decide the result and both are flags rather than findings - see SPECIES_ANI and
MIN_ALIGNED_FRACTION. Neither has been measured on this project's data yet.

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
#: ``-E`` is the sparse form: one row per pair instead of a square matrix, which is what the
#: clustering wants and is the difference between a file of N rows and one of N^2 cells on a
#: two-thousand-genome catalogue. ``--min-af`` filters at the source so a pair sharing only a
#: plasmid never reaches the clustering, and ``-s`` is skani's own screen - pairs below it are not
#: compared at all, which is where its speed comes from.
SKANI_TRIANGLE = (
    "skani triangle -E -l {genome_list} -o {output} -t {threads} "
    "-s {screen} --min-af {min_af}{preset}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: ANI at which two bins are called the same species. 95% is the convention the large MAG
#: catalogues are built on, which is the only reason to start here: it makes a species count
#: comparable with a published one. It is not a measurement of ours, and the right value depends on
#: what the catalogue is for - 99 separates strains and is what an inStrain-style analysis wants.
SPECIES_ANI = 95.0
#: How much of a genome has to align before an ANI estimate is allowed to speak for the whole
#: genome, as a percentage. Deliberately permissive: at 95% identity the ANI estimate is doing the
#: species call, and this filter exists only to stop a short shared region - a prophage, a plasmid,
#: a conserved operon - standing in for a genome-wide comparison. Raising it makes clusters smaller
#: and the catalogue larger, which is a real choice and not a correctness fix.
MIN_ALIGNED_FRACTION = 10.0
#: skani's pre-filter, in approximate percent identity. Anything it screens out cannot be in a
#: cluster, so it has to sit below SPECIES_ANI with room for the approximation; skani's own default
#: is 80, which is far enough below 95 to be safe and is left alone.
SKANI_SCREEN = 80
#: skani presets trade compression factor for accuracy on fragmented assemblies. MAGs are
#: fragmented, which argues for ``medium``, but skani's help puts that at an N50 below 10 kb and
#: costs 2x for it, so which one a MAG catalogue actually needs is a measurement nobody here has
#: made. Left at skani's default until it is one.
SKANI_PRESET = "default"
SKANI_PRESETS = ("default", "fast", "medium", "slow")
#: Where the kept genomes and the record of what was collapsed into them are written.
CATALOGUE_DIR = "dereplicated_genomes"
CLUSTERS_TABLE = "clusters" + TSV_SUFFIX
ANI_TABLE = "ani" + TSV_SUFFIX
#: The cluster table's columns. ``ani_to_representative`` is in there on purpose - it is what makes
#: an over-merged cluster visible rather than a number nobody can question.
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

        The larger, not the smaller: a 60% complete bin of the same organism aligns almost
        entirely to its complete neighbour while the neighbour aligns only 60% to it, and taking
        the smaller would refuse to collapse exactly the duplicates this module exists to collapse.
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
        help="ANI%% at which two bins are the same genome - the species boundary; "
        "99 gives strain-level clusters instead",
    )
    p.add_argument(
        "--min-aligned-fraction",
        type=percent(),
        default=MIN_ALIGNED_FRACTION,
        help="minimum %% of either genome that must align for an ANI value to be believed",
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

    Bins are numbered within one binner's output, so ``bin_001.fasta`` exists in every sample and
    the names cannot be used as they are - skani would be handed several files with the same label
    and the clustering could not say which sample a genome came from. The source directory's name
    becomes a prefix, which is also what makes the catalogue's filenames say where each kept genome
    was found.
    """
    genomes: List[Genome] = []
    seen: Set[str] = set()
    for directory in bin_dirs:
        source = os.path.basename(os.path.normpath(directory))
        for filename in list_bins(directory, FASTA_EXTENSIONS):
            name = "%s__%s" % (source, bin_name_from_filename(filename))
            # Two input directories can share a basename (one per sample, all called
            # "metawrap_70_10_bins"). Disambiguating here rather than refusing keeps the common
            # cohort layout usable without asking the caller to rename anything.
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

    The columns are named in the header rather than assumed by position, because skani's ``--ci``
    and ``--detailed`` flags insert columns and a positional read would silently take the wrong
    one if anybody added either to the command above.
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

    Single linkage, which is a choice with a known cost: identity is not transitive, so a chain of
    genomes each 95% identical to the next can be pulled into one cluster whose ends are not 95%
    identical to each other. The alternative is to require every pair inside a cluster to clear the
    threshold, which needs the full matrix and makes the answer depend on the order genomes are
    considered in. Single linkage is taken here because it never splits a genuine duplicate - which
    is the failure this module exists to prevent - and because what it over-merges can be seen: the
    cluster table records every member's ANI to its representative, so a chained cluster shows up
    as a member well below the threshold.

    Genomes with no edge to anything are their own cluster, so the output always accounts for every
    name given. Returned in input order, each cluster's members in input order, so two runs over
    the same bins give the same catalogue.
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
            # Lower index wins, so a cluster's root is its earliest member and the ordering below
            # follows the input rather than the order edges happened to arrive in.
            parent[max(root_a, root_b)] = min(root_a, root_b)

    groups: Dict[int, List[str]] = {}
    for i, name in enumerate(names):
        groups.setdefault(find(i), []).append(name)
    return [groups[root] for root in sorted(groups)]


def read_quality(bin_dirs: Sequence[str], genomes: Sequence[Genome]) -> Dict[str, float]:
    """Each genome's quality score, from the CheckM table beside its bin directory if there is one.

    ``<bins>.stats.tsv`` is what this project's CheckM runs leave next to the bins they scored, so
    a cohort that has been through bin_refinement already has the numbers and nothing needs to be
    re-run. Genomes with no row are absent from the result rather than scored zero, which lets the
    caller tell "scored badly" from "never scored" - they pick different representatives.
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
    """Total bases in a fasta, used to break ties and to rank genomes nobody scored."""
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

    A scored genome always beats an unscored one, however large: size says a bin is big, not that
    it is one organism, and a 6 Mb chimera would otherwise represent a cluster containing the clean
    3 Mb version of the same genome. Ties fall to the input order, which keeps the choice
    reproducible.
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

    Recorded because it is how an over-merged cluster is spotted: single linkage can chain, and a
    member whose ANI to the representative is below the threshold - or absent entirely - was pulled
    in through somebody else.
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
                "no CheckM table beside any input folder, so representatives are being chosen on\n"
                "genome size alone. Size says a bin is large, not that it is one organism, so a\n"
                "chimera can win its cluster. Run bin_refinement (or CheckM) over these bins first\n"
                "and this module will use the scores it leaves in <bins>%s." % STATS_SUFFIX
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
