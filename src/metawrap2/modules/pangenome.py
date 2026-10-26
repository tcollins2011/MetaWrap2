"""MetaWrap2 pangenome: one pangenome per species cluster, from dereplication's clusters.

`dereplicate` keeps one representative per organism. The copies of a species recovered from
different samples differ in gene content, and this module builds a pangenome from each cluster's
members to describe that variation.

A pangenome is defined over one species, so the input is `dereplicate`'s `clusters.tsv` plus the
bin folders it was built from, and the whole set is never treated as one group. `--genomes` takes a
single folder known to hold one species.

With few samples most clusters are singletons (almost all of them on three samples), so many
cohorts produce no pangenome; this is reported. See :mod:`metawrap2.pangenome` for the minimum
group size and why the partitioning needs about 15 genomes to be reliable.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import List

from .. import pangenome as _pan
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..progress import bar
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    require_file,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["pangenome"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# `workflow` runs annotation, clustering, the graph, partitioning and flat-file output in one call.
#
# --fasta, not --anno: PPanGGOLiN calls genes with Prodigal consistently across the group, which
# gene-family clustering needs.
PPANGGOLIN = "ppanggolin workflow --fasta {genome_list} -o {outdir} -c {threads} --force"
# Core-gene alignment for `phylogeny`, written only with --msa because it is a second pass.
# --phylo writes one concatenated alignment of the core families. --single_copy excludes families
# present more than once in a genome, since paralogs would align different genes to each other.
PPANGGOLIN_MSA = (
    "ppanggolin msa -p {pangenome} -o {outdir} -c {threads} --partition core --source dna"
    " --phylo --single_copy --force"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per species group: how many genomes, and the gene-content partitioning.
SUMMARY_TABLE = "pangenomes" + TSV_SUFFIX
#: Where each group's PPanGGOLiN output tree lives.
GROUPS_DIR = "groups"
#: Its pangenome database, which `ppanggolin msa` reads.
PANGENOME_FILE = "pangenome.h5"
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 pangenome",
        usage="metawrap2 pangenome [options] -o out -c clusters.tsv -b bins [bins ...]",
        description="A pangenome per species in a cohort, from dereplication's clusters.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-c",
        "--clusters",
        help="`dereplicate`'s clusters.tsv. Each cluster with enough members becomes one "
        "pangenome (one per species, never one over the whole catalogue)",
    )
    p.add_argument(
        "-b",
        "--bins",
        nargs="+",
        help="the bin folder(s) the clusters were built from - every genome, not the "
        "dereplicated representatives",
    )
    p.add_argument(
        "-g",
        "--genomes",
        help="a single folder of genomes already known to be one species, instead of --clusters",
    )
    p.add_argument(
        "--min-genomes",
        type=int,
        default=_pan.MIN_GENOMES,
        metavar="N",
        # Concatenated, not %-formatted, because argparse interpolates %(default)s itself.
        help="smallest group to build a pangenome for (default: %(default)s). Below "
        + str(_pan.RELIABLE_GENOMES)
        + " genomes the persistent/shell/cloud partitioning is marked provisional in the "
        "output, whatever this is set to",
    )
    p.add_argument(
        "--msa",
        action="store_true",
        help="also write each group's core-gene alignment, which `metawrap2 phylogeny` takes",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _groups(args) -> List[_pan.GenomeGroup]:
    """The species groups to build, from clusters or from a single folder."""
    if args.genomes:
        paths = _pan.genome_paths([args.genomes])
        if not paths:
            error("no genomes in %s" % args.genomes)
        group = _pan.GenomeGroup(name=os.path.basename(args.genomes.rstrip("/")), genomes=paths)
        if not group.is_a_pangenome:
            error(
                "%d genome(s) in %s: a pangenome needs at least %d, since one or two genomes "
                "cannot be partitioned into persistent, shell and cloud genes"
                % (group.size, args.genomes, _pan.MIN_GENOMES)
            )
        return [group]

    clusters = _pan.read_clusters(args.clusters)
    paths = _pan.genome_paths(args.bins or [])
    if not paths:
        error("no genomes under %s" % ", ".join(args.bins or ["(nothing given)"]))
    groups, census = _pan.groups_from_clusters(clusters, paths, min_genomes=args.min_genomes)
    for line in _pan.census_lines(census):
        comm(line)
    if not groups:
        comm(_pan.NO_GROUPS_CAVEAT)
    return groups


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    if not args.clusters and not args.genomes:
        error("nothing to build: pass --clusters with --bins, or --genomes")
    if args.clusters and not args.bins:
        error("--clusters needs --bins: the cluster table lists genome names, not the genomes")
    if args.clusters:
        validate_inputs("pangenome", [(require_file, args.clusters, "cluster table (-c)")])

    env = env_for("pangenome", settings)
    rec = start_run("pangenome", args, env, settings, inputs=[args.clusters or args.genomes])
    try:
        announcement("BEGIN PANGENOME PIPELINE!")
        ensure_dir(args.output)
        groups = _groups(args)

        if dry_run():
            comm("would build %d pangenome(s)" % len(groups))
            announcement("PANGENOME PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0
        if not groups:
            announcement("PANGENOME PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        work = os.path.join(args.output, "work_files")
        ensure_dir(work)
        # PPanGGOLiN creates its output directory with a non-recursive mkdir, so the parent must
        # exist.
        ensure_dir(os.path.join(args.output, GROUPS_DIR))
        summaries: List[_pan.PangenomeSummary] = []
        for group in bar(groups, desc="pangenomes", unit=" group", total=len(groups)):
            outdir = os.path.join(args.output, GROUPS_DIR, group.name)
            listing = _pan.write_genome_list(
                group, os.path.join(work, "%s_genomes%s" % (group.name, TSV_SUFFIX))
            )
            rc = run(
                PPANGGOLIN.format(genome_list=listing, outdir=outdir, threads=args.threads),
                env=env,
                tool="ppanggolin workflow",
                check=False,
            )
            if rc != 0:
                warning("PPanGGOLiN failed on %s; that group has no pangenome" % group.name)
                continue
            summary = _pan.parse_partitions(outdir, group.name, group.is_reliable)
            if summary is None:
                warning(
                    "PPanGGOLiN wrote no partition files for %s, so its partitioning is unknown"
                    % group.name
                )
                continue
            summaries.append(summary)
            if args.msa:
                rc = run(
                    PPANGGOLIN_MSA.format(
                        pangenome=os.path.join(outdir, PANGENOME_FILE),
                        outdir=os.path.join(outdir, "msa"),
                        threads=args.threads,
                    ),
                    env=env,
                    tool="ppanggolin msa",
                    check=False,
                )
                if rc != 0:
                    warning("could not write a core-gene alignment for %s" % group.name)

        if not summaries:
            error("no group produced a pangenome")

        path = os.path.join(args.output, SUMMARY_TABLE)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_pan.SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for summary in sorted(summaries, key=lambda s: (-s.genomes, s.group)):
                writer.writerow(_pan.summary_row(summary))
        comm("%d pangenome(s) -> %s" % (len(summaries), path))
        for summary in sorted(summaries, key=lambda s: -s.genomes)[:10]:
            shares = summary.shares
            comm(
                "  %s: %d genomes, %d families (%.0f%% persistent, %.0f%% shell, %.0f%% cloud)%s"
                % (
                    summary.group,
                    summary.genomes,
                    summary.total_families,
                    shares.get(_pan.PERSISTENT, 0.0),
                    shares.get(_pan.SHELL, 0.0),
                    shares.get(_pan.CLOUD, 0.0),
                    "" if summary.reliable else "  [provisional]",
                )
            )
        if any(not summary.reliable for summary in summaries):
            comm(_pan.SMALL_GROUP_CAVEAT)
        suspect = _pan.suspect_groups(summaries)
        if suspect:
            warning(
                "%d group(s) share less than %g%% of their gene families: %s"
                % (
                    len(suspect),
                    _pan.SUSPECT_PERSISTENT,
                    ", ".join(
                        "%s (%.0f%% persistent, %d families over %.0f genes per genome)"
                        % (
                            s.group,
                            s.shares.get(_pan.PERSISTENT, 0.0),
                            s.total_families,
                            s.mean_genes or 0.0,
                        )
                        for s in suspect
                    ),
                )
            )
            comm(_pan.SUSPECT_CLUSTER_CAVEAT)
        announcement("PANGENOME PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
