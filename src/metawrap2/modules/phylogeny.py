"""MetaWrap2 phylogeny: a tree from an alignment the pipeline already produced.

This module does no alignment. It builds a tree with support values from one of two alignments
the pipeline already writes:

* `classify_bins --gtdbtk` writes a concatenated marker alignment (120 bacterial or 53 archaeal
  single-copy genes), which places a catalogue's genomes relative to each other.
* `pangenome --msa` writes one species group's core-gene alignment, which resolves strains within
  that species.

The summary records which kind was given. Universal markers cannot resolve strains, and one
species' core genes do not place that species.

IQ-TREE is the default, with ModelFinder choosing the substitution model and ultrafast bootstrap
for support. `--fast` switches to FastTree, suited to a few hundred genomes, where model
selection dominates the cost.

See :mod:`metawrap2.phylogeny` on interpreting support values on a MAG phylogeny: they measure
the alignment and do not reflect incompleteness, transfer or contamination.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import List, Optional, Tuple

from .. import phylogeny as _phylo
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
    resolve_threads,
    run,
    start_run,
    threads_arg,
    warning,
)

CONDA_ENV = MODULE_ENVS["phylogeny"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# -m MFP is ModelFinder Plus: it picks the substitution model from the alignment, which on a marker
# concatenation changes branch lengths materially.
#
# -B is the ultrafast bootstrap; 1000 is IQ-TREE's minimum. --alrt adds the SH-like approximate
# likelihood ratio test, which is cheap and often disagrees with UFBoot; a branch supported by only
# one of them is doubtful.
#
# -T AUTO lets IQ-TREE's scaling test pick the thread count, since over-threading a small alignment
# slows it down. --threads-max caps it.
# The binary is `iqtree`: bioconda's IQ-TREE 3 installs `iqtree` and `iqtree3` but no `iqtree2`.
# The flags used here are the same in versions 2 and 3.
IQTREE = (
    "iqtree -s {alignment} --prefix {prefix} -m MFP -B {replicates} --alrt 1000"
    " -T AUTO --threads-max {threads} --seed {seed} -redo"
)
# Two templates, because FastTree's data type is set by a flag: `-nt` declares nucleotide, and
# without it FastTree reads protein. `-nt` on a protein alignment treats every residue other than
# A, C, G and T as ambiguous and returns a tree built from almost nothing, with no error. The
# upstream alignments (pangenome core genes, GTDB-Tk's MSA) are protein.
# :func:`metawrap2.phylogeny.read_alignment` detects the type, and the summary records the model.
# FastTree's SH-like supports are not bootstrap values; the summary records which tool built the
# tree.
FASTTREE_NT = "fasttree -nt -gtr -log {log} {alignment}"
#: -lg instead of FastTree's default JTT: LG fits concatenated single-copy markers better, and
#: FastTree has supported it since 2.1.5.
FASTTREE_AA = "fasttree -lg -log {log} {alignment}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per tree.
SUMMARY_TABLE = "trees" + TSV_SUFFIX
#: Where each tree's own files are kept.
TREES_DIR = "trees"
#: Fixed so two runs over one alignment give the same tree. IQ-TREE's search is stochastic, and an
#: unseeded rerun with a different topology would look like a real disagreement.
SEED = 1337
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 phylogeny",
        usage="metawrap2 phylogeny [options] -o out (-a alignment.fasta | -d search_dir)",
        description="Build a tree with support values from a marker or core-gene alignment.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-a",
        "--alignment",
        nargs="+",
        help="one or more alignments in FASTA (gzip accepted), such as GTDB-Tk's "
        "`*.user_msa.fasta` or the core-gene alignment from `pangenome --msa`",
    )
    source.add_argument(
        "-d",
        "--search",
        help="a directory to search for alignments: a GTDB-Tk output tree or a `pangenome` output "
        "directory",
    )
    p.add_argument(
        "--fast",
        action="store_true",
        help="use FastTree instead of IQ-TREE. Orders of magnitude faster, and suited to a few "
        "hundred genomes. Its SH-like values are not bootstrap support; the table labels them",
    )
    p.add_argument(
        "--replicates",
        type=int,
        default=_phylo.BOOTSTRAP_REPLICATES,
        metavar="N",
        help="ultrafast bootstrap replicates (default: %(default)s, IQ-TREE's minimum)",
    )
    p.add_argument(
        "--support-threshold",
        type=float,
        default=_phylo.SUPPORTED_BRANCH,
        metavar="PERCENT",
        help="support at or above which a branch is counted as resolved (default: %(default)s). "
        "For UFBoot this is about as strict as 70 for a standard bootstrap",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _alignments(args) -> List[str]:
    if args.alignment:
        for path in args.alignment:
            if not os.path.isfile(path):
                error("no alignment at %s" % path)
        return list(args.alignment)
    found = _phylo.find_alignments(args.search)
    if not found:
        error(
            "no alignment found under %s. This module reads GTDB-Tk's `*.user_msa.fasta` or the "
            "core-gene alignment written by `pangenome --msa`" % args.search
        )
    comm("%d alignment(s) found under %s" % (len(found), args.search))
    return found


def _tree_name(path: str) -> str:
    """A short, unique-enough name for a tree, from the alignment's own path.

    The parent directory is included because every group's core-gene alignment has the same
    filename.
    """
    base = os.path.basename(path)
    for extension in (".gz", ".fasta", ".fa", ".aln", ".faa", ".fna"):
        base = base.removesuffix(extension)
    parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
    return "%s__%s" % (parent, base) if parent else base


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)

    env = env_for("phylogeny", settings)
    rec = start_run("phylogeny", args, env, settings, inputs=list(args.alignment or []))
    try:
        announcement("BEGIN PHYLOGENY PIPELINE!")
        ensure_dir(args.output)
        paths = _alignments(args)

        if dry_run():
            comm(
                "would build %d tree(s) with %s"
                % (len(paths), "FastTree" if args.fast else "IQ-TREE")
            )
            announcement("PHYLOGENY PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        rows: List[Tuple[str, dict]] = []
        for path in bar(paths, desc="trees", unit=" alignment", total=len(paths)):
            name = _tree_name(path)
            try:
                alignment = _phylo.read_alignment(path)
            except (OSError, ValueError) as problem:
                warning("could not read %s: %s" % (path, problem))
                continue
            reason = _phylo.too_few_taxa(alignment)
            if reason:
                warning(reason)
                continue
            comm(
                "%s: %d taxa, %d columns, mean occupancy %.2f"
                % (
                    name,
                    alignment.size,
                    alignment.columns,
                    alignment.mean_occupancy or 0.0,
                )
            )
            if alignment.sparse_taxa:
                warning(
                    "%s: %d taxon/taxa fill less than %.0f%% of the alignment (%s)"
                    % (
                        name,
                        len(alignment.sparse_taxa),
                        _phylo.OCCUPANCY_FLOOR * 100,
                        ", ".join(alignment.sparse_taxa[:5]),
                    )
                )
                comm(_phylo.OCCUPANCY_CAVEAT)

            treedir = os.path.join(args.output, TREES_DIR)
            ensure_dir(treedir)
            prefix = os.path.join(treedir, name)
            model: Optional[str] = None
            if args.fast:
                newick = prefix + ".treefile"
                nucleotide = alignment.kind == "nucleotide"
                template = FASTTREE_NT if nucleotide else FASTTREE_AA
                comm(
                    "%s appears to be %s, so FastTree runs with %s"
                    % (name, alignment.kind, "-nt -gtr" if nucleotide else "-lg")
                )
                rc = run(
                    template.format(log=prefix + ".log", alignment=path),
                    env=env,
                    tool="fasttree",
                    check=False,
                    stdout_path=newick,
                )
                model = "GTR (FastTree, fixed)" if nucleotide else "LG (FastTree, fixed)"
            else:
                rc = run(
                    IQTREE.format(
                        alignment=path,
                        prefix=prefix,
                        replicates=args.replicates,
                        threads=args.threads,
                        seed=SEED,
                    ),
                    env=env,
                    tool="iqtree2",
                    check=False,
                )
                newick = prefix + ".treefile"
                model = _phylo.model_from_report(prefix + ".iqtree")
            if rc != 0 or not os.path.isfile(newick):
                warning("tree building failed for %s" % name)
                continue
            support = _phylo.support_summary(newick, threshold=args.support_threshold)
            rows.append(
                (
                    name,
                    _phylo.summary_row(
                        name,
                        alignment,
                        model,
                        support,
                        replicates=0 if args.fast else args.replicates,
                    ),
                )
            )
            fraction = support.fraction_supported
            comm(
                "  %s: %s, %d of %d internal branches at or above %g%s"
                % (
                    name,
                    model or "model not recorded",
                    support.supported,
                    support.branches,
                    args.support_threshold,
                    "" if fraction is None else " (%.0f%%)" % (100.0 * fraction),
                )
            )

        if not rows:
            error("no tree was built")

        path = os.path.join(args.output, SUMMARY_TABLE)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_phylo.SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for _name, row in sorted(rows):
                writer.writerow(row)
        comm("%d tree(s) -> %s" % (len(rows), path))
        comm(_phylo.SUPPORT_CAVEAT)
        if any("user_msa" in row["alignment"] for _name, row in rows):
            comm(_phylo.SCOPE_CAVEAT)
        announcement("PHYLOGENY PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
