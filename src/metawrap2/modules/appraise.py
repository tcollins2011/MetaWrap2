"""MetaWrap2 appraise: what share of the community the bins represent, and what has no genome at all.

A bin set reports what it found. The question it cannot answer about itself is what it missed, and
the binned fraction of the assembly - which `bin_refinement` already reports - is not that answer:
it is a fraction of *assembled base pairs*, and what assembles well and what is abundant are
different things. A bin set can hold 60% of the assembly and stand for 30% of the cells.

SingleM answers the per-cell version from the reads directly. It finds single-copy marker genes in
the unassembled reads, finds the same markers in the genomes, and compares: each marker sequence in
the reads is either present in a genome (binned), present in the assembly but no genome, or absent
from both. Weighted by read coverage, that gives the share of the community the genomes represent -
and, from the leftover, a list of lineages present in the sample with no genomic representation at
all. The second output is the more useful one for deciding what to do next.

Three steps, all of them SingleM's: profile the reads, profile the genomes, appraise one against
the other. The module runs all three because the intermediate OTU tables are not interesting on
their own, and keeps them anyway - they are expensive and SingleM can do more with them.

See :mod:`metawrap2.recovery` for the two caveats this prints, and for the assembled-bp fraction
that answers the other question.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import List, Optional

from .. import recovery as _recovery
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    require_file,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["appraise"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# Markers in the reads. Forward reads only: SingleM treats a pair as two observations of the same
# cell unless told the library layout, and a marker read and its mate would then count twice.
SINGLEM_PIPE_READS = (
    "singlem pipe -1 {reads} --otu-table {out} --metapackage {metapackage} --threads {threads}"
)
# The same markers in the genomes. -d/-x rather than a file list because a bin folder is what every
# other module in this pipeline passes around; SingleM raises its own --min-taxon-coverage and
# --min-orf-length defaults for genome input, which is why this is not just --forward on a fasta.
SINGLEM_PIPE_GENOMES = (
    "singlem pipe -d {bins} -x {extension} --otu-table {out} --metapackage {metapackage}"
    " --threads {threads}"
)
# The comparison. --output-unaccounted-for-otu-table is the output worth having: the markers present
# in the reads that no genome holds, which is the list of lineages this bin set says nothing about.
SINGLEM_APPRAISE = (
    "singlem appraise --metagenome-otu-tables {metagenome} --genome-otu-tables {genomes}"
    " --metapackage {metapackage} --output-unaccounted-for-otu-table {unaccounted}"
    " --output-binned-otu-table {binned} --output-found-in"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per sample per domain: what share of the community the genomes stand for.
APPRAISAL_TABLE = "community_represented" + TSV_SUFFIX
#: Lineages in the reads with no genome representing them, most abundant first.
UNREPRESENTED_TABLE = "unrepresented_lineages" + TSV_SUFFIX
#: SingleM's own tables, kept because they cost minutes and SingleM can do more with them than this.
OTU_DIR = "singlem"
#: How many lineages to print. The table holds all of them.
PRINTED_LINEAGES = 10
#: SingleM compares a fixed 60 bp window of each marker gene, so a read shorter than that cannot
#: cover one and contributes nothing whatever its quality. Measured, not assumed: the OTU tables
#: SingleM writes carry 60 bp sequences.
#:
#: This is the difference between the right input and the wrong one. One real gut sample's *trimmed*
#: reads (mean 34 bp) yielded 10 marker sequences from 12 million reads, which would appraise as a
#: bin set representing almost nothing; the same sample's raw 75 bp reads are all long enough. So
#: this module wants reads that are quality-filtered but not length-trimmed below the window.
SINGLEM_WINDOW = 60
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 appraise",
        usage="metawrap2 appraise [options] -o out -b bin_folder -r reads_1.fastq [reads_1.fastq ...]",
        description="What fraction of each sample's community the bins represent, per cell.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-b",
        "--bins",
        required=True,
        help="the genomes to appraise - a refined bin set, or a dereplicated catalogue",
    )
    p.add_argument(
        "-r",
        "--reads",
        required=True,
        nargs="+",
        metavar="FILE",
        help="one read file per sample, forward mate only. SingleM counts a marker read and its "
        "mate as two observations of one cell, so passing both inflates coverage",
    )
    p.add_argument(
        "--metapackage",
        help="path to the SingleM metapackage (default: SINGLEM_METAPACKAGE from metawrap2.toml)",
    )
    p.add_argument(
        "--extension",
        default=None,
        help="bin file extension (default: inferred from the folder)",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _singlem_label(path: str) -> str:
    """The sample name SingleM will put in its OTU table for *path*.

    It labels rows with the file's basename minus its extensions, so a forward read file becomes
    ``ERR011347_1`` where every other module in this pipeline calls that sample ``ERR011347``. The
    mapping is needed, not cosmetic: the appraisal table is joined against bin and function tables
    on the sample name, and a row labelled ``ERR011347_1`` joins to nothing.
    """
    name = os.path.basename(path)
    for extension in (".gz", ".bz2"):
        name = name.removesuffix(extension)
    return os.path.splitext(name)[0]


def _sample_name(path: str) -> str:
    name = os.path.basename(path)
    for extension in (".gz", ".bz2"):
        name = name.removesuffix(extension)
    name = os.path.splitext(name)[0]
    for suffix in ("_final_pure_reads_1", "_1", "_R1", ".1"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _infer_extension(folder: str, given: Optional[str]) -> str:
    """The bin extension SingleM should look for, from the folder itself.

    Asked of the folder rather than defaulted, because SingleM's own default is ``fna`` and
    MetaWrap2 writes ``.fasta`` - so the default silently appraises zero genomes and the appraisal
    then reports that none of the community is represented, which looks like a result.
    """
    if given:
        return given.lstrip(".")
    counts: dict = {}
    for name in sorted(os.listdir(folder)):
        _stem, dot, extension = name.rpartition(".")
        if dot and extension:
            counts[extension] = counts.get(extension, 0) + 1
    if not counts:
        error("no files with an extension in %s, so there are no genomes to appraise" % folder)
    best = max(counts, key=lambda extension: counts[extension])
    if len(counts) > 1:
        warning(
            "%s holds more than one extension (%s); appraising the %d .%s file(s). Pass "
            "--extension to choose"
            % (
                folder,
                ", ".join(sorted(counts)),
                counts[best],
                best,
            )
        )
    return best


def _check_read_lengths(reads: List[str]) -> None:
    """Warn when the reads are too short for SingleM's window, which is silent otherwise.

    The failure this prevents does not look like a failure: short reads produce an almost empty
    marker profile, and an almost empty profile appraises as "these genomes represent none of the
    community". Reusing nonpareil's scanner rather than writing a second one - it reads the head of
    the file and answers the same question.
    """
    from .. import nonpareil as _np

    for path in reads:
        lengths = _np.scan_read_lengths(path, k=SINGLEM_WINDOW, limit=100_000)
        if not lengths.scanned:
            continue
        comm(
            "%s: %d reads scanned, mean length %.0f bp"
            % (os.path.basename(path), lengths.scanned, lengths.mean)
        )
        if lengths.fraction_short > 0.5:
            warning(
                "%s: %.0f%% of its reads are shorter than SingleM's %d bp window and cannot carry "
                "a marker. The appraisal below will understate what the genomes represent, "
                "possibly to nothing. Use reads that were quality-filtered but not trimmed below "
                "%d bp - on one real sample, trimmed reads (mean 34 bp) gave 10 marker sequences "
                "from 12 million reads where the raw 75 bp reads gave thousands"
                % (
                    os.path.basename(path),
                    lengths.fraction_short * 100,
                    SINGLEM_WINDOW,
                    SINGLEM_WINDOW,
                )
            )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    absolutize_paths(args)
    validate_inputs(
        "appraise",
        [(require_nonempty_dir, args.bins, "bin folder (-b)")]
        + [(require_file, path, "read file (-r)") for path in args.reads],
    )
    resolve_threads(args, settings)

    metapackage = args.metapackage or settings.db("SINGLEM_METAPACKAGE")
    if not metapackage:
        error(
            "no SingleM metapackage: pass --metapackage or set SINGLEM_METAPACKAGE in "
            "metawrap2.toml. `metawrap2 install-db SINGLEM_METAPACKAGE` downloads it (2.2 GB "
            "compressed, 9.7 GB on disk)"
        )

    env = env_for("appraise", settings)
    rec = start_run("appraise", args, env, settings, inputs=[args.bins] + list(args.reads))
    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN APPRAISE PIPELINE!")
        ensure_dir(args.output)
        otu_dir = os.path.join(args.output, OTU_DIR)
        ensure_dir(otu_dir)

        _check_read_lengths(args.reads)
        extension = _infer_extension(args.bins, args.extension)
        genome_table = os.path.join(otu_dir, "genomes.otu" + TSV_SUFFIX)
        metagenome_table = os.path.join(otu_dir, "metagenome.otu" + TSV_SUFFIX)

        if dry_run():
            comm(
                "would profile %d read file(s) and the .%s genomes in %s, then appraise one "
                "against the other" % (len(args.reads), extension, args.bins)
            )
            announcement("APPRAISE PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        if ckpt.todo("genomes"):
            rc = run(
                SINGLEM_PIPE_GENOMES.format(
                    bins=args.bins,
                    extension=extension,
                    out=genome_table,
                    metapackage=metapackage,
                    threads=args.threads,
                ),
                env=env,
                tool="singlem pipe (genomes)",
                check=False,
            )
            if rc != 0 or not os.path.isfile(genome_table):
                error("singlem pipe failed on the genomes; nothing to appraise against")
                return 1
            ckpt.done("genomes")
        else:
            comm("skipping the genome profile (already done; --resume)")

        # One OTU table per sample, concatenated: `pipe` labels rows with the sample name it derives
        # from the file, and appraise takes several metagenome tables, so this keeps the per-sample
        # rows distinguishable without relying on how pipe names them.
        per_sample: List[str] = []
        labels: dict = {}
        for path in args.reads:
            sample = _sample_name(path)
            labels[_singlem_label(path)] = sample
            table = os.path.join(otu_dir, "%s.otu%s" % (sample, TSV_SUFFIX))
            if ckpt.todo("reads:%s" % sample):
                comm("%s: finding marker genes in the reads" % sample)
                rc = run(
                    SINGLEM_PIPE_READS.format(
                        reads=path,
                        out=table,
                        metapackage=metapackage,
                        threads=args.threads,
                    ),
                    env=env,
                    tool="singlem pipe (reads)",
                    check=False,
                )
                if rc != 0 or not os.path.isfile(table):
                    warning("singlem pipe failed on %s; that sample is not appraised" % sample)
                    continue
                ckpt.done("reads:%s" % sample)
            per_sample.append(table)
            markers = _count_rows(table)
            if markers == 0:
                warning(
                    "%s: no marker genes found. Single-copy markers are about 0.05%% of reads, so "
                    "a shallow or subsampled file can miss them entirely - this is a depth "
                    "problem, not an empty community" % sample
                )

        if not per_sample:
            error("no sample produced a marker profile, so there is nothing to appraise")
            return 1

        with open(metagenome_table, "w") as out:
            for index, table in enumerate(per_sample):
                with open(table) as fh:
                    lines = fh.readlines()
                out.writelines(lines if index == 0 else lines[1:])

        unaccounted = os.path.join(otu_dir, "unaccounted_for.otu" + TSV_SUFFIX)
        binned = os.path.join(otu_dir, "binned.otu" + TSV_SUFFIX)
        summary = os.path.join(otu_dir, "appraise_summary" + TSV_SUFFIX)
        rc = run(
            SINGLEM_APPRAISE.format(
                metagenome=metagenome_table,
                genomes=genome_table,
                metapackage=metapackage,
                unaccounted=unaccounted,
                binned=binned,
            ),
            env=env,
            tool="singlem appraise",
            check=False,
            # appraise writes its summary table to stdout; the runner captures it to a file rather
            # than this command going through a shell for a redirect.
            stdout_path=summary,
        )
        if rc != 0 or not os.path.isfile(summary):
            error("singlem appraise failed; its own output above says why")
            return 1

        appraisals = _recovery.parse_appraisal(summary)
        table = os.path.join(args.output, APPRAISAL_TABLE)
        with open(table, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=_recovery.APPRAISAL_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for appraisal in appraisals:
                row = _recovery.appraisal_row(appraisal)
                row["sample"] = labels.get(appraisal.sample, appraisal.sample)
                writer.writerow(row)
        comm("-> %s" % table)

        for appraisal in appraisals:
            if appraisal.total == 0:
                continue
            comm(
                "%s (%s): the genomes represent %.1f%% of the community"
                % (
                    labels.get(appraisal.sample, appraisal.sample),
                    appraisal.domain or "all domains",
                    appraisal.percent_binned or 0.0,
                )
            )
        comm(_recovery.APPRAISAL_CAVEAT)

        if os.path.isfile(unaccounted):
            lineages = _recovery.unrepresented_lineages(unaccounted)
            out_path = os.path.join(args.output, UNREPRESENTED_TABLE)
            with open(out_path, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=_recovery.UNREPRESENTED_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in _recovery.unrepresented_rows(lineages):
                    writer.writerow(row)
            comm("%d unrepresented lineage(s) -> %s" % (len(lineages), out_path))
            for lineage in lineages[:PRINTED_LINEAGES]:
                comm("  %-28s coverage %.3g" % (lineage.rank_reached, lineage.coverage))
            if len(lineages) > PRINTED_LINEAGES:
                comm("  ... and %d more in the table" % (len(lineages) - PRINTED_LINEAGES))

        announcement("APPRAISE PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


def _count_rows(path: str) -> int:
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
