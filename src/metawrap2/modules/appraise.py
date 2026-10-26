"""MetaWrap2 appraise: the share of the community the bins represent, and lineages with no genome.

The binned fraction of the assembly that `bin_refinement` reports is a fraction of assembled base
pairs, and assembly and abundance differ: a bin set can hold 60% of the assembly and represent 30%
of the cells.

SingleM finds single-copy marker genes in the unassembled reads and in the genomes and compares
them. Each marker sequence in the reads is present in a genome (binned), present in the assembly
but no genome, or absent from both. Weighted by read coverage this gives the share of the
community the genomes represent, and the remainder gives a list of lineages in the sample with no
genome.

The module runs SingleM's three steps (profile the reads, profile the genomes, appraise) and keeps
the intermediate OTU tables for further SingleM use.

See :mod:`metawrap2.recovery` for the two caveats this prints and for the assembled-bp fraction.
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
# Markers in the reads. Forward reads only: SingleM counts a read and its mate as two observations
# of one cell unless told the library layout.
SINGLEM_PIPE_READS = (
    "singlem pipe -1 {reads} --otu-table {out} --metapackage {metapackage} --threads {threads}"
)
# The same markers in the genomes. -d/-x takes a bin folder, as the other modules pass around.
# SingleM uses higher --min-taxon-coverage and --min-orf-length defaults for genome input than for
# --forward on a FASTA.
SINGLEM_PIPE_GENOMES = (
    "singlem pipe -d {bins} -x {extension} --otu-table {out} --metapackage {metapackage}"
    " --threads {threads}"
)
# The comparison. --output-unaccounted-for-otu-table lists the markers in the reads that no genome
# holds: the lineages the bin set does not represent.
SINGLEM_APPRAISE = (
    "singlem appraise --metagenome-otu-tables {metagenome} --genome-otu-tables {genomes}"
    " --metapackage {metapackage} --output-unaccounted-for-otu-table {unaccounted}"
    " --output-binned-otu-table {binned} --output-found-in"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per sample per domain: the share of the community the genomes represent.
APPRAISAL_TABLE = "community_represented" + TSV_SUFFIX
#: Lineages in the reads with no genome representing them, most abundant first.
UNREPRESENTED_TABLE = "unrepresented_lineages" + TSV_SUFFIX
#: SingleM's own OTU tables, kept for further SingleM use.
OTU_DIR = "singlem"
#: How many lineages to print. The table holds all of them.
PRINTED_LINEAGES = 10
#: SingleM compares a fixed 60 bp window of each marker gene (its OTU tables carry 60 bp
#: sequences), so a shorter read contributes nothing.
#:
#: One gut sample's trimmed reads (mean 34 bp) gave 10 marker sequences from 12 million reads; the
#: raw 75 bp reads are all long enough. Input reads should be quality-filtered but not trimmed below
#: the window.
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
        help="one read file per sample, forward mate only. SingleM would count a read and its "
        "mate as two cells, so passing both mates inflates coverage",
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

    SingleM labels rows with the file's basename minus its extensions (``ERR011347_1``), while
    other modules call the sample ``ERR011347``. The appraisal table is joined to bin and function
    tables on the sample name, so the labels are mapped back.
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

    SingleM's default is ``fna`` and MetaWrap2 writes ``.fasta``; with the default, SingleM finds
    no genomes and reports none of the community as represented.
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
    """Warn when most reads are shorter than SingleM's window.

    Short reads give an almost empty marker profile, which appraises as the genomes representing
    almost none of the community. Uses nonpareil's read-length scanner on the head of each file.
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
                "a marker, so the appraisal will underestimate what the genomes represent, "
                "possibly down to zero. Use reads that are quality-filtered but not trimmed below "
                "%d bp. On one sample, trimmed reads (mean 34 bp) gave 10 marker sequences from 12 "
                "million reads, and the untrimmed 75 bp reads gave thousands"
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

        # One OTU table per sample, concatenated afterwards. `pipe` labels rows with a sample name
        # derived from the file name; `labels` maps those back to MetaWrap2 sample names.
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
                    "a shallow or subsampled file can have none. Sequence deeper or use the full "
                    "read file" % sample
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
            # appraise writes its summary table to stdout; the runner captures it to a file, so no
            # shell redirect is needed.
            stdout_path=summary,
        )
        if rc != 0 or not os.path.isfile(summary):
            error("singlem appraise failed; see its output above for the cause")
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
