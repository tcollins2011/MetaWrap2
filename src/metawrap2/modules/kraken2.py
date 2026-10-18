"""MetaWrap2 kraken2 module: taxonomic profiling of reads and/or contigs with KRAKEN2.

Runs KRAKEN2 on any number of assembly FASTA files and/or paired-end read sets, translates
the taxid output into full taxonomy lineages, collapses each into KRONA format, and draws a
single interactive kronagram of all samples with KronaTools.

With ``--bracken`` it then turns each read sample's KRAKEN2 report into an abundance estimate at
one rank. That is a separate step rather than the default because the KRAKEN2 report's read counts
are not abundances: a read whose k-mers are shared across a genus is assigned to the genus node, so
a species present at 10% of the community can show a small share of reads at its own node with the
rest stranded above it. Bracken redistributes those reads downward using the k-mer distribution of
the KRAKEN2 database, so it needs a database built from *that* KRAKEN2 database at a matching read
length. See :mod:`metawrap2.bracken` for the method, the column layout and the guards.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers/defaults in the CONSTANTS block. Change a flag there
and it takes effect — no need to follow the orchestration logic underneath.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import glob
import os
import random
from typing import Dict, List, Sequence, Tuple

from .. import bracken as _bracken
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..io.seqio import smart_open
from ..progress import bar
from ..scripts import kraken2_translate, kraken_to_krona
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
    resolve_threads,
    run,
    start_run,
    threads_arg,
    warning,
)

CONDA_ENV = MODULE_ENVS["kraken2"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
KRAKEN2_PAIRED = (
    "kraken2 --use-names --db {db} --paired --threads {threads} --output {out} {r1} {r2}"
)
KRAKEN2_SINGLE = "kraken2 --use-names --db {db} --threads {threads} --output {out} {seq}"
MEMORY_MAPPING = "--memory-mapping"  # appended when --no-preload is given (lower memory, slower)
# Appended only under --bracken. Bracken takes the aggregated *report*, not the per-read --output
# file, and refuses the latter outright: its check_report_file() exits 1 when the first character is
# C or U, which is exactly what a --output file starts with. Kept off by default so that a run
# without --bracken writes byte-identical output to one from before this flag existed.
#
# Deliberately not --report-minimizer-data: that inserts two columns before the taxid, and Bracken
# reads the rank and taxid by negative index, so the report would parse into the wrong fields.
KRAKEN2_REPORT = "--report {report}"
# -w is the second output: a Kraken2-style report with Bracken's estimates in place of the assigned
# counts, which is what KrakenTools and Pavian consume. -t is a read floor applied *before*
# redistribution, not after; see BRACKEN_THRESHOLD.
BRACKEN = (
    "bracken -d {db} -i {report} -o {out} -w {outreport} "
    "-r {read_length} -l {rank} -t {threshold}"
)
KTIMPORTTEXT = "ktImportText -o {out} {krona_files}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
FASTA_EXTS = (".fa", ".fasta", ".fa.gz", ".fasta.gz")
#: The first file of a pair. Shared by the Kraken2 loop and the Bracken stage so the two agree on
#: which inputs are reads - Bracken runs on those and only those.
READ_1_SUFFIXES = ("_1.fastq", "_1.fq", "_1.fastq.gz", "_1.fq.gz")
#: Kraken2's aggregated report, written per sample only under --bracken.
REPORT_SUFFIX = ".kreport"
#: Bracken's own abundance table, per sample.
BRACKEN_SUFFIX = ".bracken"
#: One long-format abundance table across every sample, in the four-column shape
#: `catalogue --quantify` writes and `diff_abundance` reads.
BRACKEN_TABLE = "bracken_abundance" + TSV_SUFFIX
#: Default read length. 100 because that is Bracken's own default and because the pre-built
#: distributions shipped with the Kraken2 standard database cover 50/75/100/150/200/250/300 - but it
#: must match the reads, so there is no safe default and the guard reports what the database has.
BRACKEN_READ_LENGTH = 100
#: Minimum reads at a node before Bracken will redistribute to it. 10 is Bracken's wrapper default
#: (its own -h text claims 0 and is wrong). Raising it drops rare taxa rather than reassigning them,
#: so their reads leave the table entirely instead of landing somewhere else.
BRACKEN_THRESHOLD = 10
#: Which column of Bracken's table goes into the long abundance table. The fraction, because it is
#: the compositional quantity `diff_abundance`'s ANCOM-BC expects; `_bracken.READS_UNIT` writes the
#: re-estimated counts instead, which are not comparable across libraries of different depth.
BRACKEN_UNIT = _bracken.FRACTION_UNIT
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 kraken2",
        usage="metawrap2 kraken2 [options] -o output_dir assembly.fasta reads_1.fastq reads_2.fastq ...",
        description="Taxonomic profiling of any number of assembly FASTA files (*.fa/*.fasta) "
        "and/or paired-end reads (*_1.fastq/*_2.fastq, .gz accepted).",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "-s",
        "--subsample",
        type=int,
        default=None,
        help="read subsampling number (default: all reads)",
    )
    p.add_argument(
        "--no-preload",
        action="store_true",
        help="do not pre-load the kraken2 DB into memory (slower, lower memory)",
    )
    p.add_argument(
        "--bracken",
        action="store_true",
        help="re-estimate abundance at one rank with Bracken (Lu et al., PeerJ Comput Sci 2017). "
        "kraken2 reports the reads assigned *at* each node, which is not an abundance: a read "
        "matching sequence conserved across a genus stops at the genus node, so a species at 10%% "
        "of the community can show a fraction of that at its own node. Bracken redistributes those "
        "reads down to --bracken-rank, and needs a database built from this kraken2 database at a "
        "matching read length. Runs on read samples only, not on assemblies",
    )
    p.add_argument(
        "--bracken-db",
        help="directory holding the Bracken k-mer distributions (default: the kraken2 database "
        "directory, which is where bracken-build writes them)",
    )
    p.add_argument(
        "--bracken-read-length",
        type=int,
        default=BRACKEN_READ_LENGTH,
        metavar="BP",
        help="read length the Bracken database was built for (default: %(default)s). Must match a "
        "databaseNmers.kmer_distrib in --bracken-db: the distribution depends on how many k-mers a "
        "read has, so a mismatch is not an approximation Bracken can make",
    )
    p.add_argument(
        "--bracken-rank",
        choices=list(_bracken.RANKS),
        default="S",
        help="rank to estimate abundance at (default: %(default)s for species)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument("seqs", nargs="+", help="assembly FASTA files and/or paired-end read files")
    return p.parse_args(argv)


def _is_fasta(path: str) -> bool:
    return path.endswith(FASTA_EXTS)


def _fasta_sample(path: str) -> str:
    base = os.path.basename(path)
    base = base.removesuffix(".gz")
    return os.path.splitext(base)[0]


def _fastq_sample(reads_1: str) -> str:
    base = os.path.basename(reads_1)
    for token in ("_1.fastq", "_1.fq"):
        if token in base:
            return base.split(token)[0]
    return os.path.splitext(base)[0]


def _iter_records(path):
    """Yield 4-line FASTQ records as single strings (gz/bz2-transparent)."""
    with smart_open(path, "rt") as fh:
        block = []
        for line in fh:
            block.append(line if line.endswith("\n") else line + "\n")
            if len(block) == 4:
                yield "".join(block)
                block = []
    if block:
        raise ValueError(
            "%s ends in a truncated FASTQ record (%d trailing lines)" % (path, len(block))
        )


def _subsample_pairs(r1: str, r2: str, depth: int, out1: str, out2: str) -> None:
    """Randomly subsample *depth* read pairs from r1/r2 into out1/out2 (gz-transparent input).

    Uses reservoir sampling, so peak memory is bounded by *depth* rather than by the input
    size. The previous implementation materialised every record of both files in Python
    lists first, which needs many times the file size in RAM and made ``--subsample``
    unusable on exactly the large libraries it exists to shrink.
    """
    reservoir = []
    pairs = zip(_iter_records(r1), _iter_records(r2))
    for i, pair in enumerate(bar(pairs, desc="subsampling read pairs", unit=" pair")):
        if i < depth:
            reservoir.append(pair)
        else:
            j = random.randint(0, i)
            if j < depth:
                reservoir[j] = pair
    with open(out1, "w") as o1, open(out2, "w") as o2:
        for a, b in reservoir:
            o1.write(a)
            o2.write(b)


def _kraken2(cmd: str, args, sample: str) -> str:
    if args.no_preload:
        cmd = cmd + " " + MEMORY_MAPPING
    if args.bracken:
        cmd = cmd + " " + KRAKEN2_REPORT.format(report=_report_path(args, sample))
    return cmd


def _report_path(args, sample: str) -> str:
    return os.path.join(args.output, sample + REPORT_SUFFIX)


def _read_samples(args) -> List[Tuple[str, str]]:
    """(sample, first-mate path) for the read inputs, in the order given.

    Derived from ``args.seqs`` rather than returned by the Kraken2 loop so that it is the same list
    under ``--resume``, when that loop does not run at all. Reads only: Bracken's model is the k-mer
    content of a read of a stated length, and a contig is neither, so an assembly has no entry here.
    """
    return [(_fastq_sample(path), path) for path in args.seqs if path.endswith(READ_1_SUFFIXES)]


def _run_kraken_on_inputs(args, env, db) -> None:
    for num in args.seqs:
        if num.endswith(READ_1_SUFFIXES):
            reads_1 = num
            reads_2 = num.replace("_1.fastq", "_2.fastq").replace("_1.fq", "_2.fq")
            if not os.path.isfile(reads_2):
                error("%s does not exist. Exiting..." % reads_2)
            sample = _fastq_sample(reads_1)
            comm("Now processing %s and %s with %d threads" % (reads_1, reads_2, args.threads))

            if args.subsample is not None and not dry_run():
                comm("subsampling down to %d reads..." % args.subsample)
                tmp_1 = os.path.join(args.output, "tmp_1.fastq")
                tmp_2 = os.path.join(args.output, "tmp_2.fastq")
                _subsample_pairs(reads_1, reads_2, args.subsample, tmp_1, tmp_2)
                reads_1, reads_2 = tmp_1, tmp_2
                comm("Subsampling done. Starting KRAKEN...")
                if not os.path.getsize(reads_1):
                    error("something went wrong with subsampling sequences. Exiting...")

            krak2 = os.path.join(args.output, "%s.krak2" % sample)
            if not dry_run() and os.path.isfile(krak2) and os.path.getsize(krak2):
                comm("%s already exists - skipping running kraken2" % krak2)
            else:
                cmd = _kraken2(
                    KRAKEN2_PAIRED.format(
                        db=db, threads=args.threads, out=krak2, r1=reads_1, r2=reads_2
                    ),
                    args,
                    sample,
                )
                run(cmd, env=env, tool="kraken2")
            if not dry_run() and not (os.path.isfile(krak2) and os.path.getsize(krak2)):
                error(
                    "Something went wrong with running kraken2 on %s and %s . Exiting..."
                    % (reads_1, reads_2)
                )
            if args.subsample is not None and not dry_run():
                os.remove(os.path.join(args.output, "tmp_1.fastq"))
                os.remove(os.path.join(args.output, "tmp_2.fastq"))

        elif _is_fasta(num):
            sample = _fasta_sample(num)
            comm("Now processing %s with %d threads" % (num, args.threads))
            krak2 = os.path.join(args.output, "%s.krak2" % sample)
            if not dry_run() and os.path.isfile(krak2) and os.path.getsize(krak2):
                comm("%s already exists - skipping running kraken2" % krak2)
            else:
                cmd = _kraken2(
                    KRAKEN2_SINGLE.format(db=db, threads=args.threads, out=krak2, seq=num),
                    args,
                    sample,
                )
                run(cmd, env=env, tool="kraken2")
            if not dry_run() and not (os.path.isfile(krak2) and os.path.getsize(krak2)):
                error("Something went wrong with running kraken2 on %s. Exiting..." % num)


# ─── Bracken ─────────────────────────────────────────────────────────────────────────────


def _resolve_bracken_db(args, kraken2_db: str) -> str:
    """Where the Bracken k-mer distributions are, and a usable error when they are not there.

    Read from ``--bracken-db``, defaulting to the Kraken2 database directory because that is where
    ``bracken-build`` writes ``databaseNmers.kmer_distrib`` - the files belong to that database and
    are meaningless beside another one.

    There is deliberately no ``BRACKEN_DB`` key in ``metawrap2.toml`` yet, and so no `doctor` check
    and no `install-db` entry: those live in files this change does not own. Until they land the flag
    is the only route, so the failure has to say so rather than report a missing config key the user
    cannot set.
    """
    db = args.bracken_db or kraken2_db
    problem = _bracken.check_read_length(db, args.bracken_read_length)
    if problem:
        message = (
            "--bracken cannot run: %s. Point --bracken-db at a directory holding the "
            "distributions, or build them into the Kraken2 database. MetaWrap2 has no "
            "`install-db bracken` yet, so this is not a database it can fetch for you." % problem
        )
        # A dry run is asked to print the commands it would run, and it is routinely used on a
        # machine that has neither the database nor the reads. Erroring here would make the flag
        # impossible to inspect before committing to a real run.
        if dry_run():
            warning(message)
        else:
            error(message)
    return db


def _run_bracken(args, env, db: str, samples: Sequence[Tuple[str, str]]) -> None:
    """Estimate abundance at one rank for every read sample, and write the combined long table."""
    per_sample: Dict[str, List[_bracken.BrackenRow]] = {}
    rank = args.bracken_rank
    for sample, reads_1 in samples:
        report = _report_path(args, sample)
        if not dry_run() and not (os.path.isfile(report) and os.path.getsize(report)):
            error(
                "%s was not written, so Bracken has nothing to read for %s. kraken2 only writes a "
                "report when --bracken is given; a directory from an earlier run without the flag "
                "has none, and needs --force rather than --resume" % (report, reads_1)
            )

        table = os.path.join(args.output, sample + BRACKEN_SUFFIX)
        outreport = os.path.join(
            args.output, "%s_bracken_%s%s" % (sample, _bracken.RANK_NAMES[rank], REPORT_SUFFIX)
        )
        comm("Estimating %s abundance for %s" % (_bracken.RANK_NAMES[rank], sample))
        run(
            BRACKEN.format(
                db=db,
                report=report,
                out=table,
                outreport=outreport,
                read_length=args.bracken_read_length,
                rank=rank,
                threshold=BRACKEN_THRESHOLD,
            ),
            env=env,
            tool="bracken",
            # Bracken's exit status carries no information: its wrapper ends every error path with
            # `echo "ERROR: ..."` and a bare `exit`, which returns the echo's status of 0. A missing
            # database, a missing kmer_distrib and a missing input all look like success, so the
            # table below is the only evidence that anything happened.
            check=False,
        )
        if dry_run():
            continue
        if not (os.path.isfile(table) and os.path.getsize(table)):
            error(
                "Bracken wrote no table for %s (expected %s). It exits 0 even when it fails, so "
                "its message is in run.stderr rather than in a non-zero status" % (sample, table)
            )

        rows = _bracken.read_table(table)
        _bracken.check_table_rank(rows, rank, table)
        summary = _bracken.redistribution(rows, rank, _bracken.reads_above_rank(report, rank))
        comm(summary.describe())
        suspicious = _bracken.check_redistribution(summary)
        if suspicious:
            warning("%s: %s" % (sample, suspicious))
        per_sample[sample] = rows

    if dry_run():
        return
    if not per_sample:
        warning("no read sample produced a Bracken table; no abundance table written")
        return

    combined = os.path.join(args.output, BRACKEN_TABLE)
    written = _bracken.write_long_table(combined, per_sample, BRACKEN_UNIT)
    comm(
        "%d row(s) of %s abundance (%s) across %d sample(s) -> %s"
        % (written, _bracken.RANK_NAMES[rank], BRACKEN_UNIT, len(per_sample), combined)
    )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    env = env_for("kraken2", settings)
    rec = start_run("kraken2", args, env, settings, inputs=list(args.seqs))

    ckpt = make_checkpoint(args.output)
    try:
        db = settings.db("KRAKEN2_DB")
        if not (db and os.path.isdir(db)):
            error(
                "The folder %s doesnt exist. Please consult the MetaWrap2 database guide to "
                "download and build the KRAKEN2 database" % db
            )

        announcement("RUNNING KRAKEN ON ALL FILES")
        if os.path.isdir(args.output):
            warning("%s already exists." % args.output)
        ensure_dir(args.output)

        if ckpt.todo("kraken2"):
            _run_kraken_on_inputs(args, env, db)
            ckpt.done("kraken2")
        else:
            comm("skipping kraken2 on inputs (already done; --resume)")

        krak2_files = sorted(glob.glob(os.path.join(args.output, "*.krak2")))
        if not krak2_files and not dry_run():
            error("No fasta or fastq files detected! (must be in .fastq .fa .fasta or .fq format)")

        if args.bracken:
            read_samples = _read_samples(args)
            skipped = len(args.seqs) - len(read_samples)
            if skipped:
                warning(
                    "--bracken skips %d assembly input(s). Bracken estimates how many reads of a "
                    "stated length a genome would contribute; a contig is neither a read nor of "
                    "that length, so running it on an assembly produces numbers that look like "
                    "abundances and are not" % skipped
                )
            if not read_samples:
                warning("--bracken was given but no paired-end read input was; nothing to estimate")
            else:
                bracken_db = _resolve_bracken_db(args, db)
                # On the record for the same reason every other database path is: an abundance
                # estimate is only reproducible against the k-mer distribution that produced it, and
                # this one comes from a flag rather than from metawrap2.toml, so nothing else in the
                # run would name it.
                rec.databases["BRACKEN_DB"] = bracken_db
                if ckpt.todo("bracken"):
                    announcement("ESTIMATING ABUNDANCE WITH BRACKEN")
                    _run_bracken(args, env, bracken_db, read_samples)
                    ckpt.done("bracken")
                else:
                    comm("skipping bracken (already done; --resume)")

        if ckpt.todo("translate") and not dry_run():
            announcement("RUNNING KRAKEN-TRANSLATE ON OUTPUT")
            for krak2 in krak2_files:
                comm("Translating %s" % krak2)
                translated = os.path.splitext(krak2)[0] + ".kraken2"
                # Fall back to the configured NCBI taxdump if the Kraken2 database does not
                # carry names.dmp/nodes.dmp itself.
                kraken2_translate.translate_kraken2_annotations(
                    annotation_file=krak2,
                    kraken2_db=db,
                    output=translated,
                    taxdump_dirs=(settings.db("TAXDUMP"),),
                )
            ckpt.done("translate")
        else:
            comm("skipping kraken-translate (already done; --resume)")

        if ckpt.todo("krona"):
            announcement("MAKING KRONAGRAM OF ALL FILES")
            if not dry_run():
                for translated in sorted(glob.glob(os.path.join(args.output, "*.kraken2"))):
                    krona = os.path.splitext(translated)[0] + ".krona"
                    with open(krona, "w") as out:
                        kraken_to_krona.to_krona(translated, out)
                    if not os.path.getsize(krona):
                        error(
                            "Something went wrong with making krona file from kraken file. "
                            "Exiting..."
                        )

            krona_files = sorted(glob.glob(os.path.join(args.output, "*.krona")))
            if not krona_files and dry_run():
                krona_files = [os.path.join(args.output, "<sample>.krona")]
            kronagram = os.path.join(args.output, "kronagram.html")
            run(
                KTIMPORTTEXT.format(out=kronagram, krona_files=" ".join(krona_files)),
                env=env,
                tool="ktImportText",
            )
            if not dry_run() and not (os.path.isfile(kronagram) and os.path.getsize(kronagram)):
                error("Something went wrong with running KronaTools to make kronagram. Exiting...")
            ckpt.done("krona")
        else:
            comm("skipping kronagram (already done; --resume)")

        announcement("FINISHED RUNNING KRAKEN2 PIPELINE!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
