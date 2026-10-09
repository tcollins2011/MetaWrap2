"""MetaWrap2 quant_bins module: estimate the abundance of each bin in each sample.

Four quantifiers answer that question four ways, and all of them write the same bins x samples
``bin_abundance_table.tsv`` plus, for more than one sample, a clustered Seaborn heatmap of it.

**sylph is the default**, and the reason is measured rather than assumed. On 76 bins from the CAMI2
gastrointestinal gold-standard assembly, quantified across three samples against a read-level truth
(``metawrap2 benchmark run cami2_gi quantifier``):

    quantifier   Spearman  Pearson   L1 sim   CPU-hours  peak RSS
    sylph          0.868    0.811    0.722      0.126     1.3 GB
    salmon         0.827    0.959    0.835      1.078     7.4 GB
    CoverM         0.803    0.952    0.791      2.368    11.6 GB
    depth          0.780    0.953    0.799      0.013     small

**Spearman is an ordering and Pearson is a magnitude, and no quantifier here wins both.** sylph orders
the bins best, for a tenth of salmon's cpu; salmon gets the magnitudes closest, and sylph is worst at
them by a wide margin (0.811 against 0.959). So the rule of thumb is that sylph is right when the
question is which bins are most abundant and how they change across samples, and
``--quantifier salmon`` earns its keep when the actual numbers are the point - a ratio between two
bins, say, rather than an ordering.

``depth`` is the interesting fourth. It aligns nothing: it aggregates the per-contig depth table the
binning run already wrote into a length-weighted mean per bin, so it costs **0.013 cpu-hours** - 10x
less than sylph, 80x less than salmon, 180x less than CoverM - and it matches salmon and CoverM on
magnitude (Pearson 0.953, L1 0.799, ahead of CoverM on both). What it is worst at is the ordering,
which is the metric to lean on, and that is the trade to understand before reaching for it: Pearson is
dominated by the abundant bins, Spearman weights every bin alike, so a method that is slightly noisy on
low-coverage bins reads well on one and badly on the other. Reach for it when the magnitudes matter and
the table is already on disk; not when the question is which bins rank where.

How they differ mechanically: salmon quantifies the reads against a pseudo-transcriptome of the bins'
contigs and collapses per-contig counts to a length-weighted median coverage per bin. CoverM works in
genome mode, mapping against the bins themselves and reporting per-genome coverage directly, with no
index to build and no counts to collapse. sylph skips alignment altogether - it sketches both the
bins and the reads into k-mer sketches and profiles one against the other, which is why its cost
barely grows with sequencing depth and why it needs 1.3 GB where CoverM needs 11.6. ``depth`` reads no
reads at all: it is one pass over a TSV that another module already paid for, which is why its cost is
a rounding error and why it is the only one whose numbers are on the same scale the binner clustered
on.

Everything downstream of the table is unchanged whichever ran, but **the numbers are not
interchangeable**, which is why this is a flag and not a silent substitution: salmon's cell is a
length-weighted median contig coverage, CoverM's default is relative abundance as a percentage of
the sample, and sylph's is a k-mer-derived sequence abundance, also a percentage. Rows of a CoverM
or sylph table sum to at most 100; rows of a salmon table have no such meaning.

The two percentage tables are still not the same table. A genome sylph does not detect is absent
from its output rather than present with a zero, because a sketch below the ANI threshold is no
evidence of that genome at all - so a zero in a sylph column means "not detected above
``SYLPH_MIN_ANI``", where a zero in a CoverM column means "mapped, but under
``COVERM_MIN_COVERED_FRACTION`` of the genome was covered".

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import os
import shutil
from typing import List, TextIO

from ..cliformat import make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import FASTA_EXTENSIONS_WITH_GZ, TSV_SUFFIX, bin_label
from ..io.seqio import contig_id, fasta_files, iter_fasta
from ..progress import bar
from ..pyrun import py_module
from ..scripts import split_salmon_out_into_bins, summarize_salmon_files
from ._common import (
    absolutize_paths,
    announcement,
    check_fasta,
    check_fastq,
    collect_read_pairs,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["quant_bins"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
SALMON_INDEX = "salmon index -p {threads} -t {assembly} -i {index}"
SALMON_QUANT = "salmon quant -i {index} --libType IU -1 {r1} -2 {r2} -o {out} --meta -p {threads}"
# MetaWrap2's own helper -> host interpreter (env=None); the quant_bins env has salmon,
# not python+seaborn. See metawrap2.pyrun.
MAKE_HEATMAP = py_module("metawrap2.scripts.make_heatmap", "{table}", "{png}")
# One invocation for every sample: CoverM takes repeated --coupled pairs and writes one
# genome x sample table, so there is no per-sample loop and nothing to summarise afterwards.
#
# The bins are named in a list file rather than with --genome-fasta-directory/-x. CoverM's -x takes
# one extension and its default is `fna`, so pointing it at a directory means guessing which single
# extension the folder uses - and a wrong guess yields an empty genome set rather than an error,
# while a folder holding both .fasta and .fasta.gz bins cannot be expressed at all. Naming the files
# says exactly which bins are meant, and lets a stray .stats.tsv or checkm/ subdirectory sit in a bin
# folder without becoming CoverM's problem.
COVERM_GENOME = (
    "coverm genome --genome-fasta-list {genome_list} "
    "--coupled {coupled} --methods {methods} --min-covered-fraction {min_covered} "
    "--mapper {mapper} --threads {threads} --output-file {out}"
)
# sylph works in two stages, and both of them are sketching: once over the bins (producing one
# .syldb) and once over the reads (one .sylsp per sample), after which `profile` compares sketch to
# sketch and never looks at a read again. The genome sketch takes the same list file CoverM is given.
SYLPH_SKETCH_GENOMES = "sylph sketch --gl {genome_list} -o {db} -t {threads}"
# --l1/--l2 are lists of forward and reverse read files, matched by position, so every sample is
# sketched in one invocation. The sketches land in -d as <forward read filename>.paired.sylsp.
SYLPH_SKETCH_READS = "sylph sketch --l1 {l1} --l2 {l2} -d {sketch_dir} -t {threads}"
SYLPH_PROFILE = (
    "sylph profile {db} {sketches} --minimum-ani {min_ani} -t {threads} --output-file {out}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: CoverM's genome-mode default. Kept as the one method written to the abundance table, because
#: that table has one cell per bin per sample and a second method would have to overwrite it.
COVERM_METHOD = "relative_abundance"
#: A genome must have this percentage of its length covered before CoverM reports anything but
#: zero for it. CoverM's own default is 10%, and the reason to keep it is that a handful of reads
#: landing on a conserved gene is not evidence that a genome is present.
COVERM_MIN_COVERED_FRACTION = 10
#: CoverM's own default, pinned here rather than left implicit: it changed between CoverM releases
#: (0.8 defaults to strobealign where earlier versions did not), and a mapper swapping itself out
#: underneath a user is a change in the numbers, not just in the runtime. `coverm genome --full-help`
#: lists the alternatives - ``minimap2-sr``, ``bwa-mem``, ``bwa-mem2`` - and which one earns the
#: default here is a question for the quantifier benchmark, not for taste.
COVERM_MAPPER = "strobealign"
#: A bin whose sketch matches the sample below this ANI is reported as absent rather than as rare.
#: This is the single most consequential knob in the module, and the number is measured rather than
#: inherited. On 76 CAMI2 GI bins across three samples, sylph's own default of 95 sent 55% of the
#: abundance it reported in one sample to four bins whose genomes are absent from that sample
#: entirely - relatives close enough to match at 95% ANI, which is a species-level threshold applied
#: to a strain-level question. Raising the floor fixes it, up to a point:
#:
#:     ANI    Spearman  Pearson    L1   misattributed mass
#:      95      0.758    0.678   0.832       55%
#:      97      0.823    0.744   0.728       44%
#:      99      0.868    0.811   0.555       25%      <- default
#:      99.5    0.857    0.800   0.536       25%
#:      99.9    0.696    0.840   0.448       30%   (13 of 29 present bins missed entirely)
#:
#: 99 is the turning point: past it the threshold starts discarding bins that really are there, and
#: rank accuracy collapses even as the surviving estimates get sharper. Below 95 sylph itself warns
#: "Smaller than 95 for profile will give inaccurate results", its containment statistics being
#: calibrated for same-species comparison, so 95 is a floor rather than a tuning range.
#:
#: These bins came off the pooled gold-standard assembly and are therefore near-complete and
#: low-error; bins from a fragmented real assembly match their own genome less tightly, so a strict
#: floor could cost more there than it did here. That is why this is a flag (``--min-ani``) and not
#: just a constant.
SYLPH_MIN_ANI = 99
#: Below this sylph's own containment statistics stop being calibrated, so the flag warns rather than
#: silently producing numbers that do not mean what the column header says.
SYLPH_ANI_FLOOR = 95
#: Which of sylph's two abundance columns fills the table. ``Sequence_abundance`` is the fraction of
#: the sample's *sequence* a bin accounts for, which is what salmon's coverage and CoverM's relative
#: abundance both measure, and what the CAMI2 read-to-genome truth counts. ``Taxonomic_abundance``
#: is the fraction of *organisms*, corrected for genome length - a different and equally valid
#: number, but not one the other quantifiers or the benchmark truth can be compared against.
SYLPH_ABUNDANCE_COLUMN = "Sequence_abundance"
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 quant_bins",
        usage="metawrap2 quant_bins [options] -b bins_folder -o output_dir -a assembly.fa "
        "readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]",
        description="Estimate the abundance of each bin across samples with salmon, and "
        "draw a clustered abundance heatmap (.gz reads accepted).",
    )
    p.add_argument(
        "-b", "--bins", required=True, help="folder containing draft genomes (bins) in fasta format"
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-a",
        "--assembly",
        default="",
        help="fasta file with entire metagenomic assembly (strongly recommended!)",
    )
    p.add_argument(
        "--depth",
        help="contig depth table to aggregate, for --quantifier depth. This is the table binning or "
        "coverage already wrote (work_files/contig_depth.txt, or coverage's per-sample slice)",
    )
    p.add_argument(
        "--quantifier",
        choices=tuple(QUANTIFIERS),
        default="sylph",
        help="sylph (default): sequence abundance from k-mer sketches, no alignment at all. salmon: "
        "length-weighted median contig coverage per bin. coverm: per-genome relative abundance, "
        "mapping against the bins directly. The numbers are not interchangeable - see the "
        "module docstring",
    )
    p.add_argument(
        "--min-ani",
        type=percent(),
        default=SYLPH_MIN_ANI,
        metavar="PERCENT",
        help="sylph only: report a bin as absent unless its sketch matches the sample at this "
        "average nucleotide identity (default: %(default)s, measured; see SYLPH_MIN_ANI)",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    # Not required: --quantifier depth reads a table that already encodes the reads, so demanding them
    # would be asking for files it will not open.
    p.add_argument(
        "reads",
        nargs="*",
        help="paired read files (*_1.fastq/*_2.fastq, .gz accepted); not needed with "
        "--quantifier depth, which aggregates an existing depth table",
    )
    return p.parse_args(argv)


def _write_genome_list(bin_folder: str, path: str) -> List[str]:
    """Write the bins to quantify, one absolute path per line; return the paths written.

    Every FASTA in the folder is named explicitly, compressed or not, so the tool needs no extension
    to guess at and a folder mixing ``.fasta`` with ``.fasta.gz`` bins still quantifies in one pass.
    Anything that is not a FASTA is skipped rather than passed along: a bin folder that has been
    through another module also holds ``.stats.tsv`` files and tool subdirectories, and handing those
    over as genomes is an error at best and a bogus row at worst.

    CoverM reads this as ``--genome-fasta-list`` and sylph as ``--gl``; they want the same thing, so
    they are given the same file. The returned paths are what the sylph table is pivoted against,
    since sylph reports only the bins it detected and the missing rows have to be filled in from
    somewhere.

    Sorted, so the genome list - and anything downstream that depends on the order the genomes were
    given in - is the same on every run over the same folder.
    """
    bins = fasta_files(bin_folder)
    with open(path, "w") as fh:
        fh.writelines("%s\n" % b for b in bins)
    return bins


def _coverm_table(raw: str, samples: List[str], out: TextIO) -> int:
    """Rewrite CoverM's genome table as the ``Genomic bins`` table, returning the row count.

    Two deliberate translations, both so that a table from either quantifier can be read by the same
    eyes and the same scripts:

    * **Row labels go through :func:`bin_label`**, the same function the salmon path labels its rows
      with, so the two tables key a bin identically by construction rather than by coincidence.
      CoverM already strips the extension it was given in ``-x``; this is what makes that agreement
      a guarantee instead of an observation that happened to hold on one set of bins.
    * **Columns are matched by position, not by title.** CoverM titles each column
      ``<file> Relative Abundance (%)``, derived from the read filename, and emits them in the order
      the ``--coupled`` pairs were given - which is the order of *samples*. Using our own sample
      labels keeps the two tables' headers identical instead of having one of them say
      ``sample_1_1.fastq``.

    CoverM's ``unmapped`` row is dropped: it is not a bin, and leaving it in would put a row in a
    bins x samples table that no bin folder contains. Rows come out sorted by bin, which CoverM's
    own output is not - it emitted ``bin_2`` before ``bin_1`` here, in neither file nor alphabetical
    order - so the ordering of this table is ours and stays put whatever CoverM does internally.
    """
    rows = []
    with open(raw) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        values_per_row = len(header) - 1
        if values_per_row != len(samples):
            raise RuntimeError(
                "CoverM reported %d value column(s) for %d sample(s); cannot match them up. "
                "Its table is at %s" % (values_per_row, len(samples), raw)
            )
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0] or fields[0] == "unmapped":
                continue
            rows.append((bin_label(fields[0]), fields[1:]))
    out.write("Genomic bins\t%s\n" % "\t".join(samples))
    out.writelines("%s\t%s\n" % (name, "\t".join(values)) for name, values in sorted(rows))
    return len(rows)


def _sylph_table(raw: str, samples: List[tuple], bins: List[str], out: TextIO) -> int:
    """Pivot sylph's long profile into the ``Genomic bins`` table, returning the row count.

    Three things have to be undone here, and each of them would be a silently wrong table if it were
    not:

    * **sylph's output is long, not a matrix.** One row per (sample, genome) *detection*, so it has
      to be pivoted into the bins x samples shape the rest of MetaWrap2 reads.
    * **Undetected genomes are absent, not zero.** sylph omits a genome whose sketch matched below
      ``SYLPH_MIN_ANI`` entirely, so a bin present in one sample and not another has one row, not
      two. Pivoting only what is in the file would give a ragged table and drop a bin nothing
      detected from the output altogether; the bins are therefore taken from the genome list, which
      is every bin the tool was given, and anything unreported is filled in with zero.
    * **sylph names samples and genomes by file path.** ``Sample_file`` is the *forward read file*
      the sketch came from - not the ``.sylsp`` name - which is what makes it a reliable key back to
      our own sample labels; ``Genome_file`` goes through :func:`bin_label`, the same function the
      other two quantifiers label their rows with, so the three tables key a bin identically.

    Columns are found by name rather than by position, unlike the CoverM table where the value
    columns are titled after read filenames and only their order is meaningful. sylph's header is
    made of fixed names, so a future column inserted before ``Sequence_abundance`` should not shift
    the numbers.
    """
    labels = [sample for sample, _, _ in samples]
    #: The forward read path is the key, with its basename as a fallback in case sylph ever
    #: normalises the path it was handed.
    by_read = {}
    for sample, r1, _r2 in samples:
        by_read[r1] = sample
        by_read.setdefault(os.path.basename(r1), sample)

    rows: dict = {bin_label(os.path.basename(b)): {} for b in bins}
    with open(raw) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_sample = header.index("Sample_file")
            i_genome = header.index("Genome_file")
            i_abundance = header.index(SYLPH_ABUNDANCE_COLUMN)
        except ValueError:
            raise RuntimeError(
                "sylph's profile at %s has no %s column; its header is %s"
                % (raw, SYLPH_ABUNDANCE_COLUMN, "/".join(header) or "(empty)")
            ) from None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(i_sample, i_genome, i_abundance):
                continue
            read = fields[i_sample]
            sample = by_read.get(read) or by_read.get(os.path.basename(read))
            if sample is None:
                raise RuntimeError(
                    "sylph reported sample %s, which is none of the read files it was given. Its "
                    "profile is at %s" % (read, raw)
                )
            rows.setdefault(bin_label(os.path.basename(fields[i_genome])), {})[sample] = fields[
                i_abundance
            ]
    out.write("Genomic bins\t%s\n" % "\t".join(labels))
    out.writelines(
        "%s\t%s\n" % (name, "\t".join(by_sample.get(label, "0.0") for label in labels))
        for name, by_sample in sorted(rows.items())
    )
    return len(rows)


def _concat_bins_into_assembly(bin_folder: str, assembly: str) -> None:
    """Concatenate every bin's contigs into a single assembly fasta (supports .gz bins)."""
    if os.path.isfile(assembly):
        os.remove(assembly)
    with open(assembly, "w") as out:
        for f in fasta_files(bin_folder):
            out.writelines(">%s\n%s\n" % (contig_id(header), seq) for header, seq in iter_fasta(f))


def _heatmap(table: str, output: str, samples: int) -> None:
    """Draw the clustered abundance heatmap from *table*, or say why there is none.

    Shared by both quantifiers: the table has the same shape whichever tool filled it in, so the
    picture is drawn the same way and lives under the same name.
    """
    if samples < 2:
        warning("Cannot make clustered heatmap with just one sample... Skipping heatmap")
        return
    announcement("MAKING GENOME ABUNDANCE HEATMAP WITH SEABORN")
    comm("making heatmap with Seaborn")
    heatmap = os.path.join(output, "bin_abundance_heatmap.png")
    run(
        [a.format(table=table, png=heatmap) for a in MAKE_HEATMAP],
        env=None,
        tool="make_heatmap",
        hint="Something went wrong with making the heatmap.",
    )


def _quantify_salmon(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Index a pseudo-transcriptome, quantify every sample against it, collapse the result to bins."""
    announcement("SETTING UP OUTPUT AND INDEXING ASSEMBLY")
    assembly = args.assembly
    if not assembly:
        comm("Concatenating bins into a metagenomic assembly file.")
        assembly = os.path.join(args.output, "assembly.fa")
        if not dry_run():
            _concat_bins_into_assembly(args.bins, assembly)
    elif not os.path.isfile(assembly):
        error("Assembly file %s does not exist. Exiting..." % assembly)

    index = os.path.join(args.output, "assembly_index")
    if ckpt.todo("index"):
        comm("Indexing assembly file with salmon. Ignore any warnings")
        run(
            SALMON_INDEX.format(threads=args.threads, assembly=assembly, index=index),
            env=env,
            tool="salmon index",
            hint="Something went wrong with indexing the assembly.",
        )
        ckpt.done("index")
    else:
        comm("skipping salmon index (already done; --resume)")

    announcement("ALIGNING READS FROM ALL SAMPLES BACK TO BINS WITH SALMON")
    comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

    align_dir = os.path.join(args.output, "alignment_files")
    ensure_dir(align_dir)
    if ckpt.todo("quant"):
        for sample, r1, r2 in bar(
            pairs, desc="quantifying samples", unit=" sample", total=len(pairs)
        ):
            out_quant = os.path.join(align_dir, sample + ".quant")
            comm("processing sample %s with reads %s and %s..." % (sample, r1, r2))
            run(
                SALMON_QUANT.format(index=index, r1=r1, r2=r2, out=out_quant, threads=args.threads),
                env=env,
                tool="salmon quant",
                hint="Something went wrong with aligning %s fastq files back to assembly!" % sample,
            )
        ckpt.done("quant")
    else:
        comm("skipping salmon quant (already done; --resume)")

    quant_dir = os.path.join(args.output, "quant_files")
    if ckpt.todo("summarize") and not dry_run():
        comm("summarize salmon files...")
        summarize_salmon_files.summarize(align_dir)
        ensure_dir(quant_dir)
        for f in sorted(os.listdir(align_dir)):
            if f.endswith(".quant.counts"):
                os.replace(os.path.join(align_dir, f), os.path.join(quant_dir, f))

        announcement("EXTRACTING AVERAGE ABUNDANCE OF EACH BIN")
        n = len([f for f in os.listdir(quant_dir) if "counts" in f])
        if n < 1:
            error("There were no files found in %s" % quant_dir)
        comm("There were %d samples detected. Making abundance table!" % n)
        with open(table, "w") as out:
            split_salmon_out_into_bins.build_table(quant_dir, args.bins, assembly, out)
        comm("Average bin abundance table stored in %s" % table)

        _heatmap(table, args.output, n)

        # Clean up regardless of sample count: salmon's per-sample output is large and is
        # fully summarised into quant_files/ and the abundance table by this point. It used
        # to be deleted only on the >1-sample path, so single-sample runs silently kept it.
        comm("cleaning up salmon's per-sample alignment output...")
        shutil.rmtree(align_dir, ignore_errors=True)
        ckpt.done("summarize")
    elif dry_run():
        ensure_dir(quant_dir)
        comm("(dry run) would summarize salmon output into %s" % table)
        if len(pairs) > 1:
            _heatmap(table, args.output, len(pairs))
    else:
        comm("skipping abundance summary (already done; --resume)")


def _quantify_coverm(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Map every sample against the bins with CoverM and write the genome x sample table.

    No assembly and no index are involved: CoverM maps against the bin fastas themselves, so ``-a``
    has nothing to contribute here. It is also a single command for the whole run - every sample's
    pair goes on the same command line and the complete table comes back at once - so there is no
    per-sample checkpoint to resume from, only the one covering the whole invocation. Converting
    that table is cheap, so it is redone on every resume rather than checkpointed too.
    """
    announcement("QUANTIFYING BINS ACROSS SAMPLES WITH COVERM")
    if args.assembly:
        warning(
            "--quantifier coverm maps reads against the bins themselves; the assembly (-a) is unused"
        )
    comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

    genome_list = os.path.join(args.output, "coverm_genomes.txt")
    n_bins = len(_write_genome_list(args.bins, genome_list))
    if n_bins < 1:
        error("No bins to quantify in %s. Exiting..." % args.bins)
    comm("%d bins to quantify" % n_bins)

    raw = os.path.join(args.output, "coverm_genome" + TSV_SUFFIX)
    if ckpt.todo("coverm"):
        run(
            COVERM_GENOME.format(
                genome_list=genome_list,
                coupled=" ".join("%s %s" % (r1, r2) for _, r1, r2 in pairs),
                methods=COVERM_METHOD,
                min_covered=COVERM_MIN_COVERED_FRACTION,
                mapper=COVERM_MAPPER,
                threads=args.threads,
                out=raw,
            ),
            env=env,
            tool="coverm genome",
            hint="Something went wrong quantifying the bins with CoverM.",
        )
        ckpt.done("coverm")
    else:
        comm("skipping coverm genome (already done; --resume)")

    if dry_run():
        comm("(dry run) would rewrite %s as %s" % (raw, table))
        _heatmap(table, args.output, len(pairs))
        return

    announcement("EXTRACTING AVERAGE ABUNDANCE OF EACH BIN")
    with open(table, "w") as out:
        rows = _coverm_table(raw, [sample for sample, _, _ in pairs], out)
    if rows < 1:
        error(
            "CoverM reported no genomes, though %d bin(s) were given to it in %s. Its own table is "
            "at %s" % (n_bins, genome_list, raw)
        )
    comm("Abundance table for %d bins stored in %s" % (rows, table))
    _heatmap(table, args.output, len(pairs))


def _write_read_lists(pairs: list, output: str) -> tuple:
    """Write sylph's ``--l1``/``--l2`` read lists; return their paths.

    Two files rather than one because sylph pairs the lists up by position, which is also why they
    are written from the same ``pairs`` the other quantifiers use: whatever paired the reads up, all
    three quantifiers see the same pairing.
    """
    l1 = os.path.join(output, "sylph_reads_1.txt")
    l2 = os.path.join(output, "sylph_reads_2.txt")
    for path, index in ((l1, 1), (l2, 2)):
        with open(path, "w") as fh:
            fh.writelines("%s\n" % pair[index] for pair in pairs)
    return l1, l2


def _quantify_sylph(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Sketch the bins and the reads, profile one against the other, write the abundance table.

    Nothing is aligned and nothing is indexed in the mapping sense: both sides become k-mer sketches
    and ``profile`` compares sketches, which is what makes the cost of the profiling step independent
    of how deep the samples were sequenced. The two sketching steps are checkpointed separately
    because the genome sketch is the cheap one and the read sketches are not - a run interrupted
    after sketching 20 samples' reads should not sketch them again.

    Like the CoverM path, ``-a`` has nothing to contribute: the bins are the reference.
    """
    announcement("QUANTIFYING BINS ACROSS SAMPLES WITH SYLPH")
    if args.assembly:
        warning("--quantifier sylph profiles the bins themselves; the assembly (-a) is unused")
    min_ani = getattr(args, "min_ani", SYLPH_MIN_ANI)
    if min_ani < SYLPH_ANI_FLOOR:
        warning(
            "--min-ani %g is below %g, where sylph's own containment statistics stop being "
            "calibrated; it warns about this itself and the abundances stop meaning what the column "
            "header says" % (min_ani, SYLPH_ANI_FLOOR)
        )
    comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

    genome_list = os.path.join(args.output, "sylph_genomes.txt")
    bins = _write_genome_list(args.bins, genome_list)
    if not bins:
        error("No bins to quantify in %s. Exiting..." % args.bins)
    comm("%d bins to quantify" % len(bins))
    l1, l2 = _write_read_lists(pairs, args.output)

    # `sylph sketch -o bins` writes bins.syldb, so the suffix is sylph's to add, not ours.
    db = os.path.join(args.output, "sylph_bins")
    sketch_dir = os.path.join(args.output, "sylph_read_sketches")
    if ckpt.todo("sketch_genomes"):
        comm("sketching bins with sylph")
        run(
            SYLPH_SKETCH_GENOMES.format(genome_list=genome_list, db=db, threads=args.threads),
            env=env,
            tool="sylph sketch",
            hint="Something went wrong sketching the bins with sylph.",
        )
        ckpt.done("sketch_genomes")
    else:
        comm("skipping sylph genome sketch (already done; --resume)")

    if ckpt.todo("sketch_reads"):
        comm("sketching reads from %d sample(s) with sylph" % len(pairs))
        ensure_dir(sketch_dir)
        run(
            SYLPH_SKETCH_READS.format(l1=l1, l2=l2, sketch_dir=sketch_dir, threads=args.threads),
            env=env,
            tool="sylph sketch",
            hint="Something went wrong sketching the reads with sylph.",
        )
        ckpt.done("sketch_reads")
    else:
        comm("skipping sylph read sketches (already done; --resume)")

    raw = os.path.join(args.output, "sylph_profile" + TSV_SUFFIX)
    if ckpt.todo("profile"):
        # The sketch filenames are sylph's to choose (it appends .paired.sylsp to the forward read's
        # name), so they are collected from the directory rather than predicted. Counting them is
        # worth doing: a sample sylph declined to sketch would otherwise turn into a column of
        # zeroes, which reads as a sample with nothing in it.
        sketches = sorted(
            os.path.join(sketch_dir, name)
            for name in os.listdir(sketch_dir)
            if name.endswith(".sylsp")
        )
        if len(sketches) != len(pairs):
            error(
                "sylph sketched %d sample(s) from %d pair(s) of read files. Its sketches are in %s"
                % (len(sketches), len(pairs), sketch_dir)
            )
        run(
            SYLPH_PROFILE.format(
                db=db + ".syldb",
                sketches=" ".join(sketches),
                min_ani=min_ani,
                threads=args.threads,
                out=raw,
            ),
            env=env,
            tool="sylph profile",
            hint="Something went wrong profiling the samples with sylph.",
        )
        ckpt.done("profile")
    else:
        comm("skipping sylph profile (already done; --resume)")

    if dry_run():
        comm("(dry run) would rewrite %s as %s" % (raw, table))
        _heatmap(table, args.output, len(pairs))
        return

    announcement("EXTRACTING AVERAGE ABUNDANCE OF EACH BIN")
    with open(table, "w") as out:
        rows = _sylph_table(raw, pairs, bins, out)
    comm("Abundance table for %d bins stored in %s" % (rows, table))
    _heatmap(table, args.output, len(pairs))


#: The quantifiers ``--quantifier`` can choose between, in the order it offers them. Each takes the
#: same arguments and leaves the same ``bin_abundance_table.tsv`` behind, so adding a fourth is a
#: matter of writing its function and naming it here - and the flag's choices come from this dict, so
#: there is no second list to keep in step.
def _quantify_depth(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Aggregate the per-contig depth table the pipeline already built into a per-bin abundance.

    No alignment, no index, no third-party tool, and no second opinion about what a contig's coverage
    is. `binning` and `coverage` already measure every contig's depth in every sample - that is the
    table the binner itself read - and a bin's abundance is the length-weighted mean of its contigs'
    depths. Running salmon or CoverM to find that out again aligns the same reads a second time and
    answers on a different scale: CoverM's numbers are twice a jgi-style depth on strain-rich data,
    because jgi discards reads below 97% identity and CoverM does not.

    So this is the baseline the other quantifiers should be judged against, and it is close to free:
    one pass over a TSV against the hours an alignment costs.

    Length-weighted mean rather than median, deliberately, and the difference matters for a bin made of
    one long contig and twenty short ones. The mean is what "how much of this genome was in the sample"
    means - total bases over total length - whereas an unweighted median asks a question about contigs.
    The median over contigs is reported as a second column for anyone who wants the robust version.
    """
    announcement("AGGREGATING THE CONTIG DEPTH TABLE INTO PER-BIN ABUNDANCE")
    if not args.depth:
        error(
            "--quantifier depth needs the contig depth table that binning or coverage wrote.\n"
            "Pass it with --depth, e.g. --depth BINNING_OUT/work_files/contig_depth.txt\n"
            "It is the same table the binner read, which is the point: no realignment and no second "
            "opinion about what a contig's coverage is."
        )
    if args.assembly:
        warning(
            "--quantifier depth reads the depth table and the bins; the assembly (-a) is unused"
        )

    import pandas as pd

    frame = pd.read_csv(args.depth, sep="\t")
    if frame.columns[0] != "contigName" or "contigLen" not in frame.columns:
        error(
            "%s is not a contig depth table: its first columns are %s.\n"
            "Build one with `metawrap2 coverage` or take it from a binning run's work_files/."
            % (args.depth, ", ".join(map(str, frame.columns[:3])))
        )
    samples = [
        c
        for c in frame.columns
        if c not in ("contigName", "contigLen", "totalAvgDepth") and not c.endswith("-var")
    ]
    if not samples:
        error("%s has no per-sample depth column, so there is nothing to aggregate." % args.depth)
    comm("%d contigs, %d sample column(s): %s" % (len(frame), len(samples), ", ".join(samples)))
    frame = frame.set_index("contigName")

    rows = []
    for bin_path in fasta_files(args.bins):
        name = os.path.splitext(os.path.basename(bin_path))[0]
        contigs = [contig_id(h) for h, _seq in iter_fasta(bin_path)]
        present = frame.index.intersection(contigs)
        if len(present) == 0:
            warning(
                "%s: none of its %d contigs are in the depth table, so it has no abundance. Are the "
                "bins and the table from the same assembly?" % (name, len(contigs))
            )
            continue
        if len(present) < len(contigs):
            comm("%s: %d of %d contigs found in the table" % (name, len(present), len(contigs)))
        sub = frame.loc[present]
        lengths = sub["contigLen"].to_numpy(dtype=float)
        row = {"bin": name, "length": int(lengths.sum()), "contigs": len(present)}
        for sample in samples:
            depths = sub[sample].to_numpy(dtype=float)
            row[sample] = round(float((depths * lengths).sum() / lengths.sum()), 6)
            row["%s-median" % sample] = round(float(pd.Series(depths).median()), 6)
        rows.append(row)
    if not rows:
        error("No bin had any contig in the depth table; nothing to write.")

    out = pd.DataFrame(rows).set_index("bin")
    out.to_csv(table, sep="\t")
    comm("wrote %s: %d bins x %d samples" % (table, len(out), len(samples)))


QUANTIFIERS = {
    "depth": _quantify_depth,
    "salmon": _quantify_salmon,
    "coverm": _quantify_coverm,
    "sylph": _quantify_sylph,
}


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not os.path.isdir(args.bins):
        error("%s does not exist. Exiting..." % args.bins)

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "quant_bins",
        [
            # Gzipped bins are accepted here specifically, unlike in the modules that pass a bin
            # file to prokka or CheckM1: both quantifiers read compression without noticing it -
            # salmon's path goes through iter_fasta, and CoverM reads a gzipped genome FASTA
            # directly (verified: 100% of reads mapped, to the right genome).
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS_WITH_GZ),
                args.bins,
                "bin folder (-b)",
            )
        ]
        + ([(check_fasta, args.assembly, "assembly (-a)")] if args.assembly else [])
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    # All three quantifiers live in this module's env. CoverM used to have one of its own for its
    # size; it is a single static binary and the separate solve bought nothing but a second install
    # command, so --quantifier now never needs one.
    env = env_for("quant_bins", settings)
    rec = start_run(
        "quant_bins", args, env, settings, inputs=[args.assembly, args.bins] + list(args.reads)
    )

    ckpt = make_checkpoint(args.output)
    try:
        ensure_dir(args.output)
        # Before anything expensive: a read file that cannot be paired up is worth hearing about
        # now rather than after an index or a mapping run.
        if args.quantifier == "depth":
            if args.reads:
                warning(
                    "--quantifier depth aggregates the depth table; the %d read file(s) are unused"
                    % len(args.reads)
                )
            pairs = []
        else:
            if not args.reads:
                error(
                    "--quantifier %s needs read files. Only --quantifier depth works without them, "
                    "by aggregating a depth table that binning or coverage already wrote."
                    % args.quantifier
                )
            try:
                pairs = collect_read_pairs(args.reads)
            except ValueError as e:
                error(str(e))

        table = os.path.join(args.output, "bin_abundance_table" + TSV_SUFFIX)
        QUANTIFIERS[args.quantifier](args, pairs, table, env, ckpt)

        announcement("QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
