"""MetaWrap2 quant_bins module: estimate the abundance of each bin in each sample.

Four quantifiers are available. All write the same bins x samples ``bin_abundance_table.tsv`` and,
for more than one sample, a clustered Seaborn heatmap of it.

**sylph is the default**, based on a benchmark of 76 bins from the CAMI2 gastrointestinal
gold-standard assembly, quantified across three samples against a read-level truth
(``metawrap2 benchmark run cami2_gi quantifier``):

    quantifier   Spearman  Pearson   L1 sim   CPU-hours  peak RSS
    sylph          0.868    0.811    0.722      0.126     1.3 GB
    salmon         0.827    0.959    0.835      1.078     7.4 GB
    CoverM         0.803    0.952    0.791      2.368    11.6 GB
    depth          0.780    0.953    0.799      0.013     small

Spearman measures ordering and Pearson measures magnitude; no quantifier is best at both. sylph
orders the bins best at a tenth of salmon's CPU time, but has the lowest Pearson (0.811 against
salmon's 0.959). Use sylph to rank bins and follow them across samples, and ``--quantifier salmon``
when the values themselves matter, such as a ratio between two bins.

``depth`` aligns nothing. It aggregates the per-contig depth table the binning run already wrote
into a length-weighted mean per bin, costing 0.013 CPU-hours (10x less than sylph, 80x less than
salmon, 180x less than CoverM). Its magnitudes match salmon and CoverM (Pearson 0.953, L1 0.799,
both ahead of CoverM), but its ordering is the worst. Pearson is dominated by abundant bins while
Spearman weights all bins equally, so noise on low-coverage bins lowers Spearman more. Use it when
magnitudes matter and the depth table already exists.

Methods: salmon quantifies reads against a pseudo-transcriptome of the bins' contigs and collapses
per-contig counts to a length-weighted median coverage per bin. CoverM, in genome mode, maps
against the bins and reports per-genome coverage directly. sylph builds k-mer sketches of the bins
and the reads and profiles one against the other, so its cost grows little with sequencing depth
and it needs 1.3 GB where CoverM needs 11.6. ``depth`` reads no reads: it is one pass over an
existing TSV, on the same scale the binner clustered on.

Downstream steps accept the table from any quantifier, but **the values are not interchangeable**:
salmon gives a length-weighted median contig coverage, CoverM (default method) a relative
abundance as a percentage of the sample, and sylph a k-mer-derived sequence abundance, also a
percentage. Rows of a CoverM or sylph table sum to at most 100; rows of a salmon table do not.

The two percentage tables also differ in what zero means. sylph omits genomes whose sketch matches
below the ANI threshold, so a zero in a sylph column means "not detected above ``SYLPH_MIN_ANI``".
A zero in a CoverM column means "under ``COVERM_MIN_COVERED_FRACTION`` of the genome was covered".

The commands each tool runs are in the COMMANDS block below, and the tunable numbers in the
CONSTANTS block.
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
# One invocation for all samples: CoverM takes repeated --coupled pairs and writes one
# genome x sample table.
#
# The bins are given as a list file instead of --genome-fasta-directory/-x. CoverM's -x takes a
# single extension (default `fna`); a wrong one gives an empty genome set without an error, and a
# folder mixing .fasta and .fasta.gz cannot be described. The list also keeps .stats.tsv files and
# checkm/ subdirectories in a bin folder out of CoverM's input.
COVERM_GENOME = (
    "coverm genome --genome-fasta-list {genome_list} "
    "--coupled {coupled} --methods {methods} --min-covered-fraction {min_covered} "
    "--min-read-percent-identity {min_identity} --min-read-aligned-percent {min_aligned} "
    "--mapper {mapper} --threads {threads} --output-file {out}"
)
# sylph sketches the bins (one .syldb) and the reads (one .sylsp per sample), then `profile`
# compares the sketches. The genome sketch takes the same list file CoverM is given.
SYLPH_SKETCH_GENOMES = "sylph sketch --gl {genome_list} -o {db} -t {threads}"
# --l1/--l2 are lists of forward and reverse read files, matched by position, so every sample is
# sketched in one invocation. The sketches land in -d as <forward read filename>.paired.sylsp.
SYLPH_SKETCH_READS = "sylph sketch --l1 {l1} --l2 {l2} -d {sketch_dir} -t {threads}"
SYLPH_PROFILE = (
    "sylph profile {db} {sketches} --minimum-ani {min_ani} -t {threads} --output-file {out}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: CoverM's genome-mode default. Only one method is used because the table has one cell per bin
#: per sample.
COVERM_METHOD = "relative_abundance"
#: Percentage of a genome's length that must be covered before CoverM reports a non-zero value.
#: CoverM's default of 10% is kept, so a few reads on a conserved gene do not count as presence.
COVERM_MIN_COVERED_FRACTION = 10
#: Alignment filters. CoverM's defaults are 0 and 0, so every alignment counts however poorly it
#: matches, and reads from a relative of a binned genome are counted as that genome's abundance.
#: 90% identity matches the depth table's floor (see depth.PERCENT_IDENTITY); 97 loses real coverage
#: on a MAG assembled from a population.
#:
#: ``--min-read-aligned-percent`` excludes reads that match over only part of their length, such as
#: a third of the read at 95% identity.
#:
#: Both are percentages, per `coverm genome --full-help` ("Exclude reads by overall percent identity
#: e.g. 95 for 95%"). A 0-1 scale given 90 would accept everything; a units mix-up of this kind with
#: skani's alignment fraction merged two Enterococcus species in dRep.
#:
#: The values have not yet been tuned with the quantifier benchmark.
COVERM_MIN_READ_IDENTITY = 90
COVERM_MIN_READ_ALIGNED = 75
#: CoverM's current default, pinned because it changed between releases (0.8 defaults to
#: strobealign, earlier versions did not) and a different mapper changes the numbers.
#: `coverm genome --full-help` lists ``minimap2-sr``, ``bwa-mem`` and ``bwa-mem2`` as alternatives;
#: the choice has not been benchmarked.
COVERM_MAPPER = "strobealign"
#: A bin whose sketch matches the sample below this ANI is reported as absent. This setting has the
#: largest effect in the module. On 76 CAMI2 GI bins across three samples, sylph's default of 95
#: assigned 55% of the abundance in one sample to four bins whose genomes are absent from that
#: sample; close relatives match at 95% ANI, a species-level threshold. Raising it helps up to 99:
#:
#:     ANI    Spearman  Pearson    L1   misattributed mass
#:      95      0.758    0.678   0.832       55%
#:      97      0.823    0.744   0.728       44%
#:      99      0.868    0.811   0.555       25%      <- default
#:      99.5    0.857    0.800   0.536       25%
#:      99.9    0.696    0.840   0.448       30%   (13 of 29 present bins missed entirely)
#:
#: Above 99 the threshold discards bins that are present and rank accuracy drops, though the
#: remaining estimates improve. Below 95 sylph warns "Smaller than 95 for profile will give
#: inaccurate results", since its containment statistics are calibrated for same-species
#: comparison.
#:
#: These bins come from the pooled gold-standard assembly and are near-complete and low-error. Bins
#: from a fragmented real assembly match their genome less closely, so a strict threshold may cost
#: more there; hence the ``--min-ani`` flag.
SYLPH_MIN_ANI = 99
#: Below this sylph's containment statistics are not calibrated, so --min-ani warns.
SYLPH_ANI_FLOOR = 95
#: Which of sylph's two abundance columns fills the table. ``Sequence_abundance`` is the fraction of
#: the sample's sequence a bin accounts for, which matches what salmon and CoverM measure and what
#: the CAMI2 read-to-genome truth counts. ``Taxonomic_abundance`` is the fraction of organisms,
#: corrected for genome length, and is not comparable with the other quantifiers or the truth.
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
        help="fasta file with the entire metagenomic assembly (recommended)",
    )
    p.add_argument(
        "--depth",
        help="contig depth table to aggregate, for --quantifier depth: the table binning or "
        "coverage wrote (work_files/contig_depth.txt, or coverage's per-sample slice)",
    )
    p.add_argument(
        "--quantifier",
        choices=tuple(QUANTIFIERS),
        default="sylph",
        help="sylph (default): sequence abundance from k-mer sketches, no alignment. salmon: "
        "length-weighted median contig coverage per bin. coverm: per-genome relative abundance, "
        "mapping against the bins directly. The three give different numbers; see the "
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
    # Not required: --quantifier depth reads an existing depth table and opens no reads.
    p.add_argument(
        "reads",
        nargs="*",
        help="paired read files (*_1.fastq/*_2.fastq, .gz accepted); not needed with "
        "--quantifier depth, which aggregates an existing depth table",
    )
    return p.parse_args(argv)


def _write_genome_list(bin_folder: str, path: str) -> List[str]:
    """Write the bins to quantify, one absolute path per line; return the paths written.

    Every FASTA in the folder is listed, compressed or not, so a folder mixing ``.fasta`` and
    ``.fasta.gz`` bins works in one pass. Non-FASTA files (``.stats.tsv``, tool subdirectories) are
    skipped.

    CoverM reads this as ``--genome-fasta-list`` and sylph as ``--gl``. The returned paths are used
    to fill in the bins sylph did not report. Sorted, so the order is the same on every run.
    """
    bins = fasta_files(bin_folder)
    with open(path, "w") as fh:
        fh.writelines("%s\n" % b for b in bins)
    return bins


def _coverm_table(raw: str, samples: List[str], out: TextIO) -> int:
    """Rewrite CoverM's genome table as the ``Genomic bins`` table, returning the row count.

    Two translations make it match the salmon table:

    * **Row labels go through :func:`bin_label`**, as in the salmon path, so both tables key a bin
      the same way.
    * **Columns are matched by position.** CoverM titles each column
      ``<file> Relative Abundance (%)`` after the read filename, in the order the ``--coupled``
      pairs were given, i.e. sample order. Our own sample labels are used as headers.

    CoverM's ``unmapped`` row is dropped since it is not a bin. Rows are sorted by bin, because
    CoverM's own order is neither file nor alphabetical order (it has written ``bin_2`` before
    ``bin_1``).
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

    * **Long format.** sylph writes one row per (sample, genome) detection; it is pivoted into the
      bins x samples shape the rest of MetaWrap2 reads.
    * **Undetected genomes are omitted.** sylph leaves out genomes below ``SYLPH_MIN_ANI``, so the
      bins are taken from the genome list and unreported cells are filled with zero.
    * **File-path names.** ``Sample_file`` is the forward read file the sketch came from (not the
      ``.sylsp`` name), which maps back to our sample labels. ``Genome_file`` goes through
      :func:`bin_label`, as in the other quantifiers, so all tables key a bin the same way.

    Columns are found by name, since sylph's header uses fixed names (unlike CoverM's, which are
    titled after read filenames).
    """
    labels = [sample for sample, _, _ in samples]
    #: Keyed by forward read path, with the basename as a fallback in case sylph normalises paths.
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

    Shared by all quantifiers, since the table has the same shape whichever tool filled it in.
    """
    if samples < 2:
        warning("Cannot make a clustered heatmap from one sample. Skipping the heatmap")
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
    """Index a pseudo-transcriptome, quantify every sample against it, and collapse to bins."""
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

        # Clean up regardless of sample count: salmon's per-sample output is large and is fully
        # summarised by now. It was once deleted only when there was more than one sample.
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

    CoverM maps against the bin fastas, so ``-a`` is unused. All samples go in one command, so
    there is one checkpoint for the whole invocation. Converting the table is cheap and is redone on
    every resume.
    """
    announcement("QUANTIFYING BINS ACROSS SAMPLES WITH COVERM")
    if args.assembly:
        warning("--quantifier coverm maps reads against the bins; the assembly (-a) is unused")
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
                min_identity=COVERM_MIN_READ_IDENTITY,
                min_aligned=COVERM_MIN_READ_ALIGNED,
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

    sylph pairs the two lists by position. They are written from the same ``pairs`` the other
    quantifiers use, so all quantifiers see the same pairing.
    """
    l1 = os.path.join(output, "sylph_reads_1.txt")
    l2 = os.path.join(output, "sylph_reads_2.txt")
    for path, index in ((l1, 1), (l2, 2)):
        with open(path, "w") as fh:
            fh.writelines("%s\n" % pair[index] for pair in pairs)
    return l1, l2


def _quantify_sylph(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Sketch the bins and the reads, profile one against the other, write the abundance table.

    Both sides become k-mer sketches and ``profile`` compares them, so the profiling cost does not
    depend on sequencing depth. The two sketching steps are checkpointed separately because the
    read sketches are expensive and should not be redone after an interruption.

    As in the CoverM path, ``-a`` is unused: the bins are the reference.
    """
    announcement("QUANTIFYING BINS ACROSS SAMPLES WITH SYLPH")
    if args.assembly:
        warning("--quantifier sylph profiles the bins themselves; the assembly (-a) is unused")
    min_ani = getattr(args, "min_ani", SYLPH_MIN_ANI)
    if min_ani < SYLPH_ANI_FLOOR:
        warning(
            "--min-ani %g is below %g, where sylph's containment statistics are no longer "
            "calibrated (sylph warns about this too); the abundances are then not reliable"
            % (min_ani, SYLPH_ANI_FLOOR)
        )
    comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

    genome_list = os.path.join(args.output, "sylph_genomes.txt")
    bins = _write_genome_list(args.bins, genome_list)
    if not bins:
        error("No bins to quantify in %s. Exiting..." % args.bins)
    comm("%d bins to quantify" % len(bins))
    l1, l2 = _write_read_lists(pairs, args.output)

    # `sylph sketch -o bins` writes bins.syldb; sylph adds the suffix.
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
        # sylph names the sketches (forward read name + .paired.sylsp), so collect them from the
        # directory. Check the count: a sample sylph did not sketch would otherwise become a column
        # of zeroes.
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


#: QUANTIFIERS (below) lists the ``--quantifier`` choices. Each function takes the same arguments
#: and writes the same ``bin_abundance_table.tsv``; to add one, write its function and add it to
#: the dict, which also supplies the flag's choices.
def _quantify_depth(args: argparse.Namespace, pairs: list, table: str, env, ckpt) -> None:
    """Aggregate the per-contig depth table the pipeline already built into a per-bin abundance.

    No alignment and no third-party tool. `binning` and `coverage` already measure each contig's
    depth in each sample (the table the binner read), and a bin's abundance is the length-weighted
    mean of its contigs' depths. salmon or CoverM realign the same reads and report on a different
    scale: on strain-rich data CoverM's values are twice a jgi-style depth, because jgi discards
    reads below 97% identity and CoverM does not.

    This is the baseline for the other quantifiers, and costs one pass over a TSV.

    The length-weighted mean (total bases over total length) is used because it measures how much
    of the genome was in the sample; an unweighted median describes contigs, and differs for a bin
    of one long contig and twenty short ones. The median over contigs is written as a second column.
    """
    announcement("AGGREGATING THE CONTIG DEPTH TABLE INTO PER-BIN ABUNDANCE")
    if not args.depth:
        error(
            "--quantifier depth needs the contig depth table that binning or coverage wrote.\n"
            "Pass it with --depth, e.g. --depth BINNING_OUT/work_files/contig_depth.txt\n"
            "The bins are quantified from the same coverage the binner used, with no "
            "realignment."
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
                "%s: none of its %d contigs are in the depth table, so it has no abundance. "
                "Check that the bins and the table come from the same assembly."
                % (name, len(contigs))
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
            # Gzipped bins are accepted here, unlike in modules that pass a bin to prokka or
            # CheckM1: salmon's path goes through iter_fasta, and CoverM reads a gzipped genome
            # FASTA directly (tested: 100% of reads mapped to the right genome).
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS_WITH_GZ),
                args.bins,
                "bin folder (-b)",
            )
        ]
        + ([(check_fasta, args.assembly, "assembly (-a)")] if args.assembly else [])
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    # All quantifiers are in this module's env. CoverM once had its own env; it is a single static
    # binary, so the separate env was dropped.
    env = env_for("quant_bins", settings)
    rec = start_run(
        "quant_bins", args, env, settings, inputs=[args.assembly, args.bins] + list(args.reads)
    )

    ckpt = make_checkpoint(args.output)
    try:
        ensure_dir(args.output)
        # Pair the reads before anything expensive, so a pairing error is reported early.
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
                    "--quantifier %s needs read files. Only --quantifier depth works without them; "
                    "it aggregates the depth table that binning or coverage wrote."
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
