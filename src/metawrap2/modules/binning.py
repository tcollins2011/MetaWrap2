"""MetaWrap2 binning: bin contigs with metaBAT2/MaxBin2/CONCOCT/SemiBin2/VAMB/COMEBin.

Reads back to the assembly are aligned to produce coverage, then each selected binner runs.
Run the bin_refinement module afterwards to QC, consolidate, and reassemble the bins.

``--depth <table>`` bins from an existing contig x sample depth table without aligning, so a cohort
needs one mapping pass per sample instead of one per sample per assembly. It works with metaBAT2
and MaxBin2 only; the other four binners read the alignments themselves.

Tool commands are in the COMMANDS block below and tunable numbers in the CONSTANTS block; edit
flags there.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from typing import IO, Dict, List, Optional, Sequence, Tuple

from .. import binsplit as _binsplit
from .. import checkm as _checkm
from .. import depth as _depth
from .. import limits as _limits
from .. import longreads as _long
from ..aligners import ALIGNERS, DEFAULT_ALIGNER, INDEX_SUFFIXES
from ..cliformat import bounded_int, make_parser, percent
from ..command import THREAD_LIMITS, capture_in_env
from ..config import MODULE_ENVS, load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, UNBINNED_NAME, bin_filename
from ..io.seqio import contig_id, is_binnable, iter_fasta
from ..progress import bar
from ..utils import short_tmpdir
from ._common import (
    absolutize_paths,
    add_scorer_args,
    align_read_sets,
    announcement,
    check_disk_space,
    check_fasta,
    check_fastq,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    listdir,
    make_checkpoint,
    memory_arg,
    optional_env_for,
    resolve_memory_gb,
    resolve_threads,
    run,
    scorer_notes,
    scorer_preflight,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["binning"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# The aligner-into-samtools-sort command is in metawrap2.modules._common (ALIGN_SORT,
# SORT_THREADS, align_read_sets), shared with the cohort coverage stage.
SAMTOOLS_INDEX = "samtools index -@ {threads} -b {bam}"
#: Uses the single depth table, built without ``--noIntraDepthVariance`` even for MaxBin2; see
#: :func:`_depth_table`.
METABAT2 = (
    "metabat2 -i {assembly} -a {depth} -o {out}/metabat2_bins/bin -m {metabat_len}"
    " -t {threads} --seed {seed} --unbinned"
)
MAXBIN2 = (
    "run_MaxBin.pl -contig {assembly} -markerset {markers} -thread {threads}"
    " -min_contig_length {min_len} -out {maxbin_out}/bin -abund_list {abund_list}"
    " -max_iteration {max_iteration}"
)
#: MaxBin2's EM iteration limit. 50 is MaxBin2's own default; it is passed explicitly so it can be
#: lowered with `--maxbin2-max-iteration`. Cost per iteration scales with contigs x abundance
#: columns: a user measured ~53 min/iteration on a 418k-contig water co-assembly with 36 samples,
#: so MaxBin2 alone could run ~40 hours there (the CAMI2 tier table gives 13.8 CPU-hours for the
#: whole `fast` preset).
MAXBIN2_MAX_ITERATION = 50

CONCOCT_CUTUP = "cut_up_fasta.py {assembly} -c 10000 --merge_last -b {bed} -o 0"
CONCOCT_COVERAGE = "concoct_coverage_table.py {bed} {bams}"
CONCOCT = (
    "concoct -l {min_len} -t {threads} --coverage_file {cov} --composition_file {comp}"
    " -b {concoct_out} --seed {seed}"
)
CONCOCT_MERGE = "merge_cutup_clustering.py {clustering}"
#: ``single_easy_bin`` takes one contig set and several BAMs (multi-sample coverage of one
#: assembly), which is the input unless --binsplit is given.
#: ``--compression none`` writes plain fastas, like the other binners; the default gzips each bin.
#: ``--self-supervised`` is SemiBin2's default, passed explicitly because semi-supervised training
#: needs GTDB for MMseqs2 (~100 GB, not in this env).
SEMIBIN2 = (
    "SemiBin2 single_easy_bin --input-fasta {assembly} --input-bam {bams} --output {out}"
    " --threads {threads} --min-len {min_len} --compression none --self-supervised"
    " --random-seed {seed} --tmpdir {tmp} --engine {engine}"
)
#: ``multi_easy_bin`` runs on a concatenated catalogue whose contig names carry a sample prefix
#: (--binsplit). ``-s`` is the separator between sample label and contig name (SemiBin2's default
#: is ``:``; MetaWrap2 writes ``~``, see binsplit.SEPARATOR).
#:
#: single_easy_bin fits one sample's contigs with every sample's coverage; multi_easy_bin fits all
#: samples' contigs at once, so an organism seen in several samples contributes evidence from each.
SEMIBIN2_MULTI = (
    "SemiBin2 multi_easy_bin --input-fasta {assembly} --input-bam {bams} --output {out}"
    " --threads {threads} --min-len {min_len} --compression none --self-supervised"
    " --random-seed {seed} --tmpdir {tmp} --engine {engine} -s '{separator}'"
)
#: ``--minfasta`` has no default; without it VAMB writes a cluster TSV and no FASTA files.
#:
#: ``-o`` is the binsplit separator. For a single assembly it is empty, which disables
#: binsplitting: VAMB's default ``C`` (for ``<sample>C<contig>`` names) would split bins on any
#: contig id containing a capital C, and VAMB only warns when some ids lack the separator. For a
#: concatenated catalogue (--binsplit) it is set so VAMB splits bins per sample, the mode its
#: authors report gains for.
#:
#: ``--seed`` as for SemiBin2, though VAMB does not guarantee determinism with a seed.
#:
#: ``{cuda}`` expands to ``--cuda`` or nothing. VAMB takes a boolean, SemiBin2 a word and COMEBin a
#: device string; this module's ``--device`` sets all three (see LEARNED_DEVICES).
VAMB = (
    "vamb bin default --outdir {out} --fasta {assembly} --bamdir {bamdir}"
    " -m {min_len} -p {threads} --minfasta {minfasta} --seed {seed} -o '{separator}' {cuda}"
)
#: ``run_comebin.sh`` is the only entry point (there is no ``comebin`` on PATH). Like VAMB it globs
#: ``*.bam`` from a BAM directory, so ``work_files`` is passed.
#:
#: ``-d`` comes from ``--device``. COMEBin's default is ``cuda``, which fails on a machine without a
#: GPU; see DEFAULT_LEARNED_DEVICE.
#:
#: ``-s``: ``run_comebin.sh`` forwards a seed only when given one, so by default both the network
#: and the Leiden clustering vary between runs.
COMEBIN = "run_comebin.sh -a {assembly} -o {out} -p {bamdir} -t {threads} -d {device} -s {seed}"

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are not binned
METABAT_MIN_LEN = 1500  # metaBAT2 floors the minimum contig length at 1500 bp
MAXBIN_MARKERS_BACTERIAL = 107
#: Files metaBAT2 writes under its output prefix that are not genomes; they are removed from the
#: bin set.
METABAT2_NON_BINS = ("bin.unbinned.fa", "bin.lowDepth.fa", "bin.tooShort.fa")
#: metaBAT2's ``--seed`` default is 0, which means a random seed. Three benchmark runs with the same
#: depth table and command line gave 79, 79 and 78 bins. See SEMIBIN2_SEED; metaBAT2 is also the
#: baseline the other binners are compared against.
METABAT2_SEED = 1337
#: CONCOCT's own default, passed explicitly (as COVERM_MAPPER is in quant_bins) so a change of
#: upstream default ("0 gives a random seed") cannot change results. Kept at 1, not 1337, to
#: leave existing results unchanged.
CONCOCT_SEED = 1
MAXBIN_MARKERS_UNIVERSAL = 40  # --universal (better Archaea binning)
#: SemiBin2 trains a neural network, so without a fixed seed two runs on the same input give
#: different bins, which adds noise to benchmarks and makes MAGs irreproducible. The value itself
#: is arbitrary.
SEMIBIN2_SEED = 1337
#: Where SemiBin2 writes its bins under its output directory, given the flags in SEMIBIN2
#: (reclustering on, ``--write-pre-reclustering-bins`` off). The device default for SemiBin2 and the
#: other learned binners is DEFAULT_LEARNED_DEVICE.
SEMIBIN2_BIN_DIR = "output_bins"
#: Smallest bin VAMB writes, in bases. 200 kb matches SemiBin2's minimum and roughly metaBAT2's,
#: so the binners use the same size cutoff. VAMB has no default; see VAMB.
VAMB_MIN_FASTA = 200000
#: VAMB rejects a minimum contig length below 250 and recommends 2000. Only the 250 minimum is
#: enforced, so -l is respected.
VAMB_MIN_CONTIG_LEN = 250
#: Where VAMB writes its bins, under its output directory.
VAMB_BIN_DIR = "bins"
#: VAMB does not guarantee determinism even with a seed, but a seed removes some variation. See
#: SEMIBIN2_SEED.
VAMB_SEED = 1337
#: See SEMIBIN2_SEED. COMEBin sets no seed by default.
COMEBIN_SEED = 1337
#: What ``--device`` accepts. ``cuda:N`` is also accepted and passed through, so this is checked as
#: a prefix, not as argparse ``choices``.
LEARNED_DEVICES = ("cpu", "cuda", "auto")
#: Default device for the three learned binners.
#:
#: All three can train on a GPU but default to the CPU: most users do not have a GPU, and a GPU run
#: and a CPU run of the same seeded network need not agree, so results would depend on the host.
#:
#: The GPU path is much faster: COMEBin on CPU used over 280 CPU-hours on three CAMI2 samples
#: without finishing, against 7.1 for metaBAT2 and MaxBin2 together. Use ``--device auto`` or
#: ``--device cuda``; `benchmark` measures the GPU arms separately.
DEFAULT_LEARNED_DEVICE = "cpu"
#: Where COMEBin writes its final bin set. The intermediate clusterings under
#: ``comebin_res/cluster_res`` are candidates, selected among on single-copy marker genes.
COMEBIN_BIN_DIR = os.path.join("comebin_res", "comebin_res_bins")
#: The minimum contig length COMEBin bins at. ``run_comebin.sh`` does not forward ``--contig_len``,
#: so the underlying defaults (1000/1001) always apply. When ``-l`` is higher, the assembly given to
#: COMEBin is pre-filtered; see _comebin_assembly.
COMEBIN_MIN_LEN = 1000
# ───────────────────────────────────────────────────────────────────────────────────────────

#: Every binner this module can run, in the order their bin sets are passed downstream. Each name
#: is both the flag (``--<name>``) and the output directory (``<name>_bins``). The `run` driver
#: reads this list for its -A/-B/-C wiring; a separate list there once omitted the learned binners,
#: so a sheet selecting ``--semibin2`` had its bins dropped before refinement.
BINNERS: Tuple[str, ...] = ("metabat2", "maxbin2", "concoct", "semibin2", "vamb", "comebin")

#: Where each binner's adopted bins land under the module's output directory.
BIN_DIRS: Tuple[str, ...] = tuple("%s_bins" % b for b in BINNERS)

#: Binners that read the alignments instead of a depth table, so they cannot use ``--depth``.
#: CONCOCT builds its own coverage table from the BAMs; SemiBin2, VAMB and COMEBin read the BAMs.
#: The error message names the binner that needs the alignment.
BAM_BINNERS: Tuple[str, ...] = ("concoct", "semibin2", "vamb", "comebin")

#: Binner tiers, chosen by available compute. Measured on CAMI2 gastrointestinal, three samples,
#: after bin_refinement. "Clean HQ" counts bins that are >=70% complete, <=10% contaminated and
#: non-chimeric by GUNC, joined per bin:
#:
#:     tier       binners                        HQ   clean HQ   clean/HQ   cpu-hours
#:     very-fast  metaBAT2                        37      35       94.6%       4.8
#:     fast       metaBAT2 + MaxBin2              71      59       83.1%       7.1
#:     default    metaBAT2 + SemiBin2            131     114       87.0%      86.8  (17 on a GPU)
#:     gpu        metaBAT2 + SemiBin2 + COMEBin  183     177       96.7%      52.8  (on a GPU)
#:
#: End to end on per-sample assemblies (reads in, MAGs out), the four tiers return 36, 41, 53 and 60
#: clean MAGs for 8.8, 13.8, 43.8 and 120.7 cpu-hours of binning and refinement, against metaWRAP
#: 1.3.0's 48 for 52.5. `default` is the first tier measured that beats metaWRAP 1 on both counts.
#:
#: **Why `gpu` is three binners rather than COMEBin alone.** On the gold assembly the trio and COMEBin
#: alone tie - 177 clean MAGs against 176 - which is what the ensemble used to be dropped for. End to
#: end on per-sample assemblies they do not tie: the trio returns 60 clean MAGs against COMEBin's 57.
#: The gold assembly flatters a single strong binner, because there the contigs are long and clean
#: enough that one model can place nearly everything; on a real assembly there is more that one model
#: gets wrong and another does not.
#:
#: Members were chosen by complementarity: each clean high-quality bin was assigned to the reference
#: genome owning most of its sequence, and unions were taken over pairs and trios. The binners are
#: highly redundant here (COMEBin recovers 176 of the best trio's 182 genomes, 39 found by no other
#: binner), and one dataset cannot show whether COMEBin's lead is general. Contrastive (COMEBin),
#: self-supervised (SemiBin2) and composition-plus-coverage (metaBAT2) binners fail differently, and
#: consolidation keeps the best version of each bin.
#:
#: metaBAT2 is the third member on cost. The best trio by union, COMEBin + SemiBin2 + MaxBin2,
#: reaches 182 genomes against 181, but MaxBin2 costs 24.5 cpu-hours against metaBAT2's ~1.5, and
#: only 37.7% of MaxBin2's bins pass GUNC. metaBAT2 reuses the depth table SemiBin2 already built;
#: VAMB costs 4.9 cpu-hours on a GPU (67 on a CPU) and CONCOCT 238. metaBAT2 is also the most widely
#: validated of the six.
#:
#: The tiers are not supersets of one another and "stepping up can only add genomes" does not hold:
#: `fast` runs MaxBin2, which `default` and `gpu` do not. What holds is that each tier recovers more
#: clean MAGs than the one below it.
#:
#: `gpu` needs a GPU. COMEBin on CPU reached epoch 20 of 200 after 305 cpu-hours and was stopped;
#: on one GPU the same run finishes in 40. `gpu` therefore defaults to `--device cuda` (see
#: GPU_PRESETS).
#:
#: `default` needs the learned-binner env (~1.5 GB of pytorch), for 114 clean MAGs against `fast`'s
#: 59. `fast` (the previous default) is for a quick look.
PRESETS: Dict[str, Tuple[str, ...]] = {
    "very-fast": ("metabat2",),
    "fast": ("metabat2", "maxbin2"),
    "default": ("metabat2", "semibin2"),
    # Three different algorithm families. See the note above for the margin over COMEBin alone, and
    # for why the third member is metaBAT2 and not the binner that adds one genome to the union.
    "gpu": ("metabat2", "semibin2", "comebin"),
}

#: Older preset names that are still accepted, mapped to the current name. `gpu` was called `slow`
#: before it was renamed after the hardware it needs.
PRESET_ALIASES: Dict[str, str] = {"slow": "gpu"}

#: Presets whose learned binners need a GPU to finish in reasonable time. They default to
#: ``--device cuda``; an explicit ``--device`` overrides this.
GPU_PRESETS = frozenset({"gpu"})


def resolve_preset(name: Optional[str]) -> Optional[str]:
    """The current name of preset *name*, accepting the older aliases."""
    if name is None:
        return None
    return PRESET_ALIASES.get(name, name)


#: The tier a pipeline run uses when its sheet names no binner. Named here rather than as a second
#: list of flags in the driver, so the shipped default and the tier table cannot disagree.
#:
#: `default` beats metaWRAP 1 on both counts end to end on per-sample assemblies: 53 clean
#: high-quality MAGs against 48, for 43.8 cpu-hours against 52.5. `fast`, the previous default,
#: returns 41 for 13.8.
#:
#: A plain run therefore needs metawrap2-learned (~1.5 GB of pytorch), which `quickstart` and
#: `install-env --binning-only` build. `--preset fast` restores the old behaviour.
DEFAULT_PRESET = "default"

#: The first three columns of a jgi-format depth table. Checked because metaBAT2 bins any TSV
#: without complaint.
DEPTH_HEADER: Tuple[str, ...] = ("contigName", "contigLen", "totalAvgDepth")


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 binning",
        usage="metawrap2 binning [options] -a assembly.fa -o output_dir "
        "reads_1.fastq reads_2.fastq ...",
        description="Bin assembly contigs with metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB and/or "
        "COMEBin. "
        "Provide each replicate's reads separately (.gz accepted).",
    )
    p.add_argument("-a", "--assembly", required=True, help="metagenomic assembly fasta")
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    memory_arg(p)
    p.add_argument(
        "-l",
        "--min-len",
        type=bounded_int(1),
        default=MIN_CONTIG_LEN,
        help="minimum contig length to bin; metaBAT2 never goes below 1500",
    )
    p.add_argument(
        "--preset",
        choices=tuple(PRESETS) + tuple(PRESET_ALIASES),
        help="which binners to run, chosen by how much compute you have (see PRESETS for the measured "
        "basis). Ignored if you name binners yourself. 'gpu' adds COMEBin and trains SemiBin2 and "
        "COMEBin on a GPU (--device cuda unless --device is given); on CPU COMEBin did not finish "
        "three CAMI2 samples in 280 CPU-hours. 'slow' is the old name for 'gpu'",
    )
    p.add_argument("--metabat2", action="store_true", help="bin with metaBAT2")
    p.add_argument("--maxbin2", action="store_true", help="bin with MaxBin2")
    p.add_argument("--concoct", action="store_true", help="bin with CONCOCT")
    p.add_argument(
        "--semibin2",
        action="store_true",
        help="bin with SemiBin2 (needs `metawrap2 install-env learned`, the env with pytorch "
        "shared by the three learned binners)",
    )
    p.add_argument(
        "--vamb",
        action="store_true",
        help="bin with VAMB (needs `metawrap2 install-env learned`, the env with pytorch "
        "shared by the three learned binners)",
    )
    p.add_argument(
        "--comebin",
        action="store_true",
        help="bin with COMEBin (needs `metawrap2 install-env learned`, the env with pytorch "
        "shared by the three learned binners)",
    )
    p.add_argument(
        "--device",
        default=None,
        metavar="{cpu,cuda,cuda:N,auto}",
        help="where the learned binners (SemiBin2, VAMB, COMEBin) train. 'auto' uses a GPU if "
        "torch can see one in that binner's env and the CPU otherwise; the default is cpu so that "
        "a run is reproducible on any host (see DEFAULT_LEARNED_DEVICE), except under --preset gpu, "
        "which defaults to cuda",
    )
    p.add_argument(
        "--gpu-strict",
        "--strict-gpu",
        dest="strict_gpu",
        action="store_true",
        help="refuse to start a learned binner on a GPU another process is using, instead of "
        "warning and sharing it. Use it for benchmarks: a shared card inflates wall time and the "
        "GPU trace covers both jobs. Same as `[limits] gpu_strict = true`",
    )
    p.add_argument(
        "--universal", action="store_true", help="use universal (not bacterial) markers in MaxBin2"
    )
    p.add_argument(
        "--maxbin2-max-iteration",
        type=bounded_int(1),
        default=MAXBIN2_MAX_ITERATION,
        metavar="N",
        help="maximum MaxBin2 EM iterations (default: %d, MaxBin2's own). Each iteration "
        "costs more with more contigs x samples: ~53 min/iteration on a "
        "418k-contig, 36-sample co-assembly, so 50 iterations take ~40 hours"
        % MAXBIN2_MAX_ITERATION,
    )
    p.add_argument(
        "--run-checkm",
        action="store_true",
        help="score each binner's bins with CheckM2 (CheckM1 with --checkm1)",
    )
    add_scorer_args(p, "the bins (with --run-checkm)")
    p.add_argument(
        "--long-read-type",
        choices=sorted(_long.READ_TYPES),
        default=None,
        help="the reads are long reads of this platform. Sets minimap2's preset and "
        "SemiBin2's sequencing type; long reads mapped with the short-read preset give "
        "wrong coverage",
    )
    p.add_argument(
        "--aligner",
        choices=sorted(ALIGNERS),
        default=DEFAULT_ALIGNER,
        help="aligner used to measure coverage (default %s; see ALIGNERS for the reasons)"
        % DEFAULT_ALIGNER,
    )
    p.add_argument(
        "--keep-bams",
        action="store_true",
        help="keep the per-sample BAMs in work_files/ after the depth tables are built "
        "(about 5 GB per sample; discarded by default)",
    )
    p.add_argument(
        "--bams",
        metavar="DIR",
        help="bin from existing BAMs (what `metawrap2 coverage --keep-bams` and "
        "`binning --keep-bams` leave in work_files/), with no alignment. Use it for CONCOCT and "
        "the learned binners, which read alignments and cannot use --depth. "
        "Reads are not needed. The BAMs are read in place, never deleted, and their headers are "
        "checked against the assembly first",
    )
    p.add_argument(
        "--depth",
        help="bin from a contig x sample depth table that already exists "
        "(the layout metawrap2 coverage and jgi_summarize_bam_contig_depths both write), aligning "
        "nothing. metaBAT2 and MaxBin2 only; "
        "the other binners read the alignments. Reads are not needed",
    )
    p.add_argument(
        # A required value, not nargs="?": with an optional value, `--binsplit r_1.fastq` took the
        # read file as the separator.
        "--binsplit",
        default="",
        metavar="SEP",
        # Concatenated, not %-formatted, because argparse interpolates help strings itself.
        help="treat the assembly as a concatenated multi-sample catalogue (what `metawrap2 "
        "coverage` writes), bin it as one unit, then split the bins per sample. SEP "
        "separates the sample label from the contig name; `coverage` writes '"
        + _binsplit.SEPARATOR
        + "'. This is multi-split binning: one model is fitted on all samples' contigs, "
        "reported to recover 19%% more bins and 79%% "
        "more high-quality bins. Pass the catalogue with -a and its depth table with --depth",
    )
    p.add_argument(
        "--min-identity",
        type=percent(),
        default=_depth.PERCENT_IDENTITY,
        metavar="PERCENT",
        help="minimum alignment identity for a read to count toward a contig's depth. 90 "
        "measured better than jgi's 97: 88 bins against 74, higher purity, and 47 clean "
        "high-quality MAGs against 35 (on strain-rich data 97 removes signal). "
        "Pass 97 to reproduce a jgi depth table, or 0 for no filter",
    )
    p.add_argument("--single-end", action="store_true", help="reads are single-end (*.fastq)")
    p.add_argument("--interleaved", action="store_true", help="reads are interleaved paired-end")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument(
        "reads",
        nargs="*",
        help="read files (paired *_1.fastq/*_2.fastq, or single/interleaved); "
        "omitted only with --depth",
    )
    return p.parse_args(argv)


def _split_concoct_bins(clustering_csv: str, assembly: str, out_dir: str) -> None:
    """Split contigs into per-cluster bin files from a CONCOCT merged clustering CSV.

    ``merge_cutup_clustering.py`` has already merged ``cut_up_fasta.py``'s ``<contig>.<fragment>``
    ids back into contig names, so the keys are used verbatim. The original split_concoct_bins.py
    truncated names at the first '.', which broke contig names containing a dot (metaSPAdes'
    ``NODE_1_length_5000_cov_3.5``).
    """
    contig_bin = {}
    with open(clustering_csv) as fh:
        for line in fh:
            if line.startswith("contig_id"):
                continue
            cid, _, cluster = line.strip().partition(",")
            if cid:
                contig_bin[cid] = cluster
    ensure_dir(out_dir)
    handles: Dict[str, IO[str]] = {}
    try:
        for header, seq in bar(iter_fasta(assembly), desc="splitting CONCOCT bins", unit=" contig"):
            contig = header.split()[0]
            contig_cluster = contig_bin.get(contig)
            # CONCOCT numbers clusters from 0; bins are numbered from 1, as for other binners.
            target = (
                bin_filename(int(contig_cluster) + 1)
                if contig_cluster is not None
                else UNBINNED_NAME
            )
            if target not in handles:
                # One handle per bin, closed in the finally below, to avoid reopening per contig.
                handles[target] = open(os.path.join(out_dir, target), "w")  # noqa: SIM115
            handles[target].write(">%s\n%s\n" % (contig, seq))
    finally:
        for handle in handles.values():
            handle.close()


def _align_reads(args, env, work, assembly_copy) -> None:
    """Align each read set back to the assembly, producing sorted BAMs in work/.

    Uses :func:`metawrap2.modules._common.align_read_sets`, shared with the cohort coverage stage.
    """
    align_read_sets(
        assembly_copy,
        args.reads,
        work,
        aligner=args.aligner,
        threads=args.threads,
        env=env,
        single_end=args.single_end,
        interleaved=args.interleaved,
        long_read_type=args.long_read_type,
    )


#: What `_discard_alignments` removes: the BAMs and the aligner index beside the assembly copy.
#: `_align_reads` always rebuilds the index, which is large (1.7 GB for bwa, 5.3 GB for bwa-mem2 on
#: a 228 Mbp assembly). `assembly.fa` is kept because the binners read it.
SPENT_SUFFIXES = (".bam", ".bam.bai") + INDEX_SUFFIXES


def _discard_alignments(work: str, ckpt) -> None:
    """Delete the BAMs and the aligner index once every binner has its depth table.

    Nothing after binning reads them, and they take about 5 GB per sample plus the index (~50 GB on
    a ten-sample cohort). ``--keep-bams`` keeps them.

    The align_reads and depth_table checkpoints are cleared too, so a later ``--resume`` that adds a
    binner realigns instead of running on missing BAMs. The depth table file is kept, but is rebuilt
    after realignment because the read set may differ.
    """
    spent = [f for f in listdir(work) if f.endswith(SPENT_SUFFIXES)]
    if not spent:
        return
    freed = sum(os.path.getsize(os.path.join(work, f)) for f in spent)
    for name in spent:
        os.remove(os.path.join(work, name))
    ckpt.clear("align_reads")
    ckpt.clear("depth_table")
    comm(
        "removed %d alignment file(s) from %s, freeing %.1f GB (--keep-bams to keep them)"
        % (len(spent), work, freed / 1024.0**3)
    )


def _bam_paths(work: str) -> List[str]:
    """The BAMs in *work*, sorted, as a list."""
    return sorted(os.path.join(work, f) for f in listdir(work) if f.endswith(".bam"))


def _bams(work: str) -> str:
    bams = _bam_paths(work)
    # Under --dry-run no alignment ran, so show a placeholder in the printed command.
    if not bams and dry_run():
        return os.path.join(work, "<sample>.bam")
    return " ".join(bams)


#: Name of the shared contig x sample depth table inside ``work_files``.
DEPTH_TABLE = "contig_depth.txt"


def _depth_table(args, env, work) -> str:
    """Build the contig x sample depth table that metaBAT2 and MaxBin2 both bin from.

    One table for both binners, replacing a separate jgi pass for each.

    MaxBin2's pass used ``--noIntraDepthVariance``. In jgi_summarize_bam_contig_depths 2.18 that
    code path reports wrong depths: on a 60-contig two-sample test it zeroed 59 contigs and gave
    the last a depth of 7867 (real per-sample depths 0.267 and 0.453), and 7874 on an identical
    rerun, which points to uninitialized memory. With a single BAM it returns 0 for everything.

    On the 210,632-contig CAMI2 GI assembly with three samples, compared with the same BAMs read
    without the flag:

      sample      covered contigs   reported 0   invented coverage   worst overstatement
      sample_0            50,272       10,841                   0                  1.0x
      sample_1            65,375       37,453              24,134          410,067,316x
      sample_10           56,369       32,557              25,187        1,523,050,975x

    For two of the three samples, more than half the covered contigs were reported as uncovered,
    about 25,000 contigs with no reads were given coverage, and some depths were overstated by
    hundreds of millions. MaxBin2 still produced plausible bins because it also uses composition.
    MaxBin2 now reads the mean columns of the full table.

    The table is computed in Python by `metawrap2.depth`, which implements the same algorithm (97%
    identity filter, 75 bp trimmed at each contig end, intra-contig variance) over per-base depths
    from ``samtools depth``. Over all 210,632 contigs of the CAMI2 GI assembly it agrees with jgi to
    Pearson 0.999 on the means and 0.995 on the variances (median relative difference about 5%).
    """
    depth = os.path.join(work, DEPTH_TABLE)
    comm("making contig depth file...")
    if dry_run():
        comm("(dry run) would measure depth over the BAMs in %s" % work)
        return depth
    bams = _bam_paths(work)
    _depth.write_depth_table(
        bams,
        depth,
        samples=[os.path.splitext(os.path.basename(b))[0] for b in bams],
        env=env,
        threads=args.threads,
        percent_identity=args.min_identity,
    )
    return depth


def depth_means(header: List[str]) -> List[int]:
    """Indices of the per-sample mean-depth columns in a jgi depth table header.

    The table is ``contigName contigLen totalAvgDepth <s1> <s1>-var <s2> <s2>-var ...``, so the
    means are the columns from the fourth on whose name does not end ``-var``. Selected by name so
    tables without variance columns also work.
    """
    return [i for i in range(3, len(header)) if not header[i].endswith("-var")]


def depth_table_contigs(path: str) -> List[str]:
    """The contigs a jgi depth table covers, after confirming the file is one.

    metaBAT2 bins any TSV without error, reading column 2 as a length and the rest as depths, so
    the header is checked. A table with no per-sample column is also rejected, since metaBAT2 would
    bin on composition alone.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if tuple(header[:3]) != DEPTH_HEADER:
            raise ValueError(
                "%s is not a contig depth table: its first three columns are "
                "%s, not %s.\nBuild one with `metawrap2 coverage -a <label>=<assembly> ... -o "
                "<dir> <reads>`, or drop --depth and let binning align the reads itself."
                % (path, ", ".join(header[:3]) or "(empty)", ", ".join(DEPTH_HEADER))
            )
        if not depth_means(header):
            raise ValueError(
                "%s has no per-sample depth column, so there is no coverage to bin on. In a jgi "
                "table the samples are columns 4 onwards." % path
            )
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


def _assembly_contigs(path: str) -> List[str]:
    """Contig names in *path*, read from the header lines alone.

    Avoids :func:`iter_fasta`, which loads each sequence into memory.
    """
    from ..io.seqio import smart_open

    with smart_open(path, "rt") as fh:
        return [contig_id(line[1:].rstrip("\n")) for line in fh if line.startswith(">")]


def bam_reference_names(bam: str, env: Optional[str]) -> List[str]:
    """The reference names in *bam*'s header, in order."""
    header = capture_in_env(env, ["samtools", "view", "-H", bam])
    names = []
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    names.append(field[3:])
    return names


def check_bam_references(bams: Sequence[str], assembly: str, env: Optional[str]) -> int:
    """Confirm supplied BAMs were aligned against *assembly*. Returns the overlap count.

    As with :func:`check_depth_table`: a BAM aligned against a catalogue of several assemblies is
    well formed, but the binners that read BAMs produce no bins without an error when no contigs
    match (as happened with `coverage --keep-bams` and a label). Reads one BAM header before any
    alignment or env probe.
    """
    contigs = {contig_id(header) for header, _sequence in iter_fasta(assembly)}
    if not contigs:
        raise ValueError("%s has no contigs to compare the BAMs against" % assembly)
    first = min(bams)
    names = bam_reference_names(first, env)
    if not names:
        raise ValueError("%s has no @SQ reference lines; is it a BAM?" % first)
    overlap = contigs.intersection(names)
    if not overlap:
        raise ValueError(
            "%s names none of %s's contigs: its first reference is %r and the assembly's is %r.\n"
            "A BAM aligned against a concatenated catalogue names contigs 'label<sep>contig' and "
            "the assembly names them 'contig'; reheader the BAM, or pass the catalogue as -a."
            % (first, assembly, names[0], min(contigs))
        )
    if len(overlap) < len(contigs):
        comm(
            "%d of %d assembly contig(s) are in %s; the rest have no coverage in it"
            % (len(overlap), len(contigs), os.path.basename(first))
        )
    return len(overlap)


def stage_supplied_bams(bam_dir: str, work: str) -> List[str]:
    """Symlink the BAMs in *bam_dir* into *work*, and return the links.

    Symlinked because they can be tens of gigabytes per sample and are only read. Supplied BAMs are
    never deleted: `--keep-bams` does not apply and `_discard_alignments` is skipped.
    """
    ensure_dir(work)
    staged = []
    for name in sorted(listdir(bam_dir)):
        if not name.endswith(".bam"):
            continue
        source = os.path.abspath(os.path.join(bam_dir, name))
        link = os.path.join(work, name)
        if not os.path.lexists(link):
            os.symlink(source, link)
        staged.append(link)
        index = source + ".bai"
        if os.path.isfile(index) and not os.path.lexists(link + ".bai"):
            os.symlink(index, link + ".bai")
    return staged


def check_depth_table(path: str, assembly: str) -> int:
    """Confirm a precomputed depth table describes *assembly*. Returns the number of samples.

    A table built for a different assembly is well formed (the cohort coverage stage maps every
    sample to one concatenated catalogue), so contig names are compared. No overlap is an error;
    partial overlap is reported:

    * contigs in the assembly but not the table have no coverage, and metaBAT2 and MaxBin2 leave
      them unbinned
    * contigs in the table but not the assembly are ignored, as expected for a cohort-wide table
    """
    contigs = depth_table_contigs(path)
    if not contigs:
        error("%s has a header but no contigs; there is nothing to bin." % path)
    with open(path) as fh:
        samples = len(depth_means(fh.readline().rstrip("\n").split("\t")))
    in_assembly = set(_assembly_contigs(assembly))
    shared = len(in_assembly.intersection(contigs))
    if not shared:
        error(
            "none of the %d contigs in %s are in %s; the table is for a different "
            "assembly.\nContig names must match the first word of the assembly's headers."
            % (len(contigs), path, assembly)
        )
    comm(
        "depth table: %d contig(s) over %d sample(s), %d of the assembly's %d contig(s) covered"
        % (len(contigs), samples, shared, len(in_assembly))
    )
    if shared < len(in_assembly):
        comm(
            "%d assembly contig(s) have no row in the depth table and will be left unbinned"
            % (len(in_assembly) - shared)
        )
    if shared < len(contigs):
        comm(
            "%d row(s) in the depth table name contigs that are not in this assembly; they are "
            "ignored (expected when the table covers a whole cohort's catalogue)"
            % (len(contigs) - shared)
        )
    return samples


def _run_metabat2(args, env, work, assembly_copy, depth):
    announcement("RUNNING METABAT2")
    metabat_len = max(args.min_len, METABAT_MIN_LEN)
    comm("Starting binning with metaBAT2...")
    run(
        METABAT2.format(
            assembly=assembly_copy,
            depth=depth,
            out=args.output,
            metabat_len=metabat_len,
            threads=args.threads,
            seed=METABAT2_SEED,
        ),
        env=env,
        tool="metabat2",
    )
    # metaBAT2 also writes non-genome files under the output prefix: the --unbinned leftovers,
    # contigs rejected as too short or too low-depth, and two text summaries. They are moved to
    # metabat2_extra so only bins remain.
    bins_dir = os.path.join(args.output, "metabat2_bins")
    if dry_run():
        return
    # metaBAT2 appends "total_depth=.. sample_depths=.." to contig headers; strip it so contig ids
    # match the assembly and the other binners.
    #
    # FASTA_EXTENSIONS, not BIN_EXTENSION: the files still have metaBAT2's ".fa" names until
    # renumbering below.
    for name in sorted(listdir(bins_dir)):
        if name.endswith(FASTA_EXTENSIONS) and name not in METABAT2_NON_BINS:
            path = os.path.join(bins_dir, name)
            records = list(iter_fasta(path))
            if any(" " in h or "\t" in h for h, _ in records):
                with open(path, "w") as out:
                    out.writelines(
                        ">%s\n%s\n" % (contig_id(header), seq) for header, seq in records
                    )

    extras = os.path.join(args.output, "metabat2_extra")
    for name in sorted(listdir(bins_dir)):
        is_non_bin = name.endswith((".txt", ".tsv")) or name in METABAT2_NON_BINS
        if not is_non_bin:
            continue
        ensure_dir(extras)
        os.replace(os.path.join(bins_dir, name), os.path.join(extras, name))
        comm("moved metaBAT2's %s out of the bin set (it is not a genome) into %s" % (name, extras))

    # Renumber last, after the non-bin files are moved out, so they take no bin numbers.
    renumber_bins(bins_dir)


def renumber_bins(bins_dir: str) -> int:
    """Rename every FASTA in *bins_dir* to bin_001.fasta, bin_002.fasta, ... Returns the count.

    Bins are ordered numerically by the last number in their filename, so ``bin.2.fa`` comes before
    ``bin.10.fa``. Files without a number sort last, by name.

    Files are first moved to temporary names, since an existing ``bin_002.fasta`` could otherwise
    be overwritten.
    """

    def order(name: str) -> Tuple[int, float, str]:
        digits = re.findall(r"\d+", os.path.splitext(name)[0])
        return (0, int(digits[-1]), name) if digits else (1, float("inf"), name)

    names = sorted(
        (n for n in listdir(bins_dir) if n.endswith(FASTA_EXTENSIONS)),
        key=order,
    )
    staged = []
    for i, name in enumerate(names, start=1):
        temporary = os.path.join(bins_dir, ".renumber.%d" % i)
        os.replace(os.path.join(bins_dir, name), temporary)
        staged.append((temporary, bin_filename(i)))
    for temporary, final in staged:
        os.replace(temporary, os.path.join(bins_dir, final))
    return len(staged)


def _run_maxbin2(args, env, work, assembly_copy, depth):
    announcement("RUNNING MAXBIN2")
    # One abundance file per sample (contig name and mean depth), taken from the shared table; see
    # _depth_table.
    comm("splitting the depth table into per-sample abundance files")
    abund_list = os.path.join(work, "mb2_abund_list.txt")
    if dry_run():
        markers = MAXBIN_MARKERS_UNIVERSAL if args.universal else MAXBIN_MARKERS_BACTERIAL
        maxbin_out = os.path.join(work, "maxbin2_out")
        run(
            MAXBIN2.format(
                assembly=assembly_copy,
                markers=markers,
                threads=args.threads,
                min_len=args.min_len,
                maxbin_out=maxbin_out,
                abund_list=abund_list,
                max_iteration=getattr(args, "maxbin2_max_iteration", MAXBIN2_MAX_ITERATION),
            ),
            env=env,
            tool="run_MaxBin.pl",
        )
        return
    with open(depth) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    with open(abund_list, "w") as listing:
        for i in depth_means(header):
            sample = header[i]
            per = os.path.join(work, "mb2_%s.txt" % os.path.splitext(sample)[0])
            with open(depth) as fh, open(per, "w") as out:
                for line in fh:
                    if "totalAvgDepth" in line:
                        continue
                    cols = line.rstrip("\n").split("\t")
                    out.write("%s\t%s\n" % (cols[0], cols[i]))
            listing.write(os.path.abspath(per) + "\n")

    markers = MAXBIN_MARKERS_UNIVERSAL if args.universal else MAXBIN_MARKERS_BACTERIAL
    maxbin_out = os.path.join(work, "maxbin2_out")
    ensure_dir(maxbin_out)
    comm("Starting binning with MaxBin2...")
    run(
        MAXBIN2.format(
            assembly=assembly_copy,
            markers=markers,
            threads=args.threads,
            min_len=args.min_len,
            maxbin_out=maxbin_out,
            abund_list=abund_list,
            max_iteration=getattr(args, "maxbin2_max_iteration", MAXBIN2_MAX_ITERATION),
        ),
        env=env,
        tool="run_MaxBin.pl",
    )

    out_bins = os.path.join(args.output, "maxbin2_bins")
    ensure_dir(out_bins)
    if dry_run():
        return
    n = 0
    for f in sorted(os.listdir(maxbin_out)):
        if f.startswith("bin") and f.endswith(".fasta"):
            n += 1
            os.replace(os.path.join(maxbin_out, f), os.path.join(out_bins, bin_filename(n)))
    if n == 0:
        error("MaxBin2 produced no bins. Exiting.")
    comm("MaxBin2 finished successfully, and found %d bins!" % n)


def _run_concoct(args, env, work, assembly_copy):
    announcement("RUNNING CONCOCT")
    depth = os.path.join(work, "concoct_depth.txt")
    bed = os.path.join(work, "assembly_10K.bed")
    comp = os.path.join(work, "assembly_10K.fa")
    if dry_run() or not (os.path.isfile(depth) and os.path.getsize(depth)):
        comm("indexing .bam alignment files...")
        for f in sorted(listdir(work)):
            if f.endswith(".bam"):
                run(
                    SAMTOOLS_INDEX.format(threads=args.threads, bam=os.path.join(work, f)),
                    env=env,
                    tool="samtools index",
                )
        comm("cutting up contigs into 10kb fragments for CONCOCT...")
        run(
            CONCOCT_CUTUP.format(assembly=assembly_copy, bed=bed),
            env=env,
            tool="cut_up_fasta.py",
            log_path=comp,
        )
        comm("estimating contig fragment coverage...")
        run(
            CONCOCT_COVERAGE.format(bed=bed, bams=_bams(work)),
            env=env,
            tool="concoct_coverage_table.py",
            log_path=depth,
        )

    concoct_out = os.path.join(work, "concoct_out")
    ensure_dir(concoct_out)
    comm("Starting binning with CONCOCT...")
    run(
        CONCOCT.format(
            min_len=args.min_len,
            threads=args.threads,
            cov=depth,
            comp=comp,
            concoct_out=os.path.join(concoct_out, ""),
            seed=CONCOCT_SEED,
        ).rstrip(),
        env=env,
        tool="concoct",
        # CONCOCT forks ``-t`` workers that each use BLAS, so BLAS threads are limited. See
        # THREAD_LIMITS: unlimited, this took 30 CPU-hours per sample instead of 5.6.
        extra_env=THREAD_LIMITS,
    )
    clustering = os.path.join(concoct_out, "clustering_gt%d.csv" % args.min_len)
    merged = os.path.join(concoct_out, "clustering_gt%d_merged.csv" % args.min_len)
    comm("merging 10kb fragments back into contigs")
    run(
        CONCOCT_MERGE.format(clustering=clustering),
        env=env,
        tool="merge_cutup_clustering.py",
        log_path=merged,
    )
    comm("splitting contigs into bins")
    if not dry_run():
        _split_concoct_bins(merged, assembly_copy, os.path.join(args.output, "concoct_bins"))
    comm("CONCOCT finished successfully!")


def resolve_device(device: str, env: str) -> str:
    """Turn ``--device`` into a concrete ``cpu`` or ``cuda[:N]`` for the binner about to run.

    Only ``auto`` is resolved, by asking torch in the binner's own env, since each learned binner's
    env has its own pytorch. Any failure of the probe resolves to the CPU.
    """
    if device != "auto":
        return device
    probe = "import torch; print('cuda' if torch.cuda.is_available() else 'cpu')"
    try:
        lines = capture_in_env(env, ["python", "-c", probe]).strip().splitlines()
        chosen = lines[-1].strip() if lines else ""
    except (OSError, IndexError, ValueError):
        chosen = ""
    if chosen != "cuda":
        comm("--device auto: no usable GPU in %s, training on the CPU" % env)
        return "cpu"
    comm("--device auto: training on the GPU (torch sees one in %s)" % env)
    return "cuda"


def gpu_pin(device: str) -> Optional[Dict[str, str]]:
    """``extra_env`` pinning a child to the card ``--device cuda:N`` names, else None.

    Only COMEBin takes a device index. SemiBin2 (``--engine gpu``) and VAMB (``--cuda``) take only
    the kind, so torch uses device 0 and a run asked for another card trains on the wrong one.
    ``CUDA_VISIBLE_DEVICES`` makes the requested card the child's device 0 for all three. Bare
    ``cuda`` and ``auto`` are left unpinned.

    A malformed index returns None and is passed to the tool, which reports it.
    """
    if not device.startswith("cuda:"):
        return None
    index = device.split(":", 1)[1].strip()
    if not index.isdigit():
        return None
    # base={} returns only this variable, to merge into the child's environment.
    return _limits.gpu_env(int(index), base={})


def gpu_preflight(args, device: str) -> None:
    """Report the card this step will train on, and stop if it cannot be used.

    Called once per learned binner, after :func:`resolve_device`, because for ``--device auto``
    the device is known only after the binner's env has been probed. The check itself is
    :func:`metawrap2.limits.gpu_preflight`.

    A missing card index stops the run before a step that can take hours. A card in use only
    warns, unless ``--gpu-strict`` or ``[limits] gpu_strict`` is set.

    Skipped under ``--dry-run``, which often runs on a different machine from the real run.
    """
    if dry_run():
        return
    report = _limits.gpu_preflight(device, strict=getattr(args, "strict_gpu", False))
    if not report.message:
        return
    if not report.ok:
        error(report.message)
    (warning if report.contended else comm)(report.message)


def _run_semibin2(args, env, work, assembly_copy):
    """Bin with SemiBin2, which reads the BAMs itself rather than a depth table.

    Runs in ``metawrap2-learned`` (for pytorch), not the binning env.

    SemiBin2 reads the BAMs directly. ``--depth-metabat2`` would reuse metaBAT2's table, but
    SemiBin2 accepts it only with a pretrained ``--environment`` model, not when training on the
    sample. benchmarks/configurations.toml compares the costs.
    """
    announcement("RUNNING SEMIBIN2")
    device = resolve_device(args.device, env)
    gpu_preflight(args, device)
    out = os.path.join(work, "semibin2_out")
    comm("Starting binning with SemiBin2...")
    separator = getattr(args, "binsplit", "") or ""
    fields = {
        "assembly": assembly_copy,
        "bams": _bams(work),
        "out": out,
        "threads": args.threads,
        "min_len": args.min_len,
        "seed": SEMIBIN2_SEED,
        "tmp": short_tmpdir(os.path.join(work, "tmp-semibin2")),
        # SemiBin2 calls the accelerator "gpu"; --device uses "cuda" as the other binners do.
        "engine": "gpu" if device.startswith("cuda") else "cpu",
    }
    if separator:
        comm("multi-sample mode: SemiBin2 fits one model over every sample's contigs")
        command = SEMIBIN2_MULTI.format(separator=separator, **fields)
    else:
        command = SEMIBIN2.format(**fields)
    if getattr(args, "long_read_type", None):
        # SemiBin2's default features assume short-read coverage. The type comes from
        # metawrap2.longreads, which also sets Flye's flag and minimap2's preset.
        command += " --sequencing-type=%s" % _long.resolve(args.long_read_type).semibin2_type
        comm("SemiBin2 in long-read mode")
    run(command, env=env, tool="SemiBin2", extra_env=gpu_pin(device))

    _adopt_bins(
        os.path.join(out, SEMIBIN2_BIN_DIR), os.path.join(args.output, "semibin2_bins"), "SemiBin2"
    )


#: Binners that split bins per sample themselves when given a separator (SemiBin2 via
#: multi_easy_bin -s, VAMB via -o). The others get the generic split.
NATIVE_BINSPLIT: Tuple[str, ...] = ("semibin2", "vamb")


def _post_split(args, selected: List[str]) -> None:
    """Split each bin set that the binner did not split itself, into <bins>_split/<sample>/.

    The originals are kept: the catalogue-wide bins are what a chimerism check should see, and
    consolidation needs per-sample sets matching per-sample assemblies.
    """
    if dry_run():
        return
    for binner in selected:
        if binner in NATIVE_BINSPLIT:
            comm("%s split its own bins by sample (it was given the separator)" % binner)
            continue
        bins_dir = os.path.join(args.output, "%s_bins" % binner)
        if not os.path.isdir(bins_dir):
            continue
        out_dir = os.path.join(args.output, "%s_bins_split" % binner)
        report = _binsplit.split_bin_set(bins_dir, out_dir, args.binsplit)
        # Also a flat directory of `<sample>~bin_xxx.fasta` hard links, for tools that expect a
        # folder of FASTA files, such as the benchmark scorer, which reads top-level `*_bins`
        # directories.
        flat_dir = os.path.join(args.output, "%s_split_bins" % binner)
        linked = _binsplit.flatten_split_bins(out_dir, flat_dir)
        comm(
            "%s: %d catalogue bin(s) -> %d per-sample bin(s) across %d sample(s); %d linked into %s"
            % (
                binner,
                report.bins_in,
                report.pieces_out,
                len(report.labels),
                linked,
                os.path.basename(flat_dir),
            )
        )
        if report.unprefixed_contigs:
            warning(
                "%s: %d contig(s) had no '%s' prefix and were left out of the split. The "
                "assembly is not a catalogue built by `coverage`, so --binsplit has nothing "
                "to split" % (binner, report.unprefixed_contigs, args.binsplit)
            )
        elif not report.worthwhile:
            comm(
                "%s: no bin has contigs from more than one sample; multi-split binning gave "
                "the same bins per-sample binning would" % binner
            )


def _adopt_bins(produced: str, bins_dir: str, tool: str) -> None:
    """Move *tool*'s bins out of its own output tree into *bins_dir* and give them canonical names.

    The learned binners write bins into their own subdirectory, with their own names, beside other
    files. Bins are moved, not copied, so nothing that walks the output counts them twice.

    No output is an error, because these tools exit 0 when no cluster reaches their minimum bin
    size.
    """
    ensure_dir(bins_dir)
    if dry_run():
        return
    if not os.path.isdir(produced):
        error(
            "%s finished but wrote no %s. Check its log: when no cluster reaches its minimum "
            "bin size it exits with success and writes nothing." % (tool, produced)
        )
    for name in sorted(listdir(produced)):
        if name.endswith(FASTA_EXTENSIONS):
            os.replace(os.path.join(produced, name), os.path.join(bins_dir, name))
    n = renumber_bins(bins_dir)
    if n == 0:
        error("%s produced no bins. Exiting." % tool)
    comm("%s finished successfully, and found %d bins!" % (tool, n))


def _run_vamb(args, env, work, assembly_copy):
    """Bin with VAMB, which reads the BAM directory rather than individual BAMs or a depth table.

    Runs in ``metawrap2-learned`` (for pytorch).

    VAMB refuses an existing ``--outdir``, so it is cleared first, which would otherwise break
    ``--resume`` after an interrupted VAMB. The checkpoint has already marked this stage for a
    rerun, so the directory holds no results.
    """
    announcement("RUNNING VAMB")
    device = resolve_device(args.device, env)
    gpu_preflight(args, device)
    out = os.path.join(work, "vamb_out")
    if not dry_run() and os.path.isdir(out):
        comm("clearing %s; VAMB will not write into an existing directory" % out)
        shutil.rmtree(out)
    comm("Starting binning with VAMB...")
    run(
        VAMB.format(
            out=out,
            assembly=assembly_copy,
            bamdir=work,
            min_len=max(args.min_len, VAMB_MIN_CONTIG_LEN),
            threads=args.threads,
            minfasta=VAMB_MIN_FASTA,
            seed=VAMB_SEED,
            separator=getattr(args, "binsplit", "") or "",
            cuda="--cuda" if device.startswith("cuda") else "",
        ),
        env=env,
        tool="vamb",
        # --cuda selects a GPU; the pin selects which one. See gpu_pin.
        extra_env=gpu_pin(device),
    )
    _adopt_bins(os.path.join(out, VAMB_BIN_DIR), os.path.join(args.output, "vamb_bins"), "VAMB")


def _comebin_assembly(args, work: str, assembly_copy: str) -> str:
    """Return the assembly for COMEBin, length-filtered first if ``-l`` asks for more than it is.

    COMEBin's wrapper has no ``--contig_len`` (see COMEBIN_MIN_LEN), so a higher ``-l`` is applied
    by filtering the input, keeping COMEBin consistent with the other binners. When ``-l`` is at or
    below COMEBin's minimum, the shared copy is returned unchanged.
    """
    if args.min_len <= COMEBIN_MIN_LEN:
        return assembly_copy
    filtered = os.path.join(work, "assembly_comebin.fa")
    comm(
        "filtering the assembly to >=%d bp for COMEBin, which has no minimum-length option"
        % args.min_len
    )
    if dry_run():
        return filtered
    kept = 0
    with open(filtered, "w") as out:
        for header, seq in iter_fasta(assembly_copy):
            if len(seq) >= args.min_len:
                out.write(">%s\n%s\n" % (header, seq))
                kept += 1
    if kept == 0:
        error(
            "no contig in %s is %d bp or longer, so COMEBin has nothing to bin. Lower -l."
            % (args.assembly, args.min_len)
        )
    return filtered


def _run_comebin(args, env, work, assembly_copy):
    """Bin with COMEBin, which reads the BAM directory and picks between its own candidate bin sets.

    Runs in ``metawrap2-learned`` (for pytorch).

    COMEBin produces Leiden clusterings at several resolutions, scores them on single-copy marker
    genes, and writes the best; only that set is adopted. ``run_comebin.sh`` skips stages whose
    intermediate files exist, so the output directory is not cleared, which avoids redoing the
    expensive augmentation stage.
    """
    announcement("RUNNING COMEBIN")
    device = resolve_device(args.device, env)
    # Check before the output directory is made and the augmentation pass starts. COMEBin can run
    # for hours.
    gpu_preflight(args, device)
    pin = gpu_pin(device)
    out = os.path.join(work, "comebin_out")
    ensure_dir(out)
    comm("Starting binning with COMEBin...")
    run(
        COMEBIN.format(
            assembly=_comebin_assembly(args, work, assembly_copy),
            out=out,
            bamdir=work,
            threads=args.threads,
            # COMEBin accepts an index, but is pinned like the other two and given plain `cuda`,
            # since an index other than 0 is out of range in a one-device view.
            device="cuda" if pin else device,
            seed=COMEBIN_SEED,
        ),
        env=env,
        tool="run_comebin.sh",
        extra_env=pin,
    )
    _adopt_bins(
        os.path.join(out, COMEBIN_BIN_DIR), os.path.join(args.output, "comebin_bins"), "COMEBin"
    )


def copy_assembly(assembly: str, destination: str) -> Tuple[int, int]:
    """Copy *assembly* to *destination*, leaving out records nothing can bin.

    Returns ``(empty, uninformative)``.

    **Empty records** have a header and no sequence. bwa indexes them as ``@SQ`` lines with
    ``LN:0``, and samtools then rejects the alignment at the sort ("Header includes @SQ line ...
    with invalid or no LN: tag"), hours into the run, without naming the record.

    **Uninformative records** contain no A, C, G, T or U (in practice long runs of N). They have no
    k-mer signal and no read aligns to them. metaBAT2, MaxBin2 and CONCOCT accept them; VAMB refuses
    to start.

    CAMI2's gastrointestinal gold-standard assembly has two empty records, both headed
    ``_from__to_-42_total_0``, and 20 all-N records totalling 28.9 Mbp (the largest 6.6 Mbp), which
    are masked reference genomes.
    """
    empty = uninformative = 0
    with open(destination, "w") as out:
        for header, seq in iter_fasta(assembly):
            if not seq:
                empty += 1
                continue
            if not is_binnable(seq):
                uninformative += 1
                continue
            out.write(">%s\n%s\n" % (header, seq))
    return empty, uninformative


def _count_bins(path: str) -> int:
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
    )


def _apply_preset(args) -> None:
    """Turn ``--preset`` into binner flags, or exit with an error if no binner is selected.

    If binners are named explicitly, the preset is ignored.
    """
    named = [b for b in BINNERS if getattr(args, b)]
    given = getattr(args, "preset", None)
    preset = resolve_preset(given)
    if given and given != preset:
        comm("--preset %s is the old name for --preset %s" % (given, preset))
        args.preset = preset
    if getattr(args, "device", None) is None:
        use_gpu = preset in GPU_PRESETS and not named
        args.device = "cuda" if use_gpu else DEFAULT_LEARNED_DEVICE
    if named and preset:
        comm(
            "--preset %s ignored because binners were named "
            "explicitly: %s" % (preset, ", ".join("--" + b for b in named))
        )
        return
    if preset:
        for binner in PRESETS[preset]:
            setattr(args, binner, True)
        comm("--preset %s selects %s" % (preset, " ".join("--" + b for b in PRESETS[preset])))
        if preset in GPU_PRESETS:
            comm("--preset %s trains the learned binners on --device %s" % (preset, args.device))
        return
    if not named:
        error(
            "You must select at least one binning method: --metabat2, --maxbin2, --concoct, "
            "--semibin2, --vamb, --comebin.\n"
            "Or pick a tier with --preset {very-fast,fast,default,gpu} and let it choose - see "
            "`--help`."
        )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    _apply_preset(args)

    # --depth replaces the alignment stage, so reads are not needed. A binner that reads BAMs is an
    # error without --bams; reads given with a table are ignored.
    if args.depth:
        needs_bams = [b for b in BAM_BINNERS if getattr(args, b)]
        if needs_bams and not args.bams:
            error(
                "--depth cannot be used with %s: %s read alignments, not a "
                "depth table.\nPass the alignments with --bams <dir> alongside --depth (what "
                "`coverage --split-bams` writes, and what `cohort` does), run those binners "
                "with the reads, or select only --metabat2/--maxbin2 here."
                % (
                    ", ".join("--" + b for b in needs_bams),
                    "they" if len(needs_bams) > 1 else "it",
                )
            )
        if args.reads:
            comm(
                "--depth was given, so the %d read file(s) on the command line are not used "
                "and nothing is aligned" % len(args.reads)
            )
            args.reads = []
    elif args.bams:
        if args.reads:
            comm(
                "--bams was given, so the %d read file(s) on the command line are not used "
                "and nothing is aligned" % len(args.reads)
            )
            args.reads = []
    elif not args.reads:
        error(
            "No read files given. Coverage is measured from the reads; pass them as "
            "positional arguments, or bin from an existing depth table with --depth <table>, "
            "or from existing alignments with --bams <dir>."
        )
    if args.bams and args.depth:
        comm(
            "--bams and --depth together: metaBAT2 and MaxBin2 use the table, CONCOCT and the "
            "learned binners use the BAMs, and nothing is aligned or re-measured"
        )

    settings = load_settings(args.config)
    # Combine the flag and the config key into one attribute for the GPU steps below.
    args.strict_gpu = args.strict_gpu or _limits.gpu_strict(settings.limits)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args, extra=["depth", "bams"])
    check_disk_space("binning", settings, [args.assembly] + args.reads, args.output)
    validate_inputs(
        "binning",
        [(check_fasta, args.assembly, "assembly (-a)")]
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    if args.depth:
        # Checked before the env is resolved and the run is recorded.
        try:
            check_depth_table(args.depth, args.assembly)
        except (OSError, ValueError) as exc:
            error(str(exc))
    supplied_bams: List[str] = []
    if args.bams:
        if not os.path.isdir(args.bams):
            error("--bams %s is not a directory" % args.bams)
        supplied_bams = [
            os.path.join(args.bams, f) for f in sorted(listdir(args.bams)) if f.endswith(".bam")
        ]
        if not supplied_bams:
            error("--bams %s holds no .bam files" % args.bams)
        comm("%d BAM(s) supplied with --bams; nothing will be aligned" % len(supplied_bams))
    env = env_for("binning", settings)
    # Resolved now, not when the binner runs after alignment, so a missing env fails early.
    learned_envs = {
        flag: optional_env_for(flag, settings, "--%s" % flag)
        for flag in ("semibin2", "vamb", "comebin")
        if getattr(args, flag)
    }
    # The scorer for --run-checkm: CheckM1 is in the bin_refinement env, CheckM2 in its own.
    checkm1_env = checkm2_env = None
    if args.run_checkm:
        scorer_notes(args)
        checkm2_env = scorer_preflight(args, settings)
        if not args.checkm2:
            checkm1_env = env_for("bin_refinement", settings)
    inputs = [args.assembly] + list(args.reads) + ([args.depth] if args.depth else [])
    rec = start_run("binning", args, env, settings, inputs=inputs)

    ckpt = make_checkpoint(args.output)
    try:
        announcement(
            "BINNING FROM AN EXISTING DEPTH TABLE"
            if args.depth
            else "ALIGNING READS TO MAKE COVERAGE FILES"
        )
        ensure_dir(args.output)
        work = os.path.join(args.output, "work_files")
        ensure_dir(work)
        assembly_copy = os.path.join(work, "assembly.fa")
        comm("making copy of assembly file %s" % args.assembly)
        if not dry_run():
            empty, uninformative = copy_assembly(args.assembly, assembly_copy)
            if empty:
                comm(
                    "left out %d zero-length record(s) from the assembly; bwa would index them "
                    "as '@SQ ... LN:0', which samtools cannot read" % empty
                )
            if uninformative:
                comm(
                    "left out %d record(s) with no A, C, G, T or U in them; they have no k-mer "
                    "signal and no reads align to them" % uninformative
                )

        try:
            if args.bams:
                # Check the BAMs against the assembly first; see check_bam_references.
                if not dry_run():
                    try:
                        check_bam_references(supplied_bams, assembly_copy, env)
                    except (OSError, ValueError) as exc:
                        error(str(exc))
                staged = stage_supplied_bams(args.bams, work)
                comm(
                    "linked %d supplied BAM(s) into %s; they are read in place and never "
                    "deleted" % (len(staged), work)
                )
                ckpt.done("align_reads")
            if args.depth:
                comm("using the depth table given with --depth; no reads are aligned")
            elif args.bams:
                pass  # the alignments are the supplied ones; the depth table is built from them
            elif ckpt.todo("align_reads"):
                _align_reads(args, env, work, assembly_copy)
                ckpt.done("align_reads")
            else:
                comm("skipping read alignment (already done; --resume)")
            # Built outside the metaBAT2 and MaxBin2 checkpoints, so a --resume adding MaxBin2
            # after metaBAT2 reuses it.
            depth = args.depth or os.path.join(work, DEPTH_TABLE)
            if not args.depth and (args.metabat2 or args.maxbin2) and ckpt.todo("depth_table"):
                # Also used with --bams: the table is measured from the supplied alignments.
                depth = _depth_table(args, env, work)
                ckpt.done("depth_table")
            if args.metabat2 and ckpt.todo("metabat2"):
                _run_metabat2(args, env, work, assembly_copy, depth)
                ckpt.done("metabat2")
            if args.maxbin2 and ckpt.todo("maxbin2"):
                _run_maxbin2(args, env, work, assembly_copy, depth)
                ckpt.done("maxbin2")
            if args.concoct and ckpt.todo("concoct"):
                _run_concoct(args, env, work, assembly_copy)
                ckpt.done("concoct")
            if args.semibin2 and ckpt.todo("semibin2"):
                _run_semibin2(args, learned_envs["semibin2"], work, assembly_copy)
                ckpt.done("semibin2")
            if args.vamb and ckpt.todo("vamb"):
                _run_vamb(args, learned_envs["vamb"], work, assembly_copy)
                ckpt.done("vamb")
            if args.comebin and ckpt.todo("comebin"):
                _run_comebin(args, learned_envs["comebin"], work, assembly_copy)
                ckpt.done("comebin")
            if args.binsplit:
                _post_split(args, [b for b in BINNERS if getattr(args, b, False)])
            if args.bams:
                comm("the BAMs given with --bams are not deleted")
            elif not args.keep_bams and not dry_run():
                _discard_alignments(work, ckpt)
            elif args.keep_bams and not dry_run():
                # Index kept BAMs (already sorted), since tools that read them need the index.
                for bam in _bam_paths(work):
                    if not os.path.isfile(bam + ".bai"):
                        run(
                            SAMTOOLS_INDEX.format(threads=args.threads, bam=bam),
                            env=env,
                            tool="samtools index",
                        )
        except ValueError as e:  # read-collection / validation errors
            error(str(e))

        if args.run_checkm:
            for binner_dir in BIN_DIRS:
                path = os.path.join(args.output, binner_dir)
                if _count_bins(path):
                    _checkm.run_checkm(
                        path,
                        threads=args.threads,
                        mem_gb=args.memory,
                        env=checkm1_env,
                        reduced_tree=not args.full_tree,
                        checkm2=args.checkm2,
                        checkm2_env=checkm2_env,
                    )

        announcement("BINNING PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
