"""MetaWrap2 binning: bin contigs with metaBAT2/MaxBin2/CONCOCT/SemiBin2/VAMB/COMEBin.

Reads back to the assembly are aligned to produce coverage, then each selected binner runs.
Run the bin_refinement module afterwards to QC, consolidate, and reassemble the bins.

``--depth <table>`` bins from a contig x sample depth table that already exists and aligns
nothing, which is how a cohort gets multi-sample coverage for one mapping pass per sample rather
than one per sample per assembly. metaBAT2 and MaxBin2 only: the other four binners read the
alignments themselves.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from typing import IO, Dict, List, Tuple

from .. import binsplit as _binsplit
from .. import checkm as _checkm
from .. import depth as _depth
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
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["binning"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# The alignment itself is not here: the aligner-into-samtools-sort pipeline lives in
# metawrap2.modules._common (ALIGN_SORT, SORT_THREADS, align_read_sets), because the cohort coverage
# stage maps reads the same way and two copies of that command would be two places to fix it.
SAMTOOLS_INDEX = "samtools index -@ {threads} -b {bam}"
#: The one depth table. Deliberately without ``--noIntraDepthVariance``, even for MaxBin2, which
#: does not read the variance columns: see :func:`_depth_table` for what that flag actually does.
METABAT2 = (
    "metabat2 -i {assembly} -a {depth} -o {out}/metabat2_bins/bin -m {metabat_len}"
    " -t {threads} --seed {seed} --unbinned"
)
MAXBIN2 = "run_MaxBin.pl -contig {assembly} -markerset {markers} -thread {threads} -min_contig_length {min_len} -out {maxbin_out}/bin -abund_list {abund_list}"
CONCOCT_CUTUP = "cut_up_fasta.py {assembly} -c 10000 --merge_last -b {bed} -o 0"
CONCOCT_COVERAGE = "concoct_coverage_table.py {bed} {bams}"
CONCOCT = (
    "concoct -l {min_len} -t {threads} --coverage_file {cov} --composition_file {comp}"
    " -b {concoct_out} --seed {seed}"
)
CONCOCT_MERGE = "merge_cutup_clustering.py {clustering}"
#: ``single_easy_bin`` is the right subcommand for one contig set and several BAMs, which SemiBin2
#: reads as multi-sample coverage of one assembly. That is what this module has unless --binsplit
#: says otherwise.
#: ``--compression none`` because every other binner here hands over plain fastas and the rest of
#: the pipeline reads them directly; the default gzips each bin.
#: ``--self-supervised`` is already SemiBin2's default, and is passed anyway because the alternative
#: is not merely a different model: semi-supervised training runs MMseqs2 against GTDB, which is not
#: in this env and is a ~100 GB download. Stating it keeps a future change of SemiBin2's default
#: from turning this into a missing-database error.
SEMIBIN2 = (
    "SemiBin2 single_easy_bin --input-fasta {assembly} --input-bam {bams} --output {out}"
    " --threads {threads} --min-len {min_len} --compression none --self-supervised"
    " --random-seed {seed} --tmpdir {tmp} --engine {engine}"
)
#: ``multi_easy_bin`` is the same binner over a *concatenated catalogue* whose contig names carry a
#: sample prefix - what --binsplit asks for. ``-s`` is the separator between the sample label and the
#: contig name (SemiBin2's own default is ``:``; MetaWrap2 writes ``~``, see binsplit.SEPARATOR).
#:
#: The difference from single_easy_bin is where the model is fitted. single_easy_bin fits one
#: sample's contigs with every sample's coverage; multi_easy_bin fits every sample's contigs at
#: once, so the same organism seen in ten samples contributes ten samples' worth of evidence to one
#: clustering decision. That is the whole claim behind multi-split binning.
SEMIBIN2_MULTI = (
    "SemiBin2 multi_easy_bin --input-fasta {assembly} --input-bam {bams} --output {out}"
    " --threads {threads} --min-len {min_len} --compression none --self-supervised"
    " --random-seed {seed} --tmpdir {tmp} --engine {engine} -s '{separator}'"
)
#: Three of these flags are load-bearing in ways the help text does not make obvious.
#:
#: ``--minfasta`` has no default: without it VAMB writes a cluster TSV and *no FASTA files at all*,
#: which from the outside is a binner that succeeded and found nothing.
#:
#: ``-o ''`` turns off binsplitting. VAMB's default separator is ``C``, meant for the multi-sample
#: convention where a contig is named ``<sample>C<contig>``; given contig ids that all happen to
#: contain a capital C mid-string it silently splits every bin on it, and MetaWrap2's assemblies are
#: single-assembly, so there is no sample prefix to split off. VAMB only warns and carries on when
#: *some* id lacks the separator, so this cannot be left to its leniency.
#:
#: ``--seed`` for the same reason SemiBin2 gets one, with VAMB's own caveat that it does not
#: guarantee determinism.
#:
#: ``{cuda}`` expands to ``--cuda`` or to nothing: VAMB's accelerator switch is a boolean, where
#: SemiBin2 takes a word and COMEBin takes a device string. One ``--device`` on this module covers
#: all three (see LEARNED_DEVICES), which is what keeps "run the learned binners on the GPU" a
#: single decision rather than three flags spelled differently.
#: ``-o`` is the binsplit separator, and it is a parameter rather than a constant because the right
#: value depends on what the assembly is. Empty disables binsplitting, which is correct for a single
#: assembly - there is no sample prefix to split on, and VAMB's default ``C`` would shred bins on any
#: contig id containing a capital C. Given a concatenated catalogue (--binsplit) the separator is
#: what makes VAMB do the split itself, which is the mode its own authors report the gain for.
VAMB = (
    "vamb bin default --outdir {out} --fasta {assembly} --bamdir {bamdir}"
    " -m {min_len} -p {threads} --minfasta {minfasta} --seed {seed} -o '{separator}' {cuda}"
)
#: ``run_comebin.sh`` is the only entry point; there is no ``comebin`` on PATH. Like VAMB it takes a
#: BAM *directory* and globs ``*.bam`` out of it, so ``work_files`` is handed over whole.
#:
#: ``-d`` comes from ``--device``. COMEBin's own default is plain ``cuda`` - not SemiBin2's ``auto``
#: but an outright attempt at a GPU, so on a machine without one its default does not degrade, it
#: fails. See DEFAULT_LEARNED_DEVICE for why this module asks for the CPU unless told otherwise.
#:
#: ``-s`` because COMEBin sets *no* seed unless asked: ``run_comebin.sh`` only forwards ``--seed``
#: when the flag is given, so the default is a network and a Leiden clustering that both wander
#: between runs over one input.
COMEBIN = "run_comebin.sh -a {assembly} -o {out} -p {bamdir} -t {threads} -d {device} -s {seed}"

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are not binned
METABAT_MIN_LEN = 1500  # metaBAT2 floors the minimum contig length at 1500 bp
MAXBIN_MARKERS_BACTERIAL = 107
#: Files metaBAT2 writes into its output prefix that are not genomes. Left in the bin set,
#: CheckM and bin_refinement score them as if they were bins.
METABAT2_NON_BINS = ("bin.unbinned.fa", "bin.lowDepth.fa", "bin.tooShort.fa")
#: metaBAT2's ``--seed`` default is 0, and 0 means "pick a random one" - so the binner everything
#: else in this section is compared against was the unreproducible one. Found by accident: three
#: benchmark arms ran metaBAT2 on a byte-identical depth table with a byte-identical command line,
#: and two produced 79 bins while the third produced 78. Nothing about the runs differed, which left
#: the seed as the only candidate. See SEMIBIN2_SEED for the two reasons this matters; both apply
#: here with more force, because metaBAT2 is the default binner and the baseline the others are
#: ranked against.
METABAT2_SEED = 1337
#: CONCOCT's own default is already 1, so this changes nothing today. It is passed anyway for the same
#: reason COVERM_MAPPER is pinned in quant_bins: the knob exists ("0 gives a random seed"), and a tool
#: changing its default underneath a user is a change in their MAGs, not just in their command line.
#: Not 1337 - keeping CONCOCT's own value means this commit cannot be the reason a result moved.
CONCOCT_SEED = 1
MAXBIN_MARKERS_UNIVERSAL = 40  # --universal (better Archaea binning)
#: SemiBin2 seeds a neural net, so without this two runs over identical input disagree about which
#: contigs form a bin. A benchmark comparing it against metaBAT2 then measures that noise as part of
#: the difference, and a user re-running a study cannot reproduce their own MAGs. Any fixed value
#: does; this one is only here to be fixed.
SEMIBIN2_SEED = 1337
#: SemiBin2's own default here is ``auto``, which trains on a GPU whenever torch can see one. This
#: pins it to the CPU, for two reasons that both outlive the benchmark:
#:
#: * A binner whose training silently moves to an accelerator is not measured on the same axis as
#:   metaBAT2, and worse, it would be
#:   *chosen* on hardware most people running a metagenomics pipeline do not have. The machine these
#:   comparisons ran on has several GPUs, which is exactly why this had to be noticed and nailed down.
#: * ``auto`` makes the result depend on the host. With ``--random-seed`` fixed, two CPU runs agree;
#:   a CPU run and a GPU run need not, so the same study re-analysed on a different machine produces
#:   different MAGs for no reason the user can see.
#:
#: ``--device cuda`` is there for anyone who has one and wants the speed over the agreement.
#: See DEFAULT_LEARNED_DEVICE, which is where that decision now lives for all three binners.
#: Where SemiBin2 leaves its bins under its output directory, given the flags in SEMIBIN2.
#: Reclustering is on by default and ``--write-pre-reclustering-bins`` is off, which is the
#: combination that lands here and not in ``output_recluster_bins``/``output_prerecluster_bins``.
SEMIBIN2_BIN_DIR = "output_bins"
#: Smallest bin VAMB is asked to write out, in bases. 200 kb is SemiBin2's own floor and roughly
#: metaBAT2's, so the three agree on what is too small to be a genome and none of them is being
#: credited for bins the others discarded. VAMB has no default here at all - see VAMB.
VAMB_MIN_FASTA = 200000
#: VAMB refuses a minimum contig length below 250 and recommends 2000. The floor is what is enforced
#: here: -l is the user's call, and 2000 would silently discard contigs they asked to bin.
VAMB_MIN_CONTIG_LEN = 250
#: Where VAMB puts the bins it writes, under its own output directory.
VAMB_BIN_DIR = "bins"
#: VAMB's own help says determinism is not guaranteed even with a seed, which is a reason to pass
#: one rather than not to: it removes the variation it can remove. See SEMIBIN2_SEED for why that
#: matters.
VAMB_SEED = 1337
#: Same argument as SEMIBIN2_SEED, except COMEBin passes no seed at all by default, so this is the
#: difference between a reproducible run and an unreproducible one rather than a change of value.
COMEBIN_SEED = 1337
#: What ``--device`` accepts. ``cuda:N`` is accepted too and passed through, which is why this is a
#: prefix check rather than an argparse ``choices`` list.
LEARNED_DEVICES = ("cpu", "cuda", "auto")
#: The default, and the one number in this module that is a policy rather than a measurement.
#:
#: All three learned binners can train on a GPU, and on this machine all three can see eight of them.
#: They still default to the CPU, because a default that silently moves to an accelerator is a
#: default chosen on hardware most people running a metagenomics pipeline do not have, and because a
#: GPU run and a CPU run of the same seeded network need not agree - so the same study re-analysed
#: elsewhere would produce different MAGs for no reason its author could see.
#:
#: That is an argument about the *default*, not about the capability, and the capability is worth
#: having: COMEBin on CPU spent over 280 CPU-hours on three CAMI2 samples without
#: finishing, against 7.1 for metaBAT2 and MaxBin2 together. ``--device auto`` or ``--device cuda``
#: is how you get it, and `benchmark` measures the GPU arms separately so both costs are on record.
DEFAULT_LEARNED_DEVICE = "cpu"
#: Where COMEBin leaves the bin set it finally selects. It writes several intermediate clusterings
#: under ``comebin_res/cluster_res`` and then picks between them on single-copy marker genes, so
#: this is the one directory that holds the answer rather than a candidate.
COMEBIN_BIN_DIR = os.path.join("comebin_res", "comebin_res_bins")
#: The minimum contig length COMEBin bins at. It is not a flag: ``run_comebin.sh`` never forwards
#: ``--contig_len``, so the Python defaults underneath it (1000/1001) are what every run gets. Since
#: ``-l`` is a promise the other binners keep, the assembly handed to COMEBin is pre-filtered when
#: the user asks for more than this - see _comebin_assembly.
COMEBIN_MIN_LEN = 1000
# ───────────────────────────────────────────────────────────────────────────────────────────

#: Every binner this module can run, in the order their bin sets should be handed downstream, and
#: the single place that order is written down. Each name is both the flag (``--<name>``) and the
#: output directory (``<name>_bins``), which is what lets one list serve the module's own CheckM
#: sweep and the `run` driver's -A/-B/-C wiring. Kept here rather than in the driver because adding
#: a binner is a change to this module, and a list the driver owned would be a second place to
#: forget: the three learned binners were added here and not there, so a study sheet selecting
#: ``--semibin2`` had its bins silently dropped before refinement.
BINNERS: Tuple[str, ...] = ("metabat2", "maxbin2", "concoct", "semibin2", "vamb", "comebin")

#: Where each binner's adopted bins land under the module's output directory.
BIN_DIRS: Tuple[str, ...] = tuple("%s_bins" % b for b in BINNERS)

#: Binners that read the alignments themselves rather than a contig x sample depth table, and so
#: cannot be served by ``--depth``. CONCOCT builds its own coverage table from the BAMs with its own
#: script; SemiBin2, VAMB and COMEBin each take the BAM directory and learn from it. Named here
#: rather than checked inline so the error can say which binner is forcing the alignment.
BAM_BINNERS: Tuple[str, ...] = ("concoct", "semibin2", "vamb", "comebin")

#: Which binners to run, as a tier the user picks by how much compute they have. Measured on CAMI2
#: gastrointestinal, three samples, after bin_refinement. "Clean HQ" counts bins that are *both* >=70%
#: complete / <=10% contaminated *and* non-chimeric by GUNC, joined per bin rather than estimated as a
#: product - see bench_work/clean_hq.py for why that distinction mattered:
#:
#:     tier       binners                        HQ   clean HQ   clean/HQ   cpu-hours
#:     very-fast  metaBAT2                        37      35       94.6%       4.8
#:     fast       metaBAT2 + MaxBin2              71      59       83.1%       7.1
#:     default    metaBAT2 + SemiBin2            131     114       87.0%      86.8  (17 on a GPU)
#:     slow       metaBAT2 + SemiBin2 + COMEBin  183     177       96.7%      52.8  (on a GPU)
#:
#: End to end on per-sample assemblies - reads in, MAGs out, which is the axis a user chooses on - the
#: same four tiers return 36, 41, 53 and 60 clean MAGs for 8.8, 13.8, 43.8 and 120.7 cpu-hours of
#: binning and refinement, against metaWRAP 1.3.0's 48 for 52.5. `default` is the shipped tier because
#: it is the first configuration measured here that beats metaWRAP 1 on both axes at once.
#:
#: Three things to know before reading those as a ranking.
#:
#: **Why `slow` is three binners rather than COMEBin alone.** On the gold assembly the trio and COMEBin
#: alone tie - 177 clean MAGs against 176 - which is what the ensemble used to be dropped for. End to
#: end on per-sample assemblies they do not tie: the trio returns 60 clean MAGs against COMEBin's 57.
#: The gold assembly flatters a single strong binner, because there the contigs are long and clean
#: enough that one model can place nearly everything; on a real assembly there is more that one model
#: gets wrong and another does not.
#:
#: Complementarity was measured directly to choose the members (bench_work/binner_complementarity.py:
#: assign every clean high-quality bin to the reference genome owning most of its sequence, then take
#: unions over pairs and trios). The binners are highly redundant here - COMEBin recovers 176 of the 182
#: genomes the best trio reaches, 39 of them found by nothing else - so the union says a second and
#: third binner buy little on this dataset, and the end-to-end number says they buy something anyway.
#: One dataset also cannot say whether COMEBin's 39-genome edge is general or a property of CAMI2's gut
#: community, and a top tier resting on one model has no fallback when it is wrong. Contrastive
#: (COMEBin), self-supervised (SemiBin2) and composition-plus-coverage (metaBAT2) fail differently, and
#: consolidation keeps the best version of each bin.
#:
#: metaBAT2 is the third member rather than VAMB, MaxBin2 or CONCOCT on cost. The best trio by union is
#: COMEBin + SemiBin2 + MaxBin2 at 182 genomes against metaBAT2's trio at 181 - one genome, for which
#: MaxBin2 charges 24.5 cpu-hours against metaBAT2's ~1.5, and brings a chimera tax with it: 37.7% of
#: MaxBin2's bins pass GUNC, so consolidation has more to throw away. metaBAT2 reads the depth table
#: SemiBin2 already built, which is where its ~1.5 comes from, against VAMB's 4.9 on a GPU (67 on a CPU)
#: and CONCOCT's 238. It is also the most widely validated of the six, which is worth something in the
#: member whose job is to be the fallback.
#:
#: The tiers are not supersets of one another and "stepping up can only add genomes" does not hold:
#: `fast` runs MaxBin2, which `default` and `slow` do not. What holds is that each tier recovers more
#: clean MAGs than the one below it.
#:
#: `slow` needs a GPU in practice rather than in principle. COMEBin on CPU reached epoch 20
#: of 200 after 305 cpu-hours and was stopped; on one GPU the same run finishes in 40. It is offered
#: because `--device auto` will find the GPU, not because the CPU path is slow.
#:
#: And `default` asks for the learned-binner env, which is ~1.5 GB of pytorch. That is the price of the
#: jump from 59 clean MAGs to 114, and it is a real price for someone who only wants a quick look -
#: which is what `fast` is for, and why `fast` is exactly what the default used to be.
PRESETS: Dict[str, Tuple[str, ...]] = {
    "very-fast": ("metabat2",),
    "fast": ("metabat2", "maxbin2"),
    "default": ("metabat2", "semibin2"),
    # Three different algorithm families on purpose; see the note above for the margin over COMEBin
    # alone and why metaBAT2 is the third member rather than the one that wins the union by a genome.
    "slow": ("metabat2", "semibin2", "comebin"),
}

#: The tier a pipeline run uses when its sheet names no binner. Named here rather than as a second
#: list of flags in the driver, so the shipped default and the tier table cannot disagree.
#:
#: It is `default` because that tier is the first configuration measured here that beats metaWRAP 1
#: on both axes at once: end to end, on per-sample assemblies, 53 clean high-quality MAGs against
#: metaWRAP 1's 48, for 43.8 cpu-hours against 52.5. `fast` - which this used to be - returns 41 for
#: 13.8, chosen on a quality-per-cost ratio, and losing 15% of the genomes to save compute is not a
#: default anyone asked for.
#:
#: The cost is that a plain run now wants metawrap2-learned, ~1.5 GB of pytorch, so `quickstart` and
#: `install-env --binning-only` build it. Someone who does not want that pays nothing for it by
#: asking for `--preset fast`, which is exactly the old behaviour.
DEFAULT_PRESET = "default"

#: The first three columns of a jgi-format depth table, in order. Checked rather than
#: assumed, because the failure of handing metaBAT2 some other TSV is not an error: it bins whatever
#: the columns happen to contain.
DEPTH_HEADER: Tuple[str, ...] = ("contigName", "contigLen", "totalAvgDepth")


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 binning",
        usage="metawrap2 binning [options] -a assembly.fa -o output_dir reads_1.fastq reads_2.fastq ...",
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
        help="minimum contig length to bin; metaBAT2 floors this at 1500 whatever is asked",
    )
    p.add_argument(
        "--preset",
        choices=tuple(PRESETS),
        help="which binners to run, chosen by how much compute you have (see PRESETS for the measured "
        "basis). Ignored if you name binners yourself",
    )
    p.add_argument("--metabat2", action="store_true", help="bin with metaBAT2")
    p.add_argument("--maxbin2", action="store_true", help="bin with MaxBin2")
    p.add_argument("--concoct", action="store_true", help="bin with CONCOCT")
    p.add_argument(
        "--semibin2",
        action="store_true",
        help="bin with SemiBin2 (needs `metawrap2 install-env learned`, the shared env for the "
        "three learned binners; it brings pytorch, so it is not in the binning env)",
    )
    p.add_argument(
        "--vamb",
        action="store_true",
        help="bin with VAMB (needs `metawrap2 install-env learned`, the shared env for the three "
        "learned binners; it brings pytorch, so it is not in the binning env)",
    )
    p.add_argument(
        "--comebin",
        action="store_true",
        help="bin with COMEBin (needs `metawrap2 install-env learned`, the shared env for the "
        "three learned binners; it brings pytorch, so it is not in the binning env)",
    )
    p.add_argument(
        "--device",
        default=DEFAULT_LEARNED_DEVICE,
        metavar="{cpu,cuda,cuda:N,auto}",
        help="where the learned binners (SemiBin2, VAMB, COMEBin) train. 'auto' uses a GPU if "
        "torch can see one in that binner's env and the CPU otherwise; the default is cpu so that "
        "a run is reproducible on any host (see DEFAULT_LEARNED_DEVICE)",
    )
    p.add_argument(
        "--universal", action="store_true", help="use universal (not bacterial) markers in MaxBin2"
    )
    p.add_argument(
        "--run-checkm",
        action="store_true",
        help="run CheckM on the bins (uses the bin_refinement env; needs plenty of RAM)",
    )
    p.add_argument(
        "--long-read-type",
        choices=sorted(_long.READ_TYPES),
        default=None,
        help="the reads given are long reads of this platform. Sets minimap2's preset and "
        "SemiBin2's sequencing type; mapping long reads with the short-read preset produces "
        "coverage a binner reads as a different organism",
    )
    p.add_argument(
        "--aligner",
        choices=sorted(ALIGNERS),
        default=DEFAULT_ALIGNER,
        help="aligner used to measure coverage (default %s; see ALIGNERS for why)"
        % DEFAULT_ALIGNER,
    )
    p.add_argument(
        "--keep-bams",
        action="store_true",
        help="keep the per-sample BAMs in work_files/ after the depth tables are built "
        "(about 5 GB per sample; discarded by default)",
    )
    p.add_argument(
        "--depth",
        help="bin from a contig x sample depth table that already exists "
        "(the layout metawrap2 coverage and jgi_summarize_bam_contig_depths both write), aligning "
        "nothing. metaBAT2 and MaxBin2 only - "
        "the other binners read the alignments themselves. Reads are then not needed",
    )
    p.add_argument(
        # A required value, not nargs="?". With an optional value this flag sits in front of the
        # positional read files, so `--binsplit r_1.fastq` parsed the read file as the separator and
        # then complained that no reads were given. argparse cannot tell the two apart, so the
        # ambiguity has to go rather than be documented.
        "--binsplit",
        default="",
        metavar="SEP",
        # Concatenated rather than %-formatted. argparse interpolates help strings itself (that is
        # how %(default)s works), so a string that has already been through % leaves a bare "%"
        # behind and argparse then fails with "unsupported format character". One layer of escaping,
        # argparse's, is the only safe number.
        help="treat the assembly as a concatenated multi-sample catalogue (what `metawrap2 "
        "coverage` writes), bin it as one unit, and split the bins back per sample afterwards. SEP "
        "separates the sample label from the contig name; `coverage` writes '"
        + _binsplit.SEPARATOR
        + "'. This is multi-split binning: the model is fitted on every sample's contigs at once "
        "rather than one sample at a time, which is reported to recover 19%% more bins and 79%% "
        "more high-quality bins. Pass the catalogue with -a and its depth table with --depth",
    )
    p.add_argument(
        "--min-identity",
        type=percent(),
        default=_depth.PERCENT_IDENTITY,
        metavar="PERCENT",
        help="minimum alignment identity for a read to count toward a contig's depth. 90 beats jgi's "
        "97 on every measure taken: 88 bins against 74, more pure, and 47 clean high-quality MAGs "
        "against 35 - because on strain-rich data a tight floor removes signal rather than separating "
        "strains. Pass 97 to reproduce a jgi depth table, or 0 for no filter",
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

    ``merge_cutup_clustering.py`` already folds the ``<contig>.<fragment>`` ids that
    ``cut_up_fasta.py`` produced back into whole-contig ids, so the keys here are plain
    contig names and must be used verbatim. The original split_concoct_bins.py truncated
    each name at the *first* '.', which silently mangled any assembler whose contig names
    contain a dot - metaSPAdes emits ``NODE_1_length_5000_cov_3.5`` - collapsing unrelated
    contigs onto one key and scattering them into the wrong bins.
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
            # CONCOCT numbers its clusters from 0; bins are named from 1, like every other
            # binner's, so the cluster id is shifted rather than used as the filename.
            target = (
                bin_filename(int(contig_cluster) + 1)
                if contig_cluster is not None
                else UNBINNED_NAME
            )
            if target not in handles:
                # One handle per bin, all closed in the finally below. Reopening per contig
                # instead would make this O(contigs) in syscalls.
                handles[target] = open(os.path.join(out_dir, target), "w")  # noqa: SIM115
            handles[target].write(">%s\n%s\n" % (contig, seq))
    finally:
        for handle in handles.values():
            handle.close()


def _align_reads(args, env, work, assembly_copy) -> None:
    """Align each read set back to the assembly, producing sorted BAMs in work/.

    The loop itself is :func:`metawrap2.modules._common.align_read_sets`, shared with the cohort
    coverage stage: both map every sample to one reference and want sorted BAMs out, and the only
    difference is which reference that is.
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


#: What `_discard_alignments` sweeps up: the BAMs, and whichever aligner's index is sitting next to
#: the assembly copy. Only the aligner reads an index, and `_align_reads` rebuilds one
#: unconditionally, so they are as dead as the BAMs once the depth tables exist - and nearly as
#: large (1.7 GB for bwa's, 5.3 GB for bwa-mem2's on a 228 Mbp assembly). `assembly.fa` itself
#: stays: the binners are handed its path.
SPENT_SUFFIXES = (".bam", ".bam.bai") + INDEX_SUFFIXES


def _discard_alignments(work: str, ckpt) -> None:
    """Delete the BAMs and the aligner index once every binner has its depth table.

    Nothing downstream of binning reads them: metaBAT2 and MaxBin2 consume a depth table, CONCOCT
    a coverage table, and bin_refinement and the quantifiers work from the bins and the reads. What
    is left is roughly 5 GB per sample, plus the index, whose information content is now a few
    megabytes of TSV - on a ten-sample cohort the module was sitting on ~50 GB for nothing.
    ``--keep-bams`` holds on to them for anyone wanting to inspect an alignment.

    The align_reads marker goes with them, and the depth_table marker with it. A later ``--resume``
    that adds a binner has to realign, and skipping straight to a binner that then found no BAMs
    would hand it an empty depth table rather than an error. The depth table itself survives and is
    still correct for the samples it was built from, but it is cleared anyway: the realignment it
    triggers may be over a different read set, and rebuilding one TSV is nothing against the
    alignment that has already been paid for by the time this matters.
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
    """The BAMs in *work*, sorted. The list form, for callers that are not building a command line."""
    return sorted(os.path.join(work, f) for f in listdir(work) if f.endswith(".bam"))


def _bams(work: str) -> str:
    bams = _bam_paths(work)
    # Under --dry-run no alignment ran, so show a readable placeholder rather than an empty
    # argument list that makes the printed command look malformed.
    if not bams and dry_run():
        return os.path.join(work, "<sample>.bam")
    return " ".join(bams)


#: Name of the shared contig x sample depth table inside ``work_files``.
DEPTH_TABLE = "contig_depth.txt"


def _depth_table(args, env, work) -> str:
    """Build the contig x sample depth table that metaBAT2 and MaxBin2 both bin from.

    One table for both, where there used to be one jgi pass each. That halves the work - jgi reads
    every BAM start to finish and it is the expensive part of this stage - but the reason to do it is
    correctness, not cost.

    MaxBin2's pass asked for ``--noIntraDepthVariance``, which reads as a formatting flag: MaxBin2
    wants one depth per contig per sample and has no use for the variance columns. In
    jgi_summarize_bam_contig_depths 2.18 that code path does not drop columns, it reports numbers
    that are not depths. On a 60-contig two-sample probe it zeroed 59 of the 60 contigs and gave the
    60th a depth of 7867 where the real per-sample depths were 0.267 and 0.453 - and 7867 did not
    reproduce, coming back 7874 on an identical rerun over the same BAMs, so the value comes from
    uninitialized memory rather than from some other definition of depth. Handed a single BAM it
    returns 0 for everything.

    Measured on the real thing, the 210,632-contig CAMI2 GI assembly with three samples, against the
    same BAMs read by the same tool without the flag:

      sample      covered contigs   reported 0   invented coverage   worst overstatement
      sample_0            50,272       10,841                   0                  1.0x
      sample_1            65,375       37,453              24,134          410,067,316x
      sample_10           56,369       32,557              25,187        1,523,050,975x

    So for two of the three samples more than half the contigs that have coverage were reported as
    having none, ~25,000 contigs with no reads at all were given coverage, and the overstatements run
    to depths of hundreds of millions - which in an EM over abundance is not noise, it is the signal.

    Nothing could notice. A table of wrong depths is a well-formed abundance file, and MaxBin2 bins
    on composition as well, so it still returned plausible-looking bins. Taking the mean columns out
    of the full table gives it the depths it was always supposed to have, and costs one jgi pass
    rather than two.

    That bug is also why the table is built here now rather than shelled out to a binary. It is the
    single most load-bearing number in the pipeline and it was the one thing the pipeline could not
    look inside: diagnosing the above meant reasoning about a C++ binary with no test coverage here
    and no way to check an intermediate. `metawrap2.depth` implements the same algorithm - the same
    97% identity filter, the same 75 bp of trimming at each contig end, the same intra-contig variance
    - over per-base depths from ``samtools depth``, and agrees with the tool it replaces to Pearson
    0.999 on the means and 0.995 on the variances over all 210,632 contigs of the CAMI2 GI assembly
    (median relative difference about 5%).
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
    means are every column from the fourth on whose name does not end ``-var``. Tested on the name
    rather than taken as every second column, because the stride only holds while the variance
    columns are there and the point of this function is to not care whether they are.
    """
    return [i for i in range(3, len(header)) if not header[i].endswith("-var")]


def depth_table_contigs(path: str) -> List[str]:
    """The contigs a jgi depth table covers, after confirming the file is one.

    Checked rather than trusted because the failure mode of handing metaBAT2 the wrong TSV is not
    an error message: it reads column 2 as a length and the rest as depths, bins whatever that
    means, and writes a plausible-looking set of bins. A table with no per-sample column is
    rejected for the same reason - metaBAT2 accepts it and bins on composition alone, which is a
    different algorithm than the one the user asked for.
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
                "%s has no per-sample depth column, so there is no coverage to bin on. A jgi "
                "table's fourth column onwards are the samples." % path
            )
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


def _assembly_contigs(path: str) -> List[str]:
    """Contig names in *path*, read from the header lines alone.

    Not via :func:`iter_fasta`, which assembles each sequence in memory to hand it back; nothing
    here looks at a base, and an assembly is the one file in a run big enough for that to matter.
    """
    from ..io.seqio import smart_open

    with smart_open(path, "rt") as fh:
        return [contig_id(line[1:].rstrip("\n")) for line in fh if line.startswith(">")]


def check_depth_table(path: str, assembly: str) -> int:
    """Confirm a precomputed depth table describes *assembly*. Returns the number of samples.

    The dangerous case is not a malformed table but a well-formed one built against a different
    assembly - the cohort coverage stage maps every sample against one concatenated catalogue, so a
    per-sample slice of that table is exactly as valid-looking as the wrong slice would be. No
    overlap at all is an error; a partial overlap is reported, because each side of the gap means
    something different and neither is necessarily wrong:

    * contigs in the assembly but not the table are binned with no coverage at all, which for
      metaBAT2 and MaxBin2 means they are dropped rather than mis-binned
    * contigs in the table but not the assembly are simply ignored, which is what a table built
      over a whole cohort's catalogue looks like when it is sliced correctly
    """
    contigs = depth_table_contigs(path)
    if not contigs:
        error("%s has a header but no contigs, so there is nothing to bin." % path)
    with open(path) as fh:
        samples = len(depth_means(fh.readline().rstrip("\n").split("\t")))
    in_assembly = set(_assembly_contigs(assembly))
    shared = len(in_assembly.intersection(contigs))
    if not shared:
        error(
            "none of the %d contigs in %s appear in %s, so the table describes a different "
            "assembly.\nThe names have to match the assembly's first header token exactly."
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
            "%d row(s) in the depth table name contigs that are not in this assembly and are "
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
    # metaBAT2 writes several things into the output prefix that are NOT genomes: the
    # leftovers from --unbinned, the contigs it rejected as too short or too low-depth, and
    # two plain-text summaries. Left in the bin directory, CheckM and bin_refinement score
    # them as if they were bins. Move them aside (keeping them - they are useful) so only
    # real bins remain.
    bins_dir = os.path.join(args.output, "metabat2_bins")
    if dry_run():
        return
    # metaBAT2 appends "total_depth=.. sample_depths=.." to every contig header. Strip it so
    # these bins carry the same contig ids as the assembly and as every other binner's bins -
    # otherwise nothing downstream can match them up.
    #
    # Note the extension test is FASTA_EXTENSIONS, not BIN_EXTENSION: metaBAT2 writes ".fa" and
    # is renumbered to the canonical naming further down, so at this point the files still carry
    # the tool's own names.
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

    # Only now that the non-bins are out of the way: give the real bins the same names every
    # other binner's get. Done last so a rejected-contigs file can never consume a bin number.
    renumber_bins(bins_dir)


def renumber_bins(bins_dir: str) -> int:
    """Rename every FASTA in *bins_dir* to bin_001.fasta, bin_002.fasta, ... Returns the count.

    Bins are ordered by any number already in their filename, numerically - so metaBAT2's
    ``bin.2.fa`` stays ahead of its ``bin.10.fa`` instead of being reordered by the lexicographic
    accident that makes "10" sort before "2". Files without a number sort last, by name.

    Renaming goes via a temporary name because the source and target namespaces overlap: without
    it, renaming ``bin.1.fa`` to ``bin_001.fasta`` is fine, but a binner that already wrote
    ``bin_002.fasta`` could be overwritten by another bin being renamed onto it.
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
    # One abundance file per sample, each the contig name and that sample's mean depth, taken out of
    # the shared table rather than from a second jgi pass (see _depth_table).
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
        error("MaxBin2 did not produce a single bin. Something went wrong. Exiting.")
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
        # CONCOCT forks ``-t`` workers and does array maths in each, so without this its own
        # parallelism gets multiplied by the BLAS thread count. See THREAD_LIMITS: that multiplication
        # was costing 30 CPU-hours per sample where 5.6 was enough.
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

    Only ``auto`` needs resolving, and resolving it means asking the binner's *own* env, not this
    process: the three learned binners each carry their own pytorch, so "is there a GPU" can have
    different answers in each, and the one that matters is the answer where the training will happen.
    The probe is one short python -c; a failure of any kind resolves to the CPU, because a missing or
    broken torch is exactly the case ``auto`` exists to survive.
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


def _run_semibin2(args, env, work, assembly_copy):
    """Bin with SemiBin2, which reads the BAMs itself rather than a depth table.

    It runs in ``metawrap2-learned`` because it brings pytorch, so *env* here is not the binning
    env every other binner gets.

    The BAMs are handed over directly: SemiBin2 computes both abundance and the k-mer features it
    learns from, and ``--depth-metabat2`` - which would reuse the table metaBAT2 already built, at
    no alignment cost at all - is only accepted alongside a pretrained ``--environment`` model, so
    it cannot be combined with training on this sample. That is SemiBin2's restriction, not a choice
    here; the comparison in benchmarks/configurations.toml is what decides whether the training is
    worth what it costs.
    """
    announcement("RUNNING SEMIBIN2")
    device = resolve_device(args.device, env)
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
        # SemiBin2 spells the accelerator "gpu"; --device speaks cuda like the other two.
        "engine": "gpu" if device.startswith("cuda") else "cpu",
    }
    if separator:
        comm("multi-sample mode: SemiBin2 fits one model over every sample's contigs")
        command = SEMIBIN2_MULTI.format(separator=separator, **fields)
    else:
        command = SEMIBIN2.format(**fields)
    if getattr(args, "long_read_type", None):
        # SemiBin2 distinguishes short from long and nothing finer, and the difference is real: its
        # default feature extraction assumes short-read coverage profiles. The platform name comes
        # from metawrap2.longreads so Flye's flag, minimap2's preset and this cannot disagree.
        command += " --sequencing-type=%s" % _long.resolve(args.long_read_type).semibin2_type
        comm("SemiBin2 in long-read mode")
    run(command, env=env, tool="SemiBin2")

    _adopt_bins(
        os.path.join(out, SEMIBIN2_BIN_DIR), os.path.join(args.output, "semibin2_bins"), "SemiBin2"
    )


#: Binners that do the per-sample split themselves when given a separator, so MetaWrap2 must not do
#: it again. SemiBin2 via multi_easy_bin -s and VAMB via -o; everything else gets the generic split.
NATIVE_BINSPLIT: Tuple[str, ...] = ("semibin2", "vamb")


def _post_split(args, selected: List[str]) -> None:
    """Split each bin set that the binner did not split itself, into <bins>_split/<sample>/.

    The split sets go beside the originals rather than replacing them. Two reasons: the catalogue-wide
    bins are what the binner actually produced and what a chimerism check should see, and a
    consolidation step later wants per-sample sets that line up with per-sample assemblies. Keeping
    both means neither question needs a re-run.
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
        # Also as one flat directory of `<sample>~bin_xxx.fasta`, because the nested layout is the
        # wrong shape for anything that treats a bin set as a folder of FASTA files - including this
        # project's benchmark scorer, which collects top-level `*_bins` directories and would
        # otherwise score the unsplit catalogue bins and report the split arm as identical to its
        # control. Hard-linked, so the second view costs inodes rather than sequence.
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
                "%s: %d contig(s) had no '%s' prefix and were left out of the split. That means the "
                "assembly was not a catalogue built by `coverage`, so --binsplit is binning one "
                "sample's contigs and splitting nothing"
                % (binner, report.unprefixed_contigs, args.binsplit)
            )
        elif not report.worthwhile:
            comm(
                "%s: no bin drew contigs from more than one sample, so this run got nothing from "
                "multi-split that per-sample binning would not also have found" % binner
            )


def _adopt_bins(produced: str, bins_dir: str, tool: str) -> None:
    """Move *tool*'s bins out of its own output tree into *bins_dir* and give them canonical names.

    Both learned binners need this and both fail the same way without it: each writes into a
    subdirectory of its own, under names of its own, and leaves other files beside them. Moving
    rather than copying is deliberate - a bin left in both places is scored twice by anything that
    walks the output.

    Writing nothing is treated as an error rather than as an empty bin set, because both tools exit
    0 when no cluster clears their minimum bin size. Downstream, an empty directory is
    indistinguishable from a sample with no binnable genomes in it.
    """
    ensure_dir(bins_dir)
    if dry_run():
        return
    if not os.path.isdir(produced):
        error(
            "%s finished but wrote no %s. Check its log: with no cluster over its minimum bin size "
            "it exits successfully having written nothing." % (tool, produced)
        )
    for name in sorted(listdir(produced)):
        if name.endswith(FASTA_EXTENSIONS):
            os.replace(os.path.join(produced, name), os.path.join(bins_dir, name))
    n = renumber_bins(bins_dir)
    if n == 0:
        error("%s did not produce a single bin. Something went wrong. Exiting." % tool)
    comm("%s finished successfully, and found %d bins!" % (tool, n))


def _run_vamb(args, env, work, assembly_copy):
    """Bin with VAMB, which reads the BAM directory rather than individual BAMs or a depth table.

    Runs in ``metawrap2-learned`` for the same reason SemiBin2 does: pytorch.

    VAMB refuses to start if its ``--outdir`` already exists, which makes a plain ``--resume`` after
    an interrupted VAMB fail on the leftovers rather than on anything real, so the directory is
    cleared first. That is safe here precisely because the checkpoint has already decided this stage
    has to be redone: nothing in there is a result.
    """
    announcement("RUNNING VAMB")
    device = resolve_device(args.device, env)
    out = os.path.join(work, "vamb_out")
    if not dry_run() and os.path.isdir(out):
        comm("clearing %s; VAMB refuses to run into a directory that already exists" % out)
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
    )
    _adopt_bins(os.path.join(out, VAMB_BIN_DIR), os.path.join(args.output, "vamb_bins"), "VAMB")


def _comebin_assembly(args, work: str, assembly_copy: str) -> str:
    """Return the assembly for COMEBin, length-filtered first if ``-l`` asks for more than it is.

    COMEBin has no ``--contig_len`` on its wrapper (see COMEBIN_MIN_LEN), so a user who raised
    ``-l`` would silently get COMEBin's 1000 while every other binner honoured the number asked
    for - and then a bin set built on contigs the others had excluded, which is exactly the kind of
    difference a comparison between binners must not contain. Filtering the input is the only lever
    available. When ``-l`` is at or below COMEBin's own floor there is nothing to do and the shared
    copy is used as-is, so nothing is written and nothing is duplicated in the common case.
    """
    if args.min_len <= COMEBIN_MIN_LEN:
        return assembly_copy
    filtered = os.path.join(work, "assembly_comebin.fa")
    comm(
        "filtering the assembly to >=%d bp for COMEBin, which has no minimum-length flag of its own"
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

    Runs in ``metawrap2-learned`` for pytorch, like the other two learned binners.

    COMEBin produces several Leiden clusterings at different resolutions, scores each against
    single-copy marker genes, and publishes the winner; only that winner is adopted here. Unlike
    VAMB it is resume-friendly on its own - ``run_comebin.sh`` checks for its intermediate
    manifests and skips the stages that already finished - so the output directory is left alone
    rather than cleared, which matters because the augmentation stage it would otherwise redo is the
    expensive one.
    """
    announcement("RUNNING COMEBIN")
    device = resolve_device(args.device, env)
    out = os.path.join(work, "comebin_out")
    ensure_dir(out)
    comm("Starting binning with COMEBin...")
    run(
        COMEBIN.format(
            assembly=_comebin_assembly(args, work, assembly_copy),
            out=out,
            bamdir=work,
            threads=args.threads,
            device=device,
            seed=COMEBIN_SEED,
        ),
        env=env,
        tool="run_comebin.sh",
    )
    _adopt_bins(
        os.path.join(out, COMEBIN_BIN_DIR), os.path.join(args.output, "comebin_bins"), "COMEBin"
    )


def copy_assembly(assembly: str, destination: str) -> Tuple[int, int]:
    """Copy *assembly* to *destination*, leaving out records nothing can bin.

    Returns ``(empty, uninformative)``.

    **Empty records** have a header and no sequence. Such a record is not a contig, but it is not
    harmless either: bwa indexes it as an ``@SQ`` line with ``LN:0``, and samtools then refuses the
    whole alignment with "Header includes @SQ line ... with invalid or no LN: tag" - so one empty
    record in a million-contig assembly ends the run at the sort, several hours in, with an error
    that names neither the file nor the record.

    **Uninformative records** have sequence but not one A, C, G, T or U in it - in practice a run of
    N as long as a chromosome. Every signal a binner uses is absent from such a record: its
    tetranucleotide frequency is a zero vector, and no read aligns to it, so its coverage is zero in
    every sample. metaBAT2, MaxBin2 and CONCOCT take it anyway and silently carry a contig they
    cannot place; VAMB refuses to start at all, which is the more defensible behaviour of the two and
    is how this was found. Dropping them makes the binners agree and costs nothing: there is no
    assignment of such a contig that could be right.

    Neither case is hypothetical. CAMI2's gold-standard assembly for the gastrointestinal challenge
    ships two empty records, both headed ``_from__to_-42_total_0`` - the tell-tale of a co-ordinate
    bug in whatever wrote it - and 20 all-N records totalling 28.9 Mbp, the largest of them 6.6 Mbp,
    which are whole reference genomes masked out rather than assembly artefacts.
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
    """Turn ``--preset`` into binner flags, or refuse a run that selects no binner at all.

    An explicit selection wins. A preset that quietly added COMEBin to someone's ``--metabat2`` would
    change what a benchmark arm measures and what a user's catalogue contains, and the flag they typed
    is the better evidence of what they meant.
    """
    named = [b for b in BINNERS if getattr(args, b)]
    preset = getattr(args, "preset", None)
    if named and preset:
        comm(
            "--preset %s ignored: you named %s yourself, and an explicit choice should not be "
            "silently widened" % (preset, ", ".join("--" + b for b in named))
        )
        return
    if preset:
        for binner in PRESETS[preset]:
            setattr(args, binner, True)
        comm("--preset %s selects %s" % (preset, " ".join("--" + b for b in PRESETS[preset])))
        return
    if not named:
        error(
            "You must select at least one binning method: --metabat2, --maxbin2, --concoct, "
            "--semibin2, --vamb, --comebin.\n"
            "Or pick a tier with --preset {very-fast,fast,default,slow} and let it choose - see "
            "`--help`."
        )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    _apply_preset(args)

    # --depth replaces the alignment stage, so it also replaces the reads. The two checks are
    # separate because they fail for different reasons: a binner that reads BAMs cannot be served by
    # a table at all, while reads with a table are merely work nobody would use.
    if args.depth:
        needs_bams = [b for b in BAM_BINNERS if getattr(args, b)]
        if needs_bams:
            error(
                "--depth cannot be used with %s: %s read the alignments themselves rather than a "
                "depth table.\nRun those with the reads, or select only --metabat2/--maxbin2 here."
                % (
                    ", ".join("--" + b for b in needs_bams),
                    "they" if len(needs_bams) > 1 else "it",
                )
            )
        if args.reads:
            comm(
                "--depth was given, so the %d read file(s) on the command line are not used: "
                "nothing is aligned" % len(args.reads)
            )
            args.reads = []
    elif not args.reads:
        error(
            "No read files given. Reads are what coverage is measured from; pass them as "
            "positional arguments, or bin from an existing depth table with --depth <table>."
        )

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args, extra=["depth"])
    check_disk_space("binning", settings, [args.assembly] + args.reads, args.output)
    validate_inputs(
        "binning",
        [(check_fasta, args.assembly, "assembly (-a)")]
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    if args.depth:
        # Before the env is resolved and the run recorded: a table that does not describe this
        # assembly is a mistake to catch now, not after a conda env has been probed.
        try:
            check_depth_table(args.depth, args.assembly)
        except (OSError, ValueError) as exc:
            error(str(exc))
    env = env_for("binning", settings)
    # Resolved here rather than at the point of use, because the point of use for a learned binner is
    # after every read has been aligned. A missing env found there costs the whole alignment pass.
    learned_envs = {
        flag: optional_env_for(flag, settings, "--%s" % flag)
        for flag in ("semibin2", "vamb", "comebin")
        if getattr(args, flag)
    }
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
                    "as '@SQ ... LN:0', which samtools refuses to read" % empty
                )
            if uninformative:
                comm(
                    "left out %d record(s) with no A, C, G, T or U in them; they have no k-mer "
                    "signal and no read can align to them, so no binner can place them"
                    % uninformative
                )

        try:
            if args.depth:
                comm("using the depth table given with --depth; no reads are aligned")
            elif ckpt.todo("align_reads"):
                _align_reads(args, env, work, assembly_copy)
                ckpt.done("align_reads")
            else:
                comm("skipping read alignment (already done; --resume)")
            # Built before either binner that reads it and outside both their checkpoints, so that a
            # --resume adding MaxBin2 to a finished metaBAT2 run does not have to rebuild it, and so
            # that neither binner owns a file the other depends on.
            depth = args.depth or os.path.join(work, DEPTH_TABLE)
            if not args.depth and (args.metabat2 or args.maxbin2) and ckpt.todo("depth_table"):
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
            if not args.keep_bams and not dry_run():
                _discard_alignments(work, ckpt)
        except ValueError as e:  # read-collection / validation errors
            error(str(e))

        if args.run_checkm:
            refinement_env = env_for("bin_refinement", settings)
            for binner_dir in BIN_DIRS:
                path = os.path.join(args.output, binner_dir)
                if _count_bins(path):
                    _checkm.run_checkm(
                        path,
                        threads=args.threads,
                        mem_gb=args.memory,
                        env=refinement_env,
                        quick=False,
                    )

        announcement("BINNING PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
