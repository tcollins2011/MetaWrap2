# Modules

MetaWrap2 keeps metaWRAP's modules and adds new ones. `metawrap2 <module> -h` lists every option
for a module; this page covers what each module does, its main options and defaults, and measured
behaviour where it helps choose an option. Comparisons between tools are in
[Benchmarks](benchmarks.md).

Each module runs in a conda environment grouped by workflow; see [Installation](installation.md).
`metawrap2 install-env --all` builds all of them, and `metawrap2 check <module>` verifies a module's
tools inside the environment the run will use.

## Pre-processing

### read_qc

Read trimming and host (e.g. human) read removal.

- `--trimmer` - fastp (default) or `trim-galore`. On the 10 CAMI2 GI samples fastp kept slightly
  more reads, bases and Q30 bases for 5.8% fewer CPU-hours and 29% less wall time; it runs in one pass
  and scales past cutadapt's 4-thread limit. Trim Galore peaks at 0.1 GB of memory, fastp at 1.3 GB.
  See [Benchmarks](benchmarks.md#7-smaller-changes); re-measure with
  `metawrap2 benchmark run cami2_gi trimmer`.
- `--fastp-quality-trim` - adds `--cut_front --cut_tail --cut_mean_quality 20` so that fastp trims.
  Off by default because the benchmark set has no low-quality reads to measure it on.

Retention is reported from the trimmer's own report, e.g. `12,621,104 pairs in ->
5,243,141 out (41.5%): 58.5% removed by quality filtering`, and the run warns when under 70% of reads
are kept.

With default settings fastp does no trimming and discards a read whole when more than 40% of its
bases are below Q15. On one library of 2010-era 75 bp reads this removed 58.5% of the reads, where
trim-galore (which trims the low-quality tail and keeps the shorter read) kept 95%, and `nonpareil`
on the fastp output estimated five times as much additional sequencing. For old or low-quality
libraries use `--trimmer trim-galore` or `--fastp-quality-trim`.

### assembly

Metagenomic assembly with MEGAHIT (default) or metaSPAdes, plus QC.

- `--metaspades` - on three CAMI2 GI samples co-assembled, a 38% longer N50 (12.9 kb vs 9.3 kb) for
  3.2x the CPU-hours, 9.7x the peak memory (140 GB vs 15 GB) and 14% fewer bases assembled. Use it
  when contiguity matters and the memory is available. Re-measure with
  `metawrap2 benchmark run cami2_gi assembler`.
- `--megahit --metaspades` - hybrid mode: MEGAHIT assembles the reads that did not map to the
  metaSPAdes scaffolds, and the two are merged. The alignment is piped into the filter (nothing
  written to disk) and the scaffold index is deleted once the FASTQ is written.
- `--aligner` - the aligner for hybrid mode, from the registry `binning` and `blobology` use.
  Default bwa. On a 273k-scaffold metaSPAdes assembly and its 3-sample reads, bwa took 5.45
  CPU-hours and passed 2,867,680 unmapped reads to MEGAHIT; minimap2 took 0.62 and passed 4,113,486
  (44% more; 1,249,798 reads mapped by bwa were unmapped by minimap2, 3,992 the other way, Jaccard
  0.695). MEGAHIT assembles the extra reads a second time, and the merge cannot tell a duplicate
  contig from a novel one. (In `binning` the aligner only feeds a depth estimate, so minimap2 is the
  default there.)
- `--keep-graph` - keeps each assembler's graph as `assembly_graph.<assembler>.{gfa,fastg}`, for
  plasmid-aware binning, graph-aware binners and Hi-C placement. Off by default. metaSPAdes and
  metaFlye GFAs are hard-linked (121 MB and 2.3 MB measured, no extra disk). MEGAHIT's is built with
  `megahit_toolkit contig2fastg` (235 MB, two seconds on one core) at the largest k MEGAHIT wrote,
  which can differ from the configured one. It is a separate resumable step, so
  `--resume --keep-graph` adds a graph to a finished assembly. Without the flag nothing is deleted.

Hybrid mode cannot be smoke-tested on the `tiny` fixture: metaSPAdes' error correction exits with
*Failed to determine offset!* on its synthetic quality strings. The remapping stage runs on the
fixture when given scaffolds.

### assembly (long reads)

`-L reads.fastq --read-type {ont,ont-hq,pacbio,pacbio-hifi}` assembles with metaFlye. Several `-L`
files are co-assembled as one dataset. `--meta` is always passed.

The platform setting (`metawrap2.longreads`) maps to Flye's `--nano-hq`, minimap2's `-x map-ont` and
SemiBin2's `--sequencing-type=long_read`. The reads are checked against the declared platform
(this flags a short-read FASTQ passed to `-L`, or HiFi declared for nanopore-accuracy data), and
N50, median, the fraction over 10 kb and the count under 1 kb go to `long_read_report.tsv`.

On 25× simulated R10 nanopore reads over three CAMI2 MAGs (2.50 Mbp): 72 contigs, 2.41 Mbp, N50
57.9 kb. Bin with `binning --long-read-type ont-hq`, which maps with `map-ont`; with the short-read
preset coverage is low and uneven. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#85-long-reads).

Not implemented: short-read polishing (unmeasured risk of changing one strain's consensus to a
relative's), and merging long- and short-read assemblies (the duplicate-contig problem above).

### kraken2

Taxonomy profiling of reads or contigs, with a kronagram.

- `--bracken` - adds a Bracken abundance estimate at one rank (`--bracken-rank`, default species).
  kraken2 counts reads at the lowest node all their k-mers agree on, so a species' reads are split
  between its own node and its genus; Bracken redistributes them downward.
- `--bracken-db` - default: the kraken2 database directory, where `bracken-build` writes.
- `--bracken-read-length` - must match a `databaseNmers.kmer_distrib` file in that directory. The
  standard kraken2 database ships 50, 75, 100, 150, 200, 250 and 300 bp.

Bracken files must be built by `bracken-build` from the same kraken2 database. A run takes seconds
per sample on one core; building costs one kraken2 classification of the reference set per read
length. Outputs: a `.bracken` table per sample and `bracken_abundance.tsv`, in the long format
`catalogue` writes and `diff_abundance` reads. Assemblies are skipped with a warning. Without
`--bracken` the kraken2 output is unchanged.

### nonpareil

Diversity coverage: how much of each community the reads cover, and whether samples can be
compared.

Inputs: forward reads only (`-r`), QC'd (Nonpareil assumes error well below 5%).

Per sample it reports C (abundance-weighted average coverage), Nd (the Nonpareil diversity index,
comparable across samples) and the effort to reach `--target-coverage` (default 95) as a multiple of
what was sequenced. Model fitting is a Python port of Nonpareil's R package (`metawrap2.nonpareil`)
that reproduces R's C, Nd, LRstar, model r and both gamma parameters to six significant figures on
Nonpareil's example datasets.

- `--kernel {auto,kmer,alignment}` - `auto` (default) uses the faster k-mer kernel unless a sample's
  mean read length is under twice k, in which case the whole cohort uses the alignment kernel. On
  ERR011347 (mean 34 bp; 11 k-mers per read at k=24) the k-mer kernel gave 50.5% coverage and 47x
  more sequencing, the alignment kernel 83.0% and 4.4x; about 4 and 16 minutes per sample.
- `--seed` (default 1337) - passed to the alignment kernel as `-r`. Unseeded, two runs gave 83.0% and
  84.3%. The k-mer kernel has no seed, so its subsampling is not reproducible.
- `--no-filter` - skip removing reads shorter than k, which the k-mer kernel refuses (3.4% of one
  real sample). The filtered count is reported.

Under the k-mer kernel effort is in k-mers (`AL` is L-k+1 per read); each row gives `effort_unit` and
the sequenced base pairs.

On the ERR011347-9 test samples QC'd with `--trimmer trim-galore`
(32-49 bp, alignment kernel): coverage 84.4%, 74.6% and 89.4%, Nd 17.3-18.6, 2.5x to 8x more
sequencing to reach 95%. With fastp's defaults the 75 bp reads survive, `auto` picks the k-mer
kernel, coverage is near 0.65 and the effort projection three times higher; Nd stays at 17.3-18.7.
Report the QC and kernel with the coverage. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#5-diversity-coverage).

## Bins

### binning

Initial binning with metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB and/or COMEBin. `metawrap2 binning`
has no default binner and needs at least one binner flag or a preset; `metawrap2 run` uses metaBAT2
and SemiBin2 (as of 2026-10-09), which needs the `metawrap2-learned` environment (~1.5 GB of
pytorch, built by `quickstart`).

`--preset {very-fast,fast,default,gpu}` picks a binner set by compute. Clean high-quality MAGs
(≥70/≤10 and non-chimeric by GUNC) from binning CAMI2 GI's gold assembly:

| tier | binners | HQ | clean HQ | clean/HQ | CPU-h |
|---|---|---|---|---|---|
| `very-fast` | metaBAT2 | 37 | 35 | 94.6% | 4.8 |
| `fast` | metaBAT2 + MaxBin2 | 71 | 59 | 83.1% | 7.1 |
| `default` | metaBAT2 + SemiBin2 | 131 | 114 | 87.0% | 86.8 (17 on a GPU) |
| `gpu` | metaBAT2 + SemiBin2 + COMEBin | 183 | 177 | 96.7% | 52.8 on a GPU |

The tiers are not nested, and a preset never overrides a binner you named. On real per-sample
assemblies the tiers are closer together; end to end, `default` gives 53 clean MAGs against metaWRAP
1's 48 for 43.8 CPU-hours against 52.5, and `gpu` 60 against COMEBin alone's 57. `very-fast`'s 4.8
CPU-hours is with minimap2 (22.9 under bwa). See [Benchmarks](benchmarks.md#3-binning) for both
designs and the choice of third binner.

Binner flags:

- `--metabat2`, `--maxbin2` - read the shared depth table `work_files/contig_depth.txt`. MaxBin2's
  per-sample abundance files take its mean columns. `--maxbin2-max-iteration` sets the EM iteration
  limit (default 50, MaxBin2's own); on a 418k-contig co-assembly of 36 samples one iteration took
  about 53 minutes, so 50 take about 40 hours. Preset costs in Benchmarks (3 samples, 17-63k
  contigs) do not scale to that size.
- `--concoct` - opt-in (since 2026-10-02). On `cami2_gi` it used 192 of the 213 binning CPU-hours of
  a three-binner run, for 74 clean MAGs against the classical pair's 64. Coverage comes from
  `concoct_coverage_table.py`; its depths run ~1.73x the jgi table's because `samtools bedcov`
  applies no identity filter.
- `--semibin2` - reads the BAMs (`--input-bam`). `--random-seed` is fixed. Needs contigs of at least
  4000 bp and drops bins under 200 kb regardless of `--min-len`.
- `--vamb` - variational autoencoder over composition and abundance (`--bamdir work_files`). Set
  away from VAMB's defaults: `--minfasta` 200 kb (otherwise no FASTA files), `-o ''` (no binsplitting
  on the `C` separator) and `-m 250` (so contigs `-l` asked for are not discarded). `--seed` is
  fixed, though VAMB does not guarantee determinism.
- `--comebin` - contrastive learning over six augmented views, Leiden clustering at several
  resolutions, chosen by single-copy markers. Runs `run_comebin.sh`; its env includes FragGeneScan,
  HMMER and pplacer. `-d` follows `--device`, `-s` is passed, and `-l` is applied by pre-filtering
  the assembly when `--min-len` is above COMEBin's built-in 1000 bp. Completed stages are kept across
  `--resume`. On CPU, over three CAMI2 samples, it reached epoch 20 of 200 after 305 CPU-hours.

SemiBin2, VAMB and COMEBin share `metawrap2-learned` (`metawrap2 install-env learned`), which ships a
CUDA-capable pytorch. All three apply a 200 kb minimum bin size. Benchmark arms: `semibin2`,
`metabat2_semibin2`, `vamb`, `metabat2_vamb`, `comebin`, `metabat2_comebin`, and the `semibin2_gpu`
/ `vamb_gpu` / `comebin_gpu` GPU arms.

Other options:

- `--device {cpu,cuda,cuda:N,auto}` - where the learned binners train. Default cpu, because GPU and
  CPU runs of the same seeded network can give different MAGs. `auto` uses a GPU if that binner's
  torch sees one and falls back to the CPU otherwise. `--preset gpu` defaults to cuda; an explicit
  `--device` overrides it. Peak VRAM is recorded: COMEBin 5.71 GB (17k and 63k contigs), SemiBin2
  0.41 GB on 16,936 contigs and 2.75 GB on 63,376, VAMB 0.56 GB on 16,527, so an 8 GB card runs all
  three on that data. See [Hardware](https://github.com/ursky/MetaWrap2#hardware).

  `cuda:N` pins the step to card N with `CUDA_VISIBLE_DEVICES`. Before each GPU step the card is
  checked with `nvidia-smi` and a line gives the card, its utilisation, its free memory and any
  process holding it. An index above the device count ends the run (the pin would hide every card,
  and the binner would fail at its first allocation or train on the CPU). A card in use (another
  process on it, or 25% utilisation or more) gives a warning and the step shares it. Bare `--device cuda` pins nothing: the binner's framework picks a card, usually
  cuda:0, and does not fall back to the CPU when every card is busy; the preflight lists which
  cards are free. With no `nvidia-smi`, or a driver that does not answer, the preflight says so and
  the step starts. `--device cpu` probes nothing, and `--dry-run` skips the check.
- `--gpu-strict` (also `--strict-gpu`) - refuse to start a learned binner on a card another
  process is using, instead of warning and sharing it. Use it for benchmark runs: a shared card
  inflates the step's wall time and its GPU trace includes the other job. Same as the global
  `metawrap2 --gpu-strict` or `[limits] gpu_strict = true`.
- `--min-identity PERCENT` - minimum alignment identity for a read to count toward depth (default
  90). `jgi_summarize_bam_contig_depths` uses 97, which on CAMI2 GI discards half the mapped bases.
  `--min-identity 97` reproduces a jgi table (for comparison with metaWRAP 1 or metaBAT2's published
  validation); `0` disables the filter. One set of alignments, metaBAT2 on each table:

  | floor | mean depth | bins | AMBER F1 | purity | HQ | clean HQ |
  |---|---|---|---|---|---|---|
  | none | 15.71 | 90 | 0.6607 | 0.8374 | 48 | 45 |
  | ≥90 (default) | 15.70 | 88 | 0.6624 | 0.8375 | 48 | 47 |
  | ≥95 | 13.46 | 86 | 0.6574 | 0.8264 | 43 | 40 |
  | ≥97 *(jgi's)* | 8.35 | 74 | 0.6371 | 0.8175 | 37 | 35 |
  | ≥99 | 1.48 | 21 | 0.2452 | 0.7904 | - | - |

  Measured on one dataset with one binner; see [Benchmarks](benchmarks.md#4-coverage-estimation).
- `--aligner {minimap2,bwa,bwa-mem2}` - default minimap2. On one CAMI2 GI sample, minimap2 and bwa
  gave bins of the same quality (AMBER F1 0.4519 against 0.4493) for 1.37 CPU-hours against 7.77.
  bwa-mem2 matches bwa's depths exactly at 1.5x its speed and needs 23.6 GB of RAM against 9.8.
  This was measured on CAMI2's gold assembly; aligners differ more on fragmented assemblies.
- `--depth <table>` - bin from an existing depth table, with no reads or alignment (metaBAT2 and
  MaxBin2 only; other binners are refused by name). The table is checked against the assembly: no
  overlap is an error, a partial overlap is reported both ways (contigs with no row are left
  unbinned; extra rows are expected from a cohort slice). See [coverage](#coverage).
- `--bams DIR` - bin from existing BAMs (from `coverage --split-bams`, or kept by
  `coverage --keep-bams` / `binning --keep-bams`), so alignment and binning can run on different
  machines. BAMs are symlinked, never deleted, and their headers are checked against the assembly.
  With `--depth` as well, the table serves metaBAT2 and MaxBin2, the BAMs the rest, and nothing is
  aligned.
- `--keep-bams` - keep the per-sample BAMs, `.bai` indexes and the aligner's index. By default these
  are deleted from `work_files/` once every binner has its depth (on one real sample, 0.6 GB of 1.0 GB;
  in a cohort about 5 GB of BAM per sample; minimap2's `.mmi` for a 154 MB assembly is 479 MB).
  Deleting them clears the `align_reads` checkpoint, so a later `--resume` that adds a binner
  realigns.

Depth table: `metawrap2.depth` implements `jgi_summarize_bam_contig_depths`'s algorithm (identity
filter, 75 bp trimmed off each contig end, intra-contig variance) over `samtools depth`, and agrees
with jgi to Pearson 0.9992 on means and 0.995 on variances across 210,632 CAMI2 GI contigs. MaxBin2
results from multi-sample runs before this table was used (including the `metabat2_maxbin2` and
`metabat2_maxbin2_concoct` benchmark arms) were binned on invalid values; on `cami2_gi` the fix raised
precision 0.6485 → 0.8016, ARI 0.7079 → 0.8747 and F1 0.4228 → 0.4598. See
[Benchmarks](benchmarks.md#8-defects-found-and-fixed).

Contigs with no A, C, G, T or U are dropped (VAMB refuses them; CAMI2's gold GI assembly has 20,
28.9 Mbp).

None of the learned binners can be smoke-tested end to end on the `tiny` fixture: SemiBin2 finds no
4000 bp contig, VAMB's 160 bins are all under 200 kb, and COMEBin trains (about four minutes on 8 CPU
threads) then fails with *HMMsearch completed but found zero marker hits*.

### coverage

One cross-sample depth table for a cohort, in N alignments instead of N² (20 instead of 400 for
twenty samples). The per-sample assemblies are concatenated into one catalogue, each sample's reads
are mapped to it once, and the table is split into one slice per assembly, each keeping every sample
column:

```
metawrap2 coverage -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
    -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
metawrap2 binning -a s0/final_assembly.fasta --depth cohort_coverage/depth/s0.txt \
    -o s0/binning --metabat2 --maxbin2
```

Outputs: `<out>/depth/<label>.txt` per assembly, with the original contig names.

- `-a label=path` - label each assembly (otherwise `sample_N/final_assembly.fasta` would label every
  sample `final_assembly`). Duplicate labels, and labels that match no contig, are errors.
- `--separator` - character used to prefix contigs with their label (`s0~k141_1`); default `~`,
  which SPAdes, MEGAHIT and Flye never emit. Names already containing it are refused.
- `--keep-bams`, `--keep-catalogue` - keep the alignments and catalogue (deleted by default). Kept
  BAMs are indexed. With one assembly their headers use plain contig names; with several they keep
  the catalogue names, and the run prints the `samtools` command to extract one assembly.
- `--split-bams` - also write `<out>/bams/<label>/<sample>.bam`, one indexed BAM per assembly per
  sample with that assembly's contig names, for `binning --bams`. `cohort` turns this on when a
  binner reads BAMs (the default pair does, through SemiBin2). The split costs about one pass over
  the catalogue BAMs. Other assemblies' `@SQ` lines stay, prefixed and empty. The catalogue is the
  size of all assemblies combined, and minimap2's index of it is larger.
- `--coverage fairy` - estimate depth from k-mer sketches instead of alignment, in minutes instead of
  hours. Same table format (`contigName contigLen totalAvgDepth <sample> <sample>-var ...`); the
  `-var` columns come from containment. Default `-m 95` resolves nothing below about 95% ANI.
  Sketches (`sketches/<sample>.paired.bcsp`) are kept for reuse against another catalogue;
  interleaved input is sketched as single-end. On one CAMI2 GI sample through metaBAT2, fairy scored
  AMBER F1 0.212 against minimap2's 0.452 for 0.012 CPU-hours against 1.368. Cohort comparison:
  [Benchmarks](benchmarks.md#1-cohort-support).

  ```
  metawrap2 install-env bins
  metawrap2 coverage --coverage fairy -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
      -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
  ```

Environment: binning's, `metawrap2-bins` (`metawrap2 install-env coverage` builds or finds it),
which also holds fairy.

### bin_refinement

Consolidates multiple bin sets into one refined set, written to `metawrap_<c>_<x>_bins` (named for
the `-c`/`-x` thresholds).

- `--consolidator` - how the sets are combined:
    - `metawrap` (default, metaWRAP's method) hybridises every combination of the input sets with
      Binning_refiner and asks CheckM2 which version of each bin to keep: CheckM2 on 2^n-1 sets. On the
      `metabat2_maxbin2` arm of `cami2_gi`: 10.86 CPU-hours, 70 final bins at `-c 70 -x 10`.
    - `dastool` scores input bins once on single-copy genes and aggregates in one pass. Needs
      `-a/--assembly` and its own env (`metawrap2 install-env binqc`; no database).
      `--dastool-score` is its quality floor (default 0.5); `work_files/dastool_bin_names.tsv` maps
      renumbered bins to DAS_Tool's names, which show the input set.

  Both feed the same CheckM re-scoring, `-c`/`-x` cut, plots, layout and GUNC. The `consolidator`
  comparison (`metawrap2 benchmark run cami2_gi consolidator`) prefers DAS_Tool, which is 7× cheaper,
  on a tie. See [Benchmarks](benchmarks.md#5-refinement).
- `--skip-gunc` - GUNC scores the final bins for chimerism and adds `gunc_css`/`gunc_pass` to their
  `.stats.tsv`; bins are annotated, never dropped. Skipping avoids its 6.7 GB database and 3.20
  CPU-hours on a 60-bin set (~32% on top of refinement). On the `metabat2_concoct` arm, 4 of 60 final
  bins were chimeric by GUNC but passed CheckM at `-c 70 -x 10`.
- `--score-everything` - with `metawrap`, bins outside 200 kb to 20 Mb are recorded as 0% complete /
  100% contaminated without running CheckM (the 10,738 of 18,490 scored bins passing 70/10 range
  from 0.99 to 6.66 Mb). This flag scores them anyway. Bins identical to one already scored (390 of
  1,230 in a three-binner run) reuse the cached score in `work_files`, also across `--resume`.
- `--checkm-jobs N` - with `metawrap`, how many bin sets are scored at once (default 0: as many as
  memory allows, at 24 GB per CheckM2 call or, with `--checkm1`, 16 GB per CheckM1 call on the
  reduced tree and 40 GB with `--full-tree`, and at most one per thread).
  Threads and memory are split evenly between the calls. pplacer's memory is per process, ~34 GB for
  the full tree, and it gains little past 8 threads, so sets run in parallel instead of giving one
  call more threads; pplacer is capped at 8 threads. Each bin is assigned to one set before scoring
  starts, so a shared bin is scored with the same group of bins as in a one-at-a-time run. With
  CheckM1 the scores are identical; with CheckM2 one can differ by a few hundredths of a point,
  because DIAMOND's hits for a protein depend on the rest of the query file. `--checkm-jobs 1`
  scores the sets one at a time. Timings are in [benchmarks](benchmarks.md#26-refinement-speed).
- `--checkm1` - score with CheckM1 instead of CheckM2, the default. CheckM2 needs less memory than
  CheckM1 on its full genome tree (pplacer holds about 34 GB per call), is faster, and is what all
  MetaWrap2 benchmarks used. CheckM1 needs `metawrap2 install-db checkm`, which is checked at the
  start. It runs on CheckM1's reduced tree (`--reduced_tree`); `--full-tree` uses the full genome
  tree, which can choose a different marker set. `--checkm2` is accepted and changes nothing, and
  `--quick` (the reduced tree before it became the default) prints a note. A `--resume` must use the
  scorer and tree the run started with.
- `--checkm2-ttable {4,11}` - run CheckM2 under one translation table. By default CheckM2 calls genes
  under both and keeps the one with the higher coding density. Table 4 won for 8 of 6,778 bins in 45
  CheckM2 reports, so `--checkm2-ttable 11` halves gene calling and mis-calls the rare genome that
  uses table 4.

Every input set is first filtered to bins between 50 kb and 20 Mb, and the count and share removed at
each end are reported. This filter accounts for much of the recall lost in refinement. On `cami2_gi`
the upper limit removed 7 of CONCOCT's 99 bins holding 41% of its sequence, and 3 of metaBAT2's 78
holding 39%.

### reassemble_bins

Reassembles bins to improve N50/completion and reduce contamination. Three candidates per bin are
scored with CheckM2 (strict reassembly, permissive reassembly, original bin), and the summary reports
how many bins kept their original. `reassemble_bins` has no `--checkm1`. On one real 26-bin set, 4 bins (15%) scored better before
reassembly. Before 2.1.3 the original was never scored, so a reassembly was always returned.

- `--max-coverage` (default 50, `0` to disable) - subsamples each bin's recruited reads to this
  depth, relative to the bin's assembled length, before SPAdes. SPAdes time grows with read count,
  and 50x is above the depth at which prokaryotic short-read assembly stops improving. Bins below the
  target are left as is. Reads are dropped in pairs with a fixed seed.

### dereplicate

Collapses several samples' bins into one non-redundant genome catalogue.

```
metawrap2 dereplicate -b sampleA/metawrap_70_10_bins -b sampleB/... -o catalogue/
```

Outputs: `dereplicated_genomes/` (one representative per cluster, symlinked; `--copy` to copy) and
`clusters.tsv` (one row per input bin with its cluster, source folder, representative, quality
score, size and `ani_to_representative`, empty for members added by chaining).

- `--ani` (default 95) - the species threshold large MAG catalogues use; `--ani 99` separates
  strains.
- `--min-aligned-fraction` (default 10%) - share of a genome a pair must align. Raising it gives
  smaller clusters and a larger catalogue.

Identity comes from skani (`triangle -E`). Clustering is single linkage, so A-B and B-C at 96% merge
all three even when A and C are not. The representative is the member with the best completeness
minus five times contamination from the `<bins>.stats.tsv` beside each folder, size breaking ties;
scored genomes beat unscored ones. Without a CheckM table the module warns and picks by size, which
can choose a chimera; run `bin_refinement` or `binning --run-checkm` first. Install with
`metawrap2 install-env dereplicate`.

### quant_bins

Bin abundance across samples. All quantifiers write `bin_abundance_table.tsv` and a heatmap, but the
values differ in meaning, so do not compare tables from two quantifiers.

| `--quantifier` | value | a zero means |
|---|---|---|
| `sylph` (default) | k-mer sequence abundance, % | not detected above `--min-ani` |
| `salmon` | length-weighted median contig coverage, no ceiling | |
| `coverm` | relative abundance, % of sample (row sums to at most 100) | mapped, but under `--min-covered-fraction` covered |
| `depth` | length-weighted mean depth per bin from `--depth TABLE`, with the per-contig median | |

On 76 CAMI2 GI gold-assembly bins across three samples: sylph ranked bins best (Spearman 0.868 vs
salmon 0.827, CoverM 0.803) for 0.126 CPU-hours (salmon 1.08, CoverM 2.37) and 1.3 GB peak memory
(CoverM 11.6 GB); salmon fit magnitude better (Pearson 0.959 vs 0.811). `depth` aligns nothing, costs
0.013 CPU-hours, and was second on magnitude (Pearson 0.953, L1 similarity 0.799) and last on ordering
(Spearman 0.780); use it when magnitudes matter and binning's table exists. See
[Benchmarks](benchmarks.md#6-bin-abundance); re-measure with
`metawrap2 benchmark run cami2_gi quantifier`.

- `--min-ani` (sylph only, default 99) - at sylph's default of 95, more than half the abundance in
  one sample went to four bins absent from it; at 99.9 present bins start to be discarded. Lowering
  it towards 95 is reasonable for bins from fragmented real assemblies. Below 95 MetaWrap2 warns.

sylph and salmon are in the module's env; CoverM needs `metawrap2 install-env bins`. Gzipped bins
are accepted and non-FASTA files ignored, so `-b` can point at `bin_refinement` output.

### blobology

GC-vs-abundance blobplots of contigs and bins.

- `--aligner` - binning's three choices, default minimap2. bowtie2 is no longer used (7.44
  CPU-hours, as costly as bwa, and less close to bwa than minimap2).

The alignment is streamed unsorted with no SAM on disk; BAMs and the index are removed at the end.

### classify_bins

Conservative taxonomy assignment for bins.

- Default: MEGABLAST against NCBI `nt` and taxator-tk. Deprecated (taxator-tk has had no release
  since 2015), but kept as default until compared with GTDB-Tk with both databases installed.
- `--gtdbtk` - maintained, gives GTDB taxonomy.

Both are expensive: `nt` is hundreds of GB, GTDB ~110 GB, and `classify_wf` needs a lot of RAM.

### annotate_bins

Deprecated. Use `annotate`.

## Function

Added in 2.1.2; the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md) give the
measurements behind each one.

### annotate

Calls genes once over a bin folder (`-b`) or an assembly (`-a`) and annotates them. Replaces
`annotate_bins`.

Genes are called with Pyrodigal. Each bin trains its own model (`--meta` forces metagenomic mode); an
assembly uses the metagenomic models. Every gene has `partial` and `distance_to_contig_end` (on CAMI2
bins, 8% of genes are partial in a 92-contig bin against 50-54% in bins of ~900 contigs).

- `--kofam` - KEGG orthologs with PyHMMER at KOfam's per-profile thresholds (17,890 on full-sequence
  score, 8,697 on best domain); the 1,866 profiles with no threshold are dropped by default.
- `--pfam` - Pfam domains at each family's gathering threshold. On one CAMI2 bin KEGG placed 1,620 of
  3,891 genes (41.6%) and Pfam 3,185 (81.9%), 1,578 with no KEGG ortholog; the run prints this
  increment. Output is per domain; `pfam_architectures.tsv` gives each gene's ordered domains
  (`SusD-like_3 + SusD_RagB` was most common, 41 genes). Rows carry the versioned accession and the
  bare family, which joins across releases. See the
  [2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#86-pfam-protein-domains).
- `--eggnog` - eggNOG-mapper: GO terms, EC numbers and COG categories, transferred from the
  orthologous group; `max_annot_lvl` gives the group's taxonomic level. On five bins (14,171 genes)
  it placed 94.9% in a group and gave 12,634 genes 41,075 terms. Of KO calls, 5,077 genes were called
  by both eggNOG and KOfam (4,841 the same KO), 2,295 by eggNOG only, 750 by KOfam only; both tables
  are written. Needs `metawrap2 install-db eggnog` (48 GB) and its own env (eggNOG-mapper does not
  run on Python 3.12 or later).
- `--cazy` - CAZyme families with dbCAN. Each row records how many of dbCAN's three methods agreed
  (on two CAMI2 bins, 415 CAZyme genes, 274 with at least two, dbCAN's own threshold). The per-class
  summary lists CBM apart from the catalytic classes. Needs `metawrap2 install-db dbcan` (7.4 GB,
  a pinned release).

InterProScan is not wrapped (about 50 GB of data and a four-year-old bioconda build); Pfam-A (399 MB)
is searched with PyHMMER instead.

### pathways

KEGG module completeness per genome, plus the functional hierarchy and Krona plots.

Each row has two completeness scores, `pathwise` (best alternative route) and `stepwise` (fraction
of top-level steps), and seven reliability flags (number of top-level steps, signature module,
reaction with no KO, and others).

- `--brite-db` (or config key `KEGG_BRITE`) - the KEGG BRITE ko00001 hierarchy from
  `install-db kegg`. Without it, module completeness is written and the hierarchy and Krona plots
  are skipped.
- `--skip-krona` - omit the Krona plots (`ktImportText`, installed for kraken2).

The hierarchy is gene -> chemical process -> pathway, per genome (genes weighted equally) and per
sample (genes weighted by abundance). The per-sample weight:

| flag | weights each gene by | covers | cost |
|---|---|---|---|
| `-a/--abundance` | its genome's abundance | the binned community | free, from `quant_bins` |
| `-C/--contig-depth` | its contig's depth | the whole assembly | free, from `binning`/`coverage` |
| `-G/--gene-abundance` | its own abundance | the whole assembly | one extra mapping pass |

Use `-a` when MAGs are the focus and `-C` without bins (the depth table covers the ~90% of contigs
never binned). `-G` comes from `catalogue --quantify` or `quantify_function`. A `weighting` column
names which produced each row. On the dev cohort the gene-weighted profile put 33% of the signal in
"Not Included in Pathway or Brite" against ~14% MAG-weighted. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#1-functional-annotation-and-pathways).

### report

One PDF for a finished run, and one row per genome. Reads the module output directories it is given
and computes nothing new: summary, read QC, assembly statistics, the per-genome table and its
quality tiers, abundance, `pathways` Krona plots, and row counts for resistance, viruses and strains.

- `--mag-table` - the per-genome table written beside the PDF, joining `.stats.tsv`, GUNC,
  `classify_bins` taxonomy and `quant_bins` abundance on the bin name across the spellings
  `bin_001`, `bin_001.fasta` and `bin.1`.

A blank cell means not measured, not zero. The last page before the appendix lists every module
whose output was absent and what was looked for (included even when empty). Genomes named in only one
input are listed with a count. The PDF uses matplotlib's `PdfPages`, so it needs no HTML toolchain or
network access.

### quantify_function

Sample-level gene and pathway abundance from the per-genome tables, in three files:

- relative abundance (compositional, comparable within a sample);
- copies per genome equivalent (comparable across samples; needs `--microbial-fraction`, otherwise
  skipped);
- length-normalised (comparable between gene families).

`--abundance-kind` - whether genome size is divided out: `coverage` (default, already per cell) or
`read_count` (divided by genome length). The assignable fraction is written beside the abundances.

### catalogue

One non-redundant gene catalogue across samples, for the bin-free path. MMseqs2 `linclust` at 95%
identity over 90% of the member (`--identity`, `--coverage`).

- `--quantify sample=r1,r2 ...` - maps reads to the catalogue with salmon (`--no-clip`, so poly-A
  runs in genes are kept) and writes `gene_abundance.tsv` in the long format `quantify_function` and
  `diff_abundance` read. Per-sample mapping rate goes to `gene_abundance_mapping.tsv`, with a warning
  below 50% (on the dev catalogue, 100.0% for reads from its own genes, 0.68% for unrelated gut
  reads).
- `--unit tpm` (default, comparable between genes) or `--unit reads`.

See the [2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#2-bin-free-gene-catalogue).

### amr

Resistance and virulence genes, their genomic context, and an evidence level. AMRFinderPlus runs on
nucleotides before any length filter (one measurement puts MEGAHIT's ARG-bearing contigs at 284 bp on
average). Virulence genes go to a separate table.

Each row has its distance to the nearer contig end and an evidence level on Cao et al. 2026's
four-level scale. This module reaches L1 (detection) and L2 (detection with genomic context); none
of its output supports "confers resistance". Context is `context_unknown` when the contig is too
short to establish it or no mobile-element caller has run.

- `--mobile-elements DIR` - a `metawrap2 viruses` or geNomad output directory, needed for L2.
  Proviruses give coordinates; plasmid and viral calls give whole contigs (a gene on a plasmid
  contig is at distance 0). Whole-contig calls cannot show absence.
- `--mobility-window` (default 10000) - distance for "no mobile element nearby".
- `--skip-hamronize` - skip `amr_hamronized.tsv` (hAMRonization's 36-field schema, with tool and
  database versions) and `amr_normalised.tsv` (argNorm's ARO term and drugs; on the dev set 4 calls
  map to 12 drugs). Only AMRFinderPlus rows go through argNorm, which cannot map RGI.
- `--rgi` - adds RGI against CARD (`metawrap2 install-env rgi`, `install-db card`). Rows record the
  tool, RGI rows their ARO accession, and the run prints the overlap joined on ARO. On five CAMI2
  bins AMRFinderPlus made 4 calls and RGI 8, 2 in common (`cepA` and `CepA-29` match only by ARO).
  Of RGI's 8, 1 was `Perfect` and the rest `Strict` homologs at 27.8-46.1% identity; filter on the
  `Cut_Off` tier in the `method` column. RGI missed `ant(6)-Ia` and `sat4`, found at 100% by
  AMRFinderPlus.

On one CAMI2 assembly (17 provirus regions, 344 plasmid contigs, 239 viral contigs), of 18 resistance
calls 5 were on a plasmid contig (one class-1 integron cassette), 1 had no element within 10 kb, and
12 were undetermined. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#3-antimicrobial-resistance).

### appraise

The per-cell share of each sample's community that the bins represent, and the lineages with no
genome. (`bin_refinement` reports the binned fraction of the assembly, which differs widely.) On the
test study the 26 refined bins of ERR011347 represent 11.6% of its bacterial community, with 323
unrepresented lineages, led by *Phocaeicola*, *Bacteroides* and Lachnospiraceae.

Input requirements (SingleM gives wrong results without an error otherwise):

- forward reads only (a read and its mate count as two cells);
- reads not trimmed below 60 bp, SingleM's window (trimmed 34 bp reads gave 10 marker sequences from
  12 million reads). Read lengths are checked first.

The metapackage is installed with `metawrap2 install-db singlem` (2.2 GB download, 9.7 GB on disk),
never automatically. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#6-unbinned-fraction).

### diff_abundance

Features that differ between groups of samples, by ANCOM-BC (compositional by default). Refuses
fewer than three labelled samples.

## Strains

### strains

Within-population variation, and whether two samples carry the same strain, with inStrain. Every
sample is mapped to the same genome set; `--bam` accepts existing alignments made that way.

Reports popANI between two samples per genome (same-strain threshold 99.999%, about one difference
per 100 kb) and nucleotide diversity within a sample. popANI is computed only over positions both
samples cover, so below 50% of the genome compared the verdict is `undetermined`; diversity is
filtered on `breadth_minCov`. On three gut metagenomes against 26 refined bins, 17 of 42 comparisons
reached popANI ≥ 99.999% but only 4 covered enough of the genome (the rest 0.1-46.8%).

- `--backend {instrain,tracs,both}` - TRACS estimates a transmission distance from a consensus SNP
  distance. Verdicts are reported side by side at each tool's threshold (popANI ≥ 99.999% / ≤ 10
  SNPs), never averaged. TRACS skips references under 25% covered. On two gut samples the tools
  agreed: inStrain called bin.1, bin.2 and bin.4 the same strain; TRACS called bin.1 (5 raw SNPs, 0
  filtered, 2.0 M sites) and bin.4 (2 → 0, 190 k sites) and skipped bin.2 for coverage. See the
  [2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#82-strains).

## Comparative genomics

### pangenome

One PPanGGOLiN pangenome per species, from `dereplicate`'s clusters (`-c`), with gene families
partitioned into persistent, shell and cloud.

- Groups under 15 genomes are flagged `partitioning_reliable=no` (on a three-sample cohort, 46
  clusters gave 7 pangenomes and 32 singletons).
- Groups under 50% persistent are flagged, with a pointer to the cluster's ANI and members'
  contamination. Of seven real groups, six were 68-96% persistent and one 21%.
- `--msa` - write each group's core-gene alignment, for `phylogeny`.

See the [2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#83-pangenomes).

### phylogeny

A tree from GTDB-Tk's 120-marker concatenation (species) or `pangenome --msa`'s core-gene alignment
(strains within a species). The output records which was used; universal markers cannot resolve
strains.

- Default: IQ-TREE with ModelFinder and ultrafast bootstrap (`--replicates`, default 1000).
- `--fast` - FastTree, for hundreds of genomes; `-lg` for protein, `-nt -gtr` for nucleotide,
  detected from the alignment.
- `--support-threshold` (default 95) - for UFBoot about as strict as 70 for a standard bootstrap.

Each tree is reported with the alignment's occupancy beside its support. Bootstrap support does not
detect missing markers, horizontal transfer or contamination. See the
[2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#84-phylogenies).

## Viruses

### viruses

Viral contigs (geNomad), their completeness, with host flanks trimmed off proviruses (CheckV), and
their populations. Quality tiers are CheckV's MIUViG determination; lower-bound completeness is
flagged (51 of 161 contigs on one real assembly).

- `-a` - run on the assembly. On one CAMI2 sample the bins held zero viral contigs and the assembly
  161. Use `-b` for proviruses in binned genomes. The plasmid calls feed `amr --mobile-elements`.
- `--crispr-hosts BIN_FOLDER` - host prediction from CCTyper spacers in the bins, exact match in both
  orientations, trusted arrays only (`--untrusted-arrays` for the rest). On five CAMI2 bins, 54
  spacers matched none of 161 viral contigs; contigs with no match can still have a host in the
  sample. See the [2.1.2 release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/2.1.2.md#4-viruses).
