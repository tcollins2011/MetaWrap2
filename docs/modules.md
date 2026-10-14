# Modules

MetaWrap2 keeps metaWRAP's module set and adds to it. Every module prints its own options with
`metawrap2 <module> -h`, which is the authoritative list; this page is the orientation.

## Pre-processing
- **read_qc** - read trimming and host (e.g. human) read removal. Trimming uses fastp by default.
  On the 10 CAMI2 GI samples it kept marginally more reads, bases and Q30 bases than Trim Galore for
  5.8% fewer CPU-hours and 29% less wall time, because it is single-pass and threads past the
  4-thread cutadapt ceiling. `--trimmer trim-galore` switches back, which helps if 1.3 GB of peak
  memory is too much: Trim Galore peaks at 0.1 GB. Re-measure either choice with
  `metawrap2 benchmark run cami2_gi trimmer`. See [Benchmarks](benchmarks.md#7-smaller-changes).
- **assembly** - metagenomic assembly with metaSPAdes or MEGAHIT, plus QC. MEGAHIT is the default.
  On three CAMI2 GI samples co-assembled, metaSPAdes returned a 38% longer N50 (12.9 kb vs 9.3 kb)
  but used 3.2x the CPU-hours and 9.7x the peak memory (140 GB vs 15 GB), and assembled 14% fewer
  bases in total. Use `--metaspades` when contiguity matters and the memory is available;
  `--megahit --metaspades` runs both and merges them. Re-measure with
  `metawrap2 benchmark run cami2_gi assembler`.

  The hybrid path feeds MEGAHIT the reads metaSPAdes could not place. The aligner streams straight
  into the filter, so no alignment is written to disk, and the scaffold index is deleted once the
  FASTQ exists.

  `--aligner` picks that aligner from the same registry `binning` and `blobology` use, but its
  default here is bwa. In binning the aligner only feeds a depth estimate, and minimap2 and bwa give
  the same bins (Spearman 0.984 on contig depth) for 1.37 CPU-hours against 7.77. Here the aligner is
  a filter: the reads that did not map are what MEGAHIT assembles.

  On the 273k-scaffold metaSPAdes assembly from the `assembler` arm and its 3-sample co-assembly
  reads, bwa cost 5.45 CPU-hours and recruited 2,867,680 reads; minimap2 cost 0.62 (8.8x less) and
  recruited 4,113,486, 44% more. The disagreement is one-sided: minimap2 calls 1,249,798 reads
  unmapped that bwa placed on a scaffold, against 3,992 the other way (Jaccard 0.695), because
  its `sr` preset is less sensitive on a fragmented real assembly. Those reads reach MEGAHIT labelled
  unassembled and are assembled a second time, and merging the two assemblies cannot tell a duplicate from a novel contig. `--aligner minimap2` is available if you accept that.

  The hybrid path cannot be smoke-tested on the `tiny` fixture. metaSPAdes' error-correction stage
  exits with *Failed to determine offset!* on the fixture's synthetic quality strings, before any
  assembly happens. The remapping stage itself runs on the fixture once it has scaffolds to map to.
- **kraken2** - taxonomy profiling of reads or contigs, with a kronagram.
- **assembly (long reads)** - `-L reads.fastq --read-type {ont,ont-hq,pacbio,pacbio-hifi}`
  assembles with metaFlye. Several `-L` files are one dataset to Flye, so a co-assembly needs only
  the flag. `--meta` is always passed, since Flye's default assumes one genome at even coverage.

  The platform setting lives in one place (`metawrap2.longreads`) and maps to Flye's `--nano-hq`,
  minimap2's `-x map-ont` and SemiBin2's `--sequencing-type=long_read`. The run checks the reads
  against the declared platform and writes N50, median, the fraction over 10 kb and the count under
  1 kb to `long_read_report.tsv`. This catches two mistakes: a path pointing at the short-read FASTQ
  (metaFlye assembles it badly instead of failing) and HiFi declared over nanopore-accuracy data.

  On 25× simulated R10 nanopore reads over three CAMI2 MAGs (2.50 Mbp) it produced 72 contigs,
  2.41 Mbp, N50 57.9 kb, and the scanner's read N50 of 12,326 bp matches Flye's own log. `binning
  --long-read-type ont-hq` then maps with `map-ont` for both index and alignment; the short-read
  preset gives low, uneven coverage that a binner reads as a different organism. See
  [2.1.2 release notes](releases/2.1.2.md#85-long-reads).

  Two things it does not do:
  - No short-read polishing. A polisher given reads from a close relative can rewrite one strain's
    consensus into another's, and that needs its own measurement.
  - No merging of a long-read and a short-read assembly, for the same duplicate-contig reason as above.
- **nonpareil** - diversity coverage: how much of each community the reads cover, and whether the
  samples can be compared. Pass the forward mate only (`-r`), since redundancy between a read and its
  own mate measures the library. Pass QC'd reads: Nonpareil assumes error well below 5%.

  It reports three numbers per sample: C, abundance-weighted average coverage; Nd, the Nonpareil
  diversity index, comparable across samples; and the projected effort to reach `--target-coverage`
  (default 95), as a multiple of what was sequenced. On the three ERR011347-9 gut samples in this
  repo's test study (2010-era reads, mean length 32-49 bp): coverage 84.4%, 74.6% and 89.4%,
  Nd 17.3-18.6, and 2.5x to 8x more sequencing needed to reach 95%.

  Model fitting is in Python (`metawrap2.nonpareil`), ported from Nonpareil's R package. On
  Nonpareil's own three example datasets it reproduces R's C, Nd, LRstar, model r and both gamma
  parameters to six significant figures.

  `--kernel auto` (the default) picks one kernel for the whole cohort. The two kernels disagree on
  short reads: on ERR011347 (mean read length 34 bp) the k-mer kernel reports 50.5% coverage and
  the alignment kernel 83.0%, projecting 47x against 4.4x more sequencing. At 34 bp and k=24 a read
  carries 11 k-mers. So `auto` uses the faster k-mer kernel unless a sample's mean read length is
  under twice k, in which case the whole cohort gets the alignment kernel. Mixing kernels would make
  the cohort's figures incomparable. The alignment kernel cost about 16 minutes per sample against 4
  for the k-mer kernel on these reads.

  Caveats:
  - The alignment kernel seeds from the clock unless given a seed; two runs over the same reads
    returned 83.0% and 84.3%. `-r 1337` is passed (`--seed`). The k-mer kernel has no seed, so its
    subsampling is not reproducible.
  - The k-mer kernel refuses any read shorter than k, and trimmed reads often contain them (3.4% of
    one real sample here). Those reads are filtered out first (`--no-filter` to skip) and the count
    is reported.
  - Under the k-mer kernel the effort axis is in k-mers: Nonpareil's `AL` is k-mers per read
    (L-k+1). Every row states its `effort_unit` and carries the sequenced base pairs separately.

  Details: [2.1.2 release notes](releases/2.1.2.md#5-diversity-coverage-can-these-samples-be-compared).

## Bins
- **binning** - initial binning with metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB and/or COMEBin.

  `metawrap2 run` uses metaBAT2 and SemiBin2 by default (as of 2026-10-09). End to end, that pair
  gives 53 clean MAGs against metaWRAP 1's 48, for 43.8 CPU-hours against 52.5. It needs the
  `metawrap2-learned` environment (~1.5 GB of pytorch, built by `quickstart`); `--preset fast`
  gives metaBAT2 + MaxBin2 instead.

  CONCOCT is opt-in (since 2026-10-02). On `cami2_gi` it took 192 of the 213 CPU-hours a
  three-binner run spent binning, for 74 clean MAGs against the classical pair's 64, while the
  learned binners reach 105 and more for a fraction of its cost. Add `--concoct` when the catalogue
  matters more than the run time.

  `metawrap2 binning` itself has no default and needs at least one binner flag or a preset.

  `--preset {very-fast,fast,default,slow}` picks a binner set by available compute.
  Measured on CAMI2 GI as clean high-quality MAGs (≥70/≤10 and non-chimeric by GUNC, counted per bin):

  | tier | binners | HQ | clean HQ | clean/HQ | CPU-h |
  |---|---|---|---|---|---|
  | `very-fast` | metaBAT2 | 37 | 35 | 94.6% | 4.8 |
  | `fast` | metaBAT2 + MaxBin2 | 71 | 59 | 83.1% | 7.1 |
  | `default` *(shipped)* | metaBAT2 + SemiBin2 | 131 | 114 | 87.0% | 86.8 (17 on a GPU) |
  | `slow` | metaBAT2 + SemiBin2 + COMEBin | 183 | 177 | 96.7% | 52.8 on a GPU |

  Each tier costs more than the one below and recovers more clean MAGs. The tiers are not nested,
  and a preset never overrides a binner you named. These counts come from binning CAMI2's gold
  assembly, which isolates the binner; on real per-sample assemblies the tiers sit closer together.
  See [Benchmarks](benchmarks.md#3-binning) for both designs side by side.

  `slow` runs three binners instead of COMEBin alone. On the gold assembly they tie (177 clean MAGs
  against 176); end to end on per-sample assemblies the trio returns 60 against COMEBin's 57.
  metaBAT2 is the third member because it reads the depth table SemiBin2 already built and costs
  ~1.5 CPU-hours; the MaxBin2 trio recovers one more genome for 24.5. In `default`, SemiBin2 alone on
  a GPU recovers 116 clean MAGs against the pair's 114 for a fifth of the CPU; metaBAT2 is kept for
  ~1.3 CPU-hours. `--comebin` alone still works. `very-fast`'s 4.8 CPU-hours is with minimap2; it was
  22.9 under bwa, nearly all of it alignment.

  `--min-identity PERCENT` (default 90) is the minimum alignment identity for a read to count
  toward a contig's depth. It is the one place MetaWrap2 departs from
  `jgi_summarize_bam_contig_depths`' parameters (97). On CAMI2 GI, 97 discards half the mapped bases
  against the true reference genomes. One set of alignments, six depth tables, metaBAT2 on each:

  | floor | mean depth | bins | AMBER F1 | purity | HQ | clean HQ |
  |---|---|---|---|---|---|---|
  | none | 15.71 | 90 | 0.6607 | 0.8374 | 48 | 45 |
  | ≥90 (default) | 15.70 | 88 | 0.6624 | 0.8375 | 48 | 47 |
  | ≥95 | 13.46 | 86 | 0.6574 | 0.8264 | 43 | 40 |
  | ≥97 *(jgi's)* | 8.35 | 74 | 0.6371 | 0.8175 | 37 | 35 |
  | ≥99 | 1.48 | 21 | 0.2452 | 0.7904 | — | — |

  At 90 the bins are more complete, more pure and less often chimeric than at 97: 47 clean
  high-quality MAGs against 35, a 34% gain. A looser floor did not cause the binner to merge strains
  here. Pass `--min-identity 97` to reproduce a jgi depth table (needed to compare against metaWRAP 1
  or metaBAT2's published validation), or `0` for no filter. This was measured on one dataset with
  one binner. See [Benchmarks](benchmarks.md#4-coverage-estimation).

  The binners read coverage. metaBAT2, MaxBin2 and SemiBin2 get it from a depth table
  MetaWrap2 computes (`metawrap2.depth`), CONCOCT through `concoct_coverage_table.py`. The BAM is a
  temporary file. `--aligner` picks the aligner: minimap2 by default, with `bwa` and `bwa-mem2`.
  Mapping one CAMI2 GI sample to a 210k-contig assembly, minimap2 and bwa gave bins of
  indistinguishable quality (AMBER F1 0.4519 against 0.4493, recall 0.3118 against 0.3103,
  precision 0.8204 against 0.8133, all within noise) for 1.37 CPU-hours against bwa's 7.77.
  minimap2 ranks contigs almost identically to bwa (Spearman 0.984) but differs on magnitude
  (Pearson 0.946), which is why the comparison was carried through to the bins. bwa-mem2 is
  bit-identical to bwa on the depth vector (Pearson 1.0000) at 1.5x its speed and needs 23.6 GB of
  RAM against bwa's 9.8. bowtie2 cost what bwa cost and agreed with it less closely than minimap2,
  so it is not offered. Caveat: this was one sample against CAMI2's pooled gold-standard assembly,
  which is far less fragmented than a real metaSPAdes assembly, and fragmentation is what separates
  aligners.

  Two cheaper routes were measured and rejected:
  - `strobealign --aemb` emits per-contig depth without a BAM, but not the per-base coverage that
    metaBAT2's intra-contig variance needs.
  - `fairy` sketches k-mers instead of aligning. It is 648x cheaper than bwa and scored highest on
    raw F1 per CPU-hour (17.68 against minimap2's 0.33). Like `--aemb` it is per-contig only, and
    CONCOCT needs per-fragment coverage from an indexed BAM, so neither can feed it.

  minimap2 keeps the BAM and works for all three classical binners. Single-sample fairy lost, but
  fairy's claim is multi-sample, so it is available as `coverage --coverage fairy`.

  `--depth <table>` bins from an existing depth table, aligning nothing and needing no reads.
  metaBAT2 and MaxBin2 only: CONCOCT builds its own coverage table from an indexed BAM and the three
  learned binners read the BAM directory, so `--depth` with one of those is refused by name. This is
  how multi-sample coverage costs N mappings instead of N per assembly (see `coverage` below). The
  table is checked against the assembly first. No overlap is an error. A partial overlap is reported
  from both sides: assembly contigs with no row are left unbinned, while rows naming contigs this
  assembly lacks are expected when slicing a cohort catalogue.

  metaBAT2 and MaxBin2 share one depth table, `work_files/contig_depth.txt`, built
  in its own checkpointed stage. MaxBin2's per-sample abundance files are split out of it, taking the
  mean columns and skipping the variances. `jgi_summarize_bam_contig_depths` 2.18 with
  `--noIntraDepthVariance` reports values that are not depths, and MaxBin2 still returns plausible
  bins on them. Re-binning `cami2_gi` on a correct table: precision 0.6485 → 0.8016, ARI 0.7079 →
  0.8747, F1 0.4228 → 0.4598. Any multi-sample MaxBin2 result produced before this fix is suspect,
  including the `metabat2_maxbin2` and `metabat2_maxbin2_concoct` benchmark arms. See
  [Benchmarks](benchmarks.md#8-defects-found-and-fixed).

  The table is computed by `metawrap2.depth` instead of `jgi_summarize_bam_contig_depths`, using the
  same algorithm: an identity filter, 75 bp trimmed off each contig end before averaging, and the
  intra-contig depth variance, over per-base depths from `samtools depth`. It agrees with jgi to
  Pearson 0.9992 on per-sample means and 0.995 on variances across all 210,632 contigs of the CAMI2
  GI assembly. Owning the code makes the table testable.

  CONCOCT's coverage path (`concoct_coverage_table.py` over a cut-up BED and indexed BAMs) was
  checked the same way and agrees with the jgi table: all 210,632 contigs present in both, sample
  columns in the same order, no contig covered in one and not the other. CONCOCT reads ~1.73x higher
  depths throughout (median 1.742 / 1.737 / 1.729 across the three samples), because jgi applies
  `--percentIdentity 97` and `samtools bedcov` counts every aligned base. This is close to a
  per-contig scale factor (within-contig spread of the ratio, median 1.22) with a tail where it is
  not (p99 3.4). `bench_work/check_concoct_coverage.py` is the check.

  Contigs with no usable base are dropped. A record with sequence but no A, C, G, T or U has a
  zero tetranucleotide vector and no read can align to it. VAMB refuses to start on one; metaBAT2,
  MaxBin2 and CONCOCT accept it and carry a contig they cannot place. CAMI2's gold-standard
  gastrointestinal assembly has 20 of them, 28.9 Mbp in total and the largest 6.6 Mbp (masked-out
  reference genomes).

  Once every requested binner has its depth table, `work_files/` drops the
  per-sample BAMs, their `.bai` indexes and the aligner's index of the assembly, keeping the depth
  tables, the assembly copy and CONCOCT's cut-up fasta. On one real sample that is 0.6 GB of 1.0 GB;
  the BAM alone is about 5 GB per sample in a cohort. minimap2's `.mmi` for a 154 MB assembly is
  479 MB, larger than the BAM beside it. `--keep-bams` keeps everything. Deleting them also clears
  the `align_reads` checkpoint, so a later `--resume` that adds a binner realigns.

  SemiBin2, VAMB and COMEBin share one environment (`metawrap2 install-env
  learned`), since all three are pytorch applications that solve together. The benchmark arms
  (`semibin2`, `metabat2_semibin2`, `vamb`, `metabat2_vamb`, `comebin`, `metabat2_comebin`) decide
  how they are used, on recovered clean MAGs with cost reported beside them. All three apply a
  200 kb minimum bin size, so none is credited for bins the others discard.

  `--semibin2` adds SemiBin2, the one binner that reads the BAMs itself: it computes its own
  abundance features alongside the k-mer embedding it trains, so it is given `--input-bam`. It lives
  in the learned-binner env (`metawrap2 install-env learned`), about 1.5 GB of pytorch. Two settings
  are pinned for reproducibility:
  - `--random-seed` is fixed. Without a seed, two runs over identical input disagree about which
    contigs form a bin.
  - The training device comes from `--device` (default cpu). SemiBin2's own default is `auto`.

  SemiBin2 needs contigs of at least 4000 bp to generate must-link pairs and refuses to run without
  them. It drops any bin under 200 kb, independent of `--min-len`.

  `--vamb` adds VAMB, a variational autoencoder over composition and abundance. It reads the BAM
  directory (`--bamdir work_files`). Three flags are set against VAMB's defaults, because each
  default fails without an error:
  - `--minfasta` has no default; without it VAMB writes a cluster TSV and no FASTA files. It is set
    to 200 kb, matching SemiBin2's floor and roughly metaBAT2's.
  - `-o ''` turns binsplitting off. VAMB's separator defaults to `C` (for `<sample>C<contig>`
    names); MetaWrap2 assembles one sample at a time, so a contig id containing a capital C would
    have its bin split on it. VAMB only warns when some identifier lacks the separator.
  - `-m` is floored at 250, VAMB's hard minimum, instead of its recommended 2000, so it does not
    discard contigs `-l` asked to bin.

  `--seed` is also fixed, though VAMB notes it does not guarantee determinism. VAMB refuses to run
  into an existing `--outdir`, so the module clears its working output first; the checkpoint has
  already decided the stage is being redone.

  `--comebin` adds COMEBin, the most expensive learned binner. It trains contrastively over six
  augmented views of the assembly, produces several Leiden clusterings at different resolutions, and
  picks one by scoring against single-copy marker genes. Its entry point is `run_comebin.sh`, which
  takes a BAM directory and globs `*.bam`. Its env also carries FragGeneScan, HMMER and pplacer.
  Three settings differ from its defaults:
  - `-d` comes from `--device` (default cpu). COMEBin's own default is `cuda`, which fails without a
    GPU.
  - `-s` is passed, because `run_comebin.sh` forwards `--seed` only when given; otherwise both the
    network and the Leiden clustering are unseeded.
  - `-l` is applied by filtering the input: `run_comebin.sh` never forwards `--contig_len`, so
    COMEBin bins at its built-in 1000 bp. When `--min-len` is higher, the assembly handed to COMEBin
    is pre-filtered so it sees the same contigs as the other binners. At or below 1000 the shared copy
    is used as-is.

  COMEBin needs no clearing between attempts: `run_comebin.sh` skips stages whose manifests exist,
  so the augmentation pass survives a `--resume`.

  None of the three learned binners can be smoke-tested end to end on the `tiny` fixture, each for a
  different reason. SemiBin2 refuses (no contig of 4000 bp). VAMB clustered 257 contigs into 160
  bins, all under the 200 kb floor, so it ends with *did not produce a single bin*. COMEBin completes
  augmentation and training (about four minutes on 8 CPU threads for 257 contigs) and then fails in
  clustering with *HMMsearch completed but found zero marker hits*, since the synthetic genomes carry
  no recognisable single-copy markers.

  `--device {cpu,cuda,cuda:N,auto}` sets where all three learned binners train. `auto` asks that
  binner's environment whether torch can see a device and falls back to the CPU otherwise, including
  when torch is missing or broken. The default is cpu for reproducibility: a GPU run and a CPU run of
  the same seeded network need not agree, so the same study re-analysed elsewhere could produce
  different MAGs.

  The CPU default is expensive for COMEBin. On CPU over three CAMI2 samples it reached
  epoch 20 of 200 after 305 CPU-hours, a projected ~1,485 for one arm, against 7.1 for metaBAT2 and
  MaxBin2 together. On the CPU, COMEBin does not finish. `metawrap2-learned` ships a CUDA-capable
  pytorch, so `--device cuda` needs no second install; the `semibin2_gpu` / `vamb_gpu` /
  `comebin_gpu` benchmark arms measure the gap.
- **coverage** - one cross-sample depth table for a whole cohort, in N alignments instead of N².

  Bins improve when every contig's depth is known in every sample. Mapping all N samples against
  each of N per-sample assemblies costs N² alignments: four hundred for twenty samples. This module
  produces the same twenty tables from twenty alignments. It concatenates the per-sample assemblies
  into one contig catalogue, maps each sample's reads once against it, builds one depth table, and
  splits that table back into one slice per assembly. Each sample is then binned on its own contigs
  with every sample's coverage:

  ```
  metawrap2 coverage -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
      -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
  metawrap2 binning -a s0/final_assembly.fasta --depth cohort_coverage/depth/s0.txt \
      -o s0/binning --metabat2 --maxbin2
  ```

  Naming makes the split safe. Contig names collide across assemblies, so each contig is prefixed
  with its assembly's label (`s0~k141_1`) before concatenation, and the split demultiplexes on that
  prefix. A catalogue where an original name contains the separator is refused; `--separator` takes
  any character the names do not use (`~` by default, which SPAdes, MEGAHIT and Flye never emit).
  Two assemblies sharing a label are also refused. That is why `-a label=path` exists: a cohort laid
  out as `sample_N/final_assembly.fasta` would otherwise label every sample `final_assembly`.

  The slices carry the original contig names, so they match the original assembly that `binning
  --depth` checks against. Every slice keeps all the sample columns. A label that matched no contig
  in the catalogue is an error.

  Only metaBAT2 and MaxBin2 can use the result; the other binners read the alignments themselves,
  and `binning` refuses them by name. The BAMs and the catalogue are deleted once the tables are
  written (`--keep-bams`, `--keep-catalogue` to keep them). The catalogue is the size of every
  assembly together, and minimap2's index of it is larger.

  It uses binning's conda environment, since it runs the same aligner, samtools and
  `jgi_summarize_bam_contig_depths`: `metawrap2 install-env coverage` builds (or finds)
  `metawrap2-bins`.

  `--coverage fairy` replaces the alignment. fairy sketches each sample's reads once, then estimates
  every contig's depth in every sample from k-mer containment. It writes the same
  `contigName contigLen totalAvgDepth <sample> <sample>-var ...` table
  `jgi_summarize_bam_contig_depths` writes, so everything downstream is the same:

  ```
  metawrap2 install-env bins
  metawrap2 coverage --coverage fairy -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
      -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
  ```

  It takes minutes where alignment takes hours, at a cost in resolution. fairy's default `-m 95`
  resolves nothing below roughly 95% ANI, so reads from a strain that diverges from its contig are
  credited to that contig or to nothing. A sketch cannot distinguish positions within a contig, so
  the `-var` columns are filled from containment, not per-base depth.

  Alignment is the default because it measured better. On one CAMI2 GI sample through metaBAT2,
  single-sample fairy scored AMBER F1 0.212 against minimap2's 0.452 for 0.012 CPU-hours against
  1.368, and raising its ANI threshold made agreement with the aligners worse. fairy's published
  claim is multi-sample (98.5% of MAGs recovered); the cohort comparison is in
  [Benchmarks](benchmarks.md#1-cohort-support).

  Two operational differences from the aligned route:
  - The sketches (`sketches/<sample>.paired.bcsp`) are kept. They are small and independent of the
    catalogue, so they serve a later run against a different one.
  - Interleaved input is sketched as single-end, since a sketch does not record which mate a k-mer
    came from.

  fairy is installed in `metawrap2-bins`.
- **bin_refinement** - consolidate multiple bin sets into a refined set. GUNC scores the final bins
  for chimerism and adds `gunc_css`/`gunc_pass` to their `.stats.tsv`; bins are annotated, never
  dropped. `--skip-gunc` opts out (and out of its 6.7 GB database).

  GUNC is on by default because it catches what CheckM cannot: on the `cami2_gi` `metabat2_concoct`
  arm, 4 of the 60 final bins are chimeric by GUNC despite passing CheckM at `-c 70 -x 10`. A bin
  built from halves of two genomes has a complete, uncontaminated-looking marker set. GUNC costs
  3.20 CPU-hours on a 60-bin set, ~32% on top of refinement. The same run shows consolidation working:
  the two input sets were 12.0% chimeric (metaBAT2, 9/75) and 36.0% (CONCOCT, 32/89), and the
  consolidated set is 6.7% (4/60). The final bins land in `metawrap_<c>_<x>_bins`, named for the
  thresholds they were cut at.

  `--consolidator` chooses how the sets become one set:

  - `metawrap` (default) hybridises every combination of the input sets with Binning_refiner and
    asks CheckM which version of each bin to keep. With *n* input sets that is CheckM on 2^n-1 sets;
    on a 3-binner arm, four Binning_refiner runs and four CheckM runs. On the `metabat2_maxbin2` arm
    of `cami2_gi`: 10.86 CPU-hours, 70 final bins at `-c 70 -x 10`.
  - `dastool` scores the input bins once on single-copy genes and aggregates in one pass, with no
    hybridisation and no per-set CheckM. It needs `-a/--assembly` and its own env (`metawrap2
    install-env binqc`), which brings R, ruby, blast, diamond, prodigal and pullseq. No database.
    `--dastool-score` is its quality floor (default 0.5); `work_files/dastool_bin_names.tsv` maps the
    renumbered bins back to DAS_Tool's names, which say which input set each came from.

  Both leave one set in `binsO`, and everything after that (CheckM re-scoring, the `-c`/`-x` cut,
  plots, final layout, GUNC) is shared code, so only the route to the one set differs. `metawrap` is
  the default as the incumbent. The `consolidator` comparison
  (`metawrap2 benchmark run cami2_gi consolidator`) decides it on NC MAGs recovered; a tie favours
  DAS_Tool, which is 7× cheaper. See [Benchmarks](benchmarks.md#5-refinement).

  Before either runs, every input bin set is filtered to bins between 50 kb and 20 Mb. This can
  remove a lot. A bin over 20 Mb holds several genomes and cannot become a MAG; one under 50 kb is a
  fragment. On the `cami2_gi` CONCOCT set the ceiling alone removed 7 of 99 bins carrying 41% of the
  set's sequence (the largest 131.8 MB), and 3 of metaBAT2's 78 carrying 39%. Refinement reports
  both ends with a count and a share, since much of any recall drop across refinement happens here
  and not in consolidation.

  With `metawrap`, two shortcuts reduce the cost of the 2^n-1 CheckM runs without changing a
  reported number:
  - A bin whose contigs match one already scored takes that score. Across the seven combinations of
    a three-binner run, 390 of 1,230 bins are repeats. Bins are keyed by a hash of their sorted contig
    ids, and the cache lives in `work_files` so `--resume` keeps it.
  - A bin outside 200 kb to 20 Mb is recorded as 0% complete / 100% contaminated without scoring.
    The bounds are calibrated on 18,490 scored bins, of which the 10,738 passing 70/10 run from
    0.99 Mb to 6.66 Mb. `--score-everything` turns this off. In refinement it saves little, because
    the 50 kb–20 Mb input filter has already run; it matters when the benchmark scores raw binner
    output, where 0.7% of bins hold 11.8% of all the sequence.

  An N50 floor was tested and dropped: the worst N50 among passing bins is 2,025 bp, so a safe floor
  skips nothing the size bounds do not, and a higher one discards real genomes.
- **reassemble_bins** - reassemble bins to improve N50/completion and reduce contamination.

  `--max-coverage` (default 50, `0` to disable) thins each bin's recruited reads before SPAdes.
  Recruitment is per bin, so an abundant organism hands SPAdes hundreds of times the coverage it can
  use, and SPAdes time grows with read count, not genome size. 50x is above where prokaryotic
  short-read assembly stops improving.

  The target is measured against that bin's own assembled length, so a 1.2 Mbp bin and a 6 Mbp bin
  get different read counts for the same depth. A bin already below the target is not touched,
  which is the common case. Reads are dropped in pairs with a fixed seed, so a re-run keeps the same
  reads and SPAdes never gets a library with one mate missing.
- **dereplicate** - collapse several samples' bins into one non-redundant genome catalogue.

  Binning twenty samples separately returns the same abundant organism twenty times. A MAG count
  over the pooled folders then counts bins, and a quantifier given near-identical genomes splits one
  organism's reads between them. `metawrap2 dereplicate -b sampleA/metawrap_70_10_bins -b
  sampleB/... -o catalogue/` takes any number of bin folders and returns one representative per
  cluster in `dereplicated_genomes/`, symlinked (or copied with `--copy`).

  Identity comes from skani (`triangle -E`), which sketches and scales to the thousands of pairs a
  cohort produces. Clustering and representative choice are in Python and unit-tested.
  - `--ani 95` (default) is the species threshold large MAG catalogues use, so counts are comparable
    with published ones. `--ani 99` separates strains.
  - `--min-aligned-fraction` (default 10%, permissive) is the share of a genome a pair must align:
    99% identity over 5% of two genomes is a shared prophage. Raising it gives smaller clusters and a
    larger catalogue.

  The representative is the best-scoring member, read from the `<bins>.stats.tsv` that
  `bin_refinement` leaves beside each folder (completeness minus five times contamination, size
  breaking ties). A scored genome beats an unscored one of any size, so a 6 Mb chimera cannot stand
  in for the clean 3 Mb copy of the same genome. With no CheckM table the module says so and falls
  back to size, where a chimera can win its cluster; run `bin_refinement`/`binning --run-checkm`
  first.

  Clustering is single linkage, so it chains: A-B and B-C at 96% merge all three even when A and C
  are not 96% identical. This never splits a real duplicate. `clusters.tsv` has a row per input bin,
  kept or not, with its cluster, source folder, representative, quality score, size and
  `ani_to_representative`; a member pulled in by chaining has no value there. One env, once:
  `metawrap2 install-env dereplicate`.
- **quant_bins** - estimate bin abundance across samples, with four quantifiers.

  sylph is the default. It does not align: it sketches the bins and the reads and profiles one
  against the other. On 76 bins from the CAMI2 GI gold-standard assembly across three samples,
  against a read-level truth, it ordered the bins best (Spearman 0.868 vs salmon's 0.827 and CoverM's
  0.803) for 0.19 CPU-hours against salmon's 1.08 and CoverM's 2.37, and peaked at 1.3 GB where
  CoverM peaked at 11.6 GB. It is weaker on magnitude: salmon correlates better linearly (Pearson
  0.959 vs 0.811), so use `--quantifier salmon` when the values matter more than the ranking.
  `--quantifier coverm` maps against the bins directly and was the slowest and most memory-hungry.
  Re-measure with `metawrap2 benchmark run cami2_gi quantifier`. See
  [Benchmarks](benchmarks.md#6-bin-abundance).

  `--quantifier depth --depth TABLE` is the cheapest and aligns nothing again. It aggregates
  binning's per-contig depth table into a length-weighted mean per bin (total bases over total
  length), with the per-contig median beside it: 0.013 CPU-hours, against 0.126 for sylph and 2.37
  for CoverM. Against the same truth, Pearson 0.953 and L1 similarity 0.799 put it second of four on
  magnitude, ahead of CoverM on both; Spearman 0.780 puts it last on ordering. Pearson is dominated
  by abundant bins and Spearman weights every bin alike, so noise on low-coverage bins hurts only the
  latter. Use it when the magnitudes matter and the table is already on disk.

  `--min-ani` (sylph only, default 99) moved the numbers most. At sylph's own default of 95%, more
  than half the abundance reported in one sample went to four bins whose genomes are absent from
  that sample. At 99.9 the floor starts discarding bins that are present. Bins from a fragmented real
  assembly match their own genome less tightly than these gold-standard ones, so lowering it towards
  95 is reasonable there; below 95 sylph's statistics are not calibrated and MetaWrap2 says so.

  All write `bin_abundance_table.tsv` and the same heatmap, but the cells mean different things.
  Do not compare tables from two quantifiers.
  - salmon: length-weighted median contig coverage, no ceiling.
  - CoverM: relative abundance as a percentage of the sample; a row sums to at most 100. A zero means
    "mapped, but under `--min-covered-fraction` of the genome was covered".
  - sylph: k-mer sequence abundance, also a percentage. A zero means "not detected above
    `--min-ani`".

  sylph and salmon ship in the module's own env; CoverM needs its own, once:
  `metawrap2 install-env bins`. All accept gzipped bins and ignore non-FASTA files in the bin folder,
  so a `-b` straight from `bin_refinement` works.
- **blobology** - GC-vs-abundance blobplots of contigs and bins.

  `--aligner` offers the same three choices as binning, from the shared registry in
  `metawrap2.aligners`, defaulting to minimap2. bowtie2 is no longer used: it cost what bwa cost
  (7.44 CPU-hours, 2.58 of it building the index) and agreed with bwa less closely than minimap2.

  The BAM is not sorted, since `gc_cov_annotate` sums aligned span per contig with a plain
  `samtools view` and record order does not matter. The alignment streams without a SAM on disk, and
  the cleanup step removes the BAMs and the aligner's index at the end of the run.
- **classify_bins** - conservative taxonomy assignment for bins.

  The default route is MEGABLAST against NCBI `nt` followed by taxator-tk, and it is deprecated:
  taxator-tk has had no upstream release since 2015. `--gtdbtk` is the maintained route and gives
  GTDB taxonomy, which the MAG literature uses. The default stays until a comparison with `nt` and
  GTDB present at once has been run, and because switching would break installs without GTDB's
  reference data. Neither route is cheap: `nt` is hundreds of GB, GTDB is ~110 GB and `classify_wf`
  needs a lot of RAM.
- **annotate_bins** - deprecated. Use `annotate`.

## Function

These modules are new. [plan_functional_annotation.md](plan_functional_annotation.md) explains the
design choices.

- **annotate** - calls genes once over a bin folder or a bare assembly, and annotates them. It
  replaces `annotate_bins`, whose per-bin Prokka output could not be joined across genomes and did
  not cover unbinned sequence.

  Genes come from Pyrodigal in-process. A bin trains its own coding model (`--meta` forces
  metagenomic mode); an assembly uses Prodigal's metagenomic models, since one model trained on many
  organisms is an average of them. Every gene carries `partial` and `distance_to_contig_end`: on
  CAMI2 bins, 8% of genes are partial in a bin of 92 contigs against 50-54% in bins of ~900.

  `--kofam` assigns KEGG orthologs with PyHMMER against KOfam's curated per-profile thresholds,
  applying each to the right score: 17,890 profiles calibrate on the full-sequence score and 8,697
  on the best domain. The 1,866 profiles with no threshold are dropped by default.

  `--pfam` assigns Pfam domains, covering genes KEGG cannot place. On one real CAMI2 bin, KEGG placed
  1,620 of 3,891 genes (41.6%), Pfam placed 3,185 (81.9%), and 1,578 of those have no KEGG ortholog.
  The run prints that increment. Domains are called at each family's own gathering threshold, since
  one e-value across 21,000 families suits long conserved families and not short ones.

  Pfam output is per domain, since a protein can carry several, or the same one twice: 54.7% of
  annotated genes on that bin have more than one. `pfam_architectures.tsv` gives the ordered list
  (`SusD-like_3 + SusD_RagB` is the commonest there, 41 genes). Every row has the versioned accession,
  which dates the call, and the bare family, which joins across releases. See
  [2.1.2 release notes](releases/2.1.2.md#86-protein-domains-for-the-genes-kegg-cannot-place).

  InterProScan is not wrapped: it is a Java wrapper around a dozen member databases, with a bioconda
  build four years old and about 50 GB of data. Pfam-A is 399 MB and searched with PyHMMER like the
  other HMMs. The other member databases are a gap, listed in the plan.

  `--cazy` adds CAZyme families with dbCAN over the same proteins. It answers which glycosidic bond
  an enzyme acts on, which a KO table answers poorly; the two annotations barely overlap for
  polysaccharide-degrading enzymes. dbCAN runs three independent methods and every row says how many
  agreed: on two CAMI2 bins, 415 CAZyme genes, 274 of them with at least two methods (dbCAN's own
  bar). The per-class summary keeps CBM apart from the catalytic classes, since a
  carbohydrate-binding module need not cut a polysaccharide. Genomes differ widely: bin_003 has 194
  glycoside hydrolases in 73 families, bin_006 has 24 in 18. Needs `metawrap2 install-db dbcan`
  (7.4 GB, a pinned release instead of dbCAN's moving `db_current` snapshot).
- **pathways** - KEGG module completeness per genome, plus the functional hierarchy and Krona plots.

  Every row has two completeness scores: `pathwise` (the best alternative route through the
  definition) and `stepwise` (the fraction of top-level steps satisfied). Five published tools use
  five incompatible definitions, and they disagree by more than the biological signal. Seven flags
  per row say how far to trust it: how many top-level steps the module has, whether it is a signature
  module, whether one of its reactions has no KO, and so on.

  The hierarchy is gene -> chemical process -> pathway, from KEGG BRITE: per genome with every gene
  weighted equally, and per sample with every gene weighted by how much of it is in the sample. The
  per-sample profile is not the average of the per-genome ones. Krona plots use `ktImportText`,
  already installed for kraken2 (`--skip-krona` to omit them).

  A sample profile needs a weight per gene, and the right weight depends on whether MAGs are the unit
  of analysis:

  | flag | weights each gene by | covers | cost |
  |---|---|---|---|
  | `-a/--abundance` | its genome's abundance | the binned community | free, from `quant_bins` |
  | `-C/--contig-depth` | its contig's depth | the whole assembly | free, from `binning`/`coverage` |
  | `-G/--gene-abundance` | its own abundance | the whole assembly | one extra mapping pass |

  Use `-a` when MAGs are the focus. Use `-C` when no bins are passed: the depth table from `binning`
  and `coverage` covers every contig, including the ~90% that never reach a MAG. `-G` answers per
  gene instead of per contig and costs a mapping pass (`catalogue --quantify` or
  `quantify_function`). Every hierarchy row carries a `weighting` column naming which of the four
  produced it. The choice changes the answer: on the dev cohort the gene-weighted profile puts 33% of
  the signal in "Not Included in Pathway or Brite" against ~14% for the MAG-weighted one, because
  unbinned sequence is enriched for genes BRITE has no place for. See
  [2.1.2 release notes](releases/2.1.2.md#1-function-from-gene-names-to-pathway-completeness).
- **quantify_function** - sample-level gene and pathway abundance from the per-genome tables.

  Three units, each in its own file:
  - relative abundance: compositional, comparable only within a sample;
  - copies per genome equivalent: comparable across samples (needs `--microbial-fraction`, otherwise
    skipped);
  - length-normalised: comparable between gene families.

  Whether genome size is divided out is set by `--abundance-kind` (default coverage, which is already
  per-cell; `read_count` divides by genome length). The assignable fraction is written beside the
  abundances, since a MAG-based profile describes only the binned part of a sample.
- **catalogue** - one non-redundant gene catalogue across any number of samples, for the bin-free
  path. Most assembled sequence never reaches a MAG, and this covers the rest. MMseqs2 `linclust` at
  95% identity over 90% of the member (`--identity`, `--coverage`).

  `--quantify sample=r1,r2 ...` maps reads against the finished catalogue with salmon and writes
  `gene_abundance.tsv` in the long format `quantify_function` and `diff_abundance` read, so the
  bin-free path reaches the same tables as the binned one. `--unit tpm` (default) is
  length-normalised and comparable between genes; `--unit reads` is not.

  Two details specific to a gene catalogue:
  - Salmon trims a trailing run of ten or more A's from every reference as a poly-A tail. In a gene
    catalogue that is sequence, and length is the denominator of every salmon unit, so `--no-clip`
    is passed.
  - The mapping rate per sample is written to `gene_abundance_mapping.tsv`, read from salmon's own
    record (a read mapping to several genes appears in several rows of the table). It measures how
    well the catalogue represents the sample, and the run warns below 50%. On the dev catalogue,
    reads cut from its own genes map at 100.0% and unrelated gut reads at 0.68%.

  See [2.1.2 release notes](releases/2.1.2.md#2-the-bin-free-path-now-quantified).
- **amr** - resistance and virulence genes, their genomic context, and the evidence level each call
  supports.

  AMRFinderPlus runs on nucleotides before any length filter, since assemblies break around these
  genes: one measurement puts MEGAHIT's ARG-bearing contigs at 284 bp on average. Every row carries
  its distance to the nearer contig end. "No mobile element within 10 kb" is unsupported for a gene
  2 kb from the contig edge, so the context verdict has a `context_unknown` state, used whenever no
  mobile-element caller has run. `--mobility-window` (default 10000) sets the distance.

  Every row carries an evidence level from Cao et al. 2026's four-level ladder. This module reaches
  L1 and L2: detection, and detection with genomic context. Nothing it produces supports "confers
  resistance".

  L2 needs a mobile-element caller, supplied with `--mobile-elements DIR`: a `metawrap2 viruses`
  output directory or a geNomad run directory. geNomad's provirus table gives elements with
  coordinates on the host contig; its plasmid and viral calls give whole contigs. A gene on a plasmid
  contig is reported at distance 0 with the kind named. Whole-contig calls cannot support the
  negative: a gene on a contig not called a plasmid stays `context_unknown`.

  On one CAMI2 assembly (17 provirus regions, 344 plasmid contigs, 239 viral contigs), of 18
  resistance calls, 5 sit on a plasmid contig, 1 has no element within 10 kb, and 12 are undetermined
  because their contig is too short. The five are one finding: `ant(2'')-Ia`, `qacEdelta1` and `sul1`
  on one 8.3 kb plasmid contig is a class-1 integron cassette, `qacEdelta1` being its
  3'-conserved-segment marker. Virulence goes to its own table.

  Two more tables come from the same run (`--skip-hamronize` turns both off):
  - `amr_hamronized.tsv`: the calls in hAMRonization's 36-field cross-tool schema, with the tool and
    database versions in the rows, since an ARG call depends on the database date.
  - `amr_normalised.tsv`: adds argNorm's ARO term and the drugs each gene confers resistance to. This
    lets two backends be combined (one gene has a different symbol in every database). On the dev
    set, 4 calls map to 12 drugs, because one aminoglycoside phosphotransferase covers nine of them.

  `--rgi` adds RGI against CARD as a second backend (`metawrap2 install-env rgi`, `install-db card`).
  Every row says which tool found it, and RGI's rows carry their ARO accession directly. The run
  prints what the two backends agree on, joined on ARO where known: on five CAMI2 bins, AMRFinderPlus
  makes 4 calls and RGI 8, with 2 genes in common. A symbol join would find one, since AMRFinderPlus'
  `cepA` and RGI's `CepA-29` are the same beta-lactamase at the same 99.67% identity.

  Read RGI's extra calls with its `Cut_Off` tier, in the `method` column. Of its 8 calls here, 1 is
  `Perfect` and the rest are `Strict` homologs at 27.8-46.1% identity, mostly vancomycin-cluster
  genes. RGI misses `ant(6)-Ia` and `sat4`, which AMRFinderPlus finds at 100%. Neither backend is a
  superset of the other.

  argNorm cannot map RGI rows (it accepts argsoap, deeparg, resfinder, amrfinderplus, groot and
  abricate) and raises on a hAMRonization table containing anything else. So only AMRFinderPlus rows
  go to argNorm; RGI's already carry their terms. See
  [2.1.2 release notes](releases/2.1.2.md#3-resistance-detection-context-and-a-shared-vocabulary).
- **appraise** - what share of each sample's community the bins represent, per cell, and which
  lineages have no genome. `bin_refinement` reports the binned fraction of the assembly; this
  reports the binned fraction of the community, and the two are far apart. On this repo's test
  study, the 26 refined bins of ERR011347 represent 11.6% of its bacterial community by coverage,
  with 323 unrepresented lineages, the most abundant being *Phocaeicola*, *Bacteroides* and
  Lachnospiraceae.

  SingleM finds markers in the reads and in the genomes and compares them. Two input requirements,
  both of which give a wrong answer instead of an error if missed:
  - Pass the forward reads only. SingleM counts a marker read and its mate as two observations of
    one cell.
  - Pass reads that were quality-filtered but not trimmed below 60 bp, SingleM's window length. On
    this sample the trimmed reads (mean 34 bp) yielded 10 marker sequences from 12 million reads,
    while the raw 75 bp reads give thousands. The module checks read lengths before running.

  The metapackage is `metawrap2 install-db singlem` (2.2 GB download, 9.7 GB on disk). Nothing
  fetches it implicitly. See [2.1.2 release notes](releases/2.1.2.md#6-what-the-bins-leave-out).
- **diff_abundance** - which features differ between groups of samples, by ANCOM-BC. It is
  compositional by default because a relative abundance table is sum-constrained, and a t-test on it
  tests the constraint. It refuses to test fewer than three labelled samples.

## Strains

- **strains** - within-population variation, and whether two samples carry the same strain. A MAG is
  a population consensus: two samples whose MAGs are 99.9% identical can carry different strains.
  inStrain reads the alignments, so it has every position's allele frequencies.

  The module maps every sample to the same genome set, since popANI between samples mapped to
  different references is not meaningful. `--bam` takes existing alignments, with the same
  requirement.

  It reports popANI between two samples for one genome (inStrain's "same strain" threshold is
  99.999%, about one difference per 100 kb) and nucleotide diversity within one sample.

  Both are gated on breadth. popANI is computed only over positions both samples covered deeply
  enough, so two samples sharing 2% of a genome at 100% identity give popANI 100%. Below 50% of the
  genome compared, the verdict is `undetermined`. On three real gut metagenomes against this repo's
  26 refined bins, 17 of 42 comparisons reach popANI ≥ 99.999%, and only 4 do so over enough of the
  genome; the other 13 sit at 0.1-46.8% compared. Diversity is gated on `breadth_minCov` instead of
  bare breadth: a genome can be 91% covered with a quarter of that deep enough to call a variant.

  `--backend tracs` adds a second caller, and `both` runs each. inStrain compares allele populations
  and reports an identity. TRACS estimates a transmission distance from a consensus SNP distance,
  accounting for within-host diversity and a clock; in Shaw & Yu's evaluation most tools
  overestimate SNP distances, TRACS least. Verdicts are reported side by side at each tool's
  published threshold (popANI ≥ 99.999% / ≤ 10 SNPs) and never averaged.

  On the same two gut samples they agree. inStrain called bin.1, bin.2 and bin.4 the same strain;
  TRACS called bin.1 (5 raw SNPs, 0 after filtering, over 2.0 M sites) and bin.4 (2 → 0 over 190 k
  sites) the same strain and skipped bin.2 for coverage. TRACS refuses any reference with under 25%
  of its genome covered, so a shallow sample produces no row. See
  [2.1.2 release notes](releases/2.1.2.md#82-strains-is-this-the-same-strain-in-two-samples).

## Comparative genomics

- **pangenome** - one pangenome per species in a cohort, built from `dereplicate`'s clusters (`-c`).
  A pangenome is defined over a species; over a whole catalogue it would describe a community.
  PPanGGOLiN partitions gene families into persistent, shell and cloud. The shell holds the
  strain-level gene content that a single representative genome does not show.

  Expect few groups: on a three-sample cohort, 46 clusters gave 7 pangenomes and 32 singletons.
  Below 15 genomes the partitioning is fitted to too little data, and those rows are flagged
  `partitioning_reliable=no`.

  The persistent share also checks the dereplication. Of seven real groups, six sat at 68-96%
  persistent and one at 21%: 7,361 families over a mean of 3,214 genes per genome, so those three
  "genomes" barely share gene content. Any group below 50% is flagged with a pointer to the cluster's
  ANI and its members' contamination.

  `--msa` writes each group's core-gene alignment, which `phylogeny` takes. See
  [2.1.2 release notes](releases/2.1.2.md#83-pangenomes-and-what-they-say-about-the-clustering).
- **phylogeny** - a tree from an alignment the pipeline already wrote: GTDB-Tk's 120-marker
  concatenation (species relative to each other) or `pangenome --msa`'s core-gene alignment (strains
  within one species). The output says which of the two it was given, because a tree from universal
  markers cannot resolve strains.

  IQ-TREE with ModelFinder and the ultrafast bootstrap by default (`--replicates`, default 1000);
  `--fast` uses FastTree for hundreds of genomes. Every tree is reported with the alignment's
  occupancy beside its support. Ultrafast bootstrap measures whether the alignment supports a branch;
  it cannot see markers missing from an incomplete MAG, horizontal transfer, or contamination placing
  another organism's marker in a genome's row. Read support at 95 (`--support-threshold`), which for
  UFBoot is about as strict as 70 for a standard bootstrap. See
  [2.1.2 release notes](releases/2.1.2.md#84-phylogenies-from-alignments-that-already-exist).

## Viruses

- **viruses** - viral contigs, their completeness, and the populations they form. geNomad finds
  them; CheckV estimates how much of a genome each one is and trims host flanks off proviruses.

  Point it at the assembly (`-a`). On one CAMI2 sample, screening the bins finds zero viral contigs
  and screening the assembly finds 161, since a phage contig is usually in no prokaryotic bin. Use
  bins (`-b`) for proviruses integrated in a binned genome. geNomad's plasmid calls from the same run
  feed `amr --mobile-elements`.

  Quality tiers are CheckV's own MIUViG determination. A completeness CheckV reports as a lower bound
  is flagged as one: on a real assembly 51 of 161 contigs are not reliable estimates, and treating
  them as measurements under-reports novel viruses.

  `--crispr-hosts BIN_FOLDER` predicts hosts from bins you already have. CCTyper reads CRISPR arrays
  off the prokaryotic bins; a spacer from bin X found in a viral contig is evidence that X is a host
  of it. Matching is exact, in both orientations, and uses only spacers from arrays CCTyper trusts
  (`--untrusted-arrays` to include the rest). Expect a low yield: on five CAMI2 bins, three carried
  arrays (I-B, I-C) with 54 spacers between them, and none matched the sample's 161 viral contigs.
  A link needs the host to carry a CRISPR system that has recorded this virus, so no link is not
  evidence of a missing host. See [2.1.2 release notes](releases/2.1.2.md#4-viruses).

Each module runs in a conda environment grouped by workflow; see
[Installation](installation.md). `metawrap2 install-env --all` builds every one of them, and
`metawrap2 check <module>` verifies a module's tools inside the environment the run will use.
