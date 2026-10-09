# Modules

MetaWrap2 keeps metaWRAP's module set and adds three. Every module prints its own options with
`metawrap2 <module> -h`, which is the authoritative list; this page is the orientation.

## Pre-processing
- **read_qc** - read trimming and host (e.g. human) read removal. Trimming is fastp by default,
  measured on the 10 CAMI2 GI samples: it kept marginally more reads, bases and Q30 bases than
  Trim Galore for 5.8% fewer CPU-hours and 29% less wall time, since it is single-pass and threads
  past the four cutadapt stops at. `--trimmer trim-galore` switches back - worth it if 1.3 GB of
  peak memory is too much, as Trim Galore peaks at 0.1 GB. Re-measure either choice with
  `metawrap2 benchmark run cami2_gi trimmer`.
- **assembly** - metagenomic assembly with metaSPAdes or MEGAHIT, plus QC. MEGAHIT by default,
  measured on three CAMI2 GI samples co-assembled: metaSPAdes returned a 38% longer N50 (12.9 kb
  vs 9.3 kb) but wanted 3.2x the CPU-hours and 9.7x the peak memory (140 GB vs 15 GB), and
  assembled 14% fewer bases in total. `--metaspades` is there for when contiguity is the point and
  the memory is available; `--megahit --metaspades` runs both and merges them. Re-measure with
  `metawrap2 benchmark run cami2_gi assembler`. The hybrid path, which feeds MEGAHIT the reads
  metaSPAdes could not place, used to write the aligner's full SAM to disk, filter it into a second
  SAM and delete both - roughly 25 GB of text at peak to produce one FASTQ. It now streams the
  aligner straight into the filter, so no alignment is written at all, and the scaffold index is
  deleted once the FASTQ exists.

  `--aligner` picks that aligner from the same registry `binning` and `blobology` use, but its
  default here is **bwa**, not the registry's minimap2 - and that is a measured result rather than
  caution. Everywhere else the aligner is an implementation detail of a coverage estimator, and
  minimap2 and bwa agree on contig depth closely enough (Spearman 0.984) that the bins come out the
  same for 1.37 CPU-hours against 7.77. Here the aligner is a *filter*: the only thing read out of
  the alignment is whether a read mapped, and the reads that did not are what MEGAHIT assembles.

  Measured on the 273k-scaffold metaSPAdes assembly from the `assembler` arm and the 3-sample
  co-assembly reads it was built from, bwa cost **5.45 CPU-hours** and recruited 2,867,680 reads;
  minimap2 cost **0.62** - 8.8x less - and recruited 4,113,486, 44% more. The disagreement is
  one-sided: minimap2 calls 1,249,798 reads unmapped that bwa placed on a scaffold, against 3,992 the
  other way (Jaccard 0.695). Its `sr` preset is less sensitive on a fragmented real assembly, so
  those reads reach MEGAHIT labelled unassembled and get assembled a second time - and combining the
  two assemblies cannot tell a duplicate from a novel contig. The cheaper aligner changes what the
  assembly contains, in the bad direction, which is the one thing a cost saving cannot buy.
  `--aligner minimap2` is there for anyone who wants the 8.8x knowing that. Worth noting this is the
  *opposite* conclusion to binning's, from the same two tools, because the property read out of the
  alignment is different.

  The hybrid path cannot be smoke-tested on the `tiny` fixture, for a reason that is nothing to do
  with the aligner: metaSPAdes' error-correction stage exits with *Failed to determine offset!* on the
  fixture's synthetic quality strings, before any assembly happens. The remapping stage itself runs on
  the fixture fine once it has scaffolds to map to.
- **kraken2** - taxonomy profiling of reads or contigs, with a kronagram.

## Bins
- **binning** - initial binning with metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB and/or COMEBin.

  `metawrap2 run` uses **metaBAT2 and SemiBin2** by default, as of 2026-10-09, because that pair is
  the first configuration measured here that beats metaWRAP 1 on both quality and cost at once: end
  to end, 53 clean MAGs against metaWRAP 1's 48 for 43.8 CPU-hours against 52.5. It costs an extra
  environment — `metawrap2-learned`, ~1.5 GB of pytorch, which `quickstart` builds — and
  `--preset fast` is how to decline that and get the previous metaBAT2 + MaxBin2 default.

  CONCOCT has been opt-in since 2026-10-02, dropped on cost: on `cami2_gi` it was 192 of the 213
  CPU-hours a three-binner run spent binning, and it bought 10 extra clean MAGs for that — 74
  against the classical pair's 64. What makes it hard to justify now is not its accuracy, which is
  good, but that the learned binners reach 105 and more for a fraction of its cost. Add `--concoct`
  back when the catalogue matters more than the clock.

  `metawrap2 binning` itself has no default and still wants at least one binner flag, or a preset.

  **`--preset {very-fast,fast,default,slow}`** picks a binner set by how much compute you have, since
  "which binners to support" and "which to run" are different questions. Measured on CAMI2 GI as clean
  high-quality MAGs — bins that are both ≥70/≤10 *and* non-chimeric by GUNC, counted per bin:

  | tier | binners | HQ | clean HQ | clean/HQ | CPU-h |
  |---|---|---|---|---|---|
  | `very-fast` | metaBAT2 | 37 | 35 | 94.6% | 4.8 |
  | `fast` | metaBAT2 + MaxBin2 | 71 | 59 | 83.1% | 7.1 |
  | `default` *(shipped)* | metaBAT2 + SemiBin2 | 131 | 114 | 87.0% | 86.8 (17 on a GPU) |
  | `slow` | COMEBin | 178 | **176** | **98.9%** | 39.6 on a GPU |

  Each tier costs more than the one below and recovers more clean MAGs; they are **not** nested. A preset
  never overrides a binner you named explicitly. These counts come from binning CAMI2's gold assembly,
  which isolates the binner; on a real per-sample assembly the tiers sit closer together — see
  [Benchmarks](benchmarks.md) for both designs side by side.

  `slow` was metaBAT2 + SemiBin2 + COMEBin until the ensemble was compared against a correctly counted
  COMEBin. The earlier figure for COMEBin alone — 158 clean MAGs — came from multiplying its 178 HQ bins by
  an 88.8% GUNC pass rate that GUNC had measured over all **376** bins COMEBin emits, 198 of which never
  pass 70/10 and were never MAG candidates. Counted per bin, 176 of its 178 HQ bins pass GUNC, so the
  ensemble's margin is one bin for 33% more CPU. The same thing happens a tier down: SemiBin2 alone on a GPU
  recovers 116 clean MAGs against the pair's 114 for a fifth of the CPU. Consolidating a strong learned
  binner with weaker classical ones is not reliably worth its cost on this data. `default` keeps the pair
  anyway, because metaBAT2 reads the same depth table SemiBin2 does and costs ~1.3 CPU-hours beside it.
  `--metabat2 --semibin2 --comebin` still works for anyone who wants the ensemble. `very-fast`'s figure was
  22.9 CPU-hours until that arm was re-run under today's default aligner; essentially all of it was bwa's
  alignment, and it is 4.8 with minimap2.

  **`--min-identity PERCENT`** (default **90**) is the minimum alignment identity for a read to count
  toward a contig's depth, and it is the one place MetaWrap2 departs from
  `jgi_summarize_bam_contig_depths`' parameters. jgi's 97 was inherited and never evaluated, and it turns
  out to be a 2× lever: on CAMI2 GI it discards half the mapped bases *against the true reference
  genomes*. Measured — one set of alignments, six depth tables from it, metaBAT2 on each:

  | floor | mean depth | bins | AMBER F1 | purity | HQ | clean HQ |
  |---|---|---|---|---|---|---|
  | none | 15.71 | 90 | 0.6607 | 0.8374 | 48 | 45 |
  | **≥90 (default)** | 15.70 | 88 | 0.6624 | 0.8375 | **48** | **47** |
  | ≥95 | 13.46 | 86 | 0.6574 | 0.8264 | 43 | 40 |
  | ≥97 *(jgi's)* | 8.35 | 74 | 0.6371 | 0.8175 | 37 | **35** |
  | ≥99 | 1.48 | 21 | 0.2452 | 0.7904 | — | — |

  97 loses 14 bins of 88 and gains nothing. The objection to loosening the floor is that a relative's
  reads then inflate a contig's depth, two strains come to look like one organism, and the binner merges
  them — so it should buy completeness at purity's expense. Both AMBER purity and GUNC say otherwise: at
  90 the bins are more complete, more pure, *and* less often chimeric, for **47 clean high-quality MAGs
  against 35** — a 34% gain. On strain-rich
  data the filter was not separating strains, it was removing signal the binner was using. **Pass
  `--min-identity 97` to reproduce a jgi depth table** — which is what comparing against metaWRAP 1, or
  against metaBAT2's own published validation, needs — or `0` for no filter. Measured on one dataset with
  one binner, so it is a flag and not a buried constant.

  **`--subsample FRACTION`** streams every Nth read pair to the aligner instead of all of them, with no
  subsampled copy written to disk. A regular stride rather than a random draw, so both mates keep the same
  records without coordinating a seed. It is a **trade, not a saving**: cost scales with the fraction and
  so does the MAG count — 0.5 saves 1.9× and costs 4 of 37 high-quality MAGs, 0.25 saves 3.9× and costs
  12, 0.1 saves 7.8× and costs 25. Purity holds while completeness falls, which says what happens:
  thinning does not corrupt a well-covered contig's depth, it pushes marginal contigs below the point
  where their cross-sample pattern is informative. Not a default, and the wrong lever for a whole study
  anyway — assembly is 86% of the CPU time.

  None of the binners reads an alignment. They read coverage - metaBAT2, MaxBin2 and SemiBin2 from a
  depth table MetaWrap2 computes itself (`metawrap2.depth`, which replaced
  `jgi_summarize_bam_contig_depths` and agrees with it to Pearson 0.9992 on means across all 210,632
  CAMI2 GI contigs), CONCOCT through `concoct_coverage_table.py` - so the BAM is a
  temporary file on the way to a depth vector, and the aligner is an implementation detail of a
  coverage estimator. `--aligner` picks it: **minimap2** by default, with `bwa` and `bwa-mem2`
  selectable. Mapping one CAMI2 GI sample back to a 210k-contig assembly, minimap2 and bwa produced
  bins of indistinguishable quality - AMBER F1 0.4519 against 0.4493, recall 0.3118 against 0.3103,
  precision 0.8204 against 0.8133, every one fractionally in minimap2's favour and every one inside
  noise - for **1.37 CPU-hours against bwa's 7.77**. Worth knowing that minimap2 is not the arm that
  agrees most closely with bwa's depth *numbers*: it ranks contigs almost identically (Spearman
  0.984) but differs on magnitude (Pearson 0.946). Agreement with bwa is not the goal - bwa is the
  incumbent, not the truth - which is why the comparison was carried through to the bins instead of
  stopping at the depth correlation. bwa-mem2 *is* bit-identical to bwa on the depth vector (Pearson
  1.0000) at 1.5x its speed, not the advertised 2-3x, and wants 23.6 GB of RAM against bwa's 9.8,
  which is a trade to opt into rather than to impose. bowtie2 was listed as a binning dependency and
  never actually invoked; measured alongside the rest it cost what bwa cost and agreed with it less
  closely than minimap2 did, so it is gone rather than promoted. One caveat on all of it: this was a
  single sample against CAMI2's pooled gold-standard assembly, which is far less fragmented than a
  real metaSPAdes assembly, and fragmentation is exactly what separates aligners.

  Two cheaper routes were measured and rejected, both worth knowing about. `strobealign --aemb`
  skips the BAM and emits per-contig depth directly, but per-contig depth is all it can emit: the
  intra-contig depth *variance* metaBAT2 uses as a binning signal needs per-base coverage, which is
  precisely what `--aemb` declines to compute. `fairy` does not align at all, sketching k-mers
  instead, and is 648x cheaper than bwa - and it scored highest of everything measured on raw F1 per
  CPU-hour (17.68 against minimap2's 0.33), which is exactly what a decision rule of accuracy-per-cost
  says to choose. It is still the wrong answer, and that is the useful part: both routes are
  per-contig only, and CONCOCT cuts contigs into 10 kb fragments and needs per-fragment coverage out
  of an indexed BAM, so neither can feed it at all. A ratio needs a fidelity floor underneath it or
  it will eventually recommend measuring nothing. minimap2 keeps the BAM and is therefore a drop-in
  for all three binners.

  `--aemb` is closed by that finding - a route that cannot produce metaBAT2's variance column cannot
  be fixed by wiring it differently. fairy is not: what lost above was *single-sample* fairy, and
  fairy's claim is a multi-sample one. It is wired as `coverage --coverage fairy` for that comparison
  and nothing else, and the aligned route stays the default meanwhile.

  `--depth <table>` bins from a depth table that already exists, aligning nothing and needing no
  reads. metaBAT2 and MaxBin2 only, for the reason above: CONCOCT builds its own coverage table from
  an indexed BAM and the three learned binners read the BAM directory, so asking for `--depth` with
  one of those is refused by name rather than aligning the reads anyway and ignoring the table. This
  is the piece multi-sample coverage needs to cost N mappings instead of N per assembly — map each
  sample once against the concatenated per-sample assemblies, then hand each sample's binning its own
  slice of the one table. The table is checked against the assembly first, because the input that
  does damage is not a malformed table but a well-formed one built over different contigs: no overlap
  at all is an error, and a partial overlap is reported from both sides, since assembly contigs with
  no row are left unbinned while rows naming contigs this assembly lacks are the expected result of
  slicing a cohort catalogue.

  metaBAT2 and MaxBin2 share **one** depth table, `work_files/contig_depth.txt`, built in its own
  checkpointed stage before either of them runs. They used to build one each, and MaxBin2's pass asked
  for `--noIntraDepthVariance` - which reads as a formatting flag, since MaxBin2 wants one depth per
  contig per sample and has no use for the variance columns. In `jgi_summarize_bam_contig_depths` 2.18
  that code path does not drop columns; it reports numbers that are not depths. Running the same BAMs
  through the same tool with and without the flag, on a 60-contig probe it zeroed 59 of the 60 contigs
  and gave the 60th a depth of 7867 where the real per-sample depths were 0.267 and 0.453 - and 7867
  came back as 7874 on an identical rerun, so the value is uninitialized memory. On the real
  210,632-contig CAMI2 GI assembly with three samples it reported no coverage for **57% of the contigs
  that have coverage** in two of the three samples, gave ~25,000 contigs with no reads a coverage
  anyway, and overstated others by up to 1.5e9x. Nothing failed: a table of wrong depths is a
  well-formed abundance file, and MaxBin2 bins on composition too, so it kept returning plausible
  bins. MaxBin2's per-sample abundance files are now split out of the shared table, taking the mean
  columns and skipping the variances next to them - which also removes a second full pass over every
  BAM, the expensive part of this stage. **Any multi-sample MaxBin2 result produced before this is
  suspect**, including the `metabat2_maxbin2` and `metabat2_maxbin2_concoct` benchmark arms.

  Re-binning `cami2_gi` on a real table and scoring both sets together puts a number on what it cost:
  precision **0.6485 → 0.8016**, ARI **0.7079 → 0.8747**, F1 0.4228 → 0.4598. The fixed run returns 129
  bins over 349.8 Mbp where the buggy one returned 188 over 384.0, and recall rises anyway (0.3137 →
  0.3223) - the bins and the sequence that went were not genomes, they were contigs the invented depths
  had placed in the wrong bin.

  CONCOCT's coverage comes from a different code path - `concoct_coverage_table.py` over a cut-up BED
  and indexed BAMs - so it got the same question asked of it, and it answers cleanly. Aggregating the
  fragment coverages back to contigs by length-weighted mean and comparing against the jgi table built
  from the same three BAMs: all 210,632 contigs present in both, sample columns in the same order, no
  contig with coverage in the jgi table and none in CONCOCT's. CONCOCT reads **~1.73x higher** depths
  throughout (median 1.742 / 1.737 / 1.729 across the three samples), which is `jgi`'s
  `--percentIdentity 97` - jgi counts only reads aligning end to end at 97% identity, and
  `samtools bedcov` underneath CONCOCT counts every aligned base. On CAMI2 the reads carry strain
  variation the gold-standard assembly does not, so that filter drops a large share of them. It is
  close to a per-contig scale factor rather than a distortion of the abundance profile the binners
  actually use - for contigs covered in all three samples the within-contig spread of that ratio has a
  median of 1.22 - with a tail where it is not (p99 3.4). Both numbers are what their tools document;
  the point of recording them is that the two abundance paths agree, which is what the MaxBin2 table
  did not. `bench_work/check_concoct_coverage.py` is the check.

  `--semibin2` adds **SemiBin2**, which is the one binner here that does read the BAMs - it computes
  its own abundance features alongside the k-mer embedding it trains, so it is handed
  `--input-bam` rather than a depth table. It lives in its own env (`metawrap2 install-env
  binning-semibin2`) because it brings pytorch, roughly 2-3 GB that the three default binners have no
  use for. Two of its settings are pinned rather than left at SemiBin2's defaults, and both are about
  getting the same answer twice:

  - `--random-seed` is fixed. SemiBin2 trains a network, so without a seed two runs over identical
    input disagree about which contigs form a bin - which a benchmark would measure as part of the
    difference between binners, and which means a user cannot reproduce their own MAGs.
  - Its training device comes from `--device`, which covers all three learned binners and defaults to
    **cpu**. SemiBin2's own default is `auto`, meaning "a GPU if torch can see one"; that was found by
    running it on a machine with several GPUs and noticing *Running with GPU* in the log. See
    **`--device`** below for the whole argument - it is about reproducibility, and what it costs is
    now measured rather than assumed.

  Two things about it are worth knowing before reaching for it. It needs contigs of at least 4000 bp
  to generate must-link pairs and refuses to run without them, so it cannot be smoke-tested on the
  `tiny` fixture. And it drops any bin under 200 kb, which is its own floor rather than `--min-len`.

  `--vamb` adds **VAMB**, which trains a variational autoencoder over composition and abundance
  together. It reads the BAM *directory* - `--bamdir work_files` - rather than individual BAMs or a
  depth table. It shares SemiBin2's env (`metawrap2 install-env learned`): all three
  learned binners are pytorch applications and they solve together, so three environments were three
  copies of the same ~1.5 GB dependency graph. Three of its flags are set against VAMB's
  defaults, and each of them is a silent failure if left alone:

  - `--minfasta` **has no default**. Without it VAMB writes a cluster TSV and no FASTA files at all,
    which from outside the tool is a binner that succeeded and found nothing. It is set to 200 kb,
    matching SemiBin2's floor and roughly metaBAT2's, so none of the three gets credit for bins the
    others discarded.
  - `-o ''` turns binsplitting off. VAMB's separator defaults to `C`, for the multi-sample
    convention where a contig is named `<sample>C<contig>`; MetaWrap2 assembles one sample at a
    time, so there is no prefix to split off, and a contig id that happens to contain a capital C
    gets its bin split on it. VAMB only warns when *some* identifier lacks the separator, so this
    cannot be left to its leniency.
  - `-m` is floored at **250**, VAMB's own hard minimum, rather than at its recommended 2000. `-l`
    is the user's call, and 2000 would quietly discard contigs they asked to bin.

  `--seed` is fixed as well, with VAMB's own caveat that it does not guarantee determinism - which is
  a reason to pass one rather than not to. And VAMB refuses to run into an `--outdir` that already
  exists, so the module clears its working output first; that is safe because the checkpoint has
  already decided the stage is being redone, so nothing in there is a result.

  None of the three learned binners can be smoke-tested end to end on the `tiny` fixture, and each
  fails it for a different reason, which is worth knowing before concluding one of them is broken.
  SemiBin2 refuses outright, needing a contig of at least 4000 bp. VAMB gets furthest - it trained and
  clustered 257 contigs into 160 bins - then writes none of them, every one being under the 200 kb
  floor, so the run ends with *did not produce a single bin* rather than an empty bin set. COMEBin
  completes its augmentation and contrastive training (about four minutes on 8 CPU threads for 257
  contigs, which is its own signal about cost) and then fails in clustering with *HMMsearch completed
  but found zero marker hits*: the fixture's synthetic genomes carry no recognisable single-copy
  markers for it to seed on.

  `--comebin` adds **COMEBin**, the third and most expensive of the learned binners: it trains
  contrastively over six augmented views of the assembly, then produces several Leiden clusterings at
  different resolutions and picks between them by scoring each against single-copy marker genes. That
  last step is close to what `bin_refinement` does across binners, done internally across its own
  candidates. Its entry point is a shell script, `run_comebin.sh`, and like VAMB it takes a BAM
  directory and globs `*.bam` out of it. Its env (`metawrap2 install-env learned`) also
  carries FragGeneScan, HMMER and pplacer for the marker-gene step. Three things are set against its
  defaults:

  - `-d` comes from `--device`, which defaults to **cpu**; COMEBin's own default is plain **`cuda`**.
    This is a stronger version of SemiBin2's problem: `auto` at least falls back, whereas `cuda` asks
    torch for a GPU and fails without one. **COMEBin is also the binner for which the CPU default is
    least affordable** - see `--device`.
  - `-s` is passed, because `run_comebin.sh` forwards `--seed` **only when the flag is given**. So
    COMEBin's out-of-the-box behaviour is an unseeded network *and* an unseeded Leiden clustering.
  - `-l` is honoured by filtering the input, not by a flag: `run_comebin.sh` never forwards
    `--contig_len`, so COMEBin bins at its built-in 1000 bp whatever you ask for. When `--min-len`
    is higher, the assembly handed to COMEBin is pre-filtered to it, so its bins are built on the
    same contigs every other binner saw. At or below 1000 nothing is written and the shared copy is
    used as-is.

  COMEBin needs no clearing between attempts, unlike VAMB: `run_comebin.sh` checks for its own
  intermediate manifests and skips the stages that already completed, so the expensive augmentation
  pass survives a `--resume`.

  None of the three learned binners is a default, and none will become one on its paper's claims. The
  benchmark arms (`semibin2`, `metabat2_semibin2`, `vamb`, `metabat2_vamb`, `comebin`,
  `metabat2_comebin`) are what decide it, on recovered clean MAGs with cost reported beside them, like
  everything else. All three
  happen to share a 200 kb minimum bin size, so none of them is credited for bins the others threw
  away.

  **`--device {cpu,cuda,cuda:N,auto}`** - where the learned binners train. One switch for all three,
  replacing `--semibin2-engine` and `--comebin-device` and reaching VAMB, which had no way to ask at
  all. `auto` asks that binner's own environment whether torch can see a device and falls back to the
  CPU if it cannot; a missing or broken torch resolves to the CPU rather than raising, because that is
  the case `auto` exists to survive.

  The default is **cpu**, and that is a reproducibility decision rather than a claim about hardware: a
  GPU run and a CPU run of the same seeded network need not agree, so the same study re-analysed
  elsewhere would produce different MAGs for no reason its author could see. A default that silently
  moves to an accelerator is also a default chosen on hardware most people running a metagenomics
  pipeline do not have.

  What that default costs is measured, not assumed, and for COMEBin it is severe: on CPU
  over three CAMI2 samples it reached **epoch 20 of 200 after 305 CPU-hours** - a projected ~1,485 for
  one arm, against 7.1 for metaBAT2 and MaxBin2 together. On the CPU, COMEBin is not slow; it does not
  finish. `metawrap2-learned` ships a CUDA-capable pytorch, so `--device cuda` needs no second install,
  and the `semibin2_gpu` / `vamb_gpu` / `comebin_gpu` benchmark arms exist to put a number on the gap.

  **The depth table is built here**, by `metawrap2.depth`, rather than by shelling out to
  `jgi_summarize_bam_contig_depths`. Same algorithm - the same 97% alignment-identity filter, the same
  75 bp trimmed off each contig end before averaging, the same intra-contig depth variance - computed
  over per-base depths from `samtools depth`. It agrees with the tool it replaces to Pearson 0.999 on
  the per-sample means and 0.995 on the variances across all 210,632 contigs of the CAMI2 GI assembly.

  The reason is not that the tool was wrong. It is that this table is the most load-bearing number in
  the pipeline and was the one thing the pipeline could not look inside: when MaxBin2 was found to be
  binning on depths jgi had invented, the only way to see it was to run a C++ binary twice and notice
  that the answer changed. One table, one pass, and code with tests.

  **Contigs with no usable base are left out.** A record with sequence but not one A, C, G, T or U in
  it has a zero tetranucleotide vector and no read can align to it, so there is no assignment of it
  that could be right. VAMB refuses to start on one; metaBAT2, MaxBin2 and CONCOCT take it and
  silently carry a contig they cannot place. CAMI2's gold-standard gastrointestinal assembly ships 20
  of them, 28.9 Mbp in total and the largest 6.6 Mbp - masked-out reference genomes, not assembly
  artefacts - which is how this was found.

  The BAM is scratch, and the module now treats it that way: once every requested binner has taken
  its depth table, `work_files/` drops the per-sample BAMs, their `.bai` indexes, and the aligner's
  index of the assembly, keeping the depth tables, the assembly copy and CONCOCT's cut-up fasta.
  On one real sample that is 0.6 GB of 1.0 GB, and the saving scales with the cohort - about 5 GB
  per sample for the BAM alone, against a few megabytes of TSV that is all any binner ever wanted
  out of it. The aligner index is worth deleting too and is easy to overlook: minimap2's `.mmi` for
  a 154 MB assembly is 479 MB, larger than the BAM beside it. `--keep-bams` keeps the lot for
  anyone who wants to inspect an alignment afterwards. Deleting them also clears the `align_reads`
  checkpoint, so a later `--resume` that adds a binner realigns rather than resuming into a stage
  whose inputs are gone.
- **coverage** - one cross-sample depth table for a whole cohort, in N alignments instead of N².

  A cohort gets better bins when every contig's depth is known in every sample: a contig that rises
  and falls with its neighbours across twenty samples is evidence no single sample can give. The
  obvious way to collect that costs N² alignments - map all N samples' reads against each of the N
  per-sample assemblies - which for twenty samples is four hundred alignments to produce twenty depth
  tables.

  This module gets the same twenty tables for twenty alignments, by the multi-split route:
  concatenate the per-sample assemblies into one contig catalogue, map each sample's reads once
  against the catalogue, build one depth table over it, then split that table back into one slice per
  assembly. Each sample's binning then runs on its own contigs with every sample's coverage:

  ```
  metawrap2 coverage -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
      -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
  metawrap2 binning -a s0/final_assembly.fasta --depth cohort_coverage/depth/s0.txt \
      -o s0/binning --metabat2 --maxbin2
  ```

  What makes the split safe is the naming, which is the only part of this with a way to go quietly
  wrong. Contig names collide across assemblies - every assembler numbers from one - so each contig
  is prefixed with its assembly's label (`s0~k141_1`) before concatenation, and the prefix is what
  the split demultiplexes on. That is unambiguous only if no original name contains the separator, so
  a catalogue where one does is refused rather than guessed at later; `--separator` takes any
  character the names do not use (`~` by default, which SPAdes, MEGAHIT and Flye never emit). Two
  assemblies sharing a label are refused for the same reason, which is why `-a label=path` exists at
  all: a cohort laid out as `sample_N/final_assembly.fasta` would otherwise label every sample
  `final_assembly`.

  The slices carry the original contig names, not the prefixed ones, because they are handed to a
  binning run over the original assembly. A table of prefixed names against an unprefixed assembly is
  exactly the well-formed-but-wrong input that `binning --depth` checks for, so leaving them in would
  turn a correct run into a hard stop one command later. Every slice keeps all of the sample columns
  - that is the point of the exercise - and a label that matched no contig in the catalogue is an
  error, not an empty file to hand a binner.

  Only metaBAT2 and MaxBin2 can use the result: the other binners read the alignments themselves, so
  there is no table to hand them and `binning` says so by name rather than silently aligning anyway.
  The BAMs and the catalogue are swept once the tables are written (`--keep-bams`,
  `--keep-catalogue`), which matters more here than in `binning` - the catalogue is the size of every
  assembly together and minimap2's index of it is larger again.

  It shares binning's conda environment rather than having one of its own, because it runs the same
  aligner, samtools and `jgi_summarize_bam_contig_depths`: `metawrap2 install-env coverage` builds
  (or finds) `metawrap2-bins`.

  `--coverage fairy` replaces the alignment entirely. fairy sketches each sample's reads once, then
  estimates every contig's depth in every sample from k-mer containment, writing the same
  `contigName contigLen totalAvgDepth <sample> <sample>-var ...` table
  `jgi_summarize_bam_contig_depths` writes - so the split, the slices and everything downstream are
  byte-for-byte the same shape, and the flag is the only difference between the two routes:

  ```
  metawrap2 install-env bins
  metawrap2 coverage --coverage fairy -a s0=s0/final_assembly.fasta -a s1=s1/final_assembly.fasta \
      -o cohort_coverage/ s0_1.fastq s0_2.fastq s1_1.fastq s1_2.fastq
  ```

  It is minutes where the aligned route is hours, and the cost is resolution: fairy's default
  `-m 95` means nothing below roughly 95% ANI is resolved, so reads from a strain that diverges from
  the contig it assembled into are credited to that contig or to nothing, and a sketch cannot
  distinguish positions within a contig. metaBAT2's intra-contig depth variance therefore is not a
  real measurement on this route - the `-var` columns are filled, but from containment rather than
  per-base depth.

  The aligned route is the default because it measured better, not merely because it is older. On one
  CAMI2 GI sample through the same metaBAT2, single-sample fairy scored **AMBER F1 0.212 against
  minimap2's 0.452** - half the bins - for 0.012 CPU-hours against 1.368, and raising its ANI
  threshold made agreement with the aligners *worse* rather than better. It is 100× cheaper and
  recovers half the bins, and a ranking that divided one by the other would have made it the winner —
  which is the clearest argument in this repo against ranking on a ratio at all. What that run could
  not test is fairy's own claim, which is multi-sample (98.5% of MAGs
  recovered) and needs exactly the cross-sample table this module builds; that is benchmark B5's
  coverage arm, and the default does not move before it runs.

  Two operational differences. The sketches (`sketches/<sample>.paired.bcsp`) are kept rather than
  swept, unlike the BAMs: they are small, and they do not depend on the catalogue, so the same
  sketches serve a later run against a different one. And an interleaved input is sketched as
  single-end, because a sketch records which k-mers a sample contains and not which mate they came
  from, so splitting the file would cost a pass over the reads to produce the same sketch. fairy gets
  its own environment (`metawrap2-bins`) rather than a line in binning's, so a user on the
  default route never waits for it to solve.
- **bin_refinement** - consolidate multiple bin sets into a superior refined set. GUNC also
  scores the final bins for chimerism and adds `gunc_css`/`gunc_pass` to their `.stats.tsv`;
  bins are annotated, never dropped. `--skip-gunc` opts out (and out of its 6.7 GB database).
  It is on by default because CheckM cannot see what it sees: on the `cami2_gi` `metabat2_concoct`
  arm, **4 of the 60 final bins are chimeric by GUNC despite passing CheckM at `-c 70 -x 10`** - a bin
  built from halves of two genomes has a complete, uncontaminated-looking marker set. It costs 3.20
  CPU-hours on a 60-bin set, ~32% on top of the refinement it checks. The same run is also the evidence
  that consolidation works: the two input sets were 12.0% chimeric (metaBAT2, 9/75) and 36.0% (CONCOCT,
  32/89), and the consolidated set is 6.7% (4/60).
  The final bins land in `metawrap_<c>_<x>_bins`, named for the thresholds they were cut at, so a
  result carries the thresholds that produced it.

  `--consolidator` chooses *how* the sets become one set, and the two ways are a cost/accuracy
  trade rather than two implementations of one idea:

  - **`metawrap`** (default) hybridises every combination of the input sets with Binning_refiner
    and asks CheckM which version of each bin to keep. With *n* input sets that is CheckM on
    2^n-1 sets; on a 3-binner arm, four Binning_refiner runs and four CheckM runs. Measured on
    the `metabat2_maxbin2` arm of `cami2_gi`: **10.86 CPU-hours**, 70 final bins at `-c 70 -x 10`.
  - **`dastool`** scores the input bins once on single-copy genes and aggregates in one pass -
    no hybridisation and no per-set CheckM at all. Needs `-a/--assembly`, because it reads the
    contigs itself, and its own env (`metawrap2 install-env binqc`), which
    brings R, ruby, blast, diamond, prodigal and pullseq. No database. `--dastool-score` is its
    quality floor on the bins it keeps (default 0.5); `work_files/dastool_bin_names.tsv` maps the
    renumbered bins back to DAS_Tool's own names, which say which input set each came from.

  Before either runs, every input bin set is filtered to bins between **50 kb and 20 Mb**, and it is
  worth knowing this can be a large filter rather than a formality. No prokaryotic genome is 20 Mb, so
  a bin over the ceiling is a catch-all holding several, and it could not have become a MAG under any
  completeness/contamination threshold; under the floor is a fragment. On the `cami2_gi` CONCOCT set
  the ceiling alone removed 7 of 99 bins carrying **41% of the set's sequence** (the largest 131.8 MB),
  and 3 of metaBAT2's 78 carrying 39%. Refinement reports both ends with a count and a share, because a
  binner that puts most of its sequence in catch-all bins looks like a high-recall binner right up to
  this point - and because a recall drop across refinement is otherwise easy to attribute to the
  consolidation step, which is not where most of it happens.

  With `metawrap`, two things keep the 2^n-1 CheckM runs from costing what they look like they cost, and
  neither changes a reported number. **A bin whose contigs match one already scored takes that score** —
  across the seven combinations of a three-binner run, 390 of 1,230 bins are the same bin again, some
  appearing in all seven. A bin is keyed by a hash of its sorted contig ids, and the cache lives in
  `work_files` so `--resume` keeps it. And **a bin outside 200 kb to 20 Mb is recorded as 0% complete /
  100% contaminated without being scored**; the bounds are calibrated against 18,490 bins that have real
  scores, of which the 10,738 passing 70/10 run from 0.99 Mb to 6.66 Mb, so nothing that passes is near
  them. `--score-everything` turns the filter off. In refinement the filter saves almost nothing, because
  the 50 kb–20 Mb input filter above has already run — it earns its keep when the benchmark scores bin
  sets as a binner emitted them, where 0.7% of bins hold 11.8% of all the sequence. An N50 floor was
  implemented, measured and dropped: the worst N50 among passing bins is 2,025 bp, so a floor low enough
  to be safe skips nothing the size bounds had not, and one high enough to skip anything discards real
  genomes.

  Both leave one set in `binsO`, and everything after that - re-scoring with CheckM, the `-c`/`-x`
  cut, the plots, the final layout, GUNC - is the same code either way. That is what makes them
  comparable: only the route to the one set differs, not what counts as a final bin. `metawrap`
  is the default because it is the incumbent, and the `consolidator` comparison
  (`metawrap2 benchmark run cami2_gi consolidator`) is what would change that, on NC MAGs recovered.
  A tie on those is a win for DAS_Tool, which is 7× cheaper.
- **reassemble_bins** - reassemble bins to improve N50/completion and reduce contamination.

  **`--max-coverage`** (default **50**, `0` to disable) thins each bin's recruited reads before SPAdes
  sees them. Read recruitment is per bin, so an abundant organism hands SPAdes hundreds of times the
  coverage short-read assembly can use - and SPAdes time grows with the number of reads, not with
  genome size, so a handful of high-abundance outliers decide how long the whole step takes. 50x is
  comfortably above where prokaryotic short-read assembly stops improving, which makes the thinning
  close to free in assembly quality and not close to free in time.

  The target is measured against **that bin's own assembled length**, not an assumed genome size: a
  1.2 Mbp bin and a 6 Mbp bin want very different read counts for the same depth. A bin already below
  the target is not touched, which is the common case. Reads are dropped in **pairs**, with a fixed
  seed, so a re-run keeps the same reads and SPAdes is never handed a library with one mate missing.
- **dereplicate** - collapse several samples' bins into one non-redundant genome catalogue.

  Bin twenty samples separately and the same abundant organism comes back twenty times. Nothing
  downstream knows that: a MAG count over the pooled folders counts bins, not organisms, and a
  quantifier handed near-identical genomes splits one organism's reads between them, so every copy
  looks rarer than the organism is. `metawrap2 dereplicate -b sampleA/metawrap_70_10_bins -b
  sampleB/... -o catalogue/` takes any number of bin folders and returns one representative per
  cluster in `dereplicated_genomes/`, symlinked (or copied with `--copy`).

  Identity comes from **skani** (`triangle -E`), which sketches rather than aligns and so scales to
  the thousands of pairs a cohort produces; everything after that - the clustering and the choice of
  representative - is Python, so it is unit-tested rather than inferred from a tool's behaviour.
  `--ani 95` is the default, the threshold large MAG catalogues use for a species, which is a
  convention that makes the count comparable with a published one rather than a measurement of ours.
  `--ani 99` separates strains instead. A pair also has to share `--min-aligned-fraction` of a genome
  (10% by default, deliberately permissive): 99% identity over 5% of two genomes is a shared prophage,
  not one organism. Raising it makes clusters smaller and the catalogue larger, which is a real choice
  and not a correctness fix.

  The representative is the **best-scoring** member, not the largest, read from the
  `<bins>.stats.tsv` that `bin_refinement` already leaves beside each folder (completeness minus five
  times contamination, size breaking ties). A scored genome beats an unscored one however large,
  because size says a bin is big and not that it is one organism - a 6 Mb chimera must not stand in
  for the clean 3 Mb copy of the same genome. With no CheckM table anywhere the module says so and
  falls back to size, which is when a chimera can win its cluster; running
  `bin_refinement`/`binning --run-checkm` first is the fix.

  Clustering is **single linkage**, so it chains: A-B and B-C at 96% merge all three even when A and
  C are not 96% identical. That is accepted because it never splits a genuine duplicate, which is the
  failure dereplication exists to prevent, and because `clusters.tsv` records each member's
  `ani_to_representative` - a member pulled in by chaining has no value there, which is what an
  over-merged cluster looks like from the outside. `clusters.tsv` has a row per input bin, kept or
  not, with its cluster, source folder, representative, quality score and size. One env, once:
  `metawrap2 install-env dereplicate`.
- **quant_bins** - estimate bin abundance across samples, with four quantifiers to choose from.
  sylph by default, which never aligns anything: it sketches the bins and the reads and profiles one
  against the other. Measured on 76 bins from the CAMI2 GI gold-standard assembly across three
  samples, against a read-level truth, it ordered the bins better than either aligner (Spearman
  0.868 vs salmon's 0.827 and CoverM's 0.803) for **0.19 CPU-hours against salmon's 1.08 and CoverM's
  2.37**, and peaked at 1.3 GB where CoverM peaked at 11.6 GB. What it does not do as well is
  magnitude: salmon's estimates correlate much better linearly (Pearson 0.959 vs 0.811), so
  `--quantifier salmon` is the one to reach for when the actual numbers matter rather than the
  ranking. `--quantifier coverm` maps against the bins directly and was both the slowest and the
  hungriest. Re-measure any of this with `metawrap2 benchmark run cami2_gi quantifier`.

  **`--quantifier depth --depth TABLE`** is the cheapest of the four and the only one that does not align
  anything a second time. The other three answer a question the binning run already answered, and answer
  it on a different scale; this aggregates binning's own per-contig depth table into a length-weighted
  mean per bin — total bases over total length — with the per-contig median beside it. **0.013 CPU-hours**,
  against 0.126 for sylph and 2.37 for CoverM.

  Measured against the same read-level truth, it splits the two metrics: Pearson **0.953** and L1
  similarity **0.799** put it second of four on magnitude, ahead of CoverM on both; Spearman **0.780**
  puts it last on ordering, which is the metric to lean on. That split is coherent rather than noise —
  Pearson is dominated by the abundant bins while Spearman weights every bin alike, so a method slightly
  noisy on low-coverage bins reads well on one and badly on the other. Reach for it when the magnitudes
  matter and the table is already on disk; `--quantifier sylph` still orders bins best, and
  `--quantifier salmon` is still the most accurate on magnitude outright.

  `--min-ani` is sylph's one knob worth knowing about, and it is the setting that moved the numbers
  most. At sylph's own default of 95%, more than half the abundance reported in one sample went to
  four bins whose genomes are absent from that sample - relatives close enough to match at a
  species-level threshold. The default here is 99; at 99.9 the floor starts discarding bins that
  really are present. Bins off a fragmented real assembly match their own genome less tightly than
  these gold-standard ones did, so lowering it towards 95 is reasonable there - below 95 sylph's own
  statistics stop being calibrated and MetaWrap2 says so.

  All three write `bin_abundance_table.tsv` and the same heatmap, but **the cells do not mean the
  same thing**: salmon's is a length-weighted median contig coverage, CoverM's is relative abundance
  as a percentage of the sample, sylph's is a k-mer sequence abundance, also a percentage. A CoverM
  or sylph row sums to at most 100 and a salmon row has no such ceiling. Don't compare tables from
  two different quantifiers. A zero differs too: CoverM's means "mapped, but under
  `--min-covered-fraction` of the genome was covered", sylph's means "not detected above `--min-ani`".
  sylph and salmon both ship in the module's own env; CoverM needs its own, once:
  `metawrap2 install-env bins`. Any of the three accepts gzipped bins, and ignores
  non-FASTA files in the bin folder, so a `-b` straight from `bin_refinement` works.
- **blobology** - GC-vs-abundance blobplots of contigs and bins.

  bowtie2 is gone. It was the module's only aligner, and blobology maps reads back for exactly the
  reason binning does - per-contig abundance, here for the plot's y-axis - so the same measurement
  settles it: bowtie2 cost what bwa cost (7.44 CPU-hours, 2.58 of it building the index) and agreed
  with bwa less closely than minimap2 did. `--aligner` now offers the same three choices binning
  offers, from the same registry in `metawrap2.aligners`, defaulting to **minimap2**. One registry
  rather than two means the modules cannot drift apart on a question they share.

  No sort, deliberately: `gc_cov_annotate` reads the BAM with a plain `samtools view` and sums
  aligned span per contig, so record order is irrelevant and sorting would be work spent on nothing.
  The alignment also used to go through a SAM on disk that nothing ever deleted, leaving a full
  uncompressed alignment per sample behind permanently; it now streams, and the cleanup step sweeps
  the BAMs and the aligner's index at the end of the run.
- **classify_bins** - conservative taxonomy assignment for bins.
- **annotate_bins** - functional annotation of draft genomes.

Each module runs in its own conda environment `metawrap2-<module>` (see
[Installation](installation.md)).
