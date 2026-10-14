# MetaWrap2: new capabilities, measured against metaWRAP 1

MetaWrap2 is a Python-native rewrite of metaWRAP. This document covers what it can do that metaWRAP 1
cannot, and what each capability is worth when measured. metaWRAP 1.3.0 is the baseline throughout.

| | metaWRAP 1.3.0 | MetaWrap2 |
|---|---|---|
| cohort runs | one sample at a time, by hand | `run` / `cohort` from a sheet, resumable per step |
| cross-sample coverage | none | `coverage`: one depth table, N alignments instead of N² |
| dereplication | none | `dereplicate`: a cohort's bins collapsed to organisms |
| binners | metaBAT2, MaxBin2, CONCOCT | those three plus SemiBin2, VAMB, COMEBin, GPU-capable |
| which binners run | one fixed set | four `--preset` tiers, 35 to 176 clean MAGs |
| depth table | `jgi_summarize_bam_contig_depths` | own implementation; identity floor measured and moved |
| chimera detection | none | GUNC on every final bin set, on by default |
| consolidation | CheckM on every combination | same, minus repeat and ruled-out calls |
| bin abundance | salmon | sylph default, four quantifiers measured |
| read capping for reassembly | none | `--max-coverage`, default 50× |
| defaults chosen by | the author | `metawrap2 benchmark`, results committed |
| install | one env, manual | 9 workflow envs, `install-env` / `install-db` / `doctor` |

---

# Methods

Data. CAMI2 human microbiome challenge, gastrointestinal, short read. Three samples (0, 1, 10) unless
stated. CAMI2 ships a pooled gold-standard assembly and a contig-to-genome truth table, which makes both
designs below possible.

Two designs. Results from one are not comparable with the other, and no figure here mixes them.

- End-to-end arms give each pipeline the same raw reads and let it assemble its own. This answers "what
  does a user get", and it is how metaWRAP 1 is compared. AMBER cannot score these: its truth is written
  against CAMI2's pooled gold assembly, and self-assembled contig names do not appear in it.
- Binner arms bin the gold assembly directly. This isolates the binner from assembler variance and
  allows AMBER.

Quality. CheckM2 scores bins into MIMAG tiers (near-complete ≥90% complete / ≤5% contaminated,
high-quality ≥70/≤10) and GUNC scores chimerism. The headline metric is clean high-quality MAGs: bins
that are both ≥70/≤10 and non-chimeric, counted by joining the two reports per bin. It is not the
product of an HQ count and a pass fraction. GUNC runs per bin set, so a pass rate over raw binner output
is diluted by bins that were never candidates; multiplying by it understated one arm by 18 genomes.

The same CheckM2 and GUNC score every arm of a comparison. metaWRAP 1 refines with CheckM 1.0.18, which
reports different completeness for the same genome, so letting each arm report its own numbers would
compare the scorers.

Where AMBER applies, it adds bp-weighted F1, purity and completeness against the truth table. AMBER
purity is bp-weighted, so a bin that is 95% one genome and 5% another barely registers. It is reported
alongside GUNC, never instead of it.

Cost. CPU-hours, measured by reaping the child with `os.wait4`. Wall time is not used, because it
rewards a tool for being given more cores. A module run outside the harness produces bins and no cost,
so it is not a measurement. Cost is reported beside quality and never divided into it: a
quality-per-cost ratio is dominated by the cheapest arm, and ranking on one produced a default that
recovers fewer genomes than metaWRAP 1.

Harness. `metawrap2 benchmark list|fetch|run|score`. Comparisons are declared as data in
`benchmarks/configurations.toml`, so adding an arm is a TOML entry. Results land in a committed TSV, so a
moved default arrives in review next to its numbers. `fetch` is checksum-verified, never implicit, and
refuses to cross a disk ceiling.

Baseline. metaWRAP 1.3.0 from its own conda environment at the versions its recipe pins: metaBAT2
2.12.1, MaxBin2 2.2.5, CheckM 1.0.18, bwa 0.7.18, Python 2.7.15.

---

# Results

## 1. Cohort support

metaWRAP 1 processes one sample at a time, by hand. MetaWrap2 adds three pieces that make a cohort a
first-class input. This is the largest functional gap between the two.

`metawrap2 run` / `cohort` execute a study sheet: samples, read paths, steps, per-module options. State
is a per-step manifest fingerprinted on inputs, so `--resume` reruns only the steps whose inputs changed.
A 90-CPU-hour study that dies in refinement resumes at refinement.

`metawrap2 coverage` gives each contig's depth in every sample, which is what separates two genomes at
similar abundance in any one sample. It needs N alignments instead of N² and writes one table every
binner reads.

![Cohort designs](figures/cohort.png)

The top two bars are the comparison; both give the binner one depth column per sample. Co-assembly with
every sample's coverage is metaWRAP 1's standard cohort workflow, so it is the baseline.

| | clean MAGs | near-complete | genomes | CPU-h |
|---|---|---|---|---|
| co-assembly, 3 coverages *(metaWRAP 1's workflow)* | 34 | 18 | 43 | 85.1 |
| per-sample assembly, 3 coverages *(`coverage`)* | 34 | **23** | **46** | 90.9 |

Clean high-quality MAGs are level. Per-sample assembly gives 23 near-complete against 18 and 46
organisms against 43, for 7% more CPU, and keeps each sample's bins attributable to that sample.

Scaling is untested. Co-assembly's assembler sees 3× the reads here and would see 20× at twenty samples,
where it is expected to fail on memory. For small cohorts co-assembly remains reasonable; set
`coassemble = true` in a sheet.

The bottom bar is an ablation with one coverage column instead of three: 29 clean MAGs against 34, and
19 near-complete against 23. Building the table costs 2.26 CPU-hours and returns 1.2, because binning
then reads a table instead of aligning.

`--coverage fairy` (k-mer sketches instead of alignments) makes the coverage step 53× cheaper and loses
13 of 34 clean MAGs. That step is 2.5% of the pipeline, so fairy is not the default.

`metawrap2 dereplicate` collapses a cohort's per-sample bin sets into one catalogue of organisms. Every
cohort figure here is post-dereplication.

## 2. Against metaWRAP 1, end to end

![metaWRAP 1 vs MetaWrap2](figures/mw1_vs_mw2.png)

| arm | high-quality | clean | CPU-h *(binning + refinement)* | GPU-h |
|---|---|---|---|---|
| `--preset very-fast` | 39 | 36 | **8.8** | — |
| `--preset fast` | 45 | 41 | 13.8 | — |
| **metaWRAP 1.3.0** | **52** | **48** | **52.5** | — |
| `--preset default` *(shipped)* | 56 | **53** | **43.8** | 0.6 |
| `--preset slow` | 64 | **60** | 120.7 | 2.4 |

Cost is binning plus refinement. The tiers share one set of assemblies, and assembly alone is 83
CPU-hours of every total; including it compresses a 6× spread into 1.8×. CPU-hours come from the driver
wrapping every step in `os.wait4` (`bench_work/preset_cost.py`). GPU-hours, the wall span of binners
given `--device cuda` (`bench_work/preset_gpu_hours.py`), are listed separately: a GPU-fed binner still
keeps ~40 host cores busy, but without a GPU `default` and `slow` are not available at any CPU price.

- `default` returns 10% more clean MAGs than metaWRAP 1 for 17% less CPU. It became the shipped tier on
  2026-10-09.
- `fast`, the previous shipped tier, gives up 15% of the clean MAGs for 74% less CPU. It was chosen on a
  quality-per-cost ratio.
- `slow` returns 25% more clean MAGs for 130% more CPU. It refines three bin sets: 7 hybrid combinations
  against `default`'s 3.

A plain `metawrap2 run` therefore needs `metawrap2-learned`, ~1.5 GB of pytorch, which `quickstart`
builds. `--preset fast` gives the previous behaviour without it.

These are the numbers for "what do I get". The larger figures in §3 come from CAMI2's gold assembly; on
a real per-sample assembly the tiers sit closer together.

![Cost per step](figures/mw1_steps.png)

Almost all of the cost difference is binning, which falls from 38.4 CPU-hours to 1.8: minimap2 instead
of bwa, and two binners instead of three. Refinement follows, since two bin sets give 3 combinations to
score with CheckM and three give 7.

![Share of CPU by step](figures/cost_landscape.png)

Assembly is 86% of a default run and does not depend on the binners, which is why the tier comparison
leaves it out.

## 3. Binning

![Binner sets](figures/binner_sets.png)

Fourteen arms on the gold assembly, every count post-refinement. metaWRAP 1's trio sits mid-table. The
top tier recovers 1.55× the shipped tier's clean MAGs here and 1.13× end to end (§2). Gold-assembly
contigs are long and clean, so this design measures what a binner can do; §2 measures what it does for a
user.

![Single binners](figures/binners_alone.png)

On agreement with the truth table, both learned binners beat all three of metaWRAP 1's. CONCOCT scores
third and costs more than everything else combined, so it is opt-in.

### `--preset`

![Presets](figures/presets.png)

| tier | binners | clean MAGs, end to end | clean MAGs, gold assembly | extra install |
|---|---|---|---|---|
| `very-fast` | metaBAT2 | 36 | 35 | — |
| `fast` | metaBAT2 + MaxBin2 | 41 | 59 | — |
| `default` *(shipped)* | metaBAT2 + SemiBin2 | **53** | 114 | `install-env learned` |
| `slow` | metaBAT2 + SemiBin2 + COMEBin | **60** | **177** | `install-env learned` |

The better binners pull away by 3× on a perfect assembly and by about 1.4× on real per-sample
assemblies. A preset never overrides a binner named explicitly.

The tiers are not nested: only `fast` runs MaxBin2. MaxBin2 is the dirtiest input measured (37.7% of its
bins pass GUNC), so the higher tiers leave it out. The quality ordering holds across all four.

In `default`, metaBAT2 reads the depth table SemiBin2 already built and costs ~1.3 CPU-hours. On the gold
assembly the pair gives 114 clean MAGs against SemiBin2's own 116; end to end it is worth keeping.

`metawrap2-learned` is one environment holding SemiBin2, VAMB and COMEBin, which are all pytorch
applications that solve together. Its solve pulls a CUDA-capable pytorch, so it serves both
`--device cpu` and `--device cuda`; there is no separate GPU install.

### Which binners complement each other, and what `slow` runs

Refinement keeps the best version of each bin, so a binner that finds 30 genomes nobody else finds is
worth more in a trio than one that finds 100 the others already had. To measure this, each clean
high-quality bin is assigned to the reference genome owning most of its sequence, and unions are taken
over pairs and trios (`bench_work/binner_complementarity.py`, on the existing gold-assembly bin sets).

| binner | clean MAGs | genomes recovered | of those, found by nothing else |
|---|---|---|---|
| COMEBin | 176 | 176 | **39** |
| SemiBin2 | 116 | 116 | 4 |
| VAMB | 80 | 80 | 0 |
| metaBAT2 | 35 | 35 | 0 |
| MaxBin2 | 33 | 33 | 1 |
| CONCOCT | 28 | 28 | 0 |

The binners are highly redundant here, and COMEBin contains most of the others. The best pair is
COMEBin + SemiBin2 at 181 genomes; no trio exceeds 182. The best trio is COMEBin + SemiBin2 + MaxBin2
(182); COMEBin + SemiBin2 + metaBAT2 reaches 181.

`slow` uses the metaBAT2 trio. MaxBin2 costs 24.5 CPU-hours here against metaBAT2's ~1.5, and only 37.7%
of its bins pass GUNC. Of the trios within one genome of the best, metaBAT2's is the cheapest by an order
of magnitude, against VAMB's 4.9 GPU-hours, MaxBin2's 24.5 and CONCOCT's 238.

The trio spans three algorithm families: contrastive (COMEBin), self-supervised (SemiBin2), and
composition-plus-coverage (metaBAT2). The union of 182 is an upper bound that assumes refinement picks
the right version of every bin, and one dataset cannot show whether COMEBin's 39-genome edge
generalises. Three binners that fail differently give the top tier a floor a single binner cannot.

### `--device`

![GPU arms](figures/gpu_arms.png)

One switch covers all three learned binners, including VAMB. `auto` asks each binner's environment
whether torch can see a device. The default is `cpu`, so re-analysing a study on different hardware
gives the same MAGs. For `slow`, a GPU is the difference between 40 CPU-hours and not finishing.

VAMB refuses to start on the gold assembly, correctly: it carries 20 records with no A, C, G, T or U
(masked reference genomes, 28.9 Mbp, the largest 6.6 Mbp). metaBAT2, MaxBin2 and CONCOCT accept such a
contig and cannot place it. `binning` and `coverage` drop them.

## 4. Coverage estimation

MetaWrap2 computes its own depth table instead of calling `jgi_summarize_bam_contig_depths`. The
algorithm is the same and agrees with jgi at Pearson 0.9992 on means across all 210,632 CAMI2 contigs.
The reason is auditability: MaxBin2 was found binning on depths jgi had invented, which showed only after
running the C++ binary twice and comparing. Owning the code made one inherited parameter testable.

![Identity floor](figures/identity_sweep.png)

jgi discards reads below 97% alignment identity. On strain-rich data that costs a quarter of the clean
MAGs: at a 90% floor bins are more complete, more pure and less often chimeric, with 47 clean MAGs
against 35. The risk in loosening it is that a relative's reads inflate a contig's depth and the binner
merges two strains; AMBER purity and GUNC show that this does not happen here.

The floor is 90, MetaWrap2's only departure from jgi's parameters. `--min-identity 97` reproduces a jgi
table, as reproducing a metaWRAP 1 result requires. The evidence is one dataset and one binner.

![Coverage fidelity](figures/coverage_fidelity.png)

The cheap alternatives were checked per contig. CoverM reads 2× bwa because it applies no identity
filter, so it is not a drop-in for a depth table.

## 5. Refinement

![Refinement savings](figures/refinement_savings.png)

Refinement hybridises every combination of its input bin sets and asks CheckM which version of each bin
to keep: 2ⁿ−1 CheckM runs. Two parts of that work are removed without changing an answer:

- A third of the calls are a bin already scored, identified by a hash of its sorted contig ids.
  Identical contigs give identical completeness.
- Bins outside 200 kb–20 Mb cannot be ≥70%-complete genomes, and CheckM's cost follows bases: metaBAT2's
  unbinnable remainder reaches 159 Mb. `--score-everything` disables this filter.

GUNC runs on every final bin set by default. A bin built from halves of two genomes has a complete,
uncontaminated-looking marker set that CheckM accepts. On one arm 4 of 60 final bins are chimeric by GUNC
despite passing CheckM at 70/10. Bins are annotated, never dropped. `--skip-gunc` opts out, along with
its 6.7 GB database. On the same arm the two input bin sets were 12.0% and 36.0% chimeric and the
consolidated set is 6.7%.

`--consolidator dastool` scores input bins once on single-copy genes instead of running CheckM on every
combination. On one sample's per-sample bins, `metawrap` keeps 13 MAGs for 4.52 CPU-hours and DAS Tool
keeps 10 for 0.62. Three genomes on one sample is not enough to move a default.

## 6. Bin abundance

![Quantifiers](figures/quantifiers.png)

metaWRAP 1 quantifies with salmon. MetaWrap2 offers four quantifiers and defaults to sylph. None is best
at both ordering and sizing bins: sylph orders best and sizes worst, salmon the reverse. Abundance tables
are usually read for which MAGs dominate and how that shifts across samples, which is an ordering.

`--quantifier depth` aligns nothing: it takes a length-weighted mean per bin from the depth table binning
already wrote. It costs 0.013 CPU-hours, 180× less than CoverM, and ranks second of four on magnitude and
last on ordering. It is the only quantifier on the scale the binner clustered on.

Values are not interchangeable across quantifiers:

- salmon reports a length-weighted median contig coverage, CoverM a percentage of the sample, sylph a
  k-mer sequence abundance.
- A CoverM zero means "mapped, but under `--min-covered-fraction` covered"; a sylph zero means "not
  detected above `--min-ani`".
- sylph's `--min-ani` defaults to 99, not sylph's own 95. At 95, more than half the abundance in one
  sample went to four bins whose genomes are absent from it.

## 7. Smaller changes

`read_qc` trims with fastp: 1.0% more bases and 0.4% more Q30 bases kept for 5.8% less CPU (cutadapt
stops scaling at 4 threads). Trim Galore remains supported and suits memory-constrained machines: it
peaks at 0.1 GB, fastp at 1.4.

`reassemble_bins --max-coverage` (default 50) thins each bin's recruited reads to 50× of that bin's
assembled length, in pairs, with a fixed seed. An abundant organism otherwise hands SPAdes hundreds of
times the coverage it can use, and SPAdes' time grows with reads. `0` keeps every read.

`blobology` no longer uses bowtie2. `--aligner` offers binning's three choices from one shared registry.

`mw2` is the same CLI under a shorter name.

Nine conda environments. Modules that run in sequence and share tools share an environment. Four stay
separate because their solutions conflict: CheckM2 has no common solution with CheckM1 + GUNC + DAS
Tool, and GTDB-Tk none with taxator-tk. GUNC, DAS Tool, CoverM, fairy and Bakta are installed by default.

Every command prints its defaults and accepted values from one formatter across two dozen entry points.
Defaults are read off the parser, and a numeric option declares its range once for both checking and
printing: `-c/--completeness` reads `[0-100] [default: 70]` and rejects `-5`.

## 8. Defects found and fixed

Each produced wrong numbers or lost runs without an error, and each has a regression test.

Pipeline:

- MaxBin2 binned on depths `jgi_summarize_bam_contig_depths` invented, affecting every bin it produced.
- An all-N contig broke VAMB and was binned by everything else.
- Four modules treated the unbinned leftovers as a candidate genome.
- `--resume` on `bin_refinement` deleted what it was resuming from and reported success; the same two
  defects were in `reassemble_bins`.
- metaBAT2 had no seed, so identical input did not reproduce.
- A crashed aligner read as a successful run.
- With one input bin set nothing applied the quality thresholds: refining a single COMEBin run left 376
  bins in `metawrap_70_10_bins`, of which 178 pass.
- `reassemble_bins` recruited reads that are not reads and scored unknowns as perfect.
- A study sheet selecting any `--preset` but `fast` failed at refinement: the driver resolved binner
  flags but not presets, and passed refinement two directories that were never created.

Harness:

- No multi-binner arm had been ranked: refinement writes into a subdirectory and the scanner only read
  the variant root.
- A re-score that changed an arm's bin set was added to the old one, so COMEBin's 178 MAGs read as 356.
- Clean-MAG counts were estimated as `HQ × GUNC pass fraction`, which understated COMEBin by 18 genomes
  and reversed a default.
- `benchmark score` run from outside the repository wrote a second results file and reported success.
- AMBER scored one bin set twice as two arms, and was skipped on any variant already through CheckM.
- The ranker summed fractions, so the variant measured on the most samples won.

Found in review of the harness after the results above were produced. Three of these move numbers in the
committed results table; the affected figures are marked and need a re-score.

* `read_qc` figures are one mate's count. The scorer wrote one row per mate file, but a measurement is
  keyed on (dataset, comparison, variant, sample, metric), so the two rows collapsed to the last mate
  walked. **`reads_kept`, `bases_kept` and `q30_bases_kept` in §7's trimmer comparison are half the
  real counts.** Both arms are halved identically, so the 1.0%/0.4% margins still hold; the absolute
  numbers do not.
* A refinement arm used only three of its bin sets, alphabetically: the runner zipped the bin
  directories onto `("-A", "-B", "-C")` and `zip` truncates. With
  `comebin/concoct/maxbin2/metabat2/semibin2/vamb` that drops metaBAT2. It also fed a previous
  refinement's output back in as an input set. **Any multi-binner refinement arm above four binners has
  to be re-run**; the arms in §5 use two.
* A crashed variant's partial output was scored and published. Scoring ignored the recorded non-zero
  exit, so bins from a variant that died after its first binner were ranked and listed with cost
  "unmeasured", indistinguishable from a finished arm whose rusage was lost. Failed variants are now
  skipped; `--score-failed` overrides.
* Abundance was scored over the intersection of the table and the truth, so a quantifier that omitted a
  bin was scored on a smaller, easier set. The salmon path writes only rows salmon emitted, while the
  CoverM path pre-seeds one per bin, so arms could be scored over different numbers of bins. Missing
  bins now count as zero and are reported as `abundance_missing_bins`. **§6's ordering result should be
  re-scored**; the gap between sylph and salmon is wide enough that it is unlikely to turn over.
* The depth stream could deadlock. `samtools view | samtools depth` left both stderr pipes unread until
  the depth stream was consumed. One warning per record fills the 64 KiB pipe buffer, samtools blocks,
  and the run hangs with no output and no timeout. Both pipes are now drained concurrently.
* The disk monitor could abort the wrong step. It interrupted the main thread without re-checking
  whether the watched step had finished, so under `run` or `cohort` the `KeyboardInterrupt` landed in
  the next step, and the watched step was reported as a disk stop it never hit.

## 9. Limitations

One dataset. Everything here is CAMI2 gastrointestinal, mostly three samples. The identity floor is the
change most worth re-testing elsewhere. A looser floor should matter less on a self-assembled contig set,
because a fragmented consensus already absorbs some of the strain variation the true genomes do not.

Open, with the diagnosis recorded:

- CONCOCT costs ~4× more under MetaWrap2 than under metaWRAP 1 for the same version, the same flags,
  and assemblies agreeing to one contig in 19,054. CONCOCT forks `-t` workers and each one's OpenBLAS
  started its own 96 threads; capping the numerical libraries to one thread per forked worker recovered
  1.3×. The rest is unexplained. The next suspect is the build: Python 3.12 and scikit-learn 1.9 against
  Python 2.7 and a much older scikit-learn. It does not affect the default path, since CONCOCT is
  opt-in.
- `metabat2_comebin`, the cheap half of the three-binner ensemble, is unmeasured.
- A dereplication comparison needs a scorer that does not exist.
- The cohort-scaling experiment that would settle co-assembly has not been run.
- MetaDecoder and binny are not started.

Decisions outstanding:

1. Should a GPU machine default differently? VAMB on a GPU beats the shipped tier on every axis
   including cost; COMEBin on a GPU is the best catalogue measured. Neither is a default, because both
   depend on hardware most users do not have.
2. Is moving the identity floor too aggressive on one dataset? It reverts in one line. Related and
   unsettled: the floor uses `qlen` as its denominator, which counts soft-clipped and inserted bases,
   where jgi divides by the aligned length alone, so the effective floor is looser than the number
   requested. It agrees with jgi to Pearson 0.9992 on the means, so the aggregate effect is small, but
   tightening it moves every depth table and needs re-measuring.

Decided since this was written: `--preset default` is the shipped tier (§2), because end to end it
beats metaWRAP 1 on both quality and cost.

State. The suite passes with ruff, black, isort and mypy clean over 91 source files, and `doctor`
reports 12/12 modules healthy.
