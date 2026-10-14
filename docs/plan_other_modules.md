# Plan: everything else

The companion to [plan_functional_annotation.md](plan_functional_annotation.md), which covers
annotation, pathways, AMR and abundance. This plan covers the other modules and the domains MetaWrap2
does not cover yet.

Each section starts with its status; §5 has the state of every phase. "No status recorded" means this
plan does not say whether the item is done. Numbers attributed to other tools are evidence about those
tools, not targets for MetaWrap2.

---

## 0. Where MetaWrap2 stands

Context.

Galaxy's training material describes Binette as "the successor of metaWRAP". Binette is faster, handles
more bin sets, and uses CheckM2, so bin refinement alone is not a differentiator.

MetaWrap2's strengths: a measurement harness (`benchmark run`/`score`, with CPU-hours from `os.wait4`),
provenance on every number, and one path from reads to MAGs to functions.

Where the peer set is ahead:

| Pipeline | What it has that we don't |
|---|---|
| nf-core/mag | long reads, multi-split binning, Nextflow resumability |
| ATLAS | Snakemake, strain tracking, three-point abundance |
| SqueezeMeta | functional abundance for reads, contigs and bins; R bindings |
| VEBA | viruses, pangenomics, BGC novelty, in one pipeline |
| Binette | faster, better-scaling refinement |
| MGnify / mettannotator | mobilome, versioned databases, validated on MAGs |

---

## 1. Things that change results

### 1.1 Multi-split binning

Status: done and measured, 2026-10-13 (§5, phase 2). It does not beat the shipped cohort path, so
`--binsplit` is a flag, not a cohort default. Open: the same experiment with the learned binners. No
status recorded for SemiBin2's native multi-sample mode.

MetaWrap2 bins one assembly at a time; a cohort gives each sample's binning run a cross-sample depth
table. Multi-split instead bins the concatenated catalogue, so the model sees every sample's composition
at once, then splits the result back out per sample. Söding's analysis reports that splitting the
embedding space by sample beats splitting the final clusters: 19% more bins and 79% more high-quality
bins. nf-core/mag has this mode.

Two results, in opposite directions:

- Splitting a catalogue bin set pays if you bin a catalogue. metaBAT2 goes from 86 bins and 16 clean
  high-quality MAGs to 166 and 37; MaxBin2 from 114 and 15 to 258 and 19, for identical cost (3.29
  against 3.28 cpu-hours; splitting is a pass over bins already written).
- Binning the catalogue does not pay against the shipped cohort path (same three assemblies, same
  cross-sample depth table, per-sample binning). Multi-split gives 47 genomes against 46 and 38
  high-quality bins against 37, but 29 clean high-quality MAGs against 34, for 9.09 cpu-hours against
  8.71. The final catalogue's GUNC pass fraction falls from 87.0% to 70.2%: one cluster over a
  concatenated catalogue can gather near-identical contigs from different samples, and splitting
  afterwards does not unmix them.

The likely reason the published gain is absent: per-sample binning on a cross-sample depth table already
has the cross-sample information. The paper compares against single-sample coverage, which is our arm b,
and that arm also returns 29 clean MAGs. `--binsplit` remains necessary for anyone binning a catalogue,
including SemiBin2's and VAMB's native multi-sample modes.

Caveats: one dataset, three samples, two classical binners. The paper's claim is about the embedding, so
the next experiment is the same comparison with the learned binners.

Machinery: `coverage.py` builds the prefixed concatenated catalogue. `binning.py` passes VAMB `-o ''` to
disable binsplitting for a single assembly, since there is no sample prefix and VAMB's default separator
`C` would shred bins on any contig id containing a capital C. The planned work was:

- A cohort binning mode that bins the concatenated catalogue and splits bins back per sample, giving
  VAMB the prefix separator `coverage.py` already writes.
- The same for SemiBin2's native multi-sample mode.
- A benchmark against the existing cohort arms (a re-bin, not a re-assembly).

### 1.2 Bin quality numbers that are wrong today

Status: open per this section; phase 1 in §5 has no recorded state. Real MIMAG tagging is halted on user
instruction and should not resume without direction.

- GUNC ≤1.0.6's `summarise` flagged every genome as passing; fixed in 1.1.1. `envs/binqc.yaml` asks for
  `gunc >=1.0.6`, so any arm scored with the floor version has meaningless chimerism numbers. Pin
  `>=1.1.1` and assert the version at runtime.
- CheckM2 ≤1.0.2 gives batch-dependent scores: the same genome scores differently depending on the rest
  of the batch. `envs/checkm2.yaml` asks for `>=1.0.2`. Pin above it and record the version per run.
- `scoring.QUALITY_TIERS` calls its tiers "MIMAG-style", and one is not. MIMAG high-quality is >90%
  complete, <5% contaminated, plus 23S, 16S and 5S rRNA and at least 18 tRNAs; medium is ≥50/<10, which
  `bins_medium_quality` matches. There is no 70/10 MIMAG tier, so `bins_high_quality` is a working
  threshold with a borrowed label, and `bins_near_complete` (90/5) is MIMAG high-quality without the
  RNA requirement.

For the labels, fix the comment and the docs, not the field name: `bins_high_quality` is the key in every
committed row of `benchmarks/results/`, and renaming it would invalidate those rows or need a migration.

Real MIMAG tagging would be nearly free: Bakta already runs tRNAscan-SE, Aragorn and Infernal/Rfam, so
the rRNA and tRNA counts are a parse of existing output.

### 1.3 `dereplicate` and the skani unit bug: checked, not affected

Status: done (checked; no change needed).

dRep v4 fixed a skani alignment-fraction percent-versus-fraction confusion that merged two distinct
*Enterococcus* species into one cluster. A factor-of-100 error there collapses real species and does
not show in the output. `MIN_ALIGNED_FRACTION = 10.0` is passed to
`skani triangle --min-af`, which takes a percentage (skani's own default is 15), and the constant's
docstring says "as a percentage".

### 1.4 CoverM's identity floors are unset

Status: open per this section; phase 1 in §5 has no recorded state.

`quant_bins.COVERM_GENOME` sets `--min-covered-fraction 10` and pins `--mapper`. It does not set
`--min-read-percent-identity` or `--min-read-aligned-percent`, whose CoverM defaults are 0, so every
alignment counts. For MAG abundance that inflates coverage with cross-mapping reads from relatives, the
failure the 90% identity floor in `depth.py` addresses for the depth table. Set explicit floors, document
and record them in provenance, and let the quantifier benchmark say what the floor is worth.

### 1.5 `classify_bins` defaults are the wrong way round

Status: taxator-tk is deprecated with a per-run warning; the default is not flipped (§5, phase 4). The
GTDB-Tk version and reference pins are open per this section. §5 refers here for why flipping the
default needs a measurement first, but this section does not give that reason.

taxator-tk is abandoned and is the default; GTDB-Tk is opt-in. Make GTDB-Tk the default and retire
taxator-tk.

`envs/gtdbtk.yaml` asks for `gtdbtk >=2.4.0`, which resolves to 2.7.2. That needs reference data R232,
about 140 GB of RAM for `classify_wf`, and has removed flags our invocation may still pass. Pin an exact
version and the reference release, state the RAM requirement in `doctor`, and check the flags against
2.7.x. An unpinned `>=` on a tool with a 140 GB floor fails for every user with a smaller machine, and
only after the install.

### 1.6 `reassemble_bins` coverage gating — already handled, and the other half would hurt

Status: thinning done; skipping deep bins dropped from the plan.

Reassembly helps at low to moderate coverage and wastes compute, or degrades the assembly, at high
coverage. The options are to thin each bin's reads to a target depth or to skip deeply covered bins.

Thinning is implemented: `--max-coverage` defaults to 50×, applied per bin against its own assembled
size, with a fixed seed so a re-run recruits the same reads.

Skipping is not worth adding. A bin at 200× thinned to 50× gets the coverage SPAdes can use and
reassembles well; skipping it would discard that improvement to save compute thinning has already saved.

---

## 2. Domains we don't cover

Self-contained modules, ordered by how often a user needs them.

### 2.1 Viruses — a `viruses` module

Status: done, including CRISPR host linking (`--crispr-hosts`, cctyper; §5, phase 5).

MetaWrap2 bins prokaryotes and discards the phage, in datasets where phage is usually the most abundant
thing present.

geNomad is the caller. It classifies contigs as viral, plasmid or chromosomal in one pass with one
database, so `amr`'s mobility work (annotation plan §6.4) shares the download. On one CAMI2 assembly:
161 viral contigs at ≥2 kb and 344 plasmid contigs, from the same run.

Screening the bins of that sample finds zero viral contigs; screening the assembly finds 161. A phage
contig is usually in no prokaryotic bin, so the assembly is the input for free phage, and bins are for
proviruses in a binned genome. VIBRANT and VirSorter2 are stale relative to geNomad.

CheckV does completeness and contamination and has no real competitor. It also trims host regions off
provirus contigs, so an integrated prophage is reported as a virus plus a flanking chunk of bacterial
chromosome, not as one chimeric contig.

The module produces, by analogy with `binning` and `bin_refinement`:

- viral contigs with geNomad's taxonomy and CheckV's completeness/contamination per contig;
- quality tiers following the MIUViG convention (high-quality ≥90% complete, genome-fragment below), an
  actual community standard, unlike the 70/10 tier in §1.2;
- vOTU clustering at the conventional 95% ANI / 85% AF, so counts mean populations. `dereplicate` does
  the same shape of work for genomes;
- abundance per vOTU per sample through the existing `coverage`/`quant_bins` path;
- host prediction: `cctyper` (§2.4) gives CRISPR spacers from the prokaryotic bins, and a spacer matching
  a vOTU is direct evidence of its host.

Caveat to print beside it: CheckV completeness is estimated against a reference database, so a novel
virus with no relative is reported as a fragment whether or not it is one.

### 2.2 Plasmids and mobile elements

Status: contig classification done (geNomad, 344 plasmid contigs on one CAMI2 assembly, feeding
`amr --mobile-elements`). PlasMAAG open. Blockers: it is not on bioconda (it ships inside the VAMB
repository), and it needs the assembly graph, which `assembly` deletes.

Covered in the annotation plan §6.4: geNomad for plasmid contigs (the same run `viruses` needs),
PlasmidFinder for replicon typing, ISEScan for insertion sequences, mobileOG-db for MGE families.

PlasMAAG (*Nat Biotechnol* 2026) does plasmid-aware assembly-graph binning and is a step change over
classifying contigs. `assembly` keeps contigs and discards SPAdes' `assembly_graph_with_scaffolds.gfa`
and MEGAHIT's intermediates, so the first change is to `assembly`, and it adds a disk cost to every run
for a feature most runs will not use. That trade is for a user to decide: make graph retention a flag,
measure its cost, then wire the binner. The geNomad path is the baseline for the comparison.

### 2.3 Strain resolution

Status: done per §5, phase 7; the original text of this section said "we have nothing".

inStrain 1.11.0 for within-sample strain variation, or TRACS, the only tool evaluated that does not
overestimate SNP distances, which matters for any transmission claim. Both are maintained. ATLAS has
this.

### 2.4 CRISPR, pangenomics, phylogenomics

Status: cctyper done (`viruses --crispr-hosts`). Pangenomics and phylogenomics done per §5, phase 7.

- cctyper for CRISPR-Cas arrays and phage-host links through spacer matching. On five CAMI2 bins, three
  carry arrays (I-B, two I-C) with 54 spacers, and none matches the sample's 161 viral contigs. All 54
  spacers are recovered from their own genomes, so the matcher works; expect low or zero yield for a
  single sample. A match is evidence, not prediction.
- PPanGGOLiN for pangenome graphs across `dereplicate`'s output.
- IQ-TREE for phylogenomics from the marker genes GTDB-Tk already extracts.

---

## 3. Reporting what MAGs leave out

Status: binned fraction and chimerism caveat done (`metawrap2.recovery`); `singlem appraise` open (§5,
phase 3).

Two findings that change what the pipeline reports, not what it computes.

MAGs are chimeric more often than their quality scores suggest. Treichel et al. 2026: among high-quality
MAGs, only 54.5–81.8% represent a single strain, with a mean of 2.3 reference genomes per high-quality
MAG. CheckM2 completeness and contamination do not capture this; GUNC captures part of it, at the
taxonomic level.

Most of the assembly never reaches a MAG. SPIRE, across 92,187 metagenomes: only about 10% of assembled
base pairs end up in medium- or high-quality MAGs, and the unbinned remainder holds an estimated 145
predicted undescribed phyla. The remainder is biased, so a MAG-only report describes a narrow,
non-representative slice of the sample.

In the binning and refinement report:

- the fraction of assembled bp binned, and the fraction binned into clean 70/10 bins;
- `singlem appraise` output: what fraction of the community the bins represent per cell, and which
  lineages have no genomic representation;
- a chimerism caveat beside every quality tier, from the existing GUNC result.

---

## 4. Harness and statistics

### 4.1 Long reads

Status: done per §5, phase 8 (metaFlye, long-read mapping presets, SemiBin2's long-read mode). Hybrid
assembly and polishing have no recorded status.

nf-core/mag has long reads. Scope: long-read assembly (metaFlye), hybrid assembly, long-read-aware
binning, polishing. The largest single piece of work in either plan, to be scoped on its own.

### 4.2 Statistics and interchange

Status: ANCOM-BC and BIOM done (`diff_abundance`, `--biom`; §5, phase 6). No status recorded for
MaAsLin 3.

The pipeline stops at feature tables, and every user writes the same analysis script.

- MaAsLin 3 for differential abundance, with its prevalence/abundance decomposition.
- ANCOM-BC is native in scikit-bio: a pure-Python compositional test with no R dependency. Use it as the
  default.
- BIOM 2.1 output for every feature table (one implementation, shared with the annotation plan).
- Treat relative abundances as compositional throughout; sum-constrained data must not go into a raw
  t-test.

### 4.3 Diversity coverage

Status: done, as `metawrap2 nonpareil`.

A separate module because it needs one read file per sample across the cohort, which fits neither
`read_qc` nor `coverage`. Findings from building it:

- The two kernels disagree on short reads: 50.5% against 83.0% coverage on the same 34 bp gut sample. The
  kernel is chosen from the measured read lengths, once for the whole cohort, so coverage figures are
  comparable across samples.
- The k-mer kernel's effort axis is in k-mers, not base pairs (`AL` = L-k+1); reporting it as bp
  overstates the sequencing done by L/AL, 2.3x on those reads.
- The alignment kernel seeds from the clock: two identical runs gave 83.0% and 84.3%.

The model fit is a Python port of the R package's, verified against it to six significant figures on
Nonpareil's three example datasets.

Diversity coverage, not sequencing depth, confounds cross-sample comparison and is ignored by most
published comparisons. It costs 4 minutes per sample under the k-mer kernel and 16 under the alignment
kernel.

### 4.4 Benchmark harness

Status: multi-split arms done (§1.1). Module completeness against MAG completeness done for one genome
(annotation plan §5.3). No status recorded for database version pinning.

- Multi-split binning arms (§1.1).
- Module completeness against MAG completeness (annotation plan §5.3).
- Pin and record every database version in `benchmarks/configurations.toml`, not only tool versions.
  GTDB R220 against R232, KOfam snapshots and AMRFinderPlus database dates all move results.

---

## 5. Order of work

| Phase | Work | Why here |
|---|---|---|
| 1 | Correctness: GUNC pin, CheckM2 pin, MIMAG labelling, CoverM floors, GTDB-Tk pin and default | Numbers currently reported wrongly. Cheap, no new dependencies |
| 2 | Multi-split binning plus benchmark arms | **done and measured**: it does not beat the shipped cohort path - 29 clean HQ MAGs against 34. `--binsplit` stays a flag; see §1.1 |
| 3 | Binned fraction, `singlem appraise`, chimerism caveat | binned fraction and chimerism caveat **done** (`metawrap2.recovery`); `singlem appraise` open |
| 4 | Retire taxator-tk | **deprecated with a per-run warning**, default not flipped - see §1.5 for why that needs a measurement first |
| 5 | `viruses`: geNomad + CheckV + vOTU clustering + CRISPR host linking | **done**, CRISPR host linking included (`--crispr-hosts`, cctyper) |
| 6 | Statistics, BIOM, Nonpareil | **done**: ANCOM-BC and BIOM (`diff_abundance`, `--biom`), Nonpareil (`metawrap2 nonpareil`) |
| 7 | Strain resolution (inStrain/TRACS), PlasMAAG, pangenomics, phylogenomics | **done** except PlasMAAG, which needs the assembly graph `assembly` currently deletes - see §2.2 |
| 8 | Long reads | **done**: metaFlye, long-read mapping presets, SemiBin2's long-read mode; verified on simulated R10 reads - 72 contigs at N50 57.9 kb from a 2.50 Mbp reference |

Phase 1 is about a week and fixes numbers that are wrong now, so it goes first whatever else is
prioritised.
