# Plan: other modules

Companion to [plan_functional_annotation.md](plan_functional_annotation.md), which covers
annotation, pathways, AMR and abundance. This plan covers the other modules and the domains MetaWrap2
does not cover yet.

Each section starts with its status; §5 has the state of every phase. "No status recorded" means this
plan does not say whether the item is done. Numbers quoted for other tools describe those tools and are
not targets for MetaWrap2.

---

## 0. Current state

Context.

Galaxy's training material describes Binette as "the successor of metaWRAP". Binette is faster,
handles more bin sets, and uses CheckM2.

MetaWrap2 has a measurement harness (`benchmark run`/`score`, with CPU-hours from `os.wait4`),
provenance on every number, and one path from reads to MAGs to functions.

What other pipelines have:

| Pipeline | What it has that we don't |
|---|---|
| nf-core/mag | long reads, multi-split binning, Nextflow resumability |
| ATLAS | Snakemake, strain tracking, three-point abundance |
| SqueezeMeta | functional abundance for reads, contigs and bins; R bindings |
| VEBA | viruses, pangenomics, BGC novelty, in one pipeline |
| Binette | faster, better-scaling refinement |
| MGnify / mettannotator | mobilome, versioned databases, validated on MAGs |

---

## 1. Changes that affect results

### 1.1 Multi-split binning

Status: done and measured, 2026-10-13 (§5, phase 2). It does not beat the shipped cohort path, so
`--binsplit` is a flag and not a cohort default. Open: the same experiment with the learned binners. No
status recorded for SemiBin2's native multi-sample mode.

MetaWrap2 bins one assembly at a time; in a cohort, each sample's binning uses a cross-sample depth
table. Multi-split bins the concatenated catalogue of all samples, then splits the bins per sample. Söding's analysis reports that splitting the
embedding space by sample gives 19% more bins and 79% more high-quality bins than splitting the final
clusters. nf-core/mag has this mode.

Results:

- When binning a catalogue, splitting the bin set helps. metaBAT2 goes from 86 bins and 16 clean
  high-quality MAGs to 166 and 37; MaxBin2 from 114 and 15 to 258 and 19, at the same cost (3.29
  against 3.28 CPU-hours; splitting is one pass over bins already written).
- Binning the catalogue does not improve on the shipped cohort path (same three assemblies, same
  cross-sample depth table, per-sample binning). Multi-split gives 47 genomes against 46 and 38
  high-quality bins against 37, but 29 clean high-quality MAGs against 34, for 9.09 CPU-hours against
  8.71. The final catalogue's GUNC pass fraction falls from 87.0% to 70.2%. One cluster over a
  concatenated catalogue can collect near-identical contigs from different samples, and splitting
  afterwards does not separate them.

Per-sample binning on a cross-sample depth table probably already captures the published gain. The
paper compares against single-sample coverage, which is our arm b, and that arm also gives 29 clean
MAGs. `--binsplit` is still needed when binning a catalogue, including with SemiBin2's and VAMB's
native multi-sample modes.

Caveats: one dataset, three samples, two classical binners. The paper's claim is about the embedding,
so the next experiment is the same comparison with the learned binners.

Implementation: `coverage.py` builds the prefixed concatenated catalogue. `binning.py` passes VAMB
`-o ''` to disable binsplitting for a single assembly, since there is no sample prefix and VAMB's
default separator `C` would split bins on any contig id containing a capital C. The planned work was:

- A cohort binning mode that bins the concatenated catalogue and splits bins back per sample, giving
  VAMB the prefix separator `coverage.py` already writes.
- The same for SemiBin2's native multi-sample mode.
- A benchmark against the existing cohort arms (a re-bin, not a re-assembly).

### 1.2 Bin quality errors

Status: open per this section; phase 1 in §5 has no recorded state. MIMAG tagging shipped in 2.1.3
([release notes](releases/2.1.3.md#4-mimag-quality-tags)).

- GUNC ≤1.0.6's `summarise` flagged every genome as passing; fixed in 1.1.1. `envs/binqc.yaml` asks for
  `gunc >=1.0.6`, so chimerism numbers from any arm scored with that version are invalid. Pin
  `>=1.1.1` and assert the version at runtime.
- CheckM2 ≤1.0.2 gives batch-dependent scores: the same genome scores differently depending on the rest
  of the batch. `envs/checkm2.yaml` asks for `>=1.0.2`. Pin above it and record the version per run.
- `scoring.QUALITY_TIERS` calls its tiers "MIMAG-style", and one is not. MIMAG high-quality is >90%
  complete, <5% contaminated, plus 23S, 16S and 5S rRNA and at least 18 tRNAs; medium is ≥50/<10, which
  `bins_medium_quality` matches. There is no 70/10 MIMAG tier, so `bins_high_quality` is a working
  threshold under a MIMAG name, and `bins_near_complete` (90/5) is MIMAG high-quality without the RNA
  requirement.

Fix the comment and the docs but keep the field name, which is the key in every recorded results
row; renaming it would need a migration.

Real MIMAG tagging would cost little: the rRNA and tRNA counts can be parsed from Bakta's tRNAscan-SE,
Aragorn and Infernal/Rfam output.

### 1.3 `dereplicate` and the skani unit bug (not affected)

Status: done (checked; no change needed).

dRep v4 fixed a skani alignment-fraction bug (percent read as fraction) that merged two distinct
*Enterococcus* species into one cluster, with no sign of the error in the output. `MIN_ALIGNED_FRACTION = 10.0` is passed to `skani triangle --min-af`, which
takes a percentage (skani's own default is 15), and the constant's docstring says "as a percentage".

### 1.4 CoverM identity thresholds

Status: open per this section; phase 1 in §5 has no recorded state.

`quant_bins.COVERM_GENOME` sets `--min-covered-fraction 10` and pins `--mapper`. It does not set
`--min-read-percent-identity` or `--min-read-aligned-percent`, whose CoverM defaults are 0, so every
alignment counts. This inflates MAG coverage with cross-mapping reads from related genomes (the 90%
identity floor in `depth.py` prevents this for the depth table). Set explicit floors, record them in
provenance, and measure their effect with the quantifier benchmark.

### 1.5 `classify_bins` defaults

Status: taxator-tk is deprecated with a per-run warning; the default is not flipped (§5, phase 4). The
GTDB-Tk version and reference pins are open per this section. §5 refers here for why flipping the
default needs a measurement first, but this section does not give that reason.

taxator-tk is abandoned and is the default; GTDB-Tk is opt-in. Make GTDB-Tk the default and retire
taxator-tk.

`envs/gtdbtk.yaml` asks for `gtdbtk >=2.4.0`, which resolves to 2.7.2. That version needs reference
data R232 and about 140 GB of RAM for `classify_wf`, and has removed flags our invocation may still
pass. Pin an exact version and the reference release, state the RAM requirement in `doctor`, and check
the flags against 2.7.x. Unpinned, it fails after install on any machine with less than 140 GB of RAM.

### 1.6 `reassemble_bins` coverage gating

Status: thinning done; skipping deep bins dropped from the plan.

Reassembly helps at low to moderate coverage and wastes compute or degrades the assembly at high
coverage. Options: thin each bin's reads to a target depth, or skip deeply covered bins.

Thinning is implemented: `--max-coverage` defaults to 50×, applied per bin against its own assembled
size, with a fixed seed so a re-run recruits the same reads.

Skipping is not planned. A bin at 200× thinned to 50× reassembles well, and thinning already saves
the compute.

---

## 2. Missing domains

Self-contained modules, ordered by how often users need them.

### 2.1 Viruses

Status: done, including CRISPR host linking (`--crispr-hosts`, cctyper; §5, phase 5).

MetaWrap2 bins prokaryotes and discards phage, usually the most abundant entity in these datasets.

geNomad is the caller. It classifies contigs as viral, plasmid or chromosomal in one pass with one
database, so `amr`'s mobility work (annotation plan §6.4) uses the same download. On one CAMI2
assembly, the same run found 161 viral contigs at ≥2 kb and 344 plasmid contigs.

Screening the bins of that sample finds zero viral contigs; screening the assembly finds 161. Phage
contigs are usually unbinned, so free phage is called from the assembly and bins are used for
proviruses. VIBRANT and VirSorter2 are less current than geNomad.

CheckV estimates completeness and contamination and has no real alternative. It also trims host
regions from provirus contigs, so an integrated prophage is reported as a virus plus flanking
bacterial chromosome instead of one chimeric contig.

Module outputs, analogous to `binning` and `bin_refinement`:

- viral contigs with geNomad's taxonomy and CheckV's completeness/contamination per contig;
- quality tiers following the MIUViG standard (high-quality ≥90% complete, genome-fragment below);
- vOTU clustering at the conventional 95% ANI / 85% AF, so counts represent populations (`dereplicate`
  does the equivalent for genomes);
- abundance per vOTU per sample through the existing `coverage`/`quant_bins` path;
- host prediction: `cctyper` (§2.4) finds CRISPR spacers in the prokaryotic bins, and a spacer matching
  a vOTU is direct evidence of its host.

Caveat to print with the results: CheckV completeness is estimated against a reference database, so a
novel virus with no relative is reported as a fragment whether or not it is one.

### 2.2 Plasmids and mobile elements

Status: contig classification done (geNomad, 344 plasmid contigs on one CAMI2 assembly, feeding
`amr --mobile-elements`). PlasMAAG open. Blockers: it is not on bioconda (it ships inside the VAMB
repository), and it needs the assembly graph, which `assembly --keep-graph` now writes to a
named path.

Covered in the annotation plan §6.4: geNomad for plasmid contigs (the same run `viruses` needs),
PlasmidFinder for replicon typing, ISEScan for insertion sequences, mobileOG-db for MGE families.

PlasMAAG (*Nat Biotechnol* 2026) does plasmid-aware assembly-graph binning, which improves
substantially on classifying contigs. The graph part is done: `assembly --keep-graph` keeps metaSPAdes'
`assembly_graph_with_scaffolds.gfa` (121 MB measured, hard-linked, so no extra disk) and builds
MEGAHIT's FASTG from its intermediates with `megahit_toolkit contig2fastg` (235 MB, two seconds on one
core). What remains is wiring up the binner. Off by default. Correction to earlier notes: `assembly`
never deleted a graph; metaSPAdes and metaFlye write one and it was kept, and MEGAHIT writes none. The
geNomad path is the baseline for the comparison.

### 2.3 Strain resolution

Status: done per §5, phase 7; the original text of this section said "we have nothing".

inStrain 1.11.0 for within-sample strain variation, or TRACS, the only tool evaluated that does not
overestimate SNP distances, which is required for any transmission claim. Both are maintained. ATLAS
has this.

### 2.4 CRISPR, pangenomics, phylogenomics

Status: cctyper done (`viruses --crispr-hosts`). Pangenomics and phylogenomics done per §5, phase 7.

- cctyper for CRISPR-Cas arrays and phage-host links through spacer matching. On five CAMI2 bins, three
  carry arrays (I-B, two I-C) with 54 spacers, and none matches the sample's 161 viral contigs. All 54
  spacers are recovered from their own genomes, so the matcher works. Expect few or no matches for a
  single sample. A match is evidence of a host link.
- PPanGGOLiN for pangenome graphs across `dereplicate`'s output.
- IQ-TREE for phylogenomics from the marker genes GTDB-Tk already extracts.

---

## 3. Unbinned fraction

Status: binned fraction and chimerism caveat done (`metawrap2.recovery`); `singlem appraise` open (§5,
phase 3).

MAGs are chimeric more often than their quality scores suggest. Treichel et al. 2026: among
high-quality MAGs, only 54.5–81.8% represent a single strain, with a mean of 2.3 reference genomes per
high-quality MAG. CheckM2 does not detect this; GUNC detects part of it, at the taxonomic level.

Most of the assembly never ends up in a MAG. SPIRE, across 92,187 metagenomes: only about 10% of
assembled base pairs end up in medium- or high-quality MAGs, and the unbinned remainder contains an
estimated 145 predicted undescribed phyla. A MAG-only report describes a small, non-representative
part of the sample.

To add to the binning and refinement report:

- the fraction of assembled bp binned, and the fraction binned into clean 70/10 bins;
- `singlem appraise` output: what fraction of the community the bins represent per cell, and which
  lineages have no genome;
- a chimerism caveat next to every quality tier, from the existing GUNC result.

---

## 4. Harness and statistics

### 4.1 Long reads

Status: done per §5, phase 8 (metaFlye, long-read mapping presets, SemiBin2's long-read mode). Hybrid
assembly and polishing have no recorded status.

Scope: long-read assembly (metaFlye), hybrid assembly, long-read-aware binning, polishing (nf-core/mag
supports long reads). The largest single item in either plan, to be scoped separately.

### 4.2 Statistics and interchange

Status: ANCOM-BC and BIOM done (`diff_abundance`, `--biom`; §5, phase 6). No status recorded for
MaAsLin 3.

The pipeline stops at feature tables, so users write their own downstream analysis.

- MaAsLin 3 for differential abundance, with its prevalence/abundance decomposition.
- ANCOM-BC as the default: a pure-Python compositional test in scikit-bio, no R dependency.
- BIOM 2.1 output for every feature table (one implementation, shared with the annotation plan).
- Treat relative abundances as compositional throughout; sum-constrained data must not go into a raw
  t-test.

### 4.3 Diversity coverage

Status: done, as `metawrap2 nonpareil`.

A separate module, because it needs one read file per sample across the cohort, which fits neither
`read_qc` nor `coverage`. Implementation notes:

- The two kernels disagree on short reads: 50.5% against 83.0% coverage on the same 34 bp gut sample.
  The kernel is chosen from the measured read lengths, once for the whole cohort, so coverage figures
  are comparable across samples.
- The k-mer kernel's effort axis is in k-mers, not base pairs (`AL` = L-k+1). Reporting it as bp
  overstates the sequencing effort by L/AL, 2.3x on those reads.
- The alignment kernel seeds from the clock: two identical runs gave 83.0% and 84.3%.

The model fit is ported from the R package and matches it to six significant figures on Nonpareil's
three example datasets. Diversity coverage, more than sequencing depth, confounds cross-sample
comparison, and most published comparisons ignore it. Runtime: 4 minutes per sample with the k-mer
kernel, 16 with the alignment kernel.

### 4.4 Benchmark harness

Status: multi-split arms done (§1.1). Module completeness against MAG completeness done for one genome
(annotation plan §5.3). No status recorded for database version pinning.

- Multi-split binning arms (§1.1).
- Module completeness against MAG completeness (annotation plan §5.3).
- Pin and record every database version in `benchmarks/configurations.toml`, as well as tool versions.
  GTDB R220 against R232, KOfam snapshots and AMRFinderPlus database dates all change results.

---

## 5. Order of work

| Phase | Work | Why here |
|---|---|---|
| 1 | Correctness: GUNC pin, CheckM2 pin, MIMAG labelling, CoverM floors, GTDB-Tk pin and default | Some reported numbers are currently wrong. Cheap, no new dependencies |
| 2 | Multi-split binning plus benchmark arms | **done and measured**: it does not beat the shipped cohort path - 29 clean HQ MAGs against 34. `--binsplit` stays a flag; see §1.1 |
| 3 | Binned fraction, `singlem appraise`, chimerism caveat | binned fraction and chimerism caveat **done** (`metawrap2.recovery`); `singlem appraise` open |
| 4 | Retire taxator-tk | **deprecated with a per-run warning**, default not flipped - see §1.5 for why that needs a measurement first |
| 5 | `viruses`: geNomad + CheckV + vOTU clustering + CRISPR host linking | **done**, CRISPR host linking included (`--crispr-hosts`, cctyper) |
| 6 | Statistics, BIOM, Nonpareil | **done**: ANCOM-BC and BIOM (`diff_abundance`, `--biom`), Nonpareil (`metawrap2 nonpareil`) |
| 7 | Strain resolution (inStrain/TRACS), PlasMAAG, pangenomics, phylogenomics | **done** except PlasMAAG itself; the graph it needs is now kept by `assembly --keep-graph` - see §2.2 |
| 8 | Long reads | **done**: metaFlye, long-read mapping presets, SemiBin2's long-read mode; verified on simulated R10 reads - 72 contigs at N50 57.9 kb from a 2.50 Mbp reference |

Phase 1 is about a week of work and fixes numbers that are wrong now, so it goes first.
