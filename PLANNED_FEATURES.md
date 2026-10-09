# Planned features

Ideas for modernizing MetaWrap2, based on the tools and methods the metagenomics community is
currently using. A living wishlist, not a commitment — an item graduates into the pipeline when
it is implemented behind a flag and verified.

Sections are organized **by area**. Every item carries a complexity tag:

- `[low]` — drop-in tool swap, a new optional flag, or an extra output. Small effort.
- `[moderate]` — a new module or alternative method; real work, but fits the current structure.
- `[major]` — architecture, new data types, ML, or scaling. Multi-release.

A complexity-ordered [roadmap](#suggested-roadmap) is at the end, and a [references](#references)
list at the very bottom. Items that have since graduated are listed under
[Already shipped](#already-shipped) and removed from their sections, so what is left here is what is
genuinely not built.

> **On citations.** References were compiled from established knowledge of these tools, not a live
> literature pull, and give author/journal/year rather than DOIs. Verify the exact citation,
> version, and current best practice before using this document formally — tool versions and
> benchmarks move fast (e.g. GTDB and CheckM2 databases, HUMAnN reference builds).

## Already shipped

Graduated out of this list; see [docs/benchmarks.md](docs/benchmarks.md) for what each was measured
to be worth.

| | What it does | How to ask for it |
|---|---|---|
| fastp | read trimming, single pass, threads past cutadapt's ceiling | default; `--trimmer trim-galore` for the old one |
| SemiBin2 | self-supervised learned binner | default (`--preset default`), or `--semibin2` |
| VAMB | variational-autoencoder binner | `--vamb` |
| COMEBin | contrastive learned binner, the best catalogue measured | `--comebin`, or `--preset slow` |
| GPU training for all three | one switch instead of per-tool flags | `--device {cpu,cuda,auto}` |
| DAS Tool | single-pass consolidation, an alternative to hybridising | `--consolidator dastool` |
| CheckM2 | neural-net completeness, no pplacer memory wall | `--checkm2` |
| GUNC | chimerism the marker-gene checks cannot see | on by default; `--skip-gunc` |
| GTDB-Tk | standardised bin taxonomy | `--gtdbtk` |
| Bakta | modern annotation | `--bakta` |
| CoverM | per-genome relative abundance | `--quantifier coverm` |
| sylph | alignment-free bin abundance from k-mer sketches | default |
| skani | ANI dereplication across a cohort | `metawrap2 dereplicate` |
| fairy | alignment-free coverage estimation | `--coverage fairy` |
| Cross-sample coverage | one depth table, N alignments rather than N² | `metawrap2 coverage` |
| Per-sample assemble, bin, then dereplicate | the strategy most catalogues are built with | `metawrap2 cohort` |
| A benchmark harness | comparisons declared as data, results committed | `metawrap2 benchmark` |

A note on dRep, which this list asked for: `metawrap2 dereplicate` implements the same idea - a
fast pre-filter, then ANI within clusters, then one representative per cluster by quality - with
skani doing both stages rather than Mash then gANI. dRep itself is not wrapped.

## Guiding principles

- **Opt-in until measured.** New tools arrive behind a flag. A default only moves on a
  measurement committed to `benchmarks/results/`, and the move is stated in the changelog and in
  the migration notes - never silently.
- **Flat and legible.** New capability is a new module or a command block in an existing module —
  not a framework to learn.
- **Reproducible.** Each tool gets a pinned conda env with a committed lockfile, a database
  fingerprint, and a line in the run provenance.
- **Standards-aware.** Prefer outputs and quality tags matching community standards — MIMAG MAG
  quality (Bowers et al., *Nat Biotechnol* 2017), GTDB taxonomy (Parks et al., *Nat Biotechnol*
  2018/2020), MIxS metadata (Yilmaz et al., *Nat Biotechnol* 2011).
- **Benchmarked, not assumed.** Declare each new method as an arm in
  `benchmarks/configurations.toml` and run it against a known-truth dataset (Sczyrba et al.,
  *Nat Methods* 2017; Meyer et al., *Nat Methods* 2022). Quality and cost are reported side by
  side, never divided into each other.

---

## Read QC and preprocessing

- `[low]` **Modern host removal** — a bowtie2/minimap2 scrub against T2T-CHM13 (Nurk et al.,
  *Science* 2022), which removes far more human read contamination than hg38; `hostile` wraps this
  with masking of microbial-homologous regions. *(Shen et al./Constantinides et al., Bioinformatics
  2023.)* An alternative to bmtagger.
- `[low]` **Long-read QC** — NanoPlot/NanoPack for ONT/PacBio read stats, chopper/filtlong for
  length/quality filtering, so QC is not Illumina-only. *(De Coster et al., Bioinformatics
  2018/2023.)*
- `[low]` **MultiQC** — aggregate per-sample FastQC/fastp reports into one study report. *(Ewels
  et al., Bioinformatics 2016.)*

## Assembly

Current defaults are metaSPAdes (Nurk et al., *Genome Res* 2017) and MEGAHIT (Li et al.,
*Bioinformatics* 2015); both stay.

- `[moderate]` **Long-read assembly** — metaFlye for ONT/PacBio (Kolmogorov et al., *Nat Methods*
  2020), and for PacBio HiFi metaMDBG (Benoit et al., *Nat Biotechnol* 2024) and hifiasm-meta (Feng
  et al., *Nat Methods* 2022). HiFi routinely yields complete, circular MAGs.
- `[moderate]` **Hybrid assembly** — OPERA-MS (Bertrand et al., *Nat Biotechnol* 2019) or
  Unicycler (Wick et al., *PLoS Comput Biol* 2017) combining short + long reads.
- `[moderate]` **Long-read polishing** — racon (Vaser et al., *Genome Res* 2017) and medaka (ONT)
  consensus after assembly.

## Binning

`bin_refinement` already consolidates any set of bins, so new binners slot in as `binning`
backends and feed refinement. Current binners: metaBAT2 (Kang et al., *PeerJ* 2019), MaxBin2 (Wu
et al., *Bioinformatics* 2016), CONCOCT (Alneberg et al., *Nat Methods* 2014), SemiBin2 (Pan et
al., *Bioinformatics* 2023), VAMB (Nissen et al., *Nat Biotechnol* 2021) and COMEBin (Wang et al.,
*Nat Commun* 2024).

- `[moderate]` **GraphMB** (Lamurias et al., *Bioinformatics* 2022, assembly-graph aware),
  **MetaDecoder** (Liu et al., *Microbiome* 2022), **binny** (Hickl et al., *NAR Genom Bioinform*
  2022) — recent binners with competitive CAMI/benchmark results, none started.
- `[moderate]` **Hi-C binning** — bin3C (DeMaere & Darling, *PeerJ* 2019) or MetaTOR (Baudry et
  al., *Front Genet* 2019), using contact frequency to resolve strains and place plasmids/mobile
  elements in the right MAG.

## Bin QC and quality standards

- `[low]` **MIMAG quality tags** — classify each MAG high/medium/low quality in `.stats.tsv` from
  completeness, contamination, and rRNA/tRNA presence (Bowers et al., *Nat Biotechnol* 2017).
- `[low]` **tRNA/rRNA detection** — tRNAscan-SE (Chan et al., *NAR* 2021) and barrnap (Seemann),
  feeding the MIMAG call and annotation.

## Taxonomic profiling and abundance

- `[low]` **Bracken** — Bayesian re-estimation of Kraken2 abundances at a chosen rank. *(Lu et al.,
  PeerJ CS 2017; Kraken2: Wood et al., Genome Biol 2019.)*
- `[moderate]` **Read-level profilers** — assembly-free community profiles: sylph (fast
  ANI/containment; Shaw & Yu, *Nat Biotechnol* 2024), MetaPhlAn 4 (Blanco-Míguez et al., *Nat
  Biotechnol* 2023), mOTUs (Ruscheweyh et al., *Microbiome* 2022), Ganon (Piro et al.,
  *Bioinformatics* 2020).

## Functional and pathway annotation

Two questions, both with standard abundance tables and plots: *what can each MAG do*, and *what is
the whole community doing and how much*. Baseline annotation is Prokka (Seemann, *Bioinformatics*
2014) or Bakta (Schwengers et al., *Microb Genom* 2021).

**Per-MAG pathway and metabolism reconstruction** `[moderate]` — genes → orthologs → pathway/module
completeness → metabolic summary.
- Ortholog assignment: eggNOG-mapper (COG/KEGG-KO/GO/EC; Cantalapiedra et al., *MBE* 2021),
  KOfamScan/KofamKOALA (KEGG KO via HMM thresholds; Aramaki et al., *Bioinformatics* 2020),
  InterProScan (Pfam/InterPro/GO; Jones et al., *Bioinformatics* 2014).
- Completeness & summaries: **DRAM** (distills KEGG/Pfam/dbCAN/MEROPS into a metabolism summary +
  product heatmap; Shaffer et al., *NAR* 2020); **METABOLIC** (biogeochemical trait profiles +
  cycling diagrams; Zhou et al., *Microbiome* 2022); **anvi'o `anvi-estimate-metabolism`** (KEGG
  module completeness; Eren et al., *Nat Microbiol* 2021); **KEGGDecoder** (module-completeness
  heatmap; Graham et al., *ISME J* 2018); **MinPath** (parsimonious pathway inference; Ye & Doak,
  *PLoS Comput Biol* 2009).
- Specialized: dbCAN (CAZymes; Zheng et al., *NAR* 2023), antiSMASH (biosynthetic gene clusters;
  Blin et al., *NAR* 2023), MEROPS (peptidases; Rawlings et al.), AMRFinderPlus (Feldgarden et
  al., *Sci Rep* 2021) / CARD-RGI (Alcock et al., *NAR* 2023).
- Output: a MAG × KEGG-module (or MetaCyc pathway) completeness matrix + heatmap.

**Community pathway and functional abundance profiling** `[moderate]` — quantify gene-family and
pathway abundance across a whole sample, stratified by contributing taxa.
- **HUMAnN 3** (standard): reads → UniRef gene families and MetaCyc pathways; emits
  `genefamilies`, `pathabundance`, and `pathcoverage` tables stratified by species, with utilities
  to normalize (CPM/relab), join samples, regroup to KO/EC/GO, and barplot. *(Beghini et al.,
  eLife 2021.)*
- **Woltka** — functional + taxonomic profiling from alignments to KO/pathway/GO against the
  WoL/GTDB reference; feature tables (TSV/BIOM). *(Zhu et al., mSystems 2022.)*
- **SqueezeMeta + SQMtools** — all-in-one producing KEGG/COG/PFAM abundance tables for reads and
  bins, with R-side tables and plots. *(Tamames & Puente-Sánchez, Front Microbiol 2019;
  Puente-Sánchez et al., BMC Bioinformatics 2020.)*
- **DIAMOND + MEGAN-LR** — fast protein alignment (Buchfink et al., *Nat Methods* 2021) → functional
  binning to KEGG/SEED/eggNOG and cross-sample comparison (Huson et al., *PLoS Comput Biol* 2016).
- Output: sample × pathway and sample × KO/gene-family abundance matrices, normalized.

**Genome-scale metabolic models** `[major]` — mechanistic pathway and cross-feeding analysis:
gapseq (Zimmermann et al., *Genome Biol* 2021), CarveMe (Machado et al., *NAR* 2018), or ModelSEED
per MAG; metaGEM end-to-end from metagenomes to community models (Zorrilla et al., *NAR* 2021).
Output: SBML models and predicted metabolic interactions.

**Tables, plots, and statistics** `[low–moderate]`
- Matrices (TSV/BIOM, joinable): MAG × module completeness; sample × pathway; sample × KO.
- Plots: module-completeness heatmap across MAGs; stratified pathway-abundance barplots; ordination.
- Normalization: relative abundance, CPM, per-genome-equivalent (MicrobeCensus; Nayfach &
  Pollard, *Genome Biol* 2015).
- Differential abundance: hand tables to MaAsLin2 (Mallick et al., *PLoS Comput Biol* 2021), ALDEx2
  (Fernandes et al., *Microbiome* 2014), or LEfSe (Segata et al., *Genome Biol* 2011).

## Strain-level analysis

- `[moderate]` **inStrain** — within-population microdiversity, popANI, and strain comparison across
  samples from read mappings. *(Olm et al., Nat Biotechnol 2021.)*
- `[moderate]` **Strainberry / StrainGE** — strain-aware assembly and tracking. *(Vicedomini et
  al., Nat Commun 2021; van Dijk et al., Genome Biol 2022.)*

## Viruses, plasmids, and mobile elements

- `[moderate]` **geNomad** — unified virus + plasmid identification with gene-based classification.
  *(Camargo et al., Nat Biotechnol 2024.)*
- `[moderate]` **CheckV** — viral genome quality/completeness, the MIMAG-analogue for viruses.
  *(Nayfach et al., Nat Biotechnol 2021.)*
- `[moderate]` **VirSorter2** — viral contig detection across dsDNA/ssDNA/RNA groups. *(Guo et al.,
  Microbiome 2021.)*
- `[moderate]` **MOB-suite / plasmidVerify** — plasmid reconstruction and typing. *(Robertson &
  Nash, Microb Genom 2018.)*

## Cohort and multi-sample analysis

How per-sample results are combined largely determines both compute cost and what genomes you
recover. MetaWrap2 should support both strategies well and make the cohort steps first-class.

**Both aggregation strategies now exist** — `metawrap2 cohort` assembles and bins per sample and
then dereplicates, and `metawrap2 run` co-assembles. Which is better is still open: measured on
three samples they tie on clean high-quality MAGs and per-sample wins at the stricter tier, and the
argument that favours co-assembly is about scaling, which is untested. The experiment that would
settle it is twenty samples, where co-assembly's assembler is expected to fail on memory rather
than merely cost more.

- `[low]` **Report which strategy a catalogue came from** in the cohort summary, since the two are
  not comparable and a catalogue outlives the command that made it.

**Cohort-level steps on the dereplicated MAG set:**
- `[low]` **MAG × sample abundance matrix** — map every sample's reads back to the non-redundant
  set (CoverM, or inStrain for strain-resolved), with a breadth-of-coverage detection threshold so
  low-coverage noise is not called present.
- `[low]` **Prevalence / presence-absence** across the cohort from that matrix.
- `[moderate]` **Differential abundance / association testing** between groups on MAG, taxon,
  pathway, and gene-family tables — MaAsLin2, ALDEx2, LEfSe.
- `[moderate]` **Community structure** — beta-diversity/ordination of taxonomic and functional
  profiles, with sample metadata (MIxS) integrated.
- `[moderate]` **Strain tracking** across samples/subjects via inStrain popANI (transmission,
  engraftment, temporal dynamics).
- `[major]` **Cohort mode** — one `metawrap2 run`-level command that orchestrates per-sample MAG
  recovery → dRep → back-mapping → the abundance/prevalence tables and plots above, resumable and
  provenanced like everything else.

## Scaling and execution

- `[major]` **Workflow-engine backend** — run the flat modules under Nextflow (Di Tommaso et al.,
  *Nat Biotechnol* 2017) or Snakemake (Mölder et al., *F1000* 2021) for HPC/cloud scheduling,
  retries, and caching; interop with nf-core/mag conventions (Krakau et al., *NAR Genom Bioinform*
  2022; Ewels et al., *Nat Biotechnol* 2020).
- `[major]` **Cloud-native execution** — Slurm / cloud batch / Kubernetes executors, object-store
  (S3) I/O, per-step containers (Docker/Apptainer already exist).
- `[major]` **GPU acceleration** — for DL binners/annotators and alignment.

## Machine learning and modern methods

- `[major]` **Protein language models** — ESM-2/ESMFold (Lin et al., *Science* 2023) or ProtT5 for
  functional inference on hypothetical proteins where homology annotation fails.
- `[major]` **Genomic/DNA language models** — for binning, contig classification, and novelty
  detection (an active research area; benchmark before adopting).
- `[major]` **Learned bin refinement** — improve on the current heuristic consolidation.

## New data types and analyses

- `[major]` **Long-read-first track** — a modern default path (HiFi/ONT) aimed at complete/circular
  MAGs.
- `[major]` **Metatranscriptomics** — expression against recovered MAGs; links to metaproteomics.
- `[major]` **Pangenome / comparative genomics** — core/accessory and gene gain/loss across MAGs
  (e.g. Panaroo; Tonkin-Hill et al., *Genome Biol* 2020).
- `[major]` **Real-time / streaming** nanopore analysis.

## Reporting, standards, and interoperability

- `[low]` **Per-MAG summary table** — quality, taxonomy, size, N50, coverage, tRNA/rRNA in one
  table.
- `[moderate]` **Interactive report / dashboard** — an anvi'o-style or rich HTML MAG/function
  browser (Eren et al., *Nat Microbiol* 2021).
- `[moderate]` **Standards & submission** — MIxS/MIMAG metadata capture, NCBI/ENA submission
  helpers, and workflow provenance in a standard form (RO-Crate) building on existing provenance.
- `[major]` **Benchmarking harness** — CAMI2 gold-standard datasets (Meyer et al., *Nat Methods*
  2022) tracking accuracy across releases.

## Databases

Installable and version-pinned via `install-db`, with checksums and a recorded fingerprint:

- Taxonomy/QC: GTDB-Tk reference, CheckM2 DB, Kraken2/Bracken, CheckV DB, geNomad DB.
- Function: KOfam HMMs, eggNOG, dbCAN, CARD, antiSMASH, MEROPS, InterPro; HUMAnN's ChocoPhlAn and
  UniRef90; MetaCyc.
- Licensing note: the full KEGG database needs a subscription — default to freely usable options
  (KOfam, eggNOG, MetaCyc via HUMAnN) and let licensed users point at KEGG.

---

## Suggested roadmap

Ordered by complexity and value:

1. **`[low]` the remaining small wirings** — MIMAG quality tags, tRNA/rRNA detection, Bracken,
   and the per-MAG summary table. Steps 1 and 2 of this list are otherwise done; see
   [Already shipped](#already-shipped).
2. **`[moderate]` functional and pathway annotation** — per-MAG metabolism (DRAM, KEGG-module
   completeness) and community pathway abundance (HUMAnN 3) with tables and plots.
3. **`[moderate]` long-read / hybrid assembly** and a long-read track.
4. **`[moderate]` viruses/plasmids and strain-level** — geNomad, CheckV, inStrain; Hi-C binning.
5. **`[major]` scaling and modern methods** — workflow-engine backend, cloud execution, ML
   annotation, dashboards, and standards.

---

## References

Compiled from established knowledge; verify exact citation, version, and DOI before formal use.
Reviews and benchmarks worth reading alongside the tool papers:

- Quince, Walker, Simpson, Loman, Segata. *Shotgun metagenomics, from sampling to analysis.* Nat
  Biotechnol, 2017.
- Sczyrba et al. *Critical Assessment of Metagenome Interpretation (CAMI).* Nat Methods, 2017; and
  Meyer et al. *CAMI II.* Nat Methods, 2022 — the standard benchmarking datasets/challenges.
- Bowers et al. *Minimum information about a metagenome-assembled genome (MIMAG).* Nat Biotechnol,
  2017 — the MAG quality standard.
- Parks et al. *A standardized bacterial taxonomy (GTDB).* Nat Biotechnol, 2018 (and subsequent
  GTDB releases) — the taxonomy backbone.
- Nayfach et al. *A genomic catalog of Earth's microbiomes.* Nat Biotechnol, 2021 — scale and QC
  practice for large MAG catalogs.

Individual tool citations are given inline next to each tool above.
