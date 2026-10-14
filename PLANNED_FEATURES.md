# Planned features

What is not built yet. Items leave this list when they ship; what each release added is in
[docs/releases/](docs/releases/). Two detailed plans cover parts of this list with versions,
licences and citations checked:
[functional annotation, pathways, AMR and abundance](docs/plan_functional_annotation.md) and
[the other modules](docs/plan_other_modules.md). Where they disagree with this page, follow the
plans.

Each item has a size tag:

- `[low]`: a tool swap, a new flag, or an extra output.
- `[moderate]`: a new module or an alternative method within the current structure.
- `[major]`: architecture, new data types, or scaling; more than one release.

> Citations give author, journal and year from memory rather than a literature search. Check the
> exact reference and current version before citing this page.

## Principles

- New tools arrive behind a flag. A default changes only on a measurement committed to
  `benchmarks/results/`, and the change is stated in the release notes and the migration table.
- New capability is a new module, or a command block in an existing one.
- Every tool gets a pinned conda environment with a committed lockfile, a database fingerprint, and
  a line in the run provenance.
- Prefer community standards: MIMAG quality (Bowers et al., *Nat Biotechnol* 2017), GTDB taxonomy
  (Parks et al., *Nat Biotechnol* 2018/2020), MIxS metadata (Yilmaz et al., *Nat Biotechnol*
  2011).
- Every new method is declared as an arm in `benchmarks/configurations.toml` and measured against a
  known-truth dataset (Sczyrba et al., *Nat Methods* 2017; Meyer et al., *Nat Methods* 2022).
  Quality and cost are reported side by side.

## Open from 2.1.2

These are started or decided, with a named blocker.

- `[low]` **MIMAG quality tags** in `.stats.tsv`, from completeness, contamination and rRNA/tRNA
  presence. Nearly free: `annotate --bakta` already runs tRNAscan-SE, Aragorn and Infernal, so the
  counts are a parse of existing output. Halted on the maintainer's instruction; resumes only on
  request.
- `[low]` **GTDB-Tk as the `classify_bins` default**, retiring taxator-tk. taxator-tk is deprecated
  with a warning; switching the default needs a measurement first.
- `[moderate]` **PlasMAAG**, plasmid-aware binning on the assembly graph (*Nat Biotechnol* 2026).
  Blocked: `assembly` deletes the graph, and PlasMAAG is not on bioconda.
- `[moderate]` **VFDB virulence calls** in `amr`. Needs a measured comparison arm before it ships.
- `[moderate]` **InterProScan member databases** beyond Pfam-A, for genes KEGG cannot place.

## Read QC

- `[low]` **Host removal against T2T-CHM13** (Nurk et al., *Science* 2022), which removes more human
  reads than hg38. `hostile` does this with masking of microbial-homologous regions (Constantinides
  et al., *Bioinformatics* 2023). An alternative to bmtagger.
- `[low]` **Long-read QC**: NanoPlot for read statistics, chopper or filtlong for filtering (De
  Coster et al., *Bioinformatics* 2018/2023).
- `[low]` **MultiQC** to combine per-sample QC reports into one (Ewels et al., *Bioinformatics*
  2016).

## Assembly

- `[moderate]` **HiFi assemblers**: metaMDBG (Benoit et al., *Nat Biotechnol* 2024) and hifiasm-meta
  (Feng et al., *Nat Methods* 2022). HiFi often gives complete, circular MAGs. metaFlye is already
  supported.
- `[moderate]` **Hybrid assembly** of short and long reads: OPERA-MS (Bertrand et al., *Nat
  Biotechnol* 2019) or Unicycler (Wick et al., *PLoS Comput Biol* 2017).
- `[moderate]` **Long-read polishing** with racon (Vaser et al., *Genome Res* 2017) and medaka.

## Binning

- `[moderate]` **More binners**: GraphMB (Lamurias et al., *Bioinformatics* 2022, uses the assembly
  graph), MetaDecoder (Liu et al., *Microbiome* 2022), binny (Hickl et al., *NAR Genom Bioinform*
  2022).
- `[moderate]` **Hi-C binning**: bin3C (DeMaere & Darling, *PeerJ* 2019) or MetaTOR (Baudry et al.,
  *Front Genet* 2019), which uses contact frequency to separate strains and place plasmids.

## Taxonomic profiling

- `[low]` **Bracken** to re-estimate Kraken2 abundances at a chosen rank (Lu et al., *PeerJ CS*
  2017).
- `[moderate]` **Read-level profilers**: MetaPhlAn 4 (Blanco-Míguez et al., *Nat Biotechnol* 2023),
  mOTUs (Ruscheweyh et al., *Microbiome* 2022), ganon (Piro et al., *Bioinformatics* 2020).

## Function

- `[moderate]` **eggNOG-mapper** for COG, GO and EC terms alongside KOfam (Cantalapiedra et al.,
  *MBE* 2021).
- `[moderate]` **Read-based pathway profiling** with HUMAnN 3 (Beghini et al., *eLife* 2021):
  UniRef gene families and MetaCyc pathways, stratified by species.
- `[moderate]` **Peptidases (MEROPS)** and **antiSMASH** as an optional deeper BGC caller next to
  GECCO.
- `[major]` **Genome-scale metabolic models**: gapseq (Zimmermann et al., *Genome Biol* 2021),
  CarveMe (Machado et al., *NAR* 2018), or metaGEM end to end (Zorrilla et al., *NAR* 2021).
- `[low]` **Per-genome-equivalent normalisation** with MicrobeCensus (Nayfach & Pollard, *Genome
  Biol* 2015).
- `[moderate]` **More differential-abundance tests** next to ANCOM-BC: MaAsLin2 (Mallick et al.,
  *PLoS Comput Biol* 2021) and ALDEx2 (Fernandes et al., *Microbiome* 2014).

## Strains, viruses and plasmids

- `[moderate]` **Strain-aware assembly and tracking**: Strainberry (Vicedomini et al., *Nat Commun*
  2021) and StrainGE (van Dijk et al., *Genome Biol* 2022).
- `[moderate]` **Plasmid reconstruction and typing** with MOB-suite (Robertson & Nash, *Microb
  Genom* 2018).

## Cohorts

- `[low]` **Record which strategy built a catalogue** (per-sample or co-assembly) in the cohort
  summary, since the two are not comparable.
- `[moderate]` **Settle per-sample versus co-assembly at scale.** On three samples they tie on clean
  MAGs. The deciding experiment is about twenty samples, where co-assembly is expected to run out
  of memory.
- `[moderate]` **Community structure**: beta diversity and ordination of taxonomic and functional
  profiles, with MIxS sample metadata.

## Execution

- `[major]` **Workflow-engine backend**: run the modules under Nextflow (Di Tommaso et al., *Nat
  Biotechnol* 2017) or Snakemake (Mölder et al., *F1000Research* 2021), following nf-core/mag
  conventions (Krakau et al., *NAR Genom Bioinform* 2022).
- `[major]` **Cloud execution**: Slurm, cloud batch or Kubernetes executors and object-store input/output. Docker
  and Apptainer images already exist.

## Methods research

- `[major]` **Protein language models** (ESM-2; Lin et al., *Science* 2023) for proteins homology
  cannot annotate.
- `[major]` **DNA language models** for binning and contig classification. Benchmark before
  adopting.
- `[major]` **Learned bin refinement** to improve on the current consolidation.

## New data types

- `[major]` **A long-read-first default path** (HiFi or ONT) aimed at complete, circular MAGs.
- `[major]` **Metatranscriptomics** against recovered MAGs.
- `[major]` **Real-time nanopore analysis.**

## Reporting and standards

- `[low]` **One per-MAG summary table**: quality, taxonomy, size, N50, coverage, tRNA and rRNA.
- `[moderate]` **An interactive HTML report** for browsing MAGs and their functions.
- `[moderate]` **Submission support**: MIxS/MIMAG metadata, NCBI/ENA submission helpers, and
  provenance as RO-Crate.

## Databases

Still to add to `install-db`: eggNOG, MEROPS, antiSMASH, the remaining InterPro member databases,
and HUMAnN's ChocoPhlAn and UniRef90.

## Background reading

- Quince, Walker, Simpson, Loman, Segata. *Shotgun metagenomics, from sampling to analysis.* Nat
  Biotechnol, 2017.
- Sczyrba et al. *CAMI.* Nat Methods, 2017; Meyer et al. *CAMI II.* Nat Methods, 2022.
- Bowers et al. *MIMAG.* Nat Biotechnol, 2017.
- Parks et al. *GTDB.* Nat Biotechnol, 2018.
- Nayfach et al. *A genomic catalog of Earth's microbiomes.* Nat Biotechnol, 2021.
