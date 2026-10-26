# Planned features

Features not yet built. Items are removed when they ship; what each release added is in
[releases/](releases/).

Each item has a size tag:

- `[low]`: a tool swap, a new flag, or an extra output.
- `[moderate]`: a new module or an alternative method within the current structure.
- `[major]`: architecture, new data types, or scaling; more than one release.

> Citations give author, journal and year from memory, without a literature search. Check the
> exact reference and current version before citing this page.

## Principles

- New tools are added behind a flag. A default changes only after a recorded benchmark
  measurement, and the change is listed in the release notes and the migration table.
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

Started or decided, each with a named blocker.

- `[low]` **GTDB-Tk as the `classify_bins` default**, replacing taxator-tk, which already prints a
  deprecation warning. Needs a measurement before the switch.
- `[moderate]` **PlasMAAG**, plasmid-aware binning on the assembly graph (*Nat Biotechnol* 2026).
  Blocked: PlasMAAG is not on bioconda. `assembly --keep-graph` (2.1.3) keeps the graph at a named
  path. On MEGAHIT the kept graph joins to no contig and 96.9% of its nodes are isolated, so a
  graph-aware tool gets no adjacency from it; see [Kept assembly graphs](#kept-assembly-graphs).
- `[moderate]` **VFDB virulence calls** in `amr`. Needs a measured comparison arm before it ships.
- `[moderate]` **InterProScan member databases** beyond Pfam-A, for genes KEGG cannot place.

## Read QC

- `[low]` **Host removal against T2T-CHM13** (Nurk et al., *Science* 2022), which removes more human
  reads than hg38. `hostile` does this and masks microbial-homologous regions (Constantinides
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
  2022). For GraphMB, see [Kept assembly graphs](#kept-assembly-graphs).
- `[moderate]` **Hi-C binning**: bin3C (DeMaere & Darling, *PeerJ* 2019) or MetaTOR (Baudry et al.,
  *Front Genet* 2019), which uses contact frequency to separate strains and place plasmids.

### Kept assembly graphs

Measured for a graph-aware refinement prototype that was later dropped. A graph that does not match
the assembly still parses but gives no result, so check both limits below.

- **MEGAHIT (the default assembler).** `megahit_toolkit contig2fastg` renames nodes to
  `NODE_<n>_length_..._ID_<n>`, which shares nothing with the `k141_<id>` contigs in
  `final_assembly.fasta`: the name join covered 0 of 7363 contigs. The kept graph is built from the
  largest k, whose contigs are already maximal paths, so only 3.1% of nodes have a neighbour (15.4%
  at k=59, 7.4% at k=99). Joining by sequence with minimap2 would fix the names and still leave
  almost no adjacency.
- **metaSPAdes.** `P` lines carry the FASTA's own scaffold names, so the join covers 5183 of 5183
  contigs, 100% of bp. Only 11.3% of contigs have a neighbour, and those linked and short enough
  to have gone unbinned hold 3.81% of assembled bp.

Both are from one single-sample assembly of 3M CAMI2 oral read pairs. Both assemblers resolve the
graph while building contigs, so the kept graph holds what the assembler could not resolve. It
does not hold the short contigs that binning missed. A graph method first needs a denser graph:
prefer metaSPAdes to MEGAHIT, and Hi-C or long reads to either.

## Taxonomic profiling

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
- `[moderate]` **Compare per-sample assembly and co-assembly at scale.** On three samples they give
  the same number of clean MAGs. The comparison needs about twenty samples, where co-assembly is
  expected to run out of memory.
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
- `[major]` **Learned bin refinement** as an alternative to the current consolidation.

## New data types

- `[major]` **A long-read-first default path** (HiFi or ONT) aimed at complete, circular MAGs.
- `[major]` **Metatranscriptomics** against recovered MAGs.
- `[major]` **Real-time nanopore analysis.**

## Reporting and standards

- `[low]` **Reference recovery report**: how many genomes from a supplied reference set (for
  example published MAGs from the same samples) a run recovered, by skani ANI and aligned fraction.
  It does not depend on CheckM, so it can check CheckM results on real samples.
- `[low]` **Run CPU-hours in `provenance.txt`**: each step's CPU time is in the run manifest; the
  run total is not in the provenance summary.
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
