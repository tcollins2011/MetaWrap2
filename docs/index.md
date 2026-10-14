# MetaWrap2: a benchmarked pipeline for genome-resolved metagenomics, from reads to MAGs and their functions

MetaWrap2 is a Python pipeline for genome-resolved metagenomics. It covers read QC, assembly,
taxonomic profiling, binning, bin refinement and reassembly, cross-sample coverage, dereplication,
quantification and classification, and then what the genomes do: functional and pathway annotation,
resistance and virulence genes, viruses, strains, pangenomes, phylogenies and differential
abundance. Each module runs on its own, so you use only the parts you need.

MetaWrap2 succeeds metaWRAP. Its first release reproduced metaWRAP's pipeline exactly, with the
same algorithms, tools and thresholds, so existing analyses could be migrated and checked. Golden
tests still pin that baseline. Since then four defaults have changed, each on a benchmark, and the
default configuration recovers more clean high-quality MAGs than metaWRAP 1 for less compute.

![MetaWrap2 against metaWRAP 1](figures/headline_mw1_vs_mw2.png)

- [Installation](installation.md) is the place to start.
- [Benchmarks](benchmarks.md) has the measurements and how they were made.
- [Release notes](releases/index.md) list, per version, what was added and which defaults changed.
- Two plans describe what is built and what remains:
  [functional annotation, pathways, AMR and abundance](plan_functional_annotation.md) and
  [the other modules](plan_other_modules.md).

## Highlights

- One command per module, `metawrap2 <module>`.
- Conda environments grouped by workflow and built with mamba, so an install is a few commands.
- The exact commands each module runs are in a `COMMANDS` block at the top of its file, where you
  can edit them, or override them in `metawrap2.toml`.
- Defaults are chosen by benchmark. `metawrap2 benchmark` runs a declared comparison against a
  dataset with a known answer and commits the result.
- Whole studies run from a sample sheet with `run` or `cohort`. `coverage` builds one depth table
  across all samples, and `dereplicate` merges per-sample bins into one catalogue of organisms.
- GUNC runs by default and flags chimeric bins: bins whose marker genes look clean but whose
  sequence comes from two genomes, which CheckM cannot detect.
- Every run records its configuration, commands, tool versions, database fingerprints, conda
  packages and MetaWrap2 commit, and carries that record through chained runs.
- `--resume` skips a step only if it finished, its outputs are present and non-empty, and its inputs
  are unchanged. `metawrap2 status` shows what a resume would do.
- `metawrap2 history` keeps a log of every run on the machine, and `--diff` compares two.
- One resource budget per run (threads, memory and GPUs), by default a share of the machine, and a
  disk check that stops a run before it fills the disk.
- `metawrap2 doctor` checks each module's environment, tools and databases; `--fix` rebuilds only
  what is broken.

Bins are named `bin_001.fasta` and statistics tables `*.stats.tsv`. Their contents match metaWRAP's.
