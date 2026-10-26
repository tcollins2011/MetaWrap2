# MetaWrap2

MetaWrap2 is a Python pipeline for genome-resolved metagenomics. It covers read QC, assembly,
taxonomic profiling, binning, bin refinement and reassembly, cross-sample coverage, dereplication,
quantification and classification, and genome annotation: functional and pathway annotation,
resistance and virulence genes, viruses, strains, pangenomes, phylogenies and differential
abundance. Each module runs on its own.

MetaWrap2 succeeds metaWRAP. Its first release reproduced metaWRAP's pipeline with the same
algorithms, tools and thresholds, so existing analyses could be migrated and checked against it.
Golden tests still pin that baseline. Five defaults have since changed, each based on a benchmark.
The default configuration recovers more near-complete bins than metaWRAP 1 on the three CAMI2 sites
and more clean MAGs on all four real datasets, but fewer near-complete bins on the CAMI3 toy dataset
and on the three-sample real dataset. `--preset gpu` recovers more on every dataset.

![Synthetic communities: CAMI2 and CAMI3](figures/synthetic_benchmarks.png)

- [Installation](installation.md)
- [Benchmarks](benchmarks.md): the measurements and methods.
- [Release notes](https://github.com/ursky/MetaWrap2/blob/main/releases/index.md): per version, what was added and which defaults changed.
- [Planned features](https://github.com/ursky/MetaWrap2/blob/main/PLANNED_FEATURES.md): what is not
  built yet.

## Highlights

- One command per module, `metawrap2 <module>`.
- Conda environments grouped by workflow and built with mamba.
- The commands each module runs are in a `COMMANDS` block at the top of its file. Edit them there or
  override them in `metawrap2.toml`.
- Defaults are chosen by benchmark. `metawrap2 benchmark` runs a declared comparison against a
  dataset with a known answer and commits the result.
- Whole studies run from a sample sheet with `run` or `cohort`. `coverage` builds one depth table
  across all samples, and `dereplicate` merges per-sample bins into one catalogue of organisms.
- GUNC runs by default and flags chimeric bins (sequence from two genomes, clean marker genes),
  which CheckM cannot detect.
- Every run records its configuration, commands, tool versions, database fingerprints, conda
  packages and MetaWrap2 commit, and carries that record through chained runs.
- `--resume` skips a step if it finished, its outputs are present and non-empty, and its inputs
  are unchanged. `metawrap2 status` shows what a resume would do.
- `metawrap2 history` logs every run on the machine. `--diff` compares two.
- One resource budget per run (threads, memory and GPUs), by default a share of the machine. A disk
  check stops a run before it fills the disk.
- `metawrap2 doctor` checks each module's environment, tools and databases. `--fix` rebuilds what
  is broken.

Bins are named `bin_001.fasta` and statistics tables `*.stats.tsv`. The contents match metaWRAP's.
