# MetaWrap2 - a flexible pipeline for genome-resolved metagenomic data analysis

> **⚠️ Experimental.** MetaWrap2 is an experimental, in-progress reimagining of metaWRAP. It
> is not yet production-ready: expect rough edges and interface changes, and validate results
> before relying on them. Feedback and contributions are very welcome.

 MetaWrap2 is an **easy-to-use metagenomic wrapper suite** that accomplishes the core tasks of metagenomic analysis from start to finish: read quality control, assembly, visualization, taxonomic profiling, extracting draft genomes (binning), and functional annotation. It takes bin extraction and analysis further with hybrid bin refinement and reassembly (see the module overview below). Each module is a standalone program, so you can use only the parts you need.

 MetaWrap2 is the **Python-native successor to metaWRAP**. It keeps the same algorithms and the same external tools - and therefore the same results - but the pipeline is now entirely Python (no Bash, Perl, or R), with modern packaging and a separate conda environment per module. If you are coming from metaWRAP 1.x, see [Coming from metaWRAP 1.x](#coming-from-metawrap-1x).

## Vision

metaWRAP became popular because it was simple to understand and easy to tweak: each step was a
readable script where you could see and edit the exact command being run. MetaWrap2 keeps that
spirit while fixing what made the original hard to maintain and install. The goals:

- **Same science, modern implementation.** Identical algorithms and tools as metaWRAP; the
  code is now Python 3, tested, and packaged properly.
- **Transparent and hackable.** Every module is one flat file whose commands and tunable
  constants sit right at the top - change a flag in seconds, no framework to learn.
- **Painless installation.** One small conda environment per module (driven by mamba) instead
  of a single monolith that refuses to solve.
- **Reproducible by default.** Every run records its full configuration, exact commands, tool
  versions, and a provenance lineage that follows your data across chained runs.
- **Robust.** Compressed inputs, clear preflight checks, resumable runs, and actionable errors.

 ![General walkthrough of MetaWrap2 modules](https://i.imgur.com/6GqRsm3.png)

## Metagenomic bin recovery

 In addition to being a tool wrapper, MetaWrap2 offers a **powerful hybrid approach** for extracting high-quality draft genomes (bins) from metagenomic data by combining several binners (e.g. metaBAT2, CONCOCT, MaxBin2) and using their individual strengths. The [bin refinement module](https://i.imgur.com/JL665Qo.png) outperforms not only individual binning approaches, but also other bin consolidation programs (Binning_refiner, DAS_Tool) in both synthetic and real datasets. Because this module is standalone, you can use your favorite binning software for the initial predictions - they do not have to come from metaBAT2, CONCOCT, and MaxBin2.

![Bin_refinement performance across microbiome types](https://i.imgur.com/KSk3l2B.jpg)

 MetaWrap2 also includes a bin **reassembly module**, which improves a set of bins by extracting the reads belonging to each bin and **reassembling them** with a more permissive, non-metagenomic assembler. This improves bin N50, modestly increases completion, and drastically reduces contamination.

 > The figures above are from the original metaWRAP benchmarking. They remain valid for MetaWrap2 because the algorithms and external tools are unchanged; only the implementation language and packaging changed.

## Overview of MetaWrap2 modules

#### Metagenomic data pre-processing modules:
	1) read_qc: read trimming and host (e.g. human) read removal
	2) assembly: metagenomic assembly and QC with metaSPAdes or MEGAHIT
	3) kraken2: taxonomy profiling and visualization of reads or contigs

#### Bin processing modules:
	1) binning: initial bin extraction with MaxBin2, metaBAT2, and/or CONCOCT
	2) bin_refinement: consolidate multiple binning predictions into a superior bin set
	3) reassemble_bins: reassemble bins to improve completion and N50, and reduce contamination
	4) quant_bins: estimate bin abundance across samples
	5) blobology: visualize the community and extracted bins with blobplots
	6) classify_bins: conservative but accurate taxonomy prediction for bins
	7) annotate_bins: functionally annotate genes in a set of bins

## System requirements
 Resource requirements vary greatly with the amount of data. Because some of the wrapped software is memory-hungry (KRAKEN2 and metaSPAdes, for example), we recommend 8+ cores and 64GB+ RAM. MetaWrap2 officially targets Linux x64, but also runs on macOS.

## Installation

MetaWrap2 has a light Python core and one small conda environment per module (holding only that module's tools), so environments resolve quickly and independently.

1. Clone the repository and install the core package:
```
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```
2. Create the conda environment(s) for the module(s) you plan to use:
```
# one module:
metawrap2 install-env binning
# or all of them:
metawrap2 install-env --all
```
Environments are created and driven with **mamba** (much faster than conda); `install-env` bootstraps mamba automatically if it isn't already present, and verifies each environment right after building it. Each environment is named `metawrap2-<module>` (e.g. `metawrap2-binning`) and is built from the matching file in [envs/](envs/). (If you manage your own environment instead, set `use_conda_envs = false` in `metawrap2.toml` and MetaWrap2 will call tools straight off your PATH.) Environment files for planned opt-in modern tools (CheckM2, GTDB-Tk, Bakta) are also included for a future release; see [envs/README.md](envs/README.md).

3. Check that everything a module needs is present before running:
```
metawrap2 check binning
```

## Configuration and databases

Copy `metawrap2.toml.example` to `~/.metawrap2/config.toml` (or pass `--config`) and set your database paths there. Databases are only needed for the modules that use them.

|    Database     | Size  |  Used in module |
|:---------------:|:---------------:|:-----:|
| CheckM DB |1.4GB| binning, bin_refinement, reassemble_bins |
| KRAKEN2 standard database|125GB | kraken2 |
| NCBI_nt |71GB |  blobology, classify_bins |
| NCBI_tax |283MB |  blobology, classify_bins |
| Indexed hg38 (host removal) | 20GB |  read_qc |

[Follow this guide for download and configuration instructions](installation/database_installation.md).

## Usage

MetaWrap2 wraps all of its modules under one command; each module can also be called on its own.
```
metawrap2 -h

MetaWrap2 v=2.0.0
Usage: metawrap2 [module] [options]

  Modules:
	read_qc          Read QC (trimming and host/contaminant removal)
	assembly         Metagenomic assembly (metaSPAdes or MEGAHIT)
	kraken2          Taxonomy profiling of reads/contigs with KRAKEN2
	blobology        GC-vs-abundance blobplots of contigs and bins
	binning          Initial binning (metaBAT2, MaxBin2, CONCOCT)
	bin_refinement   Consolidate bin sets into a refined set
	reassemble_bins  Reassemble bins from recruited reads
	quant_bins       Estimate bin abundance across samples
	classify_bins    Assign taxonomy to genomic bins
	annotate_bins    Functional annotation of draft genomes

  Utilities:
	check            Verify tools/databases (and conda envs) for a module
	install-env      Create (and verify) the conda env(s) for module(s)
	test             Test module envs + software, run unit tests, print a status map
	config           Scaffold (init) or inspect (show) the metawrap2.toml config
	completion       Print a shell completion script (bash|zsh)
```

Each module is run separately, and prints its own help. For example:
```
metawrap2 binning -h
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 --concoct \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```
Gzipped read/assembly inputs (`.gz`) are accepted throughout. Every module is also runnable directly, e.g. `python -m metawrap2.modules.binning ...`.

See the [MetaWrap2 usage tutorial](Usage_tutorial.md) for a full worked example.

### Tinkering with the commands
A design goal of MetaWrap2 is that the exact commands each module runs stay visible and easy to change. Open any module in [src/metawrap2/modules/](src/metawrap2/modules/) and you will find a `COMMANDS` block (the tool command templates) and a `CONSTANTS` block (the tunable numbers) at the top of the file - edit those to change what runs. You can also override any surfaced value without touching the source via `~/.metawrap2/config.toml`.

### Provenance and reproducibility
Every module writes provenance into its `-o` output directory, next to the data:
- **`run_config.json`** - the full run record: module, MetaWrap2 version, all resolved parameters (options *and* defaults), config file, conda env, databases, platform, start/end/duration, inputs, the exact commands run, the **versions of the tools actually used**, and the complete lineage.
- **`run_commands.txt`** - the exact commands that were run, in order. The first line is the top-level `metawrap2 <module> ...` invocation; re-running it reproduces the output.
- **`provenance.txt`** - a readable history of every MetaWrap2 step that led to this output. When an input is itself a MetaWrap2 output directory, that input's history is merged in, so provenance propagates across successive runs. Chain several modules (QC → assembly → binning → bin_refinement → reassembly) and the final output's `provenance.txt` traces the whole path; if a step consumed several MetaWrap2 inputs, all of their histories are merged (de-duplicated).
- **`run.stdout` / `run.stderr`** - the full stdout/stderr of every tool the run invoked (streamed live to your screen *and* captured here for later inspection).

Global options (placed **before** the module name) help here:
```
metawrap2 --dry-run binning -a asm.fa -o out ...   # print/record the commands without running them
metawrap2 --force   binning -a asm.fa -o out ...   # overwrite an existing MetaWrap2 output dir
metawrap2 --resume  binning -a asm.fa -o out ...   # continue a run, skipping finished steps
```

### Verifying your installation
```
metawrap2 test              # probe every module's conda env + key software, run unit tests
metawrap2 test binning      # just one module
```
`metawrap2 test` goes into each `metawrap2-<module>` environment, checks that the key tools are present and actually run (flagging anything installed-but-broken), runs the Python unit tests, and prints a per-module **installed / missing / broken** status map. `metawrap2 install-env` verifies each environment right after creating it.

## Coming from metaWRAP 1.x

MetaWrap2 keeps the same modules, tools, and results. To migrate:
- The command is now `metawrap2` instead of `metawrap`. Module names are unchanged.
- Configuration moved from the sourced `config-metawrap` file to `metawrap2.toml`.
- Each module now runs in its own conda environment (`metawrap2-<module>`); create them with `metawrap2 install-env`.
- The deprecated `kraken` (KRAKEN1) module and the `--metabat1` binner were dropped; use `kraken2` and `--metabat2`.
- Output file names and layout are unchanged.

## Repository layout

```
src/metawrap2/            the Python package
  cli.py                  `metawrap2 <module>` dispatcher
  command.py              the one place that runs an external tool (mamba env, logging, provenance)
  config.py provenance.py checkpoint.py checkm.py utils.py constants.py   shared helpers
  modules/<module>.py     one flat, editable "recipe" per pipeline module
  commands/               utility subcommands (check, install-env, test, config, completion)
  scripts/                bundled Python helper scripts
  vendor/                 bundled third-party tools (kept under their own licenses)
envs/<module>.yaml        per-module conda environment (metawrap2-<module>)
tests/                    pytest suite (unit + stub-based integration)
docs/                     mkdocs documentation site
flowcharts/               diagram sources (.xml) + make_flowcharts.py to (re)render the PNGs
containers/               Docker / Apptainer definitions
conda_pkg/                conda recipe
Module_descriptions.md, Usage_tutorial.md, installation/   user guides
```

### Error reporting
MetaWrap2 wraps other bioinformatics programs. If one of those programs fails, first troubleshoot that tool's installation/environment. Because every module is plain Python with the commands visible at the top of the file, you can see exactly how each program is invoked in [src/metawrap2/modules/](src/metawrap2/modules/). When reporting a bug, please include the full output (stdout and stderr) and the version (`metawrap2 -v`).

### Authors
MetaWrap2 is developed by Gherman Uritskiy. MetaWrap2 succeeds metaWRAP. It bundles third-party tools (including Binning_refiner and the blobology lineage) under their own licenses - please also cite the individual tools (e.g. Salmon, MaxBin2, SPAdes, Kraken2) that were integral to your analysis.
