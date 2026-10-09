# MetaWrap2

A Python pipeline for genome-resolved metagenomics: read QC, assembly, taxonomic profiling,
binning, bin refinement, reassembly, dereplication, quantification, classification, and
annotation. It is a reimplementation of metaWRAP as small, readable Python modules, with the
defaults chosen by measurement rather than by preference.

> Experimental. Interfaces may still change, and everything measured so far is measured on one
> dataset — validate your results. Coming from metaWRAP 1.x? See [Migrating](#migrating).

Each module is a single flat Python file with the tool commands at the top, so you can open
[`src/metawrap2/modules/binning.py`](src/metawrap2/modules/binning.py), find the `COMMANDS` block,
and change a flag without learning a framework.

## Measured against metaWRAP 1

On CAMI2 gastrointestinal, three samples, assembling per sample and counting bins that are both
≥70% complete / ≤10% contaminated **and** non-chimeric by GUNC:

| | clean MAGs | CPU-hours *(binning + refinement)* |
|---|---|---|
| metaWRAP 1.3.0 | 48 | 52.5 |
| MetaWrap2, `--preset fast` | 41 | 13.8 |
| **MetaWrap2, shipped default** | **53** | **43.8** |
| MetaWrap2, `--preset slow` | 57 | 80.8 |

Methods, the twelve figures behind them, the defects found along the way, and what is still open
are in [docs/benchmarks.md](docs/benchmarks.md). Every default that moved is justified there, and
the results table it draws from is committed under [benchmarks/results/](benchmarks/results/).

## Installation

```bash
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .                # or: pip install -e .   (to edit the code)
metawrap2 quickstart -t 8    # build the conda environments, then verify each
```

MetaWrap2 is a light Python core plus five conda environments — one per *workflow*, not one per
module, so getting from reads to an assembly is one install rather than four. `quickstart` builds
them and checks the tools run. `--binning-only` installs just what genome recovery needs. If you
manage tools yourself, set `use_conda_envs = false` and MetaWrap2 calls them from your `PATH`.

The default binning tier uses SemiBin2, which lives in a sixth environment with the other learned
binners because they all bring pytorch (~1.5 GB). `quickstart` builds it. To skip it, run with
`--preset fast`, which uses metaBAT2 and MaxBin2 only.

## Databases

Only for the modules that use them; nothing is downloaded automatically.

```bash
metawrap2 install-db --list                          # each database, size, and which module needs it
metawrap2 install-db checkm gunc host --small -t 3   # enough to recover genomes
```

| Database | Full size | `--small` | Used by |
|:---|:---:|:---:|:---|
| CheckM | 1.4 GB | – | `bin_refinement`, `reassemble_bins`, `binning --run-checkm` |
| GUNC | 6.7 GB | – | `bin_refinement`, unless `--skip-gunc` |
| CheckM2 | 1.7 GB | – | `bin_refinement --checkm2` |
| KRAKEN2 standard | ~90 GB | 8 GB | `kraken2` |
| NCBI nt | ~300 GB | ~2 GB | `blobology`, `classify_bins` |
| NCBI taxonomy | 550 MB | – | `blobology`, `classify_bins` |
| Indexed host genome | ~20 GB | ~9 GB | `read_qc` |

## Running

Each module is standalone and prints its options with `metawrap2 <module> --help`. A typical
genome-recovery workflow:

```bash
metawrap2 read_qc -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o READ_QC -t 24
metawrap2 assembly -1 READ_QC/final_pure_reads_1.fastq -2 READ_QC/final_pure_reads_2.fastq \
    -o ASSEMBLY -t 24 -m 64 --metaspades
metawrap2 binning -a ASSEMBLY/final_assembly.fasta -o INITIAL_BINNING -t 24 \
    --preset default reads_1.fastq reads_2.fastq
metawrap2 bin_refinement -o BIN_REFINEMENT -t 24 -c 70 -x 10 \
    -A INITIAL_BINNING/metabat2_bins -B INITIAL_BINNING/semibin2_bins
```

`.gz`/`.bz2` inputs are accepted, and paired files are recognised as `_1/_2`, `_R1/_R2`,
`_R1_001/_R2_001`, or `.1/.2`. `bin_refinement` takes bin directories from any binner (up to six),
not just the ones MetaWrap2 ships. Bins are named `bin_001.fasta` in every directory; the
`.stats.tsv` beside a bin folder scores each bin (completeness, contamination, GC, lineage, N50).

For a whole study, write a sample sheet once. `metawrap2 run` co-assembles the samples;
`metawrap2 cohort` assembles and bins each sample on its own, measures coverage across the whole
cohort in one pass, and dereplicates the per-sample bins into a single genome catalogue:

```bash
metawrap2 run --example > samples.toml   # a commented sheet to edit
metawrap2 cohort samples.toml -o study/ -t 96 -j 4
```

Both are resumable per step. See [docs/usage.md](docs/usage.md) for which to use and how a sheet is
written.

Global options go before the module name:

```bash
metawrap2 --dry-run   binning -a asm.fa -o out ...   # print the commands, run nothing
metawrap2 --resume    binning -a asm.fa -o out ...   # skip steps whose outputs are still valid
metawrap2 --force     binning -a asm.fa -o out ...   # overwrite an existing output directory
metawrap2 --max-threads 32  ...                      # a resource ceiling for the whole run
```

## Modules

| | Module | What it does |
|---|---|---|
| Pre-processing | `read_qc` | Read trimming (fastp or Trim Galore) and host read removal |
| | `assembly` | Assembly and QC with metaSPAdes or MEGAHIT |
| | `kraken2` | Taxonomic profiling of reads or contigs |
| Bins | `binning` | metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB, COMEBin |
| | `coverage` | One cross-sample depth table for a cohort, split per assembly |
| | `bin_refinement` | Consolidate several bin sets into a better one, and check for chimerism |
| | `reassemble_bins` | Reassemble each bin from its own reads |
| | `dereplicate` | Collapse several samples' bins into one non-redundant catalogue |
| | `quant_bins` | Bin abundance across samples (sylph, salmon, CoverM, or the depth table) |
| | `classify_bins` | Conservative taxonomy for bins |
| | `annotate_bins` | Functional gene annotation |
| | `blobology` | GC-vs-abundance blobplots of contigs and bins |

`bin_refinement`, `classify_bins`, and `annotate_bins` each have an opt-in alternative engine
(`--checkm2`, `--gtdbtk`, `--bakta`). `bin_refinement` combines binners to exploit their individual
strengths and `reassemble_bins` reassembles each bin from its own reads to improve N50 and
completeness — the two modules that are more than a wrapper. Both algorithms are unchanged from
metaWRAP and pinned by [golden tests](tests/test_bin_refinement_golden.py).

Picking binners by hand is optional: `--preset {very-fast,fast,default,slow}` chooses a tier by how
much compute you have, and naming a binner yourself always wins over a preset. `--device
{cpu,cuda,auto}` decides where the learned binners train; the default is `cpu` so that re-analysing
a study on different hardware does not quietly produce different MAGs.

![MetaWrap2 modules](https://i.imgur.com/6GqRsm3.png)

## Migrating

Mostly a rename; `metawrap2 config migrate <old-config>` converts an old shell config.

- `metawrap` → `metawrap2` (also installed as `mw2`); module names unchanged. The shell config
  becomes `metawrap2.toml`.
- Modules run inside conda environments created by `metawrap2 quickstart`.
- The `kraken` (KRAKEN1) module and `--metabat1` were removed; use `kraken2` and `--metabat2`.
- `.stats`/`.contigs` files are now `.stats.tsv`/`.contigs.tsv` (contents unchanged), and bins are
  `bin_001.fasta` for every binner — a scripted `*.fa` glob needs `*.fasta`.
- Four scientific defaults moved, each on a measurement in [docs/benchmarks.md](docs/benchmarks.md):
  trimming is fastp rather than Trim Galore, bin abundance is sylph rather than salmon, binning is
  metaBAT2 + SemiBin2 rather than metaBAT2 + MaxBin2, and the depth table's alignment-identity floor
  is 90% rather than jgi's 97%. To reproduce a metaWRAP 1 result, the first three are flags
  (`--trimmer trim-galore`, `--quantifier salmon`, `--preset fast`) and the fourth is
  `--min-identity 97`.
- GUNC now runs on the final bins by default and annotates chimeric ones. `--skip-gunc` opts out,
  along with its database.

## What's different

Measuring
- `metawrap2 benchmark` runs a declared comparison against a known-truth dataset and writes the
  result to a committed TSV, so a moved default arrives in review beside the numbers that justify
  it. Comparisons are data (`benchmarks/configurations.toml`), so adding an arm is a TOML entry.
- Cost is CPU-hours reaped from the child process, never wall time, which rewards a tool for being
  given more cores.

Cohorts
- `metawrap2 run` / `cohort` take a study sheet and execute it, resumable per step.
- `metawrap2 coverage` gives every binner each contig's depth in every sample from N alignments
  rather than N², because cross-sample covariance is what separates two genomes at similar
  abundance in any one sample.
- `metawrap2 dereplicate` collapses a cohort's per-sample bins into a catalogue of organisms, so
  counts mean genomes rather than bins.

Installation and environments
- Five workflow conda environments built and verified by `metawrap2 quickstart` / `install-env`,
  instead of one large environment that struggles to solve.
- `metawrap2 install-db` downloads, unpacks, indexes, and configures databases (`--small` for
  capped test-size variants).
- `metawrap2 doctor` reports an installed / missing / broken map per module and `--fix` rebuilds
  only the broken environments; it also checks each module's databases are present and complete.
- Exact-build lockfiles (`install-env --lock` / `--from-lock`, committed under `envs/locks/`) for
  reproducible environments.

Running and resuming
- Flexible read-layout detection (`_1/_2`, `_R1/_R2`, `_R1_001/_R2_001`, `.1/.2`; `.gz`/`.bz2`),
  plus `--single-end` and `--interleaved` — not just metaWRAP's `_1.fastq`/`_2.fastq`.
- A verified `--resume`: a step is skipped only if it finished, its outputs are present and
  non-empty, and its inputs are unchanged, tracked in a per-run `manifest.json` with content
  fingerprints. `metawrap2 status <dir>` shows what `--resume` would do next.
- Inputs and intermediates are checked (zero-byte, wrong format, truncated gzip) before a tool
  runs, so failures surface at their cause rather than several modules later. `--skip-validation`
  is the escape hatch if a check ever misfires.
- Ctrl-C, SIGTERM, and SIGHUP mark the running step interrupted so `--resume` redoes it; a lock
  file stops two runs from sharing one output directory.

Resources
- One resource budget for a run (`--max-threads`, `--max-memory-gb`, `--max-gpus`, or `[limits]` in
  the config), defaulting to a share of the machine rather than to a number someone has to guess.
  It is divided across parallel workers, and each step's peak memory is measured and recorded.
- Disk-space preflight warns when a step is unlikely to fit, and a run stops rather than filling a
  shared filesystem; one scratch policy (`scratch_dir` / `$METAWRAP2_SCRATCH`) covers every module.
- `metawrap2 calibrate` checks the built-in disk and memory estimates against real runs.

Provenance and records
- Per run: `run_config.json` (resolved parameters, exact commands, the version of every tool used),
  `run_commands.txt`, `provenance.txt` (merged across chained MetaWrap2 inputs),
  `run_environment.json` (the full conda package list), and `metawrap2.log`.
- `metawrap2 history` keeps a durable per-machine log of every run; `history --diff RUN RUN`
  compares two runs field by field (version, env, databases, command).
- Tool output is streamed to the screen and captured, with repeated messages and progress-bar spam
  collapsed (`--verbose-logs` to keep everything).

Code
- Python only — no `.sh`, `.pl`, or `.R` in the repository, enforced in CI — with a test suite
  including golden tests pinning the refinement algorithm's output and the pipeline's layout.
- Every command states its own defaults and accepted values, read off the parser rather than
  written into help text where it drifts.

## Notes

- Requirements scale with data; 8+ cores and 64 GB+ RAM is a reasonable floor. Linux x64 supported,
  macOS works, Python 3.10+.
- Change what a module runs by editing its `COMMANDS`/`CONSTANTS` blocks, or override config values
  in `~/.metawrap2/config.toml`.
- Development: `pip install -e ".[dev]"`; `black`, `isort`, `ruff`, `mypy`, `pytest` run in CI. See
  [CONTRIBUTING.md](CONTRIBUTING.md).
- What is planned, and why, is in [PLANNED_FEATURES.md](PLANNED_FEATURES.md).
- When reporting a bug, include `metawrap2 -v`, the full output, and the relevant `metawrap2.log`.

## Authors and citation

MetaWrap2 is developed by [Gherman Uritskiy](https://www.linkedin.com/in/gherman-uritskiy-phd-978075b8/)
(Amazon) and [Michael Schatz](http://schatz-lab.org/) (Johns Hopkins University), and succeeds metaWRAP.

There is no MetaWrap2 paper yet. Until one is out, please cite the original metaWRAP publication:
Uritskiy GV, DiRuggiero J, Taylor J. *MetaWRAP — a flexible pipeline for genome-resolved metagenomic
data analysis.* Microbiome (2018).
[doi:10.1186/s40168-018-0541-1](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1)

MetaWrap2 bundles third-party tools (including Binning_refiner and the blobology lineage) under
their own licenses. Please also cite the individual tools — metaBAT2, SemiBin2, MaxBin2, CONCOCT,
SPAdes, MEGAHIT, Kraken2, CheckM, GUNC, sylph, and others — integral to your analysis.
