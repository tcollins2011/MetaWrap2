# MetaWrap2

A Python pipeline for genome-resolved metagenomics: read QC, assembly, taxonomic profiling,
binning, bin refinement, reassembly, quantification, classification, and annotation. It is a
reimplementation of metaWRAP as small, readable Python modules.

> Experimental (phase 1). This release reproduces the metaWRAP pipeline's results — same
> algorithms, tools, and thresholds — on a rewritten, tested codebase. Expect rough edges and
> interface changes, and validate your results. Coming from metaWRAP 1.x? See [Migrating](#migrating).

Each module is a single flat Python file with the tool commands at the top, so you can open
[`src/metawrap2/modules/binning.py`](src/metawrap2/modules/binning.py), find the `COMMANDS` block,
and change a flag without learning a framework.

## Installation

```bash
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .                # or: pip install -e .   (to edit the code)
metawrap2 quickstart -t 8    # create one small conda env per module, then verify each
```

`quickstart` builds a separate conda environment per module and checks the tools run
(`--binning-only` installs just what genome recovery needs). If you manage tools yourself, set
`use_conda_envs = false` and MetaWrap2 calls them from your `PATH`.

## Databases

Optional, and only for the modules that use them; nothing is downloaded automatically.

```bash
metawrap2 install-db --list                     # each database, size, and which module needs it
metawrap2 install-db checkm host --small -t 3   # enough to recover genomes (capped test-size DBs)
```

| Database | Full size | `--small` | Used by |
|:---|:---:|:---:|:---|
| CheckM | 1.4 GB | – | `binning`, `bin_refinement`, `reassemble_bins` |
| KRAKEN2 standard | ~90 GB | 8 GB | `kraken2` |
| NCBI nt | ~300 GB | ~2 GB | `blobology`, `classify_bins` |
| NCBI taxonomy | 550 MB | – | `blobology`, `classify_bins` |
| Indexed hg38 | ~20 GB | ~9 GB | `read_qc` |

## Running

Run one module at a time. Each is standalone and prints its options with
`metawrap2 <module> --help`. A typical genome-recovery workflow:

```bash
metawrap2 read_qc -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o READ_QC -t 24
metawrap2 assembly -1 READ_QC/final_pure_reads_1.fastq -2 READ_QC/final_pure_reads_2.fastq \
    -o ASSEMBLY -t 24 -m 64 --metaspades
metawrap2 binning -a ASSEMBLY/final_assembly.fasta -o INITIAL_BINNING -t 24 \
    --metabat2 --maxbin2 --concoct reads_1.fastq reads_2.fastq
metawrap2 bin_refinement -o BIN_REFINEMENT -t 24 -c 50 -x 10 \
    -A INITIAL_BINNING/metabat2_bins -B INITIAL_BINNING/maxbin2_bins -C INITIAL_BINNING/concoct_bins
```

`.gz`/`.bz2` inputs are accepted, and paired files are recognised as `_1/_2`, `_R1/_R2`,
`_R1_001/_R2_001`, or `.1/.2`. `bin_refinement` takes bin directories from any binner (up to six),
not just the three MetaWrap2 ships. Bins are named `bin_001.fasta` in every directory; the
`.stats.tsv` beside a bin folder scores each bin (completeness, contamination, GC, lineage, N50).

Global options go before the module name:

```bash
metawrap2 --dry-run  binning -a asm.fa -o out ...   # print the commands, run nothing
metawrap2 --resume   binning -a asm.fa -o out ...   # skip steps whose outputs are still valid
metawrap2 --force    binning -a asm.fa -o out ...   # overwrite an existing output directory
```

To batch several samples in one command, `metawrap2 run` takes a sample sheet or a directory of
reads; see `metawrap2 run --help`. See the [usage tutorial](Usage_tutorial.md) for a worked example
and [Module_descriptions.md](Module_descriptions.md) for what each module does.

## Modules

| | Module | What it does |
|---|---|---|
| Pre-processing | `read_qc` | Read trimming and host (e.g. human) read removal |
| | `assembly` | Assembly and QC with metaSPAdes or MEGAHIT |
| | `kraken2` | Taxonomic profiling of reads or contigs |
| Bins | `binning` | Initial binning with metaBAT2, MaxBin2, and/or CONCOCT |
| | `bin_refinement` | Consolidate several bin sets into a better one |
| | `reassemble_bins` | Reassemble each bin from its own reads |
| | `quant_bins` | Bin abundance across samples (Salmon) |
| | `classify_bins` | Conservative taxonomy for bins |
| | `annotate_bins` | Functional gene annotation |
| | `blobology` | GC-vs-abundance blobplots of contigs and bins |

`bin_refinement`, `classify_bins`, and `annotate_bins` each have an opt-in modern alternative
selected by a flag (`--checkm2`, `--gtdbtk`, `--bakta`); the defaults keep the original tools.
`bin_refinement` combines binners to exploit their individual strengths, and `reassemble_bins`
reassembles each bin from its own reads to improve N50/completeness and cut contamination — the two
modules that are more than a wrapper. Their benchmarking still applies, since phase 1 leaves the
algorithms unchanged (checked on every commit by [golden tests](tests/test_bin_refinement_golden.py)).

![MetaWrap2 modules](https://i.imgur.com/6GqRsm3.png)

## Migrating

Mostly a rename; `metawrap2 config migrate <old-config>` converts an old shell config.

- `metawrap` → `metawrap2`; module names unchanged. The shell config becomes `metawrap2.toml`.
- Each module runs in its own conda env (`metawrap2-<module>`); create with `metawrap2 quickstart`.
- The `kraken` (KRAKEN1) module and `--metabat1` were removed; use `kraken2` and `--metabat2`.
- No scientific defaults changed. `.stats`/`.contigs` files are now `.stats.tsv`/`.contigs.tsv`
  (contents unchanged), and bins are `bin_001.fasta` for every binner — a scripted `*.fa` glob
  needs `*.fasta`.

## What's different

Installation and environments
- Per-module conda environments built and verified by `metawrap2 quickstart` / `install-env`,
  instead of one large environment that struggles to solve.
- `metawrap2 install-db` downloads, unpacks, indexes, and configures databases (`--small` for
  capped test-size variants).
- `metawrap2 doctor` reports an installed / missing / broken map per module and `--fix` rebuilds
  only the broken environments; it also checks that each module's databases are present and complete.
- Exact-build lockfiles (`install-env --lock` / `--from-lock`, committed under `envs/locks/`) for
  reproducible environments.

Running and resuming
- Flexible read-layout detection (`_1/_2`, `_R1/_R2`, `_R1_001/_R2_001`, `.1/.2`; `.gz`/`.bz2`), plus
  `--single-end` and `--interleaved` — not just metaWRAP's `_1.fastq`/`_2.fastq`.
- A verified `--resume`: a step is skipped only if it finished, its outputs are present and
  non-empty, and its inputs are unchanged, tracked in a per-run `manifest.json` with content
  fingerprints. `metawrap2 status <dir>` shows what `--resume` would do next.
- Inputs and intermediates are checked (zero-byte, wrong format, truncated gzip) before a tool runs,
  so failures surface at their cause rather than several modules later.
- Ctrl-C, SIGTERM, and SIGHUP mark the running step interrupted so `--resume` redoes it (a scheduler
  kill no longer leaves a step recorded as still running); a lock file stops two runs from sharing
  one output directory.

Resources
- Disk-space preflight warns when a step is unlikely to fit; one scratch policy (`scratch_dir` /
  `$METAWRAP2_SCRATCH`) covers every module.
- A memory budget is divided across parallel workers, each step's peak memory is measured and
  recorded, and `metawrap2 calibrate` checks the built-in estimates against it.

Provenance and records
- Per run: `run_config.json` (resolved parameters, exact commands, the version of every tool used),
  `run_commands.txt`, `provenance.txt` (merged across chained MetaWrap2 inputs), `run_environment.json`
  (the full conda package list), and `metawrap2.log`.
- `metawrap2 history` keeps a durable per-machine log of every run; `history --diff RUN RUN` compares
  two runs field by field (version, env, databases, command).
- Tool output is streamed to the screen and captured, with repeated messages and progress-bar spam
  collapsed (`--verbose-logs` to keep everything).

Code
- Python only — no `.sh`, `.pl`, or `.R` in the repository, enforced in CI — with a test suite
  including golden tests pinning the refinement algorithm's output and the pipeline's layout.

## Notes

- Requirements scale with data; 8+ cores and 64 GB+ RAM is a reasonable floor. Linux x64 supported,
  macOS works, Python 3.10+.
- Change what a module runs by editing its `COMMANDS`/`CONSTANTS` blocks, or override config values
  in `~/.metawrap2/config.toml`.
- Development: `pip install -e ".[dev]"`; `black`, `isort`, `ruff`, `mypy`, `pytest` run in CI. See
  [CONTRIBUTING.md](CONTRIBUTING.md).
- When reporting a bug, include `metawrap2 -v`, the full output, and the relevant `metawrap2.log`.

## Authors and citation

MetaWrap2 is developed by [Gherman Uritskiy](https://www.linkedin.com/in/gherman-uritskiy-phd-978075b8/)
(Amazon) and [Michael Schatz](http://schatz-lab.org/) (Johns Hopkins University), and succeeds metaWRAP.

There is no MetaWrap2 paper yet. Until one is out, please cite the original metaWRAP publication:
Uritskiy GV, DiRuggiero J, Taylor J. *MetaWRAP — a flexible pipeline for genome-resolved metagenomic
data analysis.* Microbiome (2018).
[doi:10.1186/s40168-018-0541-1](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1) It bundles third-party tools
(including Binning_refiner and the blobology lineage) under their own licenses. Please also cite the
individual tools — Salmon, MaxBin2, metaBAT2, CONCOCT, SPAdes, MEGAHIT, Kraken2, CheckM, and others
— integral to your analysis.
