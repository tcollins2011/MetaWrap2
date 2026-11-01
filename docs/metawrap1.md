# Coming from metaWRAP 1

## Migrating

Most changes are renames. `metawrap2 config migrate <old-config>` converts a metaWRAP 1 shell
config.

- The command is `metawrap2` (also installed as `mw2`). Module names are unchanged. The config file
  is `metawrap2.toml`.
- Tools run in conda environments created by `metawrap2 quickstart`.
- The `kraken` (KRAKEN1) module and `--metabat1` are gone; use `kraken2` and `--metabat2`.
- `.stats` and `.contigs` files are now `.stats.tsv` and `.contigs.tsv`, with the same contents.
  Bins are `bin_001.fasta` for every binner, so a script that globs `*.fa` needs `*.fasta`.
- Five defaults changed, each based on a measurement in [Benchmarks](benchmarks.md).
  To reproduce a metaWRAP 1 result:

  | | MetaWrap2 default | metaWRAP 1 behaviour |
  |---|---|---|
  | Trimming | fastp | `--trimmer trim-galore`, or `--fastp-quality-trim` to make fastp trim rather than filter |
  | Bin abundance | sylph | `--quantifier salmon` |
  | Binners | metaBAT2 + SemiBin2 | `--metabat2 --maxbin2 --concoct` |
  | Depth-table identity floor | 90% | `--min-identity 97` |
  | Bin scorer | CheckM2 | `--checkm1` |

- GUNC runs on the final bins and flags chimeric ones. With `--skip-gunc` it does not run and its
  database is not needed.

## Other changes

Binning
- Six binners, three of them learned (SemiBin2, VAMB, COMEBin), selected by `--preset` tiers, each
  benchmarked end to end. `--preset deep` refines three complementary bin sets.
- GUNC chimera filtering is on by default. A chimeric bin can score as complete and clean in CheckM;
  the clean counts above include only bins that pass both.
- `reassemble_bins` keeps the original bin when reassembly does not improve it, and every refined
  set carries MIMAG tags.

Benchmarking
- `metawrap2 benchmark` runs a comparison against a dataset with a known answer and writes the
  result to a committed table. Comparisons are TOML entries in `benchmarks/configurations.toml`.
- Cost is measured in CPU-hours.

Cohorts
- `metawrap2 run` and `cohort` execute a whole study from a sample sheet, step by step.
- `metawrap2 coverage` gives every binner each contig's depth in every sample, with one alignment
  per sample. Cross-sample coverage separates genomes that have the same abundance within one
  sample.
- `metawrap2 dereplicate` collapses per-sample bins into one catalogue of organisms.

Function and AMR
- `annotate` calls genes once, on bins or on the whole assembly, and later modules read the same
  table. Annotating the assembly also covers contigs that are in no MAG.
- `pathways` computes KEGG module completeness and reports it under both published definitions,
  since published tools disagree. Each score carries flags for cases where completeness can mislead:
  a reaction with no KO, KOs shared with many other modules, or a route search that was cut short.
- Sample-level function can be weighted by MAG abundance, contig depth or gene abundance, and every
  table says which was used.
- KOfam, Pfam and eggNOG run over the same proteins and are written side by side; where two
  disagree, both calls are kept and the overlap is printed.
- `amr` reports each resistance call with its genomic context and an evidence level. Its two
  backends are joined on ARO terms.

Installation
- `metawrap2 install-db` downloads, unpacks, indexes and configures databases (`--small` for test
  sizes).
- `metawrap2 doctor` reports, per module, whether its environment and tools are installed, missing
  or broken, and whether its databases are complete. `--fix` rebuilds only what is broken.
- Lockfiles (`install-env --lock` and `--from-lock`) pin every package build for reproducible
  environments.

Running and resuming
- Read layouts beyond metaWRAP's `_1.fastq`/`_2.fastq` are detected automatically, and
  `--single-end` and `--interleaved` are supported.
- `--resume` skips a step if it finished, its outputs are present and non-empty, and its inputs
  are unchanged. `metawrap2 status <dir>` shows what a resume would do.
- Inputs are checked before a tool starts (empty file, wrong format, truncated gzip).
  `--skip-validation` turns the checks off.
- Ctrl-C, SIGTERM and SIGHUP mark the running step as interrupted, so `--resume` redoes it. A lock
  file stops two runs from writing to the same directory.

Resources
- One resource budget per run (`--max-threads`, `--max-memory-gb`, `--max-gpus`, or `[limits]` in
  the config), by default a share of the machine. Peak memory is recorded for every step.
- A disk-space check before each step stops a run before it fills the disk. One scratch setting
  (`scratch_dir` or `$METAWRAP2_SCRATCH`) applies to every module.
- `metawrap2 calibrate` compares the built-in disk and memory estimates with real runs.

Provenance
- Every run writes `run_config.json` (parameters, commands and tool versions), `run_commands.txt`,
  `provenance.txt` (carried through chained runs), `run_environment.json` (the full conda package
  list) and `metawrap2.log`.
- `metawrap2 history` logs every run on the machine, and `history --diff RUN RUN` compares two.
- Tool output is shown live and saved, with progress bars and repeated messages collapsed
  (`--verbose-logs` keeps everything). A step with no output for 60 s prints its elapsed time and
  the tool's last line (`--heartbeat-seconds N`, 0 to turn it off). On Linux, a step that keeps
  few of its threads busy prints a warning while it runs.

Code
- Every command's `--help` shows its defaults and accepted values.
