# Usage

```bash
metawrap2 -h                       # list modules and utilities
metawrap2 <module> -h              # per-module options
mw2 -h                             # the same thing; mw2 is metawrap2 under a shorter name
```

`-h` shows each option's accepted values and default (`[0-100] [default: 70]`).

Each module runs separately and writes to one `-o` output directory:

```bash
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```

Gzipped inputs (`.gz`) are accepted. Modules also run as `python -m metawrap2.modules.binning ...`.
Per-module options are in [Modules](modules.md).

## Running a whole study: `run` and `cohort`

Both take a sample sheet (or a directory of read files) and run every module, with a manifest, run
lock, checkpointing, resource budget and `--resume`. Use `cohort` for more than a few samples: a
co-assembly holds every sample's reads at once, and per-sample assembly keeps each bin traceable
to its sample.

```bash
metawrap2 run    samples.toml -o study/ -t 32        # co-assemble every sample, bin once
metawrap2 cohort samples.toml -o study/ -t 32 -j 4   # assemble and bin each sample, then dereplicate
```

| | `run` | `cohort` |
|---|---|---|
| Assembly | one, from all the reads | one per sample |
| Binning | once, on the co-assembly | once per sample, from one cohort-wide alignment (a depth-table slice and per-assembly BAMs) |
| Result | one bin set | one non-redundant genome catalogue (`MAG_CATALOGUE/`) |
| Steps by default | read_qc, assembly, binning, bin_refinement, reassemble_bins, quant_bins | the same, plus `coverage` and `dereplicate`, minus `reassemble_bins` |

- `coverage` maps each sample's reads once against all samples' contigs and splits the depth table
  per assembly (N alignments instead of N² per study). `--coverage fairy` estimates the table from
  k-mer sketches in minutes, resolving nothing below ~95% ANI; alignment is the default (see
  [Modules](modules.md#coverage) and [Benchmarks](benchmarks.md#1-cohort-support)).
- `dereplicate` collapses per-sample bin sets into one genome per organism, so each organism is
  counted once and `quant_bins` does not split its reads between copies. In a cohort,
  `quant_bins` runs on the catalogue with no `-a`.
- `reassemble_bins` is off by default in a cohort because it would run once per sample. Add it to
  `steps` (or `--steps`) to include it.
- In `run`, `coverage` and `dereplicate` are dropped with a message.
- Under `metawrap2 cohort`, `coassemble = true` in a sheet is overridden with a warning naming `run`.
- `-j N` sets how many samples run a per-sample stage at once, dividing threads and memory. Concurrent
  binning steps on `--device cuda` (the default under `--preset deep`) or `auto` get one GPU each,
  round-robin over the cards the budget allows (`--max-gpus`), and the step banner names them. A
  `--device cuda:N` in the sheet keeps that card for every sample.
- The sheet is read once, when the run starts, so edits to it take effect on the next `--resume`;
  `[limits]` in `metawrap2.toml` is re-read by every step.

## Running on Slurm

`--executor slurm` submits each step of a `run` or `cohort` study as one Slurm batch job. The driver
stays on the login node: it submits a stage's jobs, waits for all of them, then submits the next
stage. It needs `sbatch`, `squeue`, `sacct` (job accounting on) and `scancel`, and the inputs,
databases, conda environments and output directory must be on a filesystem the compute nodes see.

```bash
metawrap2 cohort samples.toml -o study/ --executor slurm -j 8 --slurm-account my_lab
metawrap2 cohort samples.toml -o study/ --executor slurm -j 8 --dry-run   # print the sbatch commands
```

Job options go in `[slurm]` in `metawrap2.toml`, with a `[slurm.modules.<module>]` table per module
that needs something else:

```toml
[slurm]
partition = "general"
account = "my_lab"
time = "24:00:00"
threads = 16
memory_gb = 64

[slurm.modules.assembly]
partition = "bigmem"
memory_gb = 500
time = "72:00:00"
```

- Keys: `partition`, `account`, `qos`, `constraint`, `time`, `gres`, `gpus`, `threads`,
  `memory_gb`, `extra_args` (a list of further sbatch options) and `poll_seconds` (default 30,
  `[slurm]` only). An unset key keeps the cluster's default; unset `threads` and `memory_gb` are 8
  and 32 GB.
- The command line wins over the config: `-t` and `-m` set each job's CPUs and memory, and
  `--slurm-partition`, `--slurm-account`, `--slurm-time`, `--slurm-qos`, `--slurm-constraint`,
  `--slurm-gres` and `--slurm-gpus` set the sbatch option for every job.
- Each job gets its full `threads` and `memory_gb`; `-j` only limits how many jobs of a stage are
  queued at once. The `[limits]` budget of the login node does not cap them, and the global flags
  (`--max-disk-fraction`, `--heartbeat-seconds` and the rest) reach every job.
- A binning step on `--device cuda` (the default under `--preset deep`) requests `gres=gpu:1`
  unless `gres` or `gpus` is set for it. Slurm picks the card, so name no card (`cuda:N`) in the
  sheet.
- Each job's output goes to `LOGS/<step>.log`, after the step's command and its sbatch command. The
  batch scripts are in `.metawrap2/run_steps/slurm/`. The job runs with the driver's environment
  (`--export=ALL`) in the directory the driver was started from, and each step is given the
  config file by absolute path.
- Job IDs are in the manifest. If the driver stops (a closed login session), `--resume` reconnects
  to jobs that are still queued or running, keeps the results of jobs that finished, and submits a
  step again only if its job failed. Ctrl-C cancels the study's running jobs with `scancel`.
- A job fails unless it ends `COMPLETED` with exit code 0, so `TIMEOUT`, `OUT_OF_MEMORY`,
  `NODE_FAIL`, `CANCELLED` and a preempted job that is not requeued stop the study, as a failed
  local step does. The wall time, CPU time and peak memory from `sacct` go into the manifest and
  the study's `provenance.txt`.

## Functional annotation

The KEGG steps need two downloads, once. `install-db` writes each database's path into your config,
so the `--*-db` flags below are only needed to point at a different copy.

```bash
metawrap2 install-db kofam kegg   # KOfam profiles and thresholds; the BRITE hierarchy
```

Each module reads the previous one's table. From a refined bin set and its reads:

```bash
# genes, and KEGG orthologs for them. -a <assembly> instead of -b for the bin-free path.
metawrap2 annotate -o ANNOTATION -b BIN_REFINEMENT/metawrap_70_10_bins -t 32 \
    --kofam --kofam-db $KOFAM_DB

# module completeness, the gene -> process -> pathway hierarchy, and Krona plots
metawrap2 pathways -o PATHWAYS -k ANNOTATION/kos.tsv -a QUANT_BINS/bin_abundance_table.tsv

# how much of each function is in each sample, not just which genome has it
metawrap2 quantify_function -o FUNCTION_ABUNDANCE \
    -k ANNOTATION/kos.tsv -a QUANT_BINS/bin_abundance_table.tsv -m PATHWAYS/modules.tsv \
    --microbial-fraction SINGLEM/microbial_fraction.tsv --biom

# resistance and virulence, with genomic context and an evidence level
metawrap2 amr -o AMR -b BIN_REFINEMENT/metawrap_70_10_bins -t 32 --db $AMRFINDER_DB

# which functions differ between two groups of samples
metawrap2 diff_abundance -o STATS \
    -i FUNCTION_ABUNDANCE/ko_abundance_relative.tsv -g groups.tsv
```

`groups.tsv` has two columns with a header: `sample`, `group`. Samples not listed are not tested.

Caveats:

- Completeness depends on the MAGs. Degrading one CAMI2 genome to 70% completeness (the reporting
  threshold) lost 59% of the modules the intact genome completes, about twice the sequence lost.
  `modules.tsv` includes each genome's CheckM2 completeness.
- Relative abundance is compositional. Use `diff_abundance` or a CoDA method; a t-test on these
  columns is confounded by the sum constraint.
- A MAG-based profile covers only the binned part of the sample. `quantify_function` writes the
  assignable fraction and `bin_refinement` the binned fraction of the assembly; both are usually
  smaller than expected.

## Viruses

Run `viruses` on the assembly, since phage contigs are not in prokaryotic bins. It also writes
`plasmid_contigs.tsv`, which `amr` uses for mobility.

```bash
metawrap2 viruses -o VIRUSES -a ASSEMBLY/final_assembly.fasta -t 32 \
    --genomad-db $GENOMAD_DB --checkv-db $CHECKV_DB
```

## Community share of the bins: `appraise`

```bash
metawrap2 appraise -o APPRAISAL -b BIN_REFINEMENT/metawrap_70_10_bins \
    -r RAW_READS/*_1.fastq.gz -t 32 --metapackage $SINGLEM_METAPACKAGE
```

Outputs:

- `community_represented.tsv`: the per-cell share of each sample's community that the bins
  represent (11.6% of bacterial cells on the test study, though the bins hold most of the assembled
  base pairs).
- `unrepresented_lineages.tsv`: lineages with no genome, most abundant first.

Use forward reads not trimmed below 60 bp (SingleM's window); shorter reads give a result near
zero. The run checks read lengths and warns first.

## Sequencing coverage: `nonpareil`

Run it on every sample's QC'd forward reads before comparing samples:

```bash
metawrap2 nonpareil -o DIVERSITY -t 32 \
    -r READ_QC/*/final_pure_reads_1.fastq
```

When basenames collide, as here (`READ_QC/<sample>/final_pure_reads_1.fastq`), the row label is
taken from the directory, giving rows `ERR011347`, `ERR011348`, `ERR011349`. Files already named per
sample keep their names. A collision the directory cannot resolve is an error.

Outputs:

- `diversity_coverage.tsv`: per sample, coverage, the Nd diversity index, and how many times the
  current effort is needed to reach 95% coverage.
- `diversity_curves.tsv`: rarefaction points with the fitted value beside each, for plotting.

The run prints which samples fall below 50% coverage, and whether the spread is wide enough that
diversity differences partly reflect sequencing depth.

Caveats:

- Coverage is abundance-weighted: C = 0.95 is 95% average coverage by abundance, not 95% of
  species; the missing 5% is the rare tail.
- Under the k-mer kernel, effort is in k-mers (`effort_unit`), with `sequenced_bp` given
  separately. Reads shorter than k are filtered out first and counted.
- One kernel is chosen for the cohort from read lengths. The two differ (50.5% against 83.0%
  coverage on the same 34 bp reads); `--kernel kmer` or `--kernel alignment` forces one. Do not mix
  kernels across samples you compare.

## Global options (before the module name)

| Option | Effect |
|---|---|
| `--dry-run` | print/record the commands without running them |
| `--force`   | overwrite an existing MetaWrap2 output directory |
| `--resume`  | reuse an existing output directory, re-running anything it cannot verify |
| `--verbose-logs` | don't collapse repeated tool messages in the captured logs |
| `--heartbeat-seconds N` | seconds without output before a progress heartbeat; 0 turns it off |
| `--skip-space-check` | start even if the disk-space estimate says there is no room |

Under `run`, `cohort` and `benchmark`, these options and the resource limits below (`--max-threads`,
`--max-disk-fraction` and the rest) are passed on to every step. Each step still gets its own share
of threads and memory from the study's budget.

## Output while a step runs

A tool's output is shown live and copied to `run.stdout` / `run.stderr`. Progress-bar redraws and
repeated messages are collapsed: after 20 copies, one copy is kept every 30 s with its latest
values. `--verbose-logs` turns this off.

MetaWrap2 adds two kinds of line to the screen. Neither is written to `run.stdout` / `run.stderr`.
Under `metawrap2 run` a step's screen output is its log file, so they appear there
(`tail -f study/logs/binning.log`).

Heartbeat. COMEBin training and SemiBin2 clustering can print nothing for hours. After 60 s with no
output, one line per silent minute gives the elapsed time, the length of the silence and the
tool's last line. Before a tool has printed anything, the line says so.

```
[metawrap2] run_comebin.sh is still running, 4h12m elapsed, silent for 2m00s; last line: NN training: epoch=85/200, batches=27, loss=3.45, top1=83.06
```

Low-utilisation warning. If a step has run for five minutes and kept no more than a third of its
granted threads busy on average, it prints one warning:

```
[metawrap2] comebin has kept 9.1 of its 32 granted cores busy over 2h40m. It is working, not stuck: a step this far below its grant is in a phase that does not parallelise, so more threads will not shorten it. Those idle cores are better spent on other samples (-j) than on this step (-t).
```

The check waits five minutes because most tools start with a single-threaded phase (reading an
index, loading a model). `metawrap2 run` uses the same measurement after a step finishes to suggest
a larger `-j`. The warning reads CPU time from `/proc`, so it works on Linux only.

Both are set in `[settings]` in `~/.metawrap2/config.toml`; `--heartbeat-seconds N` overrides the
heartbeat interval for one run:

```toml
[settings]
heartbeat_seconds = 60              # 0 turns heartbeats off
low_utilisation_after_seconds = 300 # 0 turns the warning off
low_utilisation_fraction = 0.33     # the share of the grant a step must keep busy
```

## Resource limits

A `[limits]` block in `~/.metawrap2/config.toml` sets threads, memory, disk ceiling and GPUs for
every module. Defaults: half the cores, 80% of RAM, and stop when a filesystem the run writes to
passes 70% full.

```toml
[limits]
# Share of the cores. max_threads = 32 sets a count instead and wins when both are given.
max_thread_fraction = 0.5
# Share of the RAM; also what an unset -m resolves to.
max_memory_fraction = 0.8
# Stop a step when a filesystem it writes to passes this fraction full.
max_disk_fraction = 0.7
# GPUs to use. Unset means every GPU present; 0 means none.
max_gpus = 0
# true refuses a GPU step on a card another process is using.
gpu_strict = false
```

Each run prints its budget before starting:

```
Resource budget: 16/32 threads, 102.4 GB/128.0 GB memory, disk ceiling 70%, 1/1 GPUs
```

* The budget is divided: `metawrap2 run -j 4` gives each of four concurrent samples a quarter.
  With `--executor slurm` each job has its own allocation instead (see
  [Running on Slurm](#running-on-slurm)).
* `-t` and `-m` are capped to it, with a message naming the setting to raise. Unset, they resolve
  to the budget.
* The disk ceiling is a hard stop. A run will not start above it, and a run that crosses it is
  interrupted like Ctrl-C, so `--resume` continues from the last completed step.

Per-run global flags: `--max-threads`, `--max-thread-fraction`, `--max-memory-gb`,
`--max-memory-fraction`, `--max-disk-fraction`, `--max-gpus` and `--gpu-strict`.
`disk_check_interval` (seconds between disk checks inside a long step) is a `[limits]` setting;
`--disk-check-interval` overrides it. The grant is recorded in the manifest, and `metawrap2 calibrate` compares it with the memory used.

## Editing the commands

In `src/metawrap2/modules/<module>.py`, the `COMMANDS` block holds the tool command templates and
`CONSTANTS` the tunable numbers. Edit them, or override any surfaced value in
`~/.metawrap2/config.toml`.

## Utilities

| Command | Purpose |
|---|---|
| `metawrap2 run <sheet\|dir> -o study/` | run the whole pipeline for a study from one sample sheet |
| `metawrap2 cohort <sheet\|dir> -o study/` | the same, per sample, ending in one dereplicated catalogue |
| `metawrap2 check <module>` | verify tools/databases/conda env for a module |
| `metawrap2 check --for "<module command>"` | the databases and envs one command needs, then `ready: yes` or `ready: no` |
| `metawrap2 install-env <module\|--all> [-t N]` | create (and verify) conda env(s), N at a time |
| `metawrap2 install-db <name\|--all> [-t N] [--small]` | download, index, and configure databases |
| `metawrap2 install-db --for <module>` / `--for-sheet <sheet>` | install the databases a module or a sample sheet needs |
| `metawrap2 doctor [module]` | probe envs, software and databases; print a status map (alias: `test`) |
| `metawrap2 doctor --fix` | recreate only the broken environments, leaving healthy ones alone |
| `metawrap2 status <study>` | what has run, and what `--resume` would do next |
| `metawrap2 history` | every run on this machine, and when |
| `metawrap2 calibrate <study>...` | check the built-in disk/memory estimates against real runs |
| `metawrap2 score -o OUT <bin dir>...` | score any directory of genomes with CheckM2 and GUNC; see [Scoring genomes](score.md) |
| `metawrap2 benchmark list\|fetch\|run\|score` | measure tools and settings against known-truth data |
| `metawrap2 config init\|show` | scaffold or inspect the config |
| `metawrap2 completion bash\|zsh` | print a shell-completion script |

### Checking one command

`check --for` takes a full module command in quotes, as you would type it, and reports only the
databases and conda envs that command needs. Each missing one is printed with the command that
installs it. The exit code is 0 when the command is ready to run and 1 when something is missing.

```
$ metawrap2 check --for "annotate --kofam --cazy -b bins/ -o annotation/"

annotate --kofam --cazy -b bins/ -o annotation/

  OK       KOFAM_DB            /path/to/metawrap2_dbs/KOFAM
  MISSING  DBCAN_DB            -> metawrap2 install-db dbcan
  OK       metawrap2-annotate  conda env

ready: no (1 missing)
```

A leading `metawrap2` and global flags such as `--dry-run` are accepted.

### Checking on a run

`metawrap2 status` reads the run manifest and shows what `--resume` would do, using the same check.
The `--audit` check also runs at the end of every run.

```
metawrap2 status study/
metawrap2 status study/ --step binning     # inputs, outputs, command, and measured cost
metawrap2 status study/ --audit            # is every recorded output still there and non-empty?
```

### Comparing two runs

```
metawrap2 history                    # the 20 most recent runs
metawrap2 history --failed           # what didn't finish
metawrap2 history --diff RUN RUN     # version, env, databases and command line, side by side
```

Get run ids from `metawrap2 history --json`; a unique prefix is enough.

### Benchmarking on your own data

Defaults are chosen with `metawrap2 benchmark` on measured quality, with cost reported alongside.
Results are in [Benchmarks](benchmarks.md).

```
metawrap2 benchmark list                          # datasets, comparisons, and today's defaults
metawrap2 benchmark fetch cami2_gi --dest data/    # prints the size and stops; --yes to start
metawrap2 benchmark run cami2_gi binner_set --data data/
metawrap2 benchmark score cami2_gi binner_set     # ranked on quality, cost beside it
```

- `--max-samples N` or `--sample LABEL` limits an expensive comparison. Every variant gets the same
  subset, and `score` prints which samples were used; runs over different subsets are not comparable.
- `fetch` downloads only with `--yes`, and refuses if the dataset would cross the `[limits]` disk
  ceiling.
- `fetch` uses several ranged connections through `aria2c` where it is available, and one `curl`
  connection otherwise. ENA throttles per connection, so on a large read set this is 6-10x faster.
  Set `METAWRAP2_ARIA2C=` (empty) to force `curl` if a host objects to several connections, or set
  it to the path of an `aria2c` binary.
- Each file is checked against the registry's md5 before it is renamed into place, or decompressed
  in full where no md5 is published (a resume at the wrong offset can give a file of the declared
  size with a corrupt gzip stream). A file that fails is deleted, so
  a retry starts clean, and the error names the check that failed.
- Re-using a work directory needs `metawrap2 --force benchmark run ...`, since each variant is a
  full MetaWrap2 run.
- `run` gives every variant the same threads, memory and inputs.
- Binning comparisons use the dataset's gold-standard assembly, removing assembler variance.
- Results append to `benchmarks/results/<dataset>.tsv`, one long-format table per dataset.

`trimmer` is the cheapest comparison and runs on a laptop (`read_qc` per sample with trim-galore and
fastp, scored on surviving reads, bases and Q30 bases). `assembler` is the most expensive: the
metaSPAdes arm's 3-sample co-assembly took 4 hours of wall time and 140 GB of resident memory. The
committed `tiny` fixture runs in seconds and checks the wiring:

```
metawrap2 benchmark run   tiny trimmer     # the committed fixture: seconds, and proves the wiring
metawrap2 benchmark score tiny trimmer
```

Datasets and comparisons are defined in `benchmarks/datasets.toml` and
`benchmarks/configurations.toml`. See `benchmarks/README.md`.

#### Scoring bins

`score` assesses bins outside the measured run, so the ranked cost is binning alone. Each assessor
runs only when the comparison's metrics need it, and writes beside the bins it scored:

| assessor | metrics | output | needs |
|---|---|---|---|
| AMBER | `amber_f1_bp`, `amber_ari_bp`, `amber_completeness_bp`, `amber_purity_bp` | one `amber/` per variant | `metawrap2 install-env benchmark-amber` |
| CheckM | `bins_near_complete`, `bins_high_quality`, `bins_medium_quality`, `bins_total` | `<bin set>.stats.tsv` | `metawrap2 install-env checkm2`, `metawrap2 install-db checkm2` |
| GUNC | `gunc_pass_fraction` | `<bin set>.gunc/` | `metawrap2 install-env binqc`, `metawrap2 install-db gunc` (`GUNC_DB`) |

- AMBER scores against the dataset's contig-to-genome truth, without marker genes. bp-weighted
  metrics use contig lengths from the binners' assembly (run record or `--assembly`).
- CheckM gives the completeness/contamination tiers papers report. CheckM2 scores by default;
  `--checkm1` (deprecated) scores with CheckM1 on its full genome tree in the `bin_refinement` env.
  Each scorer has its own score cache beside the set, and a `.stats.tsv` written by the other
  scorer is scored again.
- GUNC gives the fraction of a set's bins that are not chimeric. CheckM and AMBER purity are
  insensitive to chimerism (on the `binner_set` arms, a 1.3-point AMBER purity gap came with a 2.8-fold difference in
  chimeric bins). It costs 3.2 to 5.7 CPU-hours per bin set at 48 threads, so it runs only when the
  comparison asks or you pass `--metric gunc_pass_fraction`.

Bin sets with existing output are skipped, so adding a metric runs only the missing tool.
`--reassess` forces a re-run when bins have changed.

For a variant with several binners, ranked metrics come from:

1. its consolidated set (`metawrap_<c>_<x>_bins`, named for the thresholds `bin_refinement` used),
   if there is one;
2. otherwise its only set;
3. otherwise per-set numbers are recorded under `metric__<bin set>`, and `score` lists which
   variants it did not rank and why (summing overlapping sets would count a genome more than once).

## Bin file names

Bins are named `bin_001.fasta` for every binner, from 1, zero-padded. Numbers are per binner:
`metabat2_bins/bin_001` and `maxbin2_bins/bin_001` are unrelated. For CheckM or GTDB-Tk, use
`-x fasta`.

## Other flags

* `binning --device {cpu,cuda,cuda:N,auto}` sets where SemiBin2, VAMB and COMEBin train. The default
  `cpu` gives the same result on any host; `auto` uses a GPU if torch sees one; `--preset deep`
  defaults to `cuda`. VAMB costs 4.9 CPU-hours on a GPU against 66.9 on the CPU, and COMEBin does not
  finish on a CPU.
* Before each GPU step the card is checked and described in one line:

  ```
  GPU preflight: cuda:1 (NVIDIA GPU), utilisation 0%, 15.2 GB of 16.0 GB free, no other process on it
  ```

  An index the host does not have ends the run. A card another process is using gives a warning
  and the step shares it; `--gpu-strict` or `[limits] gpu_strict = true` refuses instead. See [binning](modules.md#binning).
* `reassemble_bins --max-coverage` (default 50) thins each bin's recruited reads to that depth
  (relative to the bin's assembled length) before SPAdes; the bin with the most reads sets the
  reassembly cost. `--max-coverage 0` keeps every read.
