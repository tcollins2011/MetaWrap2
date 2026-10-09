# Usage

```bash
metawrap2 -h                       # list modules and utilities
metawrap2 <module> -h              # per-module options
mw2 -h                             # the same thing; mw2 is metawrap2 under a shorter name
```

`-h` is worth taking literally here. Every option states the values it accepts and the default it will
use if you say nothing - `[0-100] [default: 70]` - read off the parser rather than written into the
help text, so it cannot be out of date.

Each module is run separately and writes to a single `-o` output directory. Example:

```bash
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```

Gzipped inputs (`.gz`) are accepted throughout. Every module is also runnable directly, e.g.
`python -m metawrap2.modules.binning ...`.

## Running a whole study: `run` and `cohort`

Both take one sample sheet (or a directory of read files) and drive every module for you, with the
same manifest, run lock, checkpointing, resource budget and `--resume`. They differ in one decision,
and it is the decision that shapes the whole analysis:

```bash
metawrap2 run    samples.toml -o study/ -t 96        # co-assemble every sample, bin once
metawrap2 cohort samples.toml -o study/ -t 96 -j 4   # assemble and bin each sample, then dereplicate
```

`cohort` is `run --no-coassembly` under a name of its own, because the two produce different things:

| | `run` | `cohort` |
|---|---|---|
| Assembly | one, from all the reads | one per sample |
| Binning | once, on the co-assembly | once per sample, on a slice of one cohort-wide depth table |
| Result | one bin set | one non-redundant genome catalogue (`MAG_CATALOGUE/`) |
| Steps by default | read_qc, assembly, binning, bin_refinement, reassemble_bins, quant_bins | the same, plus `coverage` and `dereplicate`, minus `reassemble_bins` |

The two extra steps are what make the cohort path worth having rather than a loop over samples:

* **`coverage` maps each sample's reads once**, against one catalogue holding every sample's
  contigs, and splits the resulting depth table back out per assembly. A cohort binned with
  cross-sample coverage otherwise costs N alignments per assembly - N² for the study - and that,
  not the binning, is what makes people give up and co-assemble. `metawrap2 coverage --coverage
  fairy` skips the alignment entirely and estimates the same table from k-mer sketches in minutes,
  at the cost of resolving nothing below ~95% ANI. The aligned route remains the default because
  single-sample fairy has been measured and recovered half the bins; its multi-sample claim is what
  this table makes testable.
* **`dereplicate` collapses the per-sample bin sets** into one genome per organism. Without it a
  MAG count over N samples counts bins rather than organisms, and `quant_bins` divides one
  organism's reads between its near-identical copies. In a cohort the quantifier is pointed at the
  catalogue and given no `-a`, for exactly that reason.

`reassemble_bins` is left out of a cohort's defaults deliberately: it is the most expensive step in
the pipeline and here it runs once per sample. Add it to `steps` (or `--steps`) when the extra
completeness is worth paying N times for.

Mixing modes is caught rather than silently obeyed: `coverage` and `dereplicate` are dropped from a
co-assembly with a message (one assembly and one bin set leave them nothing to do), and a sheet
setting `coassemble = true` under `metawrap2 cohort` is overridden with a warning naming `run` as
the way to honour it. `-j N` sets how many samples move through a per-sample stage at once, with the
thread and memory budget divided between them.

## Global options (before the module name)

| Option | Effect |
|---|---|
| `--dry-run` | print/record the commands without running them |
| `--force`   | overwrite an existing MetaWrap2 output directory |
| `--resume`  | reuse an existing output directory, re-running anything it cannot verify |
| `--verbose-logs` | don't collapse repeated tool messages in the captured logs |
| `--skip-space-check` | start even if the disk-space estimate says there is no room |

## How much of the machine a run may use

One `[limits]` block in `~/.metawrap2/config.toml` governs every module: threads, memory, the
disk ceiling, and GPUs. Leave it out and the defaults apply - **half the cores, 80% of the RAM,
and a stop if a filesystem the run writes to passes 70% full**. They are fractions rather than
numbers so the same config works on a laptop and on a large node.

```toml
[limits]
max_thread_fraction = 0.5   # or max_threads = 32, which wins when both are given
max_memory_fraction = 0.8   # what an unset -m resolves to
max_disk_fraction = 0.7
max_gpus = 0                # unset means every GPU present; 0 means use none
```

Every run prints the budget it resolved, once, before anything starts:

```
Resource budget: 16/32 threads, 102.4 GB/128.0 GB memory, disk ceiling 70%, 1/1 GPUs
```

Three things follow from this being one budget rather than a per-module default:

* **It is divided, never multiplied.** `metawrap2 run -j 4` gives each of four concurrent samples
  a quarter of it. Four samples can never ask for four machines' worth of RAM.
* **`-t` and `-m` are capped to it**, with a message saying which setting to raise. Unset, they
  resolve to the budget rather than to a fixed constant.
* **The disk ceiling is measured, not estimated**, so it is a hard stop: a run refuses to start
  above it, and a run that crosses it while running is interrupted the same way Ctrl-C would
  interrupt it - which means `--resume` picks up from the last completed step.

To change the budget for one run without editing the config, use the global flags
(`--max-threads`, `--max-thread-fraction`, `--max-memory-gb`, `--max-memory-fraction`,
`--max-disk-fraction`, `--max-gpus`, `--disk-check-interval`); `metawrap2 -h` lists them. What
each run was granted is recorded in its manifest, and `metawrap2 calibrate` reads it back against
the memory the steps actually used - which is how you find out whether 80% was generous or tight
on your data rather than guessing.

## Tinkering with the commands

Open a module in `src/metawrap2/modules/<module>.py`: the `COMMANDS` block holds the tool
command templates and the `CONSTANTS` block the tunable numbers. Edit them to change what
runs, or override any surfaced value in `~/.metawrap2/config.toml` without touching the source.

## Utilities

| Command | Purpose |
|---|---|
| `metawrap2 run <sheet\|dir> -o study/` | run the whole pipeline for a study from one sample sheet |
| `metawrap2 cohort <sheet\|dir> -o study/` | the same, per sample, ending in one dereplicated catalogue |
| `metawrap2 check <module>` | verify tools/databases/conda env for a module |
| `metawrap2 install-env <module\|--all> [-t N]` | create (and verify) conda env(s), N at a time |
| `metawrap2 install-db <name\|--all> [-t N] [--small]` | download, index, and configure databases |
| `metawrap2 doctor [module]` | probe envs, software and databases; print a status map (alias: `test`) |
| `metawrap2 doctor --fix` | recreate exactly the broken environments, leaving healthy ones alone |
| `metawrap2 status <study>` | what has run, and what `--resume` would do next |
| `metawrap2 history` | every run on this machine, and when |
| `metawrap2 calibrate <study>...` | check the built-in disk/memory estimates against real runs |
| `metawrap2 benchmark list\|fetch\|run\|score` | measure tools and settings against known-truth data |
| `metawrap2 config init\|show` | scaffold or inspect the config |
| `metawrap2 completion bash\|zsh` | print a shell-completion script |

### Checking on a run

`metawrap2 status` reads the study's run manifest, so it can say what a resume would do without
doing it - the `--resume` column comes from the same check `--resume` itself makes:

```
metawrap2 status study/
metawrap2 status study/ --step binning     # inputs, outputs, command, and measured cost
metawrap2 status study/ --audit            # is every recorded output still there and non-empty?
```

The same audit runs automatically at the end of every run, so a step that reported success while
producing nothing is reported then rather than silently mis-skipped later.

### Comparing two runs

```
metawrap2 history                    # the 20 most recent runs
metawrap2 history --failed           # what didn't finish
metawrap2 history --diff RUN RUN     # version, env, databases and command line, side by side
```

Run ids come from `metawrap2 history --json`; a unique prefix is enough.

### Deciding a default for yourself

MetaWrap2's defaults are picked on measured quality, with cost reported beside it, and
`metawrap2 benchmark` is how that was measured - so you can re-measure it on data like yours rather
than taking the choice on trust.

```
metawrap2 benchmark list                          # datasets, comparisons, and today's defaults
metawrap2 benchmark fetch cami2_gi --dest data/    # prints the size and stops; --yes to start
metawrap2 benchmark run cami2_gi binner_set --data data/
metawrap2 benchmark score cami2_gi binner_set     # ranked on quality, cost beside it
```

Comparisons too expensive to run over a whole dataset can be scoped with `--max-samples N` or
`--sample LABEL`; every variant gets the same subset, and which samples took part is recorded and
printed by `score`, since two runs over different subsets are not comparable.

`fetch` never downloads without `--yes`, and refuses outright if the dataset would cross the
`[limits]` disk ceiling. Re-running a comparison over a work directory it has already used needs
`metawrap2 --force benchmark run ...`: each variant is a full MetaWrap2 run and will not overwrite
another one silently. `run` gives every variant of a comparison the same threads, memory and
inputs, because two variants given different grants are not a comparison; binning comparisons run
on the dataset's gold-standard assembly so the binner is measured without an assembler's variance
in front of it. Results append to `benchmarks/results/<dataset>.tsv`, which is committed - a
changed default arrives in review next to the numbers for it.

A binning comparison leaves directories of bins, which are not yet numbers, so `score` assesses them
before it ranks them - deliberately at score time, outside the measured run. The cost being ranked is
the cost of the thing being chosen: `metawrap2 binning` does not run CheckM, and CheckM's cost grows
with the number of bins, so charging it to the run would penalise a variant for being productive.
Three assessors, only run when the comparison's metrics ask for them:

- **AMBER** (`amber_f1_bp`, `amber_ari_bp`, `amber_completeness_bp`, `amber_purity_bp`) scores bins
  against the dataset's known contig-to-genome truth. Marker-gene independent, which is why it is not
  redundant with CheckM - several binners use single-copy markers internally, and scoring those with
  a marker-based tool asks them to grade their own homework. It needs its own env, once:
  `metawrap2 install-env benchmark-amber`. Its bp-weighted metrics are weighted by contig length,
  which comes from the assembly the binners were given: `score` takes that from the run record, or
  from `--assembly`.
- **CheckM** (`bins_near_complete`, `bins_high_quality`, `bins_medium_quality`, `bins_total`) gives
  the completeness/contamination tiers papers report, so a MetaWrap2 number can be put beside a
  published one. It runs in the `bin_refinement` env, the same CheckM the pipeline uses.
- **GUNC** (`gunc_pass_fraction`) is the fraction of a set's bins that are not chimeric, and it
  catches what neither of the others can. CheckM counts single-copy markers, so a bin that is 95% one
  genome and 5% another reads as barely contaminated; AMBER weights by base pairs, so the same bin
  costs almost nothing in purity. On the `binner_set` arms measured here a 1.3-point AMBER purity gap
  between two bin sets came with a 2.8-fold difference in chimeric bins, so an arm ranked on AMBER
  alone can be ranked on the wrong thing. It is also the expensive one - 3.2 to 5.7 CPU-hours per bin
  set at 48 threads, mostly diamond - so it runs only when asked for, either by the comparison or by
  `--metric gunc_pass_fraction`, and it needs the `bin_refinement-gunc` env and the `GUNC_DB`
  database (`metawrap2 install-env binqc`, `metawrap2 install-db gunc`).

Assessment is skipped per assessor for bin sets that already have its output, because scoring is the
cheap, re-runnable half of a benchmark and adding a metric to a table should not restart an hour of
CheckM or an afternoon of GUNC; `--reassess` forces it when the bins themselves have changed. Each
assessor writes beside the bins it scored - `<bin set>.stats.tsv` for CheckM, `<bin set>.gunc/` for
GUNC, one `amber/` per variant - so asking for a second metric later costs only the tool that is
missing.

A variant that ran three binners leaves three overlapping bin sets, and no single one of them is that
variant's answer: adding them up counts the same genome up to three times, which is how an F1 above 1
happens. So a variant's ranked metrics come from its consolidated set (`metawrap_<c>_<x>_bins`, named
for the thresholds `bin_refinement` cut it at) when there is one, or from its only set; otherwise the
per-set numbers are recorded under `metric__<bin set>` and `score` says which variants it would not
rank and why.

Some comparisons are cheap enough to try on a laptop. `trimmer` is the cheapest - it runs `read_qc`
once per sample with trim-galore and once with fastp, and scores both on the reads, bases and Q30
bases that survived, counted off the output FASTQ so neither tool is judged on its own report. It is
also the one that has already moved a default: fastp replaced trim-galore on the strength of it.

`assembler` is at the other end of the range - the 3-sample co-assembly behind the current numbers
took 4 hours of wall time and 140 GB of resident memory for the metaSPAdes arm alone. It confirmed
MEGAHIT rather than moving anything: metaSPAdes is the more contiguous assembler and stays available
as `--metaspades`, but not at 3.2x the CPU-hours for a default.

```
metawrap2 benchmark run   tiny trimmer     # the committed fixture: seconds, and proves the wiring
metawrap2 benchmark score tiny trimmer
```

Both the datasets and the comparisons are data, in `benchmarks/datasets.toml` and
`benchmarks/configurations.toml`, so adding a dataset or a variant to test is an edit rather than a
patch. See `benchmarks/README.md`.

## Bin file names

Bins are named `bin_001.fasta`, for every binner, counted from 1 and zero-padded. The number is a
label within one binner's output, not an identity shared across binners - `metabat2_bins/bin_001`
and `maxbin2_bins/bin_001` are unrelated. If you pass bins to CheckM or GTDB-Tk yourself, the
extension flag is `-x fasta`.

## Two flags worth knowing about

* **`binning --device {cpu,cuda,cuda:N,auto}`** - where the learned binners (SemiBin2, VAMB, COMEBin)
  train. The default is `cpu`, so a run reproduces on any host; `auto` uses a GPU if torch can see one
  in that binner's environment. It is worth passing if you have one: VAMB costs 4.9 CPU-hours on a GPU
  against 66.9 on the CPU, and COMEBin does not finish on a CPU at all.
* **`reassemble_bins --max-coverage`** (default 50) - thins each bin's recruited reads to that depth of
  coverage before SPAdes sees them, measured against that bin's own assembled length. Reassembly is the
  most expensive step in the pipeline and its cost is set by whichever bin recruited the most reads,
  which is rarely the bin you care about. `--max-coverage 0` keeps every read.
