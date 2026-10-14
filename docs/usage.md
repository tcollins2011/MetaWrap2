# Usage

```bash
metawrap2 -h                       # list modules and utilities
metawrap2 <module> -h              # per-module options
mw2 -h                             # the same thing; mw2 is metawrap2 under a shorter name
```

`-h` shows each option's accepted values and default (`[0-100] [default: 70]`), read from the
parser.

Each module runs separately and writes to one `-o` output directory:

```bash
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```

Gzipped inputs (`.gz`) are accepted. Modules also run directly, e.g.
`python -m metawrap2.modules.binning ...`.

## Running a whole study: `run` and `cohort`

Both take a sample sheet (or a directory of read files) and run every module, with a manifest,
run lock, checkpointing, resource budget and `--resume`. They differ in co-assembly:

```bash
metawrap2 run    samples.toml -o study/ -t 96        # co-assemble every sample, bin once
metawrap2 cohort samples.toml -o study/ -t 96 -j 4   # assemble and bin each sample, then dereplicate
```

`cohort` is `run --no-coassembly`.

| | `run` | `cohort` |
|---|---|---|
| Assembly | one, from all the reads | one per sample |
| Binning | once, on the co-assembly | once per sample, on a slice of one cohort-wide depth table |
| Result | one bin set | one non-redundant genome catalogue (`MAG_CATALOGUE/`) |
| Steps by default | read_qc, assembly, binning, bin_refinement, reassemble_bins, quant_bins | the same, plus `coverage` and `dereplicate`, minus `reassemble_bins` |

The two extra cohort steps:

* `coverage` maps each sample's reads once against all samples' contigs, then splits the depth
  table per assembly. Otherwise cross-sample coverage costs N alignments per assembly, N² per study.
  `--coverage fairy` estimates the table from k-mer sketches in minutes, resolving nothing below
  ~95% ANI. Alignment is the default: single-sample fairy recovered half the bins, and multi-sample
  is untested.
* `dereplicate` collapses per-sample bin sets into one genome per organism. Otherwise a MAG count
  over N samples counts bins, not organisms, and `quant_bins` splits an organism's reads between near-identical copies. In a
  cohort, `quant_bins` runs on the catalogue with no `-a`.

`reassemble_bins`, the most expensive step, is off by default in a cohort since it would run once
per sample. Add it to `steps` (or `--steps`) to include it.

- In `run`, `coverage` and `dereplicate` are dropped with a message.
- Under `metawrap2 cohort`, `coassemble = true` in a sheet is overridden with a warning naming `run`.
- `-j N` sets how many samples run a per-sample stage at once, dividing threads and memory.

## Functional annotation, end to end

Each module reads the previous one's table. The KEGG steps need two downloads, once:

```bash
metawrap2 install-db kofam kegg   # KOfam profiles and thresholds; the BRITE hierarchy
```

From a refined bin set and its reads:

```bash
# genes, and KEGG orthologs for them. -a <assembly> instead of -b for the bin-free path.
metawrap2 annotate -o ANNOTATION -b BIN_REFINEMENT/metawrap_70_10_bins -t 48 \
    --kofam --kofam-db $KOFAM_DB

# module completeness, the gene -> process -> pathway hierarchy, and Krona plots
metawrap2 pathways -o PATHWAYS -k ANNOTATION/kos.tsv -a QUANT_BINS/bin_abundance_table.tsv

# how much of each function is in each sample, not just which genome has it
metawrap2 quantify_function -o FUNCTION_ABUNDANCE \
    -k ANNOTATION/kos.tsv -a QUANT_BINS/bin_abundance_table.tsv -m PATHWAYS/modules.tsv \
    --microbial-fraction SINGLEM/microbial_fraction.tsv --biom

# resistance and virulence, with genomic context and an evidence level
metawrap2 amr -o AMR -b BIN_REFINEMENT/metawrap_70_10_bins -t 48 --db $AMRFINDER_DB

# which functions differ between two groups of samples
metawrap2 diff_abundance -o STATS \
    -i FUNCTION_ABUNDANCE/ko_abundance_relative.tsv -g groups.tsv
```

`groups.tsv` has two columns with a header: `sample`, `group`. Samples not listed are not tested.

Caveats:

- Completeness depends on your MAGs. Degrading one CAMI2 genome to 70% MAG completeness (this
  project's reporting threshold) lost 59% of the modules the intact genome completes, about
  two-fold the sequence lost. `modules.tsv` includes each genome's CheckM2 completeness.
- Relative abundance is compositional. Use `diff_abundance` or a CoDA method; a t-test on these
  columns tests the sum constraint.
- A MAG-based profile covers only the binned part of the sample. `quantify_function` writes the
  assignable fraction and `bin_refinement` the binned fraction of the assembly. Both are usually
  smaller than expected.

## Viruses

Run `viruses` on the assembly, not the bins. Phage contigs are in no prokaryotic bin.

```bash
metawrap2 viruses -o VIRUSES -a ASSEMBLY/final_assembly.fasta -t 48 \
    --genomad-db $GENOMAD_DB --checkv-db $CHECKV_DB
```

It also writes `plasmid_contigs.tsv`, which `amr` uses for mobility.

## What did the bins miss? `appraise`

```bash
metawrap2 appraise -o APPRAISAL -b BIN_REFINEMENT/metawrap_70_10_bins \
    -r RAW_READS/*_1.fastq.gz -t 48 --metapackage $SINGLEM_METAPACKAGE
```

Outputs:

- `community_represented.tsv`: the per-cell share of each sample's community that the bins
  represent.
- `unrepresented_lineages.tsv`: lineages with no genome, most abundant first.

Use forward reads, not trimmed below 60 bp. SingleM compares a 60 bp window, so shorter reads give
a result near zero. The run checks read lengths and warns first.

This is a per-cell share. On the test study, the refined bins hold most of the assembled base pairs
but represent 11.6% of the bacterial cells.

## Is this cohort comparable? `nonpareil`

Run it on every sample's QC'd forward reads before comparing samples:

```bash
metawrap2 nonpareil -o DIVERSITY -t 48 \
    -r READ_QC/*/final_pure_reads_1.fastq
```

Outputs:

- `diversity_coverage.tsv`: one row per sample with coverage, the Nd diversity index, and how many
  times the current effort would be needed to reach 95% coverage.
- `diversity_curves.tsv`: rarefaction points with the fitted value beside each, for plotting.

The run prints which samples fall below 50% coverage, and whether the spread is wide enough that
diversity differences partly reflect sequencing depth.

Caveats:

- Coverage is abundance-weighted. C = 0.95 is 95% average coverage by abundance, not 95% of
  species; the rare tail is what is missing.
- Under the k-mer kernel, effort is in k-mers (`effort_unit`), with `sequenced_bp` given
  separately. Reads shorter than k are filtered out first and counted, since the kernel rejects
  files containing them.
- One kernel is chosen for the cohort from read lengths, because the two disagree: 50.5% against
  83.0% coverage on the same 34 bp reads. `--kernel kmer` or `--kernel alignment` forces one. Do
  not mix kernels across samples you compare.

## Global options (before the module name)

| Option | Effect |
|---|---|
| `--dry-run` | print/record the commands without running them |
| `--force`   | overwrite an existing MetaWrap2 output directory |
| `--resume`  | reuse an existing output directory, re-running anything it cannot verify |
| `--verbose-logs` | don't collapse repeated tool messages in the captured logs |
| `--skip-space-check` | start even if the disk-space estimate says there is no room |

## How much of the machine a run may use

A `[limits]` block in `~/.metawrap2/config.toml` sets threads, memory, disk ceiling and GPUs for
every module. Defaults: half the cores, 80% of RAM, and stop when a filesystem the run writes to
passes 70% full. Fractions let one config work on a laptop or a large node.

```toml
[limits]
max_thread_fraction = 0.5   # or max_threads = 32, which wins when both are given
max_memory_fraction = 0.8   # what an unset -m resolves to
max_disk_fraction = 0.7
max_gpus = 0                # unset means every GPU present; 0 means use none
```

Each run prints its budget before starting:

```
Resource budget: 16/32 threads, 102.4 GB/128.0 GB memory, disk ceiling 70%, 1/1 GPUs
```

* The budget is divided: `metawrap2 run -j 4` gives each of four concurrent samples a quarter.
* `-t` and `-m` are capped to it, with a message naming the setting to raise. Unset, they resolve
  to the budget.
* The disk ceiling is a hard stop. A run will not start above it, and a run that crosses it is
  interrupted like Ctrl-C, so `--resume` continues from the last completed step.

For one run, use the global flags `--max-threads`, `--max-thread-fraction`,
`--max-memory-gb`, `--max-memory-fraction`, `--max-disk-fraction`, `--max-gpus` and
`--disk-check-interval`. `metawrap2 -h` lists all but `--disk-check-interval`. The grant is recorded in the manifest, and
`metawrap2 calibrate` compares it with the memory used, showing whether 80% suits your data.

## Tinkering with the commands

In `src/metawrap2/modules/<module>.py`, the `COMMANDS` block holds the tool command templates and
`CONSTANTS` the tunable numbers. Edit them, or override any surfaced value in
`~/.metawrap2/config.toml`.

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

`metawrap2 status` reads the run manifest and shows what `--resume` would do, using the same check:

```
metawrap2 status study/
metawrap2 status study/ --step binning     # inputs, outputs, command, and measured cost
metawrap2 status study/ --audit            # is every recorded output still there and non-empty?
```

The audit also runs at the end of every run.

### Comparing two runs

```
metawrap2 history                    # the 20 most recent runs
metawrap2 history --failed           # what didn't finish
metawrap2 history --diff RUN RUN     # version, env, databases and command line, side by side
```

Get run ids from `metawrap2 history --json`; a unique prefix is enough.

### Deciding a default for yourself

Defaults are chosen on measured quality, with cost beside it, using `metawrap2 benchmark`. You can
re-run it on your own data.

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
- Re-using a work directory needs `metawrap2 --force benchmark run ...`, since each variant is a
  full MetaWrap2 run.
- `run` gives every variant the same threads, memory and inputs.
- Binning comparisons use the dataset's gold-standard assembly, removing assembler variance.
- Results append to `benchmarks/results/<dataset>.tsv`, which is committed, so a changed default
  is reviewed with its numbers.

#### Scoring bins

`score` assesses bins before ranking them, outside the measured run, so the ranked cost is binning
alone. (CheckM's cost grows with bin count and would penalise productive variants.) Each of three
assessors runs only when the comparison's metrics need it:

- AMBER (`amber_f1_bp`, `amber_ari_bp`, `amber_completeness_bp`, `amber_purity_bp`) scores bins
  against the dataset's known contig-to-genome truth. It is marker-independent, which matters
  because several binners use single-copy markers internally. Install its env once:
  `metawrap2 install-env benchmark-amber`. bp-weighted metrics use contig lengths from the binners'
  assembly, taken from the run record or `--assembly`.
- CheckM (`bins_near_complete`, `bins_high_quality`, `bins_medium_quality`, `bins_total`) gives the
  completeness/contamination tiers papers report. It uses the `bin_refinement` env's CheckM.
- GUNC (`gunc_pass_fraction`) is the fraction of a set's bins that are not chimeric. A bin that
  is 95% one genome and 5% another looks barely contaminated to CheckM and costs AMBER almost
  nothing in purity. On the `binner_set` arms, a 1.3-point AMBER purity gap between two bin sets
  came with a 2.8-fold difference in chimeric bins. GUNC costs 3.2 to 5.7 CPU-hours per bin set at
  48 threads, mostly diamond, so it runs only when the comparison asks or you pass
  `--metric gunc_pass_fraction`. It needs the `metawrap2-binqc` env and
  the `GUNC_DB` database (`metawrap2 install-env binqc`, `metawrap2 install-db gunc`).

Each assessor writes beside the bins it scored: `<bin set>.stats.tsv` for CheckM, `<bin set>.gunc/`
for GUNC, one `amber/` per variant. Bin sets with existing output are skipped, so adding
a metric runs only the missing tool. `--reassess` forces a re-run when bins have changed.

A variant with three binners leaves three overlapping bin sets; summing them counts a genome up to
three times (F1 above 1). Ranked metrics come from:

1. its consolidated set (`metawrap_<c>_<x>_bins`, named for the thresholds `bin_refinement` used),
   if there is one;
2. otherwise its only set;
3. otherwise per-set numbers are recorded under `metric__<bin set>`, and `score` lists which
   variants it did not rank and why.

#### Example comparisons

`trimmer` is the cheapest and runs on a laptop. It runs `read_qc` per sample with trim-galore and
with fastp, scoring surviving reads, bases and Q30 bases counted from the output FASTQ. Its result
made fastp the default.

`assembler` is the most expensive: the 3-sample co-assembly took 4 hours of wall time and 140 GB of
resident memory for the metaSPAdes arm alone. MEGAHIT stays the default; metaSPAdes is more
contiguous and available as `--metaspades`, but costs 3.2x the CPU-hours.

```
metawrap2 benchmark run   tiny trimmer     # the committed fixture: seconds, and proves the wiring
metawrap2 benchmark score tiny trimmer
```

Datasets and comparisons are defined in `benchmarks/datasets.toml` and
`benchmarks/configurations.toml`. See `benchmarks/README.md`.

## Bin file names

Bins are named `bin_001.fasta` for every binner, from 1, zero-padded. Numbers are per binner:
`metabat2_bins/bin_001` and `maxbin2_bins/bin_001` are unrelated. For CheckM or GTDB-Tk, use
`-x fasta`.

## Two flags worth knowing about

* `binning --device {cpu,cuda,cuda:N,auto}` sets where the learned binners (SemiBin2, VAMB,
  COMEBin) train. The default `cpu` reproduces on any host; `auto` uses a GPU if torch sees one.
  VAMB costs 4.9 CPU-hours on a GPU against 66.9 on the CPU, and COMEBin does not finish on a CPU.
* `reassemble_bins --max-coverage` (default 50) thins each bin's recruited reads to that depth
  (relative to the bin's assembled length) before SPAdes. Reassembly cost is set by the bin with the
  most reads. `--max-coverage 0` keeps every read.
