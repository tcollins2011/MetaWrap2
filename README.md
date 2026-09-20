# MetaWrap2

**An ambitious, modern wrapper for genome-resolved metagenomics — built as small, readable Python
modules you can open and change.**

> **⚠️ Phase 1, and experimental.** MetaWrap2 is early. What exists today is the *foundation*: the
> established genome-recovery pipeline, reproduced exactly, on a harness built to carry a lot more
> than it currently does. Expect rough edges and interface changes, and validate results before
> relying on them. Feedback and contributions are very welcome.

MetaWrap2 is **not a port**. It is a new pipeline aimed at what metagenomics looks like now, and
the plan is considerably more ambitious than what is in this repository today. Getting there needs
a foundation that can be trusted, so that is what was built first.

**Where it is now (phase 1).** Ten modules take you from raw reads to annotated draft genomes:
read QC, assembly, taxonomic profiling, binning, hybrid bin refinement, bin reassembly,
quantification, classification, and annotation. Every result is **bit-for-bit what the established
pipeline produced** — same algorithms, same tools, same thresholds, verified by golden tests. That
equivalence is deliberate and it is a checkpoint, not a destination: it means you can migrate,
reproduce your existing analyses, and trust the harness before anything scientific changes
underneath you.

**Where it is going.** Real scientific work — modern tools, better methods, capabilities the
current pipeline simply does not have. See [Roadmap](#roadmap).

**What the foundation already buys you**, whether or not the science changes:

| Problem | What MetaWrap2 does |
|---|---|
| A forty-tool conda environment that will not solve | One small environment per module, built concurrently, each verified — and a one-line `metawrap2 quickstart` that does all of it |
| Five languages, five ways to report an error | Python only. No `.sh`, `.pl`, or `.R` in the repository, enforced in CI |
| Databases installed by copy-pasting from the docs | `metawrap2 install-db` downloads, unpacks, indexes, and writes the config |
| A `--resume` that trusts a marker file | `--resume` *verifies*: outputs still present and non-empty, inputs unchanged |
| A failure surfacing hours and several modules after its cause | Inputs and intermediates are checked for zero bytes, format, and truncation before a tool sees them |
| No answer to "what did I run last month?" | `metawrap2 history` — every run on this machine, with `--diff` between any two |
| Tool versions being whatever was on PATH that day | Every run records each tool's version, the full conda package list, the database fingerprints, and the MetaWrap2 commit |
| A killed job leaving no trace of why | Ctrl-C, SIGTERM and SIGHUP all mark the step interrupted, so `--resume` redoes it |
| Two runs quietly sharing one output directory | A lock file makes the second fail immediately, naming who holds the first |
| Inconsistent output names | One convention everywhere: `bin_001.fasta`, `*.stats.tsv` |
| No tests | 509 tests, including golden tests pinning the refinement algorithm's exact output and the pipeline's output layout |

And it keeps the property that made this style of pipeline worth using in the first place:
**every module is a single flat file with the commands at the top.** Open
[`src/metawrap2/modules/binning.py`](src/metawrap2/modules/binning.py), find the `COMMANDS` block,
change a flag, done. There is no framework to learn.

If you are coming from the 1.x pipeline, start with [Migrating](#migrating). Otherwise, start
here.

---

## Quick start

Everything below is the whole workflow, start to finish. Five minutes of reading.

### 1. Install

```bash
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .                # or `pip install -e .` to hack on it

metawrap2 quickstart -t 8    # builds each module's conda env, then verifies it
```

That is the installation. `quickstart` creates one small conda environment per module — `-t 8`
builds eight at a time — and checks that the tools inside them actually run. `--binning-only`
installs just the six modules genome recovery needs instead of all thirteen.

### 2. Get the databases you need

**Nothing is downloaded unless you ask.** They range from 550 MB to ~300 GB, so quickstart tells
you what each module needs and leaves the decision to you:

| You want to | You need |
|---|---|
| Recover genomes from reads (the common case) | `checkm`, plus `host` if you are removing human reads |
| Profile taxonomy of reads/contigs (`kraken2`) | `kraken2` — 8 GB capped, ~90 GB full |
| Classify bins or draw blobplots | `blast` + `taxdump` — ~2 GB capped, ~300 GB full |

```bash
metawrap2 install-db --list                    # every database, its size, and who needs it
metawrap2 install-db checkm host --small -t 3   # enough to recover genomes
metawrap2 install-db --all                     # everything, full size, when you mean it
```

`--small` fetches capped-but-real variants: right for confirming the pipeline works, not for
publishing. Paths are written into your config automatically. Check it any time with
`metawrap2 doctor`.

### 3. Describe your samples

One file lists your reads. `metawrap2 run --example` prints a commented starter:

```toml
[settings]
coassemble = true
steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]

[modules]
binning = "--metabat2 --maxbin2 --concoct"
bin_refinement = "-c 50 -x 10"

[[samples]]
name = "ERR011347"
r1 = "RAW_READS/ERR011347_1.fastq.gz"
r2 = "RAW_READS/ERR011347_2.fastq.gz"

[[samples]]
name = "ERR011348"
r1 = "RAW_READS/ERR011348_1.fastq.gz"
r2 = "RAW_READS/ERR011348_2.fastq.gz"
```

`[modules]` options are passed to that module exactly as you would type them, so anything in
`metawrap2 binning --help` is available without the driver needing to know about it. `.gz` and
`.bz2` are fine everywhere, and paired files are recognised as `_1/_2`, `_R1/_R2`, `_R1_001/_R2_001`
or `.1/.2`. Single-end and interleaved samples set `layout` and omit `r2`. If you point `run` at a
*directory* instead of a sheet, it pairs up the reads it finds.

### 4. Run it

```bash
metawrap2 run samples.toml -o study/ -t 24 -j 2
```

`-t` is the total thread budget and `-j` how many samples to process at once — the threads are
divided between them, as is `-m` if you set a memory budget. Before committing hours to it:

```bash
metawrap2 --dry-run run samples.toml -o study/    # print every command, run nothing
```

### 5. Read the output

```
study/
  READ_QC/<sample>/         trimmed, host-filtered reads + FastQC reports
  CLEAN_READS/              the QC'd reads, named for the next step
  ASSEMBLY/                 final_assembly.fasta + assembly report
  INITIAL_BINNING/          metabat2_bins/  maxbin2_bins/  concoct_bins/
  BIN_REFINEMENT/           metawrap_50_10_bins/   <- the consolidated bin set
                            metawrap_50_10_bins.stats.tsv
  BIN_REASSEMBLY/           reassembled_bins/ + reassembled_bins.stats.tsv
  QUANT_BINS/               bin_abundance_table.tsv + heatmap
  LOGS/                     one log per step
```

Bins are `bin_001.fasta`, `bin_002.fasta`, … in every directory. The `.stats.tsv` beside a bin
folder scores each bin — completeness, contamination, GC, lineage, N50, size — best first. High-quality
MAGs are the rows with completeness >90 and contamination <5:

```bash
awk -F'\t' '$2>90 && $3<5' study/BIN_REFINEMENT/metawrap_50_10_bins.stats.tsv
```

### 6. Check on it, and pick it back up

```bash
metawrap2 status study/            # what ran, how long it took, what --resume would do next
metawrap2 run samples.toml -o study/ --resume    # continue; finished steps are skipped
```

`--resume` *verifies* rather than assumes: a step is skipped only if it completed, its outputs are
still there and non-empty, and its inputs have not changed. If you interrupt a run — Ctrl-C, or a
scheduler killing it — the step in flight is recorded as interrupted and redone.

### Or do it module by module

The driver is convenience, not a requirement. Every step is an ordinary command, and each module
works standalone on anyone's data:

```bash
metawrap2 read_qc -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o READ_QC -t 24
metawrap2 assembly -1 READ_QC/final_pure_reads_1.fastq -2 READ_QC/final_pure_reads_2.fastq \
    -o ASSEMBLY -t 24 -m 64 --metaspades
metawrap2 binning -a ASSEMBLY/final_assembly.fasta -o INITIAL_BINNING -t 24 \
    --metabat2 --maxbin2 --concoct sampleA_1.fastq sampleA_2.fastq
metawrap2 bin_refinement -o BIN_REFINEMENT -t 24 -c 50 -x 10 \
    -A INITIAL_BINNING/metabat2_bins -B INITIAL_BINNING/maxbin2_bins \
    -C INITIAL_BINNING/concoct_bins
```

`bin_refinement` takes bin directories from *any* binner, up to six of them — the sets do not have
to come from MetaWrap2.

### Next

- [Usage tutorial](Usage_tutorial.md) — the same path worked through on real public data
- [Module_descriptions.md](Module_descriptions.md) — what each module does and why
- `metawrap2 <module> --help` — every option, always current
- [Changing what runs](#changing-what-runs) — editing the commands a module issues

---

## The modules

| | Module | What it does |
|---|---|---|
| **Pre-processing** | `read_qc` | Read trimming and host (e.g. human) read removal |
| | `assembly` | Assembly and assembly QC with metaSPAdes or MEGAHIT |
| | `kraken2` | Taxonomic profiling of reads or contigs, with a kronagram |
| **Bins** | `binning` | Initial bin extraction with metaBAT2, MaxBin2, and/or CONCOCT |
| | `bin_refinement` | Consolidate several bin sets into a better one |
| | `reassemble_bins` | Reassemble each bin from its own reads |
| | `quant_bins` | Bin abundance across samples (Salmon) |
| | `classify_bins` | Conservative taxonomy for bins |
| | `annotate_bins` | Functional gene annotation |
| | `blobology` | GC-vs-abundance blobplots of contigs and bins |

Three modules have an opt-in modern alternative, each in its own environment, selected by a
flag: `bin_refinement --checkm2`, `classify_bins --gtdbtk`, `annotate_bins --bakta`. The
defaults keep the original tools so existing analyses reproduce.

![General walkthrough of MetaWrap2 modules](https://i.imgur.com/6GqRsm3.png)

### The bin refinement and reassembly modules

These are the parts of MetaWrap2 that are not just a wrapper. `bin_refinement` combines several
binners and exploits their individual strengths, outperforming both the individual binners and
other consolidation programs (Binning_refiner, DAS_Tool) on synthetic and real data.
`reassemble_bins` extracts each bin's reads and reassembles them with a permissive
non-metagenomic assembler, improving N50 and completeness while cutting contamination.

![Bin_refinement performance across microbiome types](https://i.imgur.com/KSk3l2B.jpg)

Both are standalone: `bin_refinement` takes bin directories from *any* binner, not just the
three MetaWrap2 ships.

> The figures are from the original benchmarking of these methods, and still apply — in phase 1
> the algorithms and tools are unchanged. That is not taken on trust:
> [golden tests](tests/test_bin_refinement_golden.py) run the real refinement algorithm and compare
> every output file against recorded values, so the equivalence is checked on every commit.

**New in MetaWrap2:** refinement accepts **up to six bin sets**, not just two or three. Two- and
three-set runs produce byte-identical output to before, which the golden tests enforce.

---

## Roadmap

Phase 1 — this repository — is the foundation: the established pipeline reproduced exactly, on a
harness with verified resume, real provenance, a test suite, and reproducible environments. It is
finished enough to use and to migrate onto.

It is also the least interesting part. The reason for building it this way is that the work after
it changes results, and changing results is only safe on top of something that can tell you *what*
changed and *why*: golden tests that fail when output moves, provenance that records the exact
tools and databases a result came from, and `metawrap2 calibrate` and `history --diff` to compare
one run against another.

What comes next is scientific rather than structural — modern tools and methods, and capabilities
the current pipeline does not have at all. Those are not in this repository yet, and this README
will describe them when they are, rather than promising them now. Two things are already settled
about how they will land:

- **Nothing scientific changes silently.** A new method arrives behind a flag, with the previous
  behaviour reachable and the difference documented, so an analysis is never quietly not the
  analysis you ran last month.
- **The pipeline stays legible.** No framework, no plugin system, no indirection between you and
  the command being run. A module is a file you can read.

If there is something you want this to do, opening an issue now is more useful than later.

---

## Migrating

Phase 1 is deliberately a drop-in: the modules, the tools and the results are the same, so the
migration is mostly a rename. The last two bullets are the only things that can break a script.

- `metawrap` → `metawrap2`. Module names are unchanged.
- The sourced shell config becomes `metawrap2.toml`. **`metawrap2 config migrate <old-config>`
  converts it for you**, keeping your database paths and saying what happened to anything that no
  longer applies. `metawrap2 config init` writes a fresh one instead.
- Each module runs in its own env (`metawrap2-<module>`); create them with `metawrap2 install-env`.
- The `kraken` (KRAKEN1) module and the `--metabat1` binner are gone; use `kraken2` and
  `--metabat2`. Typing the old name tells you this rather than printing an argument error.
- Every other flag you already use still exists, and **no scientific default changed** —
  completeness/contamination thresholds, memory and length cut-offs, and the contig dereplication
  mode are all what they were.
- `.stats` files are now `.stats.tsv`, and `.contigs` files `.contigs.tsv`, so an extension says
  what a file is. Contents are unchanged.
- **Bins are named `bin_001.fasta`**, for every binner, counted from 1 and zero-padded. metaWRAP
  gave metaBAT2's bins `bin.1.fa`, MaxBin2's `bin.0.fa` and CONCOCT's its own cluster ids, so the
  same number meant a different thing in each folder — and unpadded numbers sort `bin.1`, `bin.10`,
  `bin.2`, which mattered because refined bins are numbered in sorted filename order. Bin
  *contents* are unchanged; scripts globbing `*.fa` need `*.fasta`.

`metawrap2 config migrate <old-config>` does the config for you. Beyond the rename, these are the
things that are actually new.

### Installation that finishes

One command does it: `metawrap2 quickstart`. Underneath, `install-env` builds the environments
concurrently (`-t`), skipping existing ones unless you pass `--force` — so re-running after an
interruption is cheap — and `install-db` downloads, unpacks and indexes each database, runs the
`checkm data setRoot` step everyone forgets, and writes the resulting paths into your config.
`--small` installs capped-but-real variants for testing. Nothing quickstart does is unavailable
separately; it runs those same commands in order and stops at the first failure with the exact
command to resume from.

**Databases are always an explicit choice.** They run from 550 MB to ~300 GB, and where they go is
a decision about your machine — so quickstart prints what is needed and asks, rather than starting
a multi-hundred-gigabyte download on your behalf. A non-interactive run never asks and never
downloads.

If you'd rather manage your own environment, set `use_conda_envs = false` and MetaWrap2 calls
tools straight off your PATH.

### `--resume` that verifies instead of assuming

Marker-file resume was wrong in two directions: a step whose outputs had been deleted was
skipped anyway, and a step whose *inputs* had changed was also skipped, silently mixing results
from two different inputs. MetaWrap2 keeps a **run manifest** (`manifest.json`) recording each
step's command, duration, and every input and output with a content fingerprint. A step is
skipped only if it completed, its outputs are all still present and non-empty, and none of its
inputs has changed.

Fingerprints sample the first and last megabyte plus the size, so the check costs the same on a
100 GB library as on a small one; `--strict-fingerprints` hashes whole files where that matters.

### Failing fast, at the cause

An empty FASTQ produces an empty assembly, which produces zero bins, which fails inside CheckM
with a message about file extensions — three modules and several hours from the actual problem.
MetaWrap2 checks inputs *and* intermediates for zero bytes, wrong format, and truncation
(including truncated gzip streams) before handing them to a tool, and reports **every** problem
at once rather than one per re-run.

### Seeing what a resume will do, before it does it

```bash
metawrap2 status study/                 # per-step: status, duration, and what --resume would do
metawrap2 status study/ --step binning  # one step's inputs, outputs, command, and measured cost
metawrap2 status study/ --audit         # is every recorded output still present and non-empty?
```

The `--resume` column comes from the same check `--resume` itself makes, so it answers "what
happens if I re-run this?" without re-running it. At the end of every run the same audit runs
automatically: a step that reported success but produced nothing, or declared no outputs at all,
is reported then rather than silently mis-skipped on the next resume.

### Knowing what you ran

- **Per run**, in the output directory: `run_config.json` (every resolved parameter, the exact
  commands, the version of every tool actually used, and a lineage that follows your data across
  chained runs), `run_commands.txt`, `provenance.txt`, `run_environment.json` (the complete conda
  package list), and `metawrap2.log`.
- **Per machine**, forever: `metawrap2 history` reads an append-only log at
  `~/.metawrap2/history.jsonl`. It survives deleting the output directory and is never pruned.

```bash
metawrap2 history                    # the 20 most recent runs
metawrap2 history --failed           # what didn't finish
metawrap2 history --module binning --commands   # exact command lines, to re-run or paste
metawrap2 history --stats            # runs and time spent, per module
```

`metawrap2 history --diff RUN RUN` compares two runs field by field — version, conda env,
databases, command line — which is usually how "it worked last month" gets explained.

`provenance.txt` merges the history of any input that is itself a MetaWrap2 output directory, so
chaining QC → assembly → binning → refinement → reassembly leaves the final output describing
the whole path.

### Reproducing a run later

`metawrap2 install-env <module> --lock` writes an exact-build lockfile
(`conda list --explicit --md5`) to [`src/metawrap2/envs/locks/`](src/metawrap2/envs/locks/), and
`--from-lock` rebuilds from it. Lockfiles for every shipped environment are committed, are
readable in a diff, and need no extra tooling to restore:

```bash
conda create -n metawrap2-binning --file src/metawrap2/envs/locks/binning.linux-64.lock
```

### Logs you can read

Tool output is streamed to your screen *and* captured. Repeated messages are collapsed and
progress-bar spam is filtered, which on a real run took the captured logs from 1.8 MB to 30 KB
without losing anything a person would want; `--verbose-logs` turns that off. Progress is shown
with `tqdm` wherever MetaWrap2 itself is doing the work.

### Surviving the things that go wrong

- **Ctrl-C, SIGTERM and SIGHUP** all mark the running step interrupted, so `--resume` redoes it
  rather than trusting it. The scheduler kill matters more than the Ctrl-C: a job stopped on a
  wall-clock or memory limit used to die outright, leaving the step recorded as still running.
- **One writer per output directory.** A second run writing the same `-o` fails immediately and
  says who holds it, instead of interleaving with the first.
- **Preflight disk estimates** warn when a step looks unlikely to fit — metaSPAdes needs roughly
  12× its input in scratch, and finding that out by filling the filesystem four hours in is
  expensive. It *warns* rather than refuses, because the estimate is a chosen multiplier rather
  than a measurement and you know things about your filesystem that it does not. Set
  `strict_space_check = true` to make it a hard stop (worth it on shared infrastructure), or
  `--skip-space-check` to turn it off.
- **One scratch policy**: `[settings] scratch_dir` (or `$METAWRAP2_SCRATCH`) decides for every
  module, instead of each tool scattering temporary files wherever its own default points.
  Unset, scratch stays inside the output directory.
- **Memory is divided across parallel workers**, so `-m 64 -j 4` gives each worker 16 GB instead
  of four workers each believing they own 64. Each step's *actual* peak memory is measured and
  recorded, and `metawrap2 calibrate` checks the built-in estimates against it.

### `metawrap2 doctor`

```bash
metawrap2 doctor              # per-module installed / missing / broken map
metawrap2 doctor --fix        # recreate exactly the broken environments, leave healthy ones alone
metawrap2 doctor --fix --from-lock
```

`doctor` also checks the databases each module needs — present, complete, and matching what the
config claims — because a healthy environment with a half-downloaded database is not a working
installation. It is the same command as `metawrap2 test`; both names work.

---

## System requirements

Requirements scale with data volume. Some wrapped tools are memory-hungry (KRAKEN2 and
metaSPAdes in particular), so 8+ cores and 64 GB+ RAM is a reasonable floor. Linux x64 is the
supported target; macOS works. Python 3.10+.

## Databases

`metawrap2 install-db` handles these for you. To configure by hand, `metawrap2 config init` and
edit `[databases]`. You only need the databases for the modules you use.

| Database | Full size | `--small` | Used by |
|:---|:---:|:---:|:---|
| CheckM | 1.4 GB | – | `binning`, `bin_refinement`, `reassemble_bins` |
| KRAKEN2 standard | ~90 GB | 8 GB | `kraken2` |
| NCBI nt | ~300 GB | ~2 GB | `blobology`, `classify_bins` |
| NCBI taxonomy | 550 MB | – | `blobology`, `classify_bins` |
| Indexed hg38 | ~20 GB | ~9 GB | `read_qc` |

If your BLAST database is not named `nt`, set `BLASTDB_NAME` alongside `BLASTDB`. Set
`[settings] threads` once and every module uses it as its default (`-t N` overrides; `-t all`
uses every core).

## Changing what runs

Every module in [`src/metawrap2/modules/`](src/metawrap2/modules/) opens with a `COMMANDS` block
(the tool command templates) and a `CONSTANTS` block (the tunable numbers). Edit those. Global
options go **before** the module name:

```bash
metawrap2 --dry-run binning -a asm.fa -o out ...   # print/record the commands, run nothing
metawrap2 --force   binning -a asm.fa -o out ...   # overwrite an existing output directory
metawrap2 --resume  binning -a asm.fa -o out ...   # continue, re-running anything unverifiable
metawrap2 --verbose-logs binning ...               # don't collapse repeated tool output
metawrap2 --skip-space-check binning ...           # start despite the disk-space estimate
```

Anything surfaced in the config can be overridden in `~/.metawrap2/config.toml` without touching
the source. Every module is also runnable directly:
`python -m metawrap2.modules.binning ...`.

## Repository layout

```
src/metawrap2/
  cli.py                  the `metawrap2 <module>` dispatcher
  command.py              the single place an external tool is run (env, logging, provenance)
  modules/<module>.py      one flat, editable file per pipeline module
  commands/                utility subcommands (run, check, doctor, install-env, install-db,
                           history, status, calibrate, config, completion)
  manifest.py              step records, fingerprints and the audit behind --resume
  validate.py              zero-byte / truncation / format checks
  history.py               the durable per-machine run log
  scratch.py               scratch-directory policy and disk preflight
  runlock.py               one writer per output directory
  interrupt.py             SIGTERM/SIGHUP handled like Ctrl-C
  usage.py                 per-step peak memory and CPU measurement
  dbcheck.py               database completeness and content fingerprints
  provenance.py logging.py config.py checkpoint.py checkm.py utils.py constants.py
  envs/<module>.yaml       per-module conda environment
  envs/locks/*.lock        committed exact-build lockfiles
  scripts/                 bundled Python helper scripts
  vendor/                  bundled third-party tools (under their own licenses)
tests/                     pytest suite, incl. golden tests for bin refinement
docs/                      mkdocs site        containers/   Docker / Apptainer
flowcharts/                diagram sources    conda_pkg/    conda recipe
```

## Development

```bash
pip install -e ".[dev]"
pre-commit install --hook-type pre-push
```

`black`, `isort`, `ruff`, `mypy`, and `pytest` must all pass before a push, and CI gates the same
five on Python 3.10–3.13. Formatting is checked, never applied, in CI. See
[CONTRIBUTING.md](CONTRIBUTING.md).

## Reporting problems

MetaWrap2 wraps other programs. If one of them fails, that tool's own installation is the first
thing to check — and because every module is plain Python with the commands at the top of the
file, you can see exactly how it was invoked. When reporting a bug, include `metawrap2 -v`, the
full output, and the relevant `metawrap2.log`.

## Authors and citation

MetaWrap2 is developed by Gherman Uritskiy and succeeds metaWRAP. It bundles third-party tools
(including Binning_refiner and the blobology lineage) under their own licenses. Please also cite
the individual tools — Salmon, MaxBin2, metaBAT2, CONCOCT, SPAdes, MEGAHIT, Kraken2, CheckM and
others — that were integral to your analysis.
