# Conda environments, one per workflow

The MetaWrap2 core is light: `python>=3.10`, `biopython`, and the plotting/table stack
(the `metawrap2` conda package installs these). The heavy third-party tools the modules drive
(SPAdes, CheckM, kraken2, ...) live in the environments described here.

There are five of them, plus four that stay separate for reasons given below. **The grouping is by
workflow, not by module and not by tool.** Modules that run one after another and share tools share an
environment: there is no sense in solving the same aligner-and-samtools graph five times, or in asking
a user to install five things to get from reads to an assembly.

Every grouping below was **solved before it was written down**. The two places the pipeline is split in
an awkward-looking way are splits the solver forced, not preferences:

* **CheckM2 cannot be solved with CheckM1 + GUNC + DAS Tool.** No common solution exists against
  conda-forge + bioconda, so it keeps its own file.
* **GTDB-Tk cannot be solved with taxator-tk.** Same situation, and GTDB-Tk is also the one tool here
  with a 100+ GB reference database, so a user who wants it is opting in anyway.

## Environments

| File | Env name | Modules | Tools |
|---|---|---|---|
| `reads.yaml` | `metawrap2-reads` | read_qc, assembly, reassemble_bins, blobology | trim-galore, fastp, fastqc, bmtagger, bowtie2, spades, megahit, quast, blast, minimap2, bwa, bwa-mem2, samtools |
| `bins.yaml` | `metawrap2-bins` | binning, coverage, dereplicate, quant_bins | metabat2, maxbin2, concoct, minimap2, bwa, bwa-mem2, samtools, skani, sylph, salmon, coverm, fairy |
| `binqc.yaml` | `metawrap2-binqc` | bin_refinement | checkm-genome (CheckM1), pplacer, gunc, das_tool |
| `taxonomy.yaml` | `metawrap2-taxonomy` | classify_bins, kraken2 | kraken2, krona, taxator-tk, blast |
| `annotate.yaml` | `metawrap2-annotate` | annotate_bins | prokka, bakta |

Four tools that used to have an environment each are now inside one of the five above, because keeping
them out only meant a second install command: **GUNC** (which runs by default unless `--skip-gunc`, so
it was never really optional), **DAS Tool** (an alternative consolidator for the step CheckM1 serves),
**CoverM** and **fairy** (single static binaries), and **Bakta** (which solves with Prokka).

## Usage

Normally you do not touch these files - `metawrap2 install-env` reads them:

```bash
metawrap2 install-env --all          # the five workflow envs, each built once
metawrap2 install-env bins binqc     # or just the ones you need
metawrap2 install-env learned        # the opt-in ones, by name
```

Directly, if you prefer:

```bash
mamba env create -f src/metawrap2/envs/bins.yaml
```

You only need the environments for the modules you use. The modules find theirs by name, so nothing
has to be activated first.

## Where these files live

They are inside the package (`src/metawrap2/envs/`) rather than at the repository root so that
they ship with a wheel: `metawrap2 install-env` needs them, and it has to work from a plain
`pip install metawrap2` as well as from a clone.

## Reproducible environments (lockfiles)

The YAMLs here use version lower bounds so they keep solving as upstream tools update. That is
what keeps them installable over time, and also what makes a run from last year hard to
reproduce - an unpinned solve is why QUAST once landed on a Python without `distutils` and
CheckM on a setuptools without `pkg_resources`.

`locks/` holds exact-build lockfiles, one per module per platform. Use them when you want the
environment you had rather than the one that solves today:

```bash
metawrap2 install-env --all --from-lock     # build from the committed lockfiles
metawrap2 install-env binning --lock        # regenerate a lockfile from your built env
```

A lockfile is `conda list --explicit --md5` output: every package's URL and checksum, no extra
tooling needed, and readable in a diff so a changed pin shows up in review. Regenerate them
deliberately when you intend to move versions, and commit the result.

Alternatively, generate cross-platform lockfiles with
[conda-lock](https://github.com/conda/conda-lock):

```bash
conda-lock -f src/metawrap2/envs/binning.yaml -p linux-64 --lockfile src/metawrap2/envs/locks/binning.conda-lock.yml
mamba create -n metawrap2-bins --file src/metawrap2/envs/locks/bins.conda-lock.yml
```

Committing the generated lockfiles pins exact builds; regenerate them when you intend to
update a tool.

## The four that stay separate

Each of these is separate for a stated reason, and the reason is never "it belongs to a different
module". The scores and labels an alternative tool produces are *not* interchangeable with the
default's, so pick one and stay with it within a study.

| File | Env name | Why it is not in one of the five | Selected by |
|---|---|---|---|
| `checkm2.yaml` | `metawrap2-checkm2` | **Forced**: no common solution with CheckM1 + GUNC + DAS Tool | `bin_refinement --checkm2` |
| `gtdbtk.yaml` | `metawrap2-gtdbtk` | **Forced**: no common solution with taxator-tk. Also needs 100+ GB of reference data | `classify_bins --gtdbtk` |
| `learned.yaml` | `metawrap2-learned` | **Size**: pytorch, ~1.5 GB, which nobody who only wants metaBAT2 should wait for | `binning --semibin2` / `--vamb` / `--comebin` |
| `benchmark-amber.yaml` | `metawrap2-benchmark-amber` | Benchmark-only, and pins matplotlib below 3.9 | `benchmark score` |

`learned.yaml` holds all three learned binners in one environment: they are all pytorch applications
and they solve together, so three environments were three copies of the same dependency graph. Its
pytorch is CUDA-capable, so `binning --device cuda` needs no second install - there is no separate GPU
environment to build.

Install one the same way as any other environment:
```bash
metawrap2 install-env checkm2
```

## Why these environments contain no Python

Only the module's *external* tools live in each environment. MetaWrap2's own helper scripts
(`metawrap2.scripts.*`, `metawrap2.vendor.*`) deliberately run in the host interpreter that
metawrap2 itself was installed into - see `metawrap2/pyrun.py`. That keeps the environments
small and avoids installing a second, divergent copy of metawrap2 and its plotting stack into
every one of them.
