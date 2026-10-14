# Installation

## The short way

```bash
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .
metawrap2 quickstart -t 8
```

`quickstart` builds and verifies every module's conda environment, then lists the databases each
module needs (550 MB to ~300 GB each). It downloads none unless asked:

- `--databases small`: capped-but-real builds, good for a first run.
- `--databases full`: the full databases.
- Or answer the prompt when running interactively.

`--binning-only` installs just the six modules genome recovery needs. The rest of this page is what
quickstart does, step by step.

## Step by step

MetaWrap2 is a Python core plus five conda environments, one per workflow. Modules that run in
sequence with the same tools share one.

```bash
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```

This installs `metawrap2` and `mw2`, which are the same command.

### Environments

`install-env` builds environments with mamba (installing it if missing) and verifies each one:

```bash
metawrap2 install-env bins            # one environment
metawrap2 install-env binning         # the same one, named by a module that uses it
metawrap2 install-env --all -t 4      # every environment, 4 at a time
```

- Name the environment (`bins`) or any module in it (`binning`, `coverage`). `--all` builds every
  environment once, opt-in ones included; `quickstart` builds the ten defaults plus
  `metawrap2-learned`.
- `-t` sets how many build at once. Each writes its solver output to a log file, printed on failure.
- Existing environments are kept unless you pass `--force`, so re-running after an interruption is
  cheap.
- Names from MetaWrap2 2.0.x, such as `bin_refinement-checkm2`, still work and map to the current
  environment.

The default environments:

| Environment | Modules |
|---|---|
| `metawrap2-reads` | read_qc, assembly, reassemble_bins, blobology, nonpareil |
| `metawrap2-bins` | binning, coverage, dereplicate, quant_bins |
| `metawrap2-binqc` | bin_refinement (CheckM1, GUNC, DAS Tool) |
| `metawrap2-taxonomy` | classify_bins, kraken2, pathways |
| `metawrap2-annotate` | annotate, catalogue, annotate_bins |
| `metawrap2-amr` | amr |
| `metawrap2-viruses` | viruses |
| `metawrap2-strains` | strains |
| `metawrap2-comparative` | pangenome, phylogeny |
| `metawrap2-singlem` | appraise |

### Opt-in environments

| Environment | Why separate | Selected by |
|---|---|---|
| `metawrap2-checkm2` | No common solution with CheckM1, GUNC and DAS Tool | `bin_refinement --checkm2` |
| `metawrap2-gtdbtk` | No common solution with taxator-tk; needs 100+ GB of reference data | `classify_bins --gtdbtk` |
| `metawrap2-learned` | pytorch, about 1.5 GB | `binning --semibin2` / `--vamb` / `--comebin` |
| `metawrap2-rgi` | Pins older dependencies than AMRFinderPlus | `amr --rgi` |
| `metawrap2-tracs` | Pins an older scientific-Python stack than inStrain | `strains --backend tracs` |
| `metawrap2-benchmark-amber` | Benchmark only; pins matplotlib below 3.9 | `benchmark score` |

```bash
metawrap2 install-env learned      # all three learned binners, one environment
metawrap2 install-env checkm2
```

The default binners include SemiBin2, so `quickstart` builds `metawrap2-learned`. Its pytorch is
CUDA-capable, so `binning --device cuda` needs no extra install. The default device is the CPU; see
[Modules](modules.md). [`src/metawrap2/envs/README.md`](https://github.com/ursky/MetaWrap2/blob/main/src/metawrap2/envs/README.md)
lists the tools in each environment.

### k-mer coverage (fairy)

`coverage --coverage fairy` estimates contig depth per sample without alignment. fairy ships in
`metawrap2-bins`. It is opt-in: single-sample, it recovered half the bins minimap2 did. Multi-sample
performance is untested.

```bash
metawrap2 install-env bins   # then: metawrap2 coverage --coverage fairy ...
```

### CheckM2

`--checkm2` needs its own database, separate from the one `install-db checkm` fetches (both are
often called "the CheckM database"):

```bash
metawrap2 install-env checkm2
metawrap2 install-db  checkm2          # 1.7 GB; CheckM2 downloads it, we record where it went
```

### GUNC

GUNC runs by default in `bin_refinement`. It detects chimeric bins, which CheckM can score as
complete and clean. It is in `metawrap2-binqc` and needs its own database:

```bash
metawrap2 install-env binqc
metawrap2 install-db gunc              # 6.7 GB reference; pass --skip-gunc to do without
```

`bin_refinement` checks for this database before CheckM starts, and if it is missing, says to
install it or pass `--skip-gunc`.

## Databases

```bash
metawrap2 install-db --list           # what's available, sizes, which module needs each
metawrap2 install-db --all -t 4       # download, index, and configure everything
metawrap2 install-db --all --small    # capped-but-real builds: fast, good for testing
```

`install-db` downloads, unpacks and indexes each database (including `checkm data setRoot`) and
writes its path to your config. It skips databases already present.

Check that a module is ready:

```bash
metawrap2 check binning            # conda env + databases, with a pass/fail exit code
metawrap2 doctor binning           # env, software and databases, with a status map
metawrap2 doctor --fix             # recreate exactly the environments that are broken
```

`doctor` is an alias for `test`. It also reports each database's content fingerprint, so a
database replaced in place shows as changed.

## Reproducing an environment later

Every shipped environment has a committed exact-build lockfile:

```bash
metawrap2 install-env bins --lock        # write one from an env you already have
metawrap2 install-env bins --from-lock   # build from the committed one
conda create -n metawrap2-bins --file src/metawrap2/envs/locks/bins.linux-64.lock
```

## Configuration

```bash
metawrap2 config init              # writes ~/.metawrap2/config.toml
metawrap2 config show              # prints the resolved settings
```

- `[settings] threads` sets every module's default. `-t N` overrides it; `-t all` uses every core.
- Database paths and settings go in `metawrap2.toml`.
- `use_conda_envs = false` uses tools from your PATH instead of the conda envs.

### Migrating an existing configuration

If you have a shell-style config of `KEY=/path` lines:

```bash
metawrap2 config migrate /path/to/old-config -o ~/.metawrap2/config.toml
```

This keeps database paths, drops obsolete settings with a reason, and flags paths missing on this
machine. The old file is read, not executed.

### mamba must be on PATH when a module runs

Tools run as `mamba run -n metawrap2-<module> ...`, so `mamba` must be on PATH at run time.
Interactive shells get it from `~/.bashrc`; `nohup`, `cron`, `sbatch` and CI jobs usually do not.
For unattended runs:

```bash
export PATH="$HOME/miniforge3/bin:$PATH"     # wherever mamba lives
```

Without it, the run stops before any tool starts and suggests this or `use_conda_envs = false`.
