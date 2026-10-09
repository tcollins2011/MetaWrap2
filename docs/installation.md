# Installation

## The short way

```bash
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .
metawrap2 quickstart -t 8
```

`quickstart` builds every module's conda environment and verifies the result. It then tells you
which databases each module needs and how large they are — **it never downloads one you did not
ask for**, because they run from 550 MB to ~300 GB. Pass `--databases small` (capped but real, good
for a first run) or `--databases full` to have it fetch them too, or answer the prompt when running
interactively. `--binning-only` installs just the six modules genome recovery needs.

Everything below is what quickstart does, step by step, if you would rather drive it yourself.

## Migrating an existing configuration

If you already have a shell-style config of `KEY=/path` lines:

```bash
metawrap2 config migrate /path/to/old-config -o ~/.metawrap2/config.toml
```

It keeps your database paths, drops anything that no longer applies and says why, and flags paths
that are not on this machine. The file is read, never executed.

MetaWrap2 has a light Python core plus five conda environments - one per *workflow*, not one per
module. Modules that run one after another and drive the same tools share an environment, so getting
from reads to an assembly is one install rather than four.

```bash
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```

That installs two commands, `metawrap2` and `mw2`. They are the same entry point - every command here
is typed as `metawrap2 <module> <flags>`, and the prefix is the part nobody needs to read.

Create the environment(s) for the module(s) you'll use. Environments are built and driven
with **mamba** (much faster than conda); `install-env` bootstraps mamba if it isn't present
and verifies each environment after building it:

```bash
metawrap2 install-env bins            # one environment
metawrap2 install-env binning         # the same one, named by a module that uses it
metawrap2 install-env --all -t 4      # all five, 4 at a time
```

The five:

| Environment | Covers | Built from |
|---|---|---|
| `metawrap2-reads` | read_qc, assembly, reassemble_bins, blobology | `envs/reads.yaml` |
| `metawrap2-bins` | binning, coverage, dereplicate, quant_bins | `envs/bins.yaml` |
| `metawrap2-binqc` | bin_refinement (CheckM1, GUNC, DAS Tool) | `envs/binqc.yaml` |
| `metawrap2-taxonomy` | classify_bins, kraken2 | `envs/taxonomy.yaml` |
| `metawrap2-annotate` | annotate_bins (Prokka, Bakta) | `envs/annotate.yaml` |

You can name either the environment (`bins`) or any module that uses it (`binning`, `coverage`); both
resolve to the same build, and `--all` builds each environment once rather than once per module that
wants it. `-t` controls how many are built concurrently; each one's solver output goes to its own log
file (printed on failure) rather than interleaving on screen. Existing environments are left
alone unless you pass `--force`, so re-running after an interruption is cheap.

Four more stay separate, each for a reason:

| Environment | Why separate | Selected by |
|---|---|---|
| `metawrap2-checkm2` | No common solution with CheckM1 + GUNC + DAS Tool | `bin_refinement --checkm2` |
| `metawrap2-gtdbtk` | No common solution with taxator-tk; also needs 100+ GB of reference data | `classify_bins --gtdbtk` |
| `metawrap2-learned` | pytorch, ~1.5 GB, which nobody who only wants metaBAT2 should wait for | `binning --semibin2` / `--vamb` / `--comebin` |
| `metawrap2-benchmark-amber` | Benchmark-only, and pins matplotlib below 3.9 | `benchmark score` |

```bash
metawrap2 install-env learned      # all three learned binners, one environment
metawrap2 install-env checkm2
```

`metawrap2-learned` ships a CUDA-capable pytorch, so `binning --device cuda` needs no second install -
there is no separate GPU environment. The default is still the CPU; see `docs/modules.md`.

`coverage-fairy` is the smallest of them - one static binary, no Python stack - and it buys the
k-mer coverage route, which estimates every contig's depth in every sample without aligning
anything. It is opt-in because it has been measured once and lost: single-sample, it recovered half
the bins minimap2 did. Its multi-sample claim is the untested one, and this is how it gets tested:

```bash
metawrap2 install-env bins   # then: metawrap2 coverage --coverage fairy ...
```

`--checkm2` needs a database of its own, which is **not** the one `install-db checkm` fetches -
they share no files, and both are commonly called "the CheckM database":

```bash
metawrap2 install-env checkm2
metawrap2 install-db  checkm2          # 1.7 GB; CheckM2 downloads it, we record where it went
```

One of those environments is **not** opt-in: `bin_refinement-gunc`. GUNC runs by default,
because a bin built from halves of two genomes has a complete and uncontaminated-looking marker
set and CheckM will call it good. It lives in its own environment only because it brings diamond
and prodigal, which have nothing to do with CheckM1 and would make that solve harder:

```bash
metawrap2 install-env binqc
metawrap2 install-db gunc              # 6.7 GB reference; pass --skip-gunc to do without
```

`bin_refinement` checks for the database before it starts rather than after CheckM has spent
hours on it, and names both ways out: install it, or pass `--skip-gunc`.

## Databases

```bash
metawrap2 install-db --list           # what's available, sizes, which module needs each
metawrap2 install-db --all -t 4       # download, index, and configure everything
metawrap2 install-db --all --small    # capped-but-real builds: fast, good for testing
```

This downloads each database, unpacks and indexes it (including the `checkm data setRoot` step
that is easy to forget), and writes the path into your config file. It skips anything already
present, so it is safe to re-run.

Check that a module is ready before running:

```bash
metawrap2 check binning            # conda env + databases, with a pass/fail exit code
metawrap2 doctor binning           # env, software and databases, with a status map
metawrap2 doctor --fix             # recreate exactly the environments that are broken
```

`metawrap2 doctor` is the same command as `metawrap2 test`; both names work. It also reports each
database's content fingerprint, so a database replaced in place shows up as a change rather than an
identical-looking path.

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

Set `[settings] threads` once and every module uses it as its default; `-t N` on a module
overrides it and `-t all` uses every core.

Set database paths and settings in `metawrap2.toml`. To manage tools yourself instead of per
-module conda envs, set `use_conda_envs = false` and MetaWrap2 will call tools off your PATH.

### mamba has to be on PATH for the run, not just for the install

MetaWrap2 runs each tool as `mamba run -n metawrap2-<module> ...`, so `mamba` itself must be
findable at *run* time. An interactive shell usually has it because the conda installer appended a
block to your `~/.bashrc` - and a `nohup`, `cron`, `sbatch` or CI job usually does not, because
those do not read it. Export it explicitly for unattended runs:

```bash
export PATH="$HOME/miniforge3/bin:$PATH"     # wherever mamba lives
```

A run that cannot find it stops before spawning anything, naming both this and
`use_conda_envs = false` as ways out. That is the check working; it used to surface as
`FileNotFoundError: 'mamba'`, which looks like the *tool* is missing rather than conda.
