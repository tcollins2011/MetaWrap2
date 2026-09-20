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

MetaWrap2 has a light Python core plus one small conda environment per module.

```bash
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```

Create the environment(s) for the module(s) you'll use. Environments are built and driven
with **mamba** (much faster than conda); `install-env` bootstraps mamba if it isn't present
and verifies each environment after building it:

```bash
metawrap2 install-env binning         # one module
metawrap2 install-env --all -t 8      # everything, 8 environments at a time
```

Each environment is named `metawrap2-<module>` and is built from `src/metawrap2/envs/<module>.yaml`.
`-t` controls how many are built concurrently; each one's solver output goes to its own log
file (printed on failure) rather than interleaving on screen. Existing environments are left
alone unless you pass `--force`, so re-running after an interruption is cheap.

Some modules have an opt-in modern alternative in its own environment, selected by a module
flag: `bin_refinement --checkm2`, `classify_bins --gtdbtk`, `annotate_bins --bakta`. Install
them the same way, e.g. `metawrap2 install-env bin_refinement-checkm2`.

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
metawrap2 install-env binning --lock        # write one from an env you already have
metawrap2 install-env binning --from-lock   # build from the committed one
conda create -n metawrap2-binning --file src/metawrap2/envs/locks/binning.linux-64.lock
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
