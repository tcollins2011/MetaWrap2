# Installing MetaWrap2

MetaWrap2 has a light Python core and one small conda environment per module (holding only that module's
tools), so environments resolve quickly and independently. This is the main fix for the install/solver
problems people hit with the old single-environment metaWRAP.

## 1. Install the core package
Clone the repository and install the `metawrap2` core with pip:
```bash
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```
The core needs `python>=3.10` plus `biopython` and the plotting stack (`numpy`, `pandas`, `matplotlib`,
`seaborn`) that MetaWrap2's bundled helper scripts use; `pip` installs these for you. Note that these
helpers deliberately run in *this* interpreter, not inside the per-module conda envs, so install the core
into an environment you are happy to keep around. After this, the `metawrap2` command is available; run
`metawrap2 -h` to see the modules, or `metawrap2 -v` to print the version.

## 2. Create the per-module conda environment(s)
Each module runs in its own conda environment named `metawrap2-<module>` (e.g. `metawrap2-binning`),
built from the matching file in [`src/metawrap2/envs/`](../src/metawrap2/envs/). Create the environment(s) for the module(s) you plan
to use with `metawrap2 install-env`:
```bash
# one module:
metawrap2 install-env binning

# several modules:
metawrap2 install-env assembly binning bin_refinement

# or everything at once:
metawrap2 install-env --all
```
You only need the environments for the modules you actually run. (If you prefer, you can create an
environment directly, e.g. `conda env create -f src/metawrap2/envs/binning.yaml`, or use `mamba` for a faster solver.)

Some modules have optional, newer engines in their own opt-in environments
(`metawrap2-bin_refinement-checkm2`, `metawrap2-classify_bins-gtdbtk`, `metawrap2-annotate_bins-bakta`).
The default environments keep the original tools and behavior, so existing runs reproduce the same
results. See [`src/metawrap2/envs/README.md`](../src/metawrap2/envs/README.md) for the full list of environments and the tools each
one holds.

## 3. Install and configure databases
`metawrap2 install-db` downloads each database, unpacks and indexes it, and writes the matching path into
your config file for you:
```bash
metawrap2 install-db --list                 # what's available, how big, which module needs it
metawrap2 install-db checkm taxdump         # just these two
metawrap2 install-db --all                  # everything (hundreds of GB, many hours)
metawrap2 install-db --all --small          # capped-but-real variants: good for a test run
```
Databases are only needed for the modules that use them, so install just what you need. `--small` picks a
genuinely smaller build where one exists (e.g. the capped 8 GB Kraken2 standard database instead of the
full ~90 GB one) - ideal for verifying the pipeline works, not for a production analysis.

This also runs the `checkm data setRoot` step for you, which is easy to forget and otherwise shows up much
later as a confusing `bin_refinement` failure.

To configure paths by hand instead, copy the example config and edit it:
```bash
metawrap2 config init          # writes a starter ~/.metawrap2/config.toml
metawrap2 config show          # print the resolved settings and which paths are missing
```
There is no more `config-metawrap` file - all configuration (thread defaults, whether to use the
per-module conda envs, and database paths) lives in `metawrap2.toml`. See the
[database installation guide](database_installation.md) for what each database is and how to build one
manually. You can also pass a config explicitly with `--config /path/to/config.toml`.

## 4. Preflight check before running
Before running a module, verify that its tools, conda environment, and any required databases are present:
```bash
metawrap2 check binning
```

## Running modules
Every module is invoked as `metawrap2 <module> [options]` and prints its own help:
```bash
metawrap2 binning -h
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 --concoct \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```
Gzipped read/assembly inputs (`.gz`) are accepted throughout. Every module is also runnable directly, e.g.
`python -m metawrap2.modules.binning ...`. See the [usage tutorial](../Usage_tutorial.md) for a full
worked example.

## Updating
To get the latest version, run `git pull` inside the repository. Because the core is installed with
`pip install -e .`, the pulled changes take effect immediately. If a module's tool set changed, re-create
its environment with `metawrap2 install-env <module>`. Your `~/.metawrap2/config.toml` is outside the
repository, so your database paths are preserved across updates.

## Tinkering with the commands
A design goal of MetaWrap2 is that the exact commands each module runs stay visible and easy to change.
The pipeline is 100% Python now (no Bash, Perl, or R). To change what a module runs, open it in
[`src/metawrap2/modules/<module>.py`](../src/metawrap2/modules/) and edit the `COMMANDS` block (the tool
command templates) and the `CONSTANTS` block (the tunable numbers) at the top of the file. You can also
override surfaced values without touching the source via `~/.metawrap2/config.toml`.
