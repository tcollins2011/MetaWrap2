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
The core only needs `python>=3.8` and `biopython`. After this, the `metawrap2` command is available; run
`metawrap2 -h` to see the modules, or `metawrap2 -v` to print the version.

## 2. Create the per-module conda environment(s)
Each module runs in its own conda environment named `metawrap2-<module>` (e.g. `metawrap2-binning`),
built from the matching file in [`envs/`](../envs/). Create the environment(s) for the module(s) you plan
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
environment directly, e.g. `conda env create -f envs/binning.yaml`, or use `mamba` for a faster solver.)

Some modules have optional, newer engines in their own opt-in environments
(`metawrap2-bin_refinement-checkm2`, `metawrap2-classify_bins-gtdbtk`, `metawrap2-annotate_bins-bakta`).
The default environments keep the original tools and behavior, so existing runs reproduce the same
results. See [`envs/README.md`](../envs/README.md) for the full list of environments and the tools each
one holds.

## 3. Configure databases
Copy the example configuration and set your database paths:
```bash
mkdir -p ~/.metawrap2
cp metawrap2.toml.example ~/.metawrap2/config.toml
# edit ~/.metawrap2/config.toml and point each database key at your download
```
There is no more `config-metawrap` file - all configuration (thread defaults, whether to use the
per-module conda envs, and database paths) lives in `metawrap2.toml`. Databases are only needed for the
modules that use them. See the [database installation guide](database_installation.md) for how to
download and configure each one. You can also pass a config explicitly with `--config /path/to/config.toml`.

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
