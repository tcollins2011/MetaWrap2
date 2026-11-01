# Installation

## Quick install

```bash
# 1. mamba, which solves the conda environments much faster than conda. Needs Python 3.10-3.14.
#    No conda on the machine: install Miniforge, which includes mamba and a current Python.
#      curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
#      bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3 && export PATH="$HOME/miniforge3/bin:$PATH"
#    Anaconda or Miniconda already installed: add mamba to its base environment.
#      conda install -n base -c conda-forge mamba
# 2. MetaWrap2 and its environments.
git clone https://github.com/ursky/MetaWrap2.git && cd MetaWrap2
pip install .
metawrap2 quickstart -t 8
```

On stock RHEL 8 or older Ubuntu, `python3` is 3.9 and `pip install .` fails at once with
`requires a different Python: 3.9.x not in '>=3.10'`. Current Miniforge provides
Python 3.14; `pip install .` took 17 seconds on a fresh install. mamba must also be on `PATH` when a
module runs, which unattended jobs do not set up; see
[mamba on PATH at run time](#mamba-on-path-at-run-time).

`quickstart` builds and verifies every module's conda environment, then lists the databases each
module needs (550 MB to ~300 GB each). It downloads databases only when asked:

- `--databases small`: reduced but real builds, for a first run.
- `--databases full`: the full databases.
- Or answer the prompt when running interactively.

For a reproducible install, add `--from-lock` (see [Lockfiles](#lockfiles)).

`--binning-only` installs the six modules genome recovery needs. `--skip-envs` builds no
environments and reports each database's status. With no `--databases` and no terminal to prompt
on, it downloads nothing and says so.

## Step by step

MetaWrap2 is a Python core plus conda environments, one per workflow. Modules that run in sequence
with the same tools share one.

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
  environment once, opt-in ones included; `quickstart` builds the eleven defaults plus
  `metawrap2-learned`.
- `-t` sets how many build at once. Each writes its solver output to a log file, printed on failure.
- Existing environments are kept unless you pass `--force`, so a re-run after an interruption
  builds only what is missing.
- Names from MetaWrap2 2.0.x, such as `bin_refinement-checkm2`, still work and map to the current
  environment.

The default environments:

| Environment | Modules |
|---|---|
| `metawrap2-reads` | read_qc, assembly, reassemble_bins, blobology, nonpareil |
| `metawrap2-bins` | binning, coverage, dereplicate, quant_bins |
| `metawrap2-binqc` | bin_refinement (GUNC, DAS Tool, and CheckM1 for `--checkm1`) |
| `metawrap2-checkm2` | bin_refinement, reassemble_bins, `binning --run-checkm` (CheckM2, the default scorer). Separate because it has no common solution with CheckM1, GUNC and DAS Tool |
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
| `metawrap2-gtdbtk` | No common solution with taxator-tk; needs 100+ GB of reference data | `classify_bins --gtdbtk` |
| `metawrap2-learned` | pytorch, about 1.5 GB | `binning --semibin2` / `--vamb` / `--comebin` |
| `metawrap2-rgi` | Pins older dependencies than AMRFinderPlus | `amr --rgi` |
| `metawrap2-tracs` | Pins an older scientific-Python stack than inStrain | `strains --backend tracs` |
| `metawrap2-eggnog` | eggNOG-mapper needs Python below 3.12 | `annotate --eggnog` |
| `metawrap2-benchmark-amber` | Benchmark only; pins matplotlib below 3.9 | `benchmark score` |

```bash
metawrap2 install-env learned      # all three learned binners, one environment
metawrap2 install-env gtdbtk
```

The default binners include SemiBin2, so `quickstart` builds `metawrap2-learned`. Its pytorch is
CUDA-capable, so `binning --device cuda` needs no extra install. The default device is the CPU
(cuda under `--preset deep`); see
[Modules](modules.md). [`src/metawrap2/envs/README.md`](https://github.com/ursky/MetaWrap2/blob/main/src/metawrap2/envs/README.md)
lists the tools in each environment.

### k-mer coverage (fairy)

`coverage --coverage fairy` estimates contig depth per sample without alignment. fairy is in
`metawrap2-bins` and is opt-in. On a single sample it recovered half as many bins as minimap2.
Multi-sample performance is untested.

```bash
metawrap2 install-env bins   # then: metawrap2 coverage --coverage fairy ...
```

### CheckM2 and CheckM1

Bins are scored with CheckM2 by default. CheckM2 needs less memory than CheckM1 on its full genome
tree (pplacer holds about 34 GB per call), is faster, and scored the final bins in the MetaWrap2
benchmarks.
`quickstart` builds its environment, and `quickstart --databases` fetches its database:

```bash
metawrap2 install-env checkm2
metawrap2 install-db  checkm2          # 1.7 GB; CheckM2 downloads it, we record where it went
```

`--checkm1` scores with CheckM1 instead, in `bin_refinement`, `binning --run-checkm` and `score`
(`reassemble_bins` uses CheckM2 only). `bin_refinement` and `binning` run it on CheckM1's reduced
tree, about 16 GB per call; `--full-tree` uses the full tree, about 40 GB. It is in
`metawrap2-binqc` and needs its own database, separate from CheckM2's. A run with `--checkm1` stops
at the start if it is not installed:

```bash
metawrap2 install-db checkm            # 1.4 GB
```

### GUNC

GUNC runs by default in `bin_refinement` and detects chimeric bins, which CheckM2 can score as
complete and clean. It is in `metawrap2-binqc` and needs its own database:

```bash
metawrap2 install-env binqc
metawrap2 install-db gunc              # 6.7 GB download, ~13 GB installed; or pass --skip-gunc
```

`bin_refinement` checks for this database, and for the scorer's, before CheckM starts. If it is missing, the error says
to install it or pass `--skip-gunc`.

## Databases

```bash
metawrap2 install-db --list           # what's available, its sizes, which module needs each
metawrap2 install-db --all -t 4       # download, index, and configure everything
metawrap2 install-db --all --small    # capped-but-real builds: fast, good for testing
metawrap2 install-db --all --dir /data/metawrap2_dbs   # a different location (default ~/metawrap2_dbs)
metawrap2 install-db --for viruses --for read_qc       # what these modules need
metawrap2 install-db --for "annotate --kofam --cazy"   # a module and its flags
metawrap2 install-db --for-sheet samples.toml          # what a study's planned steps need
```

`--for` and `--for-sheet` list the databases needed, with download and on-disk sizes and a total,
then install them as named databases are installed. For a sheet, the steps `metawrap2 run` would
plan and the sheet's `[modules]` flags decide what is needed, so `bin_refinement = "--skip-gunc"`
leaves GUNC out. A database already set up in the config is skipped. `--dry-run`, `--small`,
`--dir` and `-t` apply as usual.

`install-db` downloads, unpacks and indexes each database (including `checkm data setRoot`). It
writes each path to your config (`~/.metawrap2/config.toml`, or the file `--config` or
`$METAWRAP2_CONFIG` names) as soon as that database is installed, so later runs find it and an
interrupted install keeps what it finished. `--no-config-update` installs without touching the
config. Databases already present are skipped.

Downloads use several ranged connections when `aria2c` is available and one `curl` connection
otherwise. From a host that throttles per connection, `aria2c` is 6-10x faster. It is optional and
is taken from `PATH` or the `metawrap2-annotate` env. Both paths run the same
checks. These databases publish no checksums, so each archive is decompressed in full before it is
accepted; one that does not decompress is deleted and downloaded again, up to three attempts. A
server that cannot be reached is retried after 15 s, 60 s and 240 s, and the error gives the
downloader's exit code in words. Set `METAWRAP2_ARIA2C=` (empty) to force `curl`, for example if a
host rate-limits several connections.

Before it starts, `install-db` prints the space needed and the disk ceiling (`max_disk_fraction`
under `[limits]` in the config, or `--max-disk-fraction`). For GTDB-Tk, GUNC, eggNOG and SingleM
the archive is deleted only after unpacking, so the space needed counts the download and the
unpacked size together (220 GB for GTDB-Tk).

| Database | Size | `--small` size | Used by |
|:---|:---:|:---:|:---|
| CheckM2 (`checkm2`) | ~3 GB | – | `bin_refinement`, `reassemble_bins`, `binning --run-checkm` |
| CheckM1 (`checkm`) | ~1.4 GB | – | `bin_refinement` and `binning --run-checkm` with `--checkm1` |
| GUNC (`gunc`) | ~13 GB | – | `bin_refinement`, unless `--skip-gunc` |
| KRAKEN2 standard (`kraken2`) | ~90 GB | 8 GB | `kraken2`, and `--bracken` (the k-mer distributions ship with it) |
| NCBI nt (`blast`) | ~300 GB | ~2 GB | `blobology`, `classify_bins` |
| NCBI taxonomy (`taxdump`) | ~550 MB | – | `blobology`, `classify_bins` |
| Indexed host genome (`host`) | ~20 GB | ~9 GB | `read_qc` |
| KOfam profiles and thresholds (`kofam`) | ~1.5 GB | – | `annotate --kofam` |
| KEGG BRITE hierarchy (`kegg`) | ~4 MB | – | `pathways` hierarchy and Krona charts |
| Pfam-A (`pfam`) | 1.6 GB | – | `annotate --pfam` |
| eggNOG 5.0.2 (`eggnog`) | 48 GB | – | `annotate --eggnog` |
| dbCAN (`dbcan`) | ~7.4 GB | – | `annotate --cazy` |
| Bakta (`bakta`) | ~30 GB | 1.3 GB | `annotate --bakta` |
| AMRFinderPlus (`amrfinder`) | ~200 MB | – | `amr` |
| CARD (`card`) | ~4 MB | – | `amr --rgi` |
| geNomad (`genomad`) | ~1.5 GB | – | `viruses`, `amr --mobile-elements` |
| CheckV (`checkv`) | ~6.4 GB | – | `viruses` |
| SingleM metapackage (`singlem`) | 9.7 GB | – | `appraise` |
| GTDB-Tk (`gtdbtk`) | ~110 GB | – | `classify_bins --gtdbtk` |
| **everything** | **~644 GB** | **~172 GB** | |

`install-db --list` prints the current version of this table.

- Sizes are on disk after unpacking and indexing; downloads are smaller for some (GUNC 6.7 GB,
  SingleM 2.2 GB, eggNOG 12 GB). `~` marks an estimate, and the total is a lower bound. The
  `--small` column gives the download size of the smaller variant; its total is what
  `install-db --small` downloads for everything.
- `–` in the `--small` column means there is no smaller variant; `--small` installs the full
  database.
- SingleM's own downloader fetches it from Zenodo over one connection (`aria2c` is not used). It
  has taken from about 10 minutes to 2-4 hours (150-300 KB/s). GTDB-Tk needs about twice its size
  free while it unpacks.
- GUNC has one download source. If it is down, `bin_refinement --skip-gunc` runs without it and
  skips only the chimera check.

For the five databases fetched by a tool's own downloader (`amrfinder`, `genomad`, `checkv`,
`bakta`, `gtdbtk`), `install-db` runs that downloader in the tool's environment, records the path in
your config, and skips it on a re-run.

`kraken2 --bracken` needs no extra download. The pre-built KRAKEN2 standard database includes the
`databaseNmers.kmer_distrib` files Bracken reads, for 50 to 300 bp reads. For a custom Kraken2
database, build them with `bracken-build -d $KRAKEN2_DB -l <read length>` and point `BRACKEN_DB` at
the same directory.

Check that a module is ready with `doctor`:

```bash
metawrap2 doctor binning           # env, software and databases, with a status map
metawrap2 doctor --fix             # recreate exactly the environments that are broken
```

It also reports each database's content fingerprint, so a database replaced in place shows as
changed. `test` is the same command under its original name. In scripts, use
`metawrap2 check <module>`, which exits non-zero if the module is not ready.

## Lockfiles

The environment yamls set minimum versions, so a solve today can pick newer dependencies than the
ones tested. The lockfiles in `src/metawrap2/envs/locks/` list the exact linux-64 builds of an
environment that worked. For anything but development, build from them:

```bash
metawrap2 quickstart -t 8 --from-lock    # every default environment, from its lockfile
metawrap2 install-env --all --from-lock  # or with install-env
metawrap2 install-env bins --from-lock   # one environment
metawrap2 install-env bins --lock        # write a lockfile from an env you already have
conda create -n metawrap2-bins --file src/metawrap2/envs/locks/bins.linux-64.lock
```

Lockfiles exist for the `annotate`, `binqc`, `bins`, `checkm2`, `learned`, `reads` and `taxonomy`
environments. For any other environment, `--from-lock` says there is no lockfile and
solves its yaml instead.

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
machine. The old file is parsed and not executed.

### mamba on PATH at run time

Tools run as `mamba run -n metawrap2-<module> ...`, so `mamba` must be on PATH at run time.
Interactive shells get it from `~/.bashrc`; `nohup`, `cron`, `sbatch` and CI jobs usually do not.
For unattended runs:

```bash
export PATH="$HOME/miniforge3/bin:$PATH"     # wherever mamba lives
```

Without it, the run stops before any tool starts and suggests this or `use_conda_envs = false`.
