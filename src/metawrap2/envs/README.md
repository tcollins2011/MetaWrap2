# Per-module conda environments

The MetaWrap2 core is light: `python>=3.10`, `biopython`, and the plotting/table stack
(the `metawrap2` conda package installs these). The heavy third-party tools
each module drives (SPAdes, CheckM, salmon, kraken2, ...) live in their own
per-module environments here instead of one large shared environment.

Splitting them up means each tool's dependencies resolve on their own, so
adding or updating one tool no longer has to satisfy every other tool at the
same time. This is the main fix for the install/solver problems people have
hit with the older single-environment setup.

Each environment is named `metawrap2-<module>` and holds only that module's
tools. Duplication across envs (e.g. `bwa`/`samtools` in several) is intentional
so each one solves independently.

## Environments

| File | Env name | Module | Tools |
|---|---|---|---|
| `read_qc.yaml` | `metawrap2-read_qc` | read_qc | trim-galore, fastqc, bmtagger, bowtie2 |
| `assembly.yaml` | `metawrap2-assembly` | assembly | spades, megahit, quast |
| `binning.yaml` | `metawrap2-binning` | binning | metabat2, maxbin2, concoct, bowtie2, bwa, samtools |
| `bin_refinement.yaml` | `metawrap2-bin_refinement` | bin_refinement (default) | checkm-genome (CheckM1), pplacer |
| `reassemble_bins.yaml` | `metawrap2-reassemble_bins` | reassemble_bins | spades, bwa, samtools |
| `quant_bins.yaml` | `metawrap2-quant_bins` | quant_bins | salmon, bwa, samtools |
| `classify_bins.yaml` | `metawrap2-classify_bins` | classify_bins (default) | taxator-tk |
| `annotate_bins.yaml` | `metawrap2-annotate_bins` | annotate_bins (default) | prokka |
| `kraken2.yaml` | `metawrap2-kraken2` | kraken2 | kraken2, krona |
| `blobology.yaml` | `metawrap2-blobology` | blobology | blast, bowtie2, samtools, python, biopython, pandas, matplotlib-base, metawrap2 |

## Usage

Create the environment(s) for the module(s) you plan to run, e.g.:

```bash
conda env create -f src/metawrap2/envs/assembly.yaml
conda env create -f src/metawrap2/envs/binning.yaml
conda env create -f src/metawrap2/envs/bin_refinement.yaml
```

Or with mamba (faster solver):

```bash
mamba env create -f src/metawrap2/envs/assembly.yaml
```

You only need the environments for the modules you use. Activate the matching
environment before running a module, or let the module pick it up on the PATH.

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
mamba create -n metawrap2-binning --file src/metawrap2/envs/locks/binning.conda-lock.yml
```

Committing the generated lockfiles pins exact builds; regenerate them when you intend to
update a tool.

## Modern opt-in tools

Newer engines are offered as opt-in alternatives, each in its own environment and selected by
a module flag. The default environments keep the original tools and behaviour, so existing
analyses reproduce unchanged; the scores and labels these produce are *not* interchangeable
with the defaults, so pick one and stay with it within a study.

| File | Env name | Flag that selects it |
|---|---|---|
| `bin_refinement-checkm2.yaml` | `metawrap2-bin_refinement-checkm2` | `bin_refinement --checkm2` - CheckM2, avoids the CheckM1/pplacer memory cost |
| `classify_bins-gtdbtk.yaml` | `metawrap2-classify_bins-gtdbtk` | `classify_bins --gtdbtk` - GTDB-Tk classification (needs GTDB reference data) |
| `annotate_bins-bakta.yaml` | `metawrap2-annotate_bins-bakta` | `annotate_bins --bakta` - Bakta annotation (needs a Bakta database) |

Install one the same way as any other environment:
```bash
metawrap2 install-env bin_refinement-checkm2
```

## Why these environments contain no Python

Only the module's *external* tools live in each environment. MetaWrap2's own helper scripts
(`metawrap2.scripts.*`, `metawrap2.vendor.*`) deliberately run in the host interpreter that
metawrap2 itself was installed into - see `metawrap2/pyrun.py`. That keeps the environments
small and avoids installing a second, divergent copy of metawrap2 and its plotting stack into
every one of them.
