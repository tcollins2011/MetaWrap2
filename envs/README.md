# Per-module conda environments

The MetaWrap2 core is lightweight: it only needs `python>=3.8` and `biopython`
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
conda env create -f envs/assembly.yaml
conda env create -f envs/binning.yaml
conda env create -f envs/bin_refinement.yaml
```

Or with mamba (faster solver):

```bash
mamba env create -f envs/assembly.yaml
```

You only need the environments for the modules you use. Activate the matching
environment before running a module, or let the module pick it up on the PATH.

## Reproducible environments (lockfiles)

The YAMLs here use version lower bounds so they keep solving as upstream tools update. For
bit-for-bit reproducible environments, generate lockfiles on your target platform with
[conda-lock](https://github.com/conda/conda-lock), e.g.:

```bash
conda-lock -f envs/binning.yaml -p linux-64 --lockfile envs/locks/binning.conda-lock.yml
mamba create -n metawrap2-binning --file envs/locks/binning.conda-lock.yml
```

Committing the generated lockfiles pins exact builds; regenerate them when you intend to
update a tool.

## Modern opt-in tools (planned)

Environment files are provided for newer engines we plan to offer as opt-in
alternatives. The module flags that select them are **not wired in 2.0.0** yet -
the default environments keep the original tools and behavior. These are kept
here so the environments are ready when the flags land in a future release:

| File | Env name | Planned use |
|---|---|---|
| `bin_refinement-checkm2.yaml` | `metawrap2-bin_refinement-checkm2` | CheckM2 (avoids the CheckM1/pplacer memory cost) |
| `classify_bins-gtdbtk.yaml` | `metawrap2-classify_bins-gtdbtk` | GTDB-Tk classification |
| `annotate_bins-bakta.yaml` | `metawrap2-annotate_bins-bakta` | Bakta annotation |
