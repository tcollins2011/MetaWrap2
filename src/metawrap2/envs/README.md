# Conda environments

MetaWrap2's own code needs only Python 3.10 or later and a few Python packages, which `pip install`
provides. The third-party tools the modules run (SPAdes, CheckM, Kraken2 and so on) live in the
conda environments defined here.

Environments are grouped by workflow. Modules that run one after another and use the same tools
share one, so the common aligner and samtools stack is solved once. A tool has its own environment
when it cannot be solved together with the others, or when it is large and optional.

## Default environments

`metawrap2 quickstart` builds these eleven, plus `metawrap2-learned` from the next table.

| File | Environment | Modules | Main tools |
|---|---|---|---|
| `reads.yaml` | `metawrap2-reads` | read_qc, assembly, reassemble_bins, blobology, nonpareil | fastp, trim-galore, fastqc, bmtagger, bowtie2, spades, megahit, flye, quast, blast, minimap2, bwa, bwa-mem2, samtools, nonpareil |
| `bins.yaml` | `metawrap2-bins` | binning, coverage, dereplicate, quant_bins | metabat2, maxbin2, concoct, minimap2, bwa, bwa-mem2, samtools, skani, sylph, salmon, coverm, fairy |
| `binqc.yaml` | `metawrap2-binqc` | bin_refinement | checkm-genome (CheckM1, for `--checkm1`), pplacer, gunc, das_tool |
| `checkm2.yaml` | `metawrap2-checkm2` | bin_refinement, reassemble_bins | checkm2, the default bin scorer |
| `taxonomy.yaml` | `metawrap2-taxonomy` | classify_bins, kraken2, pathways | kraken2, bracken, krona, taxator-tk, blast |
| `annotate.yaml` | `metawrap2-annotate` | annotate, catalogue, annotate_bins | prokka, bakta, dbcan, gecco |
| `amr.yaml` | `metawrap2-amr` | amr | ncbi-amrfinderplus, hamronization, argnorm |
| `viruses.yaml` | `metawrap2-viruses` | viruses | genomad, checkv, cctyper |
| `strains.yaml` | `metawrap2-strains` | strains | instrain, minimap2, bwa-mem2, samtools |
| `comparative.yaml` | `metawrap2-comparative` | pangenome, phylogeny | ppanggolin, iqtree, fasttree, mafft |
| `singlem.yaml` | `metawrap2-singlem` | appraise | singlem |

`strains` has its own environment because inStrain requires Python 3.8. `appraise` has its own
because SingleM pins a narrow set of dependencies. CheckM2 has no common solution with CheckM1, GUNC
and DAS Tool.

## Opt-in environments

Built by name, for example `metawrap2 install-env gtdbtk`, or with all the others by
`install-env --all`.

| File | Environment | Used by | Why it is separate |
|---|---|---|---|
| `gtdbtk.yaml` | `metawrap2-gtdbtk` | `classify_bins --gtdbtk` | No common solution with taxator-tk; needs 100+ GB of reference data |
| `learned.yaml` | `metawrap2-learned` | `binning --semibin2`, `--vamb`, `--comebin` | pytorch, about 1.5 GB |
| `rgi.yaml` | `metawrap2-rgi` | `amr --rgi` | Pins older dependencies than AMRFinderPlus |
| `tracs.yaml` | `metawrap2-tracs` | `strains --backend tracs` | Pins an older scientific-Python stack than inStrain |
| `eggnog.yaml` | `metawrap2-eggnog` | `annotate --eggnog` | eggNOG-mapper imports `distutils` and needs Python below 3.12 |
| `benchmark-amber.yaml` | `metawrap2-benchmark-amber` | `benchmark score` | Benchmark only; pins matplotlib below 3.9 |

The default binners include SemiBin2, so `quickstart` also builds `metawrap2-learned`. Its pytorch
is CUDA-capable, so `binning --device cuda` needs no separate GPU environment.

Scores from different tools are not interchangeable (CheckM2 and CheckM1, for example). Use one of
them throughout a study.

## Usage

`metawrap2 install-env` reads these files; you rarely need to edit them.

```bash
metawrap2 install-env --all          # every environment, opt-in ones included
metawrap2 install-env bins binqc     # just these
metawrap2 install-env binning        # the environment a module uses, by module name
metawrap2 install-env learned        # an opt-in environment, by name
```

Or directly:

```bash
mamba env create -f src/metawrap2/envs/bins.yaml
```

Modules find their environment by name, so nothing needs to be activated first. Environment names
from MetaWrap2 2.0.x, such as `bin_refinement-checkm2` or `annotate_bins-bakta`, still work and are
mapped to the current environment.

The files are inside the package so they ship with the wheel, and `install-env` works from a
`pip install` as well as from a checkout.

## Lockfiles

The YAML files set minimum versions, so they keep solving as tools are updated, and a solve today
can differ from one a year ago. For an exact rebuild, `locks/` holds one linux-64 lockfile per
environment, the output of `conda list --explicit --md5` (every package's URL and checksum).

```bash
metawrap2 install-env --all --from-lock     # build from the committed lockfiles
metawrap2 install-env bins --lock           # write a lockfile from an environment you built
```

The Lockfiles workflow (`.github/workflows/lockfiles.yml`) re-solves every environment monthly and
pushes any changed lockfiles to the `chore/refresh-lockfiles` branch for review. Environments added
since the last refresh get a lockfile on the next run.

## No Python in the environments

The environments hold only external tools. MetaWrap2's helper scripts (`metawrap2.scripts.*`,
`metawrap2.vendor.*`) run in the Python that MetaWrap2 was installed into (see
`metawrap2/pyrun.py`), so no environment needs a copy of MetaWrap2 or its plotting stack.
