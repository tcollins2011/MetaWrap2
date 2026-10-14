# benchmarks/

The comparisons behind MetaWrap2's defaults, and their results. A default is chosen on quality
measured against a dataset with a known answer, with the compute cost reported next to it. A
variant that scores higher for four times the compute makes a good flag, not necessarily a good
default.

```
metawrap2 benchmark list                            # declared comparisons and current defaults
metawrap2 benchmark fetch cami2_gi --dest data/     # prints the download size, then stops
metawrap2 benchmark run cami2_gi binner_set --data data/
metawrap2 benchmark score cami2_gi binner_set
```

## Contents

| Path | What it is |
|---|---|
| `datasets.toml` | Dataset registry: URLs, pinned md5s, exact sizes, gold-standard layout |
| `configurations.toml` | Comparisons: which variants of which module, and the metrics to rank them on |
| `fixtures/tiny/` | A ~370 kB three-genome fixture, so the harness can be tested without a 95 GB download |
| `fixtures/make_tiny.py` | Generates the fixture (fixed seed) |
| `results/` | The results history, one TSV per dataset |
| `figure_inputs/` | Small JSON inputs for the figures (per-bin CheckM2/GUNC joins, the metaWRAP 1 arm, the identity sweep, end-to-end preset costs) |
| `make_figures.py` | Builds the figures in `docs/figures/` |

Nothing downloads on import, on `list` or on `run`. `fetch` downloads only after you confirm the size
with `--yes`.

## The results TSV

One long-format table per dataset:

```
dataset  comparison  variant  sample  metric  value  unit  source
```

A new metric adds rows, so old files never need a new column. Cost is stored as metrics too
(`wall_seconds`, `cpu_hours`, `peak_rss_bytes`), so quality and cost join on the same key. `sample`
is empty for a value already pooled across samples.

Cost is ranked by `cpu_hours`, since wall time favours whichever tool was given more cores.
`benchmark score --cost-metric wall_seconds` ranks by wall time if that is what you need.

Before reading a ranking, check which metrics were recorded. Scoring prints `N metric(s)` per
variant, and a comparison can be short of the metrics it declared if an assessor failed. Use
`--reassess` to fill the gap once the cause is fixed.

`benchmark score --checkm2` scores bin quality with CheckM2 instead of CheckM1, which helps on a
machine that cannot give CheckM1's pplacer stage about 40 GB per thread. Install it with
`metawrap2 install-env checkm2`. AMBER metrics are scored against the gold standard and are not
affected.

## Adding a dataset

Add a `[[name.files]]` block per archive to `datasets.toml`. Every remote file must declare an `md5`
and a `bytes` size, or the loader rejects it: the checksum catches a truncated download, and the
size lets the disk check run before the fetch starts. Take sizes from `curl -I` and checksums from
the publisher's manifest. Pinning the checksum here means an upstream re-publication shows up as a
checksum failure.

`reads_glob` names the read files and should be exact. CAMI2's `reads/` directory holds the 10 GB
interleaved FASTQ and a 3 MB `subsample_0.fq` in the same sample directory. If two files resolve to
one sample label, the loader stops with an error.

`keep = [...]` lists what to keep after extraction. About 95% of each CAMI2 archive is per-OTU BAM
files that no binning benchmark reads. Skipping them does not shrink the download, since a
`.tar.gz` is one stream, but it keeps them off the disk.

## Adding a comparison

Add a `[name]` block to `configurations.toml` with at least two variants, exactly one of them marked
`default = true`.

`metrics` lists what to rank on, and each must be something the scorer finds. It looks for
CheckM1/CheckM2 reports, GUNC `maxCSS_level` tables, AMBER `results.tsv`, assembly FASTAs and
`read_qc`'s `final_pure_reads_*` by file name under each variant's output directory.

`common_args = [...]` holds flags every variant gets, the parts of the command not under test.
`[trimmer]` uses it to turn off host removal and both FastQC reports, which are the same work
whichever trimmer runs.

## Running a subset

`--max-samples N` uses the first N samples, and `--sample LABEL` (repeatable) names them. A full
10-sample CAMI2 co-assembly with metaSPAdes takes days and hundreds of GB of RAM:

```
metawrap2 benchmark run cami2_gi assembler --max-samples 3
```

Every variant gets the same subset, so the comparison is fair. Results from two different subsets
are not comparable, and cohort rows carry no sample label. The samples used are therefore written to
`<work>/<comparison>.run.json` and printed by `score`; check them before comparing two runs.

A cohort assembly receives one concatenated pair of read files, built once and shared by every
variant, which is also how `metawrap2 run` co-assembles a study. It roughly doubles the reads on disk
for the selected samples.

Modules that take one sample at a time (currently `read_qc`, listed in
`runner.PER_SAMPLE_MODULES`) run once per sample into `work/<comparison>/<variant>/<sample>/`, and
their rows are summed per variant. A variant that fails on any sample records nothing, so a total
over nine samples is never ranked against one over ten. `score` skips variants whose run failed;
`--score-failed` scores them anyway.

## Datasets

`cami2_gi`: the CAMI2 Toy Human Microbiome gastrointestinal set. 10 samples, 95.4 GB, CC BY 4.0,
DOI [10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). Short Illumina reads with a
pooled gold-standard assembly and a contig-to-genome mapping, which AMBER uses to score binning.
Cite Meyer et al. 2022, *Nat. Methods* 19:429–440.

`tiny`: three synthetic genomes at 30, 50 and 68% GC, 257 contigs, 4000 read pairs. A test fixture
rather than a benchmark: it lets the harness's own tests run in CI in under a second. Regenerate it
with `python fixtures/make_tiny.py` only if the format changes.

## Rebuilding the figures

```bash
python benchmarks/make_figures.py        # writes docs/figures/*.png
```

It reads `results/<dataset>.tsv` and `figure_inputs/`. Set `$METAWRAP2_BENCH_WORK` to read the
figure inputs from a live working directory instead.

The release-note figures are built by `make_release_figures.py`, which reads the run outputs
themselves (`$MW2_RUNS`). Those are too large to commit, so these figures can only be rebuilt on a
machine that has the runs.
