# benchmarks/

The comparisons behind MetaWrap2's defaults, and their results. Defaults are chosen on quality
measured against a dataset with a known answer, with the compute cost reported alongside. A variant
that scores higher at four times the compute may be offered as a flag without becoming the
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
| `results/` | Written by `benchmark score`: one TSV per dataset (not tracked) |

Nothing is downloaded on import, `list` or `run`. `fetch` downloads only after you confirm the size
with `--yes`.

## Results TSV

One long-format table per dataset:

```
dataset  comparison  variant  sample  metric  value  unit  source
```

A new metric adds rows, so old files do not need a new column. Cost is stored as metrics too
(`wall_seconds`, `cpu_hours`, `peak_rss_bytes`), so quality and cost join on the same key. `sample`
is empty for a value pooled across samples.

Cost is ranked by `cpu_hours`, since wall time favours whichever tool was given more cores.
`benchmark score --cost-metric wall_seconds` ranks by wall time.

Scoring prints `N metric(s)` per variant. If an assessor failed, a comparison has fewer metrics than
it declared; check the count before reading a ranking. After fixing the cause, `--reassess` fills
in the missing metrics.

`benchmark score` scores bin quality with CheckM2 (`metawrap2 install-env checkm2`). `--checkm1`
(deprecated) scores with CheckM1 on its full genome tree, whose pplacer stage needs about 40 GB per
thread. AMBER metrics are scored against the gold standard and are not
affected.

## Adding a dataset

Add a `[[name.files]]` block per archive to `datasets.toml`. Every remote file must declare an `md5`
and a `bytes` size, or the loader rejects it. The checksum detects a truncated download, and the
size lets the disk check run before the fetch starts. Take sizes from `curl -I` and checksums from
the publisher's manifest. If the publisher re-publishes a file, the pinned checksum fails.

`reads_glob` names the read files and should match exactly. CAMI2's `reads/` directory holds the
10 GB interleaved FASTQ and a 3 MB `subsample_0.fq` in the same sample directory. If two files
resolve to one sample label, the loader stops with an error.

`keep = [...]` lists what to keep after extraction. About 95% of each CAMI2 archive is per-OTU BAM
files that no binning benchmark reads. The download size is the same (a `.tar.gz` is one stream),
but the BAM files are not written to disk.

## Adding a comparison

Add a `[name]` block to `configurations.toml` with at least two variants, exactly one of them marked
`default = true`.

`metrics` lists what to rank on. Each must be a metric the scorer can find. It looks for
CheckM1/CheckM2 reports, GUNC `maxCSS_level` tables, AMBER `results.tsv`, assembly FASTAs and
`read_qc`'s `final_pure_reads_*` by file name under each variant's output directory.

`common_args = [...]` holds flags every variant gets, for the parts of the command not under test.
`[trimmer]` uses it to turn off host removal and both FastQC reports, which do not depend on the
trimmer.

## Running a subset

`--max-samples N` uses the first N samples, and `--sample LABEL` (repeatable) names them. A full
10-sample CAMI2 co-assembly with metaSPAdes takes days and hundreds of GB of RAM:

```
metawrap2 benchmark run cami2_gi assembler --max-samples 3
```

Every variant gets the same subset. Results from different subsets are not comparable, and cohort
rows have no sample label, so the samples used are written to `<work>/<comparison>.run.json` and
printed by `score`. Check them before comparing two runs.

A cohort assembly gets one concatenated pair of read files, built once and shared by every variant,
the same way `metawrap2 run` co-assembles a study. This roughly doubles the reads on disk for the
selected samples.

Modules that take one sample at a time (currently `read_qc`, listed in
`runner.PER_SAMPLE_MODULES`) run once per sample into `work/<comparison>/<variant>/<sample>/`, and
their rows are summed per variant. A variant that fails on any sample records nothing, so a total
over nine samples cannot be ranked against one over ten. `score` skips variants whose run failed;
`--score-failed` scores them anyway.

## Datasets

`cami2_gi`: the CAMI2 Toy Human Microbiome gastrointestinal set. 10 samples, 95.4 GB, CC BY 4.0,
DOI [10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). Short Illumina reads with a
pooled gold-standard assembly and a contig-to-genome mapping, which AMBER uses to score binning.
Cite Meyer et al. 2022, *Nat. Methods* 19:429–440.

`tiny`: three synthetic genomes at 30, 50 and 68% GC, 257 contigs, 4000 read pairs. It is a test
fixture for the harness's own tests, which run in CI in under a second. Regenerate it
with `python fixtures/make_tiny.py` only if the format changes.
