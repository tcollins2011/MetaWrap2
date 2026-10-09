# benchmarks/

Where MetaWrap2's defaults come from.

Every default in MetaWrap2 — which assembler, which binners, which completeness cut-off — is
supposed to be chosen on **accuracy per CPU-hour**, not on accuracy alone. A variant that scores
higher for four times the compute is a better flag, not a better default. This directory is what
makes that claim checkable by someone other than the person who made it.

```
metawrap2 benchmark list                            # what is declared, and today's defaults
metawrap2 benchmark fetch cami2_gi --dest data/     # prints the size, then stops
metawrap2 benchmark run cami2_gi binner_set --data data/
metawrap2 benchmark score cami2_gi binner_set
```

## What is here

| Path | What it is |
|---|---|
| `datasets.toml` | the dataset registry: URLs, pinned md5s, exact sizes, gold-standard layout |
| `configurations.toml` | the comparisons: which variants of which module, and the metrics to rank them on |
| `fixtures/tiny/` | a ~370 kB three-genome fixture, committed, so the harness is testable without a 95 GB download |
| `fixtures/make_tiny.py` | how the fixture was generated (fixed seed, so it is reproducible rather than regenerated per run) |
| `results/` | the committed results history, one TSV per dataset |

Nothing here downloads on import, on `list`, or on `run`. A fetch happens when a dataset is named
and its size is confirmed with `--yes`, and not otherwise.

## The results TSV

One long-format table per dataset:

```
dataset  comparison  variant  sample  metric  value  unit  source
```

Long rather than wide on purpose. A new metric appends rows instead of requiring every historical
file to grow a column, and cost is stored as metrics too (`wall_seconds`, `cpu_hours`,
`peak_rss_bytes`) so accuracy and cost join on the same key — which is what makes the ratio
computable at all. `sample` is empty for a figure that is already pooled across samples.

`cpu_hours` is the denominator of the ratio, not wall time: wall time rewards a tool for being given
more cores. `metawrap2 benchmark score --cost-metric wall_seconds` will rank by wall time if that is
the question, but it is not how defaults are chosen.

Check which metrics actually landed before reading a ranking. A comparison declares the metrics it
wants and scoring records whichever ones its assessors produced, so a table can be short two metrics
without saying so in the ranking itself — `binner_set` on cami2_gi declared four and delivered two,
because CheckM failed on every variant. That particular failure is fixed (a unix-socket path length
limit; see `metawrap2.checkm.MAX_TMPDIR_LEN`), but the general lesson is not: read the `N metric(s)`
line scoring prints per variant, and `--reassess` to fill a gap once the cause is dealt with.

`metawrap2 benchmark score --checkm2` scores bin quality with CheckM2 rather than CheckM1, which is
worth reaching for on a machine where CheckM1's pplacer stage cannot have its roughly 40 GB per
thread. Install it once with `metawrap2 install-env checkm2`. The AMBER metrics are
unaffected either way: they score against the gold standard rather than against marker genes.

## Adding a dataset

Add a `[[name.files]]` block per archive to `datasets.toml`. The loader will refuse it unless every
remote file declares both an `md5` and a `bytes` size, because:

- without a checksum, nobody can tell whether a 10 GB download is intact — a truncated one looks
  exactly like a good one until something reads the end of it;
- without a size, there is no disk pre-flight, and the fetch finds out it does not fit at 80%.

Sizes and checksums are **measured**, not estimated — `curl -I` for the size, and the publisher's
own manifest for the md5, pinned here rather than read at fetch time so that a silent upstream
re-publication shows up as a checksum failure rather than as quietly different results.

`reads_glob` names the read files, and naming them precisely matters more than it looks. CAMI2's
`reads/` directory holds the 10 GB interleaved FASTQ *and* a 3 MB `subsample_0.fq` debugging
extract, both under the same `sample_N` directory and therefore under the same sample label. A
directory-wide match finds both, and whichever sorts last wins - a trimmer benchmarked on 3 MB of a
10 GB sample, which produces numbers rather than an error. Two files resolving to one label is now a
hard failure, so the worst case is a message instead of a plausible-looking table.

`keep = [...]` globs decide what survives extraction. For CAMI2 that matters a lot: ~95% of each
archive is per-OTU BAM alignments no binning benchmark reads. Skipping them cannot shrink the
download — a `.tar.gz` is one gzip stream, there is nothing to seek past — but it keeps what lands
to the reads, contigs and truth tables.

## Adding a comparison

Add a `[name]` block to `configurations.toml` with at least two variants and exactly one
`default = true`. The loader enforces both: a comparison with one variant is not comparing anything,
and one with no incumbent cannot show a challenger what it has to beat.

`metrics` lists what to rank on. They have to be metrics something actually writes — the scorer
discovers CheckM1/CheckM2 reports, GUNC `maxCSS_level` tables, AMBER `results.tsv`, assembly FASTAs
and `read_qc`'s `final_pure_reads_*` by filename under each variant's output directory.

`common_args = [...]` holds the flags every variant gets: the parts of the command that are *not*
under test. `[trimmer]` uses it to switch off host removal and both FastQC reports, which are
identical work whichever trimmer ran and would otherwise bury the difference. Declaring them once
beats repeating them per variant, where a flag forgotten in the fourth variant becomes a difference
nobody meant to measure.

## Scoping a comparison that is too expensive to run whole

`--max-samples N` uses the first N samples of the dataset, and `--sample LABEL` (repeatable) names
them. A 10-sample CAMI2 co-assembly through metaSPAdes is days of wall time and hundreds of GB of
RAM, and "which assembler is more contiguous" does not need every sample to answer:

```
metawrap2 benchmark run cami2_gi assembler --max-samples 3
```

The subset applies to every variant equally, so the comparison stays internally valid. What it is
*not* is comparable with a run over a **different** subset - and a cohort module's rows carry no
sample label, so the results TSV cannot show that on its own. The samples that took part are
therefore written into `<work>/<comparison>.run.json` and printed by `score`. Check that line before
reading two runs of one comparison against each other.

A cohort assembly is fed **one concatenated pair** of read files, built once under the work
directory and reused by every variant, because `metawrap2 assembly` takes one file per mate - which
is also how `metawrap2 run` co-assembles a study, so the benchmark measures the same thing a user's
run does. It roughly doubles the reads on disk for the samples selected.

Modules that take one sample per invocation rather than a whole cohort — currently `read_qc` — are
listed in `runner.PER_SAMPLE_MODULES`. Those run once per sample, into
`work/<comparison>/<variant>/<sample>/`, and their rows carry the sample label so they sum into one
ratio. A variant that fails on any one sample records nothing at all: a total over nine samples
compared against a total over ten would look like a cheap win.

## The datasets

**`cami2_gi`** — the CAMI2 Toy Human Microbiome gastrointestinal body site: 10 samples, 95.4 GB,
CC BY 4.0, DOI [10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). Short-read
Illumina with a pooled gold-standard assembly and a contig→genome mapping, which is what AMBER needs
to score binning against known truth. Cite Meyer et al. 2022, *Nat. Methods* 19:429–440.

**`tiny`** — three synthetic genomes at 30/50/68% GC, 257 contigs, 4000 read pairs. Not a
benchmark: a fixture, so the harness's own tests run in CI in under a second. Regenerate with
`python fixtures/make_tiny.py` only if the format needs to change; it is committed so that a results
diff is signal rather than noise.

## Rebuilding the figures

```bash
python benchmarks/make_figures.py        # writes docs/figures/*.png
```

It reads two things: `results/<dataset>.tsv`, which is what `benchmark score` appends to, and the
small JSON side-inputs under `figure_inputs/` — per-bin CheckM2/GUNC joins, the metaWRAP 1 arm's
scores, the identity sweep, the end-to-end preset costs. Those are committed because a figure
nobody else can regenerate is a claim rather than a result; the path used to be absolute and point
at the machine the benchmarks were run on. `$METAWRAP2_BENCH_WORK` overrides it to point at a live
working directory instead.
