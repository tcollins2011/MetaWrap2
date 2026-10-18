# Tutorial

A first end-to-end run on the benchmark data. Because the expected result is known, a low count
points to the installation and not the sample.

1. A smoke test on the committed fixture. Nothing is downloaded; it checks that the install works.
2. Three CAMI2 gastrointestinal samples, from download to MAGs, abundance, taxonomy and KEGG
   modules.

The CAMI2 data is simulated from known genomes, with a published contig-to-genome truth table. The
expected result is committed here ([Checking your results](#checking-your-results)). All benchmark
and tutorial datasets in this project come from CAMI for this reason.

Run commands from a clone of the repository. Install first: [Installation](installation.md).

## Part 1: smoke test on the shipped fixture

`benchmarks/fixtures/tiny` is three synthetic genomes split into a few hundred contigs, with one
pair of read files (4,000 read pairs) and a gold standard. It lets the pipeline and benchmark harness run without the 95 GB
download of Part 2. It is too small to measure binner accuracy.

```bash
metawrap2 doctor binning                      # env, tools and databases for one module
metawrap2 binning benchmarks/fixtures/tiny/reads/sample_1_1.fq.gz \
    benchmarks/fixtures/tiny/reads/sample_1_2.fq.gz \
    -a benchmarks/fixtures/tiny/assembly.fa -o /tmp/tiny_binning -t 4 --metabat2
ls /tmp/tiny_binning/metabat2_bins /tmp/tiny_binning/metabat2_extra
```

Expected result: zero bins. The fixture's 257 contigs are too few for metaBAT2 to call a cluster,
so it writes only bookkeeping files, which MetaWrap2 moves to `metabat2_extra/` (`bin.unbinned.fa`
holds the unbinned contigs). `metabat2_bins/` is empty and the run exits 0 after reporting
`0 bins (0 bases in total) formed`.

The run covers mapping, the depth table, metaBAT2 in its environment, filtering of non-genome files,
and the provenance record. A failure in any of them gives a non-zero exit and names the step. It
takes about 25 seconds at `-t 4`.

`--metabat2` runs one binner and does not need pytorch, so `metawrap2-learned` does not have to be
built yet.

The same fixture runs the benchmark harness, if you plan to run benchmarks:

```bash
metawrap2 benchmark run   tiny trimmer        # read_qc with each trimmer, on the fixture
metawrap2 benchmark score tiny trimmer        # read the counts into the results TSV and rank them
```

Cost: on the fixture, the metaBAT2 binning step took 3.7 and 4.9 seconds of wall time over two
runs, at 0.003 CPU-hours. `doctor` is not timed. It starts every tool and runs the Python unit tests
(`--no-unit-tests` skips those).

## Part 2: three CAMI2 samples

### What you need installed

| | Why | Install |
|---|---|---|
| `metawrap2-reads`, `metawrap2-bins`, `metawrap2-learned` | read_qc, assembly, binning (SemiBin2 is in `learned`) | `metawrap2 quickstart -t 8` |
| `metawrap2-binqc` | bin_refinement: CheckM, GUNC | as above |
| `metawrap2-taxonomy`, `metawrap2-annotate` | classify_bins, annotate, pathways | as above |
| `checkm` 1.4 GB, `gunc` 6.7 GB (~13 GB installed) | bin_refinement scores bins and checks chimerism | `metawrap2 install-db checkm gunc -t 2` |
| `kofam` 1.5 GB, `kegg` 4 MB | `annotate --kofam`, `pathways` | `metawrap2 install-db kofam kegg -t 2` |
| `blast` ~2 GB with `--small`, `taxdump` 550 MB | `classify_bins` without GTDB-Tk | `metawrap2 install-db blast taxdump --small -t 2` |

No host genome index is needed. The simulated community has no human reads, so
`read_qc --skip-bmtagger` below avoids a 20 GB download. On real data, drop that flag and install
`host`.

### Step 0: fetch the data

```bash
metawrap2 benchmark fetch cami2_gi            # prints the plan and the size, downloads nothing
metawrap2 benchmark fetch cami2_gi --yes      # 95.4 GB, checksum-verified against datasets.toml
```

- 95.4 GB over 13 files: ten sample archives of about 9.5 GB each, the pooled gold-standard
  assembly and two contig-to-genome truth tables. `fetch` downloads the whole dataset; this
  walkthrough uses three of the ten samples. Each archive is deleted after extraction, so peak disk
  use is about the download size.
- `fetch` does not start if the download would exceed the `[limits]` disk ceiling.
- Licence: CC BY 4.0. Cite Meyer F, et al. *Critical Assessment of Metagenome Interpretation - the
  second round of challenges.* Nat Methods 19:429-440 (2022), DOI 10.4126/FRL01-006425518. `fetch`
  prints both, and both are recorded in `benchmarks/datasets.toml`.

The samples land under `benchmarks/data/cami2_gi/extracted/`, one directory per sample, named with
the publisher's timestamp prefix (`2017.12.04_18.45.54_sample_0`). Reads are one interleaved FASTQ
per sample:

```bash
ls -d benchmarks/data/cami2_gi/extracted/*sample_0
ls benchmarks/data/cami2_gi/extracted/*sample_0/reads/anonymous_reads.fq
```

This walkthrough uses samples 0, 1 and 10, on which the committed results were measured
(`docs/benchmarks.md`, Methods).

```bash
export CAMI=benchmarks/data/cami2_gi/extracted   # every command below reads from here
```

### Step 1: read_qc

```bash
for S in 0 1 10; do
    metawrap2 read_qc --interleaved -1 $CAMI/*_sample_$S/reads/anonymous_reads.fq \
        -o READ_QC/sample_$S -t 24 --skip-bmtagger
done
```

`--interleaved` splits the mates and processes the sample as paired. Later steps read
`READ_QC/sample_$S/final_pure_reads_1.fastq` and `_2.fastq`. The names are the same with or without
host removal.

Cost: 58 to 62 seconds of wall time and 0.066 CPU-hours per sample for fastp, the default trimmer,
measured over all ten samples of this dataset.

### Step 2: assembly

```bash
for S in 0 1 10; do
    metawrap2 assembly -1 READ_QC/sample_$S/final_pure_reads_1.fastq \
        -2 READ_QC/sample_$S/final_pure_reads_2.fastq \
        -o ASSEMBLY/sample_$S -t 24 -m 64
done
```

MEGAHIT is the default and produced the committed results. `--metaspades` gives more contiguous
assemblies but was measured at 244.6 CPU-hours against MEGAHIT's 77.3 on the three-sample
co-assembly (`cami2_gi.tsv`, `assembler/*`), and at about 140 GB of resident memory.

Cost: assembly is the most expensive step, and its cost does not depend on later choices.
Assembling these three samples separately is 83 CPU-hours, about 65% of a default run
([Benchmarks](benchmarks.md) §2). In the co-assembly arm of the assembler comparison, MEGAHIT took
77.3 CPU-hours and 62 minutes of wall time (`cami2_gi.tsv`).

### Step 3: binning

```bash
for S in 0 1 10; do
    metawrap2 binning READ_QC/sample_$S/final_pure_reads_1.fastq \
        READ_QC/sample_$S/final_pure_reads_2.fastq \
        -a ASSEMBLY/sample_$S/final_assembly.fasta -o INITIAL_BINNING/sample_$S \
        -t 24 --preset default
done
```

Reads are positional and go before the flags, so no read path follows a flag that takes a value.
`--preset default` runs metaBAT2 and SemiBin2, giving `INITIAL_BINNING/sample_$S/metabat2_bins` and
`semibin2_bins`. Each binner gets one coverage column per read set passed. For a real cohort,
`metawrap2 coverage` builds one depth table across all samples with N alignments instead of N²
([Usage](usage.md)).

Cost: binning and refinement together are 43.8 CPU-hours plus 0.6 GPU-hours for `--preset default`
over these three samples ([Benchmarks](benchmarks.md) §2). The committed table does not split the
two. metaBAT2 takes about 1.3 CPU-hours of this, since it reuses SemiBin2's depth table (§3).

### Step 4: bin_refinement

```bash
for S in 0 1 10; do
    metawrap2 bin_refinement -o BIN_REFINEMENT/sample_$S -t 24 -c 70 -x 10 \
        -A INITIAL_BINNING/sample_$S/metabat2_bins \
        -B INITIAL_BINNING/sample_$S/semibin2_bins
done
```

This step determines the MAG count. It hybridises the two bin sets, uses CheckM to choose which
version of each bin to keep, then runs GUNC on the kept bins to flag chimeras of two genomes, which
CheckM scores as complete and clean. The result is
`BIN_REFINEMENT/sample_$S/metawrap_70_10_bins` with `metawrap_70_10_bins.stats.tsv` beside it: one
row per bin with completeness, contamination, GC, lineage, N50, size and the two GUNC columns.

`--checkm2` scores with CheckM2 instead of CheckM1. It needs much less memory (CheckM1's pplacer
stage needs roughly 40 GB per thread). The committed numbers were scored with CheckM2, so use it
when comparing with them. It is opt-in, with its own environment and database:
`metawrap2 install-env checkm2 && metawrap2 install-db checkm2`.

Cost: included in the 43.8 CPU-hours above. The consolidator comparison measured 4.52 CPU-hours to
refine one sample's bins with CheckM
([Benchmarks](benchmarks.md) §5).

### Step 5: quant_bins

```bash
for S in 0 1 10; do
    metawrap2 quant_bins READ_QC/sample_$S/final_pure_reads_1.fastq \
        READ_QC/sample_$S/final_pure_reads_2.fastq \
        -b BIN_REFINEMENT/sample_$S/metawrap_70_10_bins \
        -a ASSEMBLY/sample_$S/final_assembly.fasta -o QUANT_BINS/sample_$S -t 24
done
```

Output: `QUANT_BINS/sample_$S/bin_abundance_table.tsv`, plus a clustered heatmap when there is more
than one sample. The default quantifier, sylph, uses k-mer sketches and does no alignment. Numbers
differ between quantifiers. A sylph zero means "not detected above `--min-ani`"
([Benchmarks](benchmarks.md) §6).

Cost: 65.5 seconds of wall time, 0.126 CPU-hours and 1.31 GB peak memory for sylph on this dataset
(`cami2_gi.tsv`, `quantifier/sylph`). salmon took 1.078 CPU-hours and CoverM 2.368 for the same
work.

### Step 6: classify_bins

```bash
for S in 0 1 10; do
    metawrap2 classify_bins -b BIN_REFINEMENT/sample_$S/metawrap_70_10_bins \
        -o CLASSIFY/sample_$S -t 24
done
```

Output: `CLASSIFY/sample_$S/bin_taxonomy.tsv`, one consensus taxonomy per bin. With the `--small`
BLAST database (a representative-prokaryote set instead of nt) the taxonomy is coarse. `--gtdbtk`
is preferred: taxator-tk has had no upstream release since 2015, and the MAG literature uses GTDB.
It is not the default because it needs `metawrap2 install-env gtdbtk` plus
`metawrap2 install-db gtdbtk`, which fetches about 110 GB of reference data and records it as
`GTDBTK_DATA_PATH`.

Cost: not measured. No comparison in this repository covers `classify_bins`.

### Step 7: function (new in 2.1.2)

Gene calling and KEGG orthologs, then module completeness, per sample on its refined bins:

```bash
for S in 0 1 10; do
    metawrap2 annotate -o ANNOTATION/sample_$S -b BIN_REFINEMENT/sample_$S/metawrap_70_10_bins \
        -t 24 --kofam
    metawrap2 pathways -o PATHWAYS/sample_$S -k ANNOTATION/sample_$S/kos.tsv \
        -a QUANT_BINS/sample_$S/bin_abundance_table.tsv
done
```

`--kofam` reads the KOfam path `install-db` wrote into your config; `--kofam-db DIR` overrides it.
Outputs:

- `ANNOTATION/sample_$S/genes.tsv`: every called gene. `kos.tsv`: the KEGG ortholog per gene, with
  the evidence for it. On one CAMI2 sample's refined bins, `annotate` called 14,171 genes and 5,827
  took a curated KO, 2,008 of them distinct ([2.1.2 release notes](releases/2.1.2.md) §1).
- `PATHWAYS/sample_$S/modules.tsv`: KEGG module completeness per genome, under both the `pathwise`
  and `stepwise` definitions, with seven flags for known cases where completeness is misleading.
  `-a` weights the per-genome rows by abundance into a community profile, so `quant_bins` must run
  first.

Incomplete genomes lose a larger share of module completeness. In the degradation experiment
([2.1.2 release notes](releases/2.1.2.md) §1.2), a 70%-complete MAG, which passes the high-quality
tier, had lost 59% of the modules the intact genome completes. `modules.tsv` includes each genome's
completeness for this reason.

Cost: `annotate --kofam` took 76 seconds per bin on one CAMI2 sample's bins (§1).
`pathways` is not measured; it only computes over a KO table.

### All three samples with `cohort`

`metawrap2 cohort` runs all of the above from a sample sheet, which also serves as the record of
the run. Fill in the timestamped directory names from `ls`:

```toml
# samples.toml - metawrap2 cohort samples.toml -o study/ -t 32 -j 3
[settings]
coassemble = false            # per-sample assembly; `cohort` overrides this, `run` honours it

[modules]
read_qc = "--skip-bmtagger"   # simulated community: no host reads to remove
binning = "--preset default"
bin_refinement = "-c 70 -x 10"

[[samples]]
name = "sample_0"
r1 = "benchmarks/data/cami2_gi/extracted/<timestamp>_sample_0/reads/anonymous_reads.fq"
layout = "interleaved"

[[samples]]
name = "sample_1"
r1 = "benchmarks/data/cami2_gi/extracted/<timestamp>_sample_1/reads/anonymous_reads.fq"
layout = "interleaved"

[[samples]]
name = "sample_10"
r1 = "benchmarks/data/cami2_gi/extracted/<timestamp>_sample_10/reads/anonymous_reads.fq"
layout = "interleaved"
```

```bash
metawrap2 --dry-run cohort samples.toml -o study/ -t 32 -j 3   # print every module command
metawrap2 cohort samples.toml -o study/ -t 32 -j 3
```

A cohort run adds two steps. `coverage` gives every binner each contig's depth in every sample from
N alignments instead of N². `dereplicate` collapses the per-sample bin sets into
`study/MAG_CATALOGUE/dereplicated_genomes`, with one genome per organism. `--dry-run` prints the
plan and runs nothing. `metawrap2 status study/` shows what a `--resume` would redo.

## Checking your results

With the default settings, these three CAMI2 gastrointestinal samples give 53 clean high-quality
MAGs of 56 high-quality ones, end to end, for 43.8 CPU-hours of binning and refinement
([Benchmarks](benchmarks.md) §2, `--preset default` row). Count yours as the benchmark does, by
joining the CheckM2 and GUNC reports per bin. Multiplying a high-quality count by a GUNC pass
fraction undercounted one arm by 18 genomes.

```bash
# clean high-quality = completeness >=70, contamination <=10, and not chimeric by GUNC
for S in 0 1 10; do
    awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["completeness"] >= 70 && $col["contamination"] <= 10 && $col["gunc_pass"] == "yes"' \
        BIN_REFINEMENT/sample_$S/metawrap_70_10_bins.stats.tsv | wc -l
done
```

The loop prints one count per sample. The sum is comparable with 53 if:

- The samples are 0, 1 and 10 and the preset is `--preset default`. Other presets give 36
  (`very-fast`), 41 (`fast`) and 60 (`slow`) (§2).
- Bins are scored with CheckM2. The committed counts come from CheckM2; `bin_refinement` uses CheckM1
  unless you pass `--checkm2`. The two report different completeness for the same genome
  ([Benchmarks](benchmarks.md), Methods).
- The count is summed over the three samples. The same organism binned in two samples is two bins
  and one genome; in a cohort run, `dereplicate` handles this.

Then:

- Within a few MAGs of 53: normal. Bin counts vary with thread count, tool versions and which reads
  pass trimming.
- Far below, or zero: usually an environment that did not build (`metawrap2 doctor`), a missing
  GUNC or CheckM database (`bin_refinement` reports this before CheckM starts), or binning raw
  instead of QC'd reads.
- Empty `metawrap_70_10_bins`: check `metawrap2.log` and `run_commands.txt` in the output
  directory. Every run records the commands it ran and the tool versions.

## Where to go next

- [Usage](usage.md): sample sheets, `run` against `cohort`, resource budgets, and the other
  modules: viruses, AMR, strains, pangenomes, differential abundance.
- [Benchmarks](benchmarks.md): how every default here was chosen, and what is still open.
- [Modules](modules.md): what each module does, option by option.
- `metawrap2 benchmark list`: the datasets and comparisons you can re-run yourself, with the current
  defaults marked.
