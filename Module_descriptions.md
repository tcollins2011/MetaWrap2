# Detailed descriptions of each MetaWrap2 module

Every module is a standalone Python program under `metawrap2 <module>`, runs in its own conda
environment (`metawrap2-<module>`), and prints its own help with `metawrap2 <module> -h`. The exact
commands each module runs are visible at the top of its source file in
`src/metawrap2/modules/<module>.py` (the `COMMANDS` and `CONSTANTS` blocks); edit those, or override
surfaced values in `metawrap2.toml`, to change behavior.

## read_qc
The `read_qc` module pre-processes raw Illumina paired-end reads in preparation for assembly and
alignment. Reads are quality- and adapter-trimmed with Trim Galore (its default settings), leaving
only high-quality sequence. The reads are then aligned to a host genome (e.g. human) with bmtagger,
and any host reads are removed to eliminate host contamination; read pairs where only one mate maps
to the host are removed as well. FastQC is used to generate quality reports of the input and final
read sets so you can assess the improvement. Each step is optional: `--skip-trimming`,
`--skip-bmtagger`, `--skip-pre-qc-report`, and `--skip-post-qc-report` turn off the corresponding
stage. The host index prefix defaults to `hg38` and can be changed with `-x/--host`. Gzipped inputs
(`.gz`) are decompressed automatically. The cleaned reads are written to `final_pure_reads_1.fastq`
and `final_pure_reads_2.fastq`.

## assembly
The `assembly` module assembles a set of metagenomic reads with MEGAHIT and/or metaSPAdes. With no
assembler flag (or with `--megahit`) it uses MEGAHIT, which is memory-efficient, fast, and scales
well with large datasets. With `--metaspades` it runs pure metaSPAdes, which usually produces a
superior assembly. Passing both `--metaspades --megahit` runs the hybrid pipeline: the reads are
first assembled with metaSPAdes, the reads that do not map back to the long scaffolds are pulled out,
and those leftovers are assembled with MEGAHIT (better on low-coverage data); the two assemblies are
then combined. In all cases the contigs are sorted by length and renamed to resemble SPAdes naming
(contig ID, length, and coverage), short scaffolds are discarded (`-l/--min-len`, default 1000 bp),
and an assembly report is generated with QUAST. Memory is set with `-m/--memory` (GB, default 24).

## kraken2
The `kraken2` module takes any number of fastq and/or fasta files, classifies the contained
sequences with KRAKEN2, and reports the taxonomy distribution in an interactive kronagram made with
KronaTools. Read files are recognized by the `*_1.fastq`/`*_2.fastq` convention and FASTA files by
their `.fa`/`.fasta` extension (`.gz` accepted). Reads can be randomly subsampled with `-s/--subsample`
to speed up large runs, and `--no-preload` runs KRAKEN2 in memory-mapped mode (slower but lower RAM).
For each input a `.krak2` file (raw KRAKEN2 output) is produced, translated into a `.kraken2` lineage
file and a `.krona` summary, and all samples are combined into a single `kronagram.html`.

> Note: the old `kraken` (KRAKEN1) module has been removed. Use `kraken2` with a KRAKEN2 database.

## binning
The `binning` module is a convenient wrapper around three metagenomic binners: metaBAT2, MaxBin2, and
CONCOCT. First the assembly is indexed with bwa and the paired-end reads from any number of samples are
aligned to it; the alignments are sorted and compressed with samtools. metaBAT2's
`jgi_summarize_bam_contig_depths` generates the contig abundance table, which is then converted into
the input format each binner expects. You choose which binners to run with `--metabat2`, `--maxbin2`,
and/or `--concoct` (at least one is required). Bins are written to per-binner folders
(`metabat2_bins`, `maxbin2_bins`, `concoct_bins`) with formatted FASTA files. `--universal` switches
MaxBin2 to universal (rather than bacterial) markers for better Archaea binning. Reads default to
paired-end; use `--single-end` or `--interleaved` for other layouts. Optionally, `--run-checkm`
immediately runs CheckM (via the `bin_refinement` environment) to estimate the completion and
contamination of each bin.

> Note: the `--metabat1` binner has been removed. Use `--metabat2`.

## bin_refinement
The `bin_refinement` module uses a hybrid approach to take one to three bin sets (`-A`, `-B`, `-C`)
obtained with different software (or the same software with different parameters) and produces a
consolidated, improved bin set. First, Binning_refiner creates hybridized bins from every possible
combination of sets: for three sets A, B, and C it produces AB, BC, AC, and ABC. CheckM then evaluates
the completion and contamination of the bins in each candidate set (up to 3 originals + 4 hybrids).
The sets are iteratively compared, and each pair is consolidated: the same bin is identified across two
sets based on a minimum genome-length overlap, and the better bin is chosen by the score
`S = Completion - 5 x Contamination`. After all sets are incorporated, a dereplication step handles
contigs that appear in more than one bin (by default they are kept only in the best bin;
`--keep-ambiguous` and `--remove-ambiguous` change this). CheckM is run on the final set and a report
of completion, contamination, and other statistics is generated for each bin, along with completion and
contamination rank plots comparing the refined set to the originals. `-c/--completeness` (default 70)
and `-x/--contamination` (default 10) tune the target bin quality, and drive the `metawrap_<c>_<x>_bins`
output naming. `--quick` adds CheckM's reduced-tree mode to save memory, and the `--skip-refinement`,
`--skip-checkm`, and `--skip-consolidation` flags disable individual stages.

## reassemble_bins
The `reassemble_bins` module improves a set of bins by finding the reads that map to them and
reassembling them. bwa indexes the entire assembly and aligns the reads (`-1`/`-2`) back to it; reads
mapping to contigs belonging to a bin are collected (if only one mate maps, the pair is still recruited).
For each bin two read sets are stored: a strict set (at most `--strict-cut-off` mismatches, default 2)
and a permissive set (at most `--permissive-cut-off` mismatches, default 5). Each set is reassembled
with SPAdes, and short contigs are removed (`-l/--min-len`, default 500). CheckM evaluates the three
versions of each bin (original, strict, permissive), and the best version by the score
`S = Completion - 5 x Contamination` is added to the final set. The final set is re-evaluated with CheckM,
summary statistics are generated, and completion/contamination rank plots are drawn. `-c/--completeness`
(default 70) and `-x/--contamination` (default 10) set the desired bin quality. Long reads can be added
for reassembly with `--nanopore <reads>`, `--parallel` runs the per-bin SPAdes jobs concurrently
(1 thread each), and `--skip-checkm` disables CheckM.

## quant_bins
The `quant_bins` module quickly estimates the abundance of bins across any number of samples. Salmon
indexes the metagenomic assembly and aligns the reads from each sample back to it, producing per-contig
coverage estimates. The abundance of each bin in each sample is computed as the length-weighted average
of its contigs' abundances. A final bin abundance table is produced and a clustered heatmap is drawn.
Provide the bins with `-b`, the paired reads as positional arguments, and (strongly recommended) the
entire non-binned assembly with `-a` so abundances are estimated in the context of the whole community.

## blobology
The `blobology` module creates blobplots (GC vs. coverage of all contigs) of a metagenomic assembly and
annotates them with taxonomic and/or bin information. Contigs are taxonomically classified with MEGABLAST
against the NCBI nt database, the reads from any number of samples are mapped back with bowtie2 to
estimate coverage, and the two are combined into a blobplot table with the GC, per-sample coverage, and
taxonomy of each contig. If a folder of bins is supplied with `--bins`, contigs are also annotated with
the bin they belong to (contig names must match the assembly). Blobplots are drawn at the super-kingdom,
phylum, and order levels, plus a bin-membership plot. The assembly can be randomly subsampled to a set
number of contigs with `--subsample INT`.

## classify_bins
The `classify_bins` module is a conservative but accurate way to assign taxonomy to a set of bins.
The contigs from all bins are aligned to the NCBI nt database with MEGABLAST, and taxator-tk (megan-lca)
estimates the most likely taxonomy of each contig. The overall taxonomy of each bin is then derived from
its per-contig predictions: contig taxonomies are placed on a phylogenetic tree weighted by contig
length, and the tree is traversed from the root, descending a rank only if the next branch carries more
than half the weight of the current one (a minimum-confidence threshold). Once no further rank can be
resolved, the final taxonomy of the bin is reported.

## annotate_bins
The `annotate_bins` module functionally annotates a set of bins with PROKKA. PROKKA itself drives a
variety of tools (BLAST, HMMER, Aragorn, Prodigal, tbl2asn, Infernal). The annotation is parallelized
across bins and threads. For each bin the module returns an annotation file in GFF format and two FASTA
files with the untranslated and translated genes.
