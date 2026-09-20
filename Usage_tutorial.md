# A guide to analyzing metagenomic data with MetaWrap2

Note: This pipeline is only a guide. None of MetaWrap2's modules depend on each other, so if you want
to do certain steps with other software, you are free to do so. For example, if you want to try the
`reassemble_bins` module on your own bins, you do not need to run the other modules just to get to that
point. Or if you want to use a different assembler than metaSPAdes or MEGAHIT, you can do so and then
proceed through the rest of the pipeline with your own assembly.

## Before you start

Three commands get you from a fresh checkout to a working install:

```bash
pip install -e .                       # the MetaWrap2 core
metawrap2 install-env --all -t 8       # one small conda env per module, 8 at a time
metawrap2 install-db --all -t 4        # download, index, and configure every database
```

Both installers are parallel (`-t`) and resumable: they skip anything already in place, so a
re-run after an interruption picks up where it stopped. Add `--small` to `install-db` to get
capped-but-real database builds, which is much faster and enough to verify the pipeline works
(use full-size databases for real analyses). You only need the envs and databases for the
modules you actually run, so `metawrap2 install-env binning bin_refinement` is fine too.

Then confirm everything is ready:
```bash
metawrap2 check              # conda envs + databases for every module
metawrap2 test               # go into each env and probe that the tools actually run
metawrap2 config show        # which config file is in effect, and which paths are missing
```

## The short version: run the whole thing with one command

Everything below walks through the modules one at a time, which is the right way to learn the
pipeline and the right way to work when you want control over each step. If you just want the
standard genome-recovery pipeline run over a set of samples, describe them once and let
`metawrap2 run` do the plumbing:

```bash
metawrap2 run --example > samples.toml     # then edit it to point at your reads
metawrap2 run samples.toml -o study/ -t 24 -j 2
```

It creates the intermediate layout, runs the steps in order, processes independent samples
concurrently (`-j`, dividing the thread budget between them), and does the between-step work
that the rest of this guide asks you to do by hand - collecting the QC'd reads, concatenating
them for the co-assembly, handing each module the right bin directory.

```bash
metawrap2 run samples.toml -o study/ --dry-run   # print every command it would run
metawrap2 run samples.toml -o study/ --resume    # carry on after a failure or interruption
metawrap2 run reads_directory/ -o study/         # no sheet: pair up the reads in a directory
```

Each step is an ordinary `metawrap2 <module>` invocation, logged to `study/LOGS/<step>.log`, so
nothing is hidden and anything the driver does you can also do yourself.

Every module is invoked as `metawrap2 <module> ...` and prints its own help with
`metawrap2 <module> -h`.

**Read file names.** Paired files are recognised in any of the usual conventions -
`sample_1/sample_2`, `sample_R1/sample_R2`, `sample_R1_001/sample_R2_001`, `sample.1/sample.2` -
so there is no need to rename anything first. `.gz` and `.bz2` are read directly. Single-end
and interleaved data are supported with `--single-end` / `--interleaved` (interleaved input is
split into mates first and then follows the ordinary paired path).

**Threads.** Set a default once in `~/.metawrap2/config.toml` under `[settings] threads = 24`
and every module uses it; `-t N` on a module overrides it, and `-t all` uses every core. The
commands below pass `-t 24` explicitly so they work whatever your config says.

**Resuming a run.** Every module checkpoints its stages. If a run dies partway through, re-run
the same command with `--resume` and it skips the stages that already finished. Use `--force`
to overwrite an output directory instead, and `--dry-run` to print the commands a module
*would* run without running them.

**Provenance.** Every output directory gets `run_commands.txt` (the exact commands, re-runnable),
`run_config.json` (the full machine-readable record), and `provenance.txt` (a readable history
that also includes the steps of any MetaWrap2 input), so you can trace a final bin back to QC.

## Step 0: Download sample metagenomic data from the metaHIT gut survey (or use your own demultiplexed, paired-end Illumina reads).

Download data from 3 samples:
```
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_2.fastq.gz
```

MetaWrap2 accepts gzipped inputs directly, so you can leave the reads compressed. If you prefer to work
with uncompressed files, unzip them:
```
gunzip *.gz
```


Place the raw sequencing reads into a new folder
```
mkdir RAW_READS
mv *fastq RAW_READS

ls RAW_READS
ERR011347_1.fastq
ERR011347_2.fastq
ERR011348_1.fastq
ERR011348_2.fastq
ERR011349_1.fastq
ERR011349_2.fastq
```


## Step 1: Run read_qc to trim the reads and remove human contamination
Note: you will need the bmtagger hg38 index to remove the human reads - see the MetaWrap2 database
installation instructions. You may also use another host genome to filter against with the `-x` option.
Alternatively, use the `--skip-bmtagger` flag of the `read_qc` module to only do the read trimming.

Individually process each sample
```
mkdir READ_QC
metawrap2 read_qc -1 RAW_READS/ERR011347_1.fastq -2 RAW_READS/ERR011347_2.fastq -t 24 -o READ_QC/ERR011347
metawrap2 read_qc -1 RAW_READS/ERR011348_1.fastq -2 RAW_READS/ERR011348_2.fastq -t 24 -o READ_QC/ERR011348
metawrap2 read_qc -1 RAW_READS/ERR011349_1.fastq -2 RAW_READS/ERR011349_2.fastq -t 24 -o READ_QC/ERR011349
```

Alternatively, process all samples at the same time with a parallel for loop (especially if you have many samples):
```
for F in RAW_READS/*_1.fastq; do
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap2 read_qc -1 $F -2 $R -t 1 -o READ_QC/$SAMPLE &
done
```

Or as a one-liner: `for F in RAW_READS/*_1.fastq; do R=${F%_*}_2.fastq; BASE=${F##*/}; SAMPLE=${BASE%_*}; metawrap2 read_qc -1 $F -2 $R -t 1 -o READ_QC/$SAMPLE & done`


Lets have a glance at one of the output folders: `ls READ_QC/ERR011347`

These are html reports of the read quality before and after QC:
```
post-QC_report
pre-QC_report
```

Original raw reads:
![Read quality before QC](https://i.imgur.com/x8aaFWs.png)
Final QC'ed reads:
![Read quality before QC](https://i.imgur.com/drJfxC9.png)


These are the final trimmed and de-contaminated reads:
```
final_pure_reads_1.fastq
final_pure_reads_2.fastq
```

Move over the final QC'ed reads into a new folder
```
mkdir CLEAN_READS
for i in READ_QC/*; do
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done
```


## Step 2: Assembling the metagenomes with the assembly module
Note: Depending on your goals you may want to assemble each sample separately, but for the purposes of
analyzing the whole community across samples, we will be co-assembling our samples.

Concatenate the reads from all the samples:
```
cat CLEAN_READS/ERR*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/ERR*_2.fastq > CLEAN_READS/ALL_READS_2.fastq
```

Assemble the reads with the `--metaspades` flag (metaSPAdes is usually preferred over MEGAHIT unless you
have a very large data set; running with no assembler flag uses MEGAHIT, and `--metaspades --megahit`
runs the hybrid pipeline):
```
metawrap2 assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 200 -t 24 --metaspades -o ASSEMBLY
```

You will find the assembly file in `ASSEMBLY/final_assembly.fasta`, and the QUAST assembly report html in `ASSEMBLY/assembly_report.html`!

Assembly statistics:
![Assembly stats](https://i.imgur.com/RbDldGU.png)


Looking at the top 10 contigs shows we got some longer contigs (considering that we are working with just 7Gbp of data)!
```
grep ">" ASSEMBLY/final_assembly.fasta | head

>NODE_1_length_196124_cov_2.427049
>NODE_2_length_176373_cov_3.889994
>NODE_3_length_163601_cov_3.070200
>NODE_4_length_142996_cov_2.771017
>NODE_5_length_109931_cov_3.516837
>NODE_6_length_106321_cov_2.842875
>NODE_7_length_99368_cov_2.860703
>NODE_8_length_95669_cov_2.506714
>NODE_9_length_91511_cov_12.466716
>NODE_10_length_88949_cov_2.730882
```

## Step 3: Run the kraken2 module on both the reads and the assembly
Running KRAKEN2 on the reads will give us an idea of the taxonomic composition of the communities in the
three samples, while running it on the assembly will give us an idea of what taxonomic groups were
assembled better than others (the assembly process is heavily biased and should not be used to infer
overall community composition).

Note: you will need a KRAKEN2 database for this module (see the database installation instructions).

Run the `kraken2` module on all files at once, subsetting the reads to 1M reads per sample to speed up the run
```
metawrap2 kraken2 -o KRAKEN -t 24 -s 1000000 CLEAN_READS/ERR*fastq ASSEMBLY/final_assembly.fasta
```

Lets have a look at the output folder:
```
ERR011347.krak2   ERR011348.krak2   ERR011349.krak2   final_assembly.krak2
ERR011347.kraken2 ERR011348.kraken2 ERR011349.kraken2 final_assembly.kraken2
ERR011347.krona   ERR011348.krona   ERR011349.krona   final_assembly.krona
kronagram.html
```

The `.krak2` files contain the raw KRAKEN2 output, the `.kraken2` files contain the translated taxonomy
lineage of each read or contig, and the `.krona` files summarize taxonomy statistics to be fed into
KronaTools, which makes the `kronagram.html` file. The `kronagram.html` file contains all the taxonomy
information from all the samples and the co-assembly. Inspecting the kronas in a web browser will show you
what the community composition is like.

For example, here is the taxonomic composition of our first sample:
![Krona](https://i.imgur.com/jZiFPUV.png)



## Step 4: Bin the co-assembly with three different algorithms with the binning module

The initial binning process with CONCOCT, MaxBin2, and metaBAT2 will be the more time-intensive step
(especially CONCOCT and MaxBin2), so you may want to run the `binning` module with each algorithm
separately. However, MetaWrap2 supports running all three together. Our dataset is reasonably small, so
we will run all three binning predictions at the same time.

If you are used to different binning software, feel free to run them instead. The downstream refinement
process (the `bin_refinement` module) takes in up to 3 different bin sets, although you can get around
this by splitting your bin sets into groups and then recursively consolidating them.

Run the binning module with all three binners - notice how the F and R read files go at the end of the command.
```
metawrap2 binning -o INITIAL_BINNING -t 24 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/ERR*fastq
```

In the output folder, we see folders with the 3 final bin sets and a `work_files` directory with the
intermediate alignments and depth files.
```
concoct_bins	maxbin2_bins  metabat2_bins  work_files
```
Looking inside these folders reveals how many bins each binner found. But we do not know how good these
bins are yet (unless you used the `--run-checkm` flag). We will find out the quality of the bins in the
next step!


## Step 5: Consolidate bin sets with the bin_refinement module
Note: make sure you downloaded the CheckM database (`metawrap2 install-db checkm` does this and points
CheckM at it for you).

**Optional: score bins with CheckM2 instead of CheckM1.** CheckM1 uses pplacer, which needs roughly 40 GB
of RAM per thread and is the usual reason this step fails or crawls on a smaller machine. CheckM2 is fast,
has no pplacer, and needs far less memory:
```
metawrap2 install-env bin_refinement-checkm2      # one-time
metawrap2 bin_refinement --checkm2 -o BIN_REFINEMENT -t 24 -c 50 -x 10 \
    -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/
```
The `.stats.tsv` files have the same format either way, so every later step works unchanged. The scores
themselves are not identical between CheckM1 and CheckM2, so pick one and stay with it within a study.
If you keep CheckM1 and hit its memory wall, `--quick` (CheckM's `--reduced_tree`) is the other lever.

Now that you have metaBAT2, MaxBin2, and CONCOCT bins, lets consolidate them into a single, stronger bin
set! If you used your own binning software, feel free to use any 3 bin sets. If you have more than 3, you
can run them in groups. For example, if you have 5 bin sets, try consolidating 1+2+3 and 4+5, and then
consolidate again between the outputs.

When you do your refinement, put some thought into the minimum completion (`-c`) and maximum
contamination (`-x`) parameters that you enter. During refinement, MetaWrap2 has to choose the best
version of each bin among up to 7 versions. It will dynamically adjust to prioritize the bin quality that
you desire. Consider this example: bin_123 comes in four versions in terms of completion/contamination:
95/15, 90/10, 80/5, 70/5. Which one is the best version? The high-completion but high-contamination one,
or the less complete but purer bin? This is subjective and depends on what you value in a bin and on your
purposes for bin extraction.

By default, the minimum completion is 70% and maximum contamination is 10%. However, because of the
relatively poor depth of these demonstration samples, we will set minimum completion to 50% and maximum
contamination to 10%, but feel free to be much more picky. Parameters like `-c 90 -x 5` are not
unreasonable on some data (but you will get fewer bins, of course).

Run the `bin_refinement` module:
```
metawrap2 bin_refinement -o BIN_REFINEMENT -t 24 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 50 -x 10
```

In the output directory, you will see the three original bin folders we fed in, as well as the
`metawrap_50_10_bins` directory, which contains the final, consolidated bins. You will also see a
`Binning_refiner` bin set - this is an internal benchmark produced by Binning_refiner. You can ignore this
set; it will likely have low contamination and completion. You will also see `.stats.tsv` files for each one
of the bin directories.
```
concoct_bins.stats.tsv	maxbin2_bins.stats.tsv	metabat2_bins.stats.tsv	metawrap_50_10.stats.tsv	Binning_refiner.stats.tsv
concoct_bins		maxbin2_bins		metabat2_bins		metawrap_50_10_bins	Binning_refiner
```

Note that the `_50_10_` part of the naming refers to the `-c` and `-x` options you chose. You can repeat
the run with the same output directory using different options to get different results and see what works
best on your sample (using `_90_5_` for example to get more near-complete genomes). The re-calculation
will reuse the existing binning outputs, greatly reducing run time.

The `.stats.tsv` files contain useful information about each bin, including its completeness and contamination.
For example, `cat BIN_REFINEMENT/metawrap_50_10.stats.tsv`:
```
bin	completeness	contamination	GC	lineage	N50	size	binner
bin_005	100.0	1.6	0.311	Euryarchaeota	12686	1705532	binsO.checkm
bin_004	99.32	1.342	0.408	Clostridiales	58825	2083650	binsO.checkm
bin_014	86.69	5.896	0.293	Bacteria	3754	2199676	binsO.checkm
bin_006	86.22	2.348	0.371	Clostridiales	4283	2055792	binsO.checkm
bin_008	83.16	2.516	0.446	Clostridiales	2723	1467846	binsO.checkm
bin_002	80.34	0.0	0.469	Bacteria	11936	3579466	binsO.checkm
bin_009	76.57	2.648	0.425	Selenomonadales	3155	1796524	binsO.checkm
bin_013	74.82	1.710	0.435	Bacteroidales	7456	3643185	binsO.checkm
bin_003	74.53	0.377	0.284	Clostridiales	10440	1241933	binsO.checkm
bin_010	65.78	0.0	0.263	Bacteria	3045	1159966	binsO.checkm
bin_011	64.85	3.776	0.417	Bacteroidales	2086	3103352	binsO.checkm
bin_001	57.36	0.0	0.430	Bacteria	4628	2673426	binsO.checkm
bin_007	52.94	1.724	0.501	Bacteria	3614	1465011	binsO.checkm
```

To evaluate how many "good bins" (based on our >50% comp., <10% cont. metric) MetaWrap2 produced, we can run
```
cat BIN_REFINEMENT/metawrap_50_10_bins.stats.tsv | awk '$2>50 && $3<10' | wc -l
13
```

By inspecting the other files, we find that metaBAT2, MaxBin2, CONCOCT, and MetaWrap2 produced 11, 7, 10,
and 13 bins, respectively. So MetaWrap2 produced 2 more bins than the best single binner. Not bad! But this
is just the number of bins.

To more closely compare the bin sets in terms of completion and contamination, we can look at the plots in `BIN_REFINEMENT/figures/`:
![Bin_refinement](https://i.imgur.com/m6RRJxi.jpg)

Note: This graph no longer has `Binning_refiner` in it, to reduce confusion. If you want to see
Binning_refiner's performance, look at binsABC (or binsAB if you have two bin sets) in the other figure.

As you can see, the refinement process produced the best bin set in terms of both completion and
contamination. Keep in mind that these improvements are even more dramatic in more complex samples.


## Step 6: Visualize the community and the extracted bins with the blobology module
Lets use the `blobology` module to project the entire assembly onto a GC vs. abundance plane and annotate
it with taxonomy and bin information. This will not only give us an idea of what these microbial
communities are structured like, but will also show us our binning success in a more visual way.

Note: you will need the NCBI_nt and NCBI_tax databases for this module (see the database installation instructions).

NOTE: In order to annotate the blobplot with bins with the `--bins` flag, you **MUST use the
non-reassembled bins**! In other words, use the bins produced by the `bin_refinement` module, not the
`reassemble_bins` module.

```
metawrap2 blobology -a ASSEMBLY/final_assembly.fasta -t 24 -o BLOBOLOGY --bins BIN_REFINEMENT/metawrap_50_10_bins CLEAN_READS/ERR*fastq
```

You will find that the output has a number of blobplots of our communities, annotated with different levels
of taxonomy or their bin membership. Note that to help with visualizing the bins, some of the plots only
contain the contigs that were successfully binned (these files have `.binned.blobplot.` in their names).
If your assembly is very large, you can randomly subsample the contigs to plot with `--subsample INT`.

Phyla taxonomy of the entire assembled community:
![Phyla](https://i.imgur.com/VihLGWb.jpg)

Bin membership of all the contigs:
![bins](https://i.imgur.com/GDmIYe5.jpg)


## Step 7: Find the abundances of the draft genomes (bins) across the samples
We would like to know how the extracted genomes are distributed across the samples and in what abundance
each bin is present in each sample. The `quant_bins` module can give us this information. It uses Salmon -
a tool conventionally used for transcript quantitation - to estimate the abundance of each scaffold in each
sample, and then computes the average bin abundances.

NOTE: In order to run this module, it is **highly** recommended to use the non-reassembled bins (the bins
produced by the `bin_refinement` module, not the `reassemble_bins` module) and to provide the entire
non-binned assembly with the `-a` option. This will give more accurate bin abundances that are in context
of the entire community.

Lets run the `quant_bins` module:
```
metawrap2 quant_bins -b BIN_REFINEMENT/metawrap_50_10_bins -o QUANT_BINS -a ASSEMBLY/final_assembly.fasta CLEAN_READS/ERR*fastq
```

The output contains several useful files. First, there is the `bin_abundance_heatmap.png` - a quick heatmap
made to visualize the bin abundances across the samples.
![heatmap](https://i.imgur.com/K1RaPUT.png)


The raw data for this plot (as you will most likely want to make your own heatmaps to analyze) is in
`bin_abundance_table.tab`. Note that the abundances are expressed as "genome copies per million reads" and
are calculated with Salmon in a way similar to how TPM (transcripts per million) is calculated in RNAseq
analysis. As such, they are already standardized to the individual sample size.

```
Genomic bins	ERR011349	ERR011348	ERR011347
bin_009	0.113912116828	35.851964987	39.8440491514
bin_010	0.273774684856	9.52869077293	39.988244574
bin_001	7.87827599808	31.3262582417	72.4475075589
bin_004	1.11852631889	100.052540293	111.213423224
bin_002	42.0242612674	69.0094806385	80.200001212
bin_005	2.16260151787	22.06396779	43.7720962538
bin_011	64.2884105466	25.3703846834	29.5444322752
bin_006	517.890689122	0.379711918465	0.834196723864
bin_007	0.499019812767	61.2121001057	82.5953338481
bin_014	5.49635966692	14.631433905	32.98399834
bin_013	0.230760165209	56.0018529273	91.6502833521
bin_008	98.4767064505	38.0691238971	22.8857472565
bin_003	349.730007621	0.0911113402849	0.196554603409
```

Finally, you can view the abundances of all the individual contigs in all the samples in the `quant_files` folder.

## Step 8: Re-assemble the consolidated bin set with the reassemble_bins module
Now that we have our final, consolidated bin set in `BIN_REFINEMENT/metawrap_50_10_bins`, we can try to
further improve it with reassembly. The `reassemble_bins` module will collect the reads belonging to each
bin and then reassemble them separately with a "permissive" and a "strict" algorithm. Only the bins that
improved through reassembly will be altered in the final set.

Note: make sure you downloaded the CheckM database (see the MetaWrap2 database instructions).

Let us run the `reassemble_bins` module with all the reads we have:
```
metawrap2 reassemble_bins -o BIN_REASSEMBLY -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -t 24 -m 800 -c 50 -x 10 -b BIN_REFINEMENT/metawrap_50_10_bins
```

The strict and permissive read-recruitment stringency can be tuned with `--strict-cut-off` (default 2
mismatches) and `--permissive-cut-off` (default 5 mismatches), and long reads can be added with
`--nanopore <reads>`.

Looking at the output in `BIN_REASSEMBLY/reassembled_bins.stats.tsv`, we can see which bins were improved
through strict reassembly, which through permissive reassembly, and which could not be improved
(`.strict`, `.permissive`, and `.orig` bin extensions, respectively):
```
bin	completeness	contamination	GC	lineage	N50	size	binner
bin_010.orig	65.78	0.0	0.263	Bacteria	3045	1159966	NA
bin_007.strict	54.94	0.671	0.501	Clostridiales	3947	1474089	NA
bin_004.permissive	99.32	1.342	0.408	Clostridiales	72135	2088821	NA
bin_002.permissive	82.06	0.0	0.469	Bacteria	18989	3604843	NA
bin_014.strict	85.84	3.066	0.293	Bacteria	4576	2201824	NA
bin_009.permissive	76.74	2.554	0.425	Selenomonadales	3601	1802438	NA
bin_013.permissive	78.37	1.357	0.435	Bacteroidales	9887	3675176	NA
bin_011.orig	64.85	3.776	0.417	Bacteroidales	2086	3103352	NA
bin_006.permissive	88.04	1.006	0.371	Clostridiales	6288	2070146	NA
bin_005.orig	100.0	1.6	0.311	Euryarchaeota	12686	1705532	NA
bin_003.permissive	74.91	0.396	0.284	Clostridiales	16578	1243641	NA
bin_001.orig	57.36	0.0	0.430	Bacteria	4628	2673426	NA
bin_008.strict	83.89	1.342	0.446	Clostridiales	3870	1474833	NA
```

But how much did our bin set really improve? We can look at the `BIN_REASSEMBLY/reassembly_results.png`
plot to compare the original and reassembled sets. We can see that the bin reassembly modestly improved the
bin N50 and completion metrics, and significantly reduced contamination. Fantastic!
![heatmap](https://i.imgur.com/V8IosYQ.jpg)


We can also view the CheckM plot of the final bins in `BIN_REASSEMBLY/reassembled_bins.png`:
![heatmap](https://i.imgur.com/Yx00fuQ.png)



## Step 9: Determine the taxonomy of each bin with the classify_bins module
Note: you will need the NCBI_nt and NCBI_tax databases for this module
(`metawrap2 install-db blast taxdump`). If your BLAST database is not named `nt`, set `BLASTDB_NAME`
in the config alongside `BLASTDB`.

**Optional: classify against GTDB with GTDB-Tk.** GTDB-Tk replaces the whole
MEGABLAST -> taxator-tk chain with one `classify_wf` call, against a curated, consistently-named
genome taxonomy. It is usually the better answer for MAGs, and it needs neither NCBI_nt nor NCBI_tax:
```
metawrap2 install-env classify_bins-gtdbtk        # one-time; then download the GTDB reference data
metawrap2 classify_bins --gtdbtk -b BIN_REASSEMBLY/reassembled_bins -o BIN_CLASSIFICATION -t 24
```
The output contract is the same `BIN_CLASSIFICATION/bin_taxonomy.tab` (bin name, tab, lineage), so
anything you built on top of this module keeps working.

We already got an idea of the approximate taxonomy of each bin from CheckM in the `.stats.tsv` files in the
`bin_refinement` and `reassemble_bins` modules. We can do better than that, however. The `classify_bins`
module uses MEGABLAST and taxator-tk to accurately assign taxonomy to each contig and then consolidates
the results to estimate the taxonomy of the whole bin. Of course the success and accuracy of our
predictions will rely heavily on the existing database.

Estimate the taxonomy of our final, reassembled bins with the `classify_bins` module:
```
metawrap2 classify_bins -b BIN_REASSEMBLY/reassembled_bins -o BIN_CLASSIFICATION -t 24
```

We can view the final estimated taxonomy in `BIN_CLASSIFICATION/bin_taxonomy.tab`:
```
bin_001.orig.fasta	Bacteria;Firmicutes;Clostridia;Clostridiales
bin_005.orig.fasta	Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter;Methanobrevibacter smithii
bin_011.orig.fasta	Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides
bin_002.permissive.fasta	uncultured organism
bin_010.orig.fasta	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae
bin_014.strict.fasta	Bacteria;Firmicutes
bin_008.strict.fasta	Bacteria
bin_009.permissive.fasta	Bacteria
bin_006.permissive.fasta	Bacteria;Firmicutes;Clostridia;Clostridiales
bin_003.permissive.fasta	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae
bin_004.permissive.fasta	Bacteria
bin_013.permissive.fasta	Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae
bin_007.strict.fasta	Bacteria
```
As you can see, some of the bins are annotated very deeply, while others can only be classified as
"Bacteria". This method is relatively trustworthy, but it often fails to annotate organisms that are very
distant from anything in the NCBI database. For these tricky bins, manually looking at marker genes (such
as ribosomal proteins) can result in much more sensitive taxonomy assignment.


## Step 10: Functionally annotate bins with the annotate_bins module
Now that we have our finalized reassembled bins, we are ready to functionally annotate them for downstream
functional analysis. This module simply annotates the genes with PROKKA - it cannot do the actual
functional analysis for you.

Run the functional annotation module on the final, reassembled bins:
```
metawrap2 annotate_bins -o FUNCT_ANNOT -t 24 -b BIN_REASSEMBLY/reassembled_bins/
```

**Optional: annotate with Bakta instead of PROKKA.** PROKKA is no longer maintained; Bakta is its
actively-developed successor and gives richer, better-curated annotations:
```
metawrap2 install-env annotate_bins-bakta         # one-time
metawrap2 annotate_bins --bakta --bakta-db /path/to/bakta_db -o FUNCT_ANNOT -t 24 \
    -b BIN_REASSEMBLY/reassembled_bins/
```
The output folders (`bin_funct_annotations`, `bin_translated_genes`, `bin_untranslated_genes`) are laid
out identically, so downstream analysis does not change. You can also set `BAKTA_DB` in the config
instead of passing `--bakta-db` every time.

You will find the functional annotations of each bin in GFF format in the `FUNCT_ANNOT/bin_funct_annotations` folder
```
head FUNCT_ANNOT/bin_funct_annotations/bin_001.orig.gff
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	2866	3645	.	-	0	ID=HMOHEJHL_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00001;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	3642	4478	.	-	0	ID=HMOHEJHL_00002;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00002;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	4606	5859	.	-	0	ID=HMOHEJHL_00003;inference=ab initio prediction:Prodigal:2.6;locus_tag=HMOHEJHL_00003;product=hypothetical protein
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	5856	6575	.	-	0	ID=HMOHEJHL_00004;Name=ypdB;gene=ypdB;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P0AE39;locus_tag=HMOHEJHL_00004;product=Transcriptional regulatory protein YpdB
NODE_75_length_31799_cov_0.983871	Prodigal:2.6	CDS	6603	7658	.	-	0	ID=HMOHEJHL_00005;eC_number=1.1.1.261;Name=egsA;gene=egsA;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P94527;locus_tag=HMOHEJHL_00005;product=Glycerol-1-phosphate dehydrogenase [NAD(P)+]
```

You will also find the translated and untranslated predicted genes in fasta format in
`FUNCT_ANNOT/bin_translated_genes` and `FUNCT_ANNOT/bin_untranslated_genes` folders. Finally, you can find
the raw PROKKA output files in `FUNCT_ANNOT/prokka_out`.
