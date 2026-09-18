# Modules

MetaWrap2 keeps metaWRAP's module set. See the repository's
[Module_descriptions.md](https://github.com/ursky/MetaWrap2/blob/main/Module_descriptions.md)
for full details; each module also prints `metawrap2 <module> -h`.

## Pre-processing
- **read_qc** - read trimming and host (e.g. human) read removal.
- **assembly** - metagenomic assembly with metaSPAdes or MEGAHIT, plus QC.
- **kraken2** - taxonomy profiling of reads or contigs, with a kronagram.

## Bins
- **binning** - initial binning with metaBAT2, MaxBin2, and/or CONCOCT.
- **bin_refinement** - consolidate multiple bin sets into a superior refined set.
- **reassemble_bins** - reassemble bins to improve N50/completion and reduce contamination.
- **quant_bins** - estimate bin abundance across samples.
- **blobology** - GC-vs-abundance blobplots of contigs and bins.
- **classify_bins** - conservative taxonomy assignment for bins.
- **annotate_bins** - functional annotation of draft genomes.

Each module runs in its own conda environment `metawrap2-<module>` (see
[Installation](installation.md)).
