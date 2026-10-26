# Miniature finished run for the per-MAG table and the PDF report

All files are hand-written. Columns, spellings and units are copied from real output: a reference
study (READ_QC, ASSEMBLY, BIN_REFINEMENT, BIN_CLASSIFICATION, QUANT_BINS) and single-module runs
(`paths_gene`, `amr_out`, `vir3`). Real output is too large to commit (one FastQC HTML is 250 kB,
one Krona input 20 kB), and the tests are about which cells are empty.

In this run three genomes are scored, four are quantified, and bin names are spelled differently
across files. The reference study has the same pattern: 13 refined bins have CheckM scores and 26
pre-refinement bins have taxonomy.

## Files

`REFINEMENT/metawrap_70_10_bins.stats.tsv` uses `checkm.STATS_HEADER`, including CheckM's GC as a
fraction (0.412), not a percentage. Three genomes:

* `bin_001` is ordinary.
* `bin.2` is the metaWRAP-era spelling that `os.path.splitext` turns into `bin`. Its contamination
  is a measured `0.0`, which must appear in the table as `0.0`.
* `bin_003` has empty completeness and contamination: CheckM wrote the row but did not score it.
  The file tests a blank and a `0.0` in the same column.

`REFINEMENT/concoct_bins.stats.tsv` is a decoy. A `bin_refinement` output directory holds four
stats tables, one per input binner plus the consolidated set, and CONCOCT's sorts first. Its
genomes (`concoct_42`, `unbinned`) appear in no other file, so a reader that picks this table
produces rows with the wrong names.

`REFINEMENT/metawrap_70_10_bins.contigs.tsv` is the `contig -> bin` membership table. It lists
contigs for `bin_001` and `bin.2` and none for `bin_003`, so the contig count for `bin_003` must be
blank, not 0.

`REFINEMENT/metawrap_70_10_bins.gunc/GUNC.progenomes_2.1.maxCSS_level.tsv` is GUNC's maxCSS report.
The columns are GUNC 1.0.6's, and the file name includes the reference database version, which is
why `gunc.find_report` matches on prefix and suffix. Genomes are spelled `bin_001.fa` / `bin.2.fa`,
a third convention. `bin_003` is absent: a bin GUNC did not score is not a failed bin.

`CLASSIFY/bin_taxonomy.tsv` is `classify_bins`' two-column output, without a header, with genomes
spelled `bin_001.fasta`. `bin_003` has an empty lineage, meaning it was classified but not placed,
which differs from not having been classified.

`QUANT/bin_abundance_table.tsv` is `quant_bins`' matrix, with `Genomic bins` as the first header
cell. `bin_004` appears only here, and must still get a row.

`READ_QC/` has three samples:

* `sampleA` has FastQC HTML before and after QC (the `Total Sequences</td><td>N` markup of real
  reports, cut down to the Basic Statistics table);
* `sampleB` has `fastp.json` with fastp's `summary.before_filtering.total_reads` structure;
* `sampleC` has only a post-QC report, so its reads-in value is blank.

`ASSEMBLY/QUAST_out/report.tsv` uses QUAST 5.3.0's metric names, including the
`# contigs (>= 5000 bp)` and `auN` rows, which must not appear in the report.

`PATHWAYS/hierarchy/krona/` has two Krona inputs in `hierarchy.write_krona_input` format. The
second has a single zero-weight row, which `sunburst.draw` does not draw. A genome with no KO is a
valid result, and the report must name it instead of omitting its page.

`AMR/amr.tsv`, `AMR/virulence.tsv` and `VIRUSES/viral_contigs.tsv` use the columns `amr` and
`viruses` write. `virulence.tsv` has a header and no rows: the step ran and found nothing, which
the report must distinguish from a missing table.
