# A miniature finished run, for the per-MAG table and the PDF report

**Hand-written, not captured output.** Every file here was typed, with the columns, the spellings
and the units copied from real output: a reference study (READ_QC, ASSEMBLY, BIN_REFINEMENT,
BIN_CLASSIFICATION, QUANT_BINS) and single-module runs (`paths_gene`, `amr_out`, `vir3`). Real output was too
large to commit - one FastQC HTML is 250 kB and one Krona input is 20 kB - and the properties under
test are all about which cells are empty, which no amount of real data makes clearer.

The fixture is a run where **three genomes were scored, four were quantified, and nothing agrees
about how to spell a bin's name**. That is not contrived: it is what the reference study actually looks
like, where 13 refined bins have CheckM scores and 26 pre-refinement bins have taxonomy.

## What each file is modelled on, and what it is the fixture for

`REFINEMENT/metawrap_70_10_bins.stats.tsv` - `checkm.STATS_HEADER`, exactly, including CheckM's
habit of writing GC as a **fraction** (0.412) rather than a percentage. Three genomes:

* `bin_001` is ordinary.
* `bin.2` is the metaWRAP-era spelling that `os.path.splitext` turns into `bin`, and its
  contamination is **`0.0` - a real, measured zero**, which must survive into the table as `0.0`.
* `bin_003` has **empty** completeness and contamination: CheckM wrote the row but could not score
  it. Blank and `0.0` sitting in the same column of the same file is the whole point of this file.

`REFINEMENT/concoct_bins.stats.tsv` is a decoy. A `bin_refinement` output directory really does
hold four stats tables - one per input binner plus the consolidated set - and alphabetical order
puts CONCOCT's first. Its genomes (`concoct_42`, `unbinned`) appear in no other file here, so a
reader that picks the wrong table produces a table of names that are obviously not the run's.

`REFINEMENT/metawrap_70_10_bins.contigs.tsv` - the `contig -> bin` membership table. It lists
contigs for `bin_001` and `bin.2` and **none for `bin_003`**, so the contig count is the second
place a blank has to stay a blank rather than becoming a 0.

`REFINEMENT/metawrap_70_10_bins.gunc/GUNC.progenomes_2.1.maxCSS_level.tsv` - GUNC's maxCSS report.
The column list is GUNC 1.0.6's and the filename carries the reference database version, because
`gunc.find_report` matches on prefix and suffix for exactly that reason. Genomes are spelled
`bin_001.fa` / `bin.2.fa`, a **third** convention, and `bin_003` is absent - a bin GUNC did not
score has not failed.

`CLASSIFY/bin_taxonomy.tsv` - `classify_bins`' two-column output, with no header, genomes spelled
`bin_001.fasta`. `bin_003` has an **empty lineage**: classified, not placed, which is a different
statement from not having been classified.

`QUANT/bin_abundance_table.tsv` - `quant_bins`' matrix, with its literal `Genomic bins` first
header cell. `bin_004` appears **only here**, so it is the fixture for "a genome no other input
mentions still gets a row".

`READ_QC/` - three samples, three states. `sampleA` has FastQC HTML before and after (the
`Total Sequences</td><td>N` markup the real reports carry, cut down to the Basic Statistics table);
`sampleB` has `fastp.json` with fastp's real `summary.before_filtering.total_reads` shape; `sampleC`
has only a post-QC report, so its reads-in is a blank.

`ASSEMBLY/QUAST_out/report.tsv` - QUAST 5.3.0's metric names, including the `# contigs (>= 5000 bp)`
and `auN` rows that must *not* reach the page.

`PATHWAYS/hierarchy/krona/` - two Krona inputs in `hierarchy.write_krona_input` format. The second
has a single zero-weight row, which `sunburst.draw` correctly declines to draw - a genome with no KO
is a real result, and the report has to name it rather than silently produce one fewer page.

`AMR/amr.tsv`, `AMR/virulence.tsv`, `VIRUSES/viral_contigs.tsv` - the column lists `amr` and
`viruses` write. `virulence.tsv` is **header-only**: a table that ran and found nothing, which must
not be reported the same way as a table that does not exist.
