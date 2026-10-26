# Bracken output

## Provenance

Every `.bracken` file here is unedited output from Bracken 3.0.1 (bioconda `bracken=3.1`, whose
wrapper reports `Bracken v3.0.1`), run locally on 2026-10-14. The inputs to those runs are
hand-written, to avoid an 8 GB Kraken2 standard database in the test suite:

- `sample1.kreport`, `sample2.kreport`: hand-written Kraken2 reports in the six-column default
  format (`--report`, not `--report-minimizer-data`). 1,000 and 800 reads over three species in two
  phyla, with some reads assigned at the genus and phylum nodes for Bracken to redistribute.
  Modelled on the format `kraken2-inspect` emits and on the local standard database's
  `inspect.txt`: percent, clade reads, taxon reads, rank code, taxid, and a name indented two
  spaces per level.
- A hand-written minimal Bracken database: one `database100mers.kmer_distrib` in the format
  `generate_kmer_distribution.py` writes (`mapped_taxid<TAB>genome_taxid:kmers_mapped:total_kmers`,
  space-separated per genome), with entries for the two genera, the two phyla and each species'
  self-mapping. It is an input to the fixtures and is not included here.
- `sample1_species.bracken`: Bracken `-r 100 -l S -t 10` on `sample1.kreport`.
- `sample1_genus.bracken`: the same report at `-l G`.
- `sample2_species.bracken`: Bracken `-r 100 -l S -t 10` on `sample2.kreport`.
- `sample1_mismatched_db.bracken`: the same run as `sample1_species.bracken` against a database
  with its parent-node entries removed. A Bracken database built from a different Kraken2 database
  behaves the same way: every lookup misses.

Two files are hand-written, are not Bracken output, and are used to test rejection:

- `percent_not_fraction.bracken`: `sample1_species.bracken` with `fraction_total_reads` multiplied
  by 100. Bracken does not write this; it is what the column would look like as a percentage.
- `sample1_bracken_outreport.kreport`: the first lines of the six-column, header-less file Bracken
  writes to `-w`. The reader must refuse it instead of parsing it by position.

## Format details

Bracken's documentation describes the columns in prose without naming them. The header comes from
a single `o_file.write` in `est_abundance.py`:

    name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads  new_est_reads  fraction_total_reads

The following are also undocumented:

- `fraction_total_reads` is a fraction in 0-1, not a percentage. `est_abundance.py` writes
  `"%0.5f" % (new_est_reads / sum_all_reads)`, so the column sums to 1. On
  `sample1_species.bracken` it sums to 0.99894 because of the five-decimal rounding, so the
  parser's tolerance is loose.
- `taxonomy_lvl` is the rank Bracken was run at, repeated on every row, not the taxon's own rank.
  `sample1_genus.bracken` has `G` on both rows. This distinguishes a genus table from a species
  table, which `--resume` needs so it does not pick up a table written at a different
  `--bracken-rank`.
- A taxon below `-t` has no row; it is not written as zero. *S. mutans* has 0 reads in
  `sample2.kreport` and no row in `sample2_species.bracken`, but has a row in sample 1. Samples
  must be aligned by taxon, not by row.
- `added_reads` is the only column that shows whether Bracken redistributed any reads.
  `sample1_mismatched_db.bracken` has `added_reads` of 0 on every row and fractions that sum to
  1.00000. Nothing else in the file shows that the values are Kraken2's assigned counts
  renormalised, not an abundance estimate.

## Example numbers

In `sample1.kreport`, *S. pneumoniae* has 150 of 1,000 reads at its own species node (15%), and 450
reads are assigned above it at the *Streptococcus* and Firmicutes nodes. `sample1_species.bracken`
estimates 301 reads for it, 31.7%. Reading the Kraken2 species rows as composition gives about half
that.

Against the mismatched database the same sample gives 150 reads and 30.0% for this species, with
the 450 reads not redistributed. Bracken exits 0 in this case, and also for a missing database or a
missing `kmer_distrib`, because its wrapper ends each error path with `echo` followed by a bare
`exit`.
