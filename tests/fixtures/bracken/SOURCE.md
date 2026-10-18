# Bracken output, as it really looks

## Provenance

Every `.bracken` file here is **real, unedited output from Bracken 3.0.1** (bioconda `bracken=3.1`,
whose wrapper reports `Bracken v3.0.1`), run locally on 2026-10-14. The *inputs* to those runs are
hand-written, because the alternative is an 8 GB Kraken2 standard database in the test suite:

- `sample1.kreport`, `sample2.kreport` - **hand-written** Kraken2 reports in the six-column default
  format (`--report`, not `--report-minimizer-data`). 1,000 and 800 reads over three species in two
  phyla, with reads deliberately left stranded at the genus and phylum nodes so there is something
  for Bracken to redistribute. Modelled on the format `kraken2-inspect` emits and on the local
  standard database's `inspect.txt`: percent, clade reads, taxon reads, rank code, taxid, and a name
  indented two spaces per level.
- A **hand-written** minimal Bracken database, one `database100mers.kmer_distrib` in the format
  `generate_kmer_distribution.py` writes (`mapped_taxid<TAB>genome_taxid:kmers_mapped:total_kmers`,
  space-separated per genome), with entries for the two genera, the two phyla and each species'
  self-mapping. Not vendored here - it is an input to the fixture, not a fixture.
- `sample1_species.bracken` - Bracken `-r 100 -l S -t 10` on `sample1.kreport`.
- `sample1_genus.bracken` - the same report at `-l G`.
- `sample2_species.bracken` - Bracken `-r 100 -l S -t 10` on `sample2.kreport`.
- `sample1_mismatched_db.bracken` - the same run as `sample1_species.bracken` against a database
  stripped of its parent-node entries, which is what a Bracken database built from a *different*
  Kraken2 database looks like from the inside: every lookup misses.

Two files are **hand-written and are not Bracken output**, and exist only to be rejected:

- `percent_not_fraction.bracken` - `sample1_species.bracken` with `fraction_total_reads` multiplied
  by 100. Bracken never writes this; it is the shape the column would have if it were a percentage.
- `sample1_bracken_outreport.kreport` - the first lines of the six-column, header-less file Bracken
  writes to `-w`, so the reader can be checked for refusing it instead of parsing it positionally.

## Why these are fixtures

The column layout is nowhere in Bracken's documentation, which describes the columns in prose
without naming them. The header comes from a single `o_file.write` in `est_abundance.py`:

    name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads  new_est_reads  fraction_total_reads

and the rest of what these files pin down is not documented either:

- **`fraction_total_reads` is a fraction in 0-1, not a percentage.** `est_abundance.py` writes
  `"%0.5f" % (new_est_reads / sum_all_reads)`, so the column sums to 1 across the table, not to 100.
  On `sample1_species.bracken` it sums to 0.99894 - the shortfall is the five-decimal rounding, which
  is why the parser's tolerance is loose rather than exact.
- **`taxonomy_lvl` is the rank Bracken was asked for, repeated on every row**, not the taxon's own
  rank. `sample1_genus.bracken` carries `G` on both of its rows. That makes a genus table and a
  species table distinguishable, which matters because `--resume` will otherwise pick up a table
  written at a different `--bracken-rank`.
- **A taxon below `-t` is absent, not zero.** *S. mutans* has 0 reads in `sample2.kreport` and no
  row at all in `sample2_species.bracken`, while it has a row in sample 1. A parser that assumes
  every sample's table lists the same taxa will mis-align the samples.
- **`added_reads` is the only column that says whether Bracken did anything.**
  `sample1_mismatched_db.bracken` has `added_reads` of 0 on every row and fractions that still sum
  to exactly 1.00000. Nothing else in that file reveals that the estimate is Kraken2's assigned
  counts renormalised rather than an abundance estimate.

## The numbers worth remembering

The reason Bracken exists, visible in these two files. In `sample1.kreport`, *S. pneumoniae* has 150
of 1,000 reads at its own species node - 15% - while 450 reads sit above it at the *Streptococcus*
and Firmicutes nodes. `sample1_species.bracken` estimates 301 reads for it, 31.7%. Reading the
Kraken2 report's species rows as composition would have understated it by a factor of two.

Against the mismatched database the same sample yields 150 reads and 30.0% for the same species,
with 450 reads redistributed to nothing - and Bracken exits 0 in that case, as it does for a missing
database and for a missing `kmer_distrib`, because its wrapper ends each error path with `echo`
followed by a bare `exit`.
