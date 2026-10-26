# Trimmer reports for the read-retention line

All files are hand-written, following each tool's schema. A real `fastp.json` is ~180 kB of
per-cycle quality curves, and a real trimming report has a full adapter table; the retention line
reads neither.

## `fastp_low_retention.json`

Follows fastp 0.24's JSON report: `summary.before_filtering`, `summary.after_filtering` and
`filtering_result`, with key names and units taken from fastp's writer (fastp counts reads, not
pairs). Only the keys `_fastp_retention` reads are present, plus `command`, so a test fails if a
required key stops being parsed instead of a default being used.

The numbers are fastp's report for the MetaHIT tutorial cohort (ERR011347, 75 bp reads from about
2010):

    before:            25,242,208 reads
    after:             10,486,282 reads
    low_quality_reads: 14,755,926
    too_short_reads:   0

58.5% of the library is discarded as whole reads on quality and no reads are shortened, because
fastp's defaults filter reads and do not trim them.

## `fastp_high_retention.json`

The same schema with numbers typical of a modern library (99.1% of pairs kept). It tests the case
where no warning is printed.

## `trim_galore_paired_report.txt`

Follows the end of a Trim Galore 0.6.10 `*_trimming_report.txt` for a paired run: the cutadapt
`=== Summary ===` block, the per-file `RUN STATISTICS` block, and the two pair-level lines the
parser uses (`Total number of sequences analysed...` and `Number of sequence pairs removed...`).
Other sections (adapter tables, the length-distribution histogram) are left out; the parser
searches for the lines and does not depend on their order.
