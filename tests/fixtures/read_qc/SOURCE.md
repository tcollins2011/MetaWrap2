# Trimmer reports, for the read-retention line

**Hand-written, modelled on each tool's schema.** Nothing here is captured output - a real
`fastp.json` is ~180 kB of per-cycle quality curves and a real trimming report carries a full
adapter table, none of which the retention line reads.

## `fastp_low_retention.json`

Modelled on fastp 0.24's JSON report: `summary.before_filtering`, `summary.after_filtering` and
`filtering_result`, with the key names and the units (fastp counts **reads**, not pairs) taken from
fastp's own writer. Only the keys `_fastp_retention` reads are present, plus `command`, so that a
test that stops parsing a key it needs fails rather than silently reading a default.

The numbers are the MetaHIT tutorial cohort's (ERR011347, 2010-era 75 bp reads) as reported by
fastp in the user report this fixture exists for:

    before:            25,242,208 reads
    after:             10,486,282 reads
    low_quality_reads: 14,755,926
    too_short_reads:   0

That is the shape of the defect: 58.5% of the library discarded **whole** on quality, and exactly
zero reads shortened, because fastp's defaults filter and do not trim.

## `fastp_high_retention.json`

The same schema with the numbers a modern library gives - 99.1% of pairs kept - so the "does not
warn" half of the threshold is tested against a real-shaped report rather than an edited one.

## `trim_galore_paired_report.txt`

Modelled on the tail of a Trim Galore 0.6.10 `*_trimming_report.txt` for a paired run: the
cutadapt `=== Summary ===` block, the per-file `RUN STATISTICS` block, and the two pair-level lines
the parser actually uses (`Total number of sequences analysed...` and `Number of sequence pairs
removed...`). Intervening sections (adapter tables, the length-distribution histogram) are cut -
the parser searches, it does not walk the file in order.
