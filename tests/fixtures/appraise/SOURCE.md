# SingleM appraise output

SingleM 0.21.4 output for the 26 refined bins of ERR011347, appraised against that sample's raw
reads with the GTDB r232 metapackage (S6.5.0). `unaccounted_for.tsv` is the first 15 of 1,183 rows;
`appraise_summary.tsv` is complete.

Format properties tested here, none documented by SingleM:

- the summary has `total` and `average` rows beside the per-sample rows; these are SingleM's own
  summaries and are excluded when averaging;
- `coverage_assembled` is present only when an assembly was also appraised, so the parser reads
  columns by name, and a missing column means "no assembly given", not "nothing assembled";
- the unaccounted table has one row per marker gene per sequence, so one unrepresented organism can
  appear up to 59 times, and a row count is a marker count.

In this appraisal the 26 bins represent 11.6% of the bacterial community by coverage, a much lower
share than they hold of the assembled base pairs.
