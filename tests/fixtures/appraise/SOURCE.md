# singlem appraise output, as it really looks

Both files are real SingleM 0.21.4 output, from the 26 refined bins of ERR011347 appraised against
that sample's raw reads with the GTDB r232 metapackage (S6.5.0). `unaccounted_for.tsv` is the first
15 rows of 1,183; `appraise_summary.tsv` is complete.

They are fixtures because the shapes being parsed are not documented anywhere:

- the summary carries `total` and `average` rows beside the per-sample ones, which are SingleM's
  own summaries and must not be averaged back in;
- `coverage_assembled` is absent unless an assembly was appraised too, so a parser that reads
  columns positionally or assumes the column reports "nothing assembled" where the truth is
  "no assembly was given";
- the unaccounted table has one row per marker gene per sequence, so one unrepresented organism
  appears up to 59 times and counting rows counts markers.

The appraisal itself is the result worth remembering: those 26 bins represent 11.6% of the
bacterial community by coverage, against a far higher share of the assembled base pairs.
