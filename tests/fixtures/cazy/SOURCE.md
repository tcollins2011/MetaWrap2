# dbCAN output

dbCAN 5.2.9 output, unedited, over proteins `metawrap2 annotate` called from two CAMI2 bins.
`overview.tsv` is the complete table for bin_006 (65 CAZyme genes from 1,772 proteins).
`overview_multidomain.tsv` is the first rows of bin_003, a polysaccharide-degrading Bacteroides,
and holds the multi-domain cases.

The schema properties tested here are not documented by dbCAN:

- a method's call carries its coordinates (`GT119(45-389)`), and a multi-domain protein's calls are
  joined with `+` or `|`, so families are split out before methods are compared;
- `dbCAN_sub` reports subfamilies (`GH13_e82`) where the other two methods report families
  (`GH13`), so agreement is computed at family level;
- `Substrate` can be a list with unmapped members (`-;peptidoglycan`); unhandled, this yields a
  substrate named `-`;
- a lone `-` is the empty value in every column.

Counts from these files: bin_003 has 194 glycoside hydrolases in 73 families, bin_006 has 24 in 18.
