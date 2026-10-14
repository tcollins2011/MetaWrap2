# dbCAN output, as it really looks

Real dbCAN 5.2.9 output over proteins `metawrap2 annotate` called from two CAMI2 bins, unedited.
`overview.tsv` is bin_006 complete (65 CAZyme genes from 1,772 proteins); `overview_multidomain.tsv`
is the first rows of bin_003, which is a polysaccharide-degrading Bacteroides and holds the
multi-domain cases.

The schema is the reason these are fixtures. Nothing documents that:

- a method's call carries its coordinates (`GT119(45-389)`) and a multi-domain protein's calls are
  joined with `+` or `|`, so families have to be split out before methods can be compared;
- `dbCAN_sub` reports subfamilies (`GH13_e82`) where the other two report families (`GH13`), so
  agreement has to be computed on the family or every subfamily call looks like a disagreement;
- `Substrate` can be a list with unmapped members, `-;peptidoglycan`, which left alone becomes a
  substrate named `-`;
- a lone `-` is dbCAN's empty value in every column.

For reference, the counts these produce: bin_003 has 194 glycoside hydrolases in 73 families and
bin_006 has 24 in 18 - which is the Bacteroides-versus-not difference showing up in the table.
