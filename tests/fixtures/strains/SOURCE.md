# inStrain and TRACS output

inStrain 1.10.0 output, unedited, from mapping three human gut metagenomes (ERR011347-9) against
the 26 refined bins of this repo's test study.

`genome_info.tsv` is the profile of ERR011349, the lowest-coverage sample. Most of its genomes are
below the coverage floor, so it tests genomes that are present but not covered enough to measure.
`breadth` and `breadth_minCov` are separate columns; diversity is measured over `breadth_minCov`.

`genomeWide_compare.tsv` holds all 42 pairwise comparisons. It tests two format details:

- the column is `percent_compared`, not `percent_genome_compared`;
- popANI, conANI and `percent_compared` are stored as fractions, while inStrain's documentation
  and paper quote popANI as a percentage ("99.999%"). Compared against the published threshold
  without conversion, every comparison here reads as a different strain over too little of the
  genome.

`tracs_transmission_distances.csv` is TRACS 1.0.1 output, unedited, over the same two samples and
three of the same bins. It tests two format details:

- the `MSA file` column truncates the reference name at the first dot, so a run over `bin.1` and
  `bin.4` writes `combinedbin` in both rows;
- sample labels are `<sample>_<genome>` concatenated without escaping. The genome can only be
  recovered from these labels.

`tracs pipe` skips the first line of its input file with a bare `next(infile)` and documents no
header, so an input file without a header loses its first sample without a warning. With two
samples, the result is an empty distance table.
