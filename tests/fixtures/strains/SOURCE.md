# inStrain output, as it really looks

Real inStrain 1.10.0 output, unedited, from mapping three human gut metagenomes (ERR011347-9)
against the 26 refined bins of this repo's test study.

`genome_info.tsv` is ERR011349's profile - the worst-covered of the three samples, which is the
useful one: most of its genomes sit below the coverage floor, so it is the fixture for "present but
not deeply enough to measure". Note `breadth` and `breadth_minCov` are different columns and the
second is the one a diversity figure is measured over.

`genomeWide_compare.tsv` is all 42 pairwise comparisons. It is the fixture for the two traps:
the column is `percent_compared`, not `percent_genome_compared`, and popANI, conANI and
`percent_compared` are all stored as **fractions** while inStrain's documentation and paper quote
popANI as a percentage ("99.999%"). Parsed literally against the published threshold, every
comparison here is a different strain over too little of the genome - which is what the first real
run reported.

`tracs_transmission_distances.csv` is real TRACS 1.0.1 output over the same two samples and three of
the same bins, unedited. It is the fixture for two things its format gets wrong in a way that is
easy to miss: the `MSA file` column truncates the reference name at the first dot, so a run over
`bin.1` and `bin.4` writes `combinedbin` in *both* rows, and the sample labels are
`<sample>_<genome>` concatenated with no escaping. The genome can only be recovered from the second.

Note also that `tracs pipe` skips the first line of its input file with a bare `next(infile)` and
documents no header - so a file without one silently loses its first sample, which is how the first
real run produced an empty distance table from two good samples.
