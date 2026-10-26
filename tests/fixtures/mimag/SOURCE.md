# MIMAG rRNA/tRNA evidence fixtures

All five files are hand-written, not tool output. They follow the documented layout of each format
`metawrap2.mimag` parses. Producing them would need `annotate --bakta`, which needs a ~70 GB Bakta
database and an Infernal/Rfam pass over a real genome, which the test suite cannot run.

## `bakta/bin_001/bin_001.tsv`, `bakta/bin_002/bin_002.tsv`

Follow Bakta 1.9.2's `<prefix>.tsv`: five `#`-comment lines with the software and database
versions, a `#Sequence Id Type Start Stop Strand Locus Tag Gene Product DbXrefs` header, then one
tab-separated row per feature. `Type` values are `cds`, `rRNA` and `tRNA`; products are
`16S ribosomal RNA` and `tRNA-Ala(nnn)`, as Bakta names them.

The preamble has five lines because the parser finds the header as the first `#`-line with both a
`Type` and a `Product` column, and does not skip a fixed number of lines. The number of preamble
lines has changed between Bakta versions.

* `bin_001` has all three rRNAs and tRNAs for 18 distinct amino acids: MIMAG high quality.
* `bin_002` is the same annotation with the 16S removed and 17 amino acids covered. Both shortfalls
  are measured, unlike a genome that was not annotated.

Both files end with a `tRNA-Xxx(???)` row. `Xxx` is not an amino acid; the row tests that the tRNA
count is of amino acids, not tRNA features (which would give 19/18).

## `bakta/bin_001/bin_001.gff3`

The same annotation as `bin_001.tsv` in Bakta's GFF3 format: `##gff-version 3`, feature type in
column 3 (`rRNA`, `tRNA`), and `Name=`/`product=` attributes in column 9. It tests that the two
parsers agree on one genome. The GFF is the fallback, because its attribute keys vary between Bakta
versions and the TSV columns do not.

## `rna/bin_003.infernal.tblout`

Follows `cmscan --tblout --fmt 2`: two `#`-comment header lines, then whitespace-padded columns
with an index first and the Rfam family name second. Rows are `5S_rRNA` (RF00001),
`SSU_rRNA_bacteria` (RF00177) and `LSU_rRNA_bacteria` (RF02541), the real accessions for those
families.

The fourth row is `5_8S_rRNA` (RF00002), a eukaryotic family. A parser that matched "5S" as a
substring would count it as a 5S hit, so the row tests exact family matching. Coordinates, scores
and E-values are invented and not read.

## `rna/bin_003.trnascan.tsv`

Follows tRNAscan-SE 2.0's default tabular output: a two-line wrapped header, a `---` rule, then
columns with the tRNA type fifth. There are eighteen amino-acid rows, plus one `Undet` and one
`Pseudo` row. Both occur in real output and neither names an amino acid, so the count must be 18,
not 20.

`bin_003`'s two files cover the split case: Infernal for the rRNAs and tRNAscan-SE for the tRNAs,
the fallback used when `annotate --bakta` was not run.
