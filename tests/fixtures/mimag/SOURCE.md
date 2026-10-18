# MIMAG rRNA/tRNA evidence fixtures

**These five files are hand-written, not vendored.** None of them is real tool output. They are
written to the shape of the formats `metawrap2.mimag` parses, from the documented layout of each
format, and they are here because the alternative needs a Bakta run: `annotate --bakta` wants a
~70 GB Bakta database and an Infernal/Rfam pass over a real genome, which no test in this suite can
do. What each file is modelled on, and what in it is deliberate:

## `bakta/bin_001/bin_001.tsv`, `bakta/bin_002/bin_002.tsv`

Modelled on Bakta 1.9.2's `<prefix>.tsv`: five `#`-comment lines carrying the software and database
versions, then a `#Sequence Id Type Start Stop Strand Locus Tag Gene Product DbXrefs` header, then
one tab-separated row per feature. The `Type` values used are `cds`, `rRNA` and `tRNA`; products are
`16S ribosomal RNA` and `tRNA-Ala(nnn)`, which is how Bakta names them.

The preamble is five lines rather than one because the parser locates the header by finding the
first `#`-line that has both a `Type` and a `Product` column, not by skipping a fixed count - Bakta
has changed how many preamble lines it writes between versions.

* `bin_001` carries all three rRNAs and tRNAs for 18 distinct amino acids: MIMAG high quality.
* `bin_002` is the same annotation run with the 16S absent and only 17 amino acids covered. Both
  failures are *measured*, which is what makes it the counterpart to an unannotated genome.

Both files end with a `tRNA-Xxx(???)` row. `Xxx` is not an amino acid, and it is there so the tRNA
count is a count of amino acids rather than of tRNA features - with features it would read 19/18.

## `bakta/bin_001/bin_001.gff3`

The same annotation as `bin_001.tsv` in Bakta's GFF3 form: `##gff-version 3`, feature type in
column 3 (`rRNA`, `tRNA`), and `Name=`/`product=` attributes in column 9. It exists so the two
parsers can be shown to agree on one genome; the GFF is the fallback, because its attribute keys
vary with Bakta's version while the TSV's columns do not.

## `rna/bin_003.infernal.tblout`

Modelled on `cmscan --tblout --fmt 2`: two `#`-comment header lines, then whitespace-padded columns
with an index first and the Rfam family name second. Rows are `5S_rRNA` (RF00001),
`SSU_rRNA_bacteria` (RF00177) and `LSU_rRNA_bacteria` (RF02541) - the real accessions for those
families.

The fourth row is `5_8S_rRNA` (RF00002), a eukaryotic family. It is in the file on purpose: a
parser that matched "5S" as a substring would read it as a 5S hit, so this row is what keeps the
family match exact. Coordinates, scores and E-values are plausible but invented; nothing reads them.

## `rna/bin_003.trnascan.tsv`

Modelled on tRNAscan-SE 2.0's default tabular output: a two-line wrapped header, a `---` rule, then
columns of which the fifth is the tRNA type. Eighteen amino-acid rows, plus one `Undet` and one
`Pseudo` row, both of which appear in real output and neither of which names an amino acid - they
are there so the count stays at 18 rather than 20.

Together `bin_003`'s two files are the split case: Infernal measured the rRNAs and tRNAscan-SE the
tRNAs, which is the fallback a run that never called `annotate --bakta` would have.
