# reassemble_bins candidate stats fixtures

All three files are hand-written, not CheckM output. They follow `metawrap2.checkm.write_stats`:
the fixed `bin completeness contamination GC lineage N50 size` header, then one tab-separated row
per bin, sorted by completeness descending. Producing them with CheckM would need a lineage_wf pass
over 18 candidate genomes, ~40 GB of pplacer RAM and a reference tree, which the test suite cannot
run.

The numbers are modelled on a 26-bin run in which the original bin was excluded from scoring, and 4
of 26 bins (15%) were replaced with a reassembly that scored lower than the original. Here 3 of 5
chosen bins keep their original, so that each case where the original wins has its own bin.

## `candidates.stats.tsv`

The output of one CheckM pass over `reassembled_bins/`: all three candidates for each of six bins,
18 rows. This is the input to `choose_best_bin`. At the module defaults, `-c 70 -x 10`:

| bin | winner | why |
|---|---|---|
| `bin_001` | `strict` | both reassemblies score above the original; strict scores highest |
| `bin_002` | `orig` | `permissive` is 1.5 points more complete than the original but 3.4 points more contaminated, and contamination is weighted 5x, so 585.5 beats 570.0 |
| `bin_003` | `permissive` | permissive scores highest |
| `bin_004` | `orig` | both reassemblies are less complete; `strict` and `permissive` tie at 569.0 and are separated by N50, and the original scores higher than both |
| `bin_005` | `orig` | `strict` falls to 68.0% completeness and `permissive` rises to 14.0% contamination, so neither passes the thresholds and only the original remains |
| `bin_006` | none | every version is under 70% complete, so the bin is not selected and appears in none of the three counts |

In `bin_002` the reassemblies have higher N50s (41 kb and 52 kb, against 26 kb for the original).
A selection that preferred the longest contigs, or that left out the original, would pick
`permissive` and fail the test. The score is `completeness + 5 * (100 - contamination)`, the
ranking used by `metawrap2.scripts.choose_best_bin`.

## `reassembled_bins.stats.tsv`

The five winners re-scored after promotion: the output of the second CheckM run, which the summary
plots are drawn from. Values are copied from the winning rows above, since CheckM gives the same
numbers for the same contigs. The file name matters: `plot_reassembly` labels a series by its file
name and uses a different label x-position for a series whose name contains `reasse`.

## `original_bins.stats.tsv`

A header with no rows. The originals' table is built by selecting work-stats rows whose bin name
ends in `.orig`. When the `.orig` copies were written as `.fa` into a directory CheckM scanned with
`-x fasta`, there were no such rows, and the run reported `26 permissive / 26 strict / 0 orig`. A
header-only table gives an empty series, and all three panels of `plot_reassembly` failed indexing
into it at the end of the run (62m53s). The same file is still produced when every original fails
`-c`/`-x`, so the plots must handle it.
