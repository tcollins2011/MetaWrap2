# reassemble_bins candidate stats fixtures

**All three files are hand-written, not real CheckM output.** They are written to the shape of
`metawrap2.checkm.write_stats` - the frozen `bin completeness contamination GC lineage N50 size`
header, then one tab-separated row per bin, rows sorted by completeness descending - and they are
here because the real thing needs CheckM: a lineage_wf pass over 18 candidate genomes, ~40 GB of
pplacer RAM and a reference tree, which no test in this suite can do.

The numbers are modelled on the 26-bin run that reported this module's defect, where the original
bin was excluded from scoring and 4 of 26 bins (15%) were replaced with a reassembly that scored
worse than the bin they started from. The proportion here is deliberately higher - 3 of 5 chosen
bins keep their original - because every "original wins" shape worth distinguishing is given its
own bin rather than being left to chance.

## `candidates.stats.tsv`

What one CheckM pass over `reassembled_bins/` produces: all three candidates for each of six bins,
18 rows. This is `choose_best_bin`'s input. At the module defaults, `-c 70 -x 10`:

| bin | winner | why |
|---|---|---|
| `bin_001` | `strict` | the ordinary case: both reassemblies beat the original, strict by score |
| `bin_002` | **`orig`** | `permissive` is 1.5 points *more complete* than the original but 3.4 points more contaminated, and contamination is penalised 5x, so 585.5 beats 570.0 |
| `bin_003` | `permissive` | the other ordinary case |
| `bin_004` | **`orig`** | both reassemblies lose completeness outright; `strict` and `permissive` tie at 569.0 and are separated by N50, and the original wins anyway |
| `bin_005` | **`orig`** | `strict` falls to 68.0% completeness and `permissive` rises to 14.0% contamination, so neither is a candidate at all and the original is the only one left |
| `bin_006` | none | every version is under 70% complete, so the bin is not selected and appears in none of the three counts |

`bin_002` is the fixture the reporter asked for. Its reassemblies have the *higher* N50s (41 kb and
52 kb against the original's 26 kb), so a selection that preferred the longest contigs, or that
could not see the original at all, would pick `permissive` and the test would fail. The scores are
`completeness + 5 * (100 - contamination)`, which is `metawrap2.scripts.choose_best_bin`'s ranking.

## `reassembled_bins.stats.tsv`

The same five winners re-scored after promotion, which is what the second CheckM run writes and what
the summary plots are drawn from. Values are copied from the winning rows above, because CheckM is
deterministic on the same contigs - re-scoring a promoted bin cannot change its numbers. It is named
`reassembled_bins.stats.tsv` rather than anything shorter because `plot_reassembly` labels a series
by its filename and takes a different label x-position for a series whose name contains `reasse`.

## `original_bins.stats.tsv`

Header and nothing else. This is exactly the file the defect produced: the originals' table is built
by filtering the work stats for rows whose bin name ends in `.orig`, and while the `.orig` copies
were written as `.fa` in a directory CheckM was scanning with `-x fasta`, there were none - the run
had `26 permissive / 26 strict / 0 orig`. A header-only table becomes an empty series, which is what
all three panels of `plot_reassembly` indexed into and died on, 62m53s into a finished run. It is
kept as a fixture rather than fixed because the same file is still legitimately produced whenever
every original fails `-c`/`-x`.
