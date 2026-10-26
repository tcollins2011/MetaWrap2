# Scoring genomes: `score`

`metawrap2 score` runs the scoring steps of `bin_refinement` (CheckM, GUNC, MIMAG tier) on any
directory of genome FASTAs and writes the same table `bin_refinement` writes. It does not refine,
dereplicate or delete anything.

```
metawrap2 score -t 32 -o SCORED refined/metawrap_70_10_bins
metawrap2 score -t 32 -o SCORED -x .fa metawrap1_run/metawrap_70_10_bins refined/metawrap_70_10_bins
metawrap2 score -t 32 -o SCORED --checkm1 published_mags/
```

## Purpose

Use `score` to compare bin sets from different pipelines or publications with one scorer.
metaWRAP 1 refines with CheckM1 and MetaWrap2 with CheckM2, and the two disagree most in the
high-completeness tail, where the 90/5 near-complete line sits. Bin sets that carry their own
pipelines' scores cannot be compared on near-complete counts. `score` rescores every bin set with
one tool, conda env and database.

## Output

`-o` holds one `<bin set>.stats.tsv` per input directory, plus the scorers' own output
(`<bin set>.checkm/`, `<bin set>.gunc/`, a log per tool). The genomes are symlinked into
`-o/<bin set>/`, and the input directories are not written to.

The table has the columns `bin_refinement` writes, in the same order: the CheckM header, then
`gunc_css`/`gunc_pass` from GUNC, then `mimag_tag`/`mimag_rrna`/`mimag_trnas`. Anything that reads
a refined bin set's `.stats.tsv` (the benchmark scorer, the MAG table, scripts keyed on column
positions) reads a rescored one unchanged.

An arbitrary directory of genomes has no annotation output, so a genome passing 90/5 is
tagged `high_quality_unverified_rrna`: it has not been searched for a 16S. To resolve these, run
`metawrap2 annotate --bakta` over the genomes and score them again.

With several input directories, a comparison table is printed at the end. `clean near-complete` is
the number of genomes that are near-complete and not called chimeric by GUNC, counted per genome
from its own row.

## Flags

| Flag | Meaning |
|---|---|
| `-x`, `--extension` | extension of the genome FASTAs, with or without the dot (default `.fasta`) |
| `-c`, `--completeness` | completeness floor for the passing count (default 70) |
| `--contamination` | contamination ceiling for the passing count (default 10) |
| `--checkm1` | score with CheckM1 instead of CheckM2 |
| `--no-gunc` | skip the chimerism check and its reference database |
| `--score-everything` | score genomes whose size rules them out, instead of recording 0/100 |

`-x` is the FASTA extension here, as in CheckM, GUNC and GTDB-Tk, and not the contamination ceiling
as in `bin_refinement`. `-c` means completeness in both commands. The contamination ceiling has only
the long form, `--contamination`. metaWRAP 1 writes `.fa` and MetaWrap2 writes `.fasta`, so pass
`-x .fa` for metaWRAP 1 bin sets.

## Requirements

CheckM2 (the default) needs the `metawrap2-checkm2` env; CheckM1 and GUNC are in `bin_refinement`'s
`metawrap2-binqc` env. GUNC runs by default and needs `GUNC_DB`, which is checked before CheckM
starts.

```
metawrap2 install-env checkm2
metawrap2 install-env bin_refinement
metawrap2 install-db gunc
```

## Notes

- Inputs with the same basename (for example two `metawrap_70_10_bins` from different pipelines)
  are named after their parent directory as well, so their tables and rows stay separate.
- One CheckM score cache covers the whole run, keyed on each genome's contigs, so genomes shared
  between bin sets are scored once. The cache file is named for the estimator, so a `--resume` that
  changes `--checkm1` does not reuse the other estimator's scores.
- `metawrap2 report` reads at most one `.stats.tsv` per directory, so on a `score` output directory
  it shows only one of the bin sets. Use the comparison table that `score` prints to compare them.
- The unbinned catch-all (`unbinned.fasta`) is skipped, since it is not a genome.
