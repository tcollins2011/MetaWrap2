# MetaWrap2

MetaWrap2 is a Python-native pipeline for genome-resolved metagenomic data analysis: read QC,
assembly, taxonomic profiling, binning, hybrid bin refinement, bin reassembly, cross-sample
coverage, dereplication, quantification, classification, and annotation. Each module is a
standalone program, so you use only the parts you need.

It began as a reimplementation of metaWRAP and is no longer only that. The first release reproduced
metaWRAP's pipeline exactly — same algorithms, tools and thresholds — because an equivalent
baseline is what lets you migrate and reproduce existing analyses before anything changes
underneath you. That baseline is still pinned by golden tests. Since then four scientific defaults
have moved, each on a measurement, and three modules have been added that metaWRAP has no
equivalent for. The shipped configuration now recovers more clean high-quality MAGs than metaWRAP 1
does, for less compute.

What that means in numbers, how it was measured, and what is still open is in
[Benchmarks](benchmarks.md). Start with [Installation](installation.md).

## Highlights

- One command per module (`metawrap2 <module>`), each also runnable standalone.
- Five workflow conda environments driven by mamba — grouped by what runs together rather than one
  per module, so dependencies resolve quickly and an install is a few commands rather than a dozen.
- Visible, editable commands: the exact command each module runs lives in a `COMMANDS` block at the
  top of the module file. Tune it there or in `metawrap2.toml`.
- Defaults chosen by measurement. `metawrap2 benchmark` runs a declared comparison against a
  known-truth dataset and commits the result, so a default that moves arrives in review next to the
  numbers that justify it.
- Cohorts as a first-class input: `run` and `cohort` execute a study from a sample sheet, `coverage`
  builds one cross-sample depth table from N alignments rather than N², and `dereplicate` collapses
  per-sample bins into a catalogue of organisms.
- Chimerism detection on by default. GUNC annotates bins whose markers look clean but whose
  sequence comes from two genomes — which CheckM cannot see.
- Provenance built in: every run records its full config, exact commands, tool versions, database
  fingerprints, the conda package list, the MetaWrap2 commit, and a lineage that propagates across
  chained runs.
- A `--resume` that verifies: a step is skipped only if it completed, its outputs are still present
  and non-empty, and none of its inputs has changed. `metawrap2 status` shows what a resume would do
  before you run it.
- A durable run history: `metawrap2 history` records every run ever started on this machine, outside
  any output directory, and `--diff` compares two of them.
- One resource budget per run — threads, memory and GPUs — defaulting to a share of the machine
  rather than a number you have to guess, with a disk-space preflight that stops a run before it
  fills a shared filesystem.
- `metawrap2 doctor` probes each module's env, software and databases and prints an installed /
  missing / broken map; `--fix` recreates exactly what is broken.

Some output files are named so their extensions say what they are: bins are `bin_001.fasta` and
statistics tables `*.stats.tsv`, with contents unchanged from metaWRAP's.
