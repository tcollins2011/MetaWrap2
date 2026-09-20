# MetaWrap2

MetaWrap2 is an ambitious, **Python-native** pipeline for genome-resolved metagenomic data
analysis: read QC, assembly, taxonomic profiling, binning, hybrid **bin refinement**, bin
**reassembly**, quantification, classification, and annotation. Each module is a standalone
program, so you use only the parts you need.

**This is phase 1.** MetaWrap2 is not a port - it is aimed at what metagenomics looks like now, and
the plan is considerably more ambitious than what is here today. Because the work after this one
changes results, the foundation came first: the established genome-recovery pipeline reproduced
*exactly* - same algorithms, same tools, same thresholds, pinned by golden tests - on a harness
with verified resume, real provenance, reproducible environments and a test suite. That equivalence
is a checkpoint, not a destination: it is what lets you migrate, reproduce existing analyses, and
trust the harness before any science changes underneath you.

Some output files were renamed so their extensions say what they are - bins are now
`bin_001.fasta` and statistics tables `*.stats.tsv` - with unchanged contents.

## Highlights
- **One command per module** (`metawrap2 <module>`), each also runnable standalone.
- **Per-module conda environments** (`metawrap2-<module>`), driven by mamba, so dependencies
  resolve quickly and independently.
- **Visible, editable commands** - the exact command each module runs lives in a `COMMANDS`
  block at the top of the module file; tune it there or via `metawrap2.toml`.
- **Provenance built in** - every run records its full config, exact commands, tool versions,
  database fingerprints, the conda package list, the MetaWrap2 commit, and a lineage that
  propagates across chained runs.
- **`--resume` that verifies** - a step is skipped only if it completed, its outputs are still
  present and non-empty, and none of its inputs has changed. `metawrap2 status` shows what a
  resume would do before you run it.
- **A durable run history** - `metawrap2 history` records every run ever started on this machine,
  outside any output directory, and `--diff` compares two of them.
- **`metawrap2 doctor`** - probe each module's env, software and databases and get an installed /
  missing / broken status map; `--fix` recreates exactly what is broken.

See [Installation](installation.md) to get started.
