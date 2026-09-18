# MetaWrap2

MetaWrap2 is an easy-to-use, **Python-native** pipeline for genome-resolved metagenomic data
analysis: read QC, assembly, taxonomic profiling, binning, hybrid **bin refinement**, bin
**reassembly**, quantification, classification, and annotation. Each module is a standalone
program, so you use only the parts you need.

It is the successor to metaWRAP: same algorithms and external tools (same results), but the
pipeline is entirely Python (no Bash, Perl, or R), with modern packaging and a separate conda
environment per module.

## Highlights
- **One command per module** (`metawrap2 <module>`), each also runnable standalone.
- **Per-module conda environments** (`metawrap2-<module>`), driven by mamba, so dependencies
  resolve quickly and independently.
- **Visible, editable commands** - the exact command each module runs lives in a `COMMANDS`
  block at the top of the module file; tune it there or via `metawrap2.toml`.
- **Provenance built in** - every run records its full config, exact commands, tool versions,
  and a lineage that propagates across chained runs.
- **`metawrap2 test`** - probe each module's env and software and get an installed / missing /
  broken status map.

See [Installation](installation.md) to get started.
