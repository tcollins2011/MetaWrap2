# Release notes

One page per version: new features, changed defaults and the measurements behind them. The
[CHANGELOG](https://github.com/ursky/MetaWrap2/blob/main/CHANGELOG) lists every change.

| Version | Headline | Defaults that moved |
|---|---|---|
| [2.1.4](2.1.4.md) | CheckM2 as the default scorer, `metawrap2 score`, bin sets scored in parallel in `bin_refinement`, progress heartbeat, low-utilisation warning, verified downloads, GPU check, external benchmark arms, `--preset gpu` (was `slow`), aligner index reuse, real-sample benchmarks on four datasets, CAMI2 bin counts against the known answer | CheckM2 replaces CheckM1 (`--checkm1` for the old scorer, now on the reduced tree); `--preset gpu` trains on cuda; heartbeat on |
| [2.1.3](2.1.3.md) | Installation fixes, PDF report, MIMAG tags, Bracken, eggNOG-mapper, functional hierarchy images, tutorial, CAMI2 oral and skin benchmarks | none |
| [2.1.2](2.1.2.md) | Functional annotation, viruses, AMR, diversity coverage, strains, long reads: 13 new modules and flags | GUNC, CheckM2 and GTDB-Tk pins raised; CoverM identity floors set; `classify_bins` taxator-tk deprecated |
| [2.1.1](2.1.1.md) | Three-binner `--preset slow`, measured end to end | `--preset slow` = metaBAT2 + SemiBin2 + COMEBin |
| [2.1.0](2.1.0.md) | Benchmark harness, first benchmarked default, chimerism detection | default tier = metaBAT2 + SemiBin2; GUNC on by default; depth identity floor 90% |
| [2.0.1](2.0.1.md) | Pipeline driver, verified resume, run history, one-command install | none |
| [2.0.0](2.0.0.md) | Python port of metaWRAP with the same algorithms, tools and results | none |

## Definitions

A clean high-quality MAG is a bin that is ≥70% complete and ≤10% contaminated by CheckM2 and
non-chimeric by GUNC, counted by joining the two reports per bin. Cost is CPU-hours measured by
reaping the child process, reported in its own column next to quality (no quality-per-cost ratio).
End-to-end arms assemble their own contigs; binner arms bin the CAMI2 gold assembly. Results from
the two designs are not compared with each other. [Benchmarks](../docs/benchmarks.md) describes
both.
