# Release notes

One page per version: new features, changed defaults and the measurements behind them. The
[CHANGELOG](https://github.com/ursky/MetaWrap2/blob/main/CHANGELOG) lists every change.

| Version | Headline | Defaults that moved |
|---|---|---|
| [2.1.3](2.1.3.md) | Installation fixes, PDF report, MIMAG tags, Bracken, eggNOG-mapper, functional hierarchy images, tutorial, CAMI2 oral and skin benchmarks | none |
| [2.1.2](2.1.2.md) | Functional annotation, viruses, AMR, diversity coverage, strains, long reads: 13 new modules and flags | GUNC, CheckM2 and GTDB-Tk pins raised; CoverM identity floors set; `classify_bins` taxator-tk deprecated |
| [2.1.1](2.1.1.md) | Three-binner `--preset slow`, measured end to end | `--preset slow` = metaBAT2 + SemiBin2 + COMEBin |
| [2.1.0](2.1.0.md) | Benchmark harness, first benchmarked default, chimerism detection | default tier = metaBAT2 + SemiBin2; GUNC on by default; depth identity floor 90% |
| [2.0.1](2.0.1.md) | Pipeline driver, verified resume, run history, one-command install | none |
| [2.0.0](2.0.0.md) | Python port of metaWRAP with the same algorithms, tools and results | none |

**Reading the numbers.** A clean high-quality MAG is a bin that is ≥70% complete and ≤10%
contaminated by CheckM2 and non-chimeric by GUNC, counted by joining the two reports per bin. Cost
is CPU-hours measured by reaping the child process, and is reported next to quality rather than as
a ratio. There are two experimental designs and they are not mixed: end-to-end arms assemble their
own contigs, and binner arms bin the CAMI2 gold assembly. [Benchmarks](../benchmarks.md) describes
both.
