# Release notes

One page per version: what it added, which defaults moved, and the measurements behind them. The
[CHANGELOG](https://github.com/ursky/MetaWrap2/blob/main/CHANGELOG) lists every change; these pages summarise
what to expect from upgrading.

| Version | Headline | Defaults that moved |
|---|---|---|
| [2.1.2](2.1.2.md) | Function, viruses, resistance, diversity coverage, strains - 13 new modules and flags | GUNC, CheckM2 and GTDB-Tk pins raised; CoverM identity floors set; `classify_bins` taxator-tk deprecated |
| [2.1.1](2.1.1.md) | The top binning tier is an ensemble again, measured end to end | `--preset slow` = metaBAT2 + SemiBin2 + COMEBin |
| [2.1.0](2.1.0.md) | The measurement harness, the first default it moved, and chimerism detection | shipped tier = metaBAT2 + SemiBin2; GUNC on by default; depth identity floor 90% |
| [2.0.1](2.0.1.md) | A pipeline driver, verified resume, run history, one-command install | none |
| [2.0.0](2.0.0.md) | Python-native metaWRAP: same algorithms, same tools, same results | none |

**Reading the numbers.** Clean high-quality MAGs are bins that are ≥70% complete / ≤10%
contaminated by CheckM2 and non-chimeric by GUNC, counted by joining the two reports per bin. Cost
is CPU-hours measured by reaping the child process, reported beside quality and never divided into
it. Two experimental designs appear and are never mixed: end-to-end arms assemble their own contigs
(what a user gets), and binner arms bin CAMI2's gold assembly (what a binner does).
[Benchmarks](../benchmarks.md) explains both.
