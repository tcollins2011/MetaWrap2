# Provenance & reproducibility

Every module writes provenance into its `-o` output directory, next to the data:

| File | Contents |
|---|---|
| `run_config.json` | full run record: module, version, all resolved parameters (incl. defaults), config file, conda env, databases, platform, timing, inputs, exact commands, tool versions used, and the complete lineage |
| `run_commands.txt` | the exact commands run, in order; the first line is the top-level `metawrap2 <module>` reproducer |
| `provenance.txt` | readable history of every MetaWrap2 step that led to this output |
| `run.stdout` / `run.stderr` | full tool output (streamed live to the screen *and* captured here) |

## Lineage propagation

When an input is itself a MetaWrap2 output directory, that input's history is merged in, so
provenance propagates across chained runs (QC → assembly → binning → bin_refinement →
reassembly). If a step consumed several MetaWrap2 inputs, all of their histories are merged and
de-duplicated. A final output's `provenance.txt` therefore traces the whole path back to the
raw reads.

## Reproducing a run

- Re-run the top line of `run_commands.txt` to reproduce the output.
- Use `metawrap2 --dry-run <module> ...` to see/record the commands without executing them.
- Use `metawrap2 --resume <module> -o <dir> ...` to continue an interrupted run, skipping
  stages that already finished (tracked under `<dir>/.metawrap2/steps/`).
