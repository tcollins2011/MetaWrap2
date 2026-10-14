# Provenance & reproducibility

MetaWrap2 keeps one record per run, in the output directory, and one per machine.

## Per run, in the output directory

| File | Contents |
|---|---|
| `run_config.json` | full run record: module, version, all resolved parameters (incl. defaults), config file, conda env, databases, platform, timing, inputs, exact commands, tool versions used, and the complete lineage |
| `run_environment.json` | what the run used, written *before* any tool runs: the environment's full conda package list, the MetaWrap2 git commit, and each database's content fingerprint |
| `run_commands.txt` | the exact commands run, in order; the first line is the top-level `metawrap2 <module>` reproducer |
| `provenance.txt` | readable history of every MetaWrap2 step that led to this output |
| `metawrap2.log` | one timestamped log of MetaWrap2's own messages and every line of every tool (rotates at 256 MB) |
| `run.stdout` / `run.stderr` | raw tool output (streamed live to the screen *and* captured here) |

`run_environment.json` is written first, so a run that fails early still records its environment.

### Database fingerprints

A database updated in place keeps its path but not its contents. For each database,
`run_environment.json` records:

- which of its defining files are present,
- their total size,
- a fingerprint built from their names, sizes and modification times.

The fingerprint changes on rebuild, not on read, so you can find which analyses used a database
you have since replaced. `metawrap2 doctor` reports the same fingerprints.

### The exact source tree

A version such as `2.1.0` can cover many source trees, including locally edited ones. From a git
checkout, `run_environment.json` also records the commit, branch, `git describe`, and whether the
tree was dirty. From an installed wheel it records `{"kind": "installed"}`, and the version string
is the full description.

## Lineage propagation

When an input is a MetaWrap2 output directory, its history is merged in (de-duplicated across
several inputs). Provenance carries across chained runs (QC → assembly → binning → bin_refinement →
reassembly), so a final `provenance.txt` traces back to the raw reads.

## Per machine, for as long as MetaWrap2 is installed

Every module invocation appends a line to `~/.metawrap2/history.jsonl` (or `$METAWRAP2_HISTORY`).
The file survives deleting outputs, is never pruned, and is greppable. A failed history write does
not stop the run.

```bash
metawrap2 history                     # the 20 most recent runs
metawrap2 history --failed            # what didn't finish
metawrap2 history --module binning --commands   # exact command lines, to re-run or paste
metawrap2 history --stats             # runs and time spent, per module
metawrap2 history --diff RUN RUN      # what differed between two runs
```

Each run writes a `started` entry and then its outcome, so a killed run still leaves a trace.

## Verified resume

`--resume` reads the run manifest (`manifest.json`): each step's command, duration, outcome, and
fingerprinted inputs and outputs. A step is skipped only if:

- it completed,
- its outputs are all still present and non-empty, and
- none of its inputs has changed.

A step whose outputs were deleted, or whose inputs changed, is re-run.

Fingerprints hash the size plus the first and last megabyte, so a 100 GB library checks as fast
as a small one. `--strict-fingerprints` hashes whole files.

```bash
metawrap2 status study/                  # per step: status, duration, and what --resume would do
metawrap2 status study/ --step binning   # its inputs, outputs, command, and measured cost
metawrap2 status study/ --audit          # is every recorded output still present and non-empty?
```

The `--resume` column uses the same check as `--resume`. The audit also runs at the end of every
run, flagging steps that reported success but produced nothing.

## Reproducing a run

- Re-run the top line of `run_commands.txt`.
- Rebuild the environment from its committed lockfile:
  `conda create -n metawrap2-bins --file src/metawrap2/envs/locks/bins.linux-64.lock`
  (or `metawrap2 install-env bins --from-lock`).
- Compare the database fingerprints in `run_environment.json` with what you have now.
- `metawrap2 --dry-run <module> ...` prints and records the commands without executing them.
- `metawrap2 --resume <module> -o <dir> ...` continues a run, re-running anything it cannot verify.

## What each step cost

Steps are reaped with `os.wait4`, so the manifest records each step's peak memory, CPU time and
wall time. `metawrap2 calibrate <study>...` compares them with the built-in disk and memory
estimates.

Caveat: the output-to-input ratio is a lower bound on the disk multiplier, because scratch is
deleted when a step finishes and its peak is not recorded. Peak memory is exact.
