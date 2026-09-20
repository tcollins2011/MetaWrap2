# Provenance & reproducibility

Two records are kept: one **per run**, beside the data it produced, and one **per machine**, which
outlives any output directory.

## Per run, in the output directory

| File | Contents |
|---|---|
| `run_config.json` | full run record: module, version, all resolved parameters (incl. defaults), config file, conda env, databases, platform, timing, inputs, exact commands, tool versions used, and the complete lineage |
| `run_environment.json` | what the run used, written *before* any tool runs: the environment's full conda package list, the MetaWrap2 git commit, and each database's content fingerprint |
| `run_commands.txt` | the exact commands run, in order; the first line is the top-level `metawrap2 <module>` reproducer |
| `provenance.txt` | readable history of every MetaWrap2 step that led to this output |
| `metawrap2.log` | one timestamped log of MetaWrap2's own messages and every line of every tool (rotates at 256 MB) |
| `run.stdout` / `run.stderr` | raw tool output (streamed live to the screen *and* captured here) |

`run_environment.json` is written up front on purpose: a run that dies in its first tool is
exactly the run whose environment you want to inspect, and it used to leave no record of what it
was running.

### Identifying a database, not just pointing at it

A path is not an identity. Replace the contents of a database directory with a newer build — which
is what happens when a database is updated in place — and every earlier run's record names
something that no longer exists while claiming to describe what produced the data.

So `run_environment.json` records, per database, which of its defining files are present, their
total size, and a fingerprint built from their names, sizes and modification times. That changes
when the database is rebuilt and not when it is merely read, which makes "which of my analyses
used the `nt` I just replaced?" answerable. `metawrap2 doctor` reports the same fingerprints.

### Which tree, not just which version

Between releases, `2.1.0` names many different trees — including one with local edits, which is
the normal state of a package whose design invites editing the module files. When MetaWrap2 runs
from a git checkout, `run_environment.json` also records the commit, branch, `git describe`, and
whether the working tree was dirty. From an installed wheel there is no repository, and it says so
(`{"kind": "installed"}`), which is itself informative: the version string is then the whole truth.

## Lineage propagation

When an input is itself a MetaWrap2 output directory, that input's history is merged in, so
provenance propagates across chained runs (QC → assembly → binning → bin_refinement →
reassembly). If a step consumed several MetaWrap2 inputs, all of their histories are merged and
de-duplicated. A final output's `provenance.txt` therefore traces the whole path back to the
raw reads.

## Per machine, for as long as MetaWrap2 is installed

Every module invocation also appends one line to `~/.metawrap2/history.jsonl` (or
`$METAWRAP2_HISTORY`), outside any output directory. It survives deleting the output, is never
pruned by MetaWrap2, and stays greppable. A failed history write is ignored — losing a history
line must never be the reason a run fails.

```bash
metawrap2 history                     # the 20 most recent runs
metawrap2 history --failed            # what didn't finish
metawrap2 history --module binning --commands   # exact command lines, to re-run or paste
metawrap2 history --stats             # runs and time spent, per module
metawrap2 history --diff RUN RUN      # what differed between two runs
```

A run records a `started` entry and then its outcome, so a run that was killed still leaves a
trace rather than vanishing.

## Resuming: verified, not assumed

`--resume` consults the study's **run manifest** (`manifest.json`), which records each step's
command, duration, outcome, and every declared input and output with a content fingerprint. A step
is skipped only if it completed, its outputs are all still present and non-empty, and none of its
inputs has changed.

This matters because a marker file proves none of those things. The old behaviour was wrong in two
directions: a step whose outputs had been deleted was skipped anyway, and a step whose *inputs* had
changed was also skipped — silently mixing results from two different inputs.

Fingerprints hash the size plus the first and last megabyte, so the check costs the same on a
100 GB library as on a small one. `--strict-fingerprints` hashes whole files, for the case where an
edit between the sampled ends matters.

```bash
metawrap2 status study/                  # per step: status, duration, and what --resume would do
metawrap2 status study/ --step binning   # its inputs, outputs, command, and measured cost
metawrap2 status study/ --audit          # is every recorded output still present and non-empty?
```

The `--resume` column comes from the same check `--resume` itself makes, not a reimplementation of
it. The same audit runs automatically at the end of every run, so a step that reported success
while producing nothing is reported then rather than discovered on the next resume.

## Reproducing a run

- Re-run the top line of `run_commands.txt`.
- Rebuild the environment from its committed lockfile:
  `conda create -n metawrap2-binning --file src/metawrap2/envs/locks/binning.linux-64.lock`
  (or `metawrap2 install-env binning --from-lock`).
- Check the database fingerprints in `run_environment.json` against what you have now.
- `metawrap2 --dry-run <module> ...` prints and records the commands without executing them.
- `metawrap2 --resume <module> -o <dir> ...` continues a run, re-running anything it cannot verify.

## What each step cost

Steps are reaped with `os.wait4`, so the manifest records that step's exact peak memory, CPU time
and wall time. `metawrap2 calibrate <study>...` compares those measurements against MetaWrap2's
built-in disk and memory estimates.

One caveat it states plainly rather than hiding: the output-to-input ratio is a *lower bound* on
the disk multiplier, because scratch is deleted when a step finishes and no manifest can record how
large it got at its peak. Peak memory is measured exactly.
