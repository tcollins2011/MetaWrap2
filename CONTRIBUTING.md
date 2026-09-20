# Contributing to MetaWrap2

Thanks for helping improve MetaWrap2! It is a pure-Python package - no Bash, Perl, or R.

## Setup
```
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e ".[dev]"
pre-commit install                      # format/lint what you touch, per commit
pre-commit install --hook-type pre-push  # run the full gate before anything leaves your machine
```

## Running the tests
```
python -m pytest -q                 # the env-independent unit + integration tests
metawrap2 doctor                    # probe module conda envs, software and databases
```
No bioinformatics tool is needed for the test suite. Two kinds of test are worth knowing about
before you change behaviour:

- **Golden tests** (`tests/test_bin_refinement_golden.py`) run the real refinement algorithm -
  real Binning_refiner, real consolidation and dereplication, with only CheckM replaced by a
  deterministic stand-in - and compare every output file against recorded values. If your change
  is not meant to alter results and one of these fails, the failure is the notification.
- **The golden end-to-end test** (`tests/test_golden_end_to_end.py`) puts stub tools on PATH and
  runs binning into bin_refinement, quant_bins and annotate_bins, pinning the output tree and the
  `.stats.tsv` contract byte-for-byte. It catches what unit tests structurally cannot: a module
  writing where the next one does not look.

Changing a golden value is fine when the change is intended - do it in the same commit, so the
diff shows what the output used to be.

## Project layout
- `src/metawrap2/modules/<module>.py` - one flat "recipe" per pipeline module. The commands a
  module runs live in the `COMMANDS` block and the tunable numbers in the `CONSTANTS` block at
  the top of the file. Keep them there and visible.
- `src/metawrap2/command.py` - the single place that runs an external command (mamba env
  wrapping, logging, provenance, dry-run). Don't call `subprocess` directly in a module.
- `src/metawrap2/{config,provenance,logging,checkpoint,checkm,constants,utils}.py` and
  `modules/_common.py` - shared helpers.
- `src/metawrap2/{manifest,validate,history,scratch,runlock,interrupt,usage,dbcheck}.py` - the
  infrastructure behind verified resume, fail-fast input checks, the durable run history, the
  scratch policy, the output-directory lock, signal handling, per-step measurement, and database
  identity. Each module's docstring explains what problem it exists for.
- `src/metawrap2/envs/<module>.yaml` - the per-module conda env (`metawrap2-<module>`), with
  committed exact-build lockfiles in `src/metawrap2/envs/locks/`. Both ship as package data, so
  they must keep working from an installed wheel and not only from a checkout.
- `tests/` - pytest suite.

## Guidelines
- Keep modules flat and legible - the commands must stay easy to read and edit. Avoid adding
  frameworks/abstractions that hide what runs.
- Don't change scientific defaults or swap tools without discussion; MetaWrap2 keeps the same
  algorithms and software as metaWRAP.
- Run tools only through `metawrap2.command.run(...)` so logging, provenance, and env-wrapping
  stay consistent.
- Add or update tests for any change (command-template formatting, pure helpers, provenance).
- Bins are named with `constants.bin_filename()`; don't hardcode `bin_001.fasta` or an extension.
  Anything comparing or globbing bins uses `constants.BIN_EXTENSION`, and any `-x` flag passed to
  CheckM or GTDB-Tk is derived from it - a stale `-x fa` silently scores zero bins.
- `black` formats, `isort` orders imports, `ruff` lints, `mypy` type-checks. `ruff-format` is
  deliberately not enabled: two formatters fight each other over the same files.
- `black`, `isort`, `ruff`, `mypy` and `pytest` must all pass before a push (the pre-push hooks
  run them), and CI gates the same five on Python 3.10-3.13. Formatting is *checked*, never
  applied, in CI - a job that rewrites your code and then passes teaches you nothing.
- `pre-commit` also enforces the "no .sh/.pl/.R" rule; so does CI.

## Reporting bugs
Open an issue with the version (`metawrap2 -v`), the full output, and the run's
`metawrap2.log` - it holds MetaWrap2's own messages and every line of every tool, timestamped.
`run_config.json` and `run_environment.json` from the same directory say exactly which tool
versions, conda packages and databases were used.
