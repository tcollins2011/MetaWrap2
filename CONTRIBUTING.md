# Contributing to MetaWrap2

Thanks for helping improve MetaWrap2! It is a pure-Python package - no Bash, Perl, or R.

## Setup
```
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e ".[dev]"
pre-commit install        # optional but recommended
```

## Running the tests
```
python -m pytest -q                 # the env-independent unit + integration tests
metawrap2 test                      # probe module conda envs + software, then run the tests
```
The unit tests don't need any bioinformatics tools; the stub-based integration test exercises
a full module run with fake tools.

## Project layout
- `src/metawrap2/modules/<module>.py` - one flat "recipe" per pipeline module. The commands a
  module runs live in the `COMMANDS` block and the tunable numbers in the `CONSTANTS` block at
  the top of the file. Keep them there and visible.
- `src/metawrap2/command.py` - the single place that runs an external command (mamba env
  wrapping, logging, provenance, dry-run). Don't call `subprocess` directly in a module.
- `src/metawrap2/{config,provenance,checkpoint,checkm}.py` and `modules/_common.py` - shared
  helpers.
- `envs/<module>.yaml` - the per-module conda env (`metawrap2-<module>`).
- `tests/` - pytest suite.

## Guidelines
- Keep modules flat and legible - the commands must stay easy to read and edit. Avoid adding
  frameworks/abstractions that hide what runs.
- Don't change scientific defaults or swap tools without discussion; MetaWrap2 keeps the same
  algorithms and software as metaWRAP.
- Run tools only through `metawrap2.command.run(...)` so logging, provenance, and env-wrapping
  stay consistent.
- Add or update tests for any change (command-template formatting, pure helpers, provenance).
- `ruff` formats/lints; `pre-commit` enforces it plus the "no .sh/.pl/.R" rule.

## Reporting bugs
Open an issue with the full output (`run.stderr` from the output dir helps) and the version
(`metawrap2 -v`).
