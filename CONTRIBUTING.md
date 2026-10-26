# Contributing to MetaWrap2

Thanks for helping. MetaWrap2 is pure Python, with no shell, Perl or R scripts.

## Setup

```
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e ".[dev]"
pre-commit install                       # format and lint what you touch, on each commit
pre-commit install --hook-type pre-push  # run the full test gate before each push
```

## Running the tests

```
python tests/run_on_small_machine.py -q   # the whole suite, as the pre-push hook and CI run it
python -m pytest -q tests/test_x.py       # one file, while you work
metawrap2 doctor                          # check conda envs, tools and databases
```

The suite needs no bioinformatics tools and gives the same result on every machine:

- `tests/conftest.py` gives every test the same simulated machine (32 cores, 128 GB, no GPU, plenty
  of free disk) and the same `use_conda_envs` setting. If a test needs a different machine,
  monkeypatch `limits.machine_cores`, `limits.machine_memory_bytes` or `shutil.disk_usage` in that
  test. Mark a test that measures the real host with `@pytest.mark.real_machine`.
- `tests/run_on_small_machine.py` runs the suite as a 2-core, 4 GB machine with a 95%-full disk. A
  test that passes only on a large machine fails here, before the push.

Two kinds of test compare outputs with recorded results:

- Golden tests (`tests/test_bin_refinement_golden.py`) run the real refinement algorithm, with
  only CheckM replaced by a deterministic stand-in, and compare every output file with recorded
  values.
- The end-to-end golden test (`tests/test_golden_end_to_end.py`) puts stub tools on `PATH`, runs
  binning through refinement, quantification and annotation, and checks the output tree and the
  `.stats.tsv` format byte for byte.

If a change alters results, update the golden values in the same commit, so the diff shows old and
new output.

## Project layout

- `src/metawrap2/modules/<module>.py`: one file per module. The commands it runs are in the
  `COMMANDS` block and its tunable numbers in the `CONSTANTS` block, both at the top of the file.
- `src/metawrap2/command.py`: the only place that runs an external command. It handles conda env
  wrapping, logging, provenance and `--dry-run`. Modules call `metawrap2.command.run(...)`, never
  `subprocess`.
- `src/metawrap2/modules/_common.py` and the top-level modules (`config`, `provenance`, `manifest`,
  `validate`, `history`, `limits`, and so on): shared infrastructure, described in each file's
  docstring.
- `src/metawrap2/envs/<name>.yaml`: the conda environments, one per workflow, with exact-build
  lockfiles in `src/metawrap2/envs/locks/`. They ship inside the wheel, so `install-env` must work
  from an installed package as well as from a checkout.
- `benchmarks/`: the comparisons behind the defaults (`configurations.toml`), the committed
  results, and the figure scripts.
- `tests/`: the pytest suite.

## Guidelines

- Keep modules flat. Avoid abstractions that hide which commands a module runs.
- A scientific default changes only on a benchmark: add the comparison to
  `benchmarks/configurations.toml`, commit the result, and record the change in the release notes.
  Open an issue first if you are unsure.
- Add or update tests for every change.
- Name bins with `constants.bin_filename()` and match them with `constants.BIN_EXTENSION`. The `-x`
  flag passed to CheckM or GTDB-Tk is derived from it; with a wrong extension no bins are scored.
- `black` formats, `isort` orders imports, `ruff` lints and `mypy` type-checks. These four and
  `pytest` must pass before a push; the pre-push hook runs them, and CI runs them on Python 3.10 to
  3.14. CI checks formatting but does not reformat.
- The hook and CI use the same tool versions, pinned in `pyproject.toml`'s dev extras and in
  `.pre-commit-config.yaml`. A test fails if the two disagree; bump both together.
- Write documentation in the present tense. Give one reason where a choice is surprising. Avoid
  paired "X rather than Y" unless a reader would otherwise do the wrong thing.
- Every module-level import must be a declared dependency, in both `pyproject.toml` and
  `conda_pkg/meta.yaml`, which a test checks. Import an optional dependency inside the function
  that needs it.

## Releasing

A release updates two documents:

- `CHANGELOG`: one entry per change, grouped as New, Changed, Fixed, Removed and Known issue. A
  changed default states the old value and how to get it back.
- `releases/<version>.md`, for users upgrading: what the version adds, which defaults changed and
  the measurement behind each, and what is still open. Add the page to `releases/index.md` in the
  same change.

Bump the version in `pyproject.toml`, `src/metawrap2/__init__.py`, `CITATION.cff` and
`conda_pkg/meta.yaml`. A test checks they match.

Tag the release commit with an annotated tag named for the version, and push it with the commit:

```bash
git tag -a v2.1.3 -m "MetaWrap2 2.1.3"
git push origin main v2.1.3
```

## Reporting bugs

Open an issue with the version (`metawrap2 -v`), the full output, and the run's `metawrap2.log`
(MetaWrap2's messages and all tool output, timestamped). `run_config.json` and
`run_environment.json` from the same directory list the tool versions, conda packages and databases
used.
