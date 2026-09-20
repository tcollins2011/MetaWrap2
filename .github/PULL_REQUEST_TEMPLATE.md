## What this changes

## Why

## Checklist
- [ ] The full gate passes: `black --check . && isort --check-only . && ruff check . && mypy && python -m pytest -q`
- [ ] Added/updated tests for the change
- [ ] Kept modules flat: commands/constants stay in the `COMMANDS`/`CONSTANTS` blocks
- [ ] External tools are run via `metawrap2.command.run(...)` (not `subprocess` directly)
- [ ] No `.sh`/`.pl`/`.R` files added (MetaWrap2 is Python-only)
- [ ] Did not change scientific defaults / swap tools without discussion
- [ ] Golden tests still pass, or their recorded values were updated deliberately in this PR
- [ ] Bin names/extensions go through `constants.bin_filename()` / `constants.BIN_EXTENSION`
- [ ] Updated docs / CHANGELOG if user-facing
