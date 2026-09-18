## What this changes

## Why

## Checklist
- [ ] Tests pass (`python -m pytest -q`)
- [ ] Added/updated tests for the change
- [ ] Kept modules flat: commands/constants stay in the `COMMANDS`/`CONSTANTS` blocks
- [ ] External tools are run via `metawrap2.command.run(...)` (not `subprocess` directly)
- [ ] No `.sh`/`.pl`/`.R` files added (MetaWrap2 is Python-only)
- [ ] Did not change scientific defaults / swap tools without discussion
- [ ] Updated docs / CHANGELOG if user-facing
