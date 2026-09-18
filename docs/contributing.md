# Contributing

See [CONTRIBUTING.md](https://github.com/ursky/MetaWrap2/blob/main/CONTRIBUTING.md) in the
repository for the full guide. In short:

```bash
pip install -e ".[dev]"
pre-commit install
python -m pytest -q
```

- MetaWrap2 is Python-only (no `.sh`/`.pl`/`.R`).
- Keep modules flat: commands live in the `COMMANDS`/`CONSTANTS` blocks at the top of each
  `src/metawrap2/modules/<module>.py`.
- Run external tools only through `metawrap2.command.run(...)`.
- Add tests for any change; `ruff` + `pre-commit` enforce style and the no-legacy-scripts rule.
