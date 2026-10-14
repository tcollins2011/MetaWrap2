# Contributing

The full guide is [CONTRIBUTING.md](https://github.com/ursky/MetaWrap2/blob/main/CONTRIBUTING.md)
in the repository. In short:

```bash
pip install -e ".[dev]"
pre-commit install && pre-commit install --hook-type pre-push
python tests/run_on_small_machine.py -q
```

- MetaWrap2 is Python only, with no `.sh`, `.pl` or `.R` files.
- Each module is one file, with its commands in the `COMMANDS` and `CONSTANTS` blocks at the top of
  `src/metawrap2/modules/<module>.py`.
- External tools run only through `metawrap2.command.run(...)`.
- Every change comes with tests. A default changes only on a committed benchmark.
