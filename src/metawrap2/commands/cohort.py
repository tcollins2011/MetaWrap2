"""``metawrap2 cohort``: assemble and bin every sample on its own, then dereplicate.

An alias for ``metawrap2 run --no-coassembly``, using the same planner and executor, so the
sample sheet, manifest, checkpointing, run lock, resource budget, history and ``--resume`` behave
the same. Per-sample assembly is a different analysis from co-assembly, so it has its own command;
the steps specific to it (``coverage``, ``dereplicate``) are on by default.

Every option is ``run``'s; ``metawrap2 cohort -h`` lists them.
"""

from __future__ import annotations

from typing import List

from . import run as _run


def main(argv: List[str]) -> int:
    return _run.main(argv, prog="metawrap2 cohort", cohort=True)


if __name__ == "__main__":  # pragma: no cover - module entry point
    import sys

    sys.exit(main(sys.argv[1:]))
