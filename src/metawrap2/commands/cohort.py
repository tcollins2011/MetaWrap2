"""``metawrap2 cohort``: assemble and bin every sample on its own, then dereplicate.

A thin alias over ``metawrap2 run --no-coassembly``, and deliberately nothing more: the sample
sheet, manifest, checkpointing, run lock, resource budget, history and ``--resume`` all behave
identically because it is the same planner and the same executor. The separate name exists because
"co-assemble a cohort" and "bin each sample and dereplicate" are two different analyses, and the
second one should not be reachable only by negating a flag - the steps that only exist in that mode
(``coverage``, ``dereplicate``) then come on by default, which is the part people miss.

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
