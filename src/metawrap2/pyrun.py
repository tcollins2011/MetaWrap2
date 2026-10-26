"""Running MetaWrap2's own Python helper scripts.

The per-module conda envs hold only that module's external tools (metabat2, spades, checkm, ...)
and have no Python with metawrap2 and its plotting stack. MetaWrap2's own helpers
(``metawrap2.scripts.*`` and ``metawrap2.vendor.*``) therefore run in the host interpreter that
metawrap2 was installed into.

Build the argv with :func:`py_module` and pass ``env=None`` to ``run()``. ``sys.executable`` is
used instead of ``"python"`` so the helper runs on the same interpreter as the pipeline, whatever
is first on PATH.
"""

from __future__ import annotations

import shlex
import sys
from typing import List

__all__ = ["PY", "py_module"]

#: The host interpreter, shell-quoted for use in a command template.
PY = shlex.quote(sys.executable)


def py_module(module: str, *args: object) -> List[str]:
    """argv that runs ``python -m <module> <args...>`` in the host interpreter."""
    return [sys.executable, "-m", module] + [str(a) for a in args]
