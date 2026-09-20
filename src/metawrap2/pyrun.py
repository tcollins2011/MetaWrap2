"""Running MetaWrap2's own Python helper scripts.

The per-module conda envs hold only that module's *external* tools (metabat2, spades,
checkm, ...). They deliberately do not contain a Python interpreter with metawrap2 and its
plotting stack installed. So MetaWrap2's own helpers - ``metawrap2.scripts.*`` and
``metawrap2.vendor.*`` - must run in the **host** interpreter that metawrap2 itself was
installed into, not inside a module env.

Use :func:`py_module` to build the argv for such a helper, and pass ``env=None`` to
``run()``. ``sys.executable`` is used rather than the string ``"python"`` so the helper
always lands on the same interpreter (and therefore the same metawrap2) that is driving
the pipeline, whatever is first on PATH.
"""

from __future__ import annotations

import shlex
import sys
from typing import List

__all__ = ["PY", "py_module"]

#: The host interpreter, shell-quoted - safe to interpolate into a command template.
PY = shlex.quote(sys.executable)


def py_module(module: str, *args: object) -> List[str]:
    """argv that runs ``python -m <module> <args...>`` in the host interpreter."""
    return [sys.executable, "-m", module] + [str(a) for a in args]
