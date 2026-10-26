"""Progress bars for MetaWrap2's own long-running loops.

External tools print their own progress. These bars cover the long loops MetaWrap2 runs itself:
filtering reads, annotating bins one at a time, reassembling each bin with SPAdes, downloading
many files.

:func:`bar` wraps ``tqdm``:

* without tqdm installed it returns the iterable unchanged (tqdm is optional);
* bars go to stderr, because several modules redirect a command's stdout into a data file;
* bars are off when stderr is not a terminal (log file, pipe, CI), so run.stderr does not fill
  with redraw lines. ``METAWRAP2_PROGRESS`` (``1``/``0``) overrides this.
"""

from __future__ import annotations

import os
import sys
from typing import Iterable, Optional

try:
    from tqdm.auto import tqdm as _tqdm
except ImportError:  # pragma: no cover - optional dependency
    _tqdm = None

__all__ = ["bar", "enabled", "write"]


def enabled() -> bool:
    """True if progress bars should be drawn."""
    override = os.environ.get("METAWRAP2_PROGRESS")
    if override is not None:
        return override.strip().lower() in ("1", "true", "yes", "on")
    if _tqdm is None:
        return False
    try:
        return bool(sys.stderr.isatty())
    except (AttributeError, ValueError):
        return False


def bar(
    iterable: Iterable,
    *,
    desc: str,
    unit: str = "it",
    total: Optional[int] = None,
    leave: bool = False,
) -> Iterable:
    """Wrap *iterable* in a tqdm progress bar, or return it unchanged if unavailable.

    ``total`` may be omitted for a streaming source; tqdm then shows a running count and rate
    instead of a percentage.
    """
    if not enabled():
        return iterable
    return _tqdm(
        iterable,
        desc=desc,
        unit=unit,
        total=total,
        leave=leave,
        file=sys.stderr,
        dynamic_ncols=True,
        smoothing=0.1,
    )


def write(message: str) -> None:
    """Print *message* without tearing an active progress bar."""
    if enabled() and _tqdm is not None:
        _tqdm.write(message, file=sys.stderr)
    else:
        sys.stderr.write(message.rstrip("\n") + "\n")
