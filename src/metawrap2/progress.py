"""Progress bars for MetaWrap2's own long-running loops.

MetaWrap2 spends most of its wall time inside external tools, which print their own progress.
The loops worth instrumenting are the ones MetaWrap2 runs itself and that used to look like a
hang: streaming millions of reads through a filter, annotating bins one at a time, reassembling
each bin with SPAdes, downloading hundreds of files.

:func:`bar` wraps ``tqdm`` and degrades gracefully:

* if tqdm is not installed, it returns the iterable unchanged, so tqdm stays an optional
  dependency and nothing breaks without it;
* bars are written to **stderr**, never stdout, because several modules redirect a command's
  stdout into a data file and a progress bar in a FASTA would corrupt it;
* bars are disabled automatically when stderr is not a terminal (a log file, a pipe, a CI
  run), so run.stderr does not fill up with thousands of redraw lines. ``METAWRAP2_PROGRESS``
  overrides the decision either way (``1``/``0``).
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
    rather than a percentage, which is still the useful signal for "is this moving?".
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
