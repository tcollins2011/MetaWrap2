"""Collapsing the repetitive noise bioinformatics tools emit into their logs.

``run.stdout`` / ``run.stderr`` are meant to be the thing you open when something went wrong.
They stop being that when a tool repeats itself tens of thousands of times. Measured on a
real 3-sample gut metagenome run, ``INITIAL_BINNING/run.stderr`` was 1.8 MB, of which 23,812
lines were two metaBAT2 warnings repeated 11,906 times each:

    WARNING: calculated a huge mean=1.2e+03. correctedLen=456 contigDepth=789
    Please report this bam file to MetaBAT under Issue #48

The signal (which contigs, which sample, what actually failed) is buried. Two filters fix
this without losing anything a user needs:

**Progress-bar redraws.** A ``\\r``-updated progress line is one line on screen but hundreds
of thousands of bytes in a file. Only the final state of each redraw is kept.

**Repeated messages.** Lines are reduced to a *signature* with numbers masked out, so
"huge mean=1.2e+03 ... contigDepth=789" and "huge mean=9.9e+02 ... contigDepth=12" share one
signature. The first :data:`MAX_REPEATS` occurrences are written verbatim; after that they are
counted, and a single summary line is emitted at the end:

    [metawrap2] ... and 11856 more lines like: WARNING: calculated a huge mean=<N>. ...

Nothing here touches a command's *data* output: when a tool's stdout is the result
(``stdout_path=``, e.g. ``bwa mem > x.sam``) it is written straight to the file and never
passes through a filter. Filtering applies only to the log/screen stream, and
``--verbose-logs`` turns it off entirely.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

#: How many times a given message is written verbatim before it starts being counted.
MAX_REPEATS = 20

#: Numbers, hex addresses and timings carry the variation between otherwise-identical lines.
_NUMBER_RE = re.compile(
    r"""
    0x[0-9a-fA-F]+            # hex addresses
  | \d+(?:\.\d+)?(?:[eE][-+]?\d+)?   # ints, floats, scientific notation
""",
    re.VERBOSE,
)

_PREFIX = "[metawrap2]"


def signature(line: str) -> str:
    """A form of *line* with the varying numeric parts masked, for grouping repeats."""
    return _NUMBER_RE.sub("<N>", line.strip())


def collapse_redraws(text: str) -> str:
    """Keep only the final state of each ``\\r``-updated progress line."""
    if "\r" not in text:
        return text
    # A trailing \r means the writer intends to overwrite this line next; keep the last
    # non-empty segment so the final rendered state survives.
    segments = [seg for seg in text.split("\r") if seg.strip()]
    return segments[-1] if segments else ""


class LogFilter:
    """Streaming filter: pass lines through, collapsing redraws and capping repeats."""

    def __init__(self, max_repeats: int = MAX_REPEATS, enabled: bool = True):
        self.max_repeats = max_repeats
        self.enabled = enabled
        self._counts: Dict[str, int] = {}
        self._examples: Dict[str, str] = {}

    def feed(self, line: str) -> Optional[str]:
        """Return the line to write, or None to suppress it."""
        if not self.enabled:
            return line
        if "\r" in line:
            stripped = line.rstrip("\n")
            collapsed = collapse_redraws(stripped)
            if not collapsed:
                return None
            line = collapsed + ("\n" if line.endswith("\n") else "")
        if not line.strip():
            return line

        sig = signature(line)
        count = self._counts.get(sig, 0) + 1
        self._counts[sig] = count
        if count <= self.max_repeats:
            return line
        if count == self.max_repeats + 1:
            self._examples[sig] = line.strip()
            return (
                "%s further occurrences of this message are being counted, not logged: "
                "%s\n" % (_PREFIX, sig[:160])
            )
        return None

    def summary(self) -> List[str]:
        """Closing lines summarising everything that was suppressed."""
        if not self.enabled:
            return []
        out: List[str] = []
        for sig, count in sorted(self._counts.items(), key=lambda kv: -kv[1]):
            hidden = count - self.max_repeats - 1
            if hidden > 0:
                out.append(
                    "%s suppressed %d further occurrences of: %s\n" % (_PREFIX, hidden, sig[:160])
                )
        return out

    def stats(self) -> Tuple[int, int]:
        """(distinct messages, total suppressed lines) - for tests and reporting."""
        total = sum(max(0, c - self.max_repeats - 1) for c in self._counts.values())
        return len(self._counts), total
