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
signature. The first :data:`MAX_REPEATS` occurrences are written verbatim; after that the
message is *collapsed but refreshed* - one occurrence every :data:`REFRESH_SECONDS`, carrying
whatever its numbers say at that moment:

    [metawrap2] NN training: epoch=85/200, loss=2.91, top1=88.4  (collapsed; 64 similar lines
    suppressed)

The refresh is what makes the signature trick safe. Masking numbers cannot tell noise from
state: COMEBin's 200 training epochs are one log line each, and after number-stripping they
are indistinguishable from metaBAT2's warning. Suppressing them outright left a nine-minute
gap in the log of a step that was working perfectly, with no uncollapsed copy anywhere - the
only way to learn the run had reached epoch 85 of 200 was to find COMEBin's own
``work_files/comebin_out/comebin_res/training.log``. CheckM's
``Finished processing 19 of 49 (38.78%) bins.`` had the same problem. Refreshing costs one
line per 30 s per message, which does not undo the saving: the 11,906 metaBAT2 warnings and
the 5,773 ``pass <N>: searching...<N>%`` lines from a blobology run are still collapsed to
roughly one line per 30 s each.

A message whose numbers are advancing towards a total (``epoch=20/200``, ``19 of 49``) is also
recognised as progress, and its *last* value - not just a count of what was dropped - is
written when the stream ends, so the log says how far the step got even if it died seconds
after a refresh.

Nothing here touches a command's *data* output: when a tool's stdout is the result
(``stdout_path=``, e.g. ``bwa mem > x.sam``) it is written straight to the file and never
passes through a filter. Filtering applies only to the log/screen stream, and
``--verbose-logs`` turns it off entirely.
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

#: How many times a given message is written verbatim before it starts being counted.
MAX_REPEATS = 20

#: Seconds between refreshes of a collapsed message.
#:
#: The interval has to satisfy two things at once. A user watching a long step must not mistake
#: it for hung: COMEBin trains for 200 epochs at roughly 27 s each, so at 30 s every refresh
#: moves the epoch counter, and the silence between refreshes is never long enough to look like
#: a stall (the defect this replaced left nine minutes of nothing). And a pathological repeat
#: must stay cheap: 11,906 copies of one metaBAT2 warning in a 40-minute step come to about 80
#: refresh lines rather than 11,906, which is the same order of saving as suppressing them
#: entirely. Anything much shorter starts to reproduce the noise; much longer and a step that
#: emits one line a minute looks dead between refreshes.
REFRESH_SECONDS = 30.0

#: Numbers, hex addresses and timings carry the variation between otherwise-identical lines.
_NUMBER_RE = re.compile(
    r"""
    0x[0-9a-fA-F]+            # hex addresses
  | \d+(?:\.\d+)?(?:[eE][-+]?\d+)?   # ints, floats, scientific notation
""",
    re.VERBOSE,
)

#: An ``n/total`` or ``n of total`` pair - the shape a counted-progress line has in every tool
#: seen so far (COMEBin's ``epoch=20/200``, CheckM's ``Finished processing 19 of 49``).
_PROGRESS_RE = re.compile(r"(\d+)\s*(?:/|\bof\b)\s*(\d+)")

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


def counter_positions(line: str) -> Tuple[Tuple[int, int], ...]:
    """Every ``n/total`` pair in *line*, in the order they appear."""
    return tuple((int(done), int(total)) for done, total in _PROGRESS_RE.findall(line))


def _is_rising(before: Tuple[Tuple[int, int], ...], after: Tuple[Tuple[int, int], ...]) -> bool:
    """True if some counter kept its total and raised its count between two occurrences.

    Compared position by position rather than only on the first pair, because the leading
    numbers of a log line are often a slash-formatted date or a time, which never moves - and
    on ``2026/10/15 ... epoch=21/200`` the first pair would hide the one that matters. The total
    has to match so that two unrelated lines sharing a signature are not read as one counter.
    """
    return any(old[1] == new[1] and new[0] > old[0] for old, new in zip(before, after))


@dataclass
class _Repeat:
    """What is known about one message signature while the stream is being read."""

    count: int = 0
    #: The most recent occurrence, verbatim - this is what a refresh carries.
    latest: str = ""
    positions: Tuple[Tuple[int, int], ...] = ()
    progress: bool = False
    #: Occurrences dropped over the whole stream, and since the last line emitted for this
    #: signature. The first is the honest total for the closing summary; the second is what a
    #: refresh line reports, since the total at that point says nothing about the gap.
    suppressed: int = 0
    since_emit: int = 0
    last_emit: float = 0.0


class LogFilter:
    """Streaming filter: pass lines through, collapsing redraws and refreshing repeats.

    *clock* is injected so the throttle can be tested without sleeping; it must be monotonic,
    because a wall clock stepped by NTP or a DST change can make the interval negative (no
    refresh ever again) or enormous (a refresh on every line).
    """

    def __init__(
        self,
        max_repeats: int = MAX_REPEATS,
        enabled: bool = True,
        refresh_seconds: float = REFRESH_SECONDS,
        clock: Callable[[], float] = time.monotonic,
    ):
        self.max_repeats = max_repeats
        self.enabled = enabled
        self.refresh_seconds = refresh_seconds
        self._clock = clock
        self._seen: Dict[str, _Repeat] = {}

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
        state = self._seen.setdefault(sig, _Repeat())
        state.count += 1
        state.latest = line.strip()
        self._note_counter(state, line)

        if state.count <= self.max_repeats:
            return self._emit(state, line)
        if state.count == self.max_repeats + 1:
            return self._emit(
                state,
                # %g, not %d: a sub-second interval (the tests, and anything tuning this) would
                # otherwise be announced as "every 0s".
                "%s further occurrences of this message are collapsed to one line every %gs, "
                "carrying its latest values: %s\n" % (_PREFIX, self.refresh_seconds, sig[:160]),
            )
        if self._clock() - state.last_emit < self.refresh_seconds:
            state.suppressed += 1
            state.since_emit += 1
            return None
        return self._emit(
            state,
            "%s %s  (collapsed; %d similar lines suppressed)\n"
            % (_PREFIX, state.latest, state.since_emit),
        )

    def _emit(self, state: _Repeat, text: str) -> str:
        """Record that something was written for this signature, and return it."""
        state.last_emit = self._clock()
        state.since_emit = 0
        return text

    def _note_counter(self, state: _Repeat, line: str) -> None:
        """Mark this signature as progress once its counter is seen to advance."""
        positions = counter_positions(line)
        if not state.progress and _is_rising(state.positions, positions):
            state.progress = True
        state.positions = positions

    def summary(self) -> List[str]:
        """Closing lines summarising everything that was suppressed."""
        if not self.enabled:
            return []
        out: List[str] = []
        for sig, state in sorted(self._seen.items(), key=lambda kv: -kv[1].suppressed):
            if not state.suppressed:
                continue
            if state.progress:
                # For a counter, the last value is worth more than the signature: this is the
                # line that says how far the step got, and a step killed a few seconds after a
                # refresh would otherwise end its log 29 s short of what it actually reached.
                out.append(
                    "%s %s  (collapsed; %d similar lines suppressed in total)\n"
                    % (_PREFIX, state.latest, state.suppressed)
                )
            else:
                out.append(
                    "%s suppressed %d further occurrences of: %s\n"
                    % (_PREFIX, state.suppressed, sig[:160])
                )
        return out

    def stats(self) -> Tuple[int, int]:
        """(distinct messages, total suppressed lines) - for tests and reporting."""
        return len(self._seen), sum(state.suppressed for state in self._seen.values())
