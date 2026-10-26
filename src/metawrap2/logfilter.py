"""Collapse repetitive tool output in ``run.stdout`` and ``run.stderr``.

On a 3-sample gut metagenome run, ``INITIAL_BINNING/run.stderr`` was 1.8 MB, of which 23,812 lines
were two metaBAT2 warnings repeated 11,906 times each:

    WARNING: calculated a huge mean=1.2e+03. correctedLen=456 contigDepth=789
    Please report this bam file to MetaBAT under Issue #48

Two filters are applied:

**Progress-bar redraws.** Only the final state of each ``\\r``-updated line is kept.

**Repeated messages.** Lines are grouped by a *signature* with numbers masked, so
"huge mean=1.2e+03 ... contigDepth=789" and "huge mean=9.9e+02 ... contigDepth=12" match. The first
:data:`MAX_REPEATS` occurrences are written verbatim; after that one occurrence is written every
:data:`REFRESH_SECONDS`, with its current values:

    [metawrap2] NN training: epoch=85/200, loss=2.91, top1=88.4  (collapsed; 64 similar lines
    suppressed)

The refresh is needed because masking numbers also groups progress lines: COMEBin's 200 training
epochs and CheckM's ``Finished processing 19 of 49 (38.78%) bins.`` share a signature per line.
Suppressing them outright left a nine-minute gap in the log of a working step. With the refresh,
the 11,906 metaBAT2 warnings and the 5,773 ``pass <N>: searching...<N>%`` lines from a blobology
run still collapse to about one line per 30 s each.

A message whose counter rises towards a fixed total (``epoch=20/200``, ``19 of 49``) is treated as
progress, and its last value is written when the stream ends, so the log shows how far the step got.

Some tools print nothing for hours (COMEBin training, SemiBin2 clustering), and a silent screen
looks the same as a hung process. :class:`Heartbeat` prints one line per :data:`HEARTBEAT_SECONDS`
of silence with the step's elapsed time and the tool's last line. It goes to the screen only, so
the captured log is the same with or without it.

Data output is not filtered: when a tool's stdout is the result (``stdout_path=``, e.g.
``bwa mem > x.sam``) it goes straight to the file. ``--verbose-logs`` turns filtering off.
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
#: COMEBin trains at about 27 s per epoch, so at 30 s each refresh shows a new epoch. 11,906
#: copies of one metaBAT2 warning in a 40-minute step become about 80 lines. A much shorter
#: interval brings back the noise; a much longer one makes slow steps look stalled.
REFRESH_SECONDS = 30.0

#: Seconds of silence before a heartbeat line is printed. 0 turns heartbeats off.
#:
#: Twice :data:`REFRESH_SECONDS`: a collapsed message still gets a refresh line every 30 s, so 60 s
#: with nothing shown means the tool has stopped printing. It is also long enough that a step
#: printing about one line a minute (metaSPAdes between stages, CheckM between bins) does not get a
#: heartbeat after each line.
HEARTBEAT_SECONDS = 60.0

#: Numbers, hex addresses and timings carry the variation between otherwise-identical lines.
_NUMBER_RE = re.compile(
    r"""
    0x[0-9a-fA-F]+            # hex addresses
  | \d+(?:\.\d+)?(?:[eE][-+]?\d+)?   # ints, floats, scientific notation
""",
    re.VERBOSE,
)

#: An ``n/total`` or ``n of total`` pair, as in COMEBin's ``epoch=20/200`` and CheckM's
#: ``Finished processing 19 of 49``.
_PROGRESS_RE = re.compile(r"(\d+)\s*(?:/|\bof\b)\s*(\d+)")

_PREFIX = "[metawrap2]"


def signature(line: str) -> str:
    """A form of *line* with the varying numeric parts masked, for grouping repeats."""
    return _NUMBER_RE.sub("<N>", line.strip())


def duration(seconds: float) -> str:
    """An elapsed time as ``2h05m``, ``1m15s`` or ``45s``."""
    seconds = max(0.0, seconds)
    if seconds >= 3600:
        return "%dh%02dm" % (int(seconds) // 3600, (int(seconds) % 3600) // 60)
    if seconds >= 60:
        return "%dm%02ds" % (int(seconds) // 60, int(seconds) % 60)
    return "%ds" % int(seconds)


class Heartbeat:
    """One line per interval of silence, so a working step can be told apart from a hung one.

    Fed every line the tool produces (:meth:`note`) and polled from a separate thread (:meth:`due`,
    :attr:`tick_seconds`), because the readers are blocked on the tool's pipe while it is silent.

    *clock* is injectable for tests and must be monotonic, as in :class:`LogFilter`.

    Heartbeats are printed whether or not the output is a terminal, as the collapser's refresh
    lines are. Each is a new line; a ``\\r`` rewrite would be overwritten by the next line the pumps
    write.
    """

    def __init__(
        self,
        step: str = "",
        interval: float = HEARTBEAT_SECONDS,
        enabled: bool = True,
        clock: Callable[[], float] = time.monotonic,
    ):
        self.step = step
        self.interval = interval
        #: An interval of 0 turns heartbeats off.
        self.enabled = enabled and interval > 0
        self._clock = clock
        self._started = clock()
        #: The last line the tool produced, and when. The time is None until the first line, so
        #: the heartbeat can say the tool has printed nothing.
        self._last_line = ""
        self._last_line_at: Optional[float] = None
        #: Start of the current interval: the last line or the last heartbeat. Separate from
        #: ``_last_line_at`` so a heartbeat restarts the interval but the reported silence keeps
        #: growing.
        self._counting_from = self._started

    @property
    def tick_seconds(self) -> float:
        """How often a watcher should call :meth:`due`.

        At most 5 s. Polling once per interval could delay a heartbeat by up to a second interval.
        """
        return min(self.interval, 5.0) if self.enabled else 5.0

    def note(self, line: str) -> None:
        """Record a line the tool produced. Blank lines do not count as output."""
        text = line.strip()
        if not text:
            return
        self._last_line = text
        self._last_line_at = self._clock()
        self._counting_from = self._last_line_at

    def due(self) -> Optional[str]:
        """The heartbeat line to show now, or None if the step has not been silent long enough."""
        if not self.enabled:
            return None
        now = self._clock()
        if now - self._counting_from < self.interval:
            return None
        self._counting_from = now
        subject = "%s is" % self.step if self.step else "this step is"
        if self._last_line_at is None:
            return "%s %s still running, %s elapsed, and has printed nothing at all yet\n" % (
                _PREFIX,
                subject,
                duration(now - self._started),
            )
        return "%s %s still running, %s elapsed, silent for %s; last line: %s\n" % (
            _PREFIX,
            subject,
            duration(now - self._started),
            duration(now - self._last_line_at),
            self._last_line[:160],
        )


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

    All positions are compared because the first pair is often a date, as in
    ``2026/10/15 ... epoch=21/200``. The total must match so that unrelated lines with the same
    signature are not taken for one counter.
    """
    return any(old[1] == new[1] and new[0] > old[0] for old, new in zip(before, after))


@dataclass
class _Repeat:
    """What is known about one message signature while the stream is being read."""

    count: int = 0
    #: The most recent occurrence, verbatim, written by a refresh.
    latest: str = ""
    positions: Tuple[Tuple[int, int], ...] = ()
    progress: bool = False
    #: Occurrences dropped over the whole stream (for the closing summary), and since the last
    #: line emitted for this signature (for refresh lines).
    suppressed: int = 0
    since_emit: int = 0
    last_emit: float = 0.0


class LogFilter:
    """Streaming filter: pass lines through, collapsing redraws and refreshing repeats.

    *clock* is injectable for tests and must be monotonic; a wall clock stepped by NTP can make
    the interval negative or very large.

    *heartbeat* is told about every line. A command's stdout and stderr filters share one
    :class:`Heartbeat`, so output on either stream counts, and a quiet minute gives one heartbeat,
    not two.
    """

    def __init__(
        self,
        max_repeats: int = MAX_REPEATS,
        enabled: bool = True,
        refresh_seconds: float = REFRESH_SECONDS,
        clock: Callable[[], float] = time.monotonic,
        heartbeat: Optional[Heartbeat] = None,
    ):
        self.max_repeats = max_repeats
        self.enabled = enabled
        self.refresh_seconds = refresh_seconds
        self._clock = clock
        self._heartbeat = heartbeat
        self._seen: Dict[str, _Repeat] = {}

    def feed(self, line: str) -> Optional[str]:
        """Return the line to write, or None to suppress it."""
        # Before the --verbose-logs check: heartbeats are not filtering, so that flag leaves
        # them on.
        if self._heartbeat is not None:
            self._heartbeat.note(line)
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
                # %g so a sub-second interval is not shown as "every 0s".
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
                # For a counter, write the last value, which shows how far the step got.
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
        """(distinct messages, total suppressed lines), for tests and reporting."""
        return len(self._seen), sum(state.suppressed for state in self._seen.values())
