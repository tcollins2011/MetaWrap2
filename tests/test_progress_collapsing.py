"""The log filter collapses noise without hiding progress.

Both behaviours come from one heuristic (numbers masked out of the line), so a change that
fixes one of these is capable of breaking the other, and both are driven here off an injected
clock rather than real time.
"""

from itertools import pairwise

from metawrap2.logfilter import LogFilter

#: A COMEBin training line, one per epoch, 200 of them per run.
EPOCH = "2026-10-15 00:41:29,011 - NN training: epoch=%d/200, batches=27, loss=3.45, top1=83.06"

#: A blobology search line: 5,773 of these came out of one real run, and they say nothing.
NOISE = "pass %d: searching...%d%%"


class FakeClock:
    """A monotonic clock the test moves by hand, so the throttle is tested without sleeping."""

    def __init__(self) -> None:
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _feed(filt: LogFilter, lines: list) -> list:
    """Every line the filter chose to write, with the line endings stripped."""
    return [out.rstrip("\n") for out in (filt.feed(line) for line in lines) if out is not None]


def test_rising_counter_is_refreshed_not_suppressed():
    """A progress line stays visible. If this breaks, a working step looks hung for minutes."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    emitted_at = []
    for epoch in range(1, 101):
        out = filt.feed(EPOCH % epoch + "\n")
        if out is not None:
            emitted_at.append((clock.now, out))
        clock.advance(27.0)  # COMEBin's real per-epoch time
    assert any("collapsed;" in line for _at, line in emitted_at), "all epochs past 20 hidden"
    gaps = [later - earlier for (earlier, _a), (later, _b) in pairwise(emitted_at)]
    # A line that only arrives every 27 s cannot be refreshed at exactly 30 s, so two epochs is
    # the worst silence this can produce - a tenth of the nine-minute gap that was reported.
    assert max(gaps) <= 54.0
    # 20 verbatim, the notice, then a refresh every other epoch: still fewer lines than epochs.
    assert len(emitted_at) < 100


def test_refresh_carries_the_latest_values():
    """The refresh shows where the step is now. If it breaks, the line is stale and useless."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    for epoch in range(1, 86):
        filt.feed(EPOCH % epoch + "\n")
    clock.advance(31.0)
    refreshed = filt.feed(EPOCH % 86 + "\n")
    assert refreshed is not None
    assert "epoch=86/200" in refreshed
    assert "collapsed; 64 similar lines suppressed" in refreshed  # epochs 22-85


def test_noisy_message_is_still_collapsed():
    """The 5,773-line case stays collapsed. If this breaks, run.stderr is unreadable again."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    written = []
    for i in range(1, 5774):
        out = filt.feed(NOISE % (i, i % 100) + "\n")
        clock.advance(0.1)  # 5,773 lines over ~10 minutes
        if out is not None:
            written.append(out)
    # 20 verbatim + the notice + one refresh per 30 s of the 577 s stream.
    assert len(written) <= 45
    # Every occurrence is accounted for: emitted or counted, none of them lost.
    assert filt.stats() == (1, 5773 - len(written))


def test_suppressed_count_is_reported_when_the_stream_ends():
    """The closing summary still totals what was hidden, so a quiet log is never a silent one."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    for i in range(1, 501):
        filt.feed(NOISE % (i, i % 100) + "\n")
    summary = "".join(filt.summary())
    assert "suppressed 479 further occurrences of: pass <N>: searching...<N>%" in summary


def test_summary_keeps_the_last_value_of_a_progress_message():
    """A step killed mid-progress still says how far it got, not just how much was hidden."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    for epoch in range(1, 86):  # no clock advance: every one after the cap is suppressed
        filt.feed(EPOCH % epoch + "\n")
    summary = "".join(filt.summary())
    assert "epoch=85/200" in summary
    assert "64 similar lines suppressed in total" in summary


def test_refresh_interval_is_respected():
    """Nothing is emitted between refreshes. If this breaks, the filter stops filtering."""
    clock = FakeClock()
    filt = LogFilter(refresh_seconds=30.0, clock=clock)
    _feed(filt, [NOISE % (i, i) + "\n" for i in range(1, 22)])  # 20 verbatim + the notice
    assert filt.feed(NOISE % (22, 22) + "\n") is None
    clock.advance(29.9)
    assert filt.feed(NOISE % (23, 23) + "\n") is None
    clock.advance(0.2)
    assert filt.feed(NOISE % (24, 24) + "\n") is not None


def test_a_date_before_the_counter_does_not_hide_it():
    """Progress is found anywhere in the line, not only in its first n/total pair."""
    clock = FakeClock()
    filt = LogFilter(clock=clock)
    for done in range(1, 50):
        filt.feed("2026/10/15 Finished processing %d of 49 (%.2f%%) bins.\n" % (done, done * 2.04))
    summary = "".join(filt.summary())
    assert "Finished processing 49 of 49" in summary


def test_verbose_logs_passes_everything_through():
    """--verbose-logs stays an exact escape hatch: no collapsing, no refreshes, no summary."""
    filt = LogFilter(enabled=False, clock=FakeClock())
    lines = [EPOCH % epoch + "\n" for epoch in range(1, 201)]
    assert _feed(filt, lines) == [line.rstrip("\n") for line in lines]
    assert filt.summary() == []
