"""The banner box: wrap the text, then pad it - never the other way round.

Every user-facing message MetaWrap2 prints goes through `format_comment`, and the one a first
install reads last is quickstart's closing banner. Padding each line to the box width *before*
wrapping measured a multi-line segment whole, so the right border was positioned from a length
that no printed line had, and landed on the line after the break. Three of that banner's four
lines came out malformed.

The invariant these tests pin is small and checkable on any message: every line between the
borders is exactly the box width, framed by the delimiter on both sides, and no word of the
message has been cut in two.
"""

from __future__ import annotations

from metawrap2 import logging as mwlog
from metawrap2.constants import BANNER_WIDTH

#: The closing banner of the install that reported this, verbatim - the message, its embedded
#: newlines, and all. It is the regression case because its newlines are what broke the old
#: rendering: `len(segment)` counted past them and placed the border off the end of the line.
QUICKSTART_CLOSING_BANNER = (
    "Some of it is not ready - the table above says which.\n"
    "`metawrap2 doctor --fix` recreates exactly the environments that are broken and\n"
    "leaves the healthy ones alone. Everything that did install is still installed."
)

LONG_MESSAGE = (
    "Creating the conda environments is the slow part of a first install, and it is a quiet one "
    "because every solver writes to its own log file rather than to this screen. Eleven "
    "environments can take half an hour, so a long "
    "silence here is the tool working rather than the tool hung. Nothing that already exists is "
    "rebuilt, which is what makes a second attempt cheap."
)


def content_lines(rendered: str) -> list[str]:
    """The lines inside the box: no leading blank, no top or bottom border, no trailing blank."""
    lines = rendered.split("\n")
    return lines[2:-2]


def test_every_line_of_a_long_message_is_exactly_the_box_width():
    """Regression: a wrapped line used to be padded to the width the *unwrapped* segment had."""
    rendered = mwlog.format_comment(LONG_MESSAGE, "-")
    body = content_lines(rendered)
    assert len(body) > 1, "a message this long has to wrap onto several lines"
    for line in body:
        assert len(line) == BANNER_WIDTH, line


def test_every_line_is_framed_by_the_delimiter_on_both_sides():
    """The right border landing on the following line is what made the banner unreadable."""
    for line in content_lines(mwlog.format_comment(LONG_MESSAGE, "*")):
        assert line.startswith("*****"), line
        assert line.endswith("*****"), line


def test_no_word_is_split_across_lines():
    """Breaking a word would break a path or a command a user is being told to copy."""
    rendered = mwlog.format_comment(LONG_MESSAGE, "-")
    words = LONG_MESSAGE.split()
    printed = " ".join(line.strip("-").strip() for line in content_lines(rendered)).split()
    assert printed == words


def test_the_real_quickstart_closing_banner_renders_as_a_box():
    """The verbatim message from the install that reported this: three of four lines were broken.

    It carries explicit newlines, which is what the old order of operations mishandled - so this
    asserts both halves at once: the box is intact, and the author's own line breaks are kept
    rather than being wrapped away into one paragraph.
    """
    rendered = mwlog.format_comment(QUICKSTART_CLOSING_BANNER, "*")
    body = content_lines(rendered)
    assert len(body) == 3, "the message's three authored lines should stay three lines"
    for line in body:
        assert len(line) == BANNER_WIDTH, line
        assert line.startswith("*****") and line.endswith("*****"), line
    assert "says which." in body[0]
    assert "`metawrap2 doctor --fix` recreates" in body[1]
    assert "still installed." in body[2]


def test_a_blank_line_in_a_message_stays_a_full_width_line():
    """A paragraph break must not escape the box: an unpadded blank line has no borders at all."""
    for line in content_lines(mwlog.format_comment("first paragraph\n\nsecond paragraph", "#")):
        assert len(line) == BANNER_WIDTH, line
        assert line.startswith("#####") and line.endswith("#####"), line


def test_a_token_longer_than_the_box_is_not_cut_in_half():
    """Long tokens here are paths and commands; half of one pasted into a shell is worse than ugly."""
    url = "https://example.invalid/" + "x" * 150
    rendered = mwlog.format_comment("fetch it from %s now" % url, "-")
    assert url in rendered, "the token must survive whole, even if it overruns the border"
