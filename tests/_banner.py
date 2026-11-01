"""Remove banner borders and line wrapping from captured output, for substring checks.

Messages are wrapped to the banner width and framed with five delimiters on each side (see
metawrap2.logging.format_comment). Where a line breaks depends on the text, and a long tmp_path
moves the break, so a phrase can be split across two framed lines.
"""

import re

_EDGE = re.compile(r"^([-*#=])\1{4}|([-*#=])\2{4}$")


def unbanner(text: str) -> str:
    """*text* with border rows and side borders dropped and whitespace collapsed to single spaces."""
    kept = []
    for raw in text.split("\n"):
        line = raw.strip()
        if line and set(line) <= set("-*#="):
            continue  # a top or bottom border
        kept.append(_EDGE.sub("", line))
    return " ".join(" ".join(kept).split())
