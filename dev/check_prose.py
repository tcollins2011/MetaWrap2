#!/usr/bin/env python
"""Warn about writing habits that crept into this repository's documentation.

The docs were edited down once to plain, direct prose. These are the patterns that made them hard to
read, counted per file. The check only warns and always exits 0: a contrast or an em-dash is
sometimes the right choice, and a human decides.

    python dev/check_prose.py README.md docs/*.md       # or let the pre-commit hook pass the files
"""

from __future__ import annotations

import re
import sys
from typing import Dict, List, Tuple

#: (name, pattern, how many per 1,000 words before it is worth a warning)
PATTERNS: List[Tuple[str, str, float]] = [
    ("'rather than' / ', not '", r"\brather than\b|, not (?:a |an |the )?\w+", 2.0),
    ("history in present-tense docs ('used to', 'no longer')", r"\bused to\b|\bno longer\b", 1.0),
    (
        "dramatic words",
        r"\b(?:silently|quietly|deliberately|on purpose|genuinely|honest(?:ly)?|load-bearing)\b",
        1.0,
    ),
    # Real em-dashes only: " - " is also the CHANGELOG's own "New - ..." separator.
    ("em-dashes", r"—", 4.0),
    # Bold mid-sentence. A bold label opening a line or a list item is a heading, not emphasis.
    (
        "bold emphasis mid-sentence",
        r"(?m)(?<!^)(?<!^- )(?<!^\* )(?<!^- CODE )(?<![|*])\*\*[^*|\n]+\*\*(?![|*])",
        3.0,
    ),
]

#: Release notes and the changelog are where history belongs, so the history rule skips them.
HISTORY_OK = ("CHANGELOG", "docs/releases/")


def _strip_code(text: str) -> str:
    """Drop fenced code blocks and inline code, where none of this applies."""
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    # Inline code becomes a placeholder word, so a list item that opens with a code tag and then a
    # bold label ("- `[low]` **Name**") still reads as a label.
    return re.sub(r"`[^`\n]*`", "CODE", text)


def check(path: str) -> List[str]:
    with open(path, encoding="utf-8") as fh:
        prose = _strip_code(fh.read())
    words = max(1, len(prose.split()))
    warnings = []
    for name, pattern, per_thousand in PATTERNS:
        if name.startswith("history") and any(part in path for part in HISTORY_OK):
            continue
        found = len(re.findall(pattern, prose))
        rate = 1000.0 * found / words
        if rate > per_thousand:
            warnings.append(
                "  %-55s %3d  (%.1f per 1k words; aim for <= %g)"
                % (name, found, rate, per_thousand)
            )
    return warnings


def main(argv: List[str]) -> int:
    report: Dict[str, List[str]] = {}
    for path in argv:
        if path.endswith((".md", "CHANGELOG")):
            found = check(path)
            if found:
                report[path] = found
    for path, lines in report.items():
        print("%s:" % path)
        print("\n".join(lines))
    if report:
        print("\n(prose check: advisory only, nothing is blocked; see CONTRIBUTING.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
