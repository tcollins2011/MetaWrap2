"""``metawrap2 report -o REPORT.pdf RUN_DIR [RUN_DIR ...]``: one PDF for a finished run.

Reads one or more module output directories and writes a single PDF, plus the per-genome TSV that
the PDF's genome pages are a view of. Nothing is computed: see :mod:`metawrap2.report`.

    metawrap2 report -o study.pdf study/
    metawrap2 report -o study.pdf study/READ_QC study/ASSEMBLY study/BIN_REFINEMENT
    metawrap2 report -o study.pdf --title "Oral cohort, run 3" study/

Several directories are accepted because ``amr``, ``viruses`` and ``pathways`` are usually
written to their own directories next to the study.

Exits 0 whenever a PDF was written, even for a partial run. Missing modules are listed on the
report's "what was not run" page.
"""

from __future__ import annotations

import os
from typing import List

from .. import report as _report
from ..cliformat import make_parser

#: The per-genome table written beside the PDF, with the full values the PDF truncates.
MAG_TABLE_SUFFIX = ".mag_table.tsv"


def _mag_table_path(pdf_path: str) -> str:
    """``report.pdf`` -> ``report.mag_table.tsv``, beside it."""
    root, _ext = os.path.splitext(pdf_path)
    return root + MAG_TABLE_SUFFIX


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 report",
        description="Assemble one PDF report from one or more finished module output directories.",
    )
    ap.add_argument(
        "run_dirs",
        nargs="+",
        metavar="RUN_DIR",
        help="module output directories (or a study directory holding them)",
    )
    ap.add_argument("-o", "--output", required=True, help="the PDF to write")
    ap.add_argument(
        "--title", default="", help="title for the first page; defaults to a generic one"
    )
    ap.add_argument(
        "--mag-table",
        default="",
        help="where to write the per-genome TSV; defaults to the PDF's name with %s"
        % MAG_TABLE_SUFFIX,
    )
    ap.add_argument(
        "--no-mag-table",
        action="store_true",
        help="do not write the per-genome TSV (the PDF's genome pages stay truncated)",
    )
    args = ap.parse_args(argv)

    missing = [path for path in args.run_dirs if not os.path.isdir(path)]
    if missing:
        # Fail before drawing anything; a mistyped path would otherwise give a PDF that reports
        # every module as absent.
        print(
            "\nNot a directory: %s\nEach RUN_DIR must be a module output directory "
            "(the one passed to that module's -o).\n" % ", ".join(missing)
        )
        return 1

    table_path = None if args.no_mag_table else (args.mag_table or _mag_table_path(args.output))
    result = _report.build_report(
        args.run_dirs, args.output, title=args.title, mag_table_path=table_path
    )

    print("\n  %-14s %s (%d page(s))" % ("report", result.path, result.pages))
    if result.mag_table:
        print("  %-14s %s" % ("genome table", result.mag_table))
    if result.absent:
        print("\n  Not in this report, because nothing was found for it:")
        for entry in result.absent:
            print("    %-26s looked for %s" % (entry.module, entry.looked_for))
        print("\n  %s" % _report.ABSENCE_NOTE)
    print()
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
