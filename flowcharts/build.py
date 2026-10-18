#!/usr/bin/env python
"""Build every MetaWrap2 diagram: write the .drawio sources and export them with draw.io.

    python flowcharts/build.py            # all diagrams
    python flowcharts/build.py overview   # just one

Each diagram is a Python file in diagrams/. Its .drawio file is written to drawio/ and opens and
edits in draw.io; the PNG is exported to png/ by the draw.io desktop app, so it is exactly what
draw.io renders. The overview is also copied to docs/figures/, where the README and the docs site
read it.
"""

from __future__ import annotations

import importlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "diagrams"))

from drawio import export

DIAGRAMS = [
    "overview",
    "end_to_end",
    "reads_and_assembly",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "genome_catalogue",
    "function",
    "beyond_bins",
]
#: Diagrams the docs embed, copied to docs/figures/ under these names.
PUBLISHED = {"overview": "metawrap2_overview.png"}


def main(names) -> int:
    for name in names or DIAGRAMS:
        diagram = importlib.import_module(name).build()
        source = diagram.save(os.path.join(HERE, "drawio"))
        png = os.path.join(HERE, "png", name + ".png")
        export(source, png)
        print("wrote", os.path.relpath(png, os.path.dirname(HERE)))
        if name in PUBLISHED:
            shutil.copyfile(
                png, os.path.join(os.path.dirname(HERE), "docs", "figures", PUBLISHED[name])
            )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
