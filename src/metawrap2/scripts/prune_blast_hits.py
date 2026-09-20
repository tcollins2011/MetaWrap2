#!/usr/bin/env python
"""Prune MEGABLAST hits down to a single ranked taxon id per line.

Reads the NCBI ``nodes.dmp`` to learn each tax id's rank, then rewrites the raw blast
table: for each hit's semicolon-separated staxids (column 6, 0-indexed 5) it keeps only the
first id that has one of the standard ranks (species..superkingdom) and drops lines with no
ranked id. Column 6 is collapsed to that single id.

Ported faithfully from metaWRAP's ``prune_blast_hits.py``.
"""

from __future__ import annotations

import sys
from typing import Dict, Iterator, List

# "domain" is included alongside "superkingdom" because NCBI renamed that rank: a current
# taxdump has no "superkingdom" nodes at all, so a hit whose only ranked ancestor is a domain
# would otherwise be discarded.
INCLUDE = {"species", "genus", "family", "order", "class", "phylum", "superkingdom", "domain"}


def load_ranks(nodes_dmp: str) -> Dict[str, str]:
    """Map tax id -> rank from an NCBI ``nodes.dmp`` file."""
    ranks = {}
    with open(nodes_dmp) as fh:
        for line in fh:
            cut = line.split("\t")
            ranks[cut[0]] = cut[4]
    return ranks


def prune(nodes_dmp: str, raw_tab: str) -> Iterator[str]:
    """Yield pruned blast lines (tab-joined, no trailing newline)."""
    ranks = load_ranks(nodes_dmp)
    with open(raw_tab) as fh:
        for line in fh:
            cut = [c.strip() for c in line.strip().split("\t")]
            ids = cut[5]
            if len(ids.split(";")) < 1:
                continue
            ct = 0
            for tax_id in ids.split(";"):
                if tax_id not in ranks:
                    continue
                if ranks[tax_id] not in INCLUDE:
                    continue
                if ct > 0:
                    continue
                cut[5] = tax_id
                ct += 1
                yield "\t".join(cut)


def main(argv: List[str]) -> int:
    for line in prune(argv[0], argv[1]):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
