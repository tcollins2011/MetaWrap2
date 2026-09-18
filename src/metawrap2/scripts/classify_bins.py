#!/usr/bin/env python
"""Call a consensus taxonomy for each bin from its contigs' per-contig taxonomy.

For every bin, a weighted taxonomy tree is built (each contig contributes its length to its
lineage) and traversed from the root: at each level the most-supported child is chosen only
if it holds more than 50% of the weight, otherwise traversal stops. The resulting lineage is
the bin's consensus taxonomy.

Ported faithfully from metaWRAP's ``classify_bins.py`` helper. ``.gz`` bins are supported via
the shared sequence reader.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterator, List, Tuple

from ..io.seqio import iter_fasta


def add_to_tree(tree: dict, tax_list: List[str], length: int) -> dict:
    """Add a lineage *tax_list* carrying *length* weight into the nested *tree*."""
    if len(tax_list) == 0:
        return tree
    if tax_list[0] not in tree:
        tree[tax_list[0]] = [length, {}]
    else:
        tree[tax_list[0]][0] += length
    add_to_tree(tree[tax_list[0]][1], tax_list[1:], length)
    return tree


def traverse(tree: dict, taxonomy: List[str], weight: int) -> List[str]:
    """Descend *tree* along the >50%-weight path, accumulating the lineage."""
    if len(tree) == 0:
        return taxonomy
    total_score = 0
    max_score = 0
    max_class = ""
    for k in tree:
        total_score += tree[k][0]
        if tree[k][0] > max_score:
            max_score = tree[k][0]
            max_class = k
    if weight != 0:
        total_score = weight

    if 100 * max_score // total_score > 50:
        taxonomy.append(max_class)
        taxonomy = traverse(tree[max_class][1], taxonomy, tree[max_class][0])
    else:
        return taxonomy
    return taxonomy


def load_contig_taxonomy(path: str) -> Dict[str, str]:
    """Load contig -> taxonomy-path (semicolon-delimited) from a two-column table."""
    taxonomy = {}
    with open(path) as fh:
        for line in fh:
            cut = line.strip().split("\t")
            if len(cut) < 2:
                continue
            taxonomy[cut[0]] = cut[1]
    return taxonomy


def consensus_taxonomy(contig_taxonomy_tab: str, bin_folder: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(bin_filename, consensus_lineage)`` for each bin in *bin_folder*."""
    taxonomy = load_contig_taxonomy(contig_taxonomy_tab)
    for filename in os.listdir(bin_folder):
        tax_tree: dict = {}
        for header, seq in iter_fasta(os.path.join(bin_folder, filename)):
            if header in taxonomy:
                tax_tree = add_to_tree(tax_tree, taxonomy[header].split(";"), len(seq))
        consensus = traverse(tax_tree, [], 0)
        yield filename, ";".join(consensus)


def main(argv: List[str]) -> int:
    for filename, consensus in consensus_taxonomy(argv[0], argv[1]):
        print(filename + "\t" + consensus)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
