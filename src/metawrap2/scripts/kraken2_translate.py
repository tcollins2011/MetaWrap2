"""Translate KRAKEN2 taxid annotations into full taxonomy-string annotations.

Ported from metaWRAP's ``kraken2_translate.py``. Uses the NCBI ``names.dmp`` / ``nodes.dmp``
files inside the KRAKEN2 database to walk each taxid up to the root and build a
``;``-joined lineage. Run standalone
(``python -m metawrap2.scripts.kraken2_translate <db> <in.krak2> <out.kraken2>``) or import
:func:`translate_kraken2_annotations` in-process.
"""

from __future__ import annotations

import os
import sys


def get_full_name(taxid, names_map, ranks_map):
    """Generate the full taxonomic lineage string for a taxid."""
    taxid_lineage = list()
    while True:
        taxid_lineage.append(taxid)
        if taxid not in ranks_map:
            break
        new_taxid = ranks_map[taxid]
        if taxid == new_taxid:
            break
        else:
            taxid = new_taxid

    names_lineage = list()
    for taxid in taxid_lineage:
        name = names_map[taxid]
        names_lineage.append(name)

    return ";".join(reversed(names_lineage))


def load_kraken_db_metadata(kraken2_db):
    """Load NCBI name and parent-node maps needed to expand taxids into lineage strings."""
    print("Loading NCBI node names")
    names_path = os.path.join(kraken2_db, "taxonomy", "names.dmp")
    names_map = dict()
    with open(names_path) as handle:
        for line in handle:
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            name = cut[2]
            entry_type = cut[6]
            if entry_type == "scientific name":
                names_map[taxid] = name

    print("Loading NCBI taxonomic ranks")
    ranks_path = os.path.join(kraken2_db, "taxonomy", "nodes.dmp")
    ranks_map = dict()
    with open(ranks_path) as handle:
        for line in handle.readlines():
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            parent_taxid = cut[2]
            ranks_map[taxid] = parent_taxid
    return names_map, ranks_map


def translate_kraken2_annotations(annotation_file=None, kraken2_db=None, output=None):
    """Translate the kraken2 annotations in *annotation_file* from taxids to lineage strings."""
    print("Translating kraken2 annotations")
    if os.path.isfile(output):
        print("Looks like the translated kraken2 output already exists. Skipping...")
        return None

    names_map, ranks_map = load_kraken_db_metadata(kraken2_db)
    print("Writing translated taxonomy names to %s" % output)
    with open(output, "w") as out_handle:
        with open(annotation_file) as handle:
            for line in handle:
                cut = line.rstrip().split("\t")
                contig = cut[1]
                if cut[0] == "U":
                    taxonomy = ""
                else:
                    taxid = cut[2].split()[-1][:-1]
                    taxonomy = get_full_name(taxid, names_map, ranks_map)
                out_handle.write("%s\t%s\n" % (contig, taxonomy))


def main(argv) -> int:
    database_location = argv[0]
    kraken_file = argv[1]
    output_file = argv[2]
    print("Translating kraken2 annotations from %s, using metadata from the kraken2 database in %s; saving to %s"
          % (kraken_file, database_location, output_file))
    translate_kraken2_annotations(annotation_file=kraken_file, kraken2_db=database_location,
                                  output=output_file)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
