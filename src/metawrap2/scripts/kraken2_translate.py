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
    taxid_lineage = []
    while True:
        taxid_lineage.append(taxid)
        if taxid not in ranks_map:
            break
        new_taxid = ranks_map[taxid]
        if taxid == new_taxid:
            break
        else:
            taxid = new_taxid

    names_lineage = []
    for taxid in taxid_lineage:
        name = names_map[taxid]
        names_lineage.append(name)

    return ";".join(reversed(names_lineage))


def find_taxdump_file(kraken2_db, filename, extra_dirs=()):
    """Locate a taxdump file (names.dmp / nodes.dmp) for a Kraken2 database.

    ``kraken2-build`` leaves the taxonomy under ``<db>/taxonomy/``, but the *prebuilt*
    databases everyone downloads (genome-idx) ship ``names.dmp``/``nodes.dmp`` at the top
    level of the database directory instead - so hard-coding ``<db>/taxonomy/`` made the
    kraken2 module fail with a bare FileNotFoundError on any prebuilt database, after the
    (slow) classification had already run. Search both, then any caller-supplied fallback
    such as the configured TAXDUMP directory.
    """
    candidates = [
        os.path.join(kraken2_db, "taxonomy", filename),
        os.path.join(kraken2_db, filename),
    ]
    candidates += [os.path.join(d, filename) for d in extra_dirs if d]
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(
        "Could not find %s for the Kraken2 database %s. Looked in:\n  %s\n"
        "Prebuilt Kraken2 databases keep it at the top level; kraken2-build puts it under "
        "taxonomy/. If yours has neither, set TAXDUMP in metawrap2.toml to an NCBI taxdump "
        "directory (metawrap2 install-db taxdump)."
        % (filename, kraken2_db, "\n  ".join(candidates))
    )


def load_kraken_db_metadata(kraken2_db, taxdump_dirs=()):
    """Load NCBI name and parent-node maps needed to expand taxids into lineage strings."""
    print("Loading NCBI node names")
    names_path = find_taxdump_file(kraken2_db, "names.dmp", taxdump_dirs)
    names_map = {}
    with open(names_path) as handle:
        for line in handle:
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            name = cut[2]
            entry_type = cut[6]
            if entry_type == "scientific name":
                names_map[taxid] = name

    print("Loading NCBI taxonomic ranks")
    ranks_path = find_taxdump_file(kraken2_db, "nodes.dmp", taxdump_dirs)
    ranks_map = {}
    with open(ranks_path) as handle:
        for line in handle:
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            parent_taxid = cut[2]
            ranks_map[taxid] = parent_taxid
    return names_map, ranks_map


def translate_kraken2_annotations(
    annotation_file=None, kraken2_db=None, output=None, taxdump_dirs=()
):
    """Translate the kraken2 annotations in *annotation_file* from taxids to lineage strings."""
    print("Translating kraken2 annotations")
    if os.path.isfile(output):
        print("Looks like the translated kraken2 output already exists. Skipping...")
        return

    names_map, ranks_map = load_kraken_db_metadata(kraken2_db, taxdump_dirs)
    print("Writing translated taxonomy names to %s" % output)
    with open(output, "w") as out_handle, open(annotation_file) as handle:
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
    print(
        "Translating kraken2 annotations from %s, using metadata from the kraken2 database in %s; saving to %s"
        % (kraken_file, database_location, output_file)
    )
    translate_kraken2_annotations(
        annotation_file=kraken_file, kraken2_db=database_location, output=output_file
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
