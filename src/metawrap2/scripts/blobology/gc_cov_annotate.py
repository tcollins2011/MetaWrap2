#!/usr/bin/env python
# Build the blobplot data table: per-contig length, GC, coverage (from BAM files) and
# taxonomy (walked up the NCBI taxdump from each contig's blast taxid).
#
# Original tool by Sujai Kumar (blaxterlab/blobology, 2013); this is a drop-in
# replacement with the same inputs and output columns.
#
# Usage:
#   gc_cov_annotate.py --blasttaxid FILE --assembly FASTA --out FILE \
#       --taxdump DIR --bam BAM [BAM ...] --taxlist species genus family ...
import argparse
import os
import re
import subprocess
import sys

from metawrap2.io.seqio import iter_fasta  # noqa: E402

_NODE_RE = re.compile(r"^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|")
_NAME_RE = re.compile(r"^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name")
_CIGAR_RE = re.compile(r"(\d+)[MIDNP]")


def load_nodes_names(taxdump_dir):
    parent, rank, name = {}, {}, {}
    with open(os.path.join(taxdump_dir, "nodes.dmp")) as fh:
        for line in fh:
            m = _NODE_RE.match(line)
            if m:
                parent[m.group(1)] = m.group(2)
                rank[m.group(1)] = m.group(3)
    with open(os.path.join(taxdump_dir, "names.dmp")) as fh:
        for line in fh:
            m = _NAME_RE.match(line)
            if m:
                name[m.group(1)] = m.group(2)
    return parent, rank, name


def lineage_taxids(taxid, parent):
    """Walk from a taxid up to the root, returning all taxids in the lineage."""
    chain = [taxid]
    current = taxid
    seen = {taxid}
    while current in parent and parent[current] != current and parent[current] not in seen:
        current = parent[current]
        chain.insert(0, current)
        seen.add(current)
    return chain


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--blasttaxid", required=True)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--taxdump", default=".")
    ap.add_argument("--bam", nargs="*", default=[])
    ap.add_argument("--cov", nargs="*", default=[])
    ap.add_argument("--taxlist", nargs="*",
                    default=["species", "order", "phylum", "superkingdom"])
    args = ap.parse_args()

    out_file = args.out or (args.assembly + ".txt")
    wanted_ranks = set(args.taxlist)

    sys.stderr.write("Loading taxonomy from %s ...\n" % args.taxdump)
    parent, rank, name = load_nodes_names(args.taxdump)

    # taxonomy per contig, keyed by rank -> name
    contig_taxinfo = {}
    with open(args.blasttaxid) as fh:
        for line in fh:
            m = re.match(r"^(\S+)\t(\d+)", line)
            if not m:
                continue
            seqid, taxid = m.group(1), m.group(2)
            info = {}
            for tid in lineage_taxids(taxid, parent):
                r = rank.get(tid)
                if r in wanted_ranks and tid in name:
                    info[r] = name[tid]
            contig_taxinfo[seqid] = info

    # length / GC per contig (compression-transparent)
    sys.stderr.write("Loading assembly %s ...\n" % args.assembly)
    length, gccount, nonatgc, cov = {}, {}, {}, {}
    order = []
    for seqid, seq in iter_fasta(args.assembly):
        order.append(seqid)
        length[seqid] = len(seq)
        gccount[seqid] = sum(1 for b in seq if b in "gcGC")
        nonatgc[seqid] = sum(1 for b in seq if b not in "atgcATGC")
        cov[seqid] = {}

    # coverage per BAM (total aligned reference span / contig length), via samtools
    for bam in args.bam:
        sys.stderr.write("Reading %s ...\n" % bam)
        proc = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
        for sam in proc.stdout:
            if not sam or sam[0] in "@#" or sam.strip() == "":
                continue
            f = sam.rstrip("\n").split("\t")
            if len(f) < 6 or (int(f[1]) & 4) == 4:
                continue
            ref = f[2]
            if ref not in length:
                sys.stderr.write("ContigID %s in %s but not in assembly\n" % (ref, bam))
                continue
            span = sum(int(n) for n in _CIGAR_RE.findall(f[5]))
            cov[ref][bam] = cov[ref].get(bam, 0.0) + span / length[ref]
        if proc.wait() != 0:
            sys.exit("samtools view failed on %s" % bam)

    # two-column coverage files (seqid, mean depth)
    for cov_file in args.cov:
        with open(cov_file) as fh:
            for line in fh:
                m = re.match(r"^(\S+)\s+(\S+)", line)
                if m and m.group(1) in cov:
                    cov[m.group(1)][cov_file] = float(m.group(2))

    sys.stderr.write("Writing %s ...\n" % out_file)
    with open(out_file, "w") as out:
        cols = ["cov_" + b for b in args.bam] + ["cov_" + c for c in args.cov]
        out.write("seqid\tlen\tgc" + "".join("\t" + c for c in cols))
        out.write("".join("\ttaxlevel_" + t for t in args.taxlist) + "\n")
        for seqid in order:
            denom = length[seqid] - nonatgc[seqid]
            gc = gccount[seqid] / denom if denom else 0
            row = [seqid, str(length[seqid]), str(gc)]
            for b in args.bam:
                row.append(str(cov[seqid].get(b, 0)))
            for c in args.cov:
                row.append(str(cov[seqid].get(c, 0)))
            for t in args.taxlist:
                row.append(contig_taxinfo.get(seqid, {}).get(t, "Not annotated"))
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
