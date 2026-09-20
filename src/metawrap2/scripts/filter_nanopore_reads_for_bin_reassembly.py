"""Route long (nanopore) reads to the bin they align to, for reassembly.

The long-read counterpart of :mod:`metawrap2.scripts.filter_reads_for_bin_reassembly`. Reads a
SAM stream on stdin (``minimap2 -ax map-ont`` against all bins concatenated) and writes one
``<bin>.nanopore.fastq`` per bin containing the reads that aligned to that bin's contigs.

There is no strict/permissive split here: long reads carry far more mismatches than a
short-read SNP cutoff would tolerate, so any aligned read is recruited and the mismatch count
is not used. ``reassemble_bins`` passes these to SPAdes as ``--nanopore`` alongside the
short-read pairs.

Usage (``reassemble_bins`` does this for you)::

    minimap2 -ax map-ont assembly.fa nanopore.fastq \\
      | python -m metawrap2.scripts.filter_nanopore_reads_for_bin_reassembly \\
            original_bins/ reads_for_reassembly/

Rewritten from the original tab-indented script. As with the short-read version, SAM flags are
now tested with a bitwise mask instead of by indexing into ``bin(flag)`` (which the original
had to wrap in ``try/except IndexError`` precisely because it broke on small flags), and
contigs are matched on their id rather than their whole header.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterator, List, Optional, TextIO

from ..progress import bar
from .filter_reads_for_bin_reassembly import (
    FLAG_REVERSE,
    has_nm_tag,
    load_contig_bins,
    reverse_complement,
)

__all__ = ["filter_reads", "main"]


def _records(stream: Iterator[str], contig_bins: Dict[str, str]) -> Iterator[List[str]]:
    """Yield the SAM records that aligned to a contig belonging to a known bin."""
    for line in stream:
        if not line or line[0] == "@":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11 or fields[2] == "*":
            continue
        if fields[2] not in contig_bins:
            continue
        if not has_nm_tag(fields):
            continue
        yield fields


def filter_reads(sam_stream: Iterator[str], bin_folder: str, out_dir: str) -> Dict[str, int]:
    """Route long reads from *sam_stream* into per-bin FASTQ files. Returns simple counts."""
    print("loading contig to bin mappings...")
    contig_bins = load_contig_bins(bin_folder)
    print("%d contigs across %d bins" % (len(contig_bins), len(set(contig_bins.values()))))

    os.makedirs(out_dir, exist_ok=True)
    handles: Dict[str, TextIO] = {}
    counts = {"reads": 0, "recruited": 0}

    print("Parsing the sam stream and routing reads to the bin they aligned to...")
    try:
        for fields in bar(
            _records(sam_stream, contig_bins), desc="recruiting nanopore reads", unit=" read"
        ):
            counts["reads"] += 1
            bin_name = contig_bins[fields[2]]
            handle = handles.get(bin_name)
            if handle is None:
                path = os.path.join(out_dir, "%s.nanopore.fastq" % bin_name)
                handle = open(path, "w")  # noqa: SIM115 - held open for the whole stream
                handles[bin_name] = handle

            seq, qual = fields[9], fields[10]
            if int(fields[1]) & FLAG_REVERSE:
                seq = reverse_complement(seq)
                qual = qual[::-1]
            handle.write("@%s/1\n%s\n+\n%s\n" % (fields[0], seq, qual))
            counts["recruited"] += 1
    finally:
        print("closing files")
        for handle in handles.values():
            handle.close()

    print(
        "Finished splitting reads! %d reads recruited across %d bins"
        % (counts["recruited"], len(handles))
    )
    return counts


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        sys.stderr.write(
            "usage: minimap2 -ax map-ont ... | python -m metawrap2.scripts."
            "filter_nanopore_reads_for_bin_reassembly <bin folder> <output dir>\n"
        )
        return 2
    filter_reads(sys.stdin, argv[0], argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
