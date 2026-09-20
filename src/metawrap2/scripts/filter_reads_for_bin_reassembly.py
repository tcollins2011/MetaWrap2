"""Route read pairs to the bin they align to, at two stringencies, for reassembly.

Reads a SAM stream on stdin (``bwa mem`` against all bins concatenated) and writes, per bin,
four FASTQ files::

    <bin>.strict_1.fastq      <bin>.strict_2.fastq
    <bin>.permissive_1.fastq  <bin>.permissive_2.fastq

A pair is recruited to a bin when **both** mates align to contigs belonging to that same bin.
It then goes into the *strict* set if the two mates carry fewer than ``strict`` mismatches
between them (SAM ``NM`` tags summed), and into the *permissive* set under the looser
``permissive`` cutoff. ``reassemble_bins`` assembles both and keeps whichever CheckM likes
best, so the two stringencies are a way of trying both a conservative and a generous read
recruitment for every genome.

Reads that aligned in reverse are reverse-complemented (and their quality strings reversed)
so the FASTQ is in the original read orientation, which is what an assembler expects.

Usage (``reassemble_bins`` does this for you)::

    bwa mem -t 8 assembly.fa r1.fastq r2.fastq \\
      | python -m metawrap2.scripts.filter_reads_for_bin_reassembly \\
            original_bins/ reads_for_reassembly/ 2 5

This is a rewrite of the original script - same recruitment rules and same output -
restructured into functions from a tab-indented, module-level-``sys.argv`` script. Two
behavioural details were corrected in the process, both noted at their call sites: SAM flags
are now tested with bitwise masks rather than by indexing into ``bin(flag)`` (which raised
IndexError on small flag values), and contigs are matched on their id rather than their whole
header (so bins carrying tool-added header annotations are not silently unmatched).
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

from ..constants import FASTA_EXTENSIONS
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar

# SAM FLAG bits we care about.
FLAG_REVERSE = 0x10  # read aligned to the reverse strand
FLAG_FIRST = 0x40  # first mate in the pair
FLAG_SECOND = 0x80  # second mate in the pair

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

STYLES = ("strict", "permissive")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def load_contig_bins(bin_folder: str) -> Dict[str, str]:
    """Map each contig id to the name of the bin that contains it.

    Keyed on :func:`contig_id` (the header up to the first whitespace), which is also what a
    SAM reference name is - so a bin whose headers carry tool-added annotations (metaBAT2
    appends ``total_depth=..``) still matches its own alignments.
    """
    contig_bins: Dict[str, str] = {}
    for filename in sorted(os.listdir(bin_folder)):
        if not filename.endswith(FASTA_EXTENSIONS):
            continue
        bin_name = os.path.splitext(filename)[0]
        for header, _seq in iter_fasta(os.path.join(bin_folder, filename)):
            contig_bins[contig_id(header)] = bin_name
    return contig_bins


def mismatches(fields: List[str]) -> int:
    """The SAM ``NM`` edit distance for one alignment record, or 0 if it has no NM tag."""
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field.rpartition(":")[2])
    return 0


def has_nm_tag(fields: List[str]) -> bool:
    """True if this record carries an NM tag, i.e. it actually aligned."""
    return any(f.startswith("NM:i:") for f in fields[11:])


def iter_pairs(stream: Iterator[str]) -> Iterator[Tuple[List[str], List[str]]]:
    """Yield (forward_fields, reverse_fields) for each mate pair in a SAM stream.

    Header lines are skipped. Records are paired by the FLAG's first/second-mate bits, tested
    with bitwise masks: the original indexed into ``bin(flag)``, which raises IndexError for
    any flag small enough to make the binary string short.
    """
    forward: Optional[List[str]] = None
    for line in stream:
        if not line or line[0] == "@":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        try:
            flag = int(fields[1])
        except ValueError:
            continue
        if flag & FLAG_FIRST:
            forward = fields
        elif flag & FLAG_SECOND and forward is not None:
            yield forward, fields
            forward = None


def bin_for_pair(
    forward: List[str], reverse: List[str], contig_bins: Dict[str, str]
) -> Optional[str]:
    """The bin both mates belong to, or None if they disagree or are unbinned/unmapped."""
    f_ref, r_ref = forward[2], reverse[2]
    if f_ref == "*" and r_ref == "*":
        return None
    if f_ref == r_ref:
        return contig_bins.get(f_ref)
    f_bin, r_bin = contig_bins.get(f_ref), contig_bins.get(r_ref)
    if f_bin is not None and f_bin == r_bin:
        return f_bin
    return None


def as_fastq(fields: List[str], mate: int) -> str:
    """One FASTQ record for a SAM line, restored to the original read orientation."""
    seq, qual = fields[9], fields[10]
    if int(fields[1]) & FLAG_REVERSE:
        seq = reverse_complement(seq)
        qual = qual[::-1]
    return "@%s/%d\n%s\n+\n%s\n" % (fields[0], mate, seq, qual)


class BinWriters:
    """Lazily-opened output handles: four files per bin, closed together at the end.

    Files are opened on first use rather than up front, so a run only creates handles for bins
    that actually recruited reads. ``reassemble_bins`` raises the open-file limit before
    calling this, because a study with many bins needs thousands of handles at once.
    """

    def __init__(self, out_dir: str):
        self.out_dir = out_dir
        self._handles: Dict[Tuple[str, str, int], TextIO] = {}

    def write(self, bin_name: str, style: str, mate: int, record: str) -> None:
        key = (bin_name, style, mate)
        handle = self._handles.get(key)
        if handle is None:
            path = os.path.join(self.out_dir, "%s.%s_%d.fastq" % (bin_name, style, mate))
            handle = open(path, "w")  # noqa: SIM115 - closed together in close()
            self._handles[key] = handle
        handle.write(record)

    def close(self) -> None:
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    @property
    def bins(self) -> int:
        return len({key[0] for key in self._handles})


def filter_reads(
    sam_stream: Iterator[str], bin_folder: str, out_dir: str, strict: int, permissive: int
) -> Dict[str, int]:
    """Route pairs from *sam_stream* into per-bin FASTQ files. Returns simple counts."""
    print("loading contig to bin mappings...")
    contig_bins = load_contig_bins(bin_folder)
    print("%d contigs across %d bins" % (len(contig_bins), len(set(contig_bins.values()))))

    os.makedirs(out_dir, exist_ok=True)
    writers = BinWriters(out_dir)
    counts = {"pairs": 0, "recruited": 0, "strict": 0, "permissive": 0}

    print("Parsing the sam stream and routing reads to the bin they aligned to...")
    try:
        # This consumes every alignment of every read in the library off a pipe and can run
        # for a long time with no output at all; the bar makes it obvious it is still moving.
        for forward, reverse in bar(
            iter_pairs(sam_stream), desc="recruiting reads to bins", unit=" pair"
        ):
            counts["pairs"] += 1
            bin_name = bin_for_pair(forward, reverse, contig_bins)
            if bin_name is None:
                continue
            # At least one mate must have actually aligned (carry an NM tag).
            if not (has_nm_tag(forward) or has_nm_tag(reverse)):
                continue
            counts["recruited"] += 1

            total = mismatches(forward) + mismatches(reverse)
            r1, r2 = as_fastq(forward, 1), as_fastq(reverse, 2)
            if total < strict:
                writers.write(bin_name, "strict", 1, r1)
                writers.write(bin_name, "strict", 2, r2)
                counts["strict"] += 1
            if total < permissive:
                writers.write(bin_name, "permissive", 1, r1)
                writers.write(bin_name, "permissive", 2, r2)
                counts["permissive"] += 1
    finally:
        print("closing files")
        writers.close()

    print(
        "Finished splitting reads! %d pairs seen, %d recruited to a bin "
        "(%d strict, %d permissive)"
        % (counts["pairs"], counts["recruited"], counts["strict"], counts["permissive"])
    )
    return counts


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 4:
        sys.stderr.write(
            "usage: bwa mem ... | python -m metawrap2.scripts."
            "filter_reads_for_bin_reassembly <bin folder> <output dir> "
            "<strict mismatch cutoff> <permissive mismatch cutoff>\n"
        )
        return 2
    bin_folder, out_dir, strict, permissive = argv[0], argv[1], int(argv[2]), int(argv[3])
    filter_reads(sys.stdin, bin_folder, out_dir, strict, permissive)
    return 0


if __name__ == "__main__":
    sys.exit(main())
