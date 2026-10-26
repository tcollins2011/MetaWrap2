"""Route read pairs to the bin they align to, at two stringencies, for reassembly.

Reads a SAM stream on stdin (``bwa mem`` against all bins concatenated) and writes, per bin,
four FASTQ files::

    <bin>.strict_1.fastq      <bin>.strict_2.fastq
    <bin>.permissive_1.fastq  <bin>.permissive_2.fastq

A pair is recruited to a bin when both mates name contigs in that bin and at least one of them
aligned (a pair with one mate off the end of a contig can help the assembler extend it). It goes
into the strict set if the aligned mates carry fewer than ``strict`` mismatches in total (SAM
``NM`` tags summed), and into the permissive set under the looser cutoff. ``reassemble_bins``
assembles both and keeps whichever CheckM scores best.

Whether a record aligned is read from the FLAG, and the mismatch sum counts only aligned mates.
(Taking both from the ``NM`` tag scored an unmapped mate as a perfect match and put half-aligned
pairs in the strict set.)

Reads that aligned in reverse are reverse-complemented (and their quality strings reversed)
so the FASTQ is in the original read orientation, which is what an assembler expects.

Usage (``reassemble_bins`` does this for you)::

    bwa mem -t 8 assembly.fa r1.fastq r2.fastq \\
      | python -m metawrap2.scripts.filter_reads_for_bin_reassembly \\
            original_bins/ reads_for_reassembly/ 2 5

A rewrite of the original metaWRAP script with the same recruitment rules and output. Two
behaviours differ, both noted at their call sites: SAM flags are tested with bitwise masks
instead of by indexing into ``bin(flag)`` (which raised IndexError on small flag values), and
contigs are matched on their id instead of the whole header (so bins with tool-added header
annotations still match).
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

from ..constants import FASTA_EXTENSIONS
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar

# SAM FLAG bits we care about.
FLAG_UNMAPPED = 0x4  # this record did not align
FLAG_REVERSE = 0x10  # read aligned to the reverse strand
FLAG_FIRST = 0x40  # first mate in the pair
FLAG_SECOND = 0x80  # second mate in the pair
FLAG_SECONDARY = 0x100  # an additional alignment of a read that already has a primary one
FLAG_SUPPLEMENTARY = 0x800  # part of a chimeric alignment, SEQ hard-clipped to the piece

#: Records that are not the read itself, and must not be turned back into a FASTQ record.
#:
#: A secondary alignment's SEQ is ``*`` (minimap2 writes this for every one, and ``-ax map-ont``
#: emits them by default), and a supplementary alignment's SEQ is only its hard-clipped piece. As
#: FASTQ they would become a one-character read and a truncated duplicate.
#:
#: Supplementary records also carry the first/second-mate bits ``iter_pairs`` pairs on, so one
#: would replace the real forward mate and be paired with the wrong reverse mate.
EXTRA_ALIGNMENT_FLAGS = FLAG_SECONDARY | FLAG_SUPPLEMENTARY

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

STYLES = ("strict", "permissive")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def load_contig_bins(bin_folder: str) -> Dict[str, str]:
    """Map each contig id to the name of the bin that contains it.

    Keyed on :func:`contig_id` (the header up to the first whitespace), which is also what a
    SAM reference name is, so a bin whose headers carry tool-added annotations (metaBAT2
    appends ``total_depth=..``) still matches its alignments.
    """
    contig_bins: Dict[str, str] = {}
    for filename in sorted(os.listdir(bin_folder)):
        if not filename.endswith(FASTA_EXTENSIONS):
            continue
        bin_name = os.path.splitext(filename)[0]
        for header, _seq in iter_fasta(os.path.join(bin_folder, filename)):
            contig_bins[contig_id(header)] = bin_name
    return contig_bins


def mismatches(fields: List[str]) -> Optional[int]:
    """The SAM ``NM`` edit distance for one alignment record, or None if it carries no NM tag.

    None, not 0, because 0 means a perfect alignment. :func:`pair_mismatches` keeps such pairs out
    of the strict set.
    """
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field.rpartition(":")[2])
    return None


def aligned(fields: List[str]) -> bool:
    """True if this record aligned, from the FLAG's unmapped bit.

    The reference name is not enough: SAM allows an unmapped mate to be placed at its mate's
    reference, and bwa does so. Testing for an ``NM:i:`` tag works only for aligners that always
    write one, such as bwa.
    """
    try:
        return not int(fields[1]) & FLAG_UNMAPPED
    except ValueError:
        return False


def pair_mismatches(forward: List[str], reverse: List[str]) -> Optional[int]:
    """Edit distance summed over the mates that aligned, or None if one of them does not say.

    Only aligned mates are counted; an unmapped mate has no edit distance. If an aligned mate has
    no NM tag, the result is None and the pair cannot go into the strict set.
    """
    total = 0
    for fields in (forward, reverse):
        if not aligned(fields):
            continue
        nm = mismatches(fields)
        if nm is None:
            return None
        total += nm
    return total


def iter_pairs(stream: Iterator[str]) -> Iterator[Tuple[List[str], List[str]]]:
    """Yield (forward_fields, reverse_fields) for each mate pair in a SAM stream.

    Header lines are skipped, as are secondary and supplementary alignments (see
    :data:`EXTRA_ALIGNMENT_FLAGS`). Records are paired by the FLAG's first/second-mate bits, tested
    with bitwise masks (the original indexed into ``bin(flag)``, which raises IndexError for small
    flag values).
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
        if flag & EXTRA_ALIGNMENT_FLAGS:
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

    Files are opened on first use, so only bins that recruited reads get handles.
    ``reassemble_bins`` raises the open-file limit first, because many bins need thousands of
    handles at once.
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
        # Reads every alignment in the library from a pipe, which can take a long time with no
        # other output, so show a progress bar.
        for forward, reverse in bar(
            iter_pairs(sam_stream), desc="recruiting reads to bins", unit=" pair"
        ):
            counts["pairs"] += 1
            bin_name = bin_for_pair(forward, reverse, contig_bins)
            if bin_name is None:
                continue
            # At least one mate must have aligned. Not both: a pair with one mate off the end of a
            # contig can help the assembler extend it.
            if not (aligned(forward) or aligned(reverse)):
                continue
            counts["recruited"] += 1

            total = pair_mismatches(forward, reverse)
            r1, r2 = as_fastq(forward, 1), as_fastq(reverse, 2)
            if total is None:
                # Edit distance unknown: permissive set only.
                writers.write(bin_name, "permissive", 1, r1)
                writers.write(bin_name, "permissive", 2, r2)
                counts["permissive"] += 1
                continue
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
