"""Route read pairs to the bin they align to, at two stringencies, for reassembly.

Reads a SAM stream on stdin (``bwa mem`` against all bins concatenated) and writes, per bin,
four FASTQ files::

    <bin>.strict_1.fastq      <bin>.strict_2.fastq
    <bin>.permissive_1.fastq  <bin>.permissive_2.fastq

A pair is recruited to a bin when both mates name contigs belonging to that same bin and **at least
one of them actually aligned** - a pair with one mate hanging off the end of a contig is exactly the
read an assembler can use to extend it. It then goes into the *strict* set if the mates that aligned
carry fewer than ``strict`` mismatches between them (SAM ``NM`` tags summed), and into the
*permissive* set under the looser cutoff. ``reassemble_bins`` assembles both and keeps whichever
CheckM likes best, so the two stringencies are a way of trying both a conservative and a generous
read recruitment for every genome.

Whether a record aligned is read from the FLAG, and the mismatch sum counts only the mates that did.
Both used to come off the ``NM`` tag - its presence standing in for "aligned" and its absence
scoring as 0 - which made an unmapped mate contribute a *perfect* match and pushed half-aligned
pairs into the strict set, the one place they do not belong.

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
FLAG_UNMAPPED = 0x4  # this record did not align
FLAG_REVERSE = 0x10  # read aligned to the reverse strand
FLAG_FIRST = 0x40  # first mate in the pair
FLAG_SECOND = 0x80  # second mate in the pair
FLAG_SECONDARY = 0x100  # an additional alignment of a read that already has a primary one
FLAG_SUPPLEMENTARY = 0x800  # part of a chimeric alignment, SEQ hard-clipped to the piece

#: Records that are not the read itself, and must not be turned back into a FASTQ record.
#:
#: Both carry a SEQ that is not the read: a secondary alignment writes ``*`` (minimap2 does this for
#: every one it emits, and ``-ax map-ont`` emits them by default), and a supplementary alignment
#: writes only the hard-clipped piece it represents. Reconstituted as FASTQ, the first becomes a
#: one-character read and the second a truncated duplicate of one SPAdes is already being given -
#: so a reassembly was being fed junk reads, in proportion to how repetitive the genome is.
#:
#: They also break pairing: ``iter_pairs`` keys on the first/second-mate bits, which a supplementary
#: record carries too, so one would overwrite the real forward mate and be yielded against the next
#: reverse mate - a pair of two different reads, recruited on where neither of them aligned.
EXTRA_ALIGNMENT_FLAGS = FLAG_SECONDARY | FLAG_SUPPLEMENTARY

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


def mismatches(fields: List[str]) -> Optional[int]:
    """The SAM ``NM`` edit distance for one alignment record, or None if it carries no NM tag.

    None rather than 0, which is what this returned before. 0 is a *perfect* alignment, so a record
    whose edit distance is unknown was being scored as the best possible read and sorted into the
    strict set on that basis - an unknown promoted to the top of the ranking. None forces the caller
    to say what it wants to do about not knowing, and :func:`pair_mismatches` says: not strict.
    """
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field.rpartition(":")[2])
    return None


def aligned(fields: List[str]) -> bool:
    """Did this record actually align? Read from the FLAG, which is where SAM puts it.

    This used to test for the presence of an ``NM:i:`` tag, which is a claim about bwa rather than
    about SAM - bwa writes the edit distance on every aligned record, so its absence stood in for
    "unmapped". The proxy was load-bearing: an unmapped mate is not necessarily parked at RNAME
    ``*``, because the spec lets it be placed at its mate's reference (and bwa does exactly that), so
    the reference name alone cannot tell the two apart.

    Reading the flag instead matters now that ``reassemble_bins`` is not the last module hardcoded to
    one aligner. An aligner that omitted NM would have made every read look unmapped here, which
    fails by reassembling every bin from no reads rather than by stopping.
    """
    try:
        return not int(fields[1]) & FLAG_UNMAPPED
    except ValueError:
        return False


def pair_mismatches(forward: List[str], reverse: List[str]) -> Optional[int]:
    """Edit distance summed over the mates that aligned, or None if one of them does not say.

    Only aligned mates are counted, because an unmapped mate has no edit distance to contribute and
    reporting 0 for it would make a half-aligned pair look like a flawless one. None - an aligned
    mate with no NM tag - propagates, so the pair is never certified strict on a number nobody wrote.
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
    :data:`EXTRA_ALIGNMENT_FLAGS` - they are not the read, and they carry the mate bits that this
    pairs on). Records are paired by the FLAG's first/second-mate bits, tested with bitwise masks:
    the original indexed into ``bin(flag)``, which raises IndexError for any flag small enough to
    make the binary string short.
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
            # At least one mate must have actually aligned. Kept as "at least one" rather than
            # tightened to both: a pair with one mate hanging off the end of a contig is exactly the
            # read an assembler can use to extend it, which is the point of reassembling at all.
            if not (aligned(forward) or aligned(reverse)):
                continue
            counts["recruited"] += 1

            total = pair_mismatches(forward, reverse)
            r1, r2 = as_fastq(forward, 1), as_fastq(reverse, 2)
            if total is None:
                # Aligned, recruited, and its edit distance is unknown. Permissive takes it, strict
                # does not: strict is a claim about how well the read matches, and there is no
                # number here to support one.
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
