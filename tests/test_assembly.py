import io
import subprocess
import sys

from metawrap2 import aligners
from metawrap2.modules import assembly
from metawrap2.scripts import rm_short_contigs, sam_to_fastq

# --- unmapped-read extraction -------------------------------------------------------------
#
# The hybrid assembly is built from whatever these select, so the selection and the plumbing that
# carries it are both worth pinning. It used to run as two SAMs on disk and a delete; it now runs
# as one pipe, and the records it picks have to be the same ones.

SAM = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:c1\tLN:100\n"
    "mapped\t0\tc1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n"
    "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\n"
)

#: Records no aligner-specific tag can classify. The selection is a property of SAM, and the module
#: that consumes it is moving off bwa - minimap2 is MetaWrap2's measured default everywhere else -
#: so each of these is a read the old NM-tag test got wrong or could get wrong.
AWKWARD_SAM = (
    "@SQ\tSN:c1\tLN:100\n"
    # Aligned, but written by a tool that does not emit NM. The old test called this unmapped, which
    # fails by handing MEGAHIT the whole read set rather than by stopping.
    "no_nm_tag\t0\tc1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\tAS:i:4\n"
    # Unmapped mate of a mapped read: flagged 0x4 but carrying its mate's coordinates, so it sits in
    # the middle of a coordinate-sorted file rather than at the end.
    "unmapped_mate\t133\tc1\t1\t0\t*\t=\t1\t0\tGGGG\tHHHH\n"
    # An unmapped read whose quality string happens to contain the tag the old test searched for.
    "qual_says_nm\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tNM:i:\n"
    # Supplementary alignment of an already-emitted read: a duplicate, not a second read.
    "chimeric\t2048\tc1\t50\t60\t2M2S\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n"
)


def test_the_unmapped_filter_reads_the_flag_and_not_an_aligners_tag(tmp_path):
    """Each record is one the NM-tag test classified on evidence that is not about alignment."""
    sam = tmp_path / "awkward.sam"
    sam.write_text(AWKWARD_SAM)
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out, unmapped_only=True)
    names = [ln[1:] for ln in out.getvalue().splitlines()[::4]]
    assert names == ["unmapped_mate", "qual_says_nm"]

    # Unfiltered, secondary and supplementary records are more alignments of a read, not more reads.
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out)
    names = [ln[1:] for ln in out.getvalue().splitlines()[::4]]
    assert names == ["no_nm_tag", "unmapped_mate", "qual_says_nm"]

    # ``-`` is what lets the aligner stream in, so the filter has to read a pipe as well as a file.
    proc = subprocess.run(
        [sys.executable, "-m", "metawrap2.scripts.sam_to_fastq", "--unmapped", "-"],
        input=SAM,
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout == "@unmapped\nCCCC\n+\nJJJJ\n"

    # And every aligner streams into that filter under pipefail, with no SAM written to disk.
    for name, spec in sorted(aligners.ALIGNERS.items()):
        cmd = assembly.EXTRACT_UNMAPPED.format(
            align=spec["paired"].format(threads=8, index="long.fasta", r1="a_1.fq", r2="a_2.fq"),
            fastq="unused.fastq",
        )
        assert "{" not in cmd, name
        assert cmd.startswith("set -o pipefail; "), name
        assert "--unmapped -" in cmd and cmd.endswith("> unused.fastq"), name
        assert cmd.count(">") == 1, name


def test_rm_short_contigs_filters_by_actual_sequence_length(tmp_path):
    # Length comes from the sequence, not from parsing the header, so this works on any
    # assembler's naming and does not depend on the file being sorted longest-first.
    fa = tmp_path / "s.fasta"
    fa.write_text(
        ">NODE_1_length_2000_cov_5\n" + "A" * 2000 + "\n"
        ">NODE_2_length_500_cov_3\n" + "C" * 500 + "\n"
        ">NODE_3_length_3000_cov_2\n" + "G" * 3000 + "\n"
    )
    out = io.StringIO()
    rm_short_contigs.remove_short_contigs(1000, str(fa), out)
    text = out.getvalue()
    assert "NODE_1_length_2000" in text
    assert "NODE_3_length_3000" in text  # no longer lost to the early break
    assert "NODE_2_length_500" not in text
