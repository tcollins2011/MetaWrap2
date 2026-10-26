import io
import subprocess
import sys

from metawrap2 import aligners
from metawrap2.modules import assembly
from metawrap2.scripts import rm_short_contigs, sam_to_fastq

# --- unmapped-read extraction -------------------------------------------------------------
#
# The hybrid assembly is built from the reads selected here. Tests cover both the selection and
# the single pipe that carries it (no intermediate SAM files).

SAM = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:c1\tLN:100\n"
    "mapped\t0\tc1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n"
    "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\n"
)

#: Records that an aligner-specific tag (NM) cannot classify correctly. Selection must use the SAM
#: flag, since the module is not tied to bwa (minimap2 is the default elsewhere).
AWKWARD_SAM = (
    "@SQ\tSN:c1\tLN:100\n"
    # Aligned, written by a tool that does not emit NM. An NM-based test calls this unmapped.
    "no_nm_tag\t0\tc1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\tAS:i:4\n"
    # Unmapped mate of a mapped read: flagged 0x4 but with its mate's coordinates, so it sits in
    # the middle of a coordinate-sorted file.
    "unmapped_mate\t133\tc1\t1\t0\t*\t=\t1\t0\tGGGG\tHHHH\n"
    # An unmapped read whose quality string contains the text "NM:i:".
    "qual_says_nm\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tNM:i:\n"
    # Supplementary alignment of an already-emitted read: a duplicate, not a second read.
    "chimeric\t2048\tc1\t50\t60\t2M2S\t*\t0\t0\tAAAA\tIIII\tNM:i:0\n"
)


def test_the_unmapped_filter_reads_the_flag_and_not_an_aligners_tag(tmp_path):
    """Unmapped reads are selected by the SAM flag; secondary/supplementary records are dropped."""
    sam = tmp_path / "awkward.sam"
    sam.write_text(AWKWARD_SAM)
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out, unmapped_only=True)
    names = [ln[1:] for ln in out.getvalue().splitlines()[::4]]
    assert names == ["unmapped_mate", "qual_says_nm"]

    # Unfiltered output also drops secondary and supplementary records.
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out)
    names = [ln[1:] for ln in out.getvalue().splitlines()[::4]]
    assert names == ["no_nm_tag", "unmapped_mate", "qual_says_nm"]

    # ``-`` reads from stdin, so the aligner can stream in.
    proc = subprocess.run(
        [sys.executable, "-m", "metawrap2.scripts.sam_to_fastq", "--unmapped", "-"],
        input=SAM,
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.stdout == "@unmapped\nCCCC\n+\nJJJJ\n"

    # Every aligner streams into the filter under pipefail, with no SAM written to disk.
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
    # Length is measured from the sequence, not parsed from the header, and the input need not
    # be sorted longest-first.
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
    assert "NODE_3_length_3000" in text  # kept although it follows a short contig
    assert "NODE_2_length_500" not in text
