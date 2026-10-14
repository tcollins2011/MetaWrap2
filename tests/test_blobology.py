from metawrap2 import aligners
from metawrap2.modules import blobology


def test_the_aligner_streams_into_the_bam_without_writing_a_sam():
    """The module used to write a full SAM per sample and never delete it.

    pipefail is part of the contract, not decoration: ``samtools view`` converts a truncated stream
    into a valid short BAM and exits 0, so without it an aligner that died mid-run would produce a
    blobplot drawn from partial coverage with nothing in the log to say so.
    """
    align = aligners.ALIGNERS["minimap2"]["paired"].format(
        threads=4, index="asm.fa.mmi", r1="a_1.fq", r2="a_2.fq"
    )
    cmd = blobology.ALIGN_BAM.format(align=align, threads=4, bam="a.bam")
    assert "{" not in cmd
    assert cmd.startswith("set -o pipefail; ")
    assert ".sam" not in cmd
    assert "minimap2 -ax sr" in cmd and "| samtools view" in cmd and cmd.endswith("> a.bam")
    # No sort: gc_cov_annotate sums aligned span per contig with a plain `samtools view`, so
    # record order is irrelevant and sorting would be work spent on nothing.
    assert "samtools sort" not in cmd
