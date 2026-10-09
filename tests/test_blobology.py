import inspect

import pytest

from metawrap2 import aligners
from metawrap2.modules import blobology


def test_command_templates_format_cleanly():
    blastn = blobology.BLASTN.format(query="asm.fa", blastdb="/db/nt", threads=8)
    assert "{" not in blastn and "-num_threads 8" in blastn and "-db /db/nt" in blastn
    mm2 = aligners.ALIGNERS["minimap2"]["paired"].format(
        index="asm.fa.mmi", threads=4, r1="a_1.fq", r2="a_2.fq"
    )
    assert "{" not in mm2 and "-t 4" in mm2 and "a_1.fq a_2.fq" in mm2
    plot = blobology.MAKEBLOBPLOT_BASE.format(
        table="t", prop=0.005, taxlevel="bin", base="Unbinned"
    )
    assert "{" not in plot and plot.endswith("t 0.005 bin Unbinned")


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


@pytest.mark.parametrize("name", sorted(aligners.ALIGNERS))
def test_every_aligner_fits_blobologys_pipe(name):
    """blobology shares binning's registry, so a new entry has to work in both templates."""
    spec = aligners.ALIGNERS[name]
    index = aligners.index_path("asm.fa", name)
    cmd = blobology.ALIGN_BAM.format(
        align=spec["paired"].format(threads=4, index=index, r1="a_1.fq", r2="a_2.fq"),
        threads=4,
        bam="a.bam",
    )
    assert "{" not in cmd and cmd.endswith("> a.bam")
    assert spec["index_cmd"].format(assembly="asm.fa", index=index).startswith(name.split("-")[0])


def test_bowtie2_is_gone_from_the_module():
    """It was blobology's only aligner; dropping it is the point, so pin that it stays dropped."""
    source = inspect.getsource(blobology)
    assert "bowtie2" not in source.lower()
