import io

from metawrap2.modules import assembly
from metawrap2.scripts import (
    fix_megahit_contig_naming,
    rm_short_contigs,
    sam_to_fastq,
    sort_contigs,
)


def test_command_templates_format_cleanly():
    cmd = assembly.METASPADES.format(tmp="t", threads=8, mem=24, out="o", r1="a", r2="b")
    assert "{" not in cmd and "metaspades.py" in cmd and "-t 8" in cmd

    cmd = assembly.MEGAHIT_PAIRED.format(
        r1="a", r2="b", out="o", tmp="t", threads=4, mem=24_000_000_000
    )
    assert "{" not in cmd and "--continue" in cmd

    cmd = assembly.MEGAHIT_LEFTOVERS.format(
        reads="u.fastq", out="o", threads=4, mem=24_000_000_000, tmp="t"
    )
    assert "{" not in cmd and "-r u.fastq" in cmd

    cmd = assembly.QUAST.format(threads=2, out="q", quast_min=500, assembly="f.fasta")
    assert "{" not in cmd and "-m 500" in cmd


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














