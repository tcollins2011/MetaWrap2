import os

from metawrap2.modules import binning


def test_metabat2_is_given_a_seed(tmp_path, monkeypatch):
    """metaBAT2 is run with a fixed non-zero ``--seed``; its default 0 means a random seed.

    Three benchmark arms with identical depth tables and command lines gave 79, 79 and 78 bins
    with the default.
    """
    calls = []
    monkeypatch.setattr(binning, "run", lambda cmd, **kw: calls.append(cmd))
    # Returns before the bin shuffling, which wants real metaBAT2 output on disk.
    monkeypatch.setattr(binning, "dry_run", lambda: True)

    class _Args:
        output = str(tmp_path / "out")
        threads = 4
        min_len = 1000

    binning._run_metabat2(_Args(), "e", str(tmp_path / "work"), "a.fa", "d.txt")
    assert "--seed %d" % binning.METABAT2_SEED in calls[0]
    # A fixed value, so two runs of the pipeline agree.
    assert binning.METABAT2_SEED != 0


def test_a_zero_length_record_is_left_out_of_the_assembly_copy(tmp_path):
    """Zero-length records are dropped from the assembly copy.

    bwa indexes them as "@SQ SN:... LN:0" and samtools then rejects the alignment, hours into the
    run. CAMI2's gold-standard assembly contains two.
    """
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1\nACGT\n>_from__to_-42_total_0\n\n>contig_2\nTTTT\n")
    copy = tmp_path / "work" / "assembly.fa"
    os.makedirs(copy.parent)

    assert binning.copy_assembly(str(assembly), str(copy)) == (1, 0)
    text = copy.read_text()
    assert ">contig_1\nACGT\n" in text and ">contig_2\nTTTT\n" in text
    assert "_from__to_" not in text

    # All-N contigs of any length are dropped too: VAMB fails on their zero tetranucleotide vector,
    # and CAMI2 GI has 20 of them (28.9 Mbp). Contigs that only contain gaps are kept.
    assembly.write_text(">gappy\nACGTNNNNACGT\n>masked\n" + "N" * 400 + "\n>contig_2\nTTTT\n")
    assert binning.copy_assembly(str(assembly), str(copy)) == (0, 1)
    text = copy.read_text()
    assert ">masked" not in text and ">gappy\nACGTNNNNACGT\n" in text
