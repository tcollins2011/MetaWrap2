import os

from metawrap2.modules import binning


def test_metabat2_is_given_a_seed(tmp_path, monkeypatch):
    """metaBAT2's ``--seed`` default is 0, which means "pick a random one".

    So the default binner - the baseline every other binner in the benchmark is ranked against - was
    the unreproducible one, and nothing about it looked like a neural net that would need a seed.
    Three benchmark arms ran it on a byte-identical depth table with a byte-identical command line
    and two produced 79 bins to the third's 78, which is the whole evidence that this mattered.
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
    # A fixed value, not a per-run one: two runs of the pipeline have to agree with each other.
    assert binning.METABAT2_SEED != 0


def test_a_zero_length_record_is_left_out_of_the_assembly_copy(tmp_path):
    """One empty record ends the run hours later, at the sort, and names neither file nor record.

    bwa indexes it as "@SQ SN:... LN:0" and samtools refuses the whole alignment. CAMI2's
    gold-standard assembly ships two, so this is the difference between a benchmark that runs and
    one that dies five variants in a row.
    """
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1\nACGT\n>_from__to_-42_total_0\n\n>contig_2\nTTTT\n")
    copy = tmp_path / "work" / "assembly.fa"
    os.makedirs(copy.parent)

    assert binning.copy_assembly(str(assembly), str(copy)) == (1, 0)
    text = copy.read_text()
    assert ">contig_1\nACGT\n" in text and ">contig_2\nTTTT\n" in text
    assert "_from__to_" not in text

    # An all-N contig goes too, however long: VAMB refuses to start on its zero tetranucleotide
    # vector, and CAMI2 GI carries 20 of them (28.9 Mbp of masked-out reference). A contig that merely
    # has gaps in it is normal and is kept.
    assembly.write_text(">gappy\nACGTNNNNACGT\n>masked\n" + "N" * 400 + "\n>contig_2\nTTTT\n")
    assert binning.copy_assembly(str(assembly), str(copy)) == (0, 1)
    text = copy.read_text()
    assert ">masked" not in text and ">gappy\nACGTNNNNACGT\n" in text
