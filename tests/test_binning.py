import os

from metawrap2.modules import binning


def test_command_templates_format_cleanly():
    # The surfaced command templates must fill without leftover placeholders.
    cmd = binning.METABAT2.format(
        assembly="a.fa", depth="d.txt", out="out", metabat_len=1500, threads=4, seed=1
    )
    assert "{" not in cmd and "metabat2 -i a.fa" in cmd and "-t 4" in cmd

    cmd = binning.CONCOCT.format(
        min_len=1000, threads=4, cov="cov.tsv", comp="comp.fa", concoct_out="out/", seed=1
    )
    assert "{" not in cmd and "concoct -l 1000" in cmd and "--seed 1" in cmd

    cmd = binning.SEMIBIN2.format(
        assembly="a.fa",
        bams="s1.bam s2.bam",
        out="out",
        threads=4,
        min_len=1000,
        seed=1,
        tmp="t",
        engine="cpu",
    )
    assert "{" not in cmd and "SemiBin2 single_easy_bin" in cmd
    cmd = binning.VAMB.format(
        out="out",
        assembly="a.fa",
        bamdir="work",
        min_len=1000,
        threads=4,
        minfasta=200000,
        seed=1,
        cuda="",
    )
    assert "{" not in cmd and "vamb bin default" in cmd
    cmd = binning.COMEBIN.format(
        assembly="a.fa", out="out", bamdir="work", threads=4, device="cpu", seed=1
    )
    assert "{" not in cmd and "run_comebin.sh" in cmd


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


def test_an_assembly_with_nothing_wrong_with_it_is_copied_whole(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1 len=4\nACGT\n>contig_2\nTTTT\n")
    copy = tmp_path / "assembly.fa"
    assert binning.copy_assembly(str(assembly), str(copy)) == (0, 0)
    assert copy.read_text() == ">contig_1 len=4\nACGT\n>contig_2\nTTTT\n"


def test_an_all_n_contig_is_left_out_however_long_it_is(tmp_path):
    """VAMB refuses to start on a contig with a zero tetranucleotide vector; the others take it.

    CAMI2's gold-standard gastrointestinal assembly carries 20 of these, 28.9 Mbp in total and the
    largest 6.6 Mbp - masked-out reference genomes, not assembly artefacts. No read aligns to them
    and they have no k-mer signal, so there is no assignment of one that could be correct. Dropping
    them is what makes the binners agree on their input.
    """
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1\nACGT\n>masked\n" + "N" * 400 + "\n>contig_2\nTTTT\n")
    copy = tmp_path / "assembly.fa"

    assert binning.copy_assembly(str(assembly), str(copy)) == (0, 1)
    text = copy.read_text()
    assert ">masked" not in text
    assert ">contig_1\nACGT\n" in text and ">contig_2\nTTTT\n" in text


def test_a_contig_with_some_n_in_it_is_kept(tmp_path):
    """The filter is "no usable base at all", not "contains an N" - gaps are normal in a contig."""
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">gappy\nACGTNNNNNNNNNNACGT\n")
    copy = tmp_path / "assembly.fa"

    assert binning.copy_assembly(str(assembly), str(copy)) == (0, 0)
    assert "gappy" in copy.read_text()


def test_split_concoct_bins(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1 foo\nACGT\n>contig_2\nTTTT\n>contig_3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\ncontig_1,0\ncontig_2,0\ncontig_3,1\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    # CONCOCT clusters 0 and 1 become bins 1 and 2: bin names count from 1 for every binner.
    assert os.path.exists(out / "bin_001.fasta")
    assert os.path.exists(out / "bin_002.fasta")
    first = (out / "bin_001.fasta").read_text()
    assert ">contig_1" in first and ">contig_2" in first
    assert ">contig_3" in (out / "bin_002.fasta").read_text()
