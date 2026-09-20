import os

from metawrap2.modules import bin_refinement


def test_command_templates_format_cleanly():
    cmd = bin_refinement.BINNING_REFINER.format(b1="binsA", b2="binsB", out="Refined_AB")
    assert "{" not in cmd and "-1 binsA -2 binsB -o Refined_AB" in cmd
    cmd3 = bin_refinement.BINNING_REFINER_3.format(b1="binsA", b2="binsB", b3="binsC", out="R")
    assert "{" not in cmd3 and "-3 binsC" in cmd3
    plot = bin_refinement.PLOT_BINNING.format(comp=70, cont=10, stats="a.stats b.stats")
    assert "{" not in plot and plot.endswith("70 10 a.stats b.stats")


def test_fix_contig_naming(tmp_path):
    f = tmp_path / "bin.fa"
    f.write_text(">contig=1 foo\nACGT\n>contig=2\nTTTT\n")
    bin_refinement._fix_contig_naming(str(f))
    text = f.read_text()
    # '=' becomes '_', and the description after the first whitespace is dropped: contig
    # identity must be the first token everywhere, or bins from different binners cannot be
    # matched up (metaBAT2 appends "total_depth=.. sample_depths=.." to every header).
    assert ">contig_1\n" in text
    assert ">contig_2\n" in text
    assert "foo" not in text
    assert "=" not in text


def test_fix_contig_naming_strips_metabat2_header_annotations(tmp_path):
    """Real metaBAT2 2.18 output; these annotations made its bins match nothing."""
    f = tmp_path / "bin.fa"
    f.write_text(
        ">NODE_2_length_158684_cov_2.764955 total_depth=41.98 "
        "sample_depths=18.3,23.5,0.2\nACGT\n"
    )
    bin_refinement._fix_contig_naming(str(f))
    assert f.read_text() == ">NODE_2_length_158684_cov_2.764955\nACGT\n"


def test_good_bins_thresholds(tmp_path):
    stats = tmp_path / "binsA.stats"
    stats.write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin.1\t90.0\t2.0\tx\ty\t1000\t100\n"  # good
        "bin.2\t60.0\t1.0\tx\ty\t1000\t100\n"  # too incomplete
        "bin.3\t95.0\t20.0\tx\ty\t1000\t100\n"  # too contaminated
    )
    good = bin_refinement._good_bins(str(stats), comp=70, cont=10)
    assert good == ["bin.1"]


def test_count_bins(tmp_path):
    d = tmp_path / "bins"
    d.mkdir()
    (d / "bin_001.fasta").write_text(">a\nAC\n")
    (d / "bin_002.fasta").write_text(">b\nAC\n")
    (d / "notes.txt").write_text("x")
    assert bin_refinement._count_bins(str(d)) == 2
    assert bin_refinement._count_bins(str(tmp_path / "missing")) == 0


def test_bin_set_dirs_excludes_m_and_o(tmp_path):
    for name in ("binsA", "binsB", "binsAB", "binsM", "binsO"):
        os.mkdir(tmp_path / name)
    (tmp_path / "binsA.stats").write_text("")
    assert bin_refinement._bin_set_dirs(str(tmp_path)) == ["binsA", "binsAB", "binsB"]
