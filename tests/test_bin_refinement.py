from metawrap2.modules import bin_refinement


def test_command_templates_format_cleanly():
    cmd = bin_refinement.BINNING_REFINER.format(b1="binsA", b2="binsB", out="Refined_AB")
    assert "{" not in cmd and "-1 binsA -2 binsB -o Refined_AB" in cmd
    cmd3 = bin_refinement.BINNING_REFINER_3.format(b1="binsA", b2="binsB", b3="binsC", out="R")
    assert "{" not in cmd3 and "-3 binsC" in cmd3
    plot = bin_refinement.PLOT_BINNING.format(comp=70, cont=10, stats="a.stats b.stats")
    assert "{" not in plot and plot.endswith("70 10 a.stats b.stats")


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
