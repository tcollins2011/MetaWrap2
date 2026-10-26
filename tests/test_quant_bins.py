import math

from metawrap2.scripts import split_salmon_out_into_bins


def _salmon_fixture(tmp_path, counts):
    """Two 200 kb single-contig bins, plus one ``.counts`` file per sample from *counts*.

    *counts* maps sample -> {contig: count}. Returns (quant_dir, bin_folder, assembly).
    """
    bins = tmp_path / "bins"
    bins.mkdir()
    for name, contig in (("bin_1.fasta", "g1_contig1"), ("bin_2.fasta", "g2_contig1")):
        (bins / name).write_text(">%s\n%s\n" % (contig, "A" * 200_000))
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(
        "".join(">%s\n%s\n" % (c, "A" * 200_000) for c in ("g1_contig1", "g2_contig1"))
    )
    quant = tmp_path / "quant_files"
    quant.mkdir()
    for sample, per_contig in counts.items():
        lines = ["transcript\tcount"] + ["%s\t%s" % (c, v) for c, v in per_contig.items()]
        (quant / ("%s.quant.counts" % sample)).write_text("\n".join(lines) + "\n")
    return str(quant), str(bins), str(assembly)


def _read_table(path):
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    header = rows[0]
    return {r[0]: dict(zip(header[1:], (float(v) for v in r[1:]))) for r in rows[1:]}


def test_a_bins_abundance_in_one_sample_does_not_leak_into_the_next(tmp_path):
    """Each cell holds that bin's abundance in that sample only.

    Guards against accumulators being reset through a bin name left over from the inner loop, which
    cleared only one bin and let the others accumulate across samples. Each bin's contig has all its
    coverage in one sample, so a leak appears as a nonzero cell where zero is expected.
    """
    quant, bins, assembly = _salmon_fixture(
        tmp_path,
        {
            "sampX": {"g1_contig1": 1_000_000.0, "g2_contig1": 0.0},
            "sampY": {"g1_contig1": 0.0, "g2_contig1": 1_000_000.0},
        },
    )
    out = tmp_path / "table.tsv"
    with open(out, "w") as fh:
        split_salmon_out_into_bins.build_table(quant, bins, assembly, fh)

    table = _read_table(out)
    assert table == {
        "bin_1": {"sampX": 1_000_000.0, "sampY": 0.0},
        "bin_2": {"sampX": 0.0, "sampY": 1_000_000.0},
    }

    # A bin missing from one sample's counts file gets NaN there and does not shift columns. Guards
    # against the header coming from the first bin's samples while rows used their own.
    second = tmp_path / "second"
    second.mkdir()
    quant, bins, assembly = _salmon_fixture(
        second,
        {"sampX": {"g1_contig1": 10.0, "g2_contig1": 20.0}, "sampY": {"g1_contig1": 30.0}},
    )
    out = second / "table.tsv"
    with open(out, "w") as fh:
        split_salmon_out_into_bins.build_table(quant, bins, assembly, fh)
    lines = out.read_text().splitlines()
    assert lines[0] == "Genomic bins\tsampX\tsampY"
    assert all(len(line.split("\t")) == 3 for line in lines)
    table = _read_table(out)
    assert table["bin_1"] == {"sampX": 10.0, "sampY": 30.0}
    assert table["bin_2"]["sampX"] == 20.0 and math.isnan(table["bin_2"]["sampY"])
