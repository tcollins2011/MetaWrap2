import math
import sys

from metawrap2.modules import quant_bins
from metawrap2.scripts import split_salmon_out_into_bins


def test_command_templates_format_cleanly():
    idx = quant_bins.SALMON_INDEX.format(threads=4, assembly="a.fa", index="idx")
    assert "{" not in idx and "salmon index -p 4" in idx and "-i idx" in idx
    q = quant_bins.SALMON_QUANT.format(index="idx", r1="a_1.fq", r2="a_2.fq", out="o", threads=4)
    assert "{" not in q and "--libType IU" in q and "--meta -p 4" in q
    # The heatmap helper is a metawrap2 module, so it must run in the *host* interpreter -
    # the quant_bins conda env holds salmon, not python+seaborn.
    heat = [a.format(table="t.tab", png="h.png") for a in quant_bins.MAKE_HEATMAP]
    assert heat == [sys.executable, "-m", "metawrap2.scripts.make_heatmap", "t.tab", "h.png"]


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
    """Each cell is that bin in that sample alone, with nothing carried over from earlier samples.

    The accumulators used to be reset through a bin name left over from the inner loop, so only one
    bin - whichever contig came last in the counts file - was ever cleared, and it was cleared
    partway through the pass that read it. Every other bin kept accumulating across samples and that
    one bin reported NaN throughout. This fixture makes the right answer unambiguous: each bin's
    contig has all the coverage in exactly one sample and none in the other, so any leak shows up as
    a nonzero cell where a zero belongs.
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


def test_a_bin_of_only_short_contigs_is_still_nan(tmp_path):
    """Kept from the original: under 1 kb a contig weighs nothing, so the weighted list is empty.

    Worth pinning separately from the leak fix, because NaN used to have two possible causes and
    only this one is intended - a reader who sees NaN should be able to conclude "short contigs".
    """
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin_1.fasta").write_text(">short1\n%s\n" % ("A" * 900))
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">short1\n%s\n" % ("A" * 900))
    quant = tmp_path / "quant_files"
    quant.mkdir()
    (quant / "sampX.quant.counts").write_text("transcript\tcount\nshort1\t50.0\n")

    out = tmp_path / "table.tsv"
    with open(out, "w") as fh:
        split_salmon_out_into_bins.build_table(str(quant), str(bins), str(assembly), fh)
    assert math.isnan(_read_table(out)["bin_1"]["sampX"])


def test_the_columns_follow_the_sample_list_not_each_bins_own_keys(tmp_path):
    """A bin absent from one sample's counts file must shift no other bin's columns.

    The header used to be written from the first bin's sample dict while every row was written from
    its own, so a bin seen in fewer samples than the first one silently emitted a short row and a
    bin seen in different samples emitted values under the wrong headings.
    """
    quant, bins, assembly = _salmon_fixture(
        tmp_path,
        {
            "sampX": {"g1_contig1": 10.0, "g2_contig1": 20.0},
            # bin_2's contig is missing from this sample's counts file entirely.
            "sampY": {"g1_contig1": 30.0},
        },
    )
    out = tmp_path / "table.tsv"
    with open(out, "w") as fh:
        split_salmon_out_into_bins.build_table(quant, bins, assembly, fh)

    lines = out.read_text().splitlines()
    assert lines[0] == "Genomic bins\tsampX\tsampY"
    assert all(len(line.split("\t")) == 3 for line in lines)
    table = _read_table(out)
    assert table["bin_1"] == {"sampX": 10.0, "sampY": 30.0}
    assert table["bin_2"]["sampX"] == 20.0
    assert math.isnan(table["bin_2"]["sampY"])
