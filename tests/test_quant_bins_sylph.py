"""quant_bins --quantifier sylph: the profile command and the pivot of sylph's output.

sylph writes one line per detection, in no fixed order, and omits undetected genomes. These tests
check the pivot into the same ``bin_abundance_table.tsv`` the other quantifiers write, with every
bin present.
"""

from metawrap2.constants import TSV_SUFFIX
from metawrap2.modules import quant_bins

SYLPH_COLUMNS = (
    "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t"
    "Eff_cov\tANI_5-95_percentile\tEff_lambda"
)


def _pairs(*names):
    """``(sample, r1, r2)`` triples as ``collect_read_pairs`` returns them."""
    return [(n, "/reads/%s_1.fastq" % n, "/reads/%s_2.fastq" % n) for n in names]


def _profile(tmp_path, rows):
    """A sylph profile: rows are ``(r1 path, genome path, taxonomic, sequence)``."""
    raw = tmp_path / ("sylph_profile" + TSV_SUFFIX)
    lines = [SYLPH_COLUMNS]
    lines += ["%s\t%s\t%s\t%s\t99.5\t12.3\t99.1-99.9\t0.5" % row for row in rows]
    raw.write_text("\n".join(lines) + "\n")
    return str(raw)


def _table(tmp_path, raw, pairs, bins):
    out = tmp_path / "bin_abundance_table.tsv"
    with open(out, "w") as fh:
        rows = quant_bins._sylph_table(raw, pairs, bins, fh)
    return rows, [line.split("\t") for line in out.read_text().splitlines()]


def test_a_long_profile_becomes_a_bins_by_samples_table(tmp_path):
    raw = _profile(
        tmp_path,
        [
            ("/reads/s1_1.fastq", "/bins/bin_1.fasta", "70.0", "65.0"),
            ("/reads/s1_1.fastq", "/bins/bin_2.fasta", "30.0", "35.0"),
            ("/reads/s2_1.fastq", "/bins/bin_1.fasta", "10.0", "12.0"),
            ("/reads/s2_1.fastq", "/bins/bin_2.fasta", "90.0", "88.0"),
        ],
    )
    rows, table = _table(
        tmp_path, raw, _pairs("s1", "s2"), ["/bins/bin_1.fasta", "/bins/bin_2.fasta"]
    )
    assert rows == 2
    assert table == [
        ["Genomic bins", "s1", "s2"],
        ["bin_1", "65.0", "12.0"],
        ["bin_2", "35.0", "88.0"],
    ]
    # Sequence abundance, not taxonomic, to match the other quantifiers and the benchmark truth.
    assert table[1][1] != "70.0"

    # sylph omits genomes below its ANI threshold; they become zeros, so every bin has a row.
    raw = _profile(tmp_path, [("/reads/s2_1.fastq", "/bins/bin_2.fasta.gz", "100.0", "100.0")])
    rows, table = _table(
        tmp_path, raw, _pairs("s1", "s2"), ["/bins/bin_1.fasta", "/bins/bin_2.fasta.gz"]
    )
    assert rows == 2
    assert table[1:] == [["bin_1", "0.0", "0.0"], ["bin_2", "0.0", "100.0"]]

    # At 95, 55% of one sample's abundance went to bins absent from it; at 99, 25%.
    assert quant_bins.SYLPH_MIN_ANI == 99
    profile = quant_bins.SYLPH_PROFILE.format(
        db="b.syldb", sketches="a.sylsp", min_ani=quant_bins.SYLPH_MIN_ANI, threads=4, out="o.tsv"
    )
    assert "--minimum-ani 99" in profile
