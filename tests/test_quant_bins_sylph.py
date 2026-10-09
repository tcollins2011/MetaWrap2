"""quant_bins --quantifier sylph: the sketch-and-profile commands, and the pivot of sylph's output.

sylph is the one quantifier whose native output is not a table of the right shape. It reports one
line per detection, in whatever order it found them, and says nothing at all about a genome it did
not detect - so almost everything worth pinning here is about turning that into the same
``bin_abundance_table.tsv`` the other two quantifiers write, without inventing numbers and without
losing bins.
"""

import os

import pytest

from metawrap2.constants import TSV_SUFFIX
from metawrap2.modules import quant_bins

SYLPH_COLUMNS = (
    "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t"
    "Eff_cov\tANI_5-95_percentile\tEff_lambda"
)


def _pairs(*names):
    """``(sample, r1, r2)`` triples of the shape ``collect_read_pairs`` hands over."""
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


# ── the commands ──────────────────────────────────────────────────────────────────────────────


def test_the_sylph_templates_format_cleanly():
    sketch = quant_bins.SYLPH_SKETCH_GENOMES.format(
        genome_list="genomes.txt", db="/o/sylph_bins", threads=8
    )
    reads = quant_bins.SYLPH_SKETCH_READS.format(
        l1="r1.txt", l2="r2.txt", sketch_dir="/o/sk", threads=8
    )
    profile = quant_bins.SYLPH_PROFILE.format(
        db="/o/sylph_bins.syldb",
        sketches="/o/sk/a.sylsp /o/sk/b.sylsp",
        min_ani=quant_bins.SYLPH_MIN_ANI,
        threads=8,
        out="p.tsv",
    )
    for cmd in (sketch, reads, profile):
        assert "{" not in cmd
    # The genome sketch takes the same list file CoverM's --genome-fasta-list is given.
    assert "--gl genomes.txt" in sketch
    # Forward and reverse lists, matched by position, so every sample is sketched in one call.
    assert "--l1 r1.txt --l2 r2.txt" in reads
    # Raised above sylph's own 95 by measurement rather than taste - see SYLPH_MIN_ANI.
    assert "--minimum-ani 99" in profile


def test_the_read_lists_are_written_in_pair_order(tmp_path):
    l1, l2 = quant_bins._write_read_lists(_pairs("s1", "s2"), str(tmp_path))
    with open(l1) as fh:
        assert fh.read().split() == ["/reads/s1_1.fastq", "/reads/s2_1.fastq"]
    with open(l2) as fh:
        assert fh.read().split() == ["/reads/s1_2.fastq", "/reads/s2_2.fastq"]


# ── the pivot ─────────────────────────────────────────────────────────────────────────────────


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


def test_the_sequence_abundance_column_is_the_one_used_not_the_taxonomic_one(tmp_path):
    """The other quantifiers and the benchmark truth all measure sequence, not organisms.

    Taxonomic abundance divides out genome length, so for two genomes of different sizes the two
    columns disagree - and taking the wrong one would make sylph look wrong against a read-fraction
    truth for a reason that has nothing to do with sylph.
    """
    raw = _profile(tmp_path, [("/reads/s1_1.fastq", "/bins/bin_1.fasta", "70.0", "42.0")])
    _rows, table = _table(tmp_path, raw, _pairs("s1"), ["/bins/bin_1.fasta"])
    assert table[1] == ["bin_1", "42.0"]


def test_a_column_inserted_before_the_abundance_does_not_shift_the_numbers(tmp_path):
    """Columns are found by name: sylph's header is fixed names, unlike CoverM's read filenames."""
    raw = tmp_path / ("sylph_profile" + TSV_SUFFIX)
    raw.write_text(
        "Sample_file\tGenome_file\tSomething_new\tTaxonomic_abundance\tSequence_abundance\n"
        "/reads/s1_1.fastq\t/bins/bin_1.fasta\tnoise\t70.0\t42.0\n"
    )
    _rows, table = _table(tmp_path, str(raw), _pairs("s1"), ["/bins/bin_1.fasta"])
    assert table[1] == ["bin_1", "42.0"]


def test_an_undetected_bin_is_a_zero_rather_than_a_missing_row(tmp_path):
    """sylph omits a genome below its ANI threshold entirely, in every sample it is absent from.

    Pivoting only what the file contains would give a bin present in one sample a row with one cell
    and a bin present in none no row at all - a ragged table, and a bin silently dropped from a
    bins x samples table that its folder does contain.
    """
    raw = _profile(tmp_path, [("/reads/s2_1.fastq", "/bins/bin_2.fasta", "100.0", "100.0")])
    rows, table = _table(
        tmp_path,
        raw,
        _pairs("s1", "s2"),
        ["/bins/bin_1.fasta", "/bins/bin_2.fasta", "/bins/bin_3.fasta"],
    )
    assert rows == 3, "every bin sylph was given is a row, detected or not"
    assert table == [
        ["Genomic bins", "s1", "s2"],
        ["bin_1", "0.0", "0.0"],
        ["bin_2", "0.0", "100.0"],
        ["bin_3", "0.0", "0.0"],
    ]


def test_a_profile_with_no_detections_at_all_is_still_a_full_table(tmp_path):
    raw = _profile(tmp_path, [])
    rows, table = _table(tmp_path, raw, _pairs("s1"), ["/bins/bin_1.fasta"])
    assert rows == 1
    assert table == [["Genomic bins", "s1"], ["bin_1", "0.0"]]


def test_a_gzipped_bins_label_is_the_same_as_if_it_were_not(tmp_path):
    """sylph reads a gzipped genome; the row it produces must not be named ``bin_1.fasta``."""
    raw = _profile(tmp_path, [("/reads/s1_1.fastq", "/bins/bin_1.fasta.gz", "100.0", "100.0")])
    _rows, table = _table(tmp_path, raw, _pairs("s1"), ["/bins/bin_1.fasta.gz"])
    assert table[1][0] == "bin_1"


def test_the_columns_are_our_sample_labels_not_sylphs_read_filenames(tmp_path):
    raw = _profile(tmp_path, [("/reads/s1_1.fastq", "/bins/bin_1.fasta", "100.0", "100.0")])
    _rows, table = _table(tmp_path, raw, _pairs("s1"), ["/bins/bin_1.fasta"])
    assert table[0] == ["Genomic bins", "s1"]


def test_a_sample_sylph_names_by_a_path_we_do_not_recognise_is_an_error(tmp_path):
    """Guessing would put a real measurement in the wrong sample's column, which reads as a result."""
    raw = _profile(tmp_path, [("/elsewhere/other_1.fastq", "/bins/bin_1.fasta", "100.0", "100.0")])
    with pytest.raises(RuntimeError, match="none of the read files"):
        _table(tmp_path, raw, _pairs("s1"), ["/bins/bin_1.fasta"])


def test_a_profile_missing_the_abundance_column_says_so(tmp_path):
    raw = tmp_path / ("sylph_profile" + TSV_SUFFIX)
    raw.write_text(
        "Sample_file\tGenome_file\tAdjusted_ANI\n/reads/s1_1.fastq\t/b/bin_1.fasta\t99\n"
    )
    with pytest.raises(RuntimeError, match="Sequence_abundance"):
        _table(tmp_path, str(raw), _pairs("s1"), ["/b/bin_1.fasta"])


# ── wiring ────────────────────────────────────────────────────────────────────────────────────


def test_sylph_is_offered_by_the_flag_and_dispatched_to():
    assert "sylph" in quant_bins.QUANTIFIERS
    assert quant_bins.QUANTIFIERS["sylph"] is quant_bins._quantify_sylph
    args = quant_bins._parse_args(
        ["-b", "bins", "-o", "out", "--quantifier", "sylph", "r_1.fq", "r_2.fq"]
    )
    assert args.quantifier == "sylph"


def test_sylph_is_the_default_quantifier():
    """Measured on 76 CAMI2 GI bins: the best rank accuracy of the three, for 17% of salmon's CPU."""
    args = quant_bins._parse_args(["-b", "bins", "-o", "out", "r_1.fq", "r_2.fq"])
    assert args.quantifier == "sylph"


def test_sylph_ships_in_the_modules_own_env_because_it_is_the_default():
    """A default that needs a second install command is not a default, so it is not an opt-in env."""
    from metawrap2.config import OPTIONAL_ENV_TOOLS, OPTIONAL_ENVS

    assert "sylph" not in OPTIONAL_ENVS
    assert "sylph" not in OPTIONAL_ENV_TOOLS
    spec = os.path.join(os.path.dirname(quant_bins.__file__), "..", "envs", "bins.yaml")
    with open(spec) as fh:
        assert "sylph" in fh.read()


def test_min_ani_defaults_to_the_measured_value_and_is_adjustable():
    """95 sent 55% of one sample's abundance to bins whose genomes were absent from it; 99 did not."""
    assert quant_bins.SYLPH_MIN_ANI == 99
    assert quant_bins._parse_args(["-b", "b", "-o", "o", "r_1.fq", "r_2.fq"]).min_ani == 99
    args = quant_bins._parse_args(["-b", "b", "-o", "o", "--min-ani", "97.5", "r_1.fq", "r_2.fq"])
    assert args.min_ani == 97.5


def test_profile_template_takes_the_requested_ani():
    """The flag has to reach the command line, or it is a setting that silently does nothing."""
    command = quant_bins.SYLPH_PROFILE.format(
        db="b.syldb", sketches="a.sylsp", min_ani=97.5, threads=4, out="out.tsv"
    )
    assert "--minimum-ani 97.5" in command
