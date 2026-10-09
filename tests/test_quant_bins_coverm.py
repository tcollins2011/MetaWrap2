"""quant_bins --quantifier coverm: one invocation for every sample, and the table it produces.

The salmon path and the coverm path write the same ``bin_abundance_table.tsv`` and everything
downstream reads it the same way, so what is worth pinning here is the translation between CoverM's
table and that contract: row labels go through the same :func:`bin_label` the salmon path uses,
columns are matched by position against our own sample labels, and CoverM's ``unmapped`` row - which
is not a bin - is gone.
"""

import inspect
import os

import pytest

from metawrap2.constants import TSV_SUFFIX, bin_label
from metawrap2.modules import quant_bins
from metawrap2.scripts import split_salmon_out_into_bins


def test_the_coverm_template_formats_cleanly():
    cmd = quant_bins.COVERM_GENOME.format(
        genome_list="genomes.txt",
        coupled="a_1.fq a_2.fq b_1.fq b_2.fq",
        methods=quant_bins.COVERM_METHOD,
        min_covered=quant_bins.COVERM_MIN_COVERED_FRACTION,
        mapper=quant_bins.COVERM_MAPPER,
        threads=8,
        out="t.tsv",
    )
    assert "{" not in cmd
    # The bins are named in a file, so there is no single extension to guess at.
    assert "coverm genome --genome-fasta-list genomes.txt" in cmd
    assert "--genome-fasta-extension" not in cmd
    # Every sample on one command line: that is why there is no per-sample loop.
    assert "--coupled a_1.fq a_2.fq b_1.fq b_2.fq" in cmd
    assert "--methods relative_abundance" in cmd
    # A few reads on a conserved gene are not evidence that a genome is present.
    assert "--min-covered-fraction 10" in cmd


def test_the_genome_list_names_every_bin_whatever_its_extension(tmp_path):
    """A folder mixing plain and gzipped bins is one CoverM run, not a choice between them.

    ``-x`` takes a single extension, so a directory holding both could only ever have named half
    its bins - and the other half would go missing with no error, because CoverM's answer to an
    extension that matches nothing is an empty genome set.
    """
    for name in ("bin_1.fasta", "bin_2.fasta.gz", "bin_3.fa", "bin_4.fna"):
        (tmp_path / name).write_text(">c1\nACGT\n")
    listing = tmp_path / "genomes.txt"
    assert len(quant_bins._write_genome_list(str(tmp_path), str(listing))) == 4
    assert [os.path.basename(p) for p in listing.read_text().split()] == [
        "bin_1.fasta",
        "bin_2.fasta.gz",
        "bin_3.fa",
        "bin_4.fna",
    ]
    # Absolute paths, because CoverM is not run from the bin folder.
    assert all(os.path.isabs(p) for p in listing.read_text().split())


def test_the_genome_list_leaves_out_what_is_not_a_bin(tmp_path):
    """A bin folder that has been through bin_refinement also holds tables and tool output."""
    (tmp_path / "bin_1.fasta").write_text(">c1\nACGT\n")
    (tmp_path / "bin_1.stats.tsv").write_text("bin\tcompleteness\n")
    (tmp_path / "checkm").mkdir()
    listing = tmp_path / "genomes.txt"
    assert len(quant_bins._write_genome_list(str(tmp_path), str(listing))) == 1
    assert [os.path.basename(p) for p in listing.read_text().split()] == ["bin_1.fasta"]


def test_an_empty_bin_folder_yields_no_genomes(tmp_path):
    listing = tmp_path / "genomes.txt"
    assert quant_bins._write_genome_list(str(tmp_path), str(listing)) == []
    assert listing.read_text() == ""


def _coverm_output(tmp_path, header_samples, rows):
    raw = tmp_path / ("coverm_genome" + TSV_SUFFIX)
    lines = ["Genome\t" + "\t".join("%s Relative Abundance (%%)" % s for s in header_samples)]
    lines += ["\t".join(r) for r in rows]
    raw.write_text("\n".join(lines) + "\n")
    return str(raw)


def test_the_table_gets_our_sample_labels_and_salmons_row_labels(tmp_path):
    # CoverM titles its columns after the read files, and emits them in --coupled order.
    raw = _coverm_output(
        tmp_path,
        ["sampleA_1.fastq", "sampleB_1.fastq"],
        [["unmapped", "20.0", "30.0"], ["bin.1", "50.0", "10.0"], ["bin.2", "30.0", "60.0"]],
    )
    out = tmp_path / "bin_abundance_table.tsv"
    with open(out, "w") as fh:
        rows = quant_bins._coverm_table(raw, ["sampleA", "sampleB"], fh)

    assert rows == 2  # unmapped is not a bin and is not in a bins x samples table
    assert out.read_text().splitlines() == [
        "Genomic bins\tsampleA\tsampleB",
        "bin.1\t50.0\t10.0",
        "bin.2\t30.0\t60.0",
    ]


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("bin_001.fasta", "bin_001"),
        ("bin.1.fa", "bin.1"),
        ("bin_001.fna.gz", "bin_001"),
        # A gzipped bin is the case that used to disagree: stripping the last dotted component
        # leaves "bin_001.fasta", which is not what CoverM calls the same bin.
        ("bin_001.fasta.gz", "bin_001"),
        # Not a FASTA extension, so not something to guess at: E.coli must not become E.
        ("E.coli.txt", "E.coli.txt"),
        ("/a/path/bin_002.fasta", "bin_002"),
    ],
)
def test_bin_label_is_the_one_name_a_bin_goes_by(filename, expected):
    assert bin_label(filename) == expected


def test_both_quantifiers_label_a_row_with_the_same_function(tmp_path):
    """Both tables key a bin the same way, or a user cannot line the two up at all.

    Asserted on the source of the two writers rather than only on their output, because agreement
    that holds because both happen to call the same helper is a guarantee, while agreement observed
    on one set of bin names is a coincidence waiting for a gzipped bin.
    """
    assert "bin_label(" in inspect.getsource(split_salmon_out_into_bins.build_table)
    assert "bin_label(" in inspect.getsource(quant_bins._coverm_table)
    # And end to end, on the name CoverM would emit for a gzipped bin.
    raw = _coverm_output(tmp_path, ["s_1.fq"], [["bin_001", "1.0"]])
    with open(tmp_path / "t.tsv", "w") as fh:
        quant_bins._coverm_table(raw, ["s"], fh)
    assert (tmp_path / "t.tsv").read_text().splitlines()[1].split("\t")[0] == bin_label(
        "bin_001.fasta.gz"
    )


def test_rows_come_out_sorted_whatever_order_coverm_used(tmp_path):
    """CoverM's own output is not sorted - it emitted bin_2 before bin_1 on a real two-bin run."""
    raw = _coverm_output(
        tmp_path,
        ["s_1.fq"],
        [["bin_2", "40.0"], ["unmapped", "0"], ["bin_10", "10.0"], ["bin_1", "50.0"]],
    )
    with open(tmp_path / "t.tsv", "w") as fh:
        assert quant_bins._coverm_table(raw, ["s"], fh) == 3
    labels = [line.split("\t")[0] for line in (tmp_path / "t.tsv").read_text().splitlines()[1:]]
    assert labels == sorted(labels)


def test_a_column_count_that_does_not_match_the_samples_is_an_error_not_a_shifted_table(tmp_path):
    """Columns are matched by position, so a miscount has to stop rather than mislabel."""
    raw = _coverm_output(tmp_path, ["sampleA_1.fastq"], [["bin.1", "50.0"]])
    with pytest.raises(RuntimeError, match="1 value column"), open(tmp_path / "out.tsv", "w") as fh:
        quant_bins._coverm_table(raw, ["sampleA", "sampleB"], fh)


class _Checkpoint:
    """Every step still to do, recording what was marked done."""

    def __init__(self):
        self.finished = []

    def todo(self, name):
        return True

    def done(self, name):
        self.finished.append(name)


@pytest.fixture
def coverm_run(tmp_path, monkeypatch):
    """Run _quantify_coverm against a fake CoverM, returning (commands, args, table)."""
    bins = tmp_path / "bins"
    bins.mkdir()
    for name in ("bin.1.fasta", "bin.2.fasta"):
        (bins / name).write_text(">c1\nACGT\n")
    output = tmp_path / "out"
    output.mkdir()

    commands = []

    def fake_run(cmd, **kw):
        commands.append(cmd)
        if isinstance(cmd, str) and cmd.startswith("coverm"):
            dest = cmd.split("--output-file ")[1].split()[0]
            with open(dest, "w") as fh:
                fh.write("Genome\ta_1.fq Relative Abundance (%)\tb_1.fq Relative Abundance (%)\n")
                fh.write("unmapped\t5.0\t5.0\nbin.1\t60.0\t20.0\nbin.2\t35.0\t75.0\n")
        return 0

    monkeypatch.setattr(quant_bins, "run", fake_run)
    monkeypatch.setattr(quant_bins, "dry_run", lambda: False)

    args = quant_bins._parse_args(
        [
            "-b",
            str(bins),
            "-o",
            str(output),
            "--quantifier",
            "coverm",
            "-t",
            "8",
            "a_1.fastq",
            "a_2.fastq",
            "b_1.fastq",
            "b_2.fastq",
        ]
    )
    pairs = [("a", "a_1.fastq", "a_2.fastq"), ("b", "b_1.fastq", "b_2.fastq")]
    table = os.path.join(str(output), "bin_abundance_table" + TSV_SUFFIX)
    ckpt = _Checkpoint()
    quant_bins._quantify_coverm(args, pairs, table, None, ckpt)
    return commands, args, table


def test_every_sample_goes_into_a_single_coverm_invocation(coverm_run):
    commands, args, _table = coverm_run
    coverm = [c for c in commands if isinstance(c, str) and c.startswith("coverm")]
    assert len(coverm) == 1
    assert "--coupled a_1.fastq a_2.fastq b_1.fastq b_2.fastq" in coverm[0]
    # No pseudo-transcriptome is involved, so nothing is indexed.
    assert not any(isinstance(c, str) and "salmon" in c for c in commands)
    # The genome list is written next to the output and names the bins that were actually there.
    listing = os.path.join(args.output, "coverm_genomes.txt")
    assert "--genome-fasta-list %s" % listing in coverm[0]
    with open(listing) as fh:
        assert [os.path.basename(p) for p in fh.read().split()] == ["bin.1.fasta", "bin.2.fasta"]


def test_neither_quantifier_trips_over_a_subdirectory_in_the_bin_folder(tmp_path):
    """bin_refinement leaves tool subdirectories among the bins; reading one used to crash.

    ``os.listdir`` handed ``checkm/`` to the FASTA reader, which raised a bare IsADirectoryError
    naming only the path - so the salmon path died on a bin folder that a user had every reason to
    consider valid. Both the concatenated assembly and the contig->bin map now skip non-FASTA
    entries, and the CoverM genome list never names them in the first place.
    """
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin_1.fasta").write_text(">c1\nACGTACGT\n")
    (bins / "checkm").mkdir()
    (bins / "bin_1.stats.tsv").write_text("bin\tcompleteness\nbin_1\t95\n")

    assembly = tmp_path / "assembly.fa"
    quant_bins._concat_bins_into_assembly(str(bins), str(assembly))
    assert assembly.read_text() == ">c1\nACGTACGT\n"
    assert split_salmon_out_into_bins.load_bins(str(bins)) == {"c1": "bin_1.fasta"}


def test_the_abundance_table_is_written_and_the_heatmap_drawn(coverm_run):
    commands, _args, table = coverm_run
    with open(table) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "Genomic bins\ta\tb"
    assert lines[1:] == ["bin.1\t60.0\t20.0", "bin.2\t35.0\t75.0"]
    # The heatmap is the same picture from the same table, whichever quantifier filled it in.
    assert any(isinstance(c, list) and "metawrap2.scripts.make_heatmap" in c for c in commands)
