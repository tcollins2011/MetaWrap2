"""quant_bins --quantifier coverm: one invocation for every sample, and the table it produces.

The salmon path and the coverm path write the same ``bin_abundance_table.tsv`` and everything
downstream reads it the same way, so what is worth pinning here is the translation between CoverM's
table and that contract: row labels go through the same :func:`bin_label` the salmon path uses,
columns are matched by position against our own sample labels, and CoverM's ``unmapped`` row - which
is not a bin - is gone.
"""

import os

import pytest

from metawrap2.constants import TSV_SUFFIX
from metawrap2.modules import quant_bins


def _coverm_output(tmp_path, header_samples, rows):
    raw = tmp_path / ("coverm_genome" + TSV_SUFFIX)
    lines = ["Genome\t" + "\t".join("%s Relative Abundance (%%)" % s for s in header_samples)]
    lines += ["\t".join(r) for r in rows]
    raw.write_text("\n".join(lines) + "\n")
    return str(raw)


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


def test_the_abundance_table_is_written_and_the_heatmap_drawn(coverm_run, tmp_path):
    commands, args, table = coverm_run
    with open(table) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "Genomic bins\ta\tb"
    assert lines[1:] == ["bin.1\t60.0\t20.0", "bin.2\t35.0\t75.0"]
    # The heatmap is the same picture from the same table, whichever quantifier filled it in.
    assert any(isinstance(c, list) and "metawrap2.scripts.make_heatmap" in c for c in commands)

    # Every sample on one command line, against a genome list naming the bins actually there.
    (coverm,) = [c for c in commands if isinstance(c, str) and c.startswith("coverm")]
    assert "--coupled a_1.fastq a_2.fastq b_1.fastq b_2.fastq" in coverm
    listing = os.path.join(args.output, "coverm_genomes.txt")
    assert "--genome-fasta-list %s" % listing in coverm
    with open(listing) as fh:
        assert [os.path.basename(p) for p in fh.read().split()] == ["bin.1.fasta", "bin.2.fasta"]
    # CoverM's own defaults for these are 0 - every alignment counts, however poorly it matches -
    # which on a MAG catalogue reports reads from a relative as this genome's abundance.
    assert "--min-covered-fraction 10" in coverm
    assert "--min-read-percent-identity 90" in coverm
    assert "--min-read-aligned-percent 75" in coverm

    # Columns are matched by position, so a miscount has to stop rather than mislabel.
    raw = _coverm_output(tmp_path, ["sampleA_1.fastq"], [["bin.1", "50.0"]])
    with pytest.raises(RuntimeError, match="1 value column"), open(tmp_path / "x.tsv", "w") as fh:
        quant_bins._coverm_table(raw, ["sampleA", "sampleB"], fh)
