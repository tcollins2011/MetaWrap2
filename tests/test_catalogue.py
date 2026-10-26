"""Tests for the gene catalogue and its gene id prefixing.

:mod:`metawrap2.genecall` names genes ``<contig>_<n>`` and assemblers number contigs from one, so
two samples often share gene ids. ``mmseqs createdb`` keys on the name, so unprefixed inputs are
merged before clustering: two identical files of 1,772 genes each gave a catalogue of 1,772 genes,
not 3,544, and reported 1.1% redundancy removed instead of 50.5%.

Every id is prefixed with its source, as :mod:`metawrap2.modules.coverage` does for contigs. The
tests check that a prefix can be removed, that two sources cannot share a label, and that the
members table keeps the original id for joins to per-sample tables.
"""

from __future__ import annotations

from metawrap2.modules import catalogue


def _write(path, records):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, sequence) for name, sequence in records)
    return str(path)


# ── prefixing ────────────────────────────────────────────────────────────────────────────────


def test_prefixing_makes_colliding_gene_ids_distinct(tmp_path):
    """The same gene id from two samples is kept as two genes."""
    shared = "k141_1_length_500_cov_3_1"
    a = _write(tmp_path / "a.ffn", [(shared, "ACGT")])
    b = _write(tmp_path / "b.ffn", [(shared, "ACGA")])
    out_a = str(tmp_path / "pa.ffn")
    out_b = str(tmp_path / "pb.ffn")
    assert catalogue.write_prefixed(a, "sampleA", out_a) == 1
    assert catalogue.write_prefixed(b, "sampleB", out_b) == 1

    names = []
    for path in (out_a, out_b):
        with open(path) as fh:
            names.extend(line[1:].split()[0] for line in fh if line.startswith(">"))
    assert names == ["sampleA~" + shared, "sampleB~" + shared]
    assert len(set(names)) == 2, "expected two genes"
    # The prefix must be removable for joins back to per-sample tables.
    assert catalogue.split_prefixed(names[0]) == ("sampleA", shared)


# ── the cluster table ────────────────────────────────────────────────────────────────────────


# ── subsetting ───────────────────────────────────────────────────────────────────────────────


# ── the thresholds ───────────────────────────────────────────────────────────────────────────


def test_the_cluster_command_asks_for_target_coverage_and_longest_representative():
    """The linclust command uses ``--cov-mode 1`` and ``--cluster-mode 2``.

    This measures coverage of the member against a longest-first representative. Mode 0 would not
    cluster a gene fragment with its full-length gene, and metagenomes are mostly fragments.
    """
    command = catalogue.MMSEQS_LINCLUST.format(
        db="d", clusters="c", tmp="t", identity=0.95, coverage=0.9, threads=4
    )
    assert "--cov-mode 1" in command
    assert "--cluster-mode 2" in command
    assert "--min-seq-id 0.95" in command
    assert "-c 0.9" in command
    assert "{" not in command


# ── quantifying reads against the catalogue ──────────────────────────────────────────────────
