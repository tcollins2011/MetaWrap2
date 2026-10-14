"""The gene catalogue, and the name collision that silently halves it.

Almost all of this file is about one bug. :mod:`metawrap2.genecall` names a gene ``<contig>_<n>``,
and every assembler numbers contigs from one, so two samples routinely produce the same gene id.
``mmseqs createdb`` keys on the name, so unprefixed inputs are merged *before clustering runs*:
feeding two identical files of 1,772 genes each produced a catalogue built from 1,772 genes, not
3,544, and reported 1.1% redundancy removed instead of 50.5%. Nothing in the output said so.

The fix is the convention :mod:`metawrap2.modules.coverage` already uses for its contig catalogue -
prefix every id with its source - and the tests here are about the properties that makes necessary:
that a prefix can always be removed again, that two sources cannot share a label, and that the
members table carries the original id so the join back to per-sample tables survives.
"""

from __future__ import annotations

from metawrap2.modules import catalogue


def _write(path, records):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, sequence) for name, sequence in records)
    return str(path)


# ── prefixing ────────────────────────────────────────────────────────────────────────────────


def test_prefixing_makes_colliding_gene_ids_distinct(tmp_path):
    """The bug in one test. The same gene id from two samples must survive as two genes."""
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
    assert len(set(names)) == 2, "the whole point: two genes, not one"
    # A prefix can always be removed again, or every join back to a per-sample table is wrong.
    assert catalogue.split_prefixed(names[0]) == ("sampleA", shared)


# ── the cluster table ────────────────────────────────────────────────────────────────────────


# ── subsetting ───────────────────────────────────────────────────────────────────────────────


# ── the thresholds ───────────────────────────────────────────────────────────────────────────


def test_the_cluster_command_asks_for_target_coverage_and_longest_representative():
    """``--cov-mode 1`` plus ``--cluster-mode 2``: coverage of the member against a
    longest-first representative. Mode 0 would refuse to cluster a gene fragment with the
    full-length gene it is part of, and a metagenome is mostly gene fragments - so the catalogue
    would be inflated with pieces of genes it already contains.
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
