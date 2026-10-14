"""Tests for the dereplicate module's clustering, which is the part that decides the catalogue.

skani is a black box here on purpose: everything that turns its table into a set of representative
genomes is a pure function, so what is asserted is the thing that can be got wrong - which genomes
end up collapsed, and which copy of a duplicated genome is the one kept.
"""

import os

from metawrap2.modules import dereplicate


def _edge(ref, query, ani, af_ref=100.0, af_query=100.0):
    return dereplicate.Edge(ref, query, ani, af_ref, af_query)


def _bins(tmp_path, name, bins):
    """A bin folder, each entry (bin name -> sequence length in bases)."""
    directory = tmp_path / name
    directory.mkdir(parents=True)
    for bin_name, length in bins.items():
        (directory / (bin_name + ".fasta")).write_text(">c1\n" + "ACGT" * (length // 4) + "\n")
    return str(directory)


def _stats(path, rows):
    """A CheckM table of the shape this project writes, beside a bin folder."""
    with open(path + ".stats.tsv", "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
        fh.writelines(
            "%s\t%s\t%s\t0.5\tBacteria\t1000\t1000\n" % (name, completeness, contamination)
            for name, (completeness, contamination) in rows.items()
        )


# --- clustering ---------------------------------------------------------------------------


def test_the_larger_aligned_fraction_decides():
    """A half-complete bin aligns fully to its complete neighbour; the neighbour does not to it.

    Taking the smaller of the two fractions would refuse to collapse exactly the duplicates this
    module exists for, so the rule is the larger.
    """
    edges = [_edge("partial", "complete", 99.0, af_ref=98.0, af_query=51.0)]
    assert dereplicate.cluster_genomes(
        ["partial", "complete"], edges, min_aligned_fraction=90.0
    ) == [["partial", "complete"]]
    # But 99% ANI over 5% of a genome is a shared prophage, not the same organism.
    edges = [_edge("a", "b", 99.0, af_ref=5.0, af_query=4.0)]
    assert dereplicate.cluster_genomes(["a", "b"], edges, min_aligned_fraction=10.0) == [
        ["a"],
        ["b"],
    ]


# --- the whole module, with skani stubbed out ---------------------------------------------


def test_the_catalogue_keeps_one_genome_per_cluster(tmp_path, monkeypatch):
    """End to end over a fake skani table: two samples sharing one genome give three genomes."""
    a = _bins(tmp_path, "sampleA_bins", {"bin_001": 4000, "bin_002": 400})
    b = _bins(tmp_path, "sampleB_bins", {"bin_001": 4000})
    _stats(b, {"bin_001": (95.0, 1.0)})
    out = str(tmp_path / "out")
    # No conda env, so the test does not need mamba to be on PATH; skani is stubbed anyway.
    config = tmp_path / "mw2.toml"
    config.write_text("[settings]\nuse_conda_envs = false\n")

    def fake_run(cmd, **kwargs):
        # Write the table skani would have written for these three genomes.
        target = cmd.split("-o ")[1].split()[0]
        with open(target, "w") as fh:
            fh.write(
                "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\t"
                "Ref_name\tQuery_name\n"
            )
            fh.write(
                "%s\t%s\t99.5\t100.0\t100.0\tc1\tc1\n"
                % (os.path.join(a, "bin_001.fasta"), os.path.join(b, "bin_001.fasta"))
            )
        return 0

    monkeypatch.setattr(dereplicate, "run", fake_run)
    argv = ["-b", a, "-b", b, "-o", out, "-t", "2", "--config", str(config)]
    assert dereplicate.main(argv) == 0

    kept = sorted(os.listdir(os.path.join(out, dereplicate.CATALOGUE_DIR)))
    # The scored copy of the shared genome, not sample A's unscored one.
    assert kept == ["sampleA_bins__bin_002.fasta", "sampleB_bins__bin_001.fasta"]

    with open(os.path.join(out, dereplicate.CLUSTERS_TABLE)) as fh:
        rows = fh.read().splitlines()
    assert rows[0].split("\t") == list(dereplicate.CLUSTER_COLUMNS)
    # Every input bin is accounted for, not just the ones that were kept.
    assert len(rows) == 4
