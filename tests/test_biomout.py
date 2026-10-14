"""BIOM output, and the orientation that must not be got wrong.

Everything MetaWrap2 writes is a TSV of its own design, which is fine to read and bad to analyse
further: a user wanting QIIME 2, phyloseq or scikit-bio writes a converter first, slightly
differently each time. BIOM is what all of those already read.

The thing worth testing is the transpose. BIOM's convention is observations-as-features and
samples-as-samples, which is the opposite of how these tables print, and a table loaded in the wrong
orientation is not obviously wrong - it has the right numbers in the wrong places and gives nonsense
answers quietly. The unit metadata matters for the same reason: a BIOM file of relative abundances
and one of genome equivalents are the same shape and are not comparable.
"""

from __future__ import annotations

import pytest

biom = pytest.importorskip("biom", reason="the biom extra is optional")

from metawrap2 import biomout


def test_features_become_observations_and_samples_stay_samples(tmp_path):
    """The transpose, pinned. Getting it backwards gives right numbers in wrong places."""
    path = str(tmp_path / "t.biom")
    written = biomout.write_biom({"s1": {"K1": 1.0, "K2": 2.0}, "s2": {"K1": 3.0}}, path)
    assert written == 2

    table = biom.load_table(path)
    assert table.shape == (2, 2), "(features, samples)"
    assert sorted(table.ids("observation")) == ["K1", "K2"]
    assert sorted(table.ids("sample")) == ["s1", "s2"]
    assert table.get_value_by_ids("K2", "s1") == pytest.approx(2.0)
    assert table.get_value_by_ids("K1", "s2") == pytest.approx(3.0)
    assert table.get_value_by_ids("K2", "s2") == pytest.approx(0.0), "absent is zero"


# ── folding a long table ─────────────────────────────────────────────────────────────────────
