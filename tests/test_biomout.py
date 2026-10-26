"""BIOM output and its orientation.

MetaWrap2 writes its own TSVs; BIOM is read by QIIME 2, phyloseq and scikit-bio.

BIOM puts features as observations and samples as samples, the transpose of how these tables print.
A table in the wrong orientation has the right numbers in the wrong places and is hard to spot. The
unit metadata is checked because relative abundances and genome equivalents have the same shape.
"""

from __future__ import annotations

import pytest

biom = pytest.importorskip("biom", reason="the biom extra is optional")

from metawrap2 import biomout


def test_features_become_observations_and_samples_stay_samples(tmp_path):
    """Features are BIOM observations and samples are BIOM samples."""
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
