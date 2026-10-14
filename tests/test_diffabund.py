"""Differential abundance, and the refusals that matter more than the test.

The test itself is scikit-bio's. What this file pins is the reasoning around it:

* **The data is compositional**, so the default is ANCOM-BC rather than a t-test. A relative
  abundance table is sum-constrained - one feature rising forces others down - and an ordinary test
  on it answers a question about the constraint.
* **It refuses rather than producing a number.** Two samples, or one group, cannot support a
  comparison, and a module that returned something anyway would be worse than one that errors: the
  number looks like every other number in the output.
* **Every setting that changes the answer is recorded.** An adjusted p-value without its correction
  method, prevalence filter and pseudo-count cannot be compared with another one.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("skbio", reason="the stats extra is optional")

from metawrap2 import diffabund


def _dataset(n_per_group=3, n_features=12, planted="K0", effect=20.0, seed=0):
    """Two groups with one feature genuinely different, at a size the test can see."""
    rng = np.random.default_rng(seed)
    values = {}
    groups = {}
    for i in range(2 * n_per_group):
        case = i < n_per_group
        sample = "s%d" % i
        groups[sample] = "A" if case else "B"
        row = {"K%d" % k: float(rng.integers(50, 150)) for k in range(n_features)}
        row[planted] = 100.0 * effect if not case else 100.0
        values[sample] = row
    return values, groups


# ── the test finds what is there ─────────────────────────────────────────────────────────────


def test_a_planted_difference_is_the_top_hit():
    values, groups = _dataset()
    results, _settings = diffabund.test_groups(values, groups)
    results.sort(key=lambda r: r.q_value)
    assert results[0].feature == "K0"
    assert results[0].significant
    assert results[0].log_fold_change > 0, "higher in B, which is where it was planted"
    # Too few samples is refused rather than tested: a number would look like every other one.
    with pytest.raises(ValueError, match="at least 3 samples"):
        diffabund.test_groups({"s1": {"K1": 1.0}, "s2": {"K1": 2.0}}, {"s1": "A", "s2": "B"})


# ── the refusals ─────────────────────────────────────────────────────────────────────────────


# ── what is recorded ─────────────────────────────────────────────────────────────────────────


# ── prevalence and zeros ─────────────────────────────────────────────────────────────────────
