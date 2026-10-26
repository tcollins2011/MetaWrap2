"""Compositional differential abundance over a feature table, using ANCOM-BC.

A relative abundance table is compositional: values sum to one, so one feature rising forces the
others down. A t-test or Wilcoxon on such data tests the constraint as much as the biology, which is
why compositional (CoDA) methods are used instead.

The default is ANCOM-BC from scikit-bio (>= 0.7.1). The other implementations are R packages; this
one needs no R and no extra container.

Scope:

* **No design is inferred.** The caller supplies the grouping. Comparing whichever groups happen to
  be present would give a number for any input, including meaningless comparisons.
* **The multiple-testing correction is recorded per row.** Testing 573 modules or 28,000 KOs at once
  needs one, and the method changes which results pass.
* **The prevalence filter is explicit and off by default.** Rare features have little power, but
  removing them changes the correction for every other feature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Tuple

#: Multiple-testing correction. Holm (family-wise error rate) is the default instead of
#: Benjamini-Hochberg (false discovery rate) because resistome or metabolic findings are usually
#: reported one feature at a time. Both are available; the method used is in every row.
P_ADJUST = "holm"
#: Significance threshold passed to ANCOM-BC for its `Signif` column. Recorded in the output.
ALPHA = 0.05
#: Minimum fraction of samples a feature must appear in to be tested. The default 0 keeps every
#: feature (see the module docstring).
MIN_PREVALENCE = 0.0


class StatsUnavailable(RuntimeError):
    """Raised when scikit-bio is missing; the message says how to install it."""


def _require_skbio():
    try:
        from skbio.stats.composition import ancombc
    except ImportError as problem:  # pragma: no cover - depends on the install
        raise StatsUnavailable(
            "differential abundance needs scikit-bio >= 0.7.1: "
            "pip install 'metawrap2[stats]' (or pip install scikit-bio)"
        ) from problem
    from skbio.stats.composition import ancombc

    return ancombc


@dataclass(frozen=True)
class Result:
    """One feature's test result."""

    feature: str
    #: The coefficient for the term being tested: positive means higher in that group.
    log_fold_change: float
    standard_error: float
    p_value: float
    q_value: float
    significant: bool
    prevalence: float
    term: str


def prevalence(values: Mapping[str, Mapping[str, float]]) -> Dict[str, float]:
    """feature -> fraction of samples it is present in at all."""
    samples = list(values)
    if not samples:
        return {}
    counts: Dict[str, int] = {}
    for sample in samples:
        for feature, value in values[sample].items():
            if value > 0:
                counts[feature] = counts.get(feature, 0) + 1
    return {feature: n / float(len(samples)) for feature, n in counts.items()}


def to_frame(
    values: Mapping[str, Mapping[str, float]],
    min_prevalence: float = MIN_PREVALENCE,
):
    """Return a samples x features DataFrame (the orientation ANCOM-BC takes) and the pseudo-count.

    Zeros are replaced with a pseudo-count because ANCOM-BC needs strictly positive values. Half the
    smallest observed non-zero value is used, the conventional choice. It is recorded in the output
    because a different pseudo-count changes the results for the sparsest features.
    """
    import pandas as pd

    keep = {
        feature for feature, fraction in prevalence(values).items() if fraction >= min_prevalence
    }
    frame = pd.DataFrame(
        {sample: {f: v for f, v in row.items() if f in keep} for sample, row in values.items()}
    ).T.fillna(0.0)
    frame = frame.reindex(sorted(frame.columns), axis=1)
    positive = frame.to_numpy()
    smallest = positive[positive > 0].min() if (positive > 0).any() else 1.0
    return frame.replace(0.0, smallest / 2.0), smallest / 2.0


def test_groups(
    values: Mapping[str, Mapping[str, float]],
    groups: Mapping[str, str],
    *,
    term: Optional[str] = None,
    alpha: float = ALPHA,
    p_adjust: str = P_ADJUST,
    min_prevalence: float = MIN_PREVALENCE,
) -> Tuple[List[Result], Dict[str, object]]:
    """Test every feature for a difference between sample groups.

    *groups* maps sample to group label. The design is ``feature ~ group``; covariates, repeated
    measures and interactions are not supported.

    Returns the per-feature results and the settings used (correction method, prevalence filter,
    pseudo-count), which are needed to compare adjusted p-values across runs.
    """
    ancombc = _require_skbio()
    import pandas as pd

    shared = [sample for sample in values if sample in groups]
    if len(shared) < 3:
        raise ValueError(
            "need at least 3 samples with a group label to test anything; got %d" % len(shared)
        )
    labels = {groups[sample] for sample in shared}
    if len(labels) < 2:
        raise ValueError(
            "need at least two groups to compare; every sample is labelled %r" % labels.pop()
        )

    subset = {sample: values[sample] for sample in shared}
    frame, pseudo = to_frame(subset, min_prevalence)
    metadata = pd.DataFrame({"group": [groups[s] for s in frame.index]}, index=frame.index)

    output = ancombc(frame, metadata, "group", alpha=alpha, p_adjust=p_adjust)
    # scikit-bio 0.7.4 returns an ANCOMBCResult whose .result is a DataFrame indexed by
    # (FeatureID, Covariate). The covariate of a two-level factor is named after the non-reference
    # level ("group[T.B]"), so the term is read from the index. The Intercept row is dropped.
    table = output.result
    covariates = [
        name for name in table.index.get_level_values("Covariate").unique() if name != "Intercept"
    ]
    if not covariates:
        raise ValueError("ANCOM-BC returned only an intercept; there is nothing to compare")
    chosen = term or covariates[0]
    if chosen not in covariates:
        raise ValueError("no such term %r; this design has %r" % (chosen, covariates))
    rows = table.xs(chosen, level="Covariate")

    prevalences = prevalence(subset)
    results = [
        Result(
            feature=str(feature),
            log_fold_change=float(row["Log(FC)"]),
            standard_error=float(row["SE"]),
            p_value=float(row["pvalue"]),
            q_value=float(row["qvalue"]),
            significant=bool(row["Signif"]),
            prevalence=prevalences.get(str(feature), 0.0),
            term=str(chosen),
        )
        for feature, row in rows.iterrows()
    ]
    settings: Dict[str, object] = {
        "test": "ancom-bc",
        "implementation": "scikit-bio",
        "term": chosen,
        "alpha": alpha,
        "p_adjust": p_adjust,
        "min_prevalence": min_prevalence,
        "pseudocount": pseudo,
        "n_samples": len(shared),
        "n_features": int(frame.shape[1]),
        "groups": sorted(labels),
    }
    return results, settings


#: Columns of the differential-abundance table.
RESULT_COLUMNS: Tuple[str, ...] = (
    "feature",
    "term",
    "log_fold_change",
    "standard_error",
    "p_value",
    "q_value",
    "significant",
    "prevalence",
    "test",
    "p_adjust",
    "alpha",
    "pseudocount",
    "caveat",
)

#: Written on every row: why a compositional method is used instead of a t-test.
COMPOSITIONAL_CAVEAT = (
    "relative abundances are compositional: they sum to a constant, so one feature rising pushes "
    "others down and a t-test or Wilcoxon test is biased by that constraint. ANCOM-BC estimates "
    "and corrects the sampling-fraction bias"
)


def result_row(result: Result, settings: Mapping[str, object]) -> Dict[str, object]:
    return {
        "feature": result.feature,
        "term": result.term,
        "log_fold_change": "%.6g" % result.log_fold_change,
        "standard_error": "%.6g" % result.standard_error,
        "p_value": "%.4g" % result.p_value,
        "q_value": "%.4g" % result.q_value,
        "significant": int(result.significant),
        "prevalence": "%.4f" % result.prevalence,
        "test": settings.get("test", ""),
        "p_adjust": settings.get("p_adjust", ""),
        "alpha": settings.get("alpha", ""),
        "pseudocount": "%.6g" % float(str(settings.get("pseudocount", 0.0) or 0.0)),
        "caveat": COMPOSITIONAL_CAVEAT,
    }
