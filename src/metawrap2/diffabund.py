"""Differential abundance over a feature table, done compositionally because the data is.

The pipeline stops at feature tables and every user then writes the same analysis script. The script
is usually wrong in the same way: a relative abundance table is **compositional** - the values are
constrained to sum to one, so one feature rising forces every other down - and a t-test or Wilcoxon
on it answers a question about the constraint rather than about the biology. That is not a subtlety
nobody acts on; it is why the field moved to CoDA methods.

So the default here is **ANCOM-BC**, which scikit-bio ships natively as of 0.7.1. That matters more
than it sounds: the alternative implementations are R packages, and a pipeline that needs R for its
statistics step does not get a statistics step. There is no new language, no new container, and the
dependency is one already in the scientific-Python stack.

What this file deliberately does not do:

* **It does not pick a formula for you.** ANCOM-BC takes a design, and the design is the experiment.
  A module that silently compared "the first two groups it found" would produce a number for any
  input, including inputs where the comparison is meaningless.
* **It does not hide the multiple-testing correction.** 573 modules or 28,000 KOs tested at once
  need one, the choice changes which results survive, and it is recorded per row.
* **It does not drop features quietly.** A feature absent from most samples has no power and will
  not reach significance, but removing it changes every other feature's correction. The prevalence
  filter is explicit, recorded, and off by default.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Tuple

#: Multiple-testing correction. Holm rather than Benjamini-Hochberg as the default because it
#: controls the family-wise error rate rather than the false discovery rate, and a resistome or
#: metabolic claim is usually reported one feature at a time rather than as a set - which is the
#: situation FWER control is for. Both are available; which one ran is in every row.
P_ADJUST = "holm"
#: Significance threshold, passed through rather than applied here: the test reports adjusted
#: p-values and the caller decides. Recorded so a table says what it was filtered at.
ALPHA = 0.05
#: Minimum fraction of samples a feature must appear in to be tested. 0 keeps everything, which is
#: the default - see the module docstring on why filtering is not free.
MIN_PREVALENCE = 0.0


class StatsUnavailable(RuntimeError):
    """Raised when scikit-bio is missing, with what to do about it."""


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
    """A samples x features DataFrame, which is the orientation ANCOM-BC takes.

    Zeros are replaced with a pseudo-count, because ANCOM-BC needs strictly positive values and a
    log-ratio of zero is undefined. Half the smallest observed non-zero value is the conventional
    choice and is recorded in the output: a pseudo-count is an assumption, and a different one gives
    different answers for the sparsest features.
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

    *groups* maps sample to group label; the design is ``feature ~ group``, which is the simplest
    thing ANCOM-BC can be asked and the only one this function offers. Anything more - covariates,
    repeated measures, interactions - is a study design, and a pipeline that guessed at one would be
    inventing the experiment.

    Returns the per-feature results and the settings that produced them, because an adjusted
    p-value without its correction method and its prevalence filter cannot be compared with another.
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
    # (FeatureID, Covariate), with the covariate for a two-level factor named after the non-reference
    # level - "group[T.B]", which is not knowable before the call. So the term is discovered from
    # the index rather than assumed, and the Intercept row is dropped: it tests whether a feature's
    # abundance differs from zero, which is not a question anyone asked.
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

#: On every row. The reason the default is a CoDA method rather than a t-test, stated where a reader
#: of the output will see it.
COMPOSITIONAL_CAVEAT = (
    "relative abundances are compositional: the values are sum-constrained, so one feature rising "
    "forces others down and an ordinary t-test or Wilcoxon tests the constraint rather than the "
    "biology. ANCOM-BC estimates and corrects the sampling-fraction bias instead"
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
