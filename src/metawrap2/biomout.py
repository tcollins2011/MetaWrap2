"""Writing feature tables as BIOM, so they can be read by something other than this pipeline.

Every table MetaWrap2 produces is a TSV of its own design. That is fine for reading and terrible for
onward analysis: a user who wants QIIME 2, phyloseq, scikit-bio or MicrobiomeAnalyst has to write a
converter first, and will write a slightly different one each time. BIOM is the format all of those
already read, nothing in the peer set emits it, and it is little code.

Two choices worth stating:

**BIOM 2.1 (HDF5), not 1.0 (JSON).** The JSON form is a dense text representation of a sparse matrix
and gets large fast - a gene catalogue is hundreds of thousands of features by tens of samples. 2.1
is what current tools read and is written sparse.

**Observations are features and samples are samples**, which is BIOM's convention and the opposite
of how these tables are usually printed. A table with features as columns loaded into a tool
expecting the convention is silently transposed, and a transposed abundance table is not obviously
wrong at a glance - it just gives nonsense answers. So the orientation is set here once rather than
at each call site.

Metadata rides along where there is any: a KO's definition, a module's name and class, the unit a
value is in. That last one matters most - a BIOM file of relative abundances and one of genome
equivalents are the same shape and not comparable, and the unit is the only thing that says which
it is.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, Mapping, Optional

#: What MetaWrap2 writes into the BIOM table's ``generated_by`` field, so a file found later says
#: what produced it.
GENERATED_BY = "MetaWrap2"


class BiomUnavailable(RuntimeError):
    """Raised when biom-format is not installed, with what to do about it."""


def _require_biom():
    try:
        import biom  # noqa: F401
    except ImportError as problem:  # pragma: no cover - depends on the install
        raise BiomUnavailable(
            "writing BIOM needs the biom-format package: pip install 'metawrap2[biom]' "
            "(or pip install biom-format h5py)"
        ) from problem
    return __import__("biom")


def write_biom(
    values: Mapping[str, Mapping[str, float]],
    path: str,
    *,
    feature_metadata: Optional[Mapping[str, Mapping[str, object]]] = None,
    sample_metadata: Optional[Mapping[str, Mapping[str, object]]] = None,
    table_id: Optional[str] = None,
) -> int:
    """Write ``{sample: {feature: value}}`` as a BIOM 2.1 file. Returns the feature count.

    The input is keyed sample-first because that is how the pipeline computes it; the table is built
    feature-first because that is BIOM's orientation. Doing the transpose here is the point - see the
    module docstring on why a silently transposed abundance table is the worst outcome.
    """
    biom = _require_biom()
    import numpy as np

    samples = sorted(values)
    features = sorted({feature for row in values.values() for feature in row})
    if not features:
        raise ValueError("nothing to write: no features in any sample")

    matrix = np.zeros((len(features), len(samples)), dtype=float)
    feature_index = {name: i for i, name in enumerate(features)}
    for column, sample in enumerate(samples):
        for feature, value in values[sample].items():
            matrix[feature_index[feature], column] = value

    table = biom.Table(
        matrix,
        observation_ids=features,
        sample_ids=samples,
        observation_metadata=(
            [dict(feature_metadata.get(name, {})) for name in features]
            if feature_metadata
            else None
        ),
        sample_metadata=(
            [dict(sample_metadata.get(name, {})) for name in samples] if sample_metadata else None
        ),
        table_id=table_id or os.path.basename(path),
    )
    import h5py

    with h5py.File(path, "w") as handle:
        table.to_hdf5(handle, generated_by=GENERATED_BY)
    return len(features)


def long_table_to_values(
    rows: Iterable[Mapping[str, object]],
    sample_column: str = "sample",
    feature_column: str = "feature",
    value_column: str = "value",
) -> Dict[str, Dict[str, float]]:
    """Fold a long TSV's rows into ``{sample: {feature: value}}``.

    Rows whose value does not parse are skipped rather than treated as zero: a missing value and a
    measured zero are different, and BIOM's sparse representation already means absent is zero.
    """
    out: Dict[str, Dict[str, float]] = {}
    for row in rows:
        sample = row.get(sample_column)
        feature = row.get(feature_column)
        if not sample or not feature:
            continue
        try:
            value = float(str(row.get(value_column, "")))
        except (TypeError, ValueError):
            continue
        out.setdefault(str(sample), {})[str(feature)] = value
    return out


def unit_metadata(rows: Iterable[Mapping[str, object]]) -> Dict[str, Dict[str, object]]:
    """Per-sample metadata carrying the unit and what it may be compared across.

    On the sample axis rather than the feature axis because the unit is a property of the whole
    table: every value in a file produced by `quantify_function` is in one unit, and attaching it to
    each feature would imply it could vary between them.
    """
    out: Dict[str, Dict[str, object]] = {}
    for row in rows:
        sample = row.get("sample")
        if not sample or str(sample) in out:
            continue
        unit = row.get("unit")
        if not unit:
            continue
        out[str(sample)] = {
            "unit": str(unit),
            "comparable_across": str(row.get("comparable_across", "")),
            "abundance_kind": str(row.get("abundance_kind", "")),
        }
    return out
