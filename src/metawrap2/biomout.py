"""Writing feature tables as BIOM, for QIIME 2, phyloseq, scikit-bio and MicrobiomeAnalyst.

MetaWrap2's own tables are TSVs in its own layout; BIOM is the format those tools read.

**BIOM 2.1 (HDF5), not 1.0 (JSON).** The JSON form stores the sparse matrix densely as text, which
is large for a gene catalogue (hundreds of thousands of features by tens of samples). 2.1 is
sparse and is what current tools read.

**Observations are features and samples are samples**, BIOM's convention and the transpose of how
these tables are usually printed. A transposed abundance table is not obviously wrong, so the
orientation is set here once.

Metadata is attached where available: a KO's definition, a module's name and class, and the unit
of the values. The unit distinguishes, for example, relative abundances from genome equivalents,
which have the same shape.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, Mapping, Optional

#: Written into the BIOM table's ``generated_by`` field.
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

    The input is keyed sample-first, as the pipeline computes it; the table is built
    feature-first, BIOM's orientation. See the module docstring.
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

    Rows whose value does not parse are skipped, not stored as zero; BIOM's sparse representation
    already treats absent as zero.
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

    Stored on the sample axis because the unit applies to the whole table (every value from
    `quantify_function` has one unit).
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
