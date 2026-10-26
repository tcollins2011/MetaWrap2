"""Benchmarking: the cost and recovery of MetaWrap2's tool and parameter choices.

MetaWrap2 defaults are chosen on measured quality and measured cost. This package runs those
measurements:

- :mod:`~metawrap2.benchmark.registry` - which datasets exist, what they cost to fetch, and where
  their gold standards live, declared as data in ``benchmarks/datasets.toml``.
- :mod:`~metawrap2.benchmark.fetch` - getting them, checksum-verified, without overrunning the
  ``[limits]`` disk ceiling.
- :mod:`~metawrap2.benchmark.runner` - running one comparison's variants over one dataset under the
  same budget, recording wall time and peak memory per step.
- :mod:`~metawrap2.benchmark.scoring` - turning tool output into one tidy TSV of
  ``dataset x tool x mode x metric``, and into the accuracy-per-CPU-hour ratio used to choose
  defaults.

Results append to ``benchmarks/results/<dataset>.tsv``, so a change to a default can be reviewed
with its numbers.
"""

from __future__ import annotations

from .registry import Dataset, RegistryError, Remote, load_configurations, load_datasets

__all__ = [
    "Dataset",
    "RegistryError",
    "Remote",
    "load_configurations",
    "load_datasets",
]
