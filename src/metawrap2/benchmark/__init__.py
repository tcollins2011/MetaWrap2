"""Benchmarking: measuring what MetaWrap2's choices actually cost and actually recover.

Every default in MetaWrap2 is supposed to be chosen on measured quality, with measured cost beside
it. Without a
harness that is an aspiration, because nobody re-runs three binners on ten samples by hand to
settle an argument about a flag. This package is the machinery that makes the measurement routine:

- :mod:`~metawrap2.benchmark.registry` - which datasets exist, what they cost to fetch, and where
  their gold standards live, declared as data in ``benchmarks/datasets.toml``.
- :mod:`~metawrap2.benchmark.fetch` - getting them, checksum-verified, without overrunning the
  ``[limits]`` disk ceiling.
- :mod:`~metawrap2.benchmark.runner` - running one comparison's variants over one dataset under the
  same budget, recording wall time and peak memory per step.
- :mod:`~metawrap2.benchmark.scoring` - turning tool output into one tidy TSV of
  ``dataset x tool x mode x metric``, and into the accuracy-per-CPU-hour ratio the defaults are
  actually decided on.

Results append to ``benchmarks/results/<dataset>.tsv``, so a change to a default can be reviewed
alongside the numbers that justify it.
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
