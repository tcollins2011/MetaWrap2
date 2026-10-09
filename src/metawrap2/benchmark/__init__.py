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

The results are committed under ``benchmarks/results/`` so that a change to a default arrives in
review alongside the numbers that justify it, and so a regression is a diff rather than a memory.
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
