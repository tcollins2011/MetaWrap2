"""One scratch directory policy, and a disk-space check before a long step starts.

**Scratch.** Modules used to scatter their temporary space: some under the output directory,
some via a tool's own ``--tmp-dir``, some wherever the system temp happens to be. On a cluster
that matters - ``/tmp`` is frequently small, node-local, and not where you want a 200 GB SPAdes
intermediate. ``[settings] scratch_dir`` in ``metawrap2.toml`` (or ``$METAWRAP2_SCRATCH``) now
decides for every module; unset, it stays inside the run's output directory, which is the safe
default because that is the one place the user definitely chose.

**Disk space.** metaSPAdes on 34M read pairs needs far more scratch than its final output, and
discovering that by filling the filesystem several hours in is expensive. :func:`check_space`
estimates a step's need from its inputs and refuses to start if the filesystem cannot plausibly
hold it. The estimates are deliberately rough multipliers of input size - being approximately
right before the run is worth much more than being exact afterwards - and
``skip_space_check = true`` under ``[settings]`` overrides them.
"""

from __future__ import annotations

import os
import shutil
from typing import Dict, Iterable, List, Optional, Tuple

#: Rough peak scratch+output need per module, as a multiple of total input size. Derived from
#: observed runs (metaSPAdes is by far the worst: error correction writes a corrected copy of
#: the library, then the graph stages exceed it again).
SPACE_MULTIPLIERS: Dict[str, float] = {
    "assembly": 12.0,
    "read_qc": 4.0,
    "binning": 3.0,
    "reassemble_bins": 6.0,
    "bin_refinement": 2.0,
    "quant_bins": 2.0,
    "kraken2": 1.5,
    "classify_bins": 2.0,
    "annotate_bins": 1.5,
    "blobology": 3.0,
}

#: Never estimate less than this: even a tiny input needs room for tool overhead.
MINIMUM_ESTIMATE_BYTES = 2 * 1024**3  # 2 GB

_SCRATCH_ENV = "METAWRAP2_SCRATCH"


def scratch_root(settings: object, output_dir: str) -> str:
    """The directory a module should put its temporary space in.

    Resolution order: ``$METAWRAP2_SCRATCH``, then ``[settings] scratch_dir``, then a
    ``.metawrap2_tmp`` directory inside *output_dir*. The output directory is the default
    because it is the one location the user explicitly chose, and it is therefore the one
    guaranteed to have been sized for this run.
    """
    configured = os.environ.get(_SCRATCH_ENV) or getattr(settings, "scratch_dir", "") or ""
    if configured:
        return os.path.abspath(os.path.expanduser(configured))
    return os.path.join(os.path.abspath(output_dir), ".metawrap2_tmp")


def scratch_dir(settings: object, output_dir: str, name: str) -> str:
    """Create and return a named scratch subdirectory, e.g. ``metaspades``."""
    path = os.path.join(scratch_root(settings, output_dir), name)
    os.makedirs(path, exist_ok=True)
    return path


def total_size(paths: Iterable[str]) -> int:
    """Total bytes of *paths*, following into directories one level deep."""
    total = 0
    for path in paths:
        if not path:
            continue
        try:
            if os.path.isdir(path):
                for name in os.listdir(path):
                    full = os.path.join(path, name)
                    if os.path.isfile(full):
                        total += os.path.getsize(full)
            elif os.path.isfile(path):
                total += os.path.getsize(path)
        except OSError:
            continue
    return total


def estimate_bytes(module: str, inputs: Iterable[str]) -> int:
    """Rough peak disk need for *module* given *inputs*."""
    multiplier = SPACE_MULTIPLIERS.get(module, 2.0)
    return max(MINIMUM_ESTIMATE_BYTES, int(total_size(inputs) * multiplier))


def free_bytes(path: str) -> Optional[int]:
    """Free space on the filesystem holding *path* (or its nearest existing parent)."""
    probe = os.path.abspath(path)
    while probe and not os.path.exists(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            return None
        probe = parent
    try:
        return shutil.disk_usage(probe).free
    except OSError:
        return None


def human(size: Optional[float]) -> str:
    if size is None:
        return "unknown"
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return "%.1f %s" % (value, unit)
        value /= 1024
    return "%.1f PB" % value


def check_space(module: str, inputs: Iterable[str], targets: Iterable[str]) -> Tuple[bool, str]:
    """Is there plausibly room for *module*? Returns ``(ok, message)``.

    *targets* are the directories that will be written to (output and scratch, which may be on
    different filesystems). Each is checked against the estimate separately, since the run
    needs room on whichever is tightest.
    """
    needed = estimate_bytes(module, inputs)
    problems: List[str] = []
    details: List[str] = []
    checked: Dict[Tuple[int, int], str] = {}
    for target in targets:
        if not target:
            continue
        free = free_bytes(target)
        if free is None:
            continue
        try:
            probe = os.path.abspath(target)
            while probe and not os.path.exists(probe):
                probe = os.path.dirname(probe)
            key = (os.stat(probe).st_dev, 0)
        except OSError:
            key = (len(checked), 0)
        if key in checked:  # same filesystem as something already reported
            continue
        checked[key] = target
        details.append("%s: %s free" % (target, human(free)))
        if free < needed:
            problems.append(
                "%s has %s free but this step may need about %s"
                % (target, human(free), human(needed))
            )
    message = "estimated need ~%s; %s" % (human(needed), "; ".join(details) or "no target checked")
    if problems:
        return False, "\n".join(problems)
    return True, message
