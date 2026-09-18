"""Shared helpers for MetaWrap2 modules.

Deliberately tiny. Modules import the banner helpers, the conda-env lookup, and a couple of
filesystem/read utilities from here; everything else stays inline in the module so a reader
sees exactly what runs. This is not a framework — just the few things every module repeats.
"""

from __future__ import annotations

import os
import re
from typing import List, Optional

from .. import command as _cmd
from .. import logging as _log
from .. import __version__
from ..checkpoint import make_checkpoint
from ..config import MODULE_ENVS, Settings
from ..provenance import RunRecorder
from ..command import ToolError, run  # re-exported for modules

# Banner logging (same look as the original comm/announcement/warning/error).
comm = _log.comm
announcement = _log.announcement
warning = _log.warning
error = _log.error

__all__ = [
    "comm", "announcement", "warning", "error", "run", "ToolError",
    "env_for", "ensure_dir", "collect_read_pairs", "collect_reads",
    "start_run", "finish_run", "make_checkpoint",
]


def env_for(module: str, settings: Settings) -> Optional[str]:
    """Conda env name for *module*, or None if the user disabled conda envs."""
    if not settings.use_conda_envs:
        return None
    return MODULE_ENVS.get(module)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def start_run(module: str, args, env: Optional[str], settings: Settings,
              inputs: List[str]) -> RunRecorder:
    """Begin provenance recording for a module run and wire it into the runner.

    ``inputs`` are the paths this run reads (assembly/bins/reads); any that live inside a
    previous MetaWrap2 output dir let provenance propagate. Call :func:`finish_run` in a
    ``finally`` so even a failed run leaves a record.
    """
    output = getattr(args, "output", None)
    # Refuse to overwrite a previous MetaWrap2 run unless --force/--resume (or --dry-run).
    if (output and os.path.isfile(os.path.join(output, "run_config.json"))
            and not _cmd.runner.force and not _cmd.runner.resume and not _cmd.runner.dry_run):
        error("Output directory '%s' already contains a MetaWrap2 run. Re-run with "
              "--resume to continue it, --force to overwrite, or choose a new -o." % output)

    params = {k: v for k, v in vars(args).items()}
    inputs = [p for p in inputs if p]
    rec = RunRecorder(
        module=module,
        version=__version__,
        parameters=params,
        config_file=getattr(args, "config", None),
        conda_env=env,
        databases=dict(settings.databases),
        inputs=inputs,
    )
    _cmd.set_recorder(rec)

    # Tee all command output to run.stdout/run.stderr in the output dir (screen unaffected).
    if output:
        os.makedirs(output, exist_ok=True)
        _cmd.set_run_logs(os.path.join(output, "run.stdout"), os.path.join(output, "run.stderr"))
    _warn_odd_paths(inputs + ([output] if output else []))
    return rec


def finish_run(rec: RunRecorder, output_dir: str, inputs: Optional[List[str]] = None) -> None:
    """Write the provenance files (run_config.json, run_commands.txt, provenance.txt)."""
    try:
        rec.finalize(output_dir, input_paths=inputs, capture_versions=not _cmd.runner.dry_run)
    finally:
        _cmd.set_recorder(None)
        _cmd.set_run_logs(None, None)


# Characters in a path that tend to break downstream shell pipelines / tools.
_ODD_PATH = re.compile(r"[\s'\"$`|;&<>()\\!*?]")


def _warn_odd_paths(paths: List[str]) -> None:
    for p in paths:
        if p and _ODD_PATH.search(p):
            warning("Path contains characters that can break some tools (spaces or shell "
                    "metacharacters): %s. Consider renaming to [A-Za-z0-9._/-]." % p)


def collect_read_pairs(reads: List[str]) -> List[tuple]:
    """From a positional list of fastq files, pair up *_1.fastq(.gz) with *_2.fastq(.gz).

    Returns a list of (sample_name, r1, r2). Raises ValueError with a clear message if a
    mate is missing. Compression is fine — the aligners read .gz directly.
    """
    pairs = []
    for r1 in reads:
        if "_1.fastq" not in os.path.basename(r1):
            continue
        r2 = r1.replace("_1.fastq", "_2.fastq")
        if not os.path.isfile(r1):
            raise ValueError("%s does not exist" % r1)
        if not os.path.isfile(r2):
            raise ValueError("expected mate %s for %s but it does not exist" % (r2, r1))
        sample = os.path.basename(r1).split("_1.fastq")[0]
        pairs.append((sample, r1, r2))
    if not pairs:
        raise ValueError(
            "No paired read files found. Expected files named *_1.fastq / *_2.fastq "
            "(.gz allowed). For single-end/interleaved data use --single-end/--interleaved."
        )
    return pairs


def collect_reads(reads: List[str]) -> List[tuple]:
    """From a positional list of fastq files, return (sample_name, path) for each *.fastq(.gz)."""
    out = []
    for r in reads:
        base = os.path.basename(r)
        if ".fastq" not in base and ".fq" not in base:
            continue
        if not os.path.isfile(r):
            raise ValueError("%s does not exist" % r)
        sample = base.split(".fastq")[0].split(".fq")[0]
        out.append((sample, r))
    if not out:
        raise ValueError("No read files found (expected *.fastq / *.fq, .gz allowed).")
    return out
