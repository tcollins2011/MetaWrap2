"""Run provenance for MetaWrap2.

Every module writes three files into its ``-o`` output directory, alongside the data:

* ``run_config.json`` - the full, machine-readable run record: module, MetaWrap2 version,
  resolved parameters (options *and* defaults), config file, conda env, databases, platform,
  start/end/duration, inputs, the exact commands run, and the complete ``lineage``.
* ``run_commands.txt`` - the exact commands that were run, in order. The first line is the
  top-level ``metawrap2 <module> ...`` invocation (re-running it reproduces the output 1:1);
  the remaining lines are the individual tool commands that invocation expanded to.
* ``provenance.txt`` - a human-readable history of every MetaWrap2 step that led to this
  output. When an input is itself a MetaWrap2 output directory, that input's lineage is
  merged in, so provenance propagates across successive runs (and across multiple MetaWrap2
  inputs, de-duplicated by run id). This lets you trace a final bin all the way back to QC.

The runner (:mod:`metawrap2.command`) records each external command into the active recorder;
modules call :func:`start_run` at the top of ``main`` and :func:`finish_run` in a ``finally``.
"""

from __future__ import annotations

import json
import os
import platform
import shlex
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from . import dbcheck as _dbcheck
from .config import Settings as _Settings

CONFIG_NAME = "run_config.json"
ENVIRONMENT_NAME = "run_environment.json"
COMMANDS_NAME = "run_commands.txt"
PROVENANCE_NAME = "provenance.txt"

# How far up from an input path to look for a MetaWrap2 output root.
_MAX_WALK_UP = 8


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def conda_env_packages(env: Optional[str], env_manager: str = "mamba") -> List[str]:
    """Every package in conda env *env* as ``name=version=build``, or [] if unavailable.

    Recorded so a run can be reproduced exactly, and so "it worked last month" can be
    compared against "it fails today" at the level of individual package builds - which is
    where these failures actually come from (an unpinned solve picking a new Python, a new
    setuptools removing pkg_resources).
    """
    if not env:
        return []
    try:
        out = subprocess.run(
            [env_manager, "list", "-n", env, "--json"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=120,
            check=False,
        )
        if out.returncode != 0 or not out.stdout:
            return []
        parsed = json.loads(out.stdout)
    except (OSError, subprocess.SubprocessError, ValueError):
        return []
    # conda returns a bare list; mamba 2.x wraps it as {"packages": [...]}.
    packages = parsed.get("packages", []) if isinstance(parsed, dict) else parsed
    if not isinstance(packages, list):
        return []
    described = []
    for pkg in packages:
        if isinstance(pkg, dict) and pkg.get("name"):
            described.append(
                "%s=%s=%s"
                % (
                    pkg.get("name"),
                    pkg.get("version", "?"),
                    pkg.get("build_string") or pkg.get("build", "?"),
                )
            )
    return sorted(described)


def source_revision() -> Dict[str, str]:
    """Which MetaWrap2 *tree* this is, when that is knowable.

    ``metawrap2_version`` says "2.1.0", and between releases that names many different trees -
    including one with local edits, which is the normal state of a package whose whole design
    invites you to edit the module files. So when MetaWrap2 is running from a git checkout, the
    commit is recorded too, along with whether the working tree was dirty at the time.

    An installed wheel has no repository, and that is not an error: the answer is simply
    ``{"kind": "installed"}``, which is itself informative - it means the version string is the
    whole truth.
    """
    package_dir = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(package_dir))  # src/metawrap2 -> repo root
    if not os.path.isdir(os.path.join(repo, ".git")):
        return {"kind": "installed"}

    def git(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(
                ["git", "-C", repo, *args],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                timeout=10,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    commit = git("rev-parse", "HEAD")
    if not commit:
        return {"kind": "installed"}
    revision = {"kind": "git", "commit": commit, "repo": repo}
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    if branch:
        revision["branch"] = branch
    described = git("describe", "--tags", "--always", "--dirty")
    if described:
        revision["describe"] = described
    status = git("status", "--porcelain")
    # "dirty" is the part that matters: it says the commit alone does not reproduce this run.
    revision["dirty"] = "true" if status else "false"
    if status:
        revision["modified_files"] = str(len(status.splitlines()))
    return revision


@dataclass
class RunRecorder:
    """Collects everything about one module run and writes the provenance files."""

    module: str
    version: str
    parameters: Dict[str, Any]
    config_file: Optional[str]
    conda_env: Optional[str]
    databases: Dict[str, str]
    inputs: List[str]
    command_line: List[str] = field(default_factory=lambda: list(sys.argv))
    platform: str = field(default_factory=platform.platform)
    run_id: str = field(
        default_factory=lambda: time.strftime("%Y%m%dT%H%M%S-") + uuid.uuid4().hex[:8]
    )
    started: str = field(default_factory=_now)
    commands: List[str] = field(default_factory=list)
    environment: Dict[str, Any] = field(default_factory=dict)
    _t0: float = field(default_factory=time.time)

    def elapsed(self) -> float:
        """Seconds since this run started."""
        return time.time() - self._t0

    def record_command(self, cmd: str) -> None:
        """Called by the runner for each executed command (already shell-quoted)."""
        self.commands.append(cmd)

    def record_environment(self, output_dir: str) -> None:
        """Write what this run is about to use, before any of it runs.

        ``finalize`` records tool versions too, but only if the run reaches the end. A run that
        dies in its first tool used to leave no record of *what* it was running - which is
        exactly the run whose environment you want to inspect. This writes the environment
        record up front and is overwritten with the same information plus outcomes at the end.
        """
        self.environment = {
            "captured": "at start",
            "conda_env": self.conda_env,
            "conda_packages": conda_env_packages(self.conda_env),
            "platform": self.platform,
            "python": sys.version.split()[0],
            "metawrap2_version": self.version,
            "metawrap2_source": source_revision(),
            "databases": self.databases,
            # Not just the paths: what was actually at those paths. A database updated in place
            # leaves every earlier run's record naming a path whose contents have changed, with
            # nothing to show that. The fingerprint is what makes "which runs used the old nt?"
            # answerable. See metawrap2.dbcheck.
            "database_contents": _dbcheck.describe_all(_Settings(databases=dict(self.databases))),
        }
        try:
            os.makedirs(output_dir, exist_ok=True)
            with open(os.path.join(output_dir, ENVIRONMENT_NAME), "w") as fh:
                json.dump(self.environment, fh, indent=2, sort_keys=True)
        except OSError:
            pass

    # -- lineage assembly -----------------------------------------------------------------

    def _self_step(self) -> Dict[str, Any]:
        return {
            "run_id": self.run_id,
            "module": self.module,
            "metawrap2_version": self.version,
            "started": self.started,
            "finished": _now(),
            "duration_seconds": round(time.time() - self._t0, 1),
            "command_line": self.command_line,
            "commands": list(self.commands),
            "inputs": list(self.inputs),
        }

    def _lineage(self, input_paths: List[str]) -> List[Dict[str, Any]]:
        ancestors = collect_ancestor_lineage(input_paths)
        return merge_lineage(ancestors, self._self_step())

    # -- output ---------------------------------------------------------------------------

    def finalize(
        self,
        output_dir: str,
        input_paths: Optional[List[str]] = None,
        capture_versions: bool = True,
    ) -> None:
        os.makedirs(output_dir, exist_ok=True)
        input_paths = input_paths if input_paths is not None else self.inputs
        lineage = self._lineage(input_paths)
        tool_versions = capture_tool_versions(self.commands) if capture_versions else {}

        config = {
            "run_id": self.run_id,
            "metawrap2_version": self.version,
            "module": self.module,
            "command_line": self.command_line,
            "started": self.started,
            "finished": _now(),
            "duration_seconds": round(time.time() - self._t0, 1),
            "parameters": self.parameters,
            "config_file": self.config_file,
            "conda_env": self.conda_env,
            "databases": self.databases,
            "platform": self.platform,
            "tool_versions": tool_versions,
            "environment": self.environment,
            "inputs": input_paths,
            "commands": list(self.commands),
            "lineage": lineage,
        }
        with open(os.path.join(output_dir, CONFIG_NAME), "w") as fh:
            json.dump(config, fh, indent=2)

        self._write_commands(os.path.join(output_dir, COMMANDS_NAME))
        _write_provenance(os.path.join(output_dir, PROVENANCE_NAME), lineage)

    def _write_commands(self, path: str) -> None:
        with open(path, "w") as fh:
            fh.write(
                "# MetaWrap2 %s - module: %s - run: %s\n" % (self.version, self.module, self.run_id)
            )
            fh.write("# Re-running the top-level command below reproduces this output:\n")
            fh.write(shlex.join(self.command_line) + "\n\n")
            fh.write("# Individual commands executed (in order):\n")
            fh.writelines(cmd + "\n" for cmd in self.commands)


# -- lineage helpers ----------------------------------------------------------------------


def find_mw2_root(path: str) -> Optional[str]:
    """Walk up from *path* looking for a MetaWrap2 output dir (one with run_config.json)."""
    if not path:
        return None
    cur = os.path.abspath(path)
    if os.path.isfile(cur) or not os.path.isdir(cur):
        cur = os.path.dirname(cur)
    for _ in range(_MAX_WALK_UP):
        if os.path.isfile(os.path.join(cur, CONFIG_NAME)):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return None


def collect_ancestor_lineage(input_paths: List[str]) -> List[Dict[str, Any]]:
    """Gather the lineage steps from any input paths that are MetaWrap2 outputs."""
    steps: List[Dict[str, Any]] = []
    seen_roots = set()
    for path in input_paths or []:
        root = find_mw2_root(path)
        if not root or root in seen_roots:
            continue
        seen_roots.add(root)
        try:
            with open(os.path.join(root, CONFIG_NAME)) as fh:
                cfg = json.load(fh)
        except (OSError, ValueError):
            continue
        steps.extend(cfg.get("lineage", []))
    return steps


def merge_lineage(ancestors: List[Dict[str, Any]], step: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Merge ancestor steps (de-duplicated by run_id, ordered by start time) then append step."""
    by_id: Dict[str, Dict[str, Any]] = {}
    for s in ancestors:
        rid = s.get("run_id")
        if rid and rid not in by_id:
            by_id[rid] = s
    ordered = sorted(by_id.values(), key=lambda s: str(s.get("started", "")))
    if step.get("run_id") not in by_id:
        ordered.append(step)
    return ordered


def _parse_env_tool(recorded: str) -> Optional[Tuple[Optional[str], str]]:
    """From a recorded command, return (env, tool). Skips bash -c pipelines."""
    try:
        toks = shlex.split(recorded)
    except ValueError:
        return None
    if "bash" in toks[:4] and "-c" in toks[:5]:
        return None
    if "run" in toks and "-n" in toks:  # mamba run --no-capture-output -n ENV TOOL ...
        i = toks.index("-n")
        if i + 2 < len(toks):
            return toks[i + 1], toks[i + 2]
        return None
    # plain command (no env)
    for t in toks:
        if t in ("cd", "&&"):
            continue
        return None, t
    return None


def capture_tool_versions(
    commands: List[str], env_manager: str = "mamba", timeout: int = 20
) -> Dict[str, str]:
    """Best-effort: record the version of each distinct tool actually invoked.

    Runs ``<tool> --version`` (falling back to ``-version``) inside the tool's conda env.
    Silently skips anything that can't be probed (missing env, no version flag, no mamba),
    so this never fails a run - it only enriches provenance when possible.
    """
    versions: Dict[str, str] = {}
    seen = set()
    for recorded in commands:
        parsed = _parse_env_tool(recorded)
        if not parsed:
            continue
        env, tool = parsed
        if tool in seen or tool in ("cd", "bash", "python", "python3"):
            continue
        seen.add(tool)
        prefix = [env_manager, "run", "-n", env] if env else []
        for flag in ("--version", "-version"):
            try:
                p = subprocess.run(
                    prefix + [tool, flag],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout,
                    check=False,
                )
            except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
                break
            out = (p.stdout or "").strip()
            if p.returncode == 0 and out:
                versions[tool] = out.splitlines()[0][:120]
                break
    return versions


def _write_provenance(path: str, lineage: List[Dict[str, Any]]) -> None:
    lines = ["MetaWrap2 provenance", "=" * 60, ""]
    lines.append("This output was produced by the following MetaWrap2 steps, in order")
    lines.append("(steps inherited from MetaWrap2 inputs are included):")
    lines.append("")
    for i, step in enumerate(lineage, 1):
        lines.append("[%d] %s   (run %s)" % (i, step.get("module", "?"), step.get("run_id", "?")))
        lines.append("    started:  %s" % step.get("started", "?"))
        lines.append("    finished: %s" % step.get("finished", "?"))
        lines.append("    command:  %s" % shlex.join(step.get("command_line", [])))
        cmds = step.get("commands", [])
        if cmds:
            lines.append("    tools run:")
            for c in cmds:
                lines.append("      $ %s" % c)
        lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
