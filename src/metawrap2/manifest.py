"""The run manifest: what each step consumed, produced, and how long it took.

A ``<step>.done`` marker alone cannot show that a step's outputs still exist (or were written at
all, if the tool exited 0 with no output), or that its inputs have not changed since (a re-run of
an earlier step, an edited assembly). The manifest records, per step:

* the command that ran, when it started, and how long it took;
* every declared input with its size, mtime and content fingerprint;
* every declared output with the same;
* the outcome: ``completed``, ``failed``, or ``interrupted``.

With ``--resume``, a step is skipped only if it completed, its outputs are all present and
non-empty, and none of its inputs has changed. Anything else re-runs.

Fingerprints: :func:`fingerprint` hashes the size plus the first and last megabyte, which detects
a file replaced, re-downloaded, truncated or regenerated at constant cost (hashing a 100 GB FASTQ
could take longer than the step). ``--strict-fingerprints`` hashes whole files.
"""

from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

#: How much of each end of a file to hash for the cheap fingerprint.
SAMPLE_BYTES = 1024 * 1024

#: Manifest filename inside a study's state directory.
MANIFEST_NAME = "manifest.json"

#: Manifest schema version, so a manifest in a different format is ignored.
SCHEMA_VERSION = 1

COMPLETED = "completed"
FAILED = "failed"
INTERRUPTED = "interrupted"


def fingerprint(path: str, whole_file: bool = False) -> Optional[str]:
    """A cheap, stable fingerprint of *path*'s contents, or None if it cannot be read.

    Hashes the size plus the first and last :data:`SAMPLE_BYTES`. Different files of equal size
    with identical ends would collide, but that does not happen for the changes this targets
    (regenerated, re-downloaded, truncated, replaced). Pass *whole_file* to hash everything.
    """
    try:
        size = os.path.getsize(path)
    except OSError:
        return None
    digest = hashlib.sha256()
    digest.update(str(size).encode())
    try:
        with open(path, "rb") as fh:
            if whole_file or size <= 2 * SAMPLE_BYTES:
                while True:
                    block = fh.read(SAMPLE_BYTES)
                    if not block:
                        break
                    digest.update(block)
            else:
                digest.update(fh.read(SAMPLE_BYTES))
                fh.seek(-SAMPLE_BYTES, os.SEEK_END)
                digest.update(fh.read(SAMPLE_BYTES))
    except OSError:
        return None
    return digest.hexdigest()[:32]


def describe(path: str, whole_file: bool = False) -> Dict[str, Any]:
    """Record of one file or directory's current state, for later change detection."""
    entry: Dict[str, Any] = {"path": path}
    if os.path.isdir(path):
        names = sorted(n for n in os.listdir(path) if not n.startswith("."))
        entry["kind"] = "directory"
        entry["entries"] = len(names)
        # A directory is fingerprinted by its file names and sizes, which change when e.g. a
        # binner is re-run.
        listing = "\n".join(
            "%s:%d" % (n, os.path.getsize(os.path.join(path, n)))
            for n in names
            if os.path.isfile(os.path.join(path, n))
        )
        entry["fingerprint"] = hashlib.sha256(listing.encode()).hexdigest()[:32]
        return entry
    if not os.path.exists(path):
        entry["kind"] = "missing"
        return entry
    entry["kind"] = "file"
    entry["size"] = os.path.getsize(path)
    entry["mtime"] = round(os.path.getmtime(path), 3)
    entry["fingerprint"] = fingerprint(path, whole_file=whole_file)
    return entry


def unchanged(recorded: Dict[str, Any], whole_file: bool = False) -> bool:
    """True if the file/directory described by *recorded* is still in that state."""
    path = recorded.get("path")
    if not path:
        return False
    current = describe(path, whole_file=whole_file)
    if current.get("kind") != recorded.get("kind"):
        return False
    if current.get("kind") == "missing":
        return False
    return current.get("fingerprint") == recorded.get("fingerprint")


class StepRecord:
    """One step's entry in the manifest."""

    def __init__(self, name: str, data: Dict[str, Any]):
        self.name = name
        self.data = data

    @property
    def status(self) -> str:
        return self.data.get("status", FAILED)

    @property
    def outputs(self) -> List[Dict[str, Any]]:
        return self.data.get("outputs", [])

    @property
    def inputs(self) -> List[Dict[str, Any]]:
        return self.data.get("inputs", [])

    def outputs_present(self) -> Optional[str]:
        """None if every recorded output still exists and is non-empty, else why not."""
        for entry in self.outputs:
            path = entry["path"]
            if entry.get("kind") == "directory":
                if not os.path.isdir(path):
                    return "output directory %s is gone" % path
                if not any(not n.startswith(".") for n in os.listdir(path)):
                    return "output directory %s is now empty" % path
                continue
            if not os.path.isfile(path):
                return "output %s is gone" % path
            if os.path.getsize(path) == 0:
                return "output %s is now empty" % path
        return None

    def inputs_unchanged(self, whole_file: bool = False) -> Optional[str]:
        """None if every recorded input is unchanged, else which one changed."""
        for entry in self.inputs:
            if not unchanged(entry, whole_file=whole_file):
                return "input %s has changed since this step ran" % entry["path"]
        return None


class Manifest:
    """The study's manifest: load, record, query. Written atomically after every change."""

    def __init__(self, path: str):
        self.path = path
        self.data: Dict[str, Any] = {"schema": SCHEMA_VERSION, "steps": {}}
        # `cohort` and `run -j N` start and finish steps from a thread pool; each mutates `data`
        # and saves. Without a lock, threads raced on the temp file (`FileNotFoundError:
        # manifest.json.tmp`, which stopped a study after 33 of 36 assemblies) and json.dump could
        # read `data` mid-update. Re-entrant because the mutators call save() while holding it.
        self._lock = threading.RLock()
        self.load()

    # -- persistence ----------------------------------------------------------------------

    def load(self) -> None:
        if not os.path.isfile(self.path):
            return
        try:
            with open(self.path) as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            # A corrupt manifest does not stop a run; nothing is skipped.
            return
        if data.get("schema") != SCHEMA_VERSION:
            return
        self.data = data
        self.data.setdefault("steps", {})

    def save(self) -> None:
        """Write the manifest. Atomic, and safe to call from several threads at once.

        The temp file name includes the pid and thread id, since with a shared `.tmp` the second
        `os.replace` raises FileNotFoundError. The lock keeps `data` from changing during the dump.
        """
        with self._lock:
            os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
            tmp = "%s.%d.%d.tmp" % (self.path, os.getpid(), threading.get_ident())
            try:
                with open(tmp, "w") as fh:
                    json.dump(self.data, fh, indent=2, sort_keys=True)
                os.replace(tmp, self.path)
            except BaseException:
                # Remove the temp file on an interrupt or full disk; the manifest is untouched.
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise

    # -- recording ------------------------------------------------------------------------

    def start(
        self,
        name: str,
        command: Sequence[str],
        inputs: Iterable[str],
        grant: Optional[Dict[str, Any]] = None,
    ) -> float:
        """Record that *name* is starting. Returns the start time to pass back to finish()."""
        started = time.time()
        # Built before taking the lock, since `describe` reads every input and would block other
        # threads.
        entry: Dict[str, Any] = {
            "status": "running",
            "command": list(command),
            "started": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(started)),
            "inputs": [describe(p) for p in inputs],
            "outputs": [],
        }
        if grant:
            # The resources the step was granted, to compare with its measured usage (e.g. whether
            # a 40 GB peak was the whole grant or a fifth of it). Recorded at start so an
            # interrupted step (e.g. after an OOM) keeps it. See metawrap2.limits and
            # `metawrap2 calibrate`.
            entry["grant"] = dict(grant)
        with self._lock:
            self.data["steps"][name] = entry
            self.save()
        return started

    def finish(
        self,
        name: str,
        started: float,
        status: str,
        outputs: Iterable[str] = (),
        log: Optional[str] = None,
        usage: Optional[Dict[str, Any]] = None,
    ) -> None:
        # Outputs are described before taking the lock, as in start().
        described = [describe(p) for p in outputs]
        with self._lock:
            entry = self.data["steps"].setdefault(name, {})
            entry["status"] = status
            entry["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            entry["duration_seconds"] = round(time.time() - started, 1)
            entry["outputs"] = described
            if log:
                entry["log"] = log
            if usage:
                # Peak memory, CPU time and wall time, used to calibrate budgets and disk estimates.
                entry["usage"] = dict(usage)
            self.save()

    # -- querying -------------------------------------------------------------------------

    def step(self, name: str) -> Optional[StepRecord]:
        data = self.data["steps"].get(name)
        return StepRecord(name, data) if data else None

    def can_skip(self, name: str, whole_file: bool = False) -> Tuple[bool, str]:
        """Whether *name* may be skipped on a resume, and the reason either way.

        A step is skipped only if it completed, all of its outputs are still present and
        non-empty, and none of its inputs has changed.
        """
        record = self.step(name)
        if record is None:
            return False, "no record of it in the manifest"
        if record.status != COMPLETED:
            return False, "last attempt %s" % record.status
        problem = record.outputs_present()
        if problem:
            return False, problem
        problem = record.inputs_unchanged(whole_file=whole_file)
        if problem:
            return False, problem
        return True, "completed in %ss, outputs intact, inputs unchanged" % record.data.get(
            "duration_seconds", "?"
        )

    def summary(self) -> List[Dict[str, Any]]:
        """Every step, in the order it was recorded, for `metawrap2 status`-style output."""
        return [dict(name=name, **data) for name, data in self.data["steps"].items()]

    def declared_paths(self) -> List[str]:
        """Every path any step declared as an output."""
        paths: List[str] = []
        for data in self.data["steps"].values():
            paths.extend(entry["path"] for entry in data.get("outputs", []) if entry.get("path"))
        return sorted(set(paths))

    def audit(
        self, root: Optional[str] = None, ignore: Optional[Iterable[str]] = None
    ) -> Dict[str, List[str]]:
        """Check the manifest against the filesystem. Returns problems by kind.

        Reports at the end of a run what ``can_skip`` would only find on the next one:

        ``missing`` / ``empty``
            A completed step's declared output is absent or zero bytes (typically a tool that
            exited 0 without producing anything).
        ``undeclared``
            A top-level entry under *root* that no step declares. Not an error (logs and state
            directories are there), but ``--resume`` cannot track it. Only the top level is
            scanned, and names in *ignore* (the study's LOGS, CLEAN_READS, ...) are skipped.
        ``no_outputs``
            A completed step that declared no outputs, which ``--resume`` cannot verify.
        """
        problems: Dict[str, List[str]] = {
            "missing": [],
            "empty": [],
            "undeclared": [],
            "no_outputs": [],
        }
        for name, data in self.data["steps"].items():
            if data.get("status") != COMPLETED:
                continue
            outputs = data.get("outputs", [])
            if not outputs:
                problems["no_outputs"].append(name)
                continue
            for entry in outputs:
                path = entry.get("path")
                if not path:
                    continue
                if entry.get("kind") == "directory":
                    if not os.path.isdir(path):
                        problems["missing"].append("%s: %s" % (name, path))
                    elif not any(not n.startswith(".") for n in os.listdir(path)):
                        problems["empty"].append("%s: %s" % (name, path))
                elif not os.path.exists(path):
                    problems["missing"].append("%s: %s" % (name, path))
                elif os.path.getsize(path) == 0:
                    problems["empty"].append("%s: %s" % (name, path))

        if root and os.path.isdir(root):
            declared = {os.path.abspath(p) for p in self.declared_paths()}
            skip = set(ignore or ())
            for entry in sorted(os.listdir(root)):
                if entry.startswith(".") or entry in skip:
                    continue
                full = os.path.abspath(os.path.join(root, entry))
                if full not in declared:
                    problems["undeclared"].append(full)
        return problems
