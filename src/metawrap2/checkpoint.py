"""Step checkpointing so a run can resume where it left off.

Each module wraps its major stages in a checkpoint: a stage that finished writes a marker
under ``<output>/.metawrap2/steps/<stage>.done``. With ``--resume`` a re-run skips any stage
whose marker exists (and reuses the output dir); without it, a fresh run clears old markers.

A marker holds ``{"status": "done", "dry_run": false}``. ``--dry-run`` writes no markers, and a
marker that records a dry run is not honoured, since nothing ran. A marker holding the plain text
``done`` (written before 2.1.5) counts as done.

Usage in a module::

    ckpt = make_checkpoint(args.output)
    if ckpt.todo("align_reads"):
        _align_reads(...)
        ckpt.done("align_reads")

Stages are listed explicitly in each module, so the resumable stages are visible in the code.
"""

from __future__ import annotations

import json
import os
import shutil

from . import command

_STEPS_SUBDIR = os.path.join(".metawrap2", "steps")


class Checkpoint:
    def __init__(self, output_dir: str, resume: bool):
        self.dir = os.path.join(output_dir, _STEPS_SUBDIR)
        self.resume = resume

    def _marker(self, step: str) -> str:
        return os.path.join(self.dir, step + ".done")

    def finished(self, step: str) -> bool:
        """True if a real (not ``--dry-run``) run recorded *step* as finished."""
        try:
            with open(self._marker(step)) as fh:
                text = fh.read().strip()
        except OSError:
            return False
        if text == "done":
            return True
        try:
            marker = json.loads(text)
        except ValueError:
            return False
        return (
            isinstance(marker, dict)
            and marker.get("status") == "done"
            and marker.get("dry_run") is False
        )

    def todo(self, step: str) -> bool:
        """True if *step* should run now (always, unless resuming and it already finished)."""
        return not (self.resume and self.finished(step))

    def done(self, step: str) -> None:
        """Record that *step* finished. Writes nothing under ``--dry-run``, where nothing ran."""
        if command.runner.dry_run:
            return
        os.makedirs(self.dir, exist_ok=True)
        with open(self._marker(step), "w") as fh:
            json.dump({"status": "done", "dry_run": False}, fh)
            fh.write("\n")

    def clear(self, step: str) -> None:
        """Forget that *step* finished, so a later --resume runs it again.

        Used for stages whose outputs a module deletes; otherwise a resume would skip the stage
        and pass missing files to the next one. Does nothing under ``--dry-run``.
        """
        if command.runner.dry_run:
            return
        if os.path.exists(self._marker(step)):
            os.remove(self._marker(step))

    def reset(self) -> None:
        """Clear all markers (used for a fresh, non-resume run)."""
        parent = os.path.dirname(self.dir)
        if os.path.isdir(parent):
            shutil.rmtree(parent, ignore_errors=True)


def make_checkpoint(output_dir: str) -> Checkpoint:
    """Build a Checkpoint honoring the run's --resume / --force state.

    With --resume, existing markers are reused. A fresh run clears them so markers from a
    previous run cannot cause steps to be skipped.
    """
    ckpt = Checkpoint(output_dir, resume=command.runner.resume)
    if not command.runner.resume and not command.runner.dry_run:
        ckpt.reset()
    return ckpt
