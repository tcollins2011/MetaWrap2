"""Step checkpointing so a run can resume where it left off.

Each module wraps its major stages in a checkpoint: a stage that finished writes a marker
under ``<output>/.metawrap2/steps/<stage>.done``. With ``--resume`` a re-run skips any stage
whose marker exists (and reuses the output dir); without it, a fresh run clears old markers.

Usage in a module::

    ckpt = make_checkpoint(args.output)
    if ckpt.todo("align_reads"):
        _align_reads(...)
        ckpt.done("align_reads")

Stages are listed explicitly in each module, so the resumable stages are visible in the code.
"""

from __future__ import annotations

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

    def todo(self, step: str) -> bool:
        """True if *step* should run now (always, unless resuming and it already finished)."""
        return not (self.resume and os.path.exists(self._marker(step)))

    def done(self, step: str) -> None:
        """Record that *step* finished."""
        os.makedirs(self.dir, exist_ok=True)
        with open(self._marker(step), "w") as fh:
            fh.write("done\n")

    def clear(self, step: str) -> None:
        """Forget that *step* finished, so a later --resume runs it again.

        Used for stages whose outputs a module deletes; otherwise a resume would skip the stage
        and pass missing files to the next one.
        """
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
