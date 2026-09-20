import os

import pytest

from metawrap2 import command
from metawrap2.checkpoint import Checkpoint, make_checkpoint


@pytest.fixture(autouse=True)
def _reset_runner():
    command.runner.configure(resume=False, force=False, dry_run=False)
    yield
    command.runner.configure(resume=False, force=False, dry_run=False)


def test_todo_and_done(tmp_path):
    c = Checkpoint(str(tmp_path), resume=True)
    assert c.todo("align") is True
    c.done("align")
    assert os.path.exists(tmp_path / ".metawrap2" / "steps" / "align.done")
    assert c.todo("align") is False  # resuming skips a finished step
    assert c.todo("binning") is True  # unfinished step still runs


def test_no_resume_always_todo(tmp_path):
    c = Checkpoint(str(tmp_path), resume=False)
    c.done("align")
    assert c.todo("align") is True  # without resume, always run (markers still written)


def test_make_checkpoint_fresh_run_clears_markers(tmp_path):
    # a finished marker from a previous run...
    old = Checkpoint(str(tmp_path), resume=True)
    old.done("align")
    assert os.path.exists(tmp_path / ".metawrap2" / "steps" / "align.done")
    # a fresh (non-resume) run wipes markers so nothing is wrongly skipped
    command.runner.configure(resume=False, dry_run=False)
    c = make_checkpoint(str(tmp_path))
    assert not os.path.exists(tmp_path / ".metawrap2" / "steps" / "align.done")
    assert c.todo("align") is True


def test_make_checkpoint_resume_keeps_markers(tmp_path):
    old = Checkpoint(str(tmp_path), resume=True)
    old.done("align")
    command.runner.configure(resume=True)
    c = make_checkpoint(str(tmp_path))
    assert c.todo("align") is False  # resume preserves and honors the marker
