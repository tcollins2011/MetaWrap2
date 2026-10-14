import os

import pytest

from metawrap2 import command
from metawrap2.checkpoint import Checkpoint, make_checkpoint


@pytest.fixture(autouse=True)
def _reset_runner():
    command.runner.configure(resume=False, force=False, dry_run=False)
    yield
    command.runner.configure(resume=False, force=False, dry_run=False)


def test_make_checkpoint_resume_keeps_markers(tmp_path):
    old = Checkpoint(str(tmp_path), resume=True)
    old.done("align")
    command.runner.configure(resume=True)
    c = make_checkpoint(str(tmp_path))
    assert c.todo("align") is False  # resume preserves and honors the marker
    assert os.path.exists(tmp_path / ".metawrap2" / "steps" / "align.done")
    assert c.todo("binning") is True  # unfinished step still runs
