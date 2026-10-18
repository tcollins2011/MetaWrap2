"""The manifest is written from a thread pool, and twice over it was not safe to be.

`cohort` and `run -j N` start per-sample steps in a ThreadPoolExecutor, and every start and finish
writes the manifest. A user running 36 water samples at `-j 8` lost a study after 33 finished
assemblies:

    File ".../metawrap2/manifest.py", line 191, in save
        os.replace(tmp, self.path)
    FileNotFoundError: .../.metawrap2/run_steps/manifest.json.tmp

Two races, not one. The temp path was shared, so the first rename took the file and the second
found nothing; and `data` was serialised while another thread was mutating it, which can raise or
write a manifest missing a step that had already started. `--resume` recovered the run, so the cost
was hours rather than results - but the manifest is what resume trusts.
"""

import json
import os
import threading
from concurrent.futures import ThreadPoolExecutor

from metawrap2.manifest import Manifest


def _manifest(tmp_path):
    return Manifest(str(tmp_path / ".metawrap2" / "run_steps" / "manifest.json"))


def test_many_threads_recording_steps_all_survive(tmp_path):
    """The reproducer, scaled down: short steps so the starts land on top of each other."""
    manifest = _manifest(tmp_path)
    names = ["assembly:sample_%d" % i for i in range(48)]
    errors = []

    def record(name):
        try:
            started = manifest.start(name, ["mw2", "assembly"], [])
            manifest.finish(name, started, "completed", outputs=[])
        except Exception as problem:  # noqa: BLE001 - the failure under test is an exception
            errors.append(problem)

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(record, names))

    assert not errors, "manifest writes raced: %r" % errors[:3]
    # Every step is in the file on disk, not just in memory: a resume reads the file.
    with open(manifest.path) as handle:
        written = json.load(handle)
    assert sorted(written["steps"]) == sorted(names)
    assert all(written["steps"][n]["status"] == "completed" for n in names)


def test_no_temp_files_are_left_behind(tmp_path):
    manifest = _manifest(tmp_path)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda i: manifest.start("step:%d" % i, ["x"], []), range(32)))
    strays = [f for f in os.listdir(os.path.dirname(manifest.path)) if f.endswith(".tmp")]
    assert strays == []


def test_the_temp_path_is_unique_per_writer(tmp_path):
    """Directly: two writers must not name the same temp file.

    Asserted on the name rather than only through a race, because a race that happens to pass is
    not evidence, and this is the property the fix relies on.
    """
    manifest = _manifest(tmp_path)
    seen = set()
    real_replace = os.replace

    def capture(src, dst):
        seen.add(src)
        real_replace(src, dst)

    try:
        os.replace = capture  # type: ignore[assignment]
        done = threading.Barrier(4)

        def write(i):
            manifest.start("step:%d" % i, ["x"], [])
            done.wait(timeout=10)

        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(write, range(4)))
    finally:
        os.replace = real_replace  # type: ignore[assignment]

    assert len(seen) == 4, "every thread must write its own temp file, saw %r" % seen


def test_a_failed_write_leaves_the_previous_manifest_readable(tmp_path, monkeypatch):
    """A full disk mid-save must not take the manifest with it - resume reads what was there."""
    manifest = _manifest(tmp_path)
    manifest.start("assembly:sample_0", ["mw2"], [])
    with open(manifest.path) as handle:
        good = handle.read()

    def explode(*_args, **_kwargs):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(json, "dump", explode)
    try:
        manifest.start("assembly:sample_1", ["mw2"], [])
    except OSError:
        pass
    else:
        raise AssertionError("the write error must reach the caller")
    with open(manifest.path) as handle:
        assert handle.read() == good
    assert [f for f in os.listdir(os.path.dirname(manifest.path)) if f.endswith(".tmp")] == []
