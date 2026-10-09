import os
import sys

import pytest

# Allow tests to import the package from src/ without an install step.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


@pytest.fixture(autouse=True)
def _conda_env_policy(request, tmp_path_factory, monkeypatch):
    """Pin ``use_conda_envs`` for every test instead of inheriting the machine's config.

    Most tests fake the tool call and want no env wrapping at all; a handful are *about* how an
    env is resolved or required, and mark themselves ``@pytest.mark.conda_envs``. Left to the
    ambient config, both groups pass on one developer's machine and fail on another's: whichever
    way ``~/.metawrap2/config.toml`` or ``$METAWRAP2_CONFIG`` happens to be set decides.
    """
    wanted = request.node.get_closest_marker("conda_envs") is not None
    path = tmp_path_factory.mktemp("config") / "metawrap2.toml"
    path.write_text("[settings]\nuse_conda_envs = %s\n" % ("true" if wanted else "false"))
    monkeypatch.setenv("METAWRAP2_CONFIG", str(path))


#: The machine every test runs on, whatever host it is actually on. 32 cores, 128 GB, no GPU, and
#: 1.8 TB free of 2 TB. Chosen to be comfortably larger than any limit a test exercises, so a test
#: about a thread cap or a disk ceiling sees the cap rather than the host.
TEST_MACHINE_CORES = 32
TEST_MACHINE_MEMORY = 128 * 1024**3
TEST_DISK = (2 * 1024**4, 200 * 1024**3, 2 * 1024**4 - 200 * 1024**3)  # total, used, free


@pytest.fixture(autouse=True)
def _fixed_machine(request, monkeypatch):
    """Every test sees the same cores, memory, GPUs and free disk, on a laptop or on a CI runner.

    Several modules clamp to what the machine has (threads to its cores, a download to its free
    disk), so a test that passed on a 12-core laptop with a large disk failed on a 4-core CI runner
    with a small one - after the pre-push gate, which runs on the laptop, had passed it. This
    removes the host from the result: the suite now gives the same answer everywhere, so a green
    pre-push means a green CI.

    A test that measures the real host - limits.machine_cores() itself - marks itself
    ``@pytest.mark.real_machine``. A test that needs a different machine (8 cores, a full disk)
    monkeypatches on top of this, which wins because it is applied later.
    """
    if request.node.get_closest_marker("real_machine") is not None:
        return
    import collections
    import shutil

    from metawrap2 import limits

    monkeypatch.setattr(limits, "machine_cores", lambda: TEST_MACHINE_CORES)
    monkeypatch.setattr(limits, "machine_memory_bytes", lambda: TEST_MACHINE_MEMORY)
    monkeypatch.setattr(limits, "machine_gpus", lambda: 0)
    usage = collections.namedtuple("usage", "total used free")(*TEST_DISK)
    monkeypatch.setattr(shutil, "disk_usage", lambda _path: usage)
