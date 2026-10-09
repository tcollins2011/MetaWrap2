"""Run the test suite as if on a small, nearly full machine: 2 cores, 4 GB of RAM, 95% disk used.

This is how both the pre-push hook and CI run pytest. The suite is meant to be independent of the
host (tests/conftest.py pins a fixed test machine), and this is what checks that promise. A test
that quietly reads the real core count or free disk passes on a big laptop and fails on a 4-core CI
runner; under this script it fails on the laptop too, before the push, where it is cheap to fix.

    python tests/run_on_small_machine.py [pytest args...]
"""

import collections
import os
import shutil
import sys

_real_sysconf = os.sysconf


def _sysconf(name):
    # 4 GB of physical memory, at the usual 4 KiB page.
    if name == "SC_PHYS_PAGES":
        return (4 * 1024**3) // 4096
    if name == "SC_PAGE_SIZE":
        return 4096
    return _real_sysconf(name)


_DiskUsage = collections.namedtuple("usage", "total used free")

os.cpu_count = lambda: 2
os.sched_getaffinity = lambda _pid: {0, 1}
os.sysconf = _sysconf
shutil.disk_usage = lambda _path: _DiskUsage(100 * 1024**3, 95 * 1024**3, 5 * 1024**3)

if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main(sys.argv[1:] or ["-q"]))
