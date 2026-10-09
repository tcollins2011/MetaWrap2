"""Tests for metawrap2.limits: one resource budget for a whole run."""

from __future__ import annotations

import os

import pytest

from metawrap2 import limits


@pytest.fixture(autouse=True)
def _no_ambient_overrides():
    # Overrides are process-global by design (they come from flags parsed before any module runs),
    # so a test that sets one must not leak it into the next. Same for the cached GPU count.
    limits.clear_overrides()
    limits._count_gpus.cache_clear()
    yield
    limits.clear_overrides()
    limits._count_gpus.cache_clear()


def machine(cores: int = 64, memory_gb: int = 256, gpus: int = 0, monkeypatch=None):
    """Pin the measured machine so the arithmetic is testable on any host."""
    monkeypatch.setattr(limits, "machine_cores", lambda: cores)
    monkeypatch.setattr(limits, "machine_memory_bytes", lambda: memory_gb * 1024**3)
    monkeypatch.setattr(limits, "machine_gpus", lambda: gpus)


# --- measuring the machine ----------------------------------------------------------------


@pytest.mark.real_machine
def test_cores_never_zero():
    assert limits.machine_cores() >= 1


@pytest.mark.real_machine
def test_memory_never_zero():
    assert limits.machine_memory_bytes() >= 1


@pytest.mark.real_machine
def test_cgroup_quota_caps_the_affinity_count(monkeypatch):
    # The case this exists for: a container that can see the host's cores but may not use them.
    monkeypatch.setattr(os, "sched_getaffinity", lambda _pid: set(range(192)), raising=False)
    monkeypatch.setattr(limits, "_cgroup_cpu_quota", lambda: 8)
    assert limits.machine_cores() == 8


def test_cgroup_quota_rounds_up(monkeypatch, tmp_path):
    # 1.5 cores of quota is 2 threads' worth of parallelism, not 1.
    path = tmp_path / "cpu.max"
    path.write_text("150000 100000\n")
    monkeypatch.setattr("builtins.open", _only_for(str(path), "150000 100000\n"))
    assert limits._cgroup_cpu_quota() == 2


def _only_for(target, contents):
    """An ``open`` that answers cgroup reads from a fixture and passes everything else through."""
    import builtins

    real = builtins.open

    def fake(path, *args, **kwargs):
        if str(path) == "/sys/fs/cgroup/cpu.max":
            return real(target, *args, **kwargs)
        return real(path, *args, **kwargs)

    return fake


def test_unlimited_cgroup_reads_as_no_limit(tmp_path):
    path = tmp_path / "memory.max"
    path.write_text("max\n")
    assert limits._read_int(str(path)) is None


def test_missing_cgroup_file_reads_as_no_limit(tmp_path):
    assert limits._read_int(str(tmp_path / "absent")) is None


@pytest.mark.real_machine
def test_gpus_respect_cuda_visible_devices(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "2,5")
    assert limits.machine_gpus() == 2


def test_gpus_disabled_by_empty_cuda_visible_devices(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    assert limits.machine_gpus() == 0


def test_no_nvidia_smi_means_no_gpus(monkeypatch):
    # A CPU-only machine is a normal machine, not an error.
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setattr(limits.shutil, "which", lambda _name: None)
    assert limits.machine_gpus() == 0


@pytest.mark.real_machine
def test_the_gpu_probe_runs_once_per_process(monkeypatch):
    # The budget is resolved from more than one place, and the answer cannot change mid-run.
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    calls = []
    monkeypatch.setattr(limits.shutil, "which", lambda _name: calls.append(1) or "/bin/nvidia-smi")
    monkeypatch.setattr(
        limits.subprocess,
        "run",
        lambda *a, **k: type("Done", (), {"stdout": "GPU 0: NVIDIA\nGPU 1: NVIDIA\n"})(),
    )
    assert limits.machine_gpus() == 2
    assert limits.machine_gpus() == 2
    assert len(calls) == 1


# --- parsing the [limits] table -----------------------------------------------------------


def test_empty_table_gives_the_documented_defaults():
    config = limits.config_from_table({})
    assert config.max_thread_fraction == 0.5
    assert config.max_memory_fraction == 0.8
    assert config.max_disk_fraction == 0.7
    assert config.max_threads is None


def test_every_documented_key_parses():
    # Guards against LIMITS_KEYS and config_from_table drifting apart, which would advertise a
    # key in the config error message that setting it then rejects.
    table = {
        "max_threads": 8,
        "max_memory_gb": 16,
        "max_disk_fraction": 0.5,
        "max_gpus": 2,
        "disk_check_interval": 30,
        "max_thread_fraction": 0.25,
        "max_memory_fraction": 0.5,
    }
    assert set(table) == set(limits.LIMITS_KEYS)
    limits.config_from_table(table)


@pytest.mark.parametrize(
    "table",
    [
        {"max_threads": 0},
        {"max_threads": "lots"},
        {"max_thread_fraction": 0},
        {"max_thread_fraction": 1.5},
        {"max_memory_fraction": -0.1},
        {"max_memory_gb": 0},
        {"max_disk_fraction": 2},
        {"max_gpus": -1},
        {"disk_check_interval": 0},
    ],
)
def test_bad_values_are_rejected(table):
    with pytest.raises(limits.LimitsError):
        limits.config_from_table(table)


def test_zero_gpus_is_allowed():
    # "this run may not touch the GPUs" is a real request, unlike zero threads.
    assert limits.config_from_table({"max_gpus": 0}).max_gpus == 0


# --- resolving the budget -----------------------------------------------------------------


def test_fractions_are_applied_to_the_measured_machine(monkeypatch):
    machine(cores=64, memory_gb=100, monkeypatch=monkeypatch)
    budget = limits.resolve()
    assert budget.threads == 32
    assert budget.memory_gb == pytest.approx(80.0)


def test_absolutes_win_over_fractions(monkeypatch):
    machine(cores=64, memory_gb=100, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=10, max_memory_gb=20))
    assert budget.threads == 10
    assert budget.memory_gb == pytest.approx(20.0)


def test_absolutes_are_still_clamped_to_the_machine(monkeypatch):
    # Asking for more than exists is a mistake to correct, not a request to honour into an OOM.
    machine(cores=8, memory_gb=16, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=64, max_memory_gb=256))
    assert budget.threads == 8
    assert budget.memory_gb == pytest.approx(16.0)


def test_threads_never_resolve_below_one(monkeypatch):
    machine(cores=1, monkeypatch=monkeypatch)
    assert limits.resolve(limits.LimitsConfig(max_thread_fraction=0.1)).threads == 1


def test_gpus_are_capped_but_not_invented(monkeypatch):
    machine(gpus=8, monkeypatch=monkeypatch)
    assert limits.resolve(limits.LimitsConfig(max_gpus=2)).gpus == 2
    assert limits.resolve(limits.LimitsConfig(max_gpus=99)).gpus == 8
    assert limits.resolve().gpus == 8


def test_budget_records_what_it_was_measured_against(monkeypatch):
    machine(cores=192, memory_gb=728, gpus=8, monkeypatch=monkeypatch)
    budget = limits.resolve()
    assert budget.machine_cores == 192
    assert budget.machine_gpus == 8
    assert "96/192 threads" in budget.describe()
    assert "8/8 GPUs" in budget.describe()
    assert budget.as_dict()["machine_cores"] == 192


def test_describe_omits_gpus_on_a_cpu_only_machine(monkeypatch):
    machine(gpus=0, monkeypatch=monkeypatch)
    assert "GPU" not in limits.resolve().describe()


# --- dividing among concurrent workers ----------------------------------------------------


def test_per_worker_divides_rather_than_multiplying(monkeypatch):
    machine(cores=96, memory_gb=200, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=96, max_memory_gb=160))
    worker = budget.per_worker(8)
    assert worker.threads == 12
    assert worker.memory_gb == pytest.approx(20.0)
    # The whole point: eight workers add up to the budget, not eight times it.
    assert worker.threads * 8 <= budget.threads
    assert worker.memory_bytes * 8 <= budget.memory_bytes


def test_per_worker_floors_threads_at_one(monkeypatch):
    machine(cores=4, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=2))
    assert budget.per_worker(20).threads == 1


def test_per_worker_gpus_floor_at_zero(monkeypatch):
    # No such thing as a fractional device: some workers get none and fall back to CPU.
    machine(gpus=8, monkeypatch=monkeypatch)
    assert limits.resolve().per_worker(20).gpus == 0
    assert limits.resolve().per_worker(4).gpus == 2


def test_per_worker_keeps_the_machine_provenance(monkeypatch):
    machine(cores=96, monkeypatch=monkeypatch)
    assert limits.resolve().per_worker(4).machine_cores == 96


def test_per_worker_rejects_nonsense_concurrency(monkeypatch):
    machine(monkeypatch=monkeypatch)
    with pytest.raises(limits.LimitsError):
        limits.resolve().per_worker(0)


# --- command-line overrides ---------------------------------------------------------------


def test_override_beats_the_config_file(monkeypatch):
    machine(cores=64, monkeypatch=monkeypatch)
    limits.set_override("max_threads", 4)
    assert limits.resolve(limits.LimitsConfig(max_threads=32)).threads == 4


def test_a_fraction_override_displaces_a_configured_absolute(monkeypatch):
    # Without the exclusive-pair rule the configured absolute would win and the flag would
    # silently do nothing, which is exactly the failure this module exists to prevent.
    machine(cores=64, monkeypatch=monkeypatch)
    limits.set_override("max_thread_fraction", 0.25)
    assert limits.resolve(limits.LimitsConfig(max_threads=64)).threads == 16


def test_override_values_are_validated():
    with pytest.raises(limits.LimitsError):
        limits.set_override("max_threads", 0)
    with pytest.raises(limits.LimitsError):
        limits.set_override("max_disk_fraction", 3)


def test_unknown_override_is_rejected():
    with pytest.raises(limits.LimitsError):
        limits.set_override("max_pizza", 1)


def test_overrides_are_reportable(monkeypatch):
    limits.set_override("max_gpus", 0)
    assert limits.overrides() == {"max_gpus": 0}


# --- the disk ceiling ---------------------------------------------------------------------


def budget_with_disk(fraction: float, interval: int = 60) -> limits.Budget:
    return limits.Budget(
        threads=1,
        memory_bytes=1024**3,
        disk_fraction=fraction,
        gpus=0,
        disk_check_interval=interval,
        machine_cores=1,
        machine_memory_bytes=1024**3,
        machine_gpus=0,
    )


def fake_usage(monkeypatch, total: int, free: int):
    monkeypatch.setattr(
        limits.shutil,
        "disk_usage",
        lambda _path: type("Usage", (), {"total": total, "free": free, "used": total - free})(),
    )


def test_disk_state_measures_used_as_total_minus_free(tmp_path, monkeypatch):
    # Not usage.used: filesystems reserve space for root, so used + free < total, and what
    # matters is what a non-root writer can still consume.
    fake_usage(monkeypatch, total=1000, free=500)
    state = limits.disk_state(str(tmp_path), budget_with_disk(0.7))
    assert state is not None
    assert state.used_bytes == 500
    assert state.ceiling_bytes == 700
    assert state.headroom_bytes == 200
    assert not state.over


def test_disk_state_reports_over_the_ceiling(tmp_path, monkeypatch):
    fake_usage(monkeypatch, total=1000, free=100)
    state = limits.disk_state(str(tmp_path), budget_with_disk(0.7))
    assert state is not None
    assert state.over
    assert state.headroom_bytes == 0  # never negative
    assert "90% full" in state.describe()


def test_disk_state_walks_up_to_an_existing_parent(tmp_path, monkeypatch):
    # A step's target directory usually does not exist yet, but its filesystem does.
    fake_usage(monkeypatch, total=1000, free=900)
    state = limits.disk_state(str(tmp_path / "not" / "yet" / "made"), budget_with_disk(0.7))
    assert state is not None
    assert state.path == str(tmp_path)


def test_unmeasurable_filesystem_is_not_a_failure(tmp_path, monkeypatch):
    def boom(_path):
        raise OSError("no stat for you")

    monkeypatch.setattr(limits.shutil, "disk_usage", boom)
    assert limits.disk_state(str(tmp_path), budget_with_disk(0.7)) is None
    ok, message = limits.check_ceiling([str(tmp_path)], budget_with_disk(0.7))
    assert ok
    assert "could be measured" in message


def test_distinct_filesystems_collapses_one_device(tmp_path):
    (tmp_path / "out").mkdir()
    (tmp_path / "scratch").mkdir()
    paths = limits.distinct_filesystems([str(tmp_path / "out"), str(tmp_path / "scratch")])
    assert len(paths) == 1


def test_distinct_filesystems_drops_unresolvable_paths():
    assert limits.distinct_filesystems(["", None]) == []  # type: ignore[list-item]


def test_check_ceiling_fails_when_over(tmp_path, monkeypatch):
    fake_usage(monkeypatch, total=1000, free=10)
    ok, message = limits.check_ceiling([str(tmp_path)], budget_with_disk(0.7))
    assert not ok
    assert "full" in message


def test_monitor_samples_nothing_when_under_the_ceiling(tmp_path, monkeypatch):
    fake_usage(monkeypatch, total=1000, free=900)
    monitor = limits.DiskMonitor([str(tmp_path)], budget_with_disk(0.7))
    assert monitor.sample() is None


def test_monitor_sample_finds_the_offending_filesystem(tmp_path, monkeypatch):
    fake_usage(monkeypatch, total=1000, free=10)
    monitor = limits.DiskMonitor([str(tmp_path)], budget_with_disk(0.7))
    state = monitor.sample()
    assert state is not None and state.over


def test_monitor_starts_no_thread_without_paths():
    with limits.DiskMonitor([], budget_with_disk(0.7)) as monitor:
        assert monitor._thread is None


def test_monitor_interrupts_the_main_thread(tmp_path, monkeypatch):
    # The stop has to arrive as a KeyboardInterrupt in the main thread, because that is the only
    # thing that reaches code blocked in a subprocess wait - and it makes --resume work for free.
    fake_usage(monkeypatch, total=1000, free=10)
    interrupted = []
    monkeypatch.setattr(
        "_thread.interrupt_main", lambda *a: interrupted.append(True), raising=False
    )
    seen = []
    with limits.DiskMonitor(
        [str(tmp_path)], budget_with_disk(0.7, interval=1), on_stop=seen.append
    ) as monitor:
        assert monitor._thread is not None
        monitor._thread.join(timeout=5)
    assert interrupted == [True]
    assert seen and seen[0].over
    assert monitor.tripped is not None


def test_monitor_stops_even_if_reporting_raises(tmp_path, monkeypatch):
    fake_usage(monkeypatch, total=1000, free=10)
    interrupted = []
    monkeypatch.setattr(
        "_thread.interrupt_main", lambda *a: interrupted.append(True), raising=False
    )

    def boom(_state):
        raise RuntimeError("logging is broken")

    with limits.DiskMonitor(
        [str(tmp_path)], budget_with_disk(0.7, interval=1), on_stop=boom
    ) as monitor:
        monitor._thread.join(timeout=5)  # type: ignore[union-attr]
    assert interrupted == [True]


def test_disk_ceiling_reached_is_an_interrupt():
    # So every finally block and --resume path written for an interrupted run applies unchanged.
    assert issubclass(limits.DiskCeilingReached, KeyboardInterrupt)


# --- GPU arbitration ----------------------------------------------------------------------


def test_pool_hands_out_each_device_once():
    pool = limits.GpuPool(count=3)
    assert sorted([pool.acquire(), pool.acquire(), pool.acquire()]) == [0, 1, 2]
    assert pool.acquire() is None  # returns rather than blocking


def test_released_devices_come_back():
    pool = limits.GpuPool(count=1)
    device = pool.acquire()
    assert pool.acquire() is None
    pool.release(device)
    assert pool.acquire() == 0


def test_releasing_nothing_is_harmless():
    # So a caller that fell back to CPU needs no branch around release().
    pool = limits.GpuPool(count=1)
    pool.release(None)
    assert len(pool) == 1


def test_double_release_does_not_duplicate_a_device():
    pool = limits.GpuPool(count=1)
    pool.release(0)
    pool.release(0)
    assert len(pool) == 1


def test_pool_from_budget_honours_the_cap(monkeypatch):
    machine(gpus=8, monkeypatch=monkeypatch)
    pool = limits.GpuPool.from_budget(limits.resolve(limits.LimitsConfig(max_gpus=2)))
    assert len(pool) == 2


def test_empty_pool_on_a_cpu_only_machine():
    assert limits.GpuPool(count=0).acquire() is None


def test_gpu_env_pins_one_device():
    env = limits.gpu_env(3, base={"PATH": "/bin"})
    assert env["CUDA_VISIBLE_DEVICES"] == "3"
    assert env["PATH"] == "/bin"


def test_gpu_env_without_a_device_disables_cuda():
    # An empty value is how CUDA is told there are none, which turns a framework's device probe
    # into a CPU fallback instead of a crash.
    assert limits.gpu_env(None, base={})["CUDA_VISIBLE_DEVICES"] == ""


# --- wiring: the config file ---------------------------------------------------------------


def write_config(tmp_path, body: str, monkeypatch) -> str:
    path = tmp_path / "metawrap2.toml"
    path.write_text(body)
    monkeypatch.setenv("METAWRAP2_CONFIG", str(path))
    return str(path)


def test_limits_table_reaches_settings(tmp_path, monkeypatch):
    from metawrap2 import config

    write_config(tmp_path, "[limits]\nmax_threads = 8\nmax_disk_fraction = 0.5\n", monkeypatch)
    settings = config.load_settings()
    assert settings.limits.max_threads == 8
    assert settings.limits.max_disk_fraction == 0.5


def test_a_config_without_limits_gets_the_defaults(tmp_path, monkeypatch):
    from metawrap2 import config

    write_config(tmp_path, "[settings]\nthreads = 4\n", monkeypatch)
    assert config.load_settings().limits.max_thread_fraction == 0.5


def test_limits_is_a_known_table():
    from metawrap2 import config

    problems, _notes = config.check_config({"limits": {"max_threads": 4}}, "cfg.toml")
    assert problems == []


def test_a_typo_in_limits_is_a_problem_with_a_suggestion():
    # Same reasoning as [settings]: a typo here means the limit silently does not apply, which is
    # indistinguishable from never having written it.
    from metawrap2 import config

    problems, _notes = config.check_config({"limits": {"max_thread": 4}}, "cfg.toml")
    assert len(problems) == 1
    assert "max_thread" in problems[0] and "did you mean 'max_threads'" in problems[0]


def test_a_bad_limits_value_is_a_problem_too():
    from metawrap2 import config

    problems, _notes = config.check_config({"limits": {"max_disk_fraction": 5}}, "cfg.toml")
    assert len(problems) == 1
    assert "max_disk_fraction" in problems[0] and "cfg.toml" in problems[0]


def test_loading_a_bad_limits_value_raises(tmp_path, monkeypatch):
    from metawrap2 import config

    write_config(tmp_path, "[limits]\nmax_memory_fraction = 4\n", monkeypatch)
    with pytest.raises(config.ConfigError):
        config.load_settings()


def test_a_bad_limits_value_falls_back_when_not_strict(tmp_path, monkeypatch):
    from metawrap2 import config

    write_config(tmp_path, "[limits]\nmax_gpus = -3\n", monkeypatch)
    settings = config.load_settings(strict=False)
    assert settings.limits.max_gpus is None  # the default, not the rejected value


# --- wiring: the thread ceiling ------------------------------------------------------------


class Args:
    def __init__(self, threads=None):
        self.threads = threads


def test_resolve_threads_caps_to_the_budget(monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    settings = Settings(limits=limits.LimitsConfig(max_threads=8))
    assert _common.resolve_threads(Args(threads=32), settings) == 8


def test_resolve_threads_honours_a_request_inside_the_budget(monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    assert _common.resolve_threads(Args(threads=4), Settings()) == 4


def test_resolve_threads_explains_the_cap(monkeypatch, capsys):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    _common.resolve_threads(Args(threads=64), Settings())
    out = capsys.readouterr()
    assert "max_thread_fraction" in (out.out + out.err)


def test_capping_a_config_default_is_silent(monkeypatch, capsys):
    # The user did not type -t here, and the budget line start_run prints already says the number.
    # Warning every run about a default the same user configured trains people to ignore warnings.
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    assert _common.resolve_threads(Args(threads=None), Settings(threads=64)) == 32
    out = capsys.readouterr()
    assert "Asked for" not in (out.out + out.err)


def test_resolve_threads_writes_the_capped_value_back_to_args(monkeypatch):
    # Modules read args.threads directly after this, so the cap has to land there too.
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    args = Args(threads=999)
    _common.resolve_threads(args, Settings())
    assert args.threads == 32


# --- wiring: the global flags --------------------------------------------------------------


def test_global_flags_are_stripped_and_applied():
    from metawrap2 import cli

    kept = cli._apply_global_flags(["--max-threads", "4", "binning", "-o", "out"])
    assert kept == ["binning", "-o", "out"]
    assert limits.overrides() == {"max_threads": 4}


def test_equals_form_works_too():
    from metawrap2 import cli

    assert cli._apply_global_flags(["--max-gpus=0", "binning"]) == ["binning"]
    assert limits.overrides() == {"max_gpus": 0}


def test_a_flag_without_its_value_is_reported():
    from metawrap2 import cli

    with pytest.raises(limits.LimitsError):
        cli._apply_global_flags(["--max-threads"])


# --- wiring: the memory ceiling ------------------------------------------------------------


class MemArgs:
    def __init__(self, memory=None):
        self.memory = memory


def test_unset_memory_resolves_to_the_budget(monkeypatch):
    # The point of having a budget: the old fixed constants (24 GB for an assembly, 4 GB for
    # binning) left most of a large machine unused.
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(memory_gb=100, monkeypatch=monkeypatch)
    assert _common.resolve_memory_gb(MemArgs(), Settings()) == 80


def test_explicit_memory_is_honoured_inside_the_budget(monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(memory_gb=100, monkeypatch=monkeypatch)
    assert _common.resolve_memory_gb(MemArgs(memory=16), Settings()) == 16


def test_explicit_memory_above_the_budget_is_capped(monkeypatch, capsys):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(memory_gb=100, monkeypatch=monkeypatch)
    args = MemArgs(memory=500)
    assert _common.resolve_memory_gb(args, Settings()) == 80
    assert args.memory == 80  # written back, since modules read it directly
    assert "max_memory_fraction" in "".join(capsys.readouterr())


def test_zero_memory_is_refused(monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(monkeypatch=monkeypatch)
    with pytest.raises(SystemExit):
        _common.resolve_memory_gb(MemArgs(memory=0), Settings())


def test_memory_never_resolves_below_one_gb(monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(memory_gb=1, monkeypatch=monkeypatch)
    assert _common.resolve_memory_gb(MemArgs(), Settings()) >= 1


@pytest.mark.parametrize("module", ["assembly", "binning", "bin_refinement", "reassemble_bins"])
def test_every_module_with_memory_defers_to_the_budget(module, capsys):
    # Each of these used to advertise its own fixed number. They must all now say the same thing,
    # or a user reading --help still believes the old constant applies.
    import importlib

    mod = importlib.import_module("metawrap2.modules.%s" % module)
    with pytest.raises(SystemExit):
        mod._parse_args(["--help"])
    # argparse wraps help to the terminal width, so match against the unwrapped text. The
    # literal percent sign is the other half of the check: argparse %-formats help strings,
    # and an unescaped one raises instead of printing this page at all.
    rendered = " ".join(capsys.readouterr().out.split())
    assert "default: the [limits] budget, 80% of this machine's RAM" in rendered


# --- wiring: the disk ceiling --------------------------------------------------------------


def test_a_full_filesystem_stops_the_step_before_it_starts(tmp_path, monkeypatch):
    # Unlike the size estimate, which only warns, this one is measured - so there is nothing to be
    # wrong about and no reason to start a step the monitor would stop minutes later.
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    fake_usage(monkeypatch, total=1000, free=10)
    with pytest.raises(SystemExit):
        _common.check_disk_space("assembly", Settings(), [], str(tmp_path))


def test_the_ceiling_check_respects_skip_space_check(tmp_path, monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    fake_usage(monkeypatch, total=1000, free=10)
    _common.check_disk_space("assembly", Settings(skip_space_check=True), [], str(tmp_path))


def test_the_monitor_is_started_and_stopped_with_the_run(tmp_path, monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    started, stopped = [], []
    monkeypatch.setattr(_common._limits.DiskMonitor, "start", lambda self: started.append(self))
    monkeypatch.setattr(_common._limits.DiskMonitor, "stop", lambda self: stopped.append(self))
    machine(monkeypatch=monkeypatch)
    budget = limits.resolve()
    _common._start_disk_monitor("assembly", Settings(), str(tmp_path), budget)
    assert len(started) == 1 and _common._active_disk_monitor is not None
    monitor = _common._active_disk_monitor
    # finish_run is what releases it, so simulate the part of it that matters here.
    monitor.stop()
    _common._active_disk_monitor = None
    assert stopped == [monitor]


def test_the_monitor_is_not_started_when_space_checks_are_off(tmp_path, monkeypatch):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(monkeypatch=monkeypatch)
    _common._active_disk_monitor = None
    _common._start_disk_monitor(
        "assembly", Settings(skip_space_check=True), str(tmp_path), limits.resolve()
    )
    assert _common._active_disk_monitor is None


def test_overrides_are_described_as_flags():
    from metawrap2.modules._common import _describe_overrides

    described = _describe_overrides({"max_threads": 8, "max_gpus": 0})
    assert described == "--max-gpus 0, --max-threads 8"


def test_a_bad_budget_flag_exits_cleanly(capsys):
    # A user error deserves a message, not a traceback - and it must be caught before a module
    # starts, not hours in.
    from metawrap2 import cli

    assert cli.main(["--max-threads", "0", "binning"]) == 2
    assert "max_threads" in capsys.readouterr().err


def test_a_monitor_stopped_mid_sample_does_not_interrupt_the_next_step(monkeypatch):
    """stop() has to win against a watcher already inside sample().

    sample() is two statvfs calls and an os.stat, which on a stalled mount takes seconds. A step
    that finished in that window has already called stop(), so interrupting afterwards delivers a
    KeyboardInterrupt into whatever the process does next - under `run` or `cohort`, the following
    step, which is then aborted and recorded as interrupted while the step that was being watched
    is reported as a disk stop it never hit.
    """
    import threading

    interrupts = []
    monkeypatch.setattr("_thread.interrupt_main", lambda: interrupts.append(1))

    budget = limits.Budget(
        threads=1,
        memory_bytes=1,
        disk_fraction=0.7,
        gpus=0,
        disk_check_interval=0.01,
        machine_cores=1,
        machine_memory_bytes=1,
        machine_gpus=0,
    )
    monitor = limits.DiskMonitor(["/"], budget)
    sampling = threading.Event()
    release = threading.Event()

    def slow_sample():
        """Over the ceiling, but slow enough that the step finishes first."""
        sampling.set()
        release.wait(5)
        return limits.DiskState(path="/", total_bytes=100, used_bytes=99, ceiling_bytes=70)

    monkeypatch.setattr(monitor, "sample", slow_sample)
    monitor.start()
    assert sampling.wait(5), "the watcher never reached sample()"

    # The step finishing, while sample() is still in flight.
    stopped = threading.Thread(target=monitor.stop)
    stopped.start()
    release.set()
    stopped.join(10)

    assert not interrupts, "interrupted the main thread after the step it was watching had finished"
    assert monitor.tripped is None
