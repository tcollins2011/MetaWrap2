"""The single place that owns *running an external command* in MetaWrap2.

Everything about executing a tool goes through :class:`CommandRunner`:

* wrapping the command in a module's conda environment via **mamba** (``mamba run -n env``);
* honoring ``--dry-run`` (record and print the command, don't execute);
* streaming the child's stdout and stderr live to the user's screen **and** teeing them to
  ``run.stdout`` / ``run.stderr`` in the run's output directory (two reader threads, so the
  screen output is never blocked or interleaved into one stream);
* redirecting a command's stdout to a data file when that stdout *is* the output
  (``stdout_path=`` - e.g. ``bwa mem > x.sam``), keeping stderr separate;
* recording each command into the active provenance recorder;
* checking the exit code and raising a clear :class:`ToolError` on failure;
* cleaning up the child process group on Ctrl-C.

Modules never call subprocess directly - they call :func:`run` (re-exported via
``modules._common``). There is one module-level :data:`runner` instance; the CLI and
``start_run`` configure it (dry-run, env manager, per-run log files, recorder).
"""

from __future__ import annotations

import os
import shlex
import shutil
import signal
import subprocess
import sys
import threading
from dataclasses import dataclass
from typing import Dict, List, Optional, Protocol, Sequence, TextIO, Union

from .logfilter import LogFilter
from .logging import log_tool_output


class Recorder(Protocol):
    """What the runner needs from a provenance recorder.

    A Protocol rather than an import of RunRecorder: provenance imports nothing from here, and
    keeping it that way avoids a circular import while still type-checking the one call made.
    """

    def record_command(self, cmd: str) -> None: ...


__all__ = [
    "CommandRunner",
    "EnvManagerMissing",
    "ToolError",
    "run",
    "runner",
    "set_dry_run",
    "set_env_manager",
    "set_force",
    "set_recorder",
    "set_resume",
    "set_run_logs",
    "set_skip_space_check",
    "set_skip_validation",
    "set_verbose_logs",
    "tool_path_in_env",
]


@dataclass
class ToolError(RuntimeError):
    """Raised when an external tool exits non-zero."""

    tool: str
    returncode: int
    cmd: Sequence[str]
    log_tail: str = ""
    hint: str = ""

    def __str__(self) -> str:
        msg = ["%s failed (exit code %d)." % (self.tool, self.returncode)]
        msg.append("Command: " + " ".join(str(c) for c in self.cmd))
        if self.log_tail:
            msg.append("Last output:\n" + self.log_tail)
        if self.hint:
            msg.append("Hint: " + self.hint)
        return "\n".join(msg)


class EnvManagerMissing(RuntimeError):
    """Raised when a command needs a conda environment but the env manager is not on PATH.

    Without this check the failure surfaces from deep inside :func:`subprocess.Popen` as a bare
    ``FileNotFoundError: 'mamba'``, which reads like the *tool* is missing rather than the thing
    that was supposed to find it. It is raised when the argv is built rather than when it runs, so
    a long pipeline stops on its first step instead of recording a directory of failed ones.
    """

    def __init__(self, env_manager: str, env: str) -> None:
        super().__init__(
            "%s is not on PATH, so tools cannot be run inside their conda environment (%s).\n"
            "Either add it - a conda/mamba install usually needs "
            'PATH="$HOME/miniforge3/bin:$PATH" exported, which a non-login shell and most '
            "batch schedulers do not do for you - or set use_conda_envs = false under [settings] "
            "in ~/.metawrap2/config.toml to run whatever is already on PATH instead."
            % (env_manager, env)
        )
        self.env_manager = env_manager
        self.env = env


_RUN_FLAGS_CACHE: dict = {}


def _run_flags(env_manager: str) -> List[str]:
    """Extra flags for ``<env_manager> run`` so the child's output streams through live.

    ``--no-capture-output`` is a *conda* flag. mamba 2.x does not accept it and mis-parses it
    into the generated wrapper script, so every command fails with
    ``exec: --: invalid option`` - which looks like the tool is missing rather than like a
    flag problem. mamba 2.x streams output by default anyway, so when the flag is not
    supported we simply omit it. The answer is cached: this shells out once per process.
    """
    if env_manager in _RUN_FLAGS_CACHE:
        return _RUN_FLAGS_CACHE[env_manager]
    flags: List[str] = []
    try:
        helptext = (
            subprocess.run(
                [env_manager, "run", "--help"],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=60,
                check=False,
            ).stdout
            or ""
        )
        if "--no-capture-output" in helptext:
            flags = ["--no-capture-output"]
    except (OSError, subprocess.SubprocessError):
        flags = []
    _RUN_FLAGS_CACHE[env_manager] = flags
    return flags


def _pump(
    pipe,
    screen: Optional[TextIO],
    logfile,
    tail: Optional[list],
    log_filter: Optional[LogFilter] = None,
    stream_name: str = "stdout",
) -> None:
    """Read *pipe* line by line, echoing to *screen* and appending to *logfile*/*tail*.

    *log_filter* collapses progress-bar redraws and caps runaway repeated messages (see
    :mod:`metawrap2.logfilter`). The error *tail* is always fed the unfiltered line, so the
    message in a ToolError still shows exactly what the tool last said.
    """
    try:
        for line in iter(pipe.readline, ""):
            if tail is not None:
                tail.append(line)
                if len(tail) > 25:
                    tail.pop(0)
            shown = log_filter.feed(line) if log_filter is not None else line
            if shown is None:
                continue
            if screen is not None:
                screen.write(shown)
                screen.flush()
            if logfile is not None:
                logfile.write(shown)
                logfile.flush()
            # Also record it in the combined run log, so metawrap2.log has the whole story of
            # the run in order rather than only what MetaWrap2 itself said.
            log_tool_output(shown, stream_name)
    finally:
        pipe.close()


@dataclass
class CommandRunner:
    """Runs external commands with env-wrapping, logging, provenance, and dry-run."""

    env_manager: str = "mamba"  # program used for `... run -n env`
    dry_run: bool = False
    force: bool = False  # overwrite an existing MetaWrap2 output dir
    resume: bool = False  # reuse an existing output dir, skip finished steps
    recorder: Optional[Recorder] = None  # provenance RunRecorder (or None)
    run_stdout_path: Optional[str] = None
    run_stderr_path: Optional[str] = None
    verbose_logs: bool = False  # True disables log noise filtering
    skip_space_check: bool = False  # True disables the preflight disk-space estimate
    skip_validation: bool = False  # True bypasses input/intermediate file validation

    # -- configuration --------------------------------------------------------------------

    def configure(self, **kw) -> None:
        for key, value in kw.items():
            if not hasattr(self, key):
                raise AttributeError("CommandRunner has no option %r" % key)
            setattr(self, key, value)

    def env_prefix(self, env: Optional[str]) -> List[str]:
        """argv prefix that runs a command inside conda env *env* (or [] if None)."""
        if not env:
            return []
        if not self.dry_run and shutil.which(self.env_manager) is None:
            # Not under a dry run: a dry run is for seeing the plan on a machine that may not have
            # anything installed yet, and printing the command it would run costs nothing there.
            raise EnvManagerMissing(self.env_manager, env)
        return [self.env_manager, "run"] + _run_flags(self.env_manager) + ["-n", env]

    # -- execution ------------------------------------------------------------------------

    def run(
        self,
        cmd: Union[str, Sequence[str]],  # noqa: UP007 - py3.8 target
        *,
        env: Optional[str] = None,
        tool: Optional[str] = None,
        stdout_path: Optional[str] = None,
        log_path: Optional[str] = None,  # backwards-compatible alias for stdout_path
        hint: str = "",
        check: bool = True,
        cwd: Optional[str] = None,
        extra_env: Optional[Dict[str, str]] = None,
    ) -> int:
        argv = shlex.split(cmd) if isinstance(cmd, str) else list(cmd)
        tool = tool or (argv[0] if argv else "external tool")
        stdout_path = stdout_path or log_path
        full = self.env_prefix(env) + argv

        recorded = shlex.join(full)
        if cwd:
            recorded = "cd %s && %s" % (shlex.quote(cwd), recorded)
        if stdout_path:
            recorded += " > " + shlex.quote(stdout_path)
        if self.recorder is not None:
            self.recorder.record_command(recorded)

        if self.dry_run:
            sys.stderr.write("[dry-run] %s\n" % recorded)
            return 0

        return self._execute(
            full,
            tool=tool,
            stdout_path=stdout_path,
            hint=hint,
            check=check,
            cwd=cwd,
            recorded=recorded,
            isolate_user_site=bool(env),
            extra_env=extra_env,
        )

    def _open_run_log(self, path: Optional[str], header: str):
        if not path:
            return None
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fh = open(path, "a")  # noqa: SIM115 - closed in _execute's finally, spans the child
        fh.write(header)
        return fh

    def _child_env(
        self, isolate_user_site: bool, extra: Optional[Dict[str, str]] = None
    ) -> Optional[Dict[str, str]]:
        """The environment for the child process, or None to inherit this one unchanged.

        A tool running inside a conda env gets ``PYTHONNOUSERSITE=1``, because ``~/.local`` comes
        *before* the env on ``sys.path`` and silently overrides whatever the env pinned. That is not
        hypothetical: AMBER's env pins matplotlib below 3.9 for an API it uses, and a matplotlib 3.11
        in the user's own site-packages was imported instead and crashed it. Only for env-wrapped
        commands - without an env there is no pin to protect, and a tool the user pip-installed into
        ``~/.local`` is then exactly the tool they meant to run.

        *extra* is for a tool that needs a variable set for it specifically; see THREAD_LIMITS.
        """
        if not isolate_user_site and not extra:
            return None
        child = dict(os.environ)
        if isolate_user_site:
            child["PYTHONNOUSERSITE"] = "1"
        child.update(extra or {})
        return child

    def _execute(
        self,
        full,
        *,
        tool,
        stdout_path,
        hint,
        check,
        cwd,
        recorded,
        isolate_user_site=False,
        extra_env=None,
    ) -> int:
        header = "\n$ %s\n" % recorded
        run_out = self._open_run_log(self.run_stdout_path, header)
        run_err = self._open_run_log(self.run_stderr_path, header)
        # A command whose stdout *is* the result goes straight to the file, unfiltered.

        # closed in the finally below.
        data_out = open(stdout_path, "w") if stdout_path else None  # noqa: SIM115
        out_filter = LogFilter(enabled=not self.verbose_logs)
        err_filter = LogFilter(enabled=not self.verbose_logs)
        tail: list = []
        threads: List[threading.Thread] = []
        proc = None
        try:
            proc = subprocess.Popen(
                list(full),
                stdout=(data_out if data_out is not None else subprocess.PIPE),
                stderr=subprocess.PIPE,
                cwd=cwd,
                text=True,
                bufsize=1,
                env=self._child_env(isolate_user_site, extra_env),
                start_new_session=True,  # own process group, so Ctrl-C can kill the tree
            )
            # stderr always streams to the screen + run.stderr (+ tail for error messages)
            threads.append(
                threading.Thread(
                    target=_pump,
                    args=(proc.stderr, sys.stderr, run_err, tail, err_filter, "stderr"),
                    daemon=True,
                )
            )
            # stdout: to the screen + run.stdout, unless it was redirected to a data file
            if data_out is None:
                threads.append(
                    threading.Thread(
                        target=_pump,
                        args=(proc.stdout, sys.stdout, run_out, None, out_filter, "stdout"),
                        daemon=True,
                    )
                )
            for t in threads:
                t.start()
            rc = proc.wait()
            for t in threads:
                t.join()
        except KeyboardInterrupt:
            # Everything from Popen onwards, not just the wait: an interrupt arriving while the
            # output pumps were being started left the tool running with all 96 threads it had
            # been granted, orphaned by the run that was exiting - and the next thing the user
            # started then competed with it for the machine.
            if proc is not None:
                self._terminate(proc)
            raise
        finally:
            # Report what was collapsed, so a quiet log never hides that it was filtered.
            for fh, filt in ((run_out, out_filter), (run_err, err_filter)):
                if fh is None:
                    continue
                for line in filt.summary():
                    fh.write(line)
            for line in err_filter.summary():
                sys.stderr.write(line)
            for fh in (run_out, run_err, data_out):
                if fh is not None:
                    fh.close()

        if check and rc != 0:
            raise ToolError(
                tool=tool, returncode=rc, cmd=full, log_tail="".join(tail).rstrip(), hint=hint
            )
        return rc

    @staticmethod
    def _terminate(proc) -> None:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            proc.terminate()


# One shared instance the whole package uses; configured by the CLI and start_run().
runner = CommandRunner()


#: Environment that confines a tool to the parallelism it was *asked* for.
#:
#: A numerical library picks its own thread count from the core count unless told otherwise, so a tool
#: that forks N workers and does linear algebra in each gets N x cores threads fighting over cores.
#:
#: Measured, not theorised, and measured *after* the fix rather than assumed from it: capping the
#: libraries took CONCOCT binning one CAMI2 sample from 29.9 CPU-hours to 23.8, a 1.3x saving, for a
#: bin set identical down to the count. Small but free, so it stays.
#:
#: Worth stating what it does not explain. The same CONCOCT version, with byte-identical arguments and
#: input, costs 5.6 CPU-hours under a pipeline whose environment carries MKL and Python 2.7 instead of
#: a pthreads OpenBLAS and 3.12. The cap closes a fifth of that; the rest is somewhere in the compiled
#: numerical stack and is not accounted for here. Do not read this constant as having fixed it.
#:
#: Applies to any tool that forks *and* does array maths. Passed explicitly rather than set globally,
#: because a tool that is internally threaded and does not fork wants the opposite.
THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def run(cmd, **kw) -> int:
    return runner.run(cmd, **kw)


def set_recorder(recorder) -> None:
    runner.configure(recorder=recorder)


def set_dry_run(value: bool) -> None:
    runner.configure(dry_run=bool(value))


def set_force(value: bool) -> None:
    runner.configure(force=bool(value))


def set_resume(value: bool) -> None:
    runner.configure(resume=bool(value))


def set_run_logs(stdout_path: Optional[str], stderr_path: Optional[str]) -> None:
    runner.configure(run_stdout_path=stdout_path, run_stderr_path=stderr_path)


def set_env_manager(name: str) -> None:
    runner.configure(env_manager=name)


def set_verbose_logs(value: bool) -> None:
    runner.configure(verbose_logs=bool(value))


def set_skip_space_check(value: bool) -> None:
    runner.configure(skip_space_check=bool(value))


def set_skip_validation(value: bool) -> None:
    runner.configure(skip_validation=bool(value))


_TOOL_PATH_CACHE: dict = {}


def tool_path_in_env(env: Optional[str], tool: str) -> str:
    """Absolute path to *tool* inside conda env *env* (or just *tool* if not found).

    MetaWrap2's own Python helpers run in the host interpreter, but a few of them shell out
    to an external tool that only exists inside a module's conda env (blobology's
    gc_cov_annotate needs samtools). Resolving the absolute path here lets the helper call it
    directly instead of relying on a PATH that will not contain it. Cached per process.
    """
    if not env:
        return tool
    key = (env, tool)
    if key in _TOOL_PATH_CACHE:
        return _TOOL_PATH_CACHE[key]
    resolved = tool
    try:
        argv = (
            [runner.env_manager, "run"]
            + _run_flags(runner.env_manager)
            + ["-n", env, "bash", "-c", "command -v %s" % shlex.quote(tool)]
        )
        out = subprocess.run(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=120,
            check=False,
        )
        candidate = (out.stdout or "").strip().splitlines()
        if out.returncode == 0 and candidate and os.path.isabs(candidate[-1]):
            resolved = candidate[-1]
    except (OSError, subprocess.SubprocessError):
        pass
    _TOOL_PATH_CACHE[key] = resolved
    return resolved


def capture_in_env(env: Optional[str], argv: Sequence[str], timeout: int = 120) -> str:
    """Run *argv* inside conda env *env* and return its stdout, or ``""`` if anything went wrong.

    For asking an env a question rather than making it do work - "can torch see a GPU in here?" -
    where the answer decides a flag and no answer is itself an answer. Deliberately silent: every
    caller has a usable fallback, and a traceback from a probe would look like a failure of the thing
    being probed for.
    """
    try:
        if env:
            full = (
                [runner.env_manager, "run"]
                + _run_flags(runner.env_manager)
                + ["-n", env]
                + list(argv)
            )
        else:
            full = list(argv)
        out = subprocess.run(
            full,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=timeout,
            check=False,
        )
        return (out.stdout or "") if out.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        return ""
