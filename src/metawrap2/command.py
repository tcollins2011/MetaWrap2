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
import signal
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, TextIO, Union

__all__ = ["CommandRunner", "ToolError", "runner", "run", "set_recorder", "set_dry_run",
           "set_force", "set_resume", "set_run_logs", "set_env_manager"]


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


def _pump(pipe, screen: Optional[TextIO], logfile, tail: Optional[list]) -> None:
    """Read *pipe* line by line, echoing to *screen* and appending to *logfile*/*tail*."""
    try:
        for line in iter(pipe.readline, ""):
            if screen is not None:
                screen.write(line)
                screen.flush()
            if logfile is not None:
                logfile.write(line)
                logfile.flush()
            if tail is not None:
                tail.append(line)
                if len(tail) > 25:
                    tail.pop(0)
    finally:
        pipe.close()


@dataclass
class CommandRunner:
    """Runs external commands with env-wrapping, logging, provenance, and dry-run."""

    env_manager: str = "mamba"          # program used for `... run -n env`
    dry_run: bool = False
    force: bool = False                 # overwrite an existing MetaWrap2 output dir
    resume: bool = False                # reuse an existing output dir, skip finished steps
    recorder: object = None             # provenance RunRecorder (or None)
    run_stdout_path: Optional[str] = None
    run_stderr_path: Optional[str] = None

    # -- configuration --------------------------------------------------------------------

    def configure(self, **kw) -> None:
        for key, value in kw.items():
            if not hasattr(self, key):
                raise AttributeError("CommandRunner has no option %r" % key)
            setattr(self, key, value)

    def env_prefix(self, env: Optional[str]) -> List[str]:
        """argv prefix that runs a command inside conda env *env* using mamba (or [] if None)."""
        return [self.env_manager, "run", "--no-capture-output", "-n", env] if env else []

    # -- execution ------------------------------------------------------------------------

    def run(
        self,
        cmd: Union[str, Sequence[str]],
        *,
        env: Optional[str] = None,
        tool: Optional[str] = None,
        stdout_path: Optional[str] = None,
        log_path: Optional[str] = None,   # backwards-compatible alias for stdout_path
        hint: str = "",
        check: bool = True,
        cwd: Optional[str] = None,
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

        return self._execute(full, tool=tool, stdout_path=stdout_path, hint=hint,
                             check=check, cwd=cwd, recorded=recorded)

    def _open_run_log(self, path: Optional[str], header: str):
        if not path:
            return None
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fh = open(path, "a")
        fh.write(header)
        return fh

    def _execute(self, full, *, tool, stdout_path, hint, check, cwd, recorded) -> int:
        header = "\n$ %s\n" % recorded
        run_out = self._open_run_log(self.run_stdout_path, header)
        run_err = self._open_run_log(self.run_stderr_path, header)
        data_out = open(stdout_path, "w") if stdout_path else None
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
                start_new_session=True,  # own process group, so Ctrl-C can kill the tree
            )
            # stderr always streams to the screen + run.stderr (+ tail for error messages)
            threads.append(threading.Thread(
                target=_pump, args=(proc.stderr, sys.stderr, run_err, tail), daemon=True))
            # stdout: to the screen + run.stdout, unless it was redirected to a data file
            if data_out is None:
                threads.append(threading.Thread(
                    target=_pump, args=(proc.stdout, sys.stdout, run_out, None), daemon=True))
            for t in threads:
                t.start()
            try:
                rc = proc.wait()
            except KeyboardInterrupt:
                self._terminate(proc)
                raise
            for t in threads:
                t.join()
        finally:
            for fh in (run_out, run_err, data_out):
                if fh is not None:
                    fh.close()

        if check and rc != 0:
            raise ToolError(tool=tool, returncode=rc, cmd=full,
                            log_tail="".join(tail).rstrip(), hint=hint)
        return rc

    @staticmethod
    def _terminate(proc) -> None:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            proc.terminate()


# One shared instance the whole package uses; configured by the CLI and start_run().
runner = CommandRunner()


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
