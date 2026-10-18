"""What built this catalogue: the strategy `run` and `cohort` record in the study manifest.

A co-assembly catalogue and a per-sample one are different objects - different genome counts,
different redundancy, different per-genome coverage - and nothing in a finished study directory used
to say which you were holding. These tests pin the record itself, because the whole value of it is
that it is there to be read later by someone who did not run the study.
"""

from __future__ import annotations

import json
import os

from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.commands import status_cmd
from metawrap2.manifest import MANIFEST_NAME

SHEET = """
[settings]
coassemble = %s

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq.gz"
r2 = "raw/B_2.fastq.gz"
"""


def _sheet(tmp_path, coassemble: bool) -> str:
    path = tmp_path / "samples.toml"
    path.write_text(SHEET % ("true" if coassemble else "false"))
    return str(path)


def _recorded(tmp_path, monkeypatch, coassemble: bool, command=driver.main, *extra):
    """Plan a study, run nothing, and return the strategy block the driver recorded."""
    monkeypatch.setattr(driver, "execute", lambda *args, **kwargs: 0)
    out = tmp_path / "study"
    assert command([_sheet(tmp_path, coassemble), "-o", str(out), *extra]) == 0
    manifest = out / driver._STATE_DIR / MANIFEST_NAME
    return json.loads(manifest.read_text())[driver.STRATEGY_KEY], str(out)


def test_a_co_assembly_run_records_co_assembly(tmp_path, monkeypatch, capsys):
    """`metawrap2 run` says so in the manifest, and `metawrap2 status --json` reads it back.

    Regresses the case this exists for: two studies in two directories, one of each kind, and
    nothing in either to tell them apart six months later.
    """
    record, out = _recorded(tmp_path, monkeypatch, coassemble=True)
    # Dereplication is not merely off here, it has nothing to do: one assembly, one bin set.
    assert record == {"strategy": driver.CO_ASSEMBLY, "samples": 2, "dereplicated": False}

    capsys.readouterr()
    assert status_cmd.main([out, "--json"]) == 0
    assert driver.CO_ASSEMBLY in capsys.readouterr().out


def test_a_per_sample_cohort_records_per_sample(tmp_path, monkeypatch):
    """`metawrap2 cohort` records the other strategy, and that dereplication collapsed the bins.

    Regresses the alias drifting from the driver: `cohort` is `run --no-coassembly`, so the record
    has to come out of the study the planner settled on rather than out of which command was typed.
    """
    record, _out = _recorded(tmp_path, monkeypatch, coassemble=False, command=cohort_cmd.main)
    assert record == {"strategy": driver.PER_SAMPLE, "samples": 2, "dereplicated": True}


def test_the_strategy_is_recorded_even_when_dereplication_was_skipped(tmp_path, monkeypatch):
    """A cohort with `dereplicate` skipped still records what it was - and that it did not collapse.

    Regresses the worst comparison this feature is meant to prevent: without dereplication the
    per-sample bin sets are never collapsed, so the catalogue is redundant by construction and its
    genome count is not comparable with a co-assembly's or with a dereplicated cohort's.
    """
    record, out = _recorded(tmp_path, monkeypatch, False, driver.main, "--skip", "dereplicate")
    assert record == {"strategy": driver.PER_SAMPLE, "samples": 2, "dereplicated": False}
    # And recorded before any step ran, so an interrupted study says what it was too.
    with open(os.path.join(out, driver._STATE_DIR, MANIFEST_NAME)) as fh:
        assert json.load(fh)["steps"] == {}
