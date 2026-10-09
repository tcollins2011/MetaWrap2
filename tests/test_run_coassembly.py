def test_a_coassembly_without_read_qc_still_gets_its_concatenated_reads(tmp_path, monkeypatch):
    """ALL_READS was only built after a read_qc step, so a study that skips QC had no inputs.

    Every arm of a benchmark that holds read QC constant skips it, which is exactly when co-assembly
    needs the raw reads concatenated instead - and `clean_reads_for` already knew how to pick those.
    The failure was an assembly planned against files nothing created, reported as "does not exist".
    """
    from metawrap2.commands import run as run_cmd

    calls = []
    monkeypatch.setattr(
        run_cmd, "concatenate_for_coassembly", lambda layout, study, dry: calls.append(study.steps)
    )
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)

    class _Study:
        coassemble = True
        steps = ("assembly", "binning")

    class _Step:
        module = "assembly"

    run_cmd._before_group(_Step(), object(), _Study(), False)
    assert calls == [("assembly", "binning")], "the assembly step must trigger the concatenation"


def test_a_per_sample_study_does_not_concatenate_anything(tmp_path, monkeypatch):
    """Concatenating for a study that assembles each sample separately would be pure waste."""
    from metawrap2.commands import run as run_cmd

    calls = []
    monkeypatch.setattr(
        run_cmd, "concatenate_for_coassembly", lambda layout, study, dry: calls.append(1)
    )
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)

    class _Study:
        coassemble = False
        steps = ("assembly",)

    class _Step:
        module = "assembly"

    run_cmd._before_group(_Step(), object(), _Study(), False)
    assert calls == []


def test_only_the_assembly_step_triggers_it(monkeypatch):
    """A binning step must not re-concatenate 99 GB of reads it has no use for."""
    from metawrap2.commands import run as run_cmd

    calls = []
    monkeypatch.setattr(
        run_cmd, "concatenate_for_coassembly", lambda layout, study, dry: calls.append(1)
    )
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)

    class _Study:
        coassemble = True
        steps = ("assembly", "binning")

    class _Step:
        module = "binning"

    run_cmd._before_group(_Step(), object(), _Study(), False)
    assert calls == []


# --- --preset in a study sheet -----------------------------------------------------------------


def test_a_preset_in_a_sheet_decides_which_bin_sets_refinement_is_given():
    """A sheet saying --preset very-fast was handed to refinement as the *default* pair.

    selected_binners scanned for --<binner> flags and knew nothing about --preset, so every preset but
    `fast` fell through to metaBAT2 + MaxBin2 and died at refinement with "maxbin2_bins is not a
    directory" - after paying for the binning. The resolution comes from binning.PRESETS so there is one
    table, not two.
    """
    from metawrap2.commands.run import selected_binners
    from metawrap2.modules import binning

    for tier, expected in binning.PRESETS.items():
        got = selected_binners(["--preset", tier])
        assert set(got) == set(expected), "%s resolved to %s" % (tier, got)


def test_a_named_binner_still_wins_over_a_preset_in_a_sheet():
    """Same precedence the binning module applies, or the driver and the module would disagree."""
    from metawrap2.commands.run import selected_binners

    assert selected_binners(["--preset", "slow", "--metabat2"]) == ["metabat2"]


def test_an_unknown_or_truncated_preset_falls_back_to_the_default_pair():
    """A sheet is user input; a bad value must not make the driver index past the end of the list."""
    from metawrap2.commands.run import selected_binners
    from metawrap2.modules.binning import DEFAULT_PRESET, PRESETS

    fallback = list(PRESETS[DEFAULT_PRESET])
    assert selected_binners(["--preset"]) == fallback
    assert selected_binners(["--preset", "nonsense"]) == fallback
