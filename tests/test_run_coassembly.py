def test_a_coassembly_without_read_qc_still_gets_its_concatenated_reads(tmp_path, monkeypatch):
    """The co-assembly input ALL_READS is built even when the study has no read_qc step.

    Previously it was built only after read_qc, so a study that skipped QC (as benchmark arms that
    hold QC constant do) failed with "does not exist". `clean_reads_for` picks the raw reads then.
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

    # A binning step does not concatenate the reads (99 GB in one study), and a per-sample study
    # has nothing to concatenate.
    _Step.module = "binning"
    run_cmd._before_group(_Step(), object(), _Study(), False)
    _Step.module = "assembly"
    _Study.coassemble = False
    run_cmd._before_group(_Step(), object(), _Study(), False)
    assert calls == [("assembly", "binning")]


# --- --preset in a study sheet -----------------------------------------------------------------


def test_a_preset_in_a_sheet_decides_which_bin_sets_refinement_is_given():
    """selected_binners resolves --preset in a study sheet through binning.PRESETS.

    Previously it read only --<binner> flags, so every preset but `fast` resolved to metaBAT2 +
    MaxBin2 and refinement failed after binning with "maxbin2_bins is not a directory".
    """
    from metawrap2.commands.run import selected_binners
    from metawrap2.modules import binning

    for tier, expected in binning.PRESETS.items():
        got = selected_binners(["--preset", tier])
        assert set(got) == set(expected), "%s resolved to %s" % (tier, got)
    # Same precedence the binning module applies, or the driver and the module would disagree.
    assert selected_binners(["--preset", "gpu", "--metabat2"]) == ["metabat2"]
