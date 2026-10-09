"""One `--device` for all three learned binners, and what each tool is actually told.

The three spell their accelerator switch three different ways - SemiBin2 takes a word, VAMB takes a
boolean, COMEBin takes a device string - which is why there used to be a flag each. The translation is
what these pin, along with the fallback behaviour of `auto`, whose whole job is to not fail.
"""

from __future__ import annotations

from metawrap2.modules import binning


def test_the_default_is_the_cpu_and_that_is_on_purpose():
    """A GPU run and a CPU run of the same seeded network need not agree.

    Not a claim about hardware: the machine these were measured on has eight usable GPUs. It is that a
    study re-analysed elsewhere should not produce different MAGs for no visible reason.
    """
    assert binning.DEFAULT_LEARNED_DEVICE == "cpu"
    assert set(binning.LEARNED_DEVICES) == {"cpu", "cuda", "auto"}


def test_an_explicit_device_is_passed_through_without_a_probe(monkeypatch):
    """cpu and cuda are answers already; only `auto` is a question. cuda:N has to survive intact."""
    called = []
    monkeypatch.setattr(
        binning, "capture_in_env", lambda env, argv: called.append(argv) or "cuda\n"
    )
    for device in ("cpu", "cuda", "cuda:3"):
        assert binning.resolve_device(device, "metawrap2-learned") == device
    assert called == [], "an explicit device must not shell into the env to ask"


def test_auto_uses_the_gpu_when_the_binners_own_env_can_see_one(monkeypatch):
    """Asked of that binner's env, not of this process: each carries its own torch."""
    seen = {}

    def fake(env, argv):
        seen["env"] = env
        seen["argv"] = argv
        return "cuda\n"

    monkeypatch.setattr(binning, "capture_in_env", fake)
    assert binning.resolve_device("auto", "metawrap2-learned") == "cuda"
    assert seen["env"] == "metawrap2-learned"
    assert seen["argv"][0] == "python" and "torch.cuda.is_available" in seen["argv"][-1]


def test_auto_falls_back_to_the_cpu_when_there_is_no_gpu(monkeypatch):
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: "cpu\n")
    assert binning.resolve_device("auto", "metawrap2-learned") == "cpu"


def test_auto_falls_back_to_the_cpu_when_the_probe_fails(monkeypatch):
    """A broken or missing torch is exactly the case `auto` exists to survive.

    Three ways it can go wrong, all of which have to resolve rather than raise: no output at all (the
    env manager is not on PATH), a traceback on stdout, and an exception from the call itself.
    """
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: "")
    assert binning.resolve_device("auto", "metawrap2-learned") == "cpu"

    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: "ImportError: no torch\n")
    assert binning.resolve_device("auto", "metawrap2-learned") == "cpu"

    def boom(env, argv):
        raise OSError("no env manager")

    monkeypatch.setattr(binning, "capture_in_env", boom)
    assert binning.resolve_device("auto", "metawrap2-learned") == "cpu"


def test_semibin2_is_told_gpu_because_that_is_the_word_it_knows():
    """SemiBin2's --engine takes cpu/gpu; --device speaks cuda like the other two."""
    for device, expected in (
        ("cpu", "--engine cpu"),
        ("cuda", "--engine gpu"),
        ("cuda:1", "--engine gpu"),
    ):
        cmd = binning.SEMIBIN2.format(
            assembly="a.fa",
            bams="s.bam",
            out="out",
            threads=4,
            min_len=1000,
            seed=binning.SEMIBIN2_SEED,
            tmp="t",
            engine="gpu" if device.startswith("cuda") else "cpu",
        )
        assert expected in cmd


def test_vamb_is_told_with_a_boolean_flag_or_not_at_all():
    """VAMB has no device name - it has --cuda, present or absent."""
    on = binning.VAMB.format(
        out="o", assembly="a.fa", bamdir="w", min_len=250, threads=4,
        minfasta=binning.VAMB_MIN_FASTA, seed=binning.VAMB_SEED, cuda="--cuda",
    )  # fmt: skip
    off = binning.VAMB.format(
        out="o", assembly="a.fa", bamdir="w", min_len=250, threads=4,
        minfasta=binning.VAMB_MIN_FASTA, seed=binning.VAMB_SEED, cuda="",
    )  # fmt: skip
    assert "--cuda" in on
    assert "--cuda" not in off
    # And the rest of the command is unchanged by the switch, including the binsplit suppression.
    assert off.replace(" ", "").rstrip() == on.replace("--cuda", "").replace(" ", "").rstrip()


def test_comebin_is_told_the_device_string_verbatim():
    """COMEBin's -d takes cpu, cuda or cuda:N, so --device needs no translation at all."""
    for device in ("cpu", "cuda", "cuda:2"):
        cmd = binning.COMEBIN.format(
            assembly="a.fa", out="o", bamdir="w", threads=4, device=device,
            seed=binning.COMEBIN_SEED,
        )  # fmt: skip
        assert "-d %s " % device in cmd + " "


def test_the_old_per_tool_flags_are_gone_rather_than_aliased():
    """New software: a second way to say the same thing is a second thing to keep in step."""
    assert not hasattr(binning, "SEMIBIN2_ENGINE")
    assert not hasattr(binning, "COMEBIN_DEVICE")
    parsed = binning._parse_args(
        ["-a", "a.fa", "-o", "o", "--comebin", "--device", "cuda", "r_1.fastq", "r_2.fastq"]
    )
    assert parsed.device == "cuda"
    assert not hasattr(parsed, "semibin2_engine")
    assert not hasattr(parsed, "comebin_device")
