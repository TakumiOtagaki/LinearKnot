#!/usr/bin/env python3
"""Exercise MEA and gamma-centroid CLI scoring, validation and PK contracts."""
import math
import os
import re
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SCFG2_BUILD_DIR", ROOT / "build"))


def run(*args):
    return subprocess.run([str(BINARY), name, *args], cwd=ROOT,
                          capture_output=True, text=True,
                          env=dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(ROOT / "params/rna_Turner04.par")))


for name in ["mea", "centroid"]:
    BINARY = BUILD / "linearknot"
    centroid = name == "centroid"
    plain = ["--sequence", "AAAAAA", "--scaffold", "......"]
    for extra in [[], ["--exact"], ["--gamma", "0"], ["--gamma", "3", "--stats"]]:
        result = run(*plain, *extra)
        assert result.returncode == 0, result.stderr
        expected = "...... (EG=15; gamma=" if centroid else "...... (EA=6; gamma="
        assert expected in result.stdout, result.stdout
        if "--stats" in extra:
            assert "decoding_seconds=" in result.stderr and "terms=" in result.stderr

    for extra in [["--gamma", "-1"], ["--gamma", "nan"], ["--gamma", "inf"],
                  ["--gamma", "1trailing"], ["--gamma"], ["--beam", "0"],
                  ["--beam", "-1"], ["--beam", "abc"]]:
        assert run(*plain, *extra).returncode != 0, extra

    sequence, scaffold = list("A" * 24), list("." * 24)
    for k in range(3):
        sequence[3 + k] = sequence[9 + k] = "G"
        sequence[14 - k] = sequence[20 - k] = "C"
        scaffold[3 + k], scaffold[14 - k] = "(", ")"
    sequence, scaffold = "".join(sequence), "".join(scaffold)
    for mode in [[], ["--exact"]]:
        result = run("--sequence", sequence, "--scaffold", scaffold,
                     "--gamma", "1000", *mode)
        assert result.returncode == 0, result.stderr
        line = result.stdout.splitlines()[3]
        prediction = line.split()[0]
        if centroid:
            values = dict(re.findall(r"(EG|ETP|ETN|gamma)=([^;)]+)", line))
            assert math.isclose(float(values["EG"]), float(values["gamma"]) * float(values["ETP"]) + float(values["ETN"]), abs_tol=1e-8), line
        assert "[" in prediction, prediction
        assert len(prediction) == len(scaffold)
        for a, b in zip(scaffold, prediction):
            if a != ".":
                assert a == b, prediction

    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "input.txt"
        target = Path(directory) / "output.txt"
        source.write_text(">rna1\nAAAAAA\n......\n>rna2\nAAAAAA\n......\n")
        result = run("-i", str(source), "-o", str(target))
        assert result.returncode == 0 and not result.stdout, result.stderr
        text = target.read_text()
        assert ">rna1\n" in text and ">rna2\n" in text
        assert text.count("EG=15;" if centroid else "EA=6;") == 2, text

    # A fixed witness where the two objectives make strictly different choices.
    result = run("--sequence", "GGGAAACCCAGGGAAACCC", "--scaffold", "(((...)))..........",
                 "--gamma", "0.5", "--exact")
    assert result.returncode == 0, result.stderr
    prediction = result.stdout.splitlines()[3].split()[0]
    assert prediction == ("(((...))).........." if centroid else "(((...))).((....))."), prediction
    print(name + "_cli=ok")
