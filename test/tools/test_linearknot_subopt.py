#!/usr/bin/env python3
"""Exercise ranked suboptimals, band boundaries, errors and file I/O."""
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SCFG2_BUILD_DIR", ROOT / "build"))
BINARY = BUILD / "linearknot"


def run(*args):
    return subprocess.run([str(BINARY), "subopt", *args], cwd=ROOT, capture_output=True,
                          text=True, env=dict(os.environ,
                          CPARTY_FIXED_ENERGY_PARAM_FILE=str(ROOT / "params/rna_Turner04.par")))


def predictions(result):
    assert result.returncode == 0, result.stderr
    return [(structure, float(energy)) for structure, energy in
            re.findall(r"^([.()\[\]_xX]+) \(([^)]+)\)$", result.stdout, re.M)]


plain = ["--sequence", "AAAAAA", "--scaffold", "......"]
for mode in [[], ["--exact"]]:
    result = run(*plain, *mode, "--stats")
    assert predictions(result) == [("......", 0.0)], result.stdout
    assert "stop=exhausted" in result.stderr, result.stderr
    assert "scope=" + ("exact" if mode else "inside_retained_beam") in result.stderr

sequence, scaffold = "GCGCGCGCGCGCGC", "(.(......)...)"
args = ["--sequence", sequence, "--scaffold", scaffold]
for mode in [[], ["--exact"], ["--beam", "1"]]:
    all_result = run(*args, *mode, "--top-k", "10000", "--stats")
    full = predictions(all_result)
    assert len(full) > 1 and len({s for s, _ in full}) == len(full), full
    assert [e for _, e in full] == sorted(e for _, e in full), full
    assert "stop=exhausted" in all_result.stderr, all_result.stderr
    top_result = run(*args, *mode, "-k", "2", "--stats")
    assert predictions(top_result) == full[:2]
    assert "stop=top_k" in top_result.stderr, top_result.stderr
    # Inclusive boundary at an actually observed energy, including all ties.
    delta = full[1][1] - full[0][1]
    expected = [(s, e) for s, e in full if e - full[0][1] <= delta + 1e-9]
    band = run(*args, *mode, "-k", "10000", "--delta-energy", str(delta), "--stats")
    assert predictions(band) == expected, (band.stdout, expected)
    for structure, _ in full:
        assert len(structure) == len(scaffold)
        assert all(g == "." or g == s for g, s in zip(scaffold, structure))

for bad in [["-k", "0"], ["-k", "-1"], ["-k", "abc"], ["-k", "1x"],
            ["-b", "0"], ["-b", "-1"], ["-b", "1x"],
            ["-e", "-1"], ["-e", "nan"], ["-e", "inf"], ["-e", "1x"],
            ["-e"], ["--unknown"]]:
    assert run(*plain, *bad).returncode != 0, bad
assert run("--sequence", "AAAAAA", "--scaffold", "[....]").returncode != 0

with tempfile.TemporaryDirectory() as directory:
    source, target = Path(directory) / "input.txt", Path(directory) / "output.txt"
    source.write_text(">rna1\nAAAAAA\n......\n>rna2\ncuucgg\n(....)\n")
    result = run("-i", str(source), "-o", str(target), "--exact", "-k", "1")
    assert result.returncode == 0 and not result.stdout, result.stderr
    text = target.read_text()
    assert ">rna1\n" in text and ">rna2\nCUUCGG\n" in text, text
    assert "(....) (" in text, text

print("linearknot_subopt_cli=ok")
