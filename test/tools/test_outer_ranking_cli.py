#!/usr/bin/env python3
"""Exercise outer/default/baseline ranking across the public toolset."""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SCFG2_BUILD_DIR", ROOT / "build"))
ENV = dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(ROOT / "params/rna_Turner04.par"))
PROGRAMS = {
    "fold": [],
    "profile": [],
    "access": ["--access-len", "1,3", "--probability"],
    "sample": ["--count", "4", "--seed", "73"],
    "subopt": ["--top-k", "3", "--stats"],
    "mea": ["--gamma", "1.5", "--stats"],
    "centroid": ["--gamma", "1.5", "--stats"],
}
SMALL = ["--sequence", "GCGCGCGCGCGCGC", "--scaffold", "(.(......)...)", "--beam", "3"]


def run(program, *args, check=True):
    result = subprocess.run([str(BUILD / "linearknot"), *([program] if program is not None else []), *args], cwd=ROOT, env=ENV,
                            capture_output=True, text=True, timeout=180)
    if check:
        assert result.returncode == 0, (program, args, result.stdout, result.stderr)
    return result


def stable(text):
    # Raccess exposes actual elapsed time in its metadata.
    return "\n".join(line for line in text.splitlines()
                     if not line.startswith("# elapsed_seconds="))


def ensemble(text):
    return float(re.search(r"ensemble free energy: ([^ ]+)", text)[1])


def stats_root(result):
    return float(re.search(r"\blog_z=([^ ]+)", result.stderr)[1])


def access_rows(text):
    rows = {}
    for line in text.splitlines():
        if not line or not line[0].isdigit():
            continue
        pos, values = line.split("\t", 1)
        rows[int(pos)] = {int(length): float(value)
                          for token in values.split(";") if token
                          for length, value in [token.split(",")]}
    return rows


root_help = run("--help").stdout
assert root_help == run("help").stdout == run("-h").stdout
assert "LinearKnot" in root_help and "Usage: linearknot COMMAND" in root_help
assert "LinKnot" not in root_help
assert run("--version").stdout.strip() == "LinearKnot " + (ROOT / "VERSION").read_text().strip()
assert run(None, check=False).returncode != 0
for bad in ["no-such-command", "LinKnot", "linknot"]:
    result = run(bad, check=False)
    assert result.returncode != 0 and "unknown command" in result.stderr
assert run("help", "unknown", check=False).returncode != 0
assert run("help", "fold", "extra", check=False).returncode != 0

for program, extra in PROGRAMS.items():
    help_result = run(program, "--help")
    assert "Usage: linearknot " + program in help_result.stdout
    assert help_result.stdout == run("help", program).stdout
    assert "--ranking outer|baseline" in help_result.stdout
    assert "default outer" in help_result.stdout
    assert "--outer-pilot-exponent" in help_result.stdout
    default = run(program, *SMALL, *extra)
    explicit = run(program, *SMALL, *extra, "--ranking", "outer",
                   "--outer-pilot-exponent", "2/3")
    assert stable(default.stdout) == stable(explicit.stdout), program
    if "--stats" in extra:
        assert "ranking=outer" in default.stderr
        assert "outer_pilot_exponent=2/3" in default.stderr
        assert stats_root(default) == stats_root(explicit)
    run(program, *SMALL, *extra, "--ranking", "outer", "--outer-pilot-exponent", "0.5")
    baseline = run(program, *SMALL, *extra, "--ranking", "baseline")
    if "--stats" in extra:
        assert "ranking=baseline" in baseline.stderr
    exact = run(program, *SMALL, *extra, "--exact")
    if "--stats" in extra:
        assert "ranking=exact" in exact.stderr
    for bad, expected in [
        (["--ranking", "prefix"], "--ranking must be"),
        (["--ranking"], "missing value"),
        (["--outer-pilot-exponent", "0.6"], "--outer-pilot-exponent must be"),
        (["--outer-pilot-exponent", "nan"], "--outer-pilot-exponent must be"),
        (["--outer-pilot-exponent"], "missing value"),
        (["--ranking", "baseline", "--outer-pilot-exponent", "0.5"], "requires --ranking outer"),
        (["--exact", "--ranking", "outer"], "cannot be used with --exact"),
        (["--ranking", "baseline", "--exact"], "cannot be used with --exact"),
        (["--exact", "--outer-pilot-exponent", "2/3"], "cannot be used with --exact"),
    ]:
        result = run(program, *SMALL, *extra, *bad, check=False)
        assert result.returncode != 0 and expected in result.stderr, (program, bad, result.stderr)
    for beam in ["0", "-1", "1x", "9999999999999999999999999"]:
        result = run(program, *SMALL, *extra, "--beam", beam, check=False)
        assert result.returncode != 0, (program, beam)


fixture = ROOT / "test/fixtures/scfg2_outer_ranking_20261014/bpRNA_RFAM_27342.txt"
_, sequence, scaffold = fixture.read_text().splitlines()
args = ["--sequence", sequence, "--scaffold", scaffold, "--beam", "20"]
with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    default_bpp, explicit_bpp = directory / "default.bpp", directory / "explicit.bpp"
    default = run("fold", *args, "--bpp", str(default_bpp))
    explicit = run("fold", *args, "--ranking", "outer", "--bpp", str(explicit_bpp))
    baseline = run("fold", *args, "--ranking", "baseline")
    inside_only = run("fold", *args)
    assert default.stdout == explicit.stdout == inside_only.stdout
    assert default_bpp.read_bytes() == explicit_bpp.read_bytes()
    # Nonvacuous witness: the default actually changes retained PF mass.
    assert ensemble(default.stdout) < ensemble(baseline.stdout) - 1.0
    paired = [0.0] * (len(sequence) + 1)
    for line in default_bpp.read_text().splitlines():
        if not line or not line[0].isdigit():
            continue
        i, j, probability = line.split()
        i, j, probability = int(i), int(j), float(probability)
        assert i < j and 0.0 <= probability <= 1.0
        paired[i] += probability
        paired[j] += probability
    capr = run("profile", *args)
    stems = [float(p) for p in next(line for line in capr.stdout.splitlines()
                                  if line.startswith("Stem ")).split()[1:]]
    access = run("access", *args, *PROGRAMS["access"])
    assert "# beam_ranking=outer" in access.stdout
    assert "# outer_pilot_exponent=2/3" in access.stdout
    access_values = access_rows(access.stdout)
    assert len(stems) == len(sequence) == len(access_values)
    for i, stem in enumerate(stems, 1):
        assert math.isclose(stem, paired[i], abs_tol=2e-9)
        assert math.isclose(access_values[i - 1][1], 1.0 - stem, abs_tol=2e-9)
    for program in ["subopt", "mea", "centroid"]:
        result = run(program, *args, *PROGRAMS[program])
        assert math.isclose(stats_root(result), 60.715913323849691, abs_tol=1e-11)
        assert "ranking=outer" in result.stderr

    batch = directory / "batch.txt"
    batch.write_text(">first\nAAAAAA\n......\n>second\ngcgcgcgcgcgcgc\n(.(......)...)\n")
    for program, extra in PROGRAMS.items():
        batch_args = ["--input", str(batch), "--beam", "3", *extra]
        a = run(program, *batch_args)
        b = run(program, *batch_args, "--ranking", "outer")
        assert stable(a.stdout) == stable(b.stdout), program
        assert ">first\n" in a.stdout and ">second\n" in a.stdout
        if program != "fold":
            target = directory / (program + ".txt")
            saved = run(program, *batch_args, "--output", str(target))
            assert not saved.stdout and stable(target.read_text()) == stable(a.stdout), program

print("outer_ranking_cli=ok programs=7 default=outer exponents=2 batch=2")
