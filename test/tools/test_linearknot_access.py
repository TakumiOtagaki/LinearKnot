#!/usr/bin/env python3
"""Check joint loop-window probabilities, CapR agreement and CLI compatibility."""
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SCFG2_BUILD_DIR", ROOT / "build"))
CODES = ("E", "H", "B", "I", "M", "P")
LABELS = ("Exterior", "Hairpin", "Bulge", "Internal", "Multiloop", "PseudoLoop")
LENGTHS = (1, 2, 3, 5, 40, 2)  # duplicates, end boundaries, oversized window


def run(tool, *args):
    return subprocess.run([str(BUILD / "linearknot"), *([tool] if tool is not None else []), *args], cwd=ROOT,
                          capture_output=True, text=True, env=dict(os.environ,
                          CPARTY_FIXED_ENERGY_PARAM_FILE=str(ROOT / "params/rna_Turner04.par")))


def access(args, *options):
    result = run("access", *args, *options)
    assert result.returncode == 0, result.stderr
    return result.stdout


def rows(text):
    return [line.split("\t") for line in text.splitlines() if line[:1].isdigit()]


def total_column(row):
    return [(int(w), float(p)) for w, p in
            (token.split(",") for token in row[1].rstrip(";").split(";"))]


def loop_column(row):
    assert len(row) == 3 and row[2].startswith("byloop:"), row
    result = []
    for token in row[2][7:].rstrip(";").split(";"):
        length, *values = token.split(",")
        assert tuple(value.split("=")[0] for value in values) == CODES, row
        result.append((int(length), [float(value.split("=")[1]) for value in values]))
    return result


def near(a, b):
    assert math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-9), (a, b)


fixtures = [
    ("AAAAAA", "......"),
    ("GGGGGAACCAACCC", "((.....))....."),
    ("GAAGAAACC", "(..(...))"),
    ("GAAGAAACAAC", "(..(...)..)"),
    ("AAGAAGAAACAAGAAACAACAA", "..(..(...)..(...)..).."),
]
covered = set()
for mode in (["--exact"], ["--beam", "5"], ["--beam", "5", "--ranking", "baseline"]):
    for sequence, scaffold in fixtures:
        args = ["--sequence", sequence, "--scaffold", scaffold,
                "--access-len", ",".join(map(str, LENGTHS)), *mode]
        plain = access(args, "--probability")
        grouped = access(args, "--probability", "--by-loop")
        energy = access(args, "--by-loop")
        # Enabling grouping leaves the legacy columns byte-identical.
        assert [r[:2] for r in rows(grouped)] == rows(plain)
        assert "# byloop_value=joint_probability\n" in grouped
        assert "P:PseudoLoop\n" in grouped
        assert "byloop" not in plain
        capr = run("profile", "--sequence", sequence, "--scaffold", scaffold, *mode)
        assert capr.returncode == 0, capr.stderr
        context = {fields[0]: list(map(float, fields[1:])) for line in capr.stdout.splitlines()
                   if (fields := line.split()) and fields[0] in LABELS}
        rt = float(next(line.split("=", 1)[1] for line in energy.splitlines()
                        if line.startswith("# folding_rt=")))
        for start, (row, ed_row) in enumerate(zip(rows(grouped), rows(energy))):
            assert int(row[0]) == start
            expected_lengths = [w for w in LENGTHS if start + w <= len(sequence)]
            assert [w for w, _ in total_column(row)] == expected_lengths
            assert [w for w, _ in loop_column(row)] == expected_lengths
            assert row[2] == ed_row[2]  # loop values always probabilities, even in ED mode
            for (w, total), (loop_w, probabilities), (_, ed) in zip(
                    total_column(row), loop_column(row), total_column(ed_row)):
                assert w == loop_w
                assert all(0.0 <= p <= total + 1e-9 for p in probabilities)
                near(sum(probabilities), total)
                if total > 0.0:
                    near(ed, -rt * math.log(total))
                else:
                    assert ed == math.inf
                if w == 1:
                    for label, p in zip(LABELS, probabilities):
                        near(p, context[label][start])
                elif mode == ["--exact"]:
                    covered.update(c for c, p in zip(CODES, probabilities) if p > 1e-10)
                if sequence == "AAAAAA":
                    assert probabilities == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    assert total == 1.0
assert covered == set(CODES), covered

plain = ["--sequence", "AAAAAA", "--scaffold", "......", "--access-len", "2,9"]
expected = rows(access(plain, "--by-loop"))
for alias in ("--byloop", "-byloop"):
    assert rows(access(plain, alias)) == expected
assert expected[-1] == ["5"]  # no valid window, no dangling byloop column
assert rows(access([*plain, "--access-len", "100"], "--by-loop"))[-1] == ["5"]
only_large = ["--sequence", "AAAAAA", "--scaffold", "......", "--access-len", "9"]
assert rows(access(only_large, "--by-loop")) == [[str(i)] for i in range(6)]
for bad in (["--access-len", "0"], ["--access-len", "-1"], ["--access-len", "abc"],
            ["--beam", "0"], ["--exact", "--ranking", "outer"], ["--unknown"]):
    assert run("access", *plain, "--by-loop", *bad).returncode != 0, bad

with tempfile.TemporaryDirectory() as directory:
    source, target = Path(directory) / "input.txt", Path(directory) / "output.txt"
    source.write_text(">rna1\nAAAAAA\n......\n>rna2\ncuucgg\n(....)\n")
    result = run("access", "-i", str(source), "-o", str(target),
                 "--exact", "--access-len", "1,3", "--by-loop", "--probability")
    assert result.returncode == 0 and not result.stdout, result.stderr
    text = target.read_text()
    assert ">rna1\n" in text and ">rna2\n" in text
    assert len(rows(text)) == 12 and all(len(row) == 3 for row in rows(text))
    assert ",H=1," in text

print("linearknot_access_cli=ok")
