#!/usr/bin/env python3
"""Resource regressions for exact outside and zero-mass sampling records."""
import os
from pathlib import Path
import resource
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(os.environ.get("SCFG2_BUILD_DIR", ROOT / "build"))
ENV = dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(ROOT / "params/rna_Turner04.par"))
LIMIT = 256 * 1024**2


def limit_address_space():
    _, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = LIMIT if hard == resource.RLIM_INFINITY else min(LIMIT, hard)
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))


with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    args = ["--sequence", "GCAU" * 40, "--scaffold", "." * 160, "--exact"]
    outputs = []
    for name, extra in [
        ("fold", []),
        ("fold", ["--bpp", str(directory / "pairs.tsv"), "--bpp-cutoff", "0"]),
        ("profile", []),
        ("access", ["--access-len", "1,17,160,161", "--probability"]),
    ]:
        result = subprocess.run([str(BUILD / "linearknot"), name, *args, *extra], cwd=ROOT, env=ENV,
                                capture_output=True, text=True, timeout=180,
                                preexec_fn=limit_address_space)
        assert result.returncode == 0, (name, extra, result.stdout, result.stderr)
        outputs.append(result.stdout)
    assert outputs[0] == outputs[1], "inside-only and BPP CLI metadata changed"
    assert (directory / "pairs.tsv").stat().st_size > 0

# This ensemble has exactly one all-unpaired structure. Dead split candidates
# must not make an arithmetic tape exhaust the same address-space budget.
sequence = "A" * 384
scaffold = "." * len(sequence)
for count in (0, 1):
    result = subprocess.run(
        [str(BUILD / "linearknot"), "sample", "--sequence", sequence,
         "--scaffold", scaffold, "--exact", "--count", str(count), "--seed", "73"],
        cwd=ROOT, env=ENV, capture_output=True, text=True, timeout=180,
        preexec_fn=limit_address_space)
    assert result.returncode == 0, (count, result.stdout, result.stderr)
    expected = [">seq1", sequence, *([scaffold] * count), ""]
    assert result.stdout.splitlines() == expected, (count, result.stdout)
    assert not result.stderr, result.stderr

print("scfg2_exact_outside_memory=ok n=160 address_space_mib=256 consumers=4 "
      "sampling_n=384 sampling_counts=0,1")
