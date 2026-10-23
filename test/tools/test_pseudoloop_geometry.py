#!/usr/bin/env python3
"""Compare production derivation labels with the independent union-pair reference."""

import json
import os
from pathlib import Path
import re
import subprocess
import sys

from pseudoloop_context_reference import classify, pairs_from_db


ROOT = Path(__file__).resolve().parents[2]
BUILD = Path(sys.argv[1] if len(sys.argv) > 1 else
             os.environ.get("SCFG2_BUILD_DIR", ROOT / "build")).resolve()
CASES = json.loads((ROOT / "test/fixtures/pseudoloop_context_change.json").read_text())["cases"]
TARGETS = {c["name"]: c for c in CASES if "sequence" in c}
CHANGED = {"WIP_EXTEND_UNPAIRED", "WIP_BASEPAIR_V", "WIP_BASEPAIR_WMB",
           "VPL_SPLIT_VP", "VPR_SPLIT_VP_BASEPAIR", "BE_WIP_BASEPAIR", "BE_BASEPAIR_WIP"}

for case in CASES:
    assert classify(case["length"], case["pairs"])[0] == case["expected"], case["name"]

result = subprocess.run(
    [str(BUILD / "scfg2_pk_profiles_oracle_test"), "--geometry-records"], cwd=ROOT,
    input="".join(f'{c["name"]}\t{c["sequence"]}\t{c["scaffold"]}\n' for c in TARGETS.values()),
    capture_output=True, text=True, timeout=180, check=True,
)
counts, found, covered = {}, {}, {}
for line in result.stdout.splitlines():
    if not line.startswith("RECORD\t"):
        continue
    _, name, policy, sequence, scaffold, emitted, generated, emitters = line.split("\t")
    pairs = pairs_from_db(scaffold) | {
        tuple(map(int, pair.split("-"))) for pair in generated.split(",") if pair}
    want, explanation = classify(len(sequence), pairs)
    got = emitted.replace("-", "S")
    assert got == want, (name, policy, sorted(pairs), got, want, explanation)
    stats = counts.setdefault(policy, {"groups": 0, "bases": 0, "free_windows": 0})
    stats["groups"] += 1
    stats["bases"] += len(sequence)
    for run in filter(None, re.split("S+", got)):
        assert len(set(run)) == 1, (name, policy, "mixed fully unpaired run", run)
        stats["free_windows"] += len(run) * (len(run) + 1) // 2
    covered.setdefault(policy, set()).update(emitters.rstrip(",").split(","))
    if name in TARGETS:
        case = TARGETS[name]
        assert (sequence, scaffold) == (case["sequence"], case["scaffold"])
        if pairs == {tuple(pair) for pair in case["pairs"]}:
            assert got == case["expected"], name
            found.setdefault(policy, set()).add(name)

assert len(counts) == 2, counts
for policy in counts:
    assert found.get(policy, set()) == set(TARGETS), (policy, found)
    assert CHANGED <= covered[policy], (policy, CHANGED - covered[policy])
print(json.dumps({"geometry_fixtures": len(CASES), "production": counts,
                  "changed_rules_covered_per_policy": len(CHANGED), "status": "ok"}))
