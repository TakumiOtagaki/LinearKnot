#!/usr/bin/env python3
"""Compare full-precision legacy inside with scfg2 LegacyCompatible (no BPP)."""

import argparse
import json
import math
from pathlib import Path
import random
import subprocess
import sys

from gen_random_seq_struct import assign_sequence, generate_structure, validate_pairing

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "test/oracle/legacy_cparty"))
from oracle_config import digest, verify_checkout


def fixed_cases():
    # Closing neighbors are inside (i,j), including i==1 / j==n. Cover
    # all canonical orientations, four boundary placements and three G roles.
    core = "GAGGGAAAACCCAGGGAAAACCCAC"
    structure = "(.(((....))).(((....))).)"
    for pair in ("GC", "CG", "AU", "UA", "GU", "UG"):
        paired = pair[0] + core[1:-1] + pair[1]
        for left in (0, 1):
            for right in (0, 1):
                sequence = "A" * left + paired + "A" * right
                outer = "." * left + "(" + "." * (len(core) - 2) + ")" + "." * right
                full = "." * left + structure + "." * right
                for role, scaffold in (("empty", "." * len(sequence)),
                                       ("outer", outer), ("full", full)):
                    yield f"ml_closing_{pair}_pad{left}{right}_{role}", sequence, scaffold
    yield ("wip_shortest_root29", "GAGAAACAAGAAGAAAAAAAAACCAAAAC",
           "(.(...)..(............)).....")
    yield "wip_binary13", "AAGAAACAGAAAC", "..(...).(...)"
    for stem in (20, 35, 50):
        fixture = ROOT / f"test/fixtures/scfg2_wmbp_leftward_be_20261008/stem{stem}.txt"
        sequence, scaffold = fixture.read_text().splitlines()[:2]
        yield f"wmbp_be_stem{stem}", sequence, scaffold
    yield "closure14", "GCGCGCGCGCGCGC", "(.(......)...)"
    yield "vp_split25", "GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)..."
    yield ("omitted_base45", "GCAUAUUAAUAUUCUGCUAUAAAUAUAGGUGCACUUAGGGCAUGC",
           "((.(((....))(.(.....))).((....((.(....)))))))")
    yield ("long_be65", "UCUACGUCAGUAGAUCUAACCGUGGAUAGCUUCUACAGGACAUGGGAUCCCAUUGGGUAAGACUU",
           "(.(X..(.X.X.X)..X).X.X..X..XX..X.XX.X.....X......X.X..).X(.XX)...")
    for index, g in enumerate(("((....))((...))(...)", "((x...))((...))(...)",
                               "((X...))((...))(...)", "((....))((.x.))(...)",
                               "((....))((.X.))(...)")):
        yield f"vpr_guard20_{index}", "UUAUAUGGGUAUUAUUCUCA", g
    # This scaffold retains the outer/inner pairs with a seven-base hairpin.
    # n=41 tests the size-30 boundary. Pure size-31 V roots have Z=0 under
    # the unchanged MAXLOOP search window; long BE roots use long_be65 above.
    for left, right in ((16, 14),):
        g = "(" + "." * left + "(" + "." * 7 + ")" + "." * right + ")"
        sequence = list("A" * len(g))
        for i in (0, left + 1):
            sequence[i] = "G"
        for i in (left + 9, len(g) - 1):
            sequence[i] = "C"
        yield f"be_loop_{left}_{right}", "".join(sequence), g


def observe(command):
    return subprocess.check_output(command, cwd=ROOT, text=True, stderr=subprocess.PIPE)


def legacy(binary, sequence, scaffold, params):
    result = json.loads(observe([str(binary), sequence, scaffold, str(params)]))
    result.setdefault("log_z", -result["energy"] / result["rt"])
    if not all(math.isfinite(result[key]) for key in ("energy", "log_z", "rt")):
        raise ValueError("legacy inside returned non-finite values")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle-checkout", type=Path, default=ROOT / "build/legacy_cparty_oracle/src")
    parser.add_argument("--scfg2-tool", type=Path, default=ROOT / "build/scfg2_restricted_partition_function_tool")
    parser.add_argument("--param-file", type=Path, default=ROOT / "params/rna_Turner04.par")
    parser.add_argument("--previous-probe", type=Path, help="optional preserved inside-only binary from the previous oracle")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261011)
    parser.add_argument("--count", type=int, default=100, help="random cases per length")
    parser.add_argument("--lengths", default="20,30,50,80,100")
    parser.add_argument("--log-tolerance", type=float, default=1e-9)
    parser.add_argument("--energy-tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    lengths = [int(value) for value in args.lengths.split(",")]
    if (args.count < 1 or not lengths or min(lengths) < 1 or
            not math.isfinite(args.log_tolerance) or args.log_tolerance < 0 or
            not math.isfinite(args.energy_tolerance) or args.energy_tolerance < 0):
        parser.error("positive counts/lengths and finite nonnegative tolerances required")
    out, checkout = args.out_dir.resolve(), args.oracle_checkout.resolve()
    if out == checkout or checkout in out.parents:
        parser.error("out-dir must be separate from the oracle checkout")
    manifest = verify_checkout(checkout)
    out.mkdir(parents=True, exist_ok=True)
    cases = list(fixed_cases())
    _, first_seq, first_g = cases[0]
    identity = json.loads(observe([
        sys.executable, str(ROOT / "test/tools/legacy_vp_recursive_split_probe.py"),
        first_seq, first_g, "--oracle-dir", str(checkout),
        "--output-dir", str(out / "probe"), "--param-file", str(args.param_file.resolve())]))
    binary = out / "probe/inside-probe"
    previous = args.previous_probe.resolve() if args.previous_probe else None
    tool, params = args.scfg2_tool.resolve(), args.param_file.resolve()
    metadata = {
        "oracle_manifest": manifest, "probe_identity": identity,
        "oracle_manifest_sha256": digest(checkout.parent / "manifest.json"),
        "probe_sha256": digest(binary), "scfg2_tool_sha256": digest(tool),
        "previous_probe_sha256": digest(previous) if previous else None,
        "parameter_file_sha256": digest(params), "script_sha256": digest(Path(__file__)),
        "scfg2_commit": observe(["git", "rev-parse", "HEAD"]).strip(),
        "seed": args.seed, "count_per_length": args.count, "lengths": lengths,
        "log_tolerance": args.log_tolerance, "energy_tolerance": args.energy_tolerance,
        "wmb_closure_policy": "legacy-compatible-wmb", "mfe_seed": 0.0,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    rng = random.Random(args.seed)
    for n in lengths:
        for index in range(args.count):
            g = generate_structure(rng, n, 3)
            sequence = assign_sequence(rng, g)
            validate_pairing(sequence, g)
            name = f"random{n}_{index:04d}"
            cases.append((name, sequence, g))
            if index % 25 == 0:
                for marker in ("x", "X"):
                    marked = "".join(marker if ch == "." and rng.random() < 0.35 else ch for ch in g)
                    cases.append((name + "_" + marker, sequence, marked))
    rows = []
    for name, sequence, g in cases:
        row = {"case_id": name, "sequence": sequence, "scaffold": g}
        try:
            current = legacy(binary, sequence, g, params)
            raw = observe([str(tool), "--details", "--legacy-compatible", sequence, g, str(params)])
            values = dict(part.split("=", 1) for part in raw.strip().split("\t"))
            energy, log_z = float(values["pf_energy"]), float(values["log_partition_scaled"])
            if (values["wmb_closure_policy"] != "legacy-compatible-wmb" or
                    not all(math.isfinite(value) for value in (energy, log_z))):
                raise ValueError("invalid scfg2 policy or inside result")
            row.update(legacy=current, scfg2_energy=energy, scfg2_log_z=log_z,
                       energy_error=abs(energy - current["energy"]),
                       log_error=abs(log_z - current["log_z"]))
            if previous:
                old = legacy(previous, sequence, g, params)
                row.update(previous=old, previous_energy_error=abs(old["energy"] - current["energy"]),
                           previous_log_error=abs(old["log_z"] - current["log_z"]))
            row["passed"] = (row["energy_error"] <= args.energy_tolerance and
                             row["log_error"] <= args.log_tolerance and
                             row.get("previous_energy_error", 0.0) <= args.energy_tolerance and
                             row.get("previous_log_error", 0.0) <= args.log_tolerance and
                             current["bp_zero_count"] == 0)
        except (subprocess.CalledProcessError, ValueError, KeyError) as error:
            row.update(passed=False, error=str(error))
            if isinstance(error, subprocess.CalledProcessError):
                row["stderr"] = error.stderr
        rows.append(row)
        if len(rows) % 50 == 0:
            print(f"checked {len(rows)}/{len(cases)} cases", flush=True)
    (out / "cases.json").write_text(json.dumps(rows, indent=2) + "\n")
    summary = {
        "cases": len(rows), "passed": sum(row["passed"] for row in rows),
        "failed_case_ids": [row["case_id"] for row in rows if not row["passed"]],
        "bp_zero_count": sum(row.get("legacy", {}).get("bp_zero_count", 0) for row in rows),
        "bp_spans_checked": sum(len(g) * (len(g) + 1) // 2 for _, _, g in cases),
    }
    for key in ("energy_error", "log_error", "previous_energy_error", "previous_log_error"):
        summary["max_" + key] = max((row[key] for row in rows if key in row), default=None)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] == summary["cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
