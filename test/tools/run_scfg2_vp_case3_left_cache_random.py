#!/usr/bin/env python3
"""Check exact C3_LEFT KeepAll factorization on random cases."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from typing import List, Sequence


def comma_ints(text: str) -> List[int]:
    values = sorted({int(token) for token in text.split(",") if token})
    if not values or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return values


def run(command: Sequence[str]) -> str:
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr}"
        )
    return completed.stdout


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", default="build/scfg2_vp_case3_left_cache_test")
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\ttails\tleft_candidates\tsegments\tsources\twi_children"
        "\tchild_target_groups\ttarget_runs\tsignature_segments"
        "\tmax_child_fanout\tmax_same_target\tmax_target_runs"
        "\tmax_signature_segments\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {key: 0 for key in (
            "tails", "left_candidates", "segments", "sources", "wi_children",
            "child_target_groups", "target_runs", "signature_segments"
        )}
        max_child_fanout = 0
        max_same_target = 0
        max_target_runs = 0
        max_signature_segments = 0
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            failures += int(row["parity"]) != 1
            for key in totals:
                totals[key] += int(row[key])
            max_child_fanout = max(max_child_fanout, int(row["max_child_fanout"]))
            max_same_target = max(max_same_target, int(row["max_same_target"]))
            max_target_runs = max(max_target_runs, int(row["max_target_runs"]))
            max_signature_segments = max(
                max_signature_segments, int(row["max_signature_segments"])
            )
        print(
            f"{length}\t{args.count}\t{totals['tails']}\t"
            f"{totals['left_candidates']}\t{totals['segments']}\t"
            f"{totals['sources']}\t{totals['wi_children']}\t"
            f"{totals['child_target_groups']}\t{totals['target_runs']}\t"
            f"{totals['signature_segments']}\t{max_child_fanout}\t"
            f"{max_same_target}\t{max_target_runs}\t"
            f"{max_signature_segments}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
