#!/usr/bin/env python3
"""Compare exact, coordinate, and pair-tree BE tail index setup time."""

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
    parser.add_argument(
        "--tool", default="build/scfg2_be_right_tail_index_profile_tool"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--repetitions", type=int, default=100)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\trepetitions\ttails\texact_plan_ms\tcoordinate_ms"
        "\tpair_tree_ms\ttree_over_coordinate\tcoordinate_visits"
        "\tpair_tree_visits\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        exact_ms = 0.0
        coordinate_ms = 0.0
        tree_ms = 0.0
        tails = 0
        coordinate_visits = 0
        tree_visits = 0
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([
                args.tool, case["seq"], case["G"], str(args.repetitions),
                args.param_file,
            ])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            exact_ms += float(row["exact_plan_ms"])
            coordinate_ms += float(row["coordinate_ms"])
            tree_ms += float(row["pair_tree_ms"])
            tails += int(row["tails"])
            coordinate_visits += int(row["coordinate_core_visits"])
            coordinate_visits += int(row["coordinate_anchor_visits"])
            tree_visits += int(row["pair_tree_ancestor_visits"])
            tree_visits += int(row["pair_tree_anchor_visits"])
            failures += int(row["parity"]) != 1
        exact_mean = exact_ms / args.count
        coordinate_mean = coordinate_ms / args.count
        tree_mean = tree_ms / args.count
        ratio = tree_mean / coordinate_mean if coordinate_mean else 0.0
        print(
            f"{length}\t{args.count}\t{args.repetitions}\t{tails}\t"
            f"{exact_mean:.9f}\t{coordinate_mean:.9f}\t{tree_mean:.9f}\t"
            f"{ratio:.6f}\t{coordinate_visits}\t{tree_visits}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
