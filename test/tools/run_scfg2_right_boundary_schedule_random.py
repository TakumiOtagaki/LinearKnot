#!/usr/bin/env python3
"""Check exact right-boundary schedule parity on reproducible random cases."""

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
        "--tool", default="build/scfg2_right_boundary_schedule_test"
    )
    parser.add_argument(
        "--generator", default="test/tools/gen_random_seq_struct.py"
    )
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)
    if args.count <= 0:
        parser.error("--count must be positive")

    print(
        "n\tcases\titems_checked\tlegacy_excluded\t"
        "inside_failures\tviterbi_failures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run(
            [
                sys.executable,
                args.generator,
                "--seed",
                str(args.seed + length_index),
                "--count",
                str(args.count),
                "--length",
                str(length),
                "--prefix",
                f"n{length}",
            ]
        )
        items_checked = 0
        legacy_excluded = 0
        inside_failures = 0
        viterbi_failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run(
                [args.tool, case["seq"], case["G"], args.param_file]
            )
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            items_checked += int(row["items"])
            legacy_excluded += int(row["legacy_excluded"])
            inside_failures += int(row["inside_parity"]) != 1
            viterbi_failures += int(row["viterbi_parity"]) != 1
        print(
            f"{length}\t{args.count}\t{items_checked}\t{legacy_excluded}\t"
            f"{inside_failures}\t{viterbi_failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
