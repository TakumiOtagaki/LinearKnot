#!/usr/bin/env python3
"""Aggregate SCFG2 right-boundary dependency classes on random cases."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from typing import Dict, List, Sequence


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
        "--tool", default="build/scfg2_right_boundary_dependency_tool"
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

    fields = [
        "edges",
        "prior_column",
        "future_column",
        "same_column_inward",
        "same_column_outward",
        "same_span_phase",
    ]
    print("n\tcases\t" + "\t".join(fields) + "\tphase_cycle_cases")
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
        totals: Dict[str, int] = {field: 0 for field in fields}
        phase_cycle_cases = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run(
                [args.tool, case["seq"], case["G"], args.param_file]
            )
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            for field in fields:
                totals[field] += int(row[field])
            if int(row["phase_acyclic"]) != 1:
                phase_cycle_cases += 1
        print(
            f"{length}\t{args.count}\t"
            + "\t".join(str(totals[field]) for field in fields)
            + f"\t{phase_cycle_cases}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
