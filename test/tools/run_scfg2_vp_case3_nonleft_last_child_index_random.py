#!/usr/bin/env python3
"""Check exact VP CASE3 Boundary/Other readiness on random scaffolds."""

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
        "--tool", default="build/scfg2_vp_case3_nonleft_last_child_index_test"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\tboundary\tother\tzero_child\tinside_positive"
        "\tviterbi_positive\tfailures"
    )
    any_failures = False
    for length_index, length in enumerate(args.lengths):
        generated = run([
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
        ])
        totals = {key: 0 for key in (
            "boundary", "other", "zero_child", "inside_positive",
            "viterbi_positive"
        )}
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            try:
                output = run([args.tool, case["seq"], case["G"], args.param_file])
                row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
                for key in totals:
                    totals[key] += int(row[key])
            except (RuntimeError, StopIteration, KeyError, ValueError) as error:
                failures += 1
                any_failures = True
                print(f"case {case['case_id']} failed: {error}", file=sys.stderr)
        print(
            f"{length}\t{args.count}\t{totals['boundary']}\t{totals['other']}\t"
            f"{totals['zero_child']}\t{totals['inside_positive']}\t"
            f"{totals['viterbi_positive']}\t{failures}"
        )
    return 1 if any_failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
