#!/usr/bin/env python3
"""Check CASE3 tail-column schedule feasibility on reproducible random cases."""

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
        "--tool", default="build/scfg2_vp_case3_tail_schedule_tool"
    )
    parser.add_argument(
        "--generator", default="test/tools/gen_random_seq_struct.py"
    )
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--schedule", choices=("span", "right"), default="span")
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)
    if args.count <= 0:
        parser.error("--count must be positive")

    print(
        "n\tcases\tcases_with_tails\ttails\tattachments\tmaterialized"
        "\treused\tcolumns\ttopological_tails\tcommittable_columns"
        "\tnoncommittable_columns\tmissing_schedule_items\tmax_overlap"
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
        totals: Dict[str, int] = {
            "tails": 0,
            "attachments": 0,
            "materialized": 0,
            "reused": 0,
            "columns": 0,
            "topological_tails": 0,
            "committable_columns": 0,
            "noncommittable_columns": 0,
            "missing_schedule_items": 0,
        }
        cases_with_tails = 0
        max_overlap = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run(
                [args.tool, case["seq"], case["G"], args.param_file, args.schedule]
            )
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            if int(row["tails"]) > 0:
                cases_with_tails += 1
            for field in totals:
                totals[field] += int(row[field])
            max_overlap = max(max_overlap, int(row["max_overlap"]))
        print(
            f"{length}\t{args.count}\t{cases_with_tails}\t{totals['tails']}\t"
            f"{totals['attachments']}\t{totals['materialized']}\t"
            f"{totals['reused']}\t{totals['columns']}\t"
            f"{totals['topological_tails']}\t"
            f"{totals['committable_columns']}\t"
            f"{totals['noncommittable_columns']}\t"
            f"{totals['missing_schedule_items']}\t{max_overlap}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
