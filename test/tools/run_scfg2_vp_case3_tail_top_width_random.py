#!/usr/bin/env python3
"""Check finite C3_MR/LTAIL width and monotonicity on random cases."""

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
        "--tool", default="build/scfg2_vp_case3_tail_top_width_test"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\ttail_candidates\ttail_selected\ttail_pruned"
        "\tleft_candidates\trejected_sources\tright_views"
        "\tnonempty_right_views\tempty_right_views\treused_right_views"
        "\tmax_tails_per_right\tright_columns\tmulti_right_columns"
        "\tmax_right_views_per_column\tmiddle_views\tnonempty_middle_views"
        "\tempty_middle_views\treused_middle_views\tmax_tails_per_middle"
        "\tmiddle_columns\tmulti_middle_columns"
        "\tmax_middle_views_per_column\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {key: 0 for key in (
            "tail_candidates", "tail_selected", "tail_pruned",
            "left_candidates", "rejected_sources", "right_views",
            "nonempty_right_views", "empty_right_views", "reused_right_views",
            "right_columns", "multi_right_columns",
            "middle_views", "nonempty_middle_views", "empty_middle_views",
            "reused_middle_views", "middle_columns", "multi_middle_columns",
        )}
        maximum_tails_per_right = 0
        maximum_right_views_per_column = 0
        maximum_tails_per_middle = 0
        maximum_middle_views_per_column = 0
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            failures += int(row["monotone"]) != 1
            for key in totals:
                totals[key] += int(row[key])
            maximum_tails_per_right = max(
                maximum_tails_per_right, int(row["max_tails_per_right"])
            )
            maximum_right_views_per_column = max(
                maximum_right_views_per_column,
                int(row["max_right_views_per_column"]),
            )
            maximum_tails_per_middle = max(
                maximum_tails_per_middle, int(row["max_tails_per_middle"])
            )
            maximum_middle_views_per_column = max(
                maximum_middle_views_per_column,
                int(row["max_middle_views_per_column"]),
            )
        print(
            f"{length}\t{args.count}\t{totals['tail_candidates']}\t"
            f"{totals['tail_selected']}\t{totals['tail_pruned']}\t"
            f"{totals['left_candidates']}\t{totals['rejected_sources']}\t"
            f"{totals['right_views']}\t{totals['nonempty_right_views']}\t"
            f"{totals['empty_right_views']}\t{totals['reused_right_views']}\t"
            f"{maximum_tails_per_right}\t{totals['right_columns']}\t"
            f"{totals['multi_right_columns']}\t"
            f"{maximum_right_views_per_column}\t{totals['middle_views']}\t"
            f"{totals['nonempty_middle_views']}\t"
            f"{totals['empty_middle_views']}\t"
            f"{totals['reused_middle_views']}\t{maximum_tails_per_middle}\t"
            f"{totals['middle_columns']}\t{totals['multi_middle_columns']}\t"
            f"{maximum_middle_views_per_column}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
