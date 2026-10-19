#!/usr/bin/env python3
"""Check width-1 C3_LEFT support and monotonicity on random cases."""

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
        "--tool", default="build/scfg2_vp_case3_left_top_width_test"
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
        "n\tcases\tcandidates\tselected\tmaterialized\tpruned"
        "\trejected_sources\tproducer_candidates\tproducer_selected"
        "\tproducer_pruned\tproducer_rejected_sources"
        "\tproducer_endpoint_pairs\tproducer_sources\tproducer_tails"
        "\tproducer_left_candidates\tcombined_producer_candidates"
        "\tcombined_producer_selected\tcombined_producer_pruned"
        "\tcombined_endpoint_pairs\tcombined_sources"
        "\tcombined_materialized_tails\tcombined_tail_candidates"
        "\tcombined_tail_selected\tcombined_tail_pruned"
        "\tcombined_left_candidates\tcombined_left_selected"
        "\tcombined_left_pruned\tcombined_rejected_sources\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {key: 0 for key in (
            "candidates", "selected", "materialized", "pruned",
            "rejected_sources", "producer_candidates", "producer_selected",
            "producer_pruned", "producer_rejected_sources",
            "producer_endpoint_pairs", "producer_sources", "producer_tails",
            "producer_left_candidates", "combined_producer_candidates",
            "combined_producer_selected", "combined_producer_pruned",
            "combined_endpoint_pairs", "combined_sources",
            "combined_materialized_tails", "combined_tail_candidates",
            "combined_tail_selected", "combined_tail_pruned",
            "combined_left_candidates", "combined_left_selected",
            "combined_left_pruned", "combined_rejected_sources"
        )}
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            failures += int(row["monotone"]) != 1
            for key in totals:
                totals[key] += int(row[key])
        print(
            f"{length}\t{args.count}\t{totals['candidates']}\t"
            f"{totals['selected']}\t{totals['materialized']}\t"
            f"{totals['pruned']}\t"
            f"{totals['rejected_sources']}\t"
            f"{totals['producer_candidates']}\t"
            f"{totals['producer_selected']}\t"
            f"{totals['producer_pruned']}\t"
            f"{totals['producer_rejected_sources']}\t"
            f"{totals['producer_endpoint_pairs']}\t"
            f"{totals['producer_sources']}\t"
            f"{totals['producer_tails']}\t"
            f"{totals['producer_left_candidates']}\t"
            f"{totals['combined_producer_candidates']}\t"
            f"{totals['combined_producer_selected']}\t"
            f"{totals['combined_producer_pruned']}\t"
            f"{totals['combined_endpoint_pairs']}\t"
            f"{totals['combined_sources']}\t"
            f"{totals['combined_materialized_tails']}\t"
            f"{totals['combined_tail_candidates']}\t"
            f"{totals['combined_tail_selected']}\t"
            f"{totals['combined_tail_pruned']}\t"
            f"{totals['combined_left_candidates']}\t"
            f"{totals['combined_left_selected']}\t"
            f"{totals['combined_left_pruned']}\t"
            f"{totals['combined_rejected_sources']}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
