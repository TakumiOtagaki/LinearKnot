#!/usr/bin/env python3
"""Aggregate finite SCFG2 quality against the KeepAll online kernel."""

from __future__ import annotations

import argparse
import csv
import io
import math
import subprocess
import sys
from collections import defaultdict
from typing import Sequence


def run(args: Sequence[str]) -> str:
    completed = subprocess.run(
        args, check=True, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    return completed.stdout


def comma_ints(text: str) -> list[int]:
    return [int(token) for token in text.split(",") if token]


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", default="build/scfg2_approximation_quality_tool")
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--lengths", type=comma_ints,
                        default=comma_ints("20,30,50"))
    parser.add_argument("--widths", type=comma_ints,
                        default=comma_ints("1,2,5,10,20,40"))
    parser.add_argument(
        "--geometry-mode", default="finite",
        choices=("finite", "beam-only", "hairpin-only", "ancestor-only"),
        help="separate beam pruning from fixed finite geometry truncation")
    args = parser.parse_args(argv)
    width_text = ",".join(str(width) for width in args.widths)

    print("n\tgeometry\twidth\tcases\tinside_ratio_geomean\tinside_ratio_min"
          "\tinside_zero_rate\tviterbi_ratio_geomean\tviterbi_ratio_min"
          "\tviterbi_zero_rate\tfull_trace_rate"
          "\titem_coverage\tchoice_coverage\tpair_precision\tpair_recall"
          "\tinside_monotone_rate\tviterbi_monotone_rate")
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"quality{length}",
        ])
        totals = defaultdict(lambda: {
            "inside_logs": [], "viterbi_logs": [],
            "inside_min": 1.0, "viterbi_min": 1.0,
            "inside_zero": 0, "viterbi_zero": 0,
            "full": 0, "items": 0, "retained": 0,
            "inside_monotone": 0, "viterbi_monotone": 0,
            "choices": 0, "matching": 0,
            "exact_pairs": 0, "finite_pairs": 0, "intersection": 0,
        })
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([
                args.tool, case["seq"], case["G"], width_text,
                args.param_file, args.geometry_mode,
            ])
            for row in csv.DictReader(io.StringIO(output), delimiter="\t"):
                width = int(row["width"])
                out = totals[width]
                inside = float(row["inside_ratio"])
                viterbi = float(row["viterbi_ratio"])
                if inside > 0.0:
                    out["inside_logs"].append(math.log(inside))
                else:
                    out["inside_zero"] += 1
                if viterbi > 0.0:
                    out["viterbi_logs"].append(math.log(viterbi))
                else:
                    out["viterbi_zero"] += 1
                out["inside_min"] = min(out["inside_min"], inside)
                out["viterbi_min"] = min(out["viterbi_min"], viterbi)
                out["full"] += int(row["full_trace_match"])
                out["inside_monotone"] += int(row["inside_monotone"])
                out["viterbi_monotone"] += int(row["viterbi_monotone"])
                out["items"] += int(row["trace_items"])
                out["retained"] += int(row["retained_trace_items"])
                out["choices"] += int(row["trace_choices"])
                out["matching"] += int(row["matching_trace_choices"])
                out["exact_pairs"] += int(row["exact_generated_pairs"])
                out["finite_pairs"] += int(row["finite_generated_pairs"])
                out["intersection"] += int(row["pair_intersection"])
        for width in args.widths:
            out = totals[width]
            geomean_inside = 0.0 if out["inside_zero"] else math.exp(
                sum(out["inside_logs"]) / args.count)
            geomean_viterbi = 0.0 if out["viterbi_zero"] else math.exp(
                sum(out["viterbi_logs"]) / args.count)
            ratio = lambda a, b: a / b if b else 1.0
            print(
                f"{length}\t{args.geometry_mode}\t{width}\t{args.count}\t"
                f"{geomean_inside:.9g}\t"
                f"{out['inside_min']:.9g}\t{out['inside_zero']/args.count:.9g}\t"
                f"{geomean_viterbi:.9g}\t{out['viterbi_min']:.9g}\t"
                f"{out['viterbi_zero']/args.count:.9g}\t"
                f"{out['full']/args.count:.9g}\t"
                f"{ratio(out['retained'], out['items']):.9g}\t"
                f"{ratio(out['matching'], out['choices']):.9g}\t"
                f"{ratio(out['intersection'], out['finite_pairs']):.9g}\t"
                f"{ratio(out['intersection'], out['exact_pairs']):.9g}\t"
                f"{out['inside_monotone']/args.count:.9g}\t"
                f"{out['viterbi_monotone']/args.count:.9g}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
