#!/usr/bin/env python3
"""Observe finite LogMax roots across widths on random scaffold inputs."""

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
    parser.add_argument(
        "--tool", default="build/scfg2_logmax_width_profile_tool")
    parser.add_argument(
        "--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260719)
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("300,400"))
    parser.add_argument(
        "--widths", type=comma_ints, default=comma_ints("20,40,80"))
    parser.add_argument(
        "--emit-failures", action="store_true",
        help="write nonmonotone or Viterbi-bound counterexamples to stderr")
    parser.add_argument(
        "--diagnose-bound-failures", action="store_true",
        help="run the selected-DAG divergence tool when logV exceeds logZ")
    parser.add_argument(
        "--divergence-tool",
        default="build/scfg2_semiring_dag_divergence_tool")
    args = parser.parse_args(argv)
    width_text = ",".join(map(str, args.widths))

    print("n\twidth\tcases\tinside_finite_rate\tviterbi_finite_rate"
          "\tinside_newly_connected\tviterbi_newly_connected"
          "\tinside_monotone_rate\tviterbi_monotone_rate"
          "\tviterbi_le_inside_rate\tlog_inside_min\tlog_viterbi_min"
          "\tcoupled_viterbi_le_inside_rate\tcoupled_inside_matches_rate"
          "\tcoupled_log_probability_max"
          "\tinside_delta_prev_mean\tviterbi_delta_prev_mean"
          "\tlog_probability_max"
          "\telapsed_mean_ms\tinside_item_visits_mean"
          "\tviterbi_item_visits_mean")
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"logmax{length}",
        ])
        totals = defaultdict(lambda: {
            "inside_finite": 0, "viterbi_finite": 0,
            "inside_new": 0, "viterbi_new": 0,
            "inside_monotone": 0, "viterbi_monotone": 0,
            "bounded": 0, "coupled_bounded": 0, "coupled_matches": 0,
            "inside_roots": [], "viterbi_roots": [],
            "inside_deltas": [], "viterbi_deltas": [], "log_probs": [],
            "coupled_log_probs": [],
            "elapsed": 0.0, "inside_visits": 0, "viterbi_visits": 0,
        })
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([
                args.tool, case["seq"], case["G"], width_text,
                args.param_file,
            ])
            previous_inside = None
            previous_viterbi = None
            for row in csv.DictReader(io.StringIO(output), delimiter="\t"):
                width = int(row["width"])
                out = totals[width]
                inside_finite = int(row["inside_finite"])
                viterbi_finite = int(row["viterbi_finite"])
                inside_root = float(row["log_inside_root"])
                viterbi_root = float(row["log_viterbi_root"])
                out["inside_finite"] += inside_finite
                out["viterbi_finite"] += viterbi_finite
                out["inside_monotone"] += int(row["inside_monotone"])
                out["viterbi_monotone"] += int(row["viterbi_monotone"])
                out["bounded"] += int(row["viterbi_le_inside"])
                out["coupled_bounded"] += int(
                    row["coupled_viterbi_le_inside"])
                out["coupled_matches"] += int(row["coupled_inside_matches"])
                coupled_inside = float(row["coupled_log_inside_root"])
                coupled_viterbi = float(row["coupled_log_viterbi_root"])
                if math.isfinite(coupled_inside) and math.isfinite(
                        coupled_viterbi):
                    out["coupled_log_probs"].append(
                        coupled_viterbi - coupled_inside)
                if args.emit_failures and (
                        row["inside_monotone"] == "0" or
                        row["viterbi_monotone"] == "0" or
                        row["viterbi_le_inside"] == "0"):
                    print(
                        "counterexample"
                        f" case_id={case['case_id']} n={length} width={width}"
                        f" log_inside={inside_root:.17g}"
                        f" log_viterbi={viterbi_root:.17g}"
                        f" coupled_log_inside="
                        f"{float(row['coupled_log_inside_root']):.17g}"
                        f" coupled_log_viterbi="
                        f"{float(row['coupled_log_viterbi_root']):.17g}"
                        f" coupled_bounded="
                        f"{row['coupled_viterbi_le_inside']}"
                        f" coupled_inside_matches="
                        f"{row['coupled_inside_matches']}"
                        f" inside_monotone={row['inside_monotone']}"
                        f" viterbi_monotone={row['viterbi_monotone']}"
                        f" viterbi_le_inside={row['viterbi_le_inside']}",
                        file=sys.stderr)
                if (args.diagnose_bound_failures and
                        row["viterbi_le_inside"] == "0"):
                    diagnosis = run([
                        args.divergence_tool, case["seq"], case["G"],
                        str(width), args.param_file,
                    ])
                    print(diagnosis, file=sys.stderr, end="")
                if inside_finite:
                    out["inside_roots"].append(inside_root)
                if viterbi_finite:
                    out["viterbi_roots"].append(viterbi_root)
                if inside_finite and previous_inside is None:
                    out["inside_new"] += 1
                if viterbi_finite and previous_viterbi is None:
                    out["viterbi_new"] += 1
                if inside_finite and previous_inside is not None:
                    out["inside_deltas"].append(inside_root - previous_inside)
                if viterbi_finite and previous_viterbi is not None:
                    out["viterbi_deltas"].append(
                        viterbi_root - previous_viterbi)
                if inside_finite and viterbi_finite:
                    out["log_probs"].append(viterbi_root - inside_root)
                out["elapsed"] += float(row["elapsed_ms"])
                out["inside_visits"] += int(row["inside_item_visits"])
                out["viterbi_visits"] += int(row["viterbi_item_visits"])
                if inside_finite:
                    previous_inside = inside_root
                if viterbi_finite:
                    previous_viterbi = viterbi_root
        for width in args.widths:
            out = totals[width]
            inside_min = min(out["inside_roots"], default=-math.inf)
            viterbi_min = min(out["viterbi_roots"], default=-math.inf)
            inside_delta_mean = (
                sum(out["inside_deltas"]) / len(out["inside_deltas"])
                if out["inside_deltas"] else math.nan)
            viterbi_delta_mean = (
                sum(out["viterbi_deltas"]) / len(out["viterbi_deltas"])
                if out["viterbi_deltas"] else math.nan)
            log_probability_max = max(out["log_probs"], default=-math.inf)
            coupled_log_probability_max = max(
                out["coupled_log_probs"], default=-math.inf)
            print(
                f"{length}\t{width}\t{args.count}\t"
                f"{out['inside_finite']/args.count:.9g}\t"
                f"{out['viterbi_finite']/args.count:.9g}\t"
                f"{out['inside_new']}\t{out['viterbi_new']}\t"
                f"{out['inside_monotone']/args.count:.9g}\t"
                f"{out['viterbi_monotone']/args.count:.9g}\t"
                f"{out['bounded']/args.count:.9g}\t"
                f"{inside_min:.17g}\t{viterbi_min:.17g}\t"
                f"{out['coupled_bounded']/args.count:.9g}\t"
                f"{out['coupled_matches']/args.count:.9g}\t"
                f"{coupled_log_probability_max:.17g}\t"
                f"{inside_delta_mean:.9g}\t{viterbi_delta_mean:.9g}\t"
                f"{log_probability_max:.17g}\t"
                f"{out['elapsed']/args.count:.9g}\t"
                f"{out['inside_visits']/args.count:.9g}\t"
                f"{out['viterbi_visits']/args.count:.9g}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
