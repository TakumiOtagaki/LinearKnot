#!/usr/bin/env python3
"""Profile oracle-free online CASE3 retain-all versus producer width 1."""

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


def profile(tool: str, seq: str, scaffold: str, width: str,
            repetitions: int, mode: str, param_file: str,
            scope: str, tail_width: str = "all",
            left_width: str = "all") -> Dict[str, str]:
    output = run([
        tool, seq, scaffold, width, str(repetitions), mode, param_file, scope,
        tail_width, left_width,
    ])
    return next(csv.DictReader(io.StringIO(output), delimiter="\t"))


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool", default="build/scfg2_vp_case3_online_profile_tool"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--mode", choices=("inside", "max"), default="inside")
    parser.add_argument(
        "--scope", choices=("eligible", "all-wi"), default="eligible"
    )
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\tmode\tscope\tall_ms\twidth1_ms\ttime_ratio"
        "\tall_endpoint_pairs\twidth1_endpoint_pairs"
        "\tall_sources\twidth1_sources\tall_tails\twidth1_tails"
        "\tall_left\twidth1_left\tall_producers\twidth1_selected_producers"
        "\tcombined_ms\tcombined_time_ratio\tcombined_tails"
        "\tcombined_tails_selected\tcombined_left"
        "\tcombined_left_selected\tcombined_rejected_sources"
        "\tcombined_active_expiry_probes"
        "\tcombined_tail_materialization_probes"
        "\tcombined_prepared_score_probes"
        "\tcombined_max_active_segments"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        all_ms = 0.0
        width1_ms = 0.0
        combined_ms = 0.0
        totals = {key: 0 for key in (
            "all_endpoint_pairs", "width1_endpoint_pairs",
            "all_sources", "width1_sources", "all_tails", "width1_tails",
            "all_left", "width1_left", "all_producers",
            "width1_selected_producers",
            "combined_tails", "combined_tails_selected", "combined_left",
            "combined_left_selected", "combined_rejected_sources",
            "combined_active_expiry_probes",
            "combined_tail_materialization_probes",
            "combined_prepared_score_probes",
            "combined_max_active_segments",
        )}
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            retain_all = profile(
                args.tool, case["seq"], case["G"], "all",
                args.repetitions, args.mode, args.param_file, args.scope,
            )
            width_one = profile(
                args.tool, case["seq"], case["G"], "1",
                args.repetitions, args.mode, args.param_file, args.scope,
            )
            combined_one = profile(
                args.tool, case["seq"], case["G"], "1",
                args.repetitions, args.mode, args.param_file, args.scope,
                "1", "1",
            )
            all_ms += float(retain_all["dp_ms_per_run"])
            width1_ms += float(width_one["dp_ms_per_run"])
            combined_ms += float(combined_one["dp_ms_per_run"])
            for prefix, row in (("all", retain_all), ("width1", width_one)):
                totals[f"{prefix}_endpoint_pairs"] += int(row["endpoint_pairs"])
                totals[f"{prefix}_sources"] += int(row["sources"])
                totals[f"{prefix}_tails"] += int(row["tails"])
                totals[f"{prefix}_left"] += int(row["left_candidates"])
            totals["all_producers"] += int(retain_all["producer_candidates"])
            totals["width1_selected_producers"] += int(
                width_one["producer_selected"]
            )
            totals["combined_tails"] += int(combined_one["tails"])
            totals["combined_tails_selected"] += int(
                combined_one["tails_selected"]
            )
            totals["combined_left"] += int(combined_one["left_candidates"])
            totals["combined_left_selected"] += int(
                combined_one["left_selected"]
            )
            totals["combined_rejected_sources"] += int(
                combined_one["rejected_sources"]
            )
            totals["combined_active_expiry_probes"] += int(
                combined_one["active_expiry_probes"]
            )
            totals["combined_tail_materialization_probes"] += int(
                combined_one["tail_materialization_probes"]
            )
            totals["combined_prepared_score_probes"] += int(
                combined_one["prepared_score_probes"]
            )
            totals["combined_max_active_segments"] = max(
                totals["combined_max_active_segments"],
                int(combined_one["max_active_segments"]),
            )
        all_mean = all_ms / args.count
        width1_mean = width1_ms / args.count
        ratio = width1_mean / all_mean if all_mean else 0.0
        combined_mean = combined_ms / args.count
        combined_ratio = combined_mean / all_mean if all_mean else 0.0
        print(
            f"{length}\t{args.count}\t{args.mode}\t{args.scope}\t"
            f"{all_mean:.6f}\t"
            f"{width1_mean:.6f}\t{ratio:.6f}\t"
            f"{totals['all_endpoint_pairs']}\t{totals['width1_endpoint_pairs']}\t"
            f"{totals['all_sources']}\t{totals['width1_sources']}\t"
            f"{totals['all_tails']}\t{totals['width1_tails']}\t"
            f"{totals['all_left']}\t{totals['width1_left']}\t"
            f"{totals['all_producers']}\t"
            f"{totals['width1_selected_producers']}\t"
            f"{combined_mean:.6f}\t{combined_ratio:.6f}\t"
            f"{totals['combined_tails']}\t"
            f"{totals['combined_tails_selected']}\t"
            f"{totals['combined_left']}\t"
            f"{totals['combined_left_selected']}\t"
            f"{totals['combined_rejected_sources']}\t"
            f"{totals['combined_active_expiry_probes']}\t"
            f"{totals['combined_tail_materialization_probes']}\t"
            f"{totals['combined_prepared_score_probes']}\t"
            f"{totals['combined_max_active_segments']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
