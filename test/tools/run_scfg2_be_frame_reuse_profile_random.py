#!/usr/bin/env python3
"""Profile exact BE_WIP_WIP core/right-tail/frame reuse on random cases."""

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
        "--tool", default="build/scfg2_be_frame_reuse_profile_test"
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
        "n\tcases\tactive_cases\tsources\tcores\tright_tails\tframes"
        "\treused_cores\treused_right_tails\tmax_sources_per_core"
        "\tmax_right_tails_per_core\tmax_sources_per_right_tail"
        "\tmax_cores_per_column\tmax_right_tails_per_column"
        "\tmax_frames_per_column\trecursive_sources\tshared_cores"
        "\treused_shared_cores\tmax_recursive_sources_per_core"
        "\tright_shell_sources\tshared_right_tails"
        "\treused_shared_right_tails\tmax_right_shell_sources_per_tail"
        "\tright_only_frames\tleft_right_frames\tfinal_frames"
        "\tmulti_frame_tails\tmax_frames_per_tail"
        "\tmax_final_frames_per_column"
        "\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        sums = {key: 0 for key in (
            "sources", "cores", "right_tails", "frames", "reused_cores",
            "reused_right_tails",
            "recursive_sources", "shared_cores", "reused_shared_cores",
            "right_shell_sources", "shared_right_tails",
            "reused_shared_right_tails",
            "right_only_frames", "left_right_frames", "final_frames",
            "multi_frame_tails",
        )}
        maxima = {key: 0 for key in (
            "max_sources_per_core", "max_right_tails_per_core",
            "max_sources_per_right_tail", "max_cores_per_column",
            "max_right_tails_per_column", "max_frames_per_column",
            "max_recursive_sources_per_core",
            "max_right_shell_sources_per_tail",
            "max_frames_per_tail", "max_final_frames_per_column",
        )}
        active_cases = 0
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            failures += int(row["parity"]) != 1
            active_cases += int(row["sources"]) > 0
            for key in sums:
                sums[key] += int(row[key])
            for key in maxima:
                maxima[key] = max(maxima[key], int(row[key]))
        print(
            f"{length}\t{args.count}\t{active_cases}\t{sums['sources']}\t"
            f"{sums['cores']}\t{sums['right_tails']}\t{sums['frames']}\t"
            f"{sums['reused_cores']}\t{sums['reused_right_tails']}\t"
            f"{maxima['max_sources_per_core']}\t"
            f"{maxima['max_right_tails_per_core']}\t"
            f"{maxima['max_sources_per_right_tail']}\t"
            f"{maxima['max_cores_per_column']}\t"
            f"{maxima['max_right_tails_per_column']}\t"
            f"{maxima['max_frames_per_column']}\t"
            f"{sums['recursive_sources']}\t{sums['shared_cores']}\t"
            f"{sums['reused_shared_cores']}\t"
            f"{maxima['max_recursive_sources_per_core']}\t"
            f"{sums['right_shell_sources']}\t{sums['shared_right_tails']}\t"
            f"{sums['reused_shared_right_tails']}\t"
            f"{maxima['max_right_shell_sources_per_tail']}\t"
            f"{sums['right_only_frames']}\t{sums['left_right_frames']}\t"
            f"{sums['final_frames']}\t{sums['multi_frame_tails']}\t"
            f"{maxima['max_frames_per_tail']}\t"
            f"{maxima['max_final_frames_per_column']}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
