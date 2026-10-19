#!/usr/bin/env python3
"""Benchmark coordinate and scaffold-indexed WMBP split scans."""

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
        "--tool", default="build/scfg2_wmbp_be_vp_split_index_profile_tool"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\tcoordinate_visits\tcoordinate_geometry_visits"
        "\tindexed_visits\tbucket_visits"
        "\tmean_build_ms\tmean_coordinate_ms\tmean_indexed_ms"
        "\tindexed_over_coordinate\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        visits = {field: 0 for field in (
            "coordinate_visits", "coordinate_geometry_visits",
            "indexed_visits", "bucket_visits",
        )}
        times = {field: 0.0 for field in (
            "index_build_ms", "coordinate_ms", "indexed_ms",
        )}
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([
                args.tool, case["seq"], case["G"], args.param_file,
                str(args.repeats),
            ])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            for field in visits:
                visits[field] += int(row[field])
            for field in times:
                times[field] += float(row[field])
            failures += int(row["parity"]) != 1
        mean = {field: value / args.count for field, value in times.items()}
        ratio = mean["indexed_ms"] / mean["coordinate_ms"]
        print(
            f"{length}\t{args.count}\t{visits['coordinate_visits']}\t"
            f"{visits['coordinate_geometry_visits']}\t"
            f"{visits['indexed_visits']}\t{visits['bucket_visits']}\t"
            f"{mean['index_build_ms']:.6f}\t{mean['coordinate_ms']:.6f}\t"
            f"{mean['indexed_ms']:.6f}\t{ratio:.6f}\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
