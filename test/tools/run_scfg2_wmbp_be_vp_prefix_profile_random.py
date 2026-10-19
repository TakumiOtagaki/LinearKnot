#!/usr/bin/env python3
"""Aggregate exact WMBP BE-prefix / right-VP reuse on random cases."""

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
        "--tool", default="build/scfg2_wmbp_be_vp_prefix_profile_test"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    sum_fields = (
        "sources", "wmbp", "wmbw", "wi_empty", "wi", "prefixes",
        "reused_prefixes", "right_vps", "reused_right_vps",
        "multi_source_prefixes", "readiness_checks",
    )
    print(
        "n\tcases\t" + "\t".join(sum_fields)
        + "\tmax_sources_per_prefix\tmax_prefixes_per_right_vp"
        + "\tlocal_variant_prefixes\tmax_locals_per_prefix"
        + "\tmax_sources_per_column\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {field: 0 for field in sum_fields}
        maxima = {
            "max_sources_per_prefix": 0,
            "max_prefixes_per_right_vp": 0,
            "local_variant_prefixes": 0,
            "max_locals_per_prefix": 0,
            "max_sources_per_column": 0,
        }
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            for field in totals:
                totals[field] += int(row[field])
            for field in maxima:
                maxima[field] = max(maxima[field], int(row[field]))
            failures += int(row["parity"]) != 1
        print(
            f"{length}\t{args.count}\t"
            + "\t".join(str(totals[field]) for field in sum_fields)
            + "\t" + "\t".join(str(maxima[field]) for field in maxima)
            + f"\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
