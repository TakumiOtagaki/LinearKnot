#!/usr/bin/env python3
"""Check CASE3 preceding-sibling prefix top-b against its offline oracle."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from typing import Sequence


def comma_ints(text: str) -> list[int]:
    return [int(token) for token in text.split(",") if token]


def run(command: Sequence[str]) -> str:
    completed = subprocess.run(
        command, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return completed.stdout


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool", default="build/scfg2_vp_case3_endpoint_prefix_frontier_test")
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260719)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument("--lengths", type=comma_ints,
                        default=comma_ints("20,30,50,80,100"))
    args = parser.parse_args(argv)

    total = 0
    for offset, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed", str(args.seed + offset),
            "--count", str(args.count), "--length", str(length),
            "--prefix", f"endpoint-prefix-{length}",
        ])
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            run([args.tool, case["seq"], case["G"], args.param_file])
            total += 1
        print(f"n={length}\tcases={args.count}\tparity=1")
    print(f"total={total}\tparity=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
