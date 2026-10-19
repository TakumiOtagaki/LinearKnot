#!/usr/bin/env python3
"""Check VP-before-WMBP phase-barrier feasibility on random cases."""

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
    parser.add_argument("--tool", default="build/scfg2_wmbp_vp_phase_barrier_test")
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    fields = (
        "items", "same_column_edges", "vp_barrier_edges",
        "barrier_wmbp_edges", "cyclic_columns", "missing_children",
    )
    mask_fields = (
        "before_mask", "after_mask", "split_mask", "prerequisite_mask",
        "wmbp_phase_mask", "be_wmbp_phase_mask",
    )
    fields = fields + (
        "vp_prerequisites", "wmbp_prerequisites",
        "wmbp_phase_prerequisites", "wmbp_phase_targets",
        "wmbp_phase_vp", "wmbp_phase_be", "wmbp_phase_wmbp",
        "wmbp_phase_missing",
        "be_wmbp_phase_prerequisites", "be_wmbp_phase_targets",
        "be_wmbp_phase_vp", "be_wmbp_phase_be", "be_wmbp_phase_wmbp",
        "be_wmbp_phase_missing",
    )
    print(
        "n\tcases\t" + "\t".join(fields) + "\t"
        + "\t".join(mask_fields) + "\tfailures"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {field: 0 for field in fields}
        masks = {field: 0 for field in mask_fields}
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            for field in fields:
                totals[field] += int(row[field])
            for field in mask_fields:
                masks[field] |= int(row[field])
            failures += int(row["feasible"]) != 1
        print(
            f"{length}\t{args.count}\t"
            + "\t".join(str(totals[field]) for field in fields)
            + "\t" + "\t".join(str(masks[field]) for field in mask_fields)
            + f"\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
