#!/usr/bin/env python3
"""Check planned LTAIL KeepAll parity on reproducible random cases."""

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
    parser.add_argument("--tool", default="build/scfg2_vp_case3_tail_cache_test")
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)

    print(
        "n\tcases\ttails\tsources\tsegments\tevent_pairs\tchild_event_pairs"
        "\tright_event_pairs"
        "\tregistry_event_pairs"
        "\tattachments\tcolumns"
        "\tnoncommittable"
    )
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {key: 0 for key in (
            "tails", "sources", "segments", "event_pairs",
            "child_event_pairs", "right_event_pairs", "registry_event_pairs",
            "attachments", "columns",
            "noncommittable"
        )}
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            if int(row["parity"]) != 1:
                raise RuntimeError(f"parity failure for {case['id']}")
            for key in totals:
                totals[key] += int(row[key])
        print(
            f"{length}\t{args.count}\t{totals['tails']}\t{totals['sources']}\t"
            f"{totals['segments']}\t{totals['event_pairs']}\t"
            f"{totals['child_event_pairs']}\t"
            f"{totals['right_event_pairs']}\t"
            f"{totals['registry_event_pairs']}\t"
            f"{totals['attachments']}\t"
            f"{totals['columns']}\t"
            f"{totals['noncommittable']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
