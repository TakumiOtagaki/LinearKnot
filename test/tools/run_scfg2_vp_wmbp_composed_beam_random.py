#!/usr/bin/env python3
"""Check composed VP/WMBP child-driven beams on random cases."""

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
        "--tool", default="build/scfg2_vp_wmbp_composed_beam_test"
    )
    parser.add_argument("--generator", default="test/tools/gen_random_seq_struct.py")
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--count", type=int, default=12)
    parser.add_argument(
        "--lengths", type=comma_ints, default=comma_ints("20,30,50,80,100")
    )
    args = parser.parse_args(argv)
    fields = (
        "vp_candidates", "vp_selected", "wmbp_candidates", "wmbp_selected",
        "rejected_vp_sources", "rejected_wmbp_sources",
        "driven_vp_sources", "driven_wmbp_sources",
        "wmbp_parent_visits", "wmbp_source_parents",
        "wmbp_split_source_parents", "wmbp_direct_source_parents",
        "wmbp_direct_only_parents",
        "sparse_wmbp_parent_visits",
        "sparse_vp_parent_visits",
        "sparse_be_parent_visits",
        "keepall_be_parent_visits",
        "v_parent_visits",
        "reverse_wmbp_parent_keys", "reverse_parent_vp_probes",
        "case3_nonleft_wi_queries", "case3_nonleft_geometry_probes",
        "case3_nonleft_right_last_probes",
        "case3_nonleft_middle_last_probes",
        "case3_nonleft_zero_child_probes",
        "case3_nonleft_emitted_sources",
        "case3_nonleft_prefix_candidates",
        "case3_nonleft_active_prefix_probes",
        "case3_nonleft_max_active_prefixes",
    )
    print("n\tcases\t" + "\t".join(fields) + "\tfailures")
    for length_index, length in enumerate(args.lengths):
        generated = run([
            sys.executable, args.generator, "--seed",
            str(args.seed + length_index), "--count", str(args.count),
            "--length", str(length), "--prefix", f"n{length}",
        ])
        totals = {field: 0 for field in fields}
        failures = 0
        for case in csv.DictReader(io.StringIO(generated), delimiter="\t"):
            output = run([args.tool, case["seq"], case["G"], args.param_file])
            row = next(csv.DictReader(io.StringIO(output), delimiter="\t"))
            for field in fields:
                if field == "case3_nonleft_max_active_prefixes":
                    totals[field] = max(totals[field], int(row[field]))
                else:
                    totals[field] += int(row[field])
            failures += int(row["monotone"]) != 1
        print(
            f"{length}\t{args.count}\t"
            + "\t".join(str(totals[field]) for field in fields)
            + f"\t{failures}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
