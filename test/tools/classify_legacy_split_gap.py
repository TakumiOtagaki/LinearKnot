#!/usr/bin/env python3
"""Classify how printed-union naive split differs from fresh legacy exact lanes."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare legacy exact lane export against naive printed-union split "
            "and input_g to classify remaining merge-less gaps."
        )
    )
    parser.add_argument("legacy_bin")
    parser.add_argument("cases_tsv")
    parser.add_argument("param_file")
    parser.add_argument("--output-tsv", default="")
    return parser.parse_args()


def pair_set(db: str) -> set[tuple[int, int]]:
    stack: list[int] = []
    out: set[tuple[int, int]] = set()
    for idx, ch in enumerate(db, start=1):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            out.add((stack.pop(), idx))
    return out


def round_lane(structure: str) -> str:
    return structure.replace("[", ".").replace("]", ".")


def square_lane(structure: str) -> str:
    return structure.replace("(", ".").replace(")", ".")


def export_legacy_replay_artifact(root: Path,
                                  legacy_bin: str,
                                  seq: str,
                                  input_g: str,
                                  param_file: str) -> dict[str, str]:
    exporter = root / "test" / "tools" / "export_legacy_replay_artifact.py"
    run = subprocess.run(
        [str(exporter), legacy_bin, param_file, seq, input_g, "--format", "kv"],
        text=True,
        capture_output=True,
        check=True,
    )
    fields: dict[str, str] = {}
    for line in run.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        fields[key] = value
    return fields


def classify(row: dict[str, str], exact_g: str, exact_gp: str) -> dict[str, str]:
    naive_g = round_lane(row["legacy_mfe_structure"])
    naive_gp = square_lane(row["legacy_mfe_structure"])
    input_g = row["input_g"]

    exact_pairs = pair_set(exact_g)
    naive_pairs = pair_set(naive_g)
    input_pairs = pair_set(input_g)

    naive_missing = sorted(exact_pairs - naive_pairs)
    naive_extra = sorted(naive_pairs - exact_pairs)
    input_missing = sorted(exact_pairs - input_pairs)
    input_extra = sorted(input_pairs - exact_pairs)

    if not naive_missing and not naive_extra:
        gap_kind = "naive_exact"
    elif naive_missing and not naive_extra and not input_missing and not input_extra:
        gap_kind = "seed_block_missing_from_union"
    elif not naive_missing and not naive_extra and (input_missing or input_extra):
        gap_kind = "legacy_added_round_pairs"
    elif naive_missing and not naive_extra and input_missing:
        gap_kind = "seed_missing_and_legacy_added"
    else:
        gap_kind = "mixed_other"

    return {
        "case_id": row["case_id"],
        "length": row["length"],
        "gap_kind": gap_kind,
        "naive_matches_exact_g": "1" if naive_g == exact_g else "0",
        "input_matches_exact_g": "1" if input_g == exact_g else "0",
        "naive_matches_exact_gp": "1" if naive_gp.replace("[", "(").replace("]", ")") == exact_gp else "0",
        "naive_missing_pairs": ",".join(f"{i}-{j}" for i, j in naive_missing),
        "naive_extra_pairs": ",".join(f"{i}-{j}" for i, j in naive_extra),
        "input_missing_pairs": ",".join(f"{i}-{j}" for i, j in input_missing),
        "input_extra_pairs": ",".join(f"{i}-{j}" for i, j in input_extra),
    }


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    rows = []
    with Path(args.cases_tsv).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            artifact = export_legacy_replay_artifact(root, args.legacy_bin, row["seq"], row["input_g"], args.param_file)
            exact_g = artifact["emitted_lane_g"]
            exact_gp = artifact["emitted_lane_gprime"].replace("[", "(").replace("]", ")")
            rows.append(classify(row, exact_g, exact_gp))

    if args.output_tsv:
        output = Path(args.output_tsv)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["gap_kind"]] = counts.get(row["gap_kind"], 0) + 1
    print(f"classified_cases={len(rows)}")
    for key in sorted(counts):
        print(f"gap_kind_{key}={counts[key]}")
    for row in rows:
        print(
            "case "
            f"case_id={row['case_id']} "
            f"gap_kind={row['gap_kind']} "
            f"naive_missing={row['naive_missing_pairs'] or 'none'} "
            f"input_missing={row['input_missing_pairs'] or 'none'} "
            f"input_extra={row['input_extra_pairs'] or 'none'}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
