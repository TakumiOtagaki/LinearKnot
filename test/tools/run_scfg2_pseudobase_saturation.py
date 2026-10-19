#!/usr/bin/env python3
"""Run staged SCFG2 saturation ladders on selected PseudoBase++ entries."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from pathlib import Path
from typing import Sequence


DEFAULT_IDS = "PKB396,PKB356,PKB357,PKB354,PKB417"


def comma_strings(text: str) -> list[str]:
    return [token.strip() for token in text.split(",") if token.strip()]


def load_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = {
            row["pkb_short"]: row
            for row in csv.DictReader(handle)
        }
    return rows


def validate_case(row: dict[str, str]) -> None:
    sequence = row["sequence"].upper()
    if row["pipeline_ready"].lower() != "true":
        raise ValueError(f"{row['pkb_short']} is not pipeline-ready")
    invalid = sorted(set(sequence) - set("ACGU"))
    if invalid:
        raise ValueError(
            f"{row['pkb_short']} contains unsupported residues: {invalid}")
    if len(sequence) != int(row["length"]):
        raise ValueError(f"{row['pkb_short']} has inconsistent length")


def scaffold_for(row: dict[str, str], mode: str) -> str:
    if mode == "empty":
        return "empty"
    structure = row["dotbracket"]
    if len(structure) != int(row["length"]):
        raise ValueError(
            f"{row['pkb_short']} has inconsistent reference structure length")
    scaffold = "".join(char if char in "()" else "." for char in structure)
    depth = 0
    for char in scaffold:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise ValueError(
                    f"{row['pkb_short']} has invalid parenthesis projection")
    if depth != 0:
        raise ValueError(
            f"{row['pkb_short']} has invalid parenthesis projection")
    if "(" not in scaffold:
        raise ValueError(
            f"{row['pkb_short']} has an empty parenthesis projection")
    return scaffold


def run_case(
    tool: Path,
    mode: str,
    row: dict[str, str],
    scaffold_mode: str,
    widths: str,
    param_file: Path,
) -> list[dict[str, str]]:
    scaffold = scaffold_for(row, scaffold_mode)
    command = [
        str(tool),
        row["sequence"].upper(),
        scaffold,
        widths,
        str(param_file),
    ]
    if mode == "inside":
        command.append("noexact")
    print(
        f"[pseudobase-saturation] start {row['pkb_short']}"
        f" n={row['length']} mode={mode} scaffold={scaffold_mode}",
        file=sys.stderr,
        flush=True,
    )
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        # Preserve live per-width progress from the C++ measurement tool.
        stderr=None,
    )
    print(
        f"[pseudobase-saturation] done {row['pkb_short']}",
        file=sys.stderr,
        flush=True,
    )
    return list(csv.DictReader(
        io.StringIO(completed.stdout), delimiter="\t"))


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--mode", choices=("inside", "outside"), required=True)
    parser.add_argument("--scaffold", choices=("empty", "round-projection"),
                        default="empty")
    parser.add_argument("--ids", default=DEFAULT_IDS)
    parser.add_argument("--widths", default="50,100,200,400")
    parser.add_argument("--param-file", type=Path,
                        default=Path("params/rna_Turner04.par"))
    parser.add_argument("--inside-tool", type=Path,
                        default=Path("build/scfg2_beam_inside_saturation_tool"))
    parser.add_argument("--outside-tool", type=Path,
                        default=Path("build/scfg2_beam_outside_saturation_tool"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    dataset = load_rows(args.dataset)
    ids = comma_strings(args.ids)
    if not ids:
        parser.error("--ids must contain at least one entry")
    missing = [case_id for case_id in ids if case_id not in dataset]
    if missing:
        parser.error(f"dataset is missing IDs: {','.join(missing)}")

    tool = args.inside_tool if args.mode == "inside" else args.outside_tool
    combined: list[dict[str, str]] = []
    tool_fields: list[str] | None = None
    for case_id in ids:
        row = dataset[case_id]
        validate_case(row)
        measurements = run_case(
            tool, args.mode, row, args.scaffold, args.widths, args.param_file)
        if not measurements:
            raise RuntimeError(f"{case_id} produced no measurements")
        if tool_fields is None:
            tool_fields = list(measurements[0])
        elif list(measurements[0]) != tool_fields:
            raise RuntimeError(f"{case_id} changed the tool output schema")
        for measurement in measurements:
            combined.append({
                "pkb_short": case_id,
                "length": row["length"],
                "rna_type": row["rna_type"],
                "mode": args.mode,
                "scaffold_mode": args.scaffold,
                **measurement,
            })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "pkb_short", "length", "rna_type", "mode", "scaffold_mode",
        *(tool_fields or []),
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(combined)
    print(
        f"[pseudobase-saturation] wrote {len(combined)} rows to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
