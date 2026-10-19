#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate random seeded legacy cases and export one replay-artifact "
            "TSV row per successful legacy run."
        )
    )
    parser.add_argument("legacy_bin")
    parser.add_argument("param_file")
    parser.add_argument("seed", type=int)
    parser.add_argument("count", type=int)
    parser.add_argument("lengths", help="comma-separated lengths, e.g. 30,50,80")
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def export_legacy_replay_artifact(root: Path,
                                  legacy_bin: str,
                                  param_file: str,
                                  seq: str,
                                  input_g: str,
                                  case_id: str) -> dict[str, str]:
    exporter = root / "test" / "tools" / "export_legacy_replay_artifact.py"
    run = subprocess.run(
        [str(exporter), legacy_bin, param_file, seq, input_g, "--case-id", case_id, "--format", "kv"],
        text=True,
        capture_output=True,
    )
    if run.returncode != 0:
        return {}
    fields: dict[str, str] = {}
    for line in run.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        fields[key] = value
    return fields


def generate_cases(root: Path, seed: int, count: int, length: int) -> list[tuple[str, str, str]]:
    generator = root / "test" / "tools" / "gen_random_seq_struct.py"
    out = subprocess.check_output(
        [
            sys.executable,
            str(generator),
            "--seed",
            str(seed),
            "--count",
            str(count),
            "--length",
            str(length),
            "--prefix",
            f"rand{length}",
        ],
        text=True,
    )
    cases: list[tuple[str, str, str]] = []
    for line in out.splitlines()[1:]:
        if not line.strip():
            continue
        case_id, seq, g = line.split("\t")
        cases.append((case_id, seq, g))
    return cases


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    lengths = [int(part) for part in args.lengths.split(",") if part]
    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    skipped = 0
    with tempfile.TemporaryDirectory() as _tmpdir:
        for length in lengths:
            for case_id, seq, g in generate_cases(root, args.seed, args.count, length):
                fields = export_legacy_replay_artifact(root, args.legacy_bin, args.param_file, seq, g, case_id)
                if not fields:
                    skipped += 1
                    continue
                fields["length"] = str(length)
                rows.append(fields)

    if not rows:
        print("artifact_rows=0")
        print(f"artifact_skipped={skipped}")
        return 0

    header = [
        "case_id",
        "length",
        "seq",
        "input_g",
        "legacy_reported_seed_g",
        "legacy_reported_mfe_r",
        "legacy_reported_mfe_energy",
        "legacy_reported_pf_r",
        "legacy_reported_pf_energy",
        "emitted_lane_g",
        "emitted_lane_gprime",
        "printed_round_projection",
        "printed_square_projection",
        "printed_split_matches_emitted",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"artifact_rows={len(rows)}")
    print(f"artifact_skipped={skipped}")
    print(f"artifact_output_tsv={output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
