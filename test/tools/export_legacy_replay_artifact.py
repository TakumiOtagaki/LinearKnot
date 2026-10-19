#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run legacy CParty once and emit the canonical replay artifact "
            "for that seeded-folding case: original seed G, legacy-emitted "
            "lane split (G, G'), printed mixed structure R, and reported MFE energy."
        )
    )
    parser.add_argument("legacy_bin")
    parser.add_argument("param_file")
    parser.add_argument("seq")
    parser.add_argument("input_g")
    parser.add_argument("--case-id", default="single_case")
    parser.add_argument(
        "--format",
        choices=("tsv", "kv"),
        default="tsv",
        help="Output format. 'tsv' prints a header plus one row; 'kv' prints key=value lines.",
    )
    return parser.parse_args()


def parse_cparty_stdout(parser_script: Path, stdout_text: str) -> dict[str, str]:
    run = subprocess.run(
        [str(parser_script)],
        input=stdout_text,
        text=True,
        capture_output=True,
        check=True,
    )
    fields = run.stdout.rstrip("\n").split("\t")
    if len(fields) != 6:
        raise RuntimeError(f"unexpected parse_cparty_stdout output: {run.stdout!r}")
    return {
        "seq": fields[0],
        "restricted": fields[1],
        "mfe_structure": fields[2],
        "mfe_energy": fields[3],
        "pf_structure": fields[4],
        "pf_energy": fields[5],
    }


def parse_legacy_lanes(stderr_text: str) -> tuple[str, str]:
    for line in stderr_text.splitlines():
        if not line.startswith("[legacy-lanes] "):
            continue
        payload = line.split(" g=", 1)
        if len(payload) != 2 or " gprime=" not in payload[1]:
            continue
        lane_g, lane_gprime = payload[1].split(" gprime=", 1)
        return lane_g, lane_gprime
    raise RuntimeError("legacy lanes were not emitted")


def round_projection(structure: str) -> str:
    return structure.replace("[", ".").replace("]", ".")


def square_projection(structure: str) -> str:
    return structure.replace("(", ".").replace(")", ".")


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    parser_script = root / "test" / "tools" / "parse_cparty_stdout.sh"

    run = subprocess.run(
        [args.legacy_bin, "-d2", "-r", args.input_g, "-P", args.param_file, args.seq],
        text=True,
        capture_output=True,
        env={**os.environ, "CPARTY_LEGACY_EMIT_LANES": "1"},
    )
    if run.returncode != 0:
        sys.stderr.write(run.stderr or run.stdout or f"legacy exited with code {run.returncode}\n")
        return run.returncode or 1

    parsed = parse_cparty_stdout(parser_script, run.stdout)
    emitted_g, emitted_gprime_round = parse_legacy_lanes(run.stderr)
    emitted_gprime = emitted_gprime_round.replace("(", "[").replace(")", "]")
    printed_r = parsed["mfe_structure"]
    printed_round_projection = round_projection(printed_r)
    printed_square_projection = square_projection(printed_r)

    fields = {
        "case_id": args.case_id,
        "seq": args.seq,
        "input_g": args.input_g,
        "legacy_reported_seed_g": parsed["restricted"],
        "legacy_reported_mfe_r": printed_r,
        "legacy_reported_mfe_energy": parsed["mfe_energy"],
        "legacy_reported_pf_r": parsed["pf_structure"],
        "legacy_reported_pf_energy": parsed["pf_energy"],
        "emitted_lane_g": emitted_g,
        "emitted_lane_gprime": emitted_gprime,
        "printed_round_projection": printed_round_projection,
        "printed_square_projection": printed_square_projection,
        "printed_split_matches_emitted": "1"
        if printed_round_projection == emitted_g and printed_square_projection == emitted_gprime
        else "0",
    }

    if args.format == "kv":
        for key, value in fields.items():
            print(f"{key}={value}")
        return 0

    header = list(fields.keys())
    print("\t".join(header))
    print("\t".join(fields[key] for key in header))
    return 0


if __name__ == "__main__":
    sys.exit(main())
