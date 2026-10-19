#!/usr/bin/env python3
"""Profile virtual SCFG2 beam buckets over reproducible random (S, G)."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


def comma_ints(text: str) -> List[int]:
    values = sorted({int(token) for token in text.split(",") if token})
    if not values or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return values


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool", default="build/scfg2_virtual_beam_profile_tool"
    )
    parser.add_argument(
        "--generator", default="test/tools/gen_random_seq_struct.py"
    )
    parser.add_argument("--param-file", default="params/rna_Turner04.par")
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument("--count", type=int, default=3, help="cases per length")
    parser.add_argument("--lengths", type=comma_ints, default=comma_ints("20,30,50"))
    parser.add_argument("--widths", type=comma_ints, default=comma_ints("10,20,50,100"))
    parser.add_argument(
        "--raw-output", help="optional path for case-level bucket rows"
    )
    return parser.parse_args(argv)


def run_checked(command: Sequence[str]) -> str:
    completed = subprocess.run(
        command, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr}"
        )
    return completed.stdout


def generate_cases(
    generator: str, seed: int, count: int, length: int
) -> Iterable[Dict[str, str]]:
    output = run_checked(
        [
            sys.executable,
            generator,
            "--seed",
            str(seed),
            "--count",
            str(count),
            "--length",
            str(length),
            "--prefix",
            f"n{length}",
        ]
    )
    yield from csv.DictReader(io.StringIO(output), delimiter="\t")


def profile_case(
    tool: str,
    sequence: str,
    scaffold: str,
    widths: Sequence[int],
    param_file: str,
) -> Iterable[Dict[str, str]]:
    output = run_checked(
        [tool, sequence, scaffold, ",".join(map(str, widths)), param_file]
    )
    yield from csv.DictReader(io.StringIO(output), delimiter="\t")


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.count <= 0:
        print("error: --count must be positive", file=sys.stderr)
        return 2

    tool = str(Path(args.tool))
    generator = str(Path(args.generator))
    raw_handle = None
    raw_writer = None
    if args.raw_output:
        raw_handle = Path(args.raw_output).open("w", encoding="utf-8", newline="")

    # Aggregates are keyed by (length, family, width). Source sums and selected
    # counts are pooled over buckets, rather than averaging per-bucket ratios.
    stats: Dict[Tuple[int, str, int], Dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    cases_with_family: Dict[Tuple[int, str, int], set[str]] = defaultdict(set)
    priority_kinds: Dict[Tuple[int, str, int], set[str]] = defaultdict(set)
    total_cases_by_length: Dict[int, int] = defaultdict(int)
    try:
        for length_index, length in enumerate(args.lengths):
            seed = args.seed + length_index
            for case in generate_cases(generator, seed, args.count, length):
                total_cases_by_length[length] += 1
                case_id = case["case_id"]
                rows = profile_case(
                    tool, case["seq"], case["G"], args.widths, args.param_file
                )
                for row in rows:
                    if raw_handle is not None and raw_writer is None:
                        raw_writer = csv.DictWriter(
                            raw_handle,
                            fieldnames=["case_id", *row.keys()],
                            delimiter="\t",
                            lineterminator="\n",
                        )
                        raw_writer.writeheader()
                    if raw_writer is not None:
                        raw_writer.writerow({"case_id": case_id, **row})

                    key = (int(row["n"]), row["family"], int(row["width"]))
                    bucket = stats[key]
                    bucket["bucket_count"] += 1
                    bucket["candidate_count"] += int(row["candidate_count"])
                    bucket["proposed_unique_count"] += int(
                        row["proposed_unique_count"]
                    )
                    bucket["proposed_unique_max"] = max(
                        bucket["proposed_unique_max"],
                        int(row["proposed_unique_count"]),
                    )
                    bucket["max_sources_per_proposed_identity"] = max(
                        bucket["max_sources_per_proposed_identity"],
                        int(row["max_sources_per_proposed_identity"]),
                    )
                    bucket["candidate_max"] = max(
                        bucket["candidate_max"], int(row["candidate_count"])
                    )
                    bucket["kept_count"] += int(row["kept_count"])
                    bucket["selected_source_count"] += int(
                        row["selected_source_count"]
                    )
                    bucket["selected_in_top"] += int(row["selected_in_top"])
                    bucket["kept_unique_count"] += int(row["kept_unique_count"])
                    bucket["selected_in_top_unique"] += int(
                        row["selected_in_top_unique"]
                    )
                    bucket["source_sum"] += float(row["source_sum"])
                    bucket["top_source_sum"] += float(row["top_source_sum"])
                    bucket["top_unique_source_sum"] += float(
                        row["top_unique_source_sum"]
                    )
                    cases_with_family[key].add(case_id)
                    priority_kinds[key].add(row["priority_kind"])
    finally:
        if raw_handle is not None:
            raw_handle.close()

    print(
        "n\tfamily\twidth\tpriority_kind\ttotal_cases\tcases_with_family\tbucket_count"
        "\tcandidate_count\tcandidate_mean\tcandidate_max"
        "\tproposed_unique_count\tproposed_unique_max\tunique_fraction"
        "\tmax_sources_per_proposed_identity\tkept_count"
        "\tselected_source_count\tselected_in_top\tselected_coverage"
        "\tkept_unique_count\tselected_in_top_unique"
        "\tselected_unique_coverage\tsource_sum\ttop_source_sum"
        "\tsource_coverage\ttop_unique_source_sum\tunique_source_coverage"
    )
    for key in sorted(stats):
        length, family, width = key
        value = stats[key]
        buckets = int(value["bucket_count"])
        selected = int(value["selected_source_count"])
        source_sum = value["source_sum"]
        print(
            f"{length}\t{family}\t{width}\t{','.join(sorted(priority_kinds[key]))}\t"
            f"{total_cases_by_length[length]}\t"
            f"{len(cases_with_family[key])}\t{buckets}\t"
            f"{int(value['candidate_count'])}\t"
            f"{value['candidate_count'] / buckets:.17g}\t"
            f"{int(value['candidate_max'])}\t"
            f"{int(value['proposed_unique_count'])}\t"
            f"{int(value['proposed_unique_max'])}\t"
            f"{value['proposed_unique_count'] / value['candidate_count']:.17g}\t"
            f"{int(value['max_sources_per_proposed_identity'])}\t"
            f"{int(value['kept_count'])}\t"
            f"{selected}\t{int(value['selected_in_top'])}\t"
            f"{(value['selected_in_top'] / selected if selected else 1.0):.17g}\t"
            f"{int(value['kept_unique_count'])}\t"
            f"{int(value['selected_in_top_unique'])}\t"
            f"{(value['selected_in_top_unique'] / selected if selected else 1.0):.17g}\t"
            f"{source_sum:.17g}\t{value['top_source_sum']:.17g}\t"
            f"{(value['top_source_sum'] / source_sum if source_sum else 0.0):.17g}\t"
            f"{value['top_unique_source_sum']:.17g}\t"
            f"{(value['top_unique_source_sum'] / source_sum if source_sum else 0.0):.17g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
