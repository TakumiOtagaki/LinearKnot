#!/usr/bin/env python3
"""Merge predictor shards and summarize the common PseudoBase++ contract."""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import tempfile
from collections import Counter
from pathlib import Path
from typing import Sequence


METRICS = (
    "all_ppv", "all_sensitivity", "all_f1",
    "crossing_ppv", "crossing_sensitivity", "crossing_f1",
    "noncrossing_ppv", "noncrossing_sensitivity", "noncrossing_f1",
)


def read_rows(paths: Sequence[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("input shards contain no rows")
    return rows


def read_ids(path: Path) -> list[str]:
    ids = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not ids:
        raise ValueError("expected IDs file is empty")
    if len(ids) != len(set(ids)):
        raise ValueError("expected IDs file contains duplicates")
    return ids


def validate_and_order(
    rows: Sequence[dict[str, str]],
    expected_ids: Sequence[str] | None,
) -> list[dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        case_id = row.get("pkb_short", "")
        if not case_id:
            raise ValueError("result row has no pkb_short")
        if case_id in by_id:
            raise ValueError(f"duplicate result ID: {case_id}")
        by_id[case_id] = row
    for field in ("dataset_sha256", "run_provenance_sha256"):
        values = {row.get(field, "") for row in rows}
        if "" in values or len(values) != 1:
            raise ValueError(f"mixed or missing provenance field: {field}")
    if expected_ids is None:
        return list(rows)
    missing = [case_id for case_id in expected_ids if case_id not in by_id]
    extra = [case_id for case_id in by_id if case_id not in set(expected_ids)]
    if missing or extra:
        raise ValueError(
            f"ID coverage mismatch: missing={','.join(missing)}"
            f" extra={','.join(extra)}")
    return [by_id[case_id] for case_id in expected_ids]


def numeric(rows: Sequence[dict[str, str]], field: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = row.get(field, "")
        if value != "":
            values.append(float(value))
    return values


def ratio(numerator: int, denominator: int, empty: float) -> float:
    return numerator / denominator if denominator else empty


def summarize(rows: Sequence[dict[str, str]]) -> dict[str, object]:
    statuses = Counter(row.get("status", "unknown") for row in rows)
    successful = [row for row in rows if row.get("status") == "ok"]
    summary: dict[str, object] = {
        "dataset_sha256": rows[0]["dataset_sha256"],
        "run_provenance_sha256": rows[0]["run_provenance_sha256"],
        "requested_cases": len(rows),
        "completed_cases": len(successful),
        "failed_cases": len(rows) - len(successful),
        "status_counts": ",".join(
            f"{status}={count}"
            for status, count in sorted(statuses.items())),
    }
    executed = [
        row for row in rows if row.get("status") != "invalid_input"]
    elapsed = numeric(executed, "execution_seconds")
    rss = numeric(executed, "peak_rss_kb")
    if elapsed:
        summary.update({
            "mean_execution_seconds": statistics.fmean(elapsed),
            "median_execution_seconds": statistics.median(elapsed),
            "maximum_execution_seconds": max(elapsed),
        })
    if rss:
        summary.update({
            "mean_peak_rss_kb": statistics.fmean(rss),
            "median_peak_rss_kb": statistics.median(rss),
            "maximum_peak_rss_kb": max(rss),
        })
    if successful and all(
        metric in row and row[metric] != ""
        for row in successful for metric in METRICS
    ):
        for metric in METRICS:
            summary[f"macro_{metric}"] = statistics.fmean(
                float(row[metric]) for row in successful)
        for family in ("all", "crossing", "noncrossing"):
            tp = sum(int(row[f"{family}_tp"]) for row in successful)
            fp = sum(int(row[f"{family}_fp"]) for row in successful)
            fn = sum(int(row[f"{family}_fn"]) for row in successful)
            summary[f"micro_{family}_tp"] = tp
            summary[f"micro_{family}_fp"] = fp
            summary[f"micro_{family}_fn"] = fn
            summary[f"micro_{family}_ppv"] = ratio(tp, tp + fp, 1.0)
            summary[f"micro_{family}_sensitivity"] = ratio(
                tp, tp + fn, 1.0)
            precision = float(summary[f"micro_{family}_ppv"])
            recall = float(summary[f"micro_{family}_sensitivity"])
            summary[f"micro_{family}_f1"] = ratio(
                2.0 * precision * recall, precision + recall, 0.0)
        summary["exact_rate"] = statistics.fmean(
            int(row["exact"]) for row in successful)
    return summary


def atomic_write_rows(
    path: Path,
    rows: Sequence[dict[str, object] | dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--expected-ids-file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args(argv)
    try:
        rows = read_rows(args.inputs)
        expected_ids = (
            read_ids(args.expected_ids_file)
            if args.expected_ids_file is not None else None)
        ordered = validate_and_order(rows, expected_ids)
        summary = summarize(ordered)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    atomic_write_rows(args.output, ordered)
    summary_path = (
        args.summary_output or args.output.with_suffix(".summary.tsv"))
    atomic_write_rows(summary_path, [summary])
    print(
        f"[pseudobase-summary] rows={len(ordered)}"
        f" completed={summary['completed_cases']}"
        f" output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
