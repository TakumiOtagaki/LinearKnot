#!/usr/bin/env python3
"""Unit tests for method-independent PseudoBase++ shard aggregation."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from summarize_pseudobase_benchmark import (  # noqa: E402
    summarize,
    validate_and_order,
)


def scored(case_id: str, status: str = "ok") -> dict[str, str]:
    row = {
        "dataset_sha256": "dataset",
        "run_provenance_sha256": "run",
        "pkb_short": case_id,
        "status": status,
        "execution_seconds": "2",
        "peak_rss_kb": "100",
        "exact": "1",
    }
    for family in ("all", "crossing", "noncrossing"):
        row[f"{family}_tp"] = "2"
        row[f"{family}_fp"] = "0"
        row[f"{family}_fn"] = "0"
        row[f"{family}_ppv"] = "1"
        row[f"{family}_sensitivity"] = "1"
        row[f"{family}_f1"] = "1"
    return row


def main() -> int:
    rows = [
        scored("PKB2"),
        {
            "dataset_sha256": "dataset",
            "run_provenance_sha256": "run",
            "pkb_short": "PKB1",
            "status": "timeout",
            "execution_seconds": "5",
            "peak_rss_kb": "200",
        },
    ]
    ordered = validate_and_order(rows, ["PKB1", "PKB2"])
    assert [row["pkb_short"] for row in ordered] == ["PKB1", "PKB2"]
    summary = summarize(ordered)
    assert summary["requested_cases"] == 2
    assert summary["completed_cases"] == 1
    assert summary["failed_cases"] == 1
    assert summary["maximum_peak_rss_kb"] == 200.0
    assert summary["macro_all_f1"] == 1.0
    assert summary["micro_all_f1"] == 1.0

    try:
        validate_and_order([scored("PKB1"), scored("PKB1")], None)
    except ValueError as error:
        assert "duplicate" in str(error)
    else:
        raise AssertionError("duplicate result ID was accepted")

    mixed = [scored("PKB1"), scored("PKB2")]
    mixed[1]["run_provenance_sha256"] = "different"
    try:
        validate_and_order(mixed, None)
    except ValueError as error:
        assert "provenance" in str(error)
    else:
        raise AssertionError("mixed run provenance was accepted")

    print("pseudobase_benchmark_summary=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
