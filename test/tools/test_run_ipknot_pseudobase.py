#!/usr/bin/env python3
"""Unit tests for the IPknot PseudoBase++ evaluator."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_ipknot_pseudobase import aggregate, parse_ipknot_output  # noqa: E402
from run_scfg2_pseudobase_viterbi import (  # noqa: E402
    load_resumable_rows,
    write_rows,
)
from scfg2_pair_accuracy import accuracy_fields, evaluate_pair_sets  # noqa: E402


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "ipknot.tsv"
        provenance = {
            "dataset_sha256": "dataset",
            "execution_sha256": "ipknot",
            "model": "LinearPartition-C",
            "beam_width": 100,
            "solver_threads": 1,
        }
        write_rows(output, [{
            **provenance,
            "pkb_short": "PKB1",
            "status": "ok",
        }])
        assert list(load_resumable_rows(
            output, provenance, ["PKB1"])) == ["PKB1"]

    sequence, structure = parse_ipknot_output(
        ">case\nACGUACGU\n((....))\n")
    assert sequence == "ACGUACGU"
    assert structure == "((....))"
    accuracy = evaluate_pair_sets({(1, 8), (2, 7)}, {(1, 8), (2, 7)})
    row = {"total_seconds": 0.25, **accuracy_fields(accuracy)}
    summary = aggregate([row])
    assert summary["macro_all_f1"] == 1.0
    assert summary["micro_all_f1"] == 1.0
    assert summary["exact_rate"] == 1.0
    print("ipknot_pseudobase_runner=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
