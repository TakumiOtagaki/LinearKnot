#!/usr/bin/env python3
"""Tests for the full-source-row PseudoBase++ audit manifest."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_pseudobase_benchmark import (  # noqa: E402
    audit_dataset,
    atomic_write_ids,
    decomposable_into_two_noncrossing_layers,
    strict_density2_status,
)
from scfg2_pair_accuracy import parse_dot_bracket  # noqa: E402


def main() -> int:
    crossing = parse_dot_bracket("(.[..).]")
    two_layers, edges = decomposable_into_two_noncrossing_layers(crossing)
    assert two_layers
    assert edges == 1

    odd_crossing_cycle = parse_dot_bracket("([{)]}")
    two_layers, edges = decomposable_into_two_noncrossing_layers(
        odd_crossing_cycle)
    assert not two_layers
    assert edges == 3

    nested_bands = {(1, 3), (2, 7), (4, 6), (5, 8)}
    two_layers, _ = decomposable_into_two_noncrossing_layers(nested_bands)
    assert two_layers
    strict, reason = strict_density2_status(nested_bands)
    assert not strict
    assert reason == "band_component_coverage_exceeds_two"

    with tempfile.TemporaryDirectory(
            prefix="pseudobase-audit-test-") as tmpdir:
        dataset = Path(tmpdir) / "dataset.csv"
        with dataset.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "pkb_short",
                    "sequence",
                    "length",
                    "dotbracket",
                    "pipeline_ready",
                    "rna_type",
                    "density2_status",
                    "density2_reason",
                ),
            )
            writer.writeheader()
            writer.writerows((
                {
                    "pkb_short": "clean",
                    "sequence": "GCAUCGCA",
                    "length": "8",
                    "dotbracket": "(.[..).]",
                    "pipeline_ready": "true",
                    "rna_type": "test",
                    "density2_status": "True",
                    "density2_reason": "bipartite",
                },
                {
                    "pkb_short": "three-layer",
                    "sequence": "GCAUGC",
                    "length": "6",
                    "dotbracket": "([{)]}",
                    "pipeline_ready": "true",
                    "rna_type": "test",
                    "density2_status": "False",
                    "density2_reason": "not_bipartite",
                },
                {
                    "pkb_short": "pk-free",
                    "sequence": "GCAUGC",
                    "length": "6",
                    "dotbracket": "(....)",
                    "pipeline_ready": "true",
                    "rna_type": "test",
                    "density2_status": "True",
                    "density2_reason": "bipartite",
                },
                {
                    "pkb_short": "ambiguous",
                    "sequence": "GCAY",
                    "length": "4",
                    "dotbracket": "....",
                    "pipeline_ready": "false",
                    "rna_type": "test",
                    "density2_status": "",
                    "density2_reason": "",
                },
                {
                    "pkb_short": "invalid-reference",
                    "sequence": "GCAU",
                    "length": "4",
                    "dotbracket": "((..",
                    "pipeline_ready": "true",
                    "rna_type": "test",
                    "density2_status": "",
                    "density2_reason": "",
                },
            ))

        manifest, summary = audit_dataset(
            dataset, minimum_loop_unpaired=3)
        by_id = {row["pkb_short"]: row for row in manifest}
        assert len(manifest) == 5
        assert by_id["clean"]["clean_pseudoknotted_eligible"] == "true"
        assert by_id["clean"]["strict_structural_eligible"] == "true"
        assert by_id["clean"]["density2_status_source"] == "True"
        assert (
            by_id["three-layer"]["two_noncrossing_layers"] == "false")
        assert by_id["three-layer"]["strict_density2"] == "false"
        assert "not_strict_density2" in (
            by_id["three-layer"]["exclusion_reasons"])
        assert by_id["pk-free"]["clean_pseudoknotted_eligible"] == "false"
        assert "no_crossing_pair" in by_id["pk-free"]["exclusion_reasons"]
        assert by_id["ambiguous"]["unsupported_residues"] == "Y"
        assert "not_pipeline_ready" in by_id["ambiguous"]["exclusion_reasons"]
        assert (
            by_id["invalid-reference"]["reference_parse_ok"] == "false")
        assert "invalid_reference" in (
            by_id["invalid-reference"]["exclusion_reasons"])
        assert summary["source_rows"] == 5
        assert summary["clean_pseudoknotted_rows"] == 2
        assert summary["strict_structural_rows"] == 1
        assert len(str(summary["source_sha256"])) == 64
        ids_path = Path(tmpdir) / "clean.ids"
        assert atomic_write_ids(
            ids_path, manifest, "clean_pseudoknotted_eligible") == 2
        assert ids_path.read_text() == "clean\nthree-layer\n"

    print("scfg2_pseudobase_audit=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
