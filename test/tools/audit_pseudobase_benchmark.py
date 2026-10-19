#!/usr/bin/env python3
"""Audit every source row in a PseudoBase++ benchmark snapshot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from collections import Counter, deque
from pathlib import Path
from typing import Iterable, Sequence


from scfg2_pair_accuracy import Pair, crossing_pairs, parse_dot_bracket


REQUIRED_COLUMNS = {
    "pkb_short",
    "sequence",
    "length",
    "dotbracket",
    "pipeline_ready",
}

MANIFEST_FIELDS = (
    "source_row",
    "pkb_short",
    "duplicate_id",
    "rna_type",
    "reference_column",
    "source_total_pairs",
    "canonical_pairs_source",
    "noncanonical_pairs_source",
    "removed_noncanonical_pairs_source",
    "declared_length",
    "sequence_length",
    "reference_length",
    "pipeline_ready",
    "unsupported_residues",
    "length_consistent",
    "reference_parse_ok",
    "reference_error",
    "reference_pair_count",
    "reference_crossing_pair_count",
    "crossing_edge_count",
    "has_crossing_pair",
    "two_noncrossing_layers",
    "strict_density2",
    "strict_density2_reason",
    "minimum_loop_unpaired_observed",
    "minimum_loop_compatible",
    "density2_status_source",
    "density2_reason_source",
    "clean_pseudoknotted_eligible",
    "strict_structural_eligible",
    "exclusion_reasons",
    "source_sha256",
)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def parse_source_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def pair_crosses(lhs: Pair, rhs: Pair) -> bool:
    i, j = lhs
    k, l = rhs
    return (i < k < j < l) or (k < i < l < j)


def pair_nested(lhs: Pair, rhs: Pair) -> bool:
    i, j = lhs
    k, l = rhs
    return (i < k < l < j) or (k < i < j < l)


def crossing_graph(pairs: Iterable[Pair]) -> tuple[list[list[int]], int]:
    ordered = sorted(pairs)
    graph = [[] for _ in ordered]
    edges = 0
    for left in range(len(ordered)):
        for right in range(left + 1, len(ordered)):
            if not pair_crosses(ordered[left], ordered[right]):
                continue
            graph[left].append(right)
            graph[right].append(left)
            edges += 1
    return graph, edges


def decomposable_into_two_noncrossing_layers(
    pairs: Iterable[Pair],
) -> tuple[bool, int]:
    graph, edges = crossing_graph(pairs)
    colors: list[int | None] = [None] * len(graph)
    for start in range(len(graph)):
        if colors[start] is not None:
            continue
        colors[start] = 0
        pending = deque([start])
        while pending:
            current = pending.popleft()
            assert colors[current] is not None
            for neighbor in graph[current]:
                if colors[neighbor] is None:
                    colors[neighbor] = 1 - colors[current]
                    pending.append(neighbor)
                elif colors[neighbor] == colors[current]:
                    return False, edges
    return True, edges


def strict_density2_status(pairs: Iterable[Pair]) -> tuple[bool, str]:
    """Apply the benchmark's band-aware strict density-2 definition."""
    ordered = sorted(pairs)
    graph, _ = crossing_graph(ordered)
    two_layers, _ = decomposable_into_two_noncrossing_layers(ordered)
    if not two_layers:
        return False, "pair_crossing_graph_not_bipartite"

    signature_groups: dict[frozenset[int], list[int]] = {}
    for index, neighbors in enumerate(graph):
        signature = frozenset(neighbors)
        if signature:
            signature_groups.setdefault(signature, []).append(index)

    bands: list[list[int]] = []
    for indices in signature_groups.values():
        nested_graph: dict[int, list[int]] = {
            index: [] for index in indices}
        for position, left in enumerate(indices):
            for right in indices[position + 1:]:
                if pair_nested(ordered[left], ordered[right]):
                    nested_graph[left].append(right)
                    nested_graph[right].append(left)
        visited: set[int] = set()
        for start in indices:
            if start in visited:
                continue
            visited.add(start)
            pending = [start]
            band: list[int] = []
            while pending:
                current = pending.pop()
                band.append(current)
                for neighbor in nested_graph[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        pending.append(neighbor)
            bands.append(band)

    band_graph: list[set[int]] = [set() for _ in bands]
    for left in range(len(bands)):
        for right in range(left + 1, len(bands)):
            if any(
                pair_crosses(ordered[left_pair], ordered[right_pair])
                for left_pair in bands[left]
                for right_pair in bands[right]
            ):
                band_graph[left].add(right)
                band_graph[right].add(left)

    visited_bands: set[int] = set()
    for start in range(len(bands)):
        if start in visited_bands:
            continue
        visited_bands.add(start)
        pending = deque([start])
        component: list[int] = []
        while pending:
            current = pending.popleft()
            component.append(current)
            for neighbor in band_graph[current]:
                if neighbor not in visited_bands:
                    visited_bands.add(neighbor)
                    pending.append(neighbor)
        for position in range(1, max((right for _, right in ordered), default=0) + 1):
            covering = sum(
                any(
                    ordered[pair_index][0] <= position
                    <= ordered[pair_index][1]
                    for pair_index in bands[band_index]
                )
                for band_index in component
            )
            if covering > 2:
                return False, "band_component_coverage_exceeds_two"
    return True, "strict_density2"


def safe_error(error: Exception) -> str:
    return " ".join(str(error).split())


def audit_row(
    row: dict[str, str],
    source_row: int,
    source_sha256: str,
    minimum_loop_unpaired: int,
) -> dict[str, object]:
    case_id = row.get("pkb_short", "").strip()
    sequence = row.get("sequence", "").strip().upper()
    reference = row.get("dotbracket", "").strip()
    pipeline_ready = parse_source_bool(row.get("pipeline_ready", ""))
    unsupported = "".join(sorted(set(sequence) - set("ACGU")))

    reasons: list[str] = []
    try:
        declared_length = int(row.get("length", ""))
        declared_length_valid = declared_length >= 0
    except ValueError:
        declared_length = -1
        declared_length_valid = False
    if not declared_length_valid:
        reasons.append("invalid_declared_length")

    length_consistent = (
        declared_length_valid
        and len(sequence) == declared_length
        and len(reference) == declared_length
    )
    if not pipeline_ready:
        reasons.append("not_pipeline_ready")
    if unsupported:
        reasons.append(f"unsupported_residues:{unsupported}")
    if declared_length_valid and len(sequence) != declared_length:
        reasons.append("sequence_length_mismatch")
    if declared_length_valid and len(reference) != declared_length:
        reasons.append("reference_length_mismatch")

    pairs: frozenset[Pair] = frozenset()
    reference_error = ""
    try:
        pairs = parse_dot_bracket(reference)
        reference_parse_ok = True
    except ValueError as error:
        reference_parse_ok = False
        reference_error = safe_error(error)
        reasons.append("invalid_reference")

    crossing_count = 0
    crossing_edge_count = 0
    two_layers = False
    strict_density2 = False
    strict_density2_reason = "reference_not_parseable"
    minimum_observed: int | str = ""
    minimum_compatible = False
    if reference_parse_ok:
        crossing_count = len(crossing_pairs(pairs))
        two_layers, crossing_edge_count = (
            decomposable_into_two_noncrossing_layers(pairs)
        )
        strict_density2, strict_density2_reason = (
            strict_density2_status(pairs)
        )
        if pairs:
            minimum_observed = min(j - i - 1 for i, j in pairs)
            minimum_compatible = minimum_observed >= minimum_loop_unpaired
        else:
            minimum_compatible = True

    has_crossing = crossing_count > 0
    basic_valid = (
        pipeline_ready
        and not unsupported
        and length_consistent
        and reference_parse_ok
    )
    clean_eligible = basic_valid and has_crossing
    strict_eligible = (
        clean_eligible and strict_density2 and minimum_compatible
    )
    if basic_valid and not has_crossing:
        reasons.append("no_crossing_pair")
    if clean_eligible and not strict_density2:
        reasons.append(f"not_strict_density2:{strict_density2_reason}")
    if clean_eligible and not minimum_compatible:
        reasons.append("minimum_loop_incompatible")

    return {
        "source_row": source_row,
        "pkb_short": case_id,
        "duplicate_id": "false",
        "rna_type": row.get("rna_type", ""),
        "reference_column": "dotbracket",
        "source_total_pairs": row.get("total_pairs", ""),
        "canonical_pairs_source": row.get("canonical_pairs", ""),
        "noncanonical_pairs_source": row.get("noncanonical_pairs", ""),
        "removed_noncanonical_pairs_source":
            row.get("removed_noncanonical_pairs", ""),
        "declared_length": row.get("length", ""),
        "sequence_length": len(sequence),
        "reference_length": len(reference),
        "pipeline_ready": bool_text(pipeline_ready),
        "unsupported_residues": unsupported,
        "length_consistent": bool_text(length_consistent),
        "reference_parse_ok": bool_text(reference_parse_ok),
        "reference_error": reference_error,
        "reference_pair_count": len(pairs),
        "reference_crossing_pair_count": crossing_count,
        "crossing_edge_count": crossing_edge_count,
        "has_crossing_pair": bool_text(has_crossing),
        "two_noncrossing_layers": bool_text(two_layers),
        "strict_density2": bool_text(strict_density2),
        "strict_density2_reason": strict_density2_reason,
        "minimum_loop_unpaired_observed": minimum_observed,
        "minimum_loop_compatible": bool_text(minimum_compatible),
        "density2_status_source": row.get("density2_status", ""),
        "density2_reason_source": row.get("density2_reason", ""),
        "clean_pseudoknotted_eligible": bool_text(clean_eligible),
        "strict_structural_eligible": bool_text(strict_eligible),
        "exclusion_reasons": ";".join(reasons),
        "source_sha256": source_sha256,
    }


def atomic_write_rows(
    path: Path,
    rows: Sequence[dict[str, object]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def atomic_write_ids(
    path: Path,
    manifest: Sequence[dict[str, object]],
    eligibility_field: str,
) -> int:
    selected = [
        str(row["pkb_short"])
        for row in manifest
        if row[eligibility_field] == "true"
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        for case_id in selected:
            handle.write(f"{case_id}\n")
    os.replace(temporary, path)
    return len(selected)


def audit_dataset(
    path: Path,
    minimum_loop_unpaired: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    payload = path.read_bytes()
    source_sha256 = hashlib.sha256(payload).hexdigest()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_COLUMNS - columns)
        if missing:
            raise ValueError(
                f"dataset is missing required columns: {','.join(missing)}")
        source_rows = list(reader)

    manifest = [
        audit_row(
            row,
            source_row=index + 2,
            source_sha256=source_sha256,
            minimum_loop_unpaired=minimum_loop_unpaired,
        )
        for index, row in enumerate(source_rows)
    ]
    ids = [str(row["pkb_short"]) for row in manifest]
    id_counts = Counter(ids)
    duplicate_ids = sorted(
        case_id for case_id, count in id_counts.items()
        if case_id and count > 1
    )
    for row in manifest:
        case_id = str(row["pkb_short"])
        added_reasons: list[str] = []
        if not case_id:
            added_reasons.append("missing_id")
        if case_id and id_counts[case_id] > 1:
            row["duplicate_id"] = "true"
            added_reasons.append("duplicate_id")
        if added_reasons:
            existing = str(row["exclusion_reasons"])
            row["exclusion_reasons"] = ";".join(
                [reason for reason in (existing, *added_reasons) if reason])
            row["clean_pseudoknotted_eligible"] = "false"
            row["strict_structural_eligible"] = "false"
    exclusion_counts = Counter(
        reason
        for row in manifest
        for reason in str(row["exclusion_reasons"]).split(";")
        if reason
    )
    summary: dict[str, object] = {
        "dataset": str(path),
        "source_bytes": len(payload),
        "source_sha256": source_sha256,
        "source_columns": ",".join(reader.fieldnames or []),
        "source_rows": len(manifest),
        "unique_nonempty_ids": len(set(ids) - {""}),
        "duplicate_id_count": len(duplicate_ids),
        "duplicate_id_rows": sum(
            row["duplicate_id"] == "true" for row in manifest),
        "duplicate_ids": ",".join(duplicate_ids),
        "minimum_loop_unpaired": minimum_loop_unpaired,
        "pipeline_ready_rows": sum(
            row["pipeline_ready"] == "true" for row in manifest),
        "reference_parseable_rows": sum(
            row["reference_parse_ok"] == "true" for row in manifest),
        "crossing_reference_rows": sum(
            row["has_crossing_pair"] == "true" for row in manifest),
        "two_layer_rows": sum(
            row["two_noncrossing_layers"] == "true" for row in manifest),
        "strict_density2_rows": sum(
            row["strict_density2"] == "true" for row in manifest),
        "clean_pseudoknotted_rows": sum(
            row["clean_pseudoknotted_eligible"] == "true"
            for row in manifest),
        "strict_structural_rows": sum(
            row["strict_structural_eligible"] == "true"
            for row in manifest),
        "exclusion_reason_counts": ",".join(
            f"{reason}={count}"
            for reason, count in sorted(exclusion_counts.items())
        ),
    }
    return manifest, summary


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--clean-ids-output", type=Path)
    parser.add_argument("--strict-ids-output", type=Path)
    parser.add_argument("--minimum-loop-unpaired", type=int, default=3)
    args = parser.parse_args(argv)
    if args.minimum_loop_unpaired < 0:
        parser.error("--minimum-loop-unpaired must be nonnegative")

    try:
        manifest, summary = audit_dataset(
            args.dataset, args.minimum_loop_unpaired)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    atomic_write_rows(args.output, manifest, MANIFEST_FIELDS)
    summary_path = (
        args.summary_output or args.output.with_suffix(".summary.tsv"))
    atomic_write_rows(summary_path, [summary], tuple(summary))
    if args.clean_ids_output is not None:
        atomic_write_ids(
            args.clean_ids_output,
            manifest,
            "clean_pseudoknotted_eligible",
        )
    if args.strict_ids_output is not None:
        atomic_write_ids(
            args.strict_ids_output,
            manifest,
            "strict_structural_eligible",
        )
    print(
        f"[pseudobase-audit] rows={summary['source_rows']}"
        f" clean={summary['clean_pseudoknotted_rows']}"
        f" strict={summary['strict_structural_rows']}"
        f" output={args.output}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
