#!/usr/bin/env python3
"""Evaluate IPknot predictions with the SCFG2 pair-set metrics."""

from __future__ import annotations

import argparse
import concurrent.futures
import shlex
import socket
import sys
import tempfile
import time
from pathlib import Path
from typing import Sequence


sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_scfg2_pseudobase_viterbi import (  # noqa: E402
    load_dataset,
    load_resumable_rows,
    provenance_sha256,
    run_captured_command,
    safe_case_name,
    select_case_ids,
    sha256_file,
    status_summary,
    validate_case,
    write_rows,
)
from scfg2_pair_accuracy import (  # noqa: E402
    PairScore,
    accuracy_fields,
    evaluate_pair_sets,
    parse_dot_bracket,
)


def parse_ipknot_output(text: str) -> tuple[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 3 or not lines[0].startswith(">"):
        raise ValueError("IPknot output is not one three-line FASTA record")
    return lines[1], lines[2]


def aggregate(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    if not rows:
        raise ValueError("cannot aggregate an empty result")
    count = len(rows)
    metric_names = (
        "all_ppv", "all_sensitivity", "all_f1",
        "crossing_ppv", "crossing_sensitivity", "crossing_f1",
        "noncrossing_ppv", "noncrossing_sensitivity", "noncrossing_f1",
    )
    out: dict[str, object] = {"cases": count}
    for name in metric_names:
        out[f"macro_{name}"] = sum(float(row[name]) for row in rows) / count
    for family in ("all", "crossing", "noncrossing"):
        score = PairScore(
            sum(int(row[f"{family}_tp"]) for row in rows),
            sum(int(row[f"{family}_fp"]) for row in rows),
            sum(int(row[f"{family}_fn"]) for row in rows),
        )
        out[f"micro_{family}_tp"] = score.true_positive
        out[f"micro_{family}_fp"] = score.false_positive
        out[f"micro_{family}_fn"] = score.false_negative
        out[f"micro_{family}_ppv"] = score.ppv
        out[f"micro_{family}_sensitivity"] = score.sensitivity
        out[f"micro_{family}_f1"] = score.f1
    out["exact_rate"] = sum(int(row["exact"]) for row in rows) / count
    out["mean_total_seconds"] = sum(
        float(row["total_seconds"]) for row in rows) / count
    return out


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--ids")
    selection.add_argument("--ids-file", type=Path)
    selection.add_argument("--all", action="store_true")
    execution = parser.add_mutually_exclusive_group(required=True)
    execution.add_argument("--binary", type=Path)
    execution.add_argument("--docker-image")
    parser.add_argument("--docker-platform", default="")
    parser.add_argument("--model", default="LinearPartition-C")
    parser.add_argument("--beam-width", type=int, default=100)
    parser.add_argument("--solver-threads", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--skip-invalid", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=3600.0)
    parser.add_argument("--artifacts-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args(argv)

    if (
        args.jobs <= 0
        or args.solver_threads <= 0
        or args.timeout_seconds < 0
    ):
        parser.error(
            "--jobs/--solver-threads must be positive and timeout nonnegative")
    artifacts_dir = (
        args.artifacts_dir
        if args.artifacts_dir is not None
        else args.output.parent / f"{args.output.stem}.artifacts"
    )
    dataset = load_dataset(args.dataset)
    try:
        ids = select_case_ids(
            dataset, args.ids, args.ids_file, args.all)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    try:
        execution_sha256 = (
            sha256_file(args.binary)
            if args.binary is not None else f"docker:{args.docker_image}")
        provenance: dict[str, object] = {
            "dataset_sha256": sha256_file(args.dataset),
            "execution_sha256": execution_sha256,
            "execution_mode": (
                "native" if args.binary is not None else "docker"),
            "docker_platform": args.docker_platform,
            "model": args.model,
            "beam_width": args.beam_width,
            "solver_threads": args.solver_threads,
            "jobs": args.jobs,
            "timeout_seconds": args.timeout_seconds,
        }
        provenance["run_provenance_sha256"] = provenance_sha256(provenance)
    except OSError as error:
        parser.error(str(error))

    valid_ids: list[str] = []
    invalid_cases: list[tuple[str, str]] = []
    for case_id in ids:
        try:
            validate_case(dataset[case_id])
        except ValueError as error:
            invalid_cases.append((case_id, str(error)))
        else:
            valid_ids.append(case_id)
    if invalid_cases and not args.skip_invalid:
        details = "; ".join(
            f"{case_id}: {reason}" for case_id, reason in invalid_cases)
        parser.error(f"invalid cases found before execution: {details}")
    for case_id, reason in invalid_cases:
        print(f"[ipknot] skip {case_id}: {reason}", file=sys.stderr)
    if not valid_ids and not args.skip_invalid:
        parser.error("no valid cases remain after input validation")

    try:
        results_by_id = (
            load_resumable_rows(args.output, provenance, ids)
            if args.resume else {}
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    reused_ok = sum(
        row.get("status") == "ok" for row in results_by_id.values())
    for case_id, reason in invalid_cases:
        row = dataset[case_id]
        results_by_id[case_id] = {
            **provenance,
            "pkb_short": case_id,
            "length": row.get("length", ""),
            "rna_type": row.get("rna_type", ""),
            "reference_structure": row.get("dotbracket", ""),
            "host": socket.gethostname(),
            "stdout_path": "",
            "stderr_path": "",
            "status": "invalid_input",
            "exit_status": "",
            "execution_seconds": 0.0,
            "peak_rss_kb": "",
            "error": reason,
        }
    pending_ids = [
        case_id for case_id in valid_ids
        if case_id not in results_by_id
        or results_by_id[case_id].get("status") != "ok"
    ]
    if reused_ok:
        print(
            f"[ipknot] resume {reused_ok} completed cases;"
            f" {len(pending_ids)} pending",
            file=sys.stderr,
            flush=True,
        )

    def run_one(case_id: str) -> dict[str, object]:
        row = dataset[case_id]
        sequence = row["sequence"].upper()
        artifact_root = artifacts_dir / safe_case_name(case_id)
        stdout_path = artifact_root / "stdout.txt"
        stderr_path = artifact_root / "stderr.txt"
        metrics_path = artifact_root / "peak_rss_kb.txt"
        common: dict[str, object] = {
            **provenance,
            "pkb_short": case_id,
            "length": row["length"],
            "rna_type": row["rna_type"],
            "reference_structure": row["dotbracket"],
            "host": socket.gethostname(),
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
        }
        with tempfile.TemporaryDirectory(prefix="ipknot-") as tmpdir:
            fasta = Path(tmpdir) / "input.fa"
            fasta.write_text(
                f">{case_id}\n{sequence}\n", encoding="utf-8")
            options = [
                "-e", args.model,
                "--beam-size", str(args.beam_width),
                "-n", str(args.solver_threads),
            ]
            if args.binary is not None:
                command = [str(args.binary), *options, str(fasta)]
            else:
                command = ["docker", "run", "--rm"]
                if args.docker_platform:
                    command.extend(("--platform", args.docker_platform))
                command.extend((
                    "-v", f"{tmpdir}:/data:ro",
                    str(args.docker_image), "ipknot", *options,
                    "/data/input.fa",
                ))
            common["command"] = shlex.join(command)
            completed = run_captured_command(
                command,
                args.timeout_seconds,
                stdout_path,
                stderr_path,
                metrics_path,
            )
        execution = {
            "exit_status": completed.exit_status,
            "execution_seconds": completed.elapsed_seconds,
            "peak_rss_kb": (
                completed.peak_rss_kb
                if completed.peak_rss_kb is not None else ""),
        }
        if completed.timed_out:
            return {
                **common,
                "status": "timeout",
                **execution,
                "error": f"exceeded {args.timeout_seconds:g} seconds",
            }
        if completed.exit_status != 0:
            return {
                **common,
                "status": "tool_failure",
                **execution,
                "error": f"predictor exited with status {completed.exit_status}",
            }
        try:
            predicted_sequence, structure = parse_ipknot_output(
                completed.stdout)
            if predicted_sequence.upper() != sequence:
                raise ValueError(
                    f"{case_id}: IPknot sequence does not match input")
            if len(structure) != len(sequence):
                raise ValueError(f"{case_id}: IPknot structure length mismatch")
            accuracy = evaluate_pair_sets(
                parse_dot_bracket(structure),
                parse_dot_bracket(row["dotbracket"]),
            )
        except (KeyError, TypeError, ValueError) as error:
            return {
                **common,
                "status": "parse_failure",
                **execution,
                "error": " ".join(str(error).split()),
            }
        print(
            f"[ipknot] done {case_id} F1={accuracy.all_pairs.f1:.6g}",
            file=sys.stderr,
            flush=True,
        )
        return {
            **common,
            "status": "ok",
            **execution,
            "error": "",
            "predicted_structure": structure,
            "total_seconds": completed.elapsed_seconds,
            **accuracy_fields(accuracy),
        }

    batch_begin = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_one, case_id): case_id
            for case_id in pending_ids
        }
        for future in concurrent.futures.as_completed(futures):
            case_id = futures[future]
            results_by_id[case_id] = future.result()
            write_rows(
                args.output,
                [results_by_id[key] for key in ids
                 if key in results_by_id],
            )
    batch_wall_seconds = time.perf_counter() - batch_begin
    results = [results_by_id[case_id] for case_id in ids]
    write_rows(args.output, results)
    summary_path = (
        args.summary_output or args.output.with_suffix(".summary.tsv"))
    successful_results = [
        row for row in results if row.get("status") == "ok"]
    summary: dict[str, object] = {
        "model": args.model,
        "beam_width": args.beam_width,
        "solver_threads": args.solver_threads,
        "jobs": args.jobs,
        "requested_cases": len(ids),
        "invalid_input_ids": ",".join(
            case_id for case_id, _ in invalid_cases),
        "batch_wall_seconds": batch_wall_seconds,
        **status_summary(results),
    }
    if successful_results:
        summary.update(aggregate(successful_results))
    write_rows(summary_path, [summary])
    print(
        f"[ipknot] wrote {len(results)} cases to {args.output}"
        f" and summary to {summary_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
