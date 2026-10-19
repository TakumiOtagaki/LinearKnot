#!/usr/bin/env python3
"""Evaluate sequence-only hotspot/LinearKnot Viterbi on PseudoBase++."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import shlex
import signal
import socket
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


sys.path.insert(0, str(Path(__file__).resolve().parent))

from scfg2_pair_accuracy import (  # noqa: E402
    PairScore,
    accuracy_fields,
    evaluate_pair_sets,
    parse_dot_bracket,
)


def comma_strings(text: str) -> list[str]:
    return [token.strip() for token in text.split(",") if token.strip()]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def provenance_sha256(values: dict[str, object]) -> str:
    payload = json.dumps(
        values, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def safe_case_name(case_id: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", case_id).strip("._")
    return safe or hashlib.sha256(case_id.encode()).hexdigest()[:16]


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


@dataclass(frozen=True)
class CommandResult:
    stdout: str
    stderr: str
    exit_status: int
    timed_out: bool
    elapsed_seconds: float
    peak_rss_kb: int | None


def run_captured_command(
    command: Sequence[str],
    timeout_seconds: float,
    stdout_path: Path,
    stderr_path: Path,
    metrics_path: Path,
) -> CommandResult:
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    timed_command = [
        "/usr/bin/time",
        "--format=%M",
        f"--output={metrics_path}",
        "--",
        *command,
    ]
    begin = time.perf_counter()
    process = subprocess.Popen(
        timed_command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    timed_out = False
    try:
        stdout, stderr = process.communicate(
            timeout=timeout_seconds if timeout_seconds > 0 else None)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
    elapsed = time.perf_counter() - begin
    atomic_write_text(stdout_path, stdout)
    atomic_write_text(stderr_path, stderr)
    peak_rss_kb: int | None = None
    try:
        metric = metrics_path.read_text(encoding="utf-8").strip()
        if metric:
            peak_rss_kb = int(metric.splitlines()[-1])
    except (FileNotFoundError, ValueError):
        pass
    return CommandResult(
        stdout=stdout,
        stderr=stderr,
        exit_status=process.returncode,
        timed_out=timed_out,
        elapsed_seconds=elapsed,
        peak_rss_kb=peak_rss_kb,
    )


def load_dataset(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {
            row["pkb_short"]: row
            for row in csv.DictReader(handle)
        }


def select_case_ids(
    dataset: dict[str, dict[str, str]],
    ids: str | None,
    ids_file: Path | None,
    select_all: bool,
) -> list[str]:
    if select_all:
        selected = list(dataset)
    elif ids_file is not None:
        selected = [
            line.strip()
            for line in ids_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    else:
        assert ids is not None
        selected = comma_strings(ids)
    if not selected:
        raise ValueError("case selection is empty")
    duplicates = sorted({
        case_id for case_id in selected if selected.count(case_id) > 1})
    if duplicates:
        raise ValueError(f"duplicate IDs selected: {','.join(duplicates)}")
    missing = [case_id for case_id in selected if case_id not in dataset]
    if missing:
        raise ValueError(f"dataset is missing IDs: {','.join(missing)}")
    return selected


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
    if len(row["dotbracket"]) != len(sequence):
        raise ValueError(
            f"{row['pkb_short']} has inconsistent reference length")
    parse_dot_bracket(row["dotbracket"])


def parse_key_values(fields: Sequence[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition("=")
        if not separator:
            raise ValueError(f"malformed key-value field: {field}")
        out[key] = value
    return out


def parse_prediction_output(text: str) -> dict[str, object]:
    summary: dict[str, str] | None = None
    record_header: list[str] | None = None
    records: dict[int, dict[str, str]] = {}
    group_header: list[str] | None = None
    groups: dict[int, dict[str, str]] = {}
    for raw in text.splitlines():
        fields = raw.split("\t")
        if fields[0] == "SUMMARY":
            summary = parse_key_values(fields[1:])
        elif fields[0] == "RECORD_HEADER":
            record_header = fields[1:]
        elif fields[0] == "RECORD":
            if record_header is None or len(fields) - 1 != len(record_header):
                raise ValueError("RECORD does not match RECORD_HEADER")
            row = dict(zip(record_header, fields[1:]))
            records[int(row["index"])] = row
        elif fields[0] == "GROUP_HEADER":
            group_header = fields[1:]
        elif fields[0] == "GROUP":
            if group_header is None or len(fields) - 1 != len(group_header):
                raise ValueError("GROUP does not match GROUP_HEADER")
            row = dict(zip(group_header, fields[1:]))
            groups[int(row["index"])] = row
    if summary is None:
        raise ValueError("prediction output has no SUMMARY")
    if summary.get("ok") != "1" or summary.get("reachable") != "1":
        raise ValueError(
            f"prediction failed: {summary.get('error', 'unknown error')}")
    selected_groups = [
        row for row in groups.values() if row["selected"] == "1"
    ]
    if len(selected_groups) != 1:
        raise ValueError("prediction has no unique selected GROUP")
    selected_group = selected_groups[0]
    selected_record_index = int(selected_group["selected_record"])
    if selected_record_index not in records:
        raise ValueError("selected GROUP refers to a missing RECORD")
    selected_record = records[selected_record_index]
    if selected_record["reachable"] != "1":
        raise ValueError("selected RECORD is unreachable")
    return {
        "summary": summary,
        "records": records,
        "groups": groups,
        "selected_group": selected_group,
        "selected_record": selected_record,
    }


def ratio(numerator: int, denominator: int, empty: float) -> float:
    return numerator / denominator if denominator else empty


def evaluate_case(
    row: dict[str, str],
    parsed: dict[str, object],
) -> dict[str, object]:
    reference = parse_dot_bracket(row["dotbracket"])
    summary = parsed["summary"]
    records = parsed["records"]
    groups = parsed["groups"]
    selected_group = parsed["selected_group"]
    selected_record = parsed["selected_record"]
    assert isinstance(summary, dict)
    assert isinstance(records, dict)
    assert isinstance(groups, dict)
    assert isinstance(selected_group, dict)
    assert isinstance(selected_record, dict)

    selected_pairs = parse_dot_bracket(selected_record["union_R"])
    selected_accuracy = evaluate_pair_sets(selected_pairs, reference)
    best_f1 = -1.0
    nonempty_scaffolds = 0
    any_nonempty_scaffold_subset = False
    best_scaffold_recall = 0.0
    best_compatible_scaffold_recall = 0.0
    reachable_records = 0
    for record in records.values():
        if record["reachable"] != "1":
            continue
        reachable_records += 1
        union_pairs = parse_dot_bracket(record["union_R"])
        best_f1 = max(
            best_f1, evaluate_pair_sets(union_pairs, reference).all_pairs.f1)
        scaffold_pairs = parse_dot_bracket(record["scaffold_G"])
        best_scaffold_recall = max(
            best_scaffold_recall,
            ratio(len(scaffold_pairs & reference), len(reference), 1.0),
        )
        if scaffold_pairs:
            nonempty_scaffolds += 1
            if scaffold_pairs.issubset(reference):
                any_nonempty_scaffold_subset = True
                best_compatible_scaffold_recall = max(
                    best_compatible_scaffold_recall,
                    ratio(len(scaffold_pairs), len(reference), 1.0),
                )

    spreads = [float(group["energy_spread"]) for group in groups.values()]
    return {
        "pkb_short": row["pkb_short"],
        "length": row["length"],
        "rna_type": row["rna_type"],
        "reference_pairs": len(reference),
        "reference_structure": row["dotbracket"],
        "selected_structure": selected_record["union_R"],
        "selected_scaffold": selected_record["scaffold_G"],
        "selected_generated": selected_record["generated_Gprime"],
        "selected_energy": selected_group["minimum_energy"],
        "selected_proposal_source":
            selected_record.get("proposal_source", "hotspot"),
        "proposal_records": len(records),
        "reachable_records": reachable_records,
        "unique_completed_structures": len(groups),
        "proposal_seconds": summary["proposal_seconds"],
        "linearfold_requested":
            int(summary.get("linearfold_requested", "0")),
        "linearfold_scaffold_added":
            int(summary.get("linearfold_scaffold_added", "0")),
        "linearfold_energy": summary.get("linearfold_energy", ""),
        "linearfold_seconds": float(summary.get("linearfold_seconds", "0")),
        "viterbi_seconds": summary["viterbi_seconds"],
        "total_seconds": summary["total_seconds"],
        # Records whose CASE3 depth cap removed candidates; blank if the
        # prediction tool predates the field.
        "case3_truncated_records": summary.get("case3_truncated_records", ""),
        "nonempty_scaffold_records": nonempty_scaffolds,
        "any_nonempty_scaffold_subset_of_reference":
            int(any_nonempty_scaffold_subset),
        "best_scaffold_pair_recall": best_scaffold_recall,
        "best_compatible_scaffold_pair_recall":
            best_compatible_scaffold_recall,
        "best_completed_f1": best_f1,
        "selection_f1_regret": best_f1 - selected_accuracy.all_pairs.f1,
        "groups_with_energy_spread": sum(spread > 0.0 for spread in spreads),
        "maximum_energy_spread": max(spreads, default=0.0),
        **accuracy_fields(selected_accuracy),
    }


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
        tp = sum(int(row[f"{family}_tp"]) for row in rows)
        fp = sum(int(row[f"{family}_fp"]) for row in rows)
        fn = sum(int(row[f"{family}_fn"]) for row in rows)
        score = PairScore(tp, fp, fn)
        out[f"micro_{family}_tp"] = tp
        out[f"micro_{family}_fp"] = fp
        out[f"micro_{family}_fn"] = fn
        out[f"micro_{family}_ppv"] = score.ppv
        out[f"micro_{family}_sensitivity"] = score.sensitivity
        out[f"micro_{family}_f1"] = score.f1
    out["exact_rate"] = sum(int(row["exact"]) for row in rows) / count
    out["mean_best_completed_f1"] = sum(
        float(row["best_completed_f1"]) for row in rows) / count
    out["mean_selection_f1_regret"] = sum(
        float(row["selection_f1_regret"]) for row in rows) / count
    out["mean_total_seconds"] = sum(
        float(row["total_seconds"]) for row in rows) / count
    out["linearfold_scaffold_added_cases"] = sum(
        int(row.get("linearfold_scaffold_added", 0)) for row in rows)
    out["linearfold_selected_cases"] = sum(
        row.get("selected_proposal_source") == "linearfold_mfe"
        for row in rows)
    return out


def write_rows(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty result")
    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(file_descriptor, "w", newline="") as handle:
            fieldnames: list[str] = []
            seen: set[str] = set()
            for row in rows:
                for field in row:
                    if field not in seen:
                        seen.add(field)
                        fieldnames.append(field)
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t")
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


def load_resumable_rows(
    path: Path,
    expected_provenance: dict[str, object],
    requested_ids: Sequence[str],
) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_id: dict[str, dict[str, object]] = {}
    requested = set(requested_ids)
    for row in rows:
        case_id = row.get("pkb_short", "")
        if not case_id or case_id in by_id:
            raise ValueError(f"resume file has missing/duplicate ID: {case_id}")
        if case_id not in requested:
            raise ValueError(
                f"resume file contains unrequested ID: {case_id}")
        for key, value in expected_provenance.items():
            if row.get(key) != str(value):
                raise ValueError(
                    f"resume provenance mismatch for {case_id}: {key}")
        by_id[case_id] = dict(row)
    return by_id


def status_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    counts = Counter(str(row.get("status", "unknown")) for row in rows)
    return {
        "completed_cases": counts["ok"],
        "timeout_cases": counts["timeout"],
        "tool_failure_cases": counts["tool_failure"],
        "parse_failure_cases": counts["parse_failure"],
        "invalid_input_cases": counts["invalid_input"],
        "failed_cases": len(rows) - counts["ok"],
        "status_counts": ",".join(
            f"{status}={count}" for status, count in sorted(counts.items())),
    }


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--ids")
    selection.add_argument("--ids-file", type=Path)
    selection.add_argument("--all", action="store_true")
    parser.add_argument("--tool", type=Path,
                        default=Path("build/scfg2_hotspot_viterbi_tool"))
    parser.add_argument(
        "--parallel-tool", type=Path,
        default=Path("test/tools/run_scfg2_parallel_hotspot_viterbi.py"))
    parser.add_argument(
        "--proposal-tool", type=Path,
        default=Path("build/scfg2_hotspot_proposal_tool"))
    parser.add_argument(
        "--worker-tool", type=Path,
        default=Path("build/scfg2_fixed_scaffold_viterbi_tool"))
    parser.add_argument(
        "--candidate-jobs", type=int, default=0,
        help="worker processes per RNA; zero uses the sequential C++ tool")
    parser.add_argument(
        "--linearfold-bin", type=Path,
        help="optionally add one pk-free LinearFold MFE scaffold")
    parser.add_argument("--linearfold-beam-width", type=int, default=100)
    parser.add_argument("--param-file", type=Path,
                        default=Path("params/rna_Turner04.par"))
    parser.add_argument("--beam-width", type=int, default=200)
    parser.add_argument("--anchor-cap", type=int, default=8)
    parser.add_argument("--hotspot-width", type=int, default=16)
    parser.add_argument("--passes", type=int, default=8)
    parser.add_argument("--max-results", type=int, default=16)
    parser.add_argument("--diversity", choices=("none", "deterministic"),
                        default="deterministic")
    parser.add_argument(
        "--jobs", type=int, default=1,
        help="independent RNA subprocesses to run concurrently")
    parser.add_argument(
        "--skip-invalid", action="store_true",
        help="report and skip dataset rows that fail input validation")
    parser.add_argument(
        "--resume", action="store_true",
        help="reuse output rows only when dataset/tool/config hashes match")
    parser.add_argument(
        "--timeout-seconds", type=float, default=3600.0,
        help="wall-time limit per RNA; zero disables the limit")
    parser.add_argument(
        "--artifacts-dir", type=Path,
        help="per-case stdout/stderr directory (default: beside output)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args(argv)

    dataset = load_dataset(args.dataset)
    try:
        ids = select_case_ids(
            dataset, args.ids, args.ids_file, args.all)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if (
        args.jobs <= 0
        or args.candidate_jobs < 0
        or args.linearfold_beam_width <= 0
        or args.timeout_seconds < 0
    ):
        parser.error(
            "--jobs must be positive; candidate jobs/timeout nonnegative")
    if args.linearfold_bin is not None and args.candidate_jobs == 0:
        parser.error("--linearfold-bin requires --candidate-jobs > 0")
    artifacts_dir = (
        args.artifacts_dir
        if args.artifacts_dir is not None
        else args.output.parent / f"{args.output.stem}.artifacts"
    )
    try:
        parallel = args.candidate_jobs > 0
        provenance: dict[str, object] = {
            "dataset_sha256": sha256_file(args.dataset),
            "tool_sha256": sha256_file(
                args.parallel_tool if parallel else args.tool),
            "param_sha256": sha256_file(args.param_file),
            "beam_width": args.beam_width,
            "anchor_cap": args.anchor_cap,
            "hotspot_width": args.hotspot_width,
            "passes": args.passes,
            "max_results": args.max_results,
            "diversity": args.diversity,
            "jobs": args.jobs,
            "candidate_jobs": args.candidate_jobs,
            "linearfold_beam_width": args.linearfold_beam_width,
            "timeout_seconds": args.timeout_seconds,
        }
        if parallel:
            provenance["proposal_tool_sha256"] = sha256_file(
                args.proposal_tool)
            provenance["worker_tool_sha256"] = sha256_file(args.worker_tool)
        if args.linearfold_bin is not None:
            provenance["linearfold_bin_sha256"] = sha256_file(
                args.linearfold_bin)
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
        print(
            f"[pseudobase-viterbi] skip {case_id}: {reason}",
            file=sys.stderr,
            flush=True,
        )
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
            f"[pseudobase-viterbi] resume {reused_ok} completed"
            f" cases; {len(pending_ids)} pending",
            file=sys.stderr,
            flush=True,
        )

    def run_one(case_id: str) -> dict[str, object]:
        row = dataset[case_id]
        artifact_root = artifacts_dir / safe_case_name(case_id)
        stdout_path = artifact_root / "stdout.txt"
        stderr_path = artifact_root / "stderr.txt"
        metrics_path = artifact_root / "peak_rss_kb.txt"
        common: dict[str, object] = {
            **provenance,
            "pkb_short": row["pkb_short"],
            "length": row["length"],
            "rna_type": row["rna_type"],
            "reference_structure": row["dotbracket"],
            "host": socket.gethostname(),
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
        }
        print(
            f"[pseudobase-viterbi] start {case_id} n={row['length']}",
            file=sys.stderr,
            flush=True,
        )
        prediction_arguments = [
            row["sequence"].upper(),
            str(args.param_file),
            str(args.beam_width),
            str(args.anchor_cap),
            str(args.hotspot_width),
            str(args.passes),
            str(args.max_results),
            args.diversity,
        ]
        if args.candidate_jobs > 0:
            command = [
                sys.executable,
                str(args.parallel_tool),
                "--proposal-tool", str(args.proposal_tool),
                "--worker-tool", str(args.worker_tool),
                "--candidate-jobs", str(args.candidate_jobs),
            ]
            if args.linearfold_bin is not None:
                command.extend([
                    "--linearfold-bin", str(args.linearfold_bin),
                    "--linearfold-beam-width",
                    str(args.linearfold_beam_width),
                ])
            command.extend(prediction_arguments)
        else:
            command = [str(args.tool), *prediction_arguments]
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
            result = evaluate_case(
                row, parse_prediction_output(completed.stdout))
        except (KeyError, TypeError, ValueError) as error:
            return {
                **common,
                "status": "parse_failure",
                **execution,
                "error": " ".join(str(error).split()),
            }
        result = {
            **common,
            "status": "ok",
            **execution,
            "error": "",
            **result,
        }
        print(
            f"[pseudobase-viterbi] done {case_id}"
            f" F1={float(result['all_f1']):.6g}",
            file=sys.stderr,
            flush=True,
        )
        return result

    batch_begin = time.perf_counter()
    if args.jobs == 1:
        for case_id in pending_ids:
            results_by_id[case_id] = run_one(case_id)
            write_rows(
                args.output,
                [results_by_id[key] for key in ids
                 if key in results_by_id],
            )
    else:
        # The C++ compatibility owner has process-global state. Each task
        # launches one isolated tool process; no BeamSession is shared across
        # threads or processes.
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.jobs) as executor:
            futures = {
                executor.submit(run_one, case_id): case_id
                for case_id in pending_ids
            }
            for future in concurrent.futures.as_completed(futures):
                case_id = futures[future]
                results_by_id[case_id] = future.result()
                # Preserve completed cases if a later subprocess fails.
                write_rows(
                    args.output,
                    [results_by_id[key] for key in ids
                     if key in results_by_id],
                )
    batch_wall_seconds = time.perf_counter() - batch_begin
    results = [results_by_id[case_id] for case_id in ids]
    write_rows(args.output, results)
    summary_path = (
        args.summary_output
        if args.summary_output is not None
        else args.output.with_suffix(".summary.tsv")
    )
    successful_results = [
        row for row in results if row.get("status") == "ok"]
    summary: dict[str, object] = {
        "beam_width": args.beam_width,
        "anchor_cap": args.anchor_cap,
        "hotspot_width": args.hotspot_width,
        "passes": args.passes,
        "max_results": args.max_results,
        "diversity": args.diversity,
        "jobs": args.jobs,
        "candidate_jobs": args.candidate_jobs,
        "linearfold_requested": int(args.linearfold_bin is not None),
        "linearfold_beam_width": args.linearfold_beam_width,
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
        f"[pseudobase-viterbi] wrote {len(results)} cases to {args.output}"
        f" and summary to {summary_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
