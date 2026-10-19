#!/usr/bin/env python3
"""Run one hotspot proposal set with fixed-scaffold Viterbi worker processes."""

from __future__ import annotations

import argparse
import concurrent.futures
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence


sys.path.insert(0, str(Path(__file__).resolve().parent))

from scfg2_pair_accuracy import parse_dot_bracket  # noqa: E402


def parse_key_values(fields: Sequence[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition("=")
        if not separator:
            raise ValueError(f"malformed key-value field: {field}")
        out[key] = value
    return out


def parse_proposals(text: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    summary: dict[str, str] | None = None
    header: list[str] | None = None
    proposals: list[dict[str, str]] = []
    for raw in text.splitlines():
        fields = raw.split("\t")
        if fields[0] == "SUMMARY":
            summary = parse_key_values(fields[1:])
        elif fields[0] == "PROPOSAL_HEADER":
            header = fields[1:]
        elif fields[0] == "PROPOSAL":
            if header is None or len(fields) - 1 != len(header):
                raise ValueError("PROPOSAL does not match PROPOSAL_HEADER")
            proposals.append(dict(zip(header, fields[1:])))
    if summary is None or summary.get("ok") != "1":
        raise ValueError(
            "proposal failed: "
            + (summary or {}).get("error", "missing SUMMARY"))
    for expected, proposal in enumerate(proposals):
        if int(proposal["index"]) != expected:
            raise ValueError("proposal indices are not contiguous")
    return summary, proposals


def parse_worker(text: str) -> dict[str, str]:
    summaries = [
        parse_key_values(line.split("\t")[1:])
        for line in text.splitlines()
        if line.startswith("SUMMARY\t")
    ]
    if len(summaries) != 1 or summaries[0].get("ok") != "1":
        raise ValueError(
            "worker failed: "
            + (summaries[0] if summaries else {}).get(
                "error", "missing SUMMARY"))
    return summaries[0]


def add_linearfold_scaffold(
    proposals: list[dict[str, str]],
    binary: Path,
    sequence: str,
    beam_width: int,
) -> tuple[bool, str, float]:
    begin = time.perf_counter()
    completed = subprocess.run(
        [str(binary), sequence, str(beam_width)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    seconds = time.perf_counter() - begin
    lines = [line.strip() for line in completed.stdout.splitlines()
             if line.strip()]
    if not lines:
        raise ValueError("LinearFold output is empty")
    scaffold = lines[0]
    if len(scaffold) != len(sequence):
        raise ValueError("LinearFold scaffold length differs from sequence")
    if any(symbol not in ".()" for symbol in scaffold):
        raise ValueError("LinearFold scaffold is not pk-free dot-bracket")
    scaffold_pairs = parse_dot_bracket(scaffold)
    existing = {
        parse_dot_bracket(proposal["scaffold_G"])
        for proposal in proposals
    }
    energy = lines[1] if len(lines) > 1 else ""
    if scaffold_pairs in existing:
        return False, energy, seconds
    proposals.append({
        "index": str(len(proposals)),
        "scaffold_pairs": str(len(scaffold_pairs)),
        # The hotspot heuristic is not defined for an external MFE seed and
        # does not participate in final energy selection.
        "heuristic_score": "0",
        "pass": str(max(
            (int(proposal["pass"]) for proposal in proposals),
            default=-1,
        ) + 1),
        "origin": "0",
        "scaffold_G": scaffold,
        "proposal_source": "linearfold_mfe",
    })
    return True, energy, seconds


def prediction_groups(
    proposals: Sequence[dict[str, str]],
    workers: Sequence[dict[str, str]],
) -> tuple[list[dict[str, object]], int | None]:
    grouped: dict[tuple[tuple[int, int], ...], list[int]] = {}
    for index, worker in enumerate(workers):
        if worker["reachable"] != "1":
            continue
        key = tuple(sorted(parse_dot_bracket(worker["union_R"])))
        grouped.setdefault(key, []).append(index)

    groups: list[dict[str, object]] = []
    for union_key in sorted(grouped):
        indices = grouped[union_key]
        selected = min(
            indices,
            key=lambda index: (
                float(workers[index]["energy"]),
                tuple(sorted(parse_dot_bracket(
                    proposals[index]["scaffold_G"]))),
                int(proposals[index]["pass"]),
                int(proposals[index]["origin"]),
                index,
            ),
        )
        energies = [float(workers[index]["energy"]) for index in indices]
        groups.append({
            "record_count": len(indices),
            "minimum_energy": min(energies),
            "maximum_energy": max(energies),
            "energy_spread": max(energies) - min(energies),
            "selected_record": selected,
            "union_key": union_key,
        })
    if not groups:
        return groups, None
    selected_group = min(
        range(len(groups)),
        key=lambda index: (
            float(groups[index]["minimum_energy"]),
            groups[index]["union_key"],
        ),
    )
    return groups, selected_group


def format_number(value: object) -> str:
    return f"{float(value):.12g}"


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposal-tool", type=Path, required=True)
    parser.add_argument("--worker-tool", type=Path, required=True)
    parser.add_argument("--candidate-jobs", type=int, required=True)
    parser.add_argument("--linearfold-bin", type=Path)
    parser.add_argument("--linearfold-beam-width", type=int, default=100)
    parser.add_argument("sequence")
    parser.add_argument("param_file")
    parser.add_argument("beam_width", type=int)
    parser.add_argument("anchor_cap", type=int)
    parser.add_argument("hotspot_width", type=int)
    parser.add_argument("passes", type=int)
    parser.add_argument("max_results", type=int)
    parser.add_argument("diversity", choices=("none", "deterministic"))
    args = parser.parse_args(argv)
    if args.candidate_jobs <= 0:
        parser.error("--candidate-jobs must be positive")
    if args.linearfold_beam_width <= 0:
        parser.error("--linearfold-beam-width must be positive")

    total_begin = time.perf_counter()
    proposal_command = [
        str(args.proposal_tool),
        args.sequence,
        str(args.anchor_cap),
        str(args.hotspot_width),
        str(args.passes),
        str(args.max_results),
        args.diversity,
    ]
    try:
        proposed = subprocess.run(
            proposal_command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proposal_summary, proposals = parse_proposals(proposed.stdout)
        for proposal in proposals:
            proposal["proposal_source"] = "hotspot"
        linearfold_added = False
        linearfold_energy = ""
        linearfold_seconds = 0.0
        if args.linearfold_bin is not None:
            (
                linearfold_added,
                linearfold_energy,
                linearfold_seconds,
            ) = add_linearfold_scaffold(
                proposals,
                args.linearfold_bin,
                args.sequence,
                args.linearfold_beam_width,
            )
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(
            "SUMMARY\tok=0\treachable=0\trecords=0\tgroups=0"
            f"\terror={str(error).replace(chr(9), ' ')}")
        return 1

    def run_worker(proposal: dict[str, str]) -> dict[str, str]:
        command = [
            str(args.worker_tool),
            args.sequence,
            args.param_file,
            str(args.beam_width),
            proposal["scaffold_G"],
        ]
        completed = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return parse_worker(completed.stdout)

    workers_by_index: dict[int, dict[str, str]] = {}
    try:
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.candidate_jobs) as executor:
            futures = {
                executor.submit(run_worker, proposal): int(proposal["index"])
                for proposal in proposals
            }
            for future in concurrent.futures.as_completed(futures):
                workers_by_index[futures[future]] = future.result()
        workers = [
            workers_by_index[index] for index in range(len(proposals))]
        groups, selected_group = prediction_groups(proposals, workers)
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(
            "SUMMARY\tok=0\treachable=0"
            f"\trecords={len(workers_by_index)}\tgroups=0"
            f"\tproposal_seconds={proposal_summary['proposal_seconds']}"
            f"\ttotal_seconds={time.perf_counter() - total_begin:.12g}"
            f"\terror={str(error).replace(chr(9), ' ')}")
        return 1

    viterbi_seconds = sum(float(worker["seconds"]) for worker in workers)
    truncated = [worker.get("case3_truncated_endpoints") for worker in workers]
    # Blank, not zero, when a worker predates the field: zero would claim the
    # CASE3 depth cap was never active.
    case3_truncated_records = (
        "" if any(value is None for value in truncated)
        else str(sum(int(value) > 0 for value in truncated)))
    total_seconds = time.perf_counter() - total_begin
    print(
        "SUMMARY\tok=1"
        f"\treachable={int(selected_group is not None)}"
        f"\tlength={len(args.sequence)}"
        f"\tbeam_width={args.beam_width}"
        f"\tanchor_cap={args.anchor_cap}"
        f"\thotspot_width={args.hotspot_width}"
        f"\tpasses={args.passes}"
        f"\tmax_results={args.max_results}"
        f"\trecords={len(proposals)}"
        f"\tgroups={len(groups)}"
        f"\tproposal_seconds={proposal_summary['proposal_seconds']}"
        f"\tlinearfold_requested={int(args.linearfold_bin is not None)}"
        f"\tlinearfold_scaffold_added={int(linearfold_added)}"
        f"\tlinearfold_energy={linearfold_energy}"
        f"\tlinearfold_seconds={linearfold_seconds:.12g}"
        f"\tviterbi_seconds={viterbi_seconds:.12g}"
        f"\ttotal_seconds={total_seconds:.12g}"
        f"\tcase3_truncated_records={case3_truncated_records}"
        "\terror=")
    print(
        "RECORD_HEADER\tindex\treachable\tscaffold_pairs"
        "\theuristic_score\tpass\torigin\tproposal_source\tenergy\tseconds"
        "\tscaffold_G\tgenerated_Gprime\tunion_R")
    for index, (proposal, worker) in enumerate(zip(proposals, workers)):
        print(
            f"RECORD\t{index}\t{worker['reachable']}"
            f"\t{proposal['scaffold_pairs']}"
            f"\t{proposal['heuristic_score']}"
            f"\t{proposal['pass']}\t{proposal['origin']}"
            f"\t{proposal['proposal_source']}"
            f"\t{worker['energy']}\t{worker['seconds']}"
            f"\t{proposal['scaffold_G']}"
            f"\t{worker['generated_Gprime']}\t{worker['union_R']}")
    print(
        "GROUP_HEADER\tindex\trecord_count\tminimum_energy"
        "\tmaximum_energy\tenergy_spread\tselected_record\tselected")
    for index, group in enumerate(groups):
        print(
            f"GROUP\t{index}\t{group['record_count']}"
            f"\t{format_number(group['minimum_energy'])}"
            f"\t{format_number(group['maximum_energy'])}"
            f"\t{format_number(group['energy_spread'])}"
            f"\t{group['selected_record']}"
            f"\t{int(index == selected_group)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
