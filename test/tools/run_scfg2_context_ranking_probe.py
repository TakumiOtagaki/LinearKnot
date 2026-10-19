#!/usr/bin/env python3
"""Run context/outer ranking diagnostics serially; exact charts evaluate ranks."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--empty-controls", action="store_true")
    parser.add_argument("--include-witnesses", action="store_true")
    parser.add_argument("--time-scalar-passes", action="store_true")
    parser.add_argument("--outer-pilot", action="store_true",
                        help="estimate max outer from finite-beam sibling alphas")
    parser.add_argument("--pilot-widths", default="20")
    parser.add_argument("--edge-caps", default="0,4,16",
                        help="outer targets per parent and target family; 0=all")
    parser.add_argument("--reserve-factors", default="",
                        help="also restrict outer states to pilot + top factor*p pre-prune candidates")
    parser.add_argument("--span-factors", default="",
                        help="also keep factor*p novel outer candidates per family/span")
    parser.add_argument("--bounded-factors", default="",
                        help="bounded split proposals and sparse reverse schedule, factor*p novel states")
    parser.add_argument("--bounded-only", action="store_true",
                        help="skip exact charts/unrestricted diagnostics; requires bounded factors")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--memory-mib", type=int, default=2048)
    parser.add_argument("--widths", default="20,100,200")
    args = parser.parse_args()
    if args.timeout <= 0 or args.memory_mib < 256:
        parser.error("positive timeout and at least 256 MiB are required")
    if args.outer_pilot and args.time_scalar_passes:
        parser.error("--time-scalar-passes applies to the exact-context probe")
    if (args.bounded_factors or args.bounded_only) and not args.outer_pilot:
        parser.error("bounded outer options require --outer-pilot")
    if args.bounded_only and not args.bounded_factors:
        parser.error("--bounded-only requires --bounded-factors")
    repo = Path(__file__).resolve().parents[2]
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    build = out / "probe-build"
    build.mkdir(exist_ok=True)
    source = repo / "test/tools" / (
        "scfg2_outer_ranking_probe.cc" if args.outer_pilot else
        "scfg2_context_ranking_probe.cc")
    snapshot = build / source.name
    shutil.copyfile(source, snapshot)
    helper_hash = None
    if args.outer_pilot:
        helper = repo / "test/tools/scfg2_bounded_outer.hh"
        shutil.copyfile(helper, build / helper.name)
        common = repo / "test/tools/scfg2_outer_probe_common.hh"
        shutil.copyfile(common, build / common.name)
        helper_hash = hashlib.sha256(helper.read_bytes()).hexdigest()
    cmake = f'''cmake_minimum_required(VERSION 3.16)
project(context_ranking_probe LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
add_executable(context_ranking_probe "{snapshot}")
target_include_directories(context_ranking_probe PRIVATE
  "{repo / 'src'}" "{repo / 'test'}" "{repo / 'build/generated'}")
target_compile_options(context_ranking_probe PRIVATE -O2 -Wall -Wextra)
target_link_libraries(context_ranking_probe PRIVATE
  "{repo / 'build/libcparty_scfg2_exact.a'}" "{repo / 'build/libRNA.a'}")
'''
    (build / "CMakeLists.txt").write_text(cmake)
    with (out / "build.log").open("w") as log:
        subprocess.run(["cmake", "-S", str(build), "-B", str(build / "build")],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
        subprocess.run(["cmake", "--build", str(build / "build"), "-j8"],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
    executable = build / "build/context_ranking_probe"
    inputs_dir = out / "inputs"
    inputs_dir.mkdir(exist_ok=True)
    cases = []
    names = set()
    for path in args.input:
        name = path.stem
        if name in names:
            parser.error(f"duplicate input name: {name}")
        names.add(name)
        copy = inputs_dir / path.name
        shutil.copyfile(path.resolve(), copy)
        cases.append((name + "_given", copy, "given"))
        if args.empty_controls:
            cases.append((name + "_empty", copy, "empty"))
    if args.include_witnesses:
        toy = inputs_dir / "toy12.txt"
        toy.write_text(">toy12\nGGAAACGAAACC\n((...)(...))\n")
        cases.extend([("toy12_given", toy, "given"),
                      ("toy12_empty", toy, "empty")])
        for depth in (8, 32, 64, 128):
            path = inputs_dir / f"nested{depth}.txt"
            path.write_text(f">nested{depth}\n" + "G" * depth + "AAAA" +
                            "C" * depth + "\n" + "(" * depth + "...." +
                            ")" * depth + "\n")
            cases.append((f"nested{depth}", path, "given"))
    if not cases:
        parser.error("provide --input or --include-witnesses")
    cpu = min(os.sched_getaffinity(0))

    def limits():
        os.sched_setaffinity(0, {cpu})
        cap = args.memory_mib * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    manifest = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                           cwd=repo, text=True).strip(),
        "source_sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        "param_sha256": hashlib.sha256(
            (repo / "params/rna_Turner04.par").read_bytes()).hexdigest(),
        "timeout_seconds": args.timeout, "memory_mib": args.memory_mib,
        "cpu": cpu, "widths": args.widths, "cases": [],
        "time_scalar_passes": args.time_scalar_passes,
        "scope": ("bounded outer only; pilot plus bounded proposals, no exact schedule/chart or full-split validation"
                  if args.bounded_only else
                  "finite-beam pilot sibling alphas; exact enumeration and rank evaluation"
                  if args.outer_pilot else "offline immutable exact chart; no forward beam changes"),
        "outer_pilot": args.outer_pilot,
        "pilot_widths": args.pilot_widths if args.outer_pilot else None,
        "edge_caps": args.edge_caps if args.outer_pilot else None,
        "reserve_factors": args.reserve_factors if args.outer_pilot else None,
        "span_factors": args.span_factors if args.outer_pilot else None,
        "bounded_factors": args.bounded_factors if args.outer_pilot else None,
        "bounded_only": args.bounded_only,
        "bounded_helper_sha256": helper_hash,
        "outer_common_sha256": (hashlib.sha256(common.read_bytes()).hexdigest()
                                 if args.outer_pilot else None),
        "shared_kernel_sha256": {
            str(path): hashlib.sha256((repo / path).read_bytes()).hexdigest()
            for path in (
                Path("src/scfg2/deductions/generate_exact_basic.hh"),
                Path("src/scfg2/deductions/split_positions.hh"),
                Path("src/scfg2/engine/basic_schedule.hh"),
            )
        },
    }
    summaries = []
    rankings = []
    failures = []
    for name, path, mode in cases:
        print(f"start {name}", flush=True)
        prefix = out / name
        command = ["/usr/bin/time", "-f", "%e\t%M", "-o",
                   str(prefix) + ".time.tsv", str(executable), str(path),
                   mode, str(prefix), args.widths]
        if args.time_scalar_passes:
            command.append("time-scalar")
        if args.outer_pilot:
            command.extend([args.pilot_widths, args.edge_caps])
            if args.reserve_factors or args.span_factors or args.bounded_factors:
                command.append(args.reserve_factors)
            if args.span_factors or args.bounded_factors:
                command.append(args.span_factors)
            if args.bounded_factors:
                command.append(args.bounded_factors)
            if args.bounded_only:
                command.append("bounded-only")
        record = {"name": name, "input_sha256": hashlib.sha256(
            path.read_bytes()).hexdigest(), "mode": mode, "command": command}
        try:
            with Path(str(prefix) + ".log").open("w") as log:
                process = subprocess.Popen(command, cwd=repo, env=environment,
                                           preexec_fn=limits, stdout=log,
                                           stderr=subprocess.STDOUT,
                                           start_new_session=True)
                try:
                    process.wait(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
            record["exit_code"] = process.returncode
            if process.returncode:
                failures.append(name)
                print(f"failed {name}; see {prefix}.log", flush=True)
            else:
                with Path(str(prefix) + ".summary.tsv").open() as stream:
                    case_summaries = list(csv.DictReader(stream, delimiter="\t"))
                timing = Path(str(prefix) + ".time.tsv").read_text().strip().split("\t")
                for summary in case_summaries:
                    summary.update(case=name, wall_seconds=timing[0], peak_rss_kib=timing[1])
                    summaries.append(summary)
                with Path(str(prefix) + ".ranking.tsv").open() as stream:
                    for row in csv.DictReader(stream, delimiter="\t"):
                        row["case"] = name
                        rankings.append(row)
                print(f"done {name}: {timing[0]} s, {int(timing[1]) / 1024:.1f} MiB",
                      flush=True)
        except subprocess.TimeoutExpired:
            record["timeout"] = True
            failures.append(name)
            print(f"timeout {name}", flush=True)
        manifest["cases"].append(record)
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for filename, rows in (("summary.tsv", summaries), ("ranking.tsv", rankings)):
        if rows:
            with (out / filename).open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
    if failures:
        raise SystemExit("failed cases: " + ", ".join(failures))


if __name__ == "__main__":
    main()
