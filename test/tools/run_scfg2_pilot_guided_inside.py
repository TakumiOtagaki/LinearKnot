#!/usr/bin/env python3
"""Compare baseline and b^q-pilot guided inside, serially on one CPU."""

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
    parser.add_argument("--widths", default="20,100")
    parser.add_argument("--outer-factor", type=int, default=4)
    parser.add_argument("--pilot-exponents", default="2/3",
                        help="comma-separated subset of 0.5,2/3 (default: 2/3)")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--empty-controls", action="store_true")
    parser.add_argument("--include-witnesses", action="store_true")
    parser.add_argument("--skip-exact", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--closure", choices=("weak", "legacy"), default="weak")
    parser.add_argument("--params", type=Path, default=Path("params/rna_Turner04.par"))
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--memory-mib", type=int, default=2048)
    args = parser.parse_args()
    if args.outer_factor < 1 or args.timeout < 1 or args.memory_mib < 256 or args.repeats < 1:
        parser.error("positive factor/repeats/timeout and at least 256 MiB are required")
    q_aliases = {"0.5": "0.5", "1/2": "0.5", "2/3": "2/3"}
    try:
        exponents = [q_aliases[q] for q in args.pilot_exponents.split(",")]
    except KeyError:
        parser.error("supported pilot exponents: 0.5,2/3")
    if len(set(exponents)) != len(exponents):
        parser.error("duplicate pilot exponent")
    repo = Path(__file__).resolve().parents[2]
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if not args.skip_build:
        with (out / "build.log").open("w") as log:
            subprocess.run(["cmake", "-S", str(repo), "-B", str(repo / "build")],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
            subprocess.run(["cmake", "--build", str(repo / "build"), "-j8", "--target",
                            "scfg2_pilot_guided_inside_probe"], check=True,
                           stdout=log, stderr=subprocess.STDOUT)
    params = args.params.resolve()
    inputs = out / "inputs"
    inputs.mkdir(exist_ok=True)
    cases, names = [], set()
    for path in args.input:
        if path.stem in names:
            parser.error(f"duplicate input: {path.stem}")
        names.add(path.stem)
        copy = inputs / path.name
        shutil.copyfile(path.resolve(), copy)
        cases.append((path.stem + "_given", copy, "given"))
        if args.empty_controls:
            cases.append((path.stem + "_empty", copy, "empty"))
    if args.include_witnesses:
        toy = inputs / "toy12.txt"
        toy.write_text(">toy12\nGGAAACGAAACC\n((...)(...))\n")
        cases.extend((("toy12_given", toy, "given"), ("toy12_empty", toy, "empty")))
        for n in (64, 128):
            path = inputs / f"gc{n}.txt"
            path.write_text(f">gc{n}\n" + "GC" * (n // 2) + "\n" + "." * n + "\n")
            cases.append((f"gc{n}", path, "empty"))
        for depth in (32, 64):
            path = inputs / f"nested{depth}.txt"
            path.write_text(f">nested{depth}\n" + "G" * depth + "AAAA" + "C" * depth +
                            "\n" + "(" * depth + "...." + ")" * depth + "\n")
            cases.append((f"nested{depth}", path, "given"))
    if not cases:
        parser.error("provide --input or --include-witnesses")
    cpu = min(os.sched_getaffinity(0))

    def limits():
        os.sched_setaffinity(0, {cpu})
        cap = args.memory_mib * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    snapshot = out / "source"
    snapshot.mkdir(exist_ok=True)
    source_paths = [
        "src/scfg2/replay/w_final_vp_wmbp_online.hh",
        "src/scfg2/deductions/generate_exact_basic.hh",
        "src/scfg2/engine/basic_schedule.hh",
        "src/scfg2/engine/legacy_span_support.hh",
        "test/tools/scfg2_pilot_guided_inside_probe.cc",
        "test/tools/scfg2_pilot_guided_inside.hh",
        "test/tools/scfg2_bounded_outer.hh",
        "test/tools/scfg2_outer_probe_common.hh",
        "test/tools/run_scfg2_pilot_guided_inside.py",
    ]
    hashes = {}
    for name in source_paths:
        path = repo / name
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        target = snapshot / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    (snapshot / "tracked_diff.patch").write_bytes(
        subprocess.check_output(["git", "diff", "HEAD"], cwd=repo))
    manifest = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "sources": hashes, "params": str(params),
        "param_sha256": hashlib.sha256(params.read_bytes()).hexdigest(),
        "cpu": cpu, "widths": args.widths, "outer_factor": args.outer_factor,
        "pilot_exponents": exponents, "repeats": args.repeats,
        "closure": args.closure, "exact_scoring": not args.skip_exact,
        "timeout_seconds": args.timeout, "memory_mib": args.memory_mib,
        "scope": "inside only; public ranking guided, private factors unchanged; shared scaffold initialization timed separately; outer scratch cleanup included; whole helper call timed",
        "cases": [],
    }
    summaries, families, failures = [], [], []
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for name, path, mode in cases:
        print(f"start {name}", flush=True)
        prefix = out / name
        command = ["/usr/bin/time", "-f", "%e\t%M", "-o", str(prefix) + ".time.tsv",
                   str(repo / "build/scfg2_pilot_guided_inside_probe"), str(path), mode,
                   str(prefix), args.widths, "no-exact" if args.skip_exact else "exact",
                   str(params), str(args.outer_factor), args.closure,
                   ",".join(exponents), str(args.repeats)]
        record = {"name": name, "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "mode": mode, "command": command}
        with Path(str(prefix) + ".log").open("w") as log:
            process = subprocess.Popen(command, cwd=repo, env=environment, preexec_fn=limits,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                record["timeout"] = True
        record["exit_code"] = process.returncode
        if process.returncode:
            failures.append(name)
            print(f"failed {name}; see {prefix}.log", flush=True)
        else:
            timing = Path(str(prefix) + ".time.tsv").read_text().strip().split("\t")
            for suffix, rows in (("summary", summaries), ("families", families)):
                with Path(str(prefix) + f".{suffix}.tsv").open() as stream:
                    for row in csv.DictReader(stream, delimiter="\t"):
                        row.update(case=name)
                        if suffix == "summary":
                            row.update(wall_seconds=timing[0], peak_rss_kib=timing[1])
                        rows.append(row)
            print(f"done {name}: {timing[0]} s, {int(timing[1]) / 1024:.1f} MiB", flush=True)
        manifest["cases"].append(record)
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for name, rows in (("summary", summaries), ("families", families)):
        if rows:
            with (out / f"{name}.tsv").open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
    if failures:
        raise SystemExit("failed cases: " + ", ".join(failures))


if __name__ == "__main__":
    main()
