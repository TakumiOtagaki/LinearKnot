#!/usr/bin/env python3
"""Measure completed-WMB closure using an external, instrumented legacy TU.

The pinned oracle checkout and its executable are never modified. The only
legacy source copied to /tmp is part_func.cc; the probe reuses existing objects.
SCFG2 selects the public kernel policy; the independent legacy gate is confined
to this external audit probe.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
PREFIX = "WMB_AUDIT_JSON "


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True, **kwargs)


def write_log(path: Path, *chunks: str) -> None:
    path.write_text("\n".join(line.rstrip() for line in "".join(chunks).splitlines()) + "\n")


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"expected one instrumentation anchor: {old!r}")
    return source.replace(old, new, 1)


def instrument(source: str) -> str:
    source = replace_once(source, '#include <algorithm>', '#include <algorithm>\n#include <cstdlib>\n#include <iomanip>')
    source = replace_once(source, 'pf_t W_final_pf::hfold_pf(sparse_tree &tree) {', '''pf_t W_final_pf::hfold_pf(sparse_tree &tree) {
    // Independent definition check: the proposed gate must not reject a
    // closed span because of a sparse-tree predicate error.
    for (cand_pos_t i = 1; i <= n; ++i) {
        for (cand_pos_t j = i; j <= n; ++j) {
            bool closed_by_scan = true;
            for (cand_pos_t k = i; k <= j; ++k) {
                int partner = tree.tree[k].pair;
                if (partner > 0 && (partner < i || partner > j))
                    closed_by_scan = false;
            }
            if (closed_by_scan != tree.weakly_closed(i, j)) {
                std::cerr << "WMB_AUDIT_CLOSURE_MISMATCH " << i << " " << j << "\\n";
                std::abort();
            }
        }
    }''')
    source = replace_once(source, '    WMB[ij] = contributions;\n}', '''    // Experimental gate on completed WMB only. Partial WMBP is untouched.
    if (std::getenv("CPARTY_AUDIT_WMB_CLOSURE") && !tree.weakly_closed(i, j))
        contributions = 0;
    WMB[ij] = contributions;
}''')
    return replace_once(source, '    pf_t energy = to_Energy(W[n], n);', '''    pf_t energy = to_Energy(W[n], n);
    std::cerr << std::setprecision(17)
              << "WMB_AUDIT_JSON {\\"logZ\\":"
              << std::log(W[n]) + n * std::log(exp_params_->pf_scale)
              << ",\\"energy\\":" << energy
              << ",\\"rt_kcal\\":" << exp_params_->kT / 1000.0
              << ",\\"pf_scale\\":" << exp_params_->pf_scale
              << ",\\"psm_kcal\\":" << PSM_penalty / 100.0
              << ",\\"closure_spans_verified\\":" << n * (n + 1) / 2 << "}\\n";
    if (std::getenv("CPARTY_AUDIT_WMB_TRACE") && n == 14) {
        for (auto span : {std::pair<int, int>{2, 7}, {8, 13}, {2, 13}}) {
            int i = span.first, j = span.second;
            std::cerr << "WMB_AUDIT_STATE i=" << i << " j=" << j
                      << " closed=" << tree.weakly_closed(i, j)
                      << " WMBP=" << get_energy_WMBP(i, j)
                      << " WMB=" << get_energy_WMB(i, j)
                      << " WM=" << get_energy_WM(i, j)
                      << " WMp=" << get_energy_WMp(i, j) << "\\n";
        }
        // Root contribution of the illegal two-branch split and the one-chain
        // branch. Both include the closing-pair scale and exterior factor.
        pf_t closing = exp_Mbloop(1, n) * exp_params_->expMLclosing
                       * scale[2] * exp_Extloop(1, n);
        pf_t split = get_energy_WM(2, 7) * get_energy_WMp(8, 13) * closing;
        pf_t chain = get_energy_WMp(2, 13) * closing;
        pf_t target_chain = get_BE(3, 10, 3, 10, tree)
                            * get_energy_WMBP(2, 7) * get_energy_VP(8, 13)
                            * pow(expPB_penalty, 2) * expPSM_penalty
                            * expb_penalty * closing;
        std::cerr << "WMB_AUDIT_WITNESS_JSON {\\"illegal_target_weight\\":"
                  << split / scale[n] << ",\\"canonical_target_weight\\":"
                  << target_chain / scale[n] << ",\\"illegal_to_canonical_ratio\\":"
                  << split / target_chain << "}\\n";
        std::cerr << "WMB_AUDIT_TERMS split=" << split << " chain=" << chain
                  << " split_energy=" << to_Energy(split, n)
                  << " chain_energy=" << to_Energy(chain, n)
                  << " target_chain_energy=" << to_Energy(target_chain, n) << "\\n";
    }''')


def build_probe(checkout: Path, work: Path, out: Path) -> tuple[Path, Path]:
    sys.path.insert(0, str(ROOT / "test/oracle/legacy_cparty"))
    from oracle_config import verify_checkout
    verify_checkout(checkout)
    source_path = checkout / "src/part_func.cc"
    original = source_path.read_text()
    modified = instrument(original)
    (work / "part_func.cc").write_text(modified)
    scfg_source = (ROOT / "test/scfg2_restricted_partition_function_tool.cc").read_text()
    scfg_source = replace_once(scfg_source,
        '    scfg2::replay::WFinalLogWeightContext ctx(pf, policy);',
        '    scfg2::replay::WFinalLogWeightContext ctx(pf, std::getenv("CPARTY_AUDIT_WMB_CLOSURE") '
        '? scfg2::WmbClosurePolicy::RequireWeakClosure : scfg2::WmbClosurePolicy::LegacyCompatible);')
    (work / "scfg2_probe.cc").write_text(scfg_source)
    (out / "legacy_instrumentation.patch").write_text("".join(difflib.unified_diff(
        original.splitlines(keepends=True), modified.splitlines(keepends=True),
        fromfile="legacy/src/part_func.cc", tofile="external-probe/part_func.cc", n=0)))
    build = checkout / "build"
    link_args = shlex.split((build / "CMakeFiles/CParty.dir/link.txt").read_text())
    objects = [build / arg for arg in link_args if arg.endswith(".o") and not arg.endswith("/part_func.cc.o")]
    if len(objects) != 12:
        raise ValueError(f"unexpected legacy object inventory: {len(objects)}")
    quote = lambda p: '"' + str(p).replace('\\', '/') + '"'
    cmake = '\n'.join([
        'cmake_minimum_required(VERSION 3.16)',
        'project(wmb_closure_audit LANGUAGES CXX)',
        'set(CMAKE_CXX_STANDARD 17)',
        'add_executable(wmb_closure_probe ' + quote(work / 'part_func.cc') + ' ' + ' '.join(map(quote, objects)) + ')',
        'target_include_directories(wmb_closure_probe PRIVATE ' + quote(checkout / 'src') + ')',
        'target_compile_options(wmb_closure_probe PRIVATE -Wall -Wextra -O3 -flto)',
        'target_link_options(wmb_closure_probe PRIVATE -flto)',
        'target_link_libraries(wmb_closure_probe PRIVATE ' + quote(build / 'libRNA.a') + ')',
        'add_executable(scfg2_closure_probe ' + quote(work / 'scfg2_probe.cc') + ')',
        'target_include_directories(scfg2_closure_probe PRIVATE ' +
        ' '.join(map(quote, [ROOT / 'src', ROOT / 'test', ROOT / 'build/generated'])) + ')',
        'target_compile_options(scfg2_closure_probe PRIVATE -Wall -Wextra -O3 -flto)',
        'target_link_options(scfg2_closure_probe PRIVATE -flto)',
        'target_link_libraries(scfg2_closure_probe PRIVATE ' + quote(ROOT / 'build/libcparty_scfg2_exact.a') + ' ' + quote(ROOT / 'build/libRNA.a') + ')',
        '',
    ])
    (work / "CMakeLists.txt").write_text(cmake)
    result = run(["cmake", "-S", str(work), "-B", str(work / "build")])
    (out / "configure.log").write_text(result.stdout + result.stderr)
    result = run(["cmake", "--build", str(work / "build"), "-j8"])
    (out / "build.log").write_text(result.stdout + result.stderr)
    metadata = {
        "oracle_manifest_sha256": hashlib.sha256((checkout.parent / "manifest.json").read_bytes()).hexdigest(),
        "oracle_commit": run(["git", "-C", str(checkout), "rev-parse", "HEAD"]).stdout.strip(),
        "oracle_diff_sha256": hashlib.sha256(run(["git", "-C", str(checkout), "diff", "HEAD"]).stdout.encode()).hexdigest(),
        "legacy_part_func_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "legacy_executable_sha256": hashlib.sha256((build / "CParty").read_bytes()).hexdigest(),
        "instrumented_part_func_sha256": hashlib.sha256(modified.encode()).hexdigest(),
        "scfg2_commit": run(["git", "-C", str(ROOT), "rev-parse", "HEAD"]).stdout.strip(),
        "scfg2_source_sha256": {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in ("src/scfg2/grammar/model_policy.hh",
                         "src/scfg2/deductions/exact_context.hh",
                         "src/scfg2/deductions/generate_exact_basic.hh",
                         "src/scfg2/replay/w_final_exact_inside.hh",
                         "test/scfg2_restricted_partition_function_tool.cc")},
        "scfg2_probe_sha256": hashlib.sha256(scfg_source.encode()).hexdigest(),
        "scfg2_library_sha256": hashlib.sha256((ROOT / 'build/libcparty_scfg2_exact.a').read_bytes()).hexdigest(),
        "parameter_file_sha256": hashlib.sha256((ROOT / 'params/rna_Turner04.par').read_bytes()).hexdigest(),
        "external_work_dir": str(work),
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return work / "build/wmb_closure_probe", work / "build/scfg2_closure_probe"


def probe(binary: Path, sequence: str, scaffold: str, closed: bool, trace: bool = False):
    env = os.environ.copy()
    for key in ("CPARTY_AUDIT_WMB_CLOSURE", "CPARTY_AUDIT_WMB_TRACE", "SCFG_DENSITY2_VIEW"):
        env.pop(key, None)
    if closed:
        env["CPARTY_AUDIT_WMB_CLOSURE"] = "1"
    if trace:
        env["CPARTY_AUDIT_WMB_TRACE"] = "1"
    args = [str(binary), "-s", "0", "-d2", "-P", str(ROOT / "params/rna_Turner04.par"), "-r", scaffold, sequence]
    result = run(args, cwd=ROOT, env=env)
    records = [json.loads(line[len(PREFIX):]) for line in result.stderr.splitlines() if line.startswith(PREFIX)]
    if len(records) != 1 or not all(math.isfinite(records[0][key]) for key in ("logZ", "energy")):
        raise ValueError(f"invalid precision record: {result.stderr}")
    return records[0], result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle-checkout", type=Path, default=ROOT / "build/legacy_cparty_oracle/src")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--count", type=int, default=100, help="cases per length")
    parser.add_argument("--lengths", default="20,30,50,80,100")
    args = parser.parse_args()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="cparty_wmb_audit_", dir="/tmp"))
    binary, scfg_binary = build_probe(args.oracle_checkout.resolve(), work, out)
    fixed = ("witness14", "GCGCGCGCGCGCGC", "(.(......)...)")
    cases = [fixed]
    for length in map(int, args.lengths.split(",")):
        generated = run([sys.executable, str(ROOT / "test/tools/gen_random_seq_struct.py"),
                         "--seed", str(args.seed), "--count", str(args.count),
                         "--length", str(length), "--prefix", f"rand{length}"]).stdout
        cases.extend(tuple(line.split("\t")) for line in generated.splitlines()[1:] if line)
    rows = []
    witness_terms = {}
    for case_id, sequence, scaffold in cases:
        before, result_before = probe(binary, sequence, scaffold, False, case_id == fixed[0])
        after, result_after = probe(binary, sequence, scaffold, True, case_id == fixed[0])
        scfg_energies = []
        for closed in (False, True):
            env = os.environ.copy()
            env.pop("CPARTY_AUDIT_WMB_CLOSURE", None)
            env.pop("SCFG_DENSITY2_VIEW", None)
            if closed:
                env["CPARTY_AUDIT_WMB_CLOSURE"] = "1"
            scfg_energies.append(float(run([str(scfg_binary), sequence, scaffold,
                str(ROOT / "params/rna_Turner04.par")], cwd=ROOT, env=env).stdout))
        if case_id == fixed[0]:
            original = run([str(args.oracle_checkout.resolve() / "build/CParty"), "-s", "0", "-d2", "-P",
                            str(ROOT / "params/rna_Turner04.par"), "-r", scaffold, sequence], cwd=ROOT)
            if original.stdout != result_before.stdout:
                raise ValueError("unguarded probe differs from original legacy stdout")
            write_log(out / "witness_legacy.log", result_before.stdout, result_before.stderr)
            write_log(out / "witness_completed_wmb_guard.log", result_after.stdout, result_after.stderr)
            term_prefix = "WMB_AUDIT_WITNESS_JSON "
            witness_terms = json.loads(next(line[len(term_prefix):] for line in result_before.stderr.splitlines()
                                           if line.startswith(term_prefix)))
        rows.append({"case_id": case_id, "n": len(sequence), "sequence": sequence, "scaffold": scaffold,
                     "legacy_logZ": before["logZ"], "guard_logZ": after["logZ"],
                     "legacy_energy": before["energy"], "guard_energy": after["energy"],
                     "energy_increase": after["energy"] - before["energy"],
                     "removed_fraction": -math.expm1(after["logZ"] - before["logZ"]),
                     "scfg2_legacy_energy": scfg_energies[0], "scfg2_guard_energy": scfg_energies[1],
                     "scfg2_legacy_abs_error": abs(scfg_energies[0] - before["energy"]),
                     "scfg2_guard_abs_error": abs(scfg_energies[1] - after["energy"]),
                     "closure_spans_verified": before["closure_spans_verified"]})
        if len(rows) % 100 == 0:
            print(f"audited {len(rows)}/{len(cases)} cases", flush=True)
    with (out / "cases.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    random_rows = rows[1:]
    worst = max(rows, key=lambda row: row["energy_increase"])
    _, worst_before = probe(binary, worst["sequence"], worst["scaffold"], False)
    _, worst_after = probe(binary, worst["sequence"], worst["scaffold"], True)
    original = run([str(args.oracle_checkout.resolve() / "build/CParty"), "-s", "0", "-d2", "-P",
                    str(ROOT / "params/rna_Turner04.par"), "-r", worst["scaffold"], worst["sequence"]], cwd=ROOT)
    if original.stdout != worst_before.stdout:
        raise ValueError("unguarded probe differs from original legacy on largest-effect case")
    write_log(out / "largest_effect_legacy.log", worst_before.stdout, worst_before.stderr)
    write_log(out / "largest_effect_completed_wmb_guard.log", worst_after.stdout, worst_after.stderr)
    summary = {"seed": args.seed, "random_cases": len(random_rows), "lengths": args.lengths,
               "cases_changed_above_1e-14_kcal": sum(row["energy_increase"] > 1e-14 for row in random_rows),
               "max_energy_increase_kcal": max((row["energy_increase"] for row in random_rows), default=0),
               "max_removed_fraction": max((row["removed_fraction"] for row in random_rows), default=0),
               "max_scfg2_legacy_abs_error_kcal": max(row["scfg2_legacy_abs_error"] for row in rows),
               "max_scfg2_guard_abs_error_kcal": max(row["scfg2_guard_abs_error"] for row in rows),
               "max_scfg2_delta_abs_error_kcal": max(abs(
                   (row["scfg2_guard_energy"] - row["scfg2_legacy_energy"]) - row["energy_increase"])
                   for row in rows),
               "cases_changed_above_1e-3_kcal": sum(row["energy_increase"] > 1e-3 for row in random_rows),
               "cases_negative_beyond_roundoff": sum(row["energy_increase"] < -1e-12 for row in rows),
               "closure_spans_verified": sum(row["closure_spans_verified"] for row in rows),
               "witness": rows[0], "witness_terms": witness_terms}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
