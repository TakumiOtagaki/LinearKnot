#!/usr/bin/env python3
"""Mark a V subtree in an isolated copy of patched legacy inside.

Uses the pinned oracle's headers and objects, replacing only part_func.cc.
Neither the primary oracle checkout nor the SCFG2 kernel is modified.
MFE=0 fixes pf_scale=1; this is a direct inside probe, not the public CLI.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys


MAIN = r'''
#include "h_globals.hh"
#include "part_func.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
  if (argc != 3) return 2;
  std::ifstream input(argv[1]);
  std::string seq, scaffold, line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '>') continue;
    if (seq.empty()) seq = line; else { scaffold = line; break; }
  }
  if (seq.empty() || seq.size() != scaffold.size()) return 2;
  if (!vrna_params_load(argv[2], VRNA_PARAMETER_FORMAT_DEFAULT)) return 2;
  sparse_tree tree(scaffold, seq.size());
  SHAPEData shape("", seq.size());
  W_final_pf inside(seq, scaffold, shape, false, false, 2, 0.0,
                    0, false, false, 1.0);
  const auto energy = inside.hfold_pf(tree);
  std::cout << std::setprecision(17) << energy << ' '
            << inside.exp_params_->kT / 1000.0 << '\n';
}
'''

MARK_V = r'''
    // External probe: each consumption of this subtree contributes one x.
    const char *mark_i = std::getenv("CPARTY_PROBE_V_I");
    const char *mark_j = std::getenv("CPARTY_PROBE_V_J");
    const char *mark_x = std::getenv("CPARTY_PROBE_V_SCALE");
    if (mark_i && mark_j && mark_x &&
        i == std::atoi(mark_i) && j == std::atoi(mark_j))
        contributions *= std::atof(mark_x);
    V[ij] = contributions;'''

TRACE_READ = r'''
                        if (std::getenv("CPARTY_PROBE_TRACE_I") &&
                            i == std::atoi(std::getenv("CPARTY_PROBE_TRACE_I")) &&
                            j == std::atoi(std::getenv("CPARTY_PROBE_TRACE_J")) &&
                            l == std::atoi(std::getenv("CPARTY_PROBE_TRACE_L")))
                          std::cerr << "WMBP(" << i << ',' << j << ") split=" << l
                                    << " BE(" << tree.tree[B_lj].pair << ',' << B_lj
                                    << ';' << tree.tree[Bp_lj].pair << ',' << Bp_lj
                                    << ") read="
                                    << get_BE(tree.tree[B_lj].pair, B_lj,
                                              tree.tree[Bp_lj].pair, Bp_lj, tree)
                                    << " m1=" << m1 << '\n';
                        contributions += m1;'''

TRACE_WRITE = r'''
    if (std::getenv("CPARTY_PROBE_TRACE_I")) {
      const int trace_j = std::atoi(std::getenv("CPARTY_PROBE_TRACE_J"));
      const int trace_l = std::atoi(std::getenv("CPARTY_PROBE_TRACE_L"));
      const int outer = tree.B(trace_l, trace_j);
      const int inner = tree.Bp(trace_l, trace_j);
      if (outer > 0 && inner > 0 && i == tree.tree[outer].pair &&
          j == outer && ip == tree.tree[inner].pair && jp == inner)
        std::cerr << "later compute_BE(" << i << ',' << j << ';'
                  << ip << ',' << jp << ")=" << contributions << '\n';
    }
    BE[iip] = contributions;'''


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise RuntimeError(f"legacy source no longer has a unique anchor: {old!r}")
    return source.replace(old, new)


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("i", type=int)
    parser.add_argument("j", type=int)
    parser.add_argument("--oracle-dir", type=Path,
                        default=root / "build/legacy_cparty_oracle/src")
    parser.add_argument("--output-dir", type=Path,
                        default=root / "build/legacy_cparty_pair_occurrence_probe")
    parser.add_argument("--param-file", type=Path,
                        default=root / "params/rna_Turner04.par")
    parser.add_argument("--trace-wmbp", nargs=3, type=int, metavar=("I", "J", "L"))
    args = parser.parse_args()
    oracle = args.oracle_dir.resolve()
    output = args.output_dir.resolve()
    if output == oracle or oracle in output.parents:
        parser.error("output-dir must be separate from the primary oracle")
    sys.path.insert(0, str(root / "test/oracle/legacy_cparty"))
    from oracle_config import verify_checkout
    manifest = verify_checkout(oracle)
    head = manifest["upstream_commit"]
    source = (oracle / "src/part_func.cc").read_text()
    source_sha = hashlib.sha256(source.encode()).hexdigest()
    source = replace_once(source, "    V[ij] = contributions;", MARK_V)
    source = replace_once(source, "                        contributions += m1;", TRACE_READ)
    source = replace_once(source, "    BE[iip] = contributions;", TRACE_WRITE)
    source = "#include <cstdlib>\n" + source
    output.mkdir(parents=True, exist_ok=True)
    (output / "part_func.cc").write_text(source)
    (output / "main.cc").write_text(MAIN)
    cxx = os.environ.get("CXX", "c++")
    common = [cxx, "-O2", "-std=gnu++17", "-I" + str(oracle / "src")]
    subprocess.run(common + ["-c", str(output / "part_func.cc"),
                            "-o", str(output / "part_func.o")], check=True)
    objects = oracle / "build/CMakeFiles/CParty.dir/src"
    deps = [str(objects / (name + ".cc.o"))
            for name in ["sparse_tree", "SHAPE", "dot_plot", "mea", "centroid"]]
    executable = output / "inside-probe"
    subprocess.run(common + [str(output / "main.cc"), str(output / "part_func.o")]
                   + deps + [str(oracle / "build/libRNA.a"), "-o", str(executable)],
                   check=True)
    energies = {}
    trace = ""
    for x in [0.0, 0.5, 1.0, 2.0, 4.0]:
        env = {k: v for k, v in os.environ.items()
               if not k.startswith("CPARTY_PROBE_")}
        env.update(CPARTY_PROBE_V_I=str(args.i), CPARTY_PROBE_V_J=str(args.j),
                   CPARTY_PROBE_V_SCALE=str(x))
        if x == 1.0 and args.trace_wmbp:
            env.update(zip(["CPARTY_PROBE_TRACE_I", "CPARTY_PROBE_TRACE_J",
                            "CPARTY_PROBE_TRACE_L"], map(str, args.trace_wmbp)))
        result = subprocess.run([str(executable), str(args.input.resolve()),
                                 str(args.param_file.resolve())], env=env,
                                text=True, capture_output=True, check=True)
        energy, rt = map(float, result.stdout.split())
        energies[x] = energy
        if x == 1.0:
            trace = result.stderr
    q = {x: math.exp(-(energy - energies[1.0]) / rt)
         for x, energy in energies.items()}
    p0 = q[0.0]
    p2 = (q[2.0] - 2.0 + p0) / 2.0
    p1 = 1.0 - p0 - p2
    residuals = {x: q[x] - (p0 + x * p1 + x * x * p2) for x in [0.5, 4.0]}
    print(json.dumps({"oracle_commit": head, "part_func_sha256": source_sha,
                      "oracle_manifest_sha256": hashlib.sha256((oracle.parent / "manifest.json").read_bytes()).hexdigest(),
                      "input": str(args.input), "target": [args.i, args.j],
                      "mfe_seed": 0.0, "pf_scale": 1.0,
                      "energy": energies[1.0], "RT": rt, "normalized_Z": q,
                      "quadratic_fit": {"p_zero": p0, "p_once": p1, "p_twice": p2,
                                        "expected_count": p1 + 2 * p2,
                                        "marker_residuals": residuals,
                                        "fit_ok": max(map(abs, residuals.values())) < 1e-8},
                      "availability_trace": trace}, indent=2))


if __name__ == "__main__":
    main()
