#!/usr/bin/env python3
"""Observe patched legacy inside and VP split ranges without changing its code."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


MAIN = r'''
#include "h_globals.hh"
#include "part_func.hh"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
  if (argc != 4 && argc != 6) return 2;
  std::string seq = argv[1], g = argv[2];
  if (seq.empty() || seq.size() != g.size()) return 2;
  if (!vrna_params_load(argv[3], VRNA_PARAMETER_FORMAT_DEFAULT)) return 2;
  sparse_tree tree(g, seq.size());
  SHAPEData shape("", seq.size());
  W_final_pf inside(seq, g, shape, false, false, 2, 0.0, 0, false, false, 1.0);
  const auto energy = inside.hfold_pf(tree);
  std::size_t bp_zero_count = 0;
  for (int l = 1; l <= static_cast<int>(seq.size()); ++l)
    for (int j = l; j <= static_cast<int>(seq.size()); ++j)
      bp_zero_count += tree.Bp(l,j) == 0;
  std::cout << std::setprecision(17) << "{\"energy\":" << energy
            << ",\"rt\":" << inside.exp_params_->kT / 1000.0
            << ",\"log_z\":" << -energy / (inside.exp_params_->kT / 1000.0)
            << ",\"bp_zero_count\":" << bp_zero_count;
  if (argc == 6) {
    const cand_pos_t i = std::atoi(argv[4]), j = std::atoi(argv[5]);
    if (i < 1 || i >= j || j > static_cast<cand_pos_t>(seq.size())) return 2;
    const cand_pos_t end_left = std::min(static_cast<cand_pos_tu>(tree.b(i,j)),
                                        static_cast<cand_pos_tu>(tree.Bp(i,j)));
    const cand_pos_t begin_right = std::max(tree.B(i,j), tree.bp(i,j)) + 1;
    std::cout << ",\"parent\":[" << i << ',' << j << "]"
              << ",\"vp\":" << inside.get_energy_VP(i,j)
              << ",\"left_range\":[" << i+1 << ',' << end_left << "]"
              << ",\"right_range\":[" << begin_right << ',' << j << "]"
              << ",\"positive_child_products\":[";
    bool first = true;
    const auto emit = [&](const char *rule, cand_pos_t k, double a, double b) {
      if (a == 0.0 || b == 0.0) return;
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"rule\":\"" << rule << "\",\"k\":" << k
                << ",\"left\":" << a << ",\"right\":" << b << '}';
    };
    for (cand_pos_t k = i+1; k < end_left; ++k) {
      emit("VP_WIP_VP_LEFT", k, inside.get_energy_WIP(i+1,k-1), inside.get_energy_VP(k,j-1));
      emit("VP_WIP_VPR", k, inside.get_energy_WIP(i+1,k-1), inside.get_energy_VPR(k,j-1));
    }
    for (cand_pos_t k = begin_right; k < j; ++k) {
      emit("VP_VP_WIP_RIGHT", k, inside.get_energy_VP(i+1,k), inside.get_energy_WIP(k+1,j-1));
      emit("VP_VPL_WIP", k, inside.get_energy_VPL(i+1,k), inside.get_energy_WIP(k+1,j-1));
    }
    std::cout << ']';
  }
  std::cout << "}\n";
}
'''


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sequence")
    parser.add_argument("scaffold")
    parser.add_argument("--parent", nargs=2, type=int, metavar=("I", "J"))
    parser.add_argument("--oracle-dir", type=Path,
                        default=root / "build/legacy_cparty_oracle/src")
    parser.add_argument("--output-dir", type=Path,
                        default=root / "build/legacy_cparty_vp_split_probe")
    parser.add_argument("--param-file", type=Path,
                        default=root / "params/rna_Turner04.par")
    args = parser.parse_args()
    if not args.sequence or len(args.sequence) != len(args.scaffold):
        parser.error("sequence and scaffold must have equal nonzero lengths")
    oracle, output = args.oracle_dir.resolve(), args.output_dir.resolve()
    if output == oracle or oracle in output.parents:
        parser.error("output-dir must be separate from the primary oracle")
    sys.path.insert(0, str(root / "test/oracle/legacy_cparty"))
    from oracle_config import verify_checkout
    manifest = verify_checkout(oracle)
    head = manifest["upstream_commit"]
    output.mkdir(parents=True, exist_ok=True)
    source = output / "main.cc"
    source.write_text(MAIN)
    objects = oracle / "build/CMakeFiles/CParty.dir/src"
    deps = [objects / (name + ".cc.o") for name in
            ("part_func", "sparse_tree", "SHAPE", "dot_plot", "mea", "centroid")]
    executable = output / "inside-probe"
    subprocess.run([os.environ.get("CXX", "c++"), "-O2", "-std=gnu++17",
                    "-I" + str(oracle / "src"), str(source), *map(str, deps),
                    str(oracle / "build/libRNA.a"), "-o", str(executable)], check=True)
    command = [str(executable), args.sequence, args.scaffold, str(args.param_file.resolve())]
    if args.parent:
        command += list(map(str, args.parent))
    result = json.loads(subprocess.check_output(command, text=True))
    result["oracle_commit"] = head
    result["oracle_manifest_sha256"] = hashlib.sha256((oracle.parent / "manifest.json").read_bytes()).hexdigest()
    result["sequence"] = args.sequence
    result["scaffold"] = args.scaffold
    result["part_func_sha256"] = hashlib.sha256((oracle / "src/part_func.cc").read_bytes()).hexdigest()
    result["mfe_seed"] = 0.0
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
