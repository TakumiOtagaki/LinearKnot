#!/usr/bin/env python3
"""Observe WIP PF/MFE cells and local MFE traceback in external legacy objects.

No oracle source is edited. MFE is a boundary regression observation, not the
semantic oracle; only legacy inside is used for numerical parity.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "test/oracle/legacy_cparty"))
from oracle_config import digest, verify_checkout

MAIN = r'''
#include "h_globals.hh"
#include "part_func.hh"
#include "W_final.hh"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

class MfeProbe : public W_final {
 public:
  using W_final::W_final;
  int wip(int i, int j) { return WMB->get_WIP(i,j); }
  std::string trace_wip(int i, int j, sparse_tree &tree) {
    if (wip(i,j) >= INF) return "";
    structure.assign(n + 1, '.');
    for (int pos = 1; pos <= n; ++pos) f[pos] = minimum_fold{};
    auto *current = new seq_interval{i, j, wip(i,j), P_WIP, nullptr};
    stack_interval = nullptr;
    while (current) {
      backtrack_restricted(current, tree);
      delete current;
      current = stack_interval;
      if (current) stack_interval = current->next;
    }
    return structure.substr(1,n);
  }
};

int main(int argc, char **argv) {
  if (argc != 7) return 2;
  std::string seq = argv[1], g = argv[2];
  const int dangles = std::atoi(argv[4]), i = std::atoi(argv[5]), j = std::atoi(argv[6]);
  if (seq.empty() || seq.size() != g.size() || i < 1 || j < i || j > int(seq.size())) return 2;
  if (!vrna_params_load(argv[3], VRNA_PARAMETER_FORMAT_DEFAULT)) return 2;
  sparse_tree tree(g, seq.size());
  SHAPEData shape("", seq.size());
  W_final_pf pf(seq, g, shape, false, false, dangles, 0.0, 0, false, false, 1.0);
  const double ensemble = pf.hfold_pf(tree);
  const double alpha = pf.get_energy_WIP(i,j);
  const double rt = pf.exp_params_->kT / 1000.0;
  MfeProbe mfe(seq, g, shape, false, false, dangles);
  mfe.hfold(tree);
  const auto trace = mfe.trace_wip(i,j,tree);
  std::cout << std::setprecision(17)
            << "{\"ensemble_energy\":" << ensemble << ",\"rt\":" << rt
            << ",\"pf_wip_weight\":" << alpha << ",\"pf_wip_energy\":";
  if (alpha > 0) std::cout << -rt * std::log(alpha); else std::cout << "null";
  std::cout << ",\"mfe_wip_dcal\":" << mfe.wip(i,j)
            << ",\"mfe_wip_trace\":\"" << trace << "\"}\n";
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sequence")
    parser.add_argument("scaffold")
    parser.add_argument("--cell", type=int, nargs=2, required=True, metavar=("I", "J"))
    parser.add_argument("--oracle-dir", type=Path, default=ROOT / "build/legacy_cparty_oracle/src")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--param-file", type=Path, default=ROOT / "params/rna_DirksPierce09.par")
    parser.add_argument("--dangles", type=int, choices=(0, 2), default=0)
    args = parser.parse_args()
    oracle, out, params = args.oracle_dir.resolve(), args.output_dir.resolve(), args.param_file.resolve()
    if out == oracle or oracle in out.parents:
        parser.error("output-dir must be separate from the oracle checkout")
    manifest = verify_checkout(oracle)
    out.mkdir(parents=True, exist_ok=True)
    source, binary = out / "main.cc", out / "wip-probe"
    source.write_text(MAIN)
    # Use the actual upstream objects, including its MFE traceback. The only
    # subclass above exposes protected observation calls; no recurrence copy.
    objects = sorted(p for p in (oracle / "build/CMakeFiles/CParty.dir/src").glob("*.cc.o")
                     if p.name != "CParty.cc.o")
    compiler = os.environ.get("CXX", "c++")
    subprocess.run([compiler, "-O2", "-std=gnu++17", "-I" + str(oracle / "src"),
                    str(source), *map(str, objects), str(oracle / "build/libRNA.a"),
                    "-o", str(binary)], check=True)
    result = json.loads(subprocess.check_output([
        str(binary), args.sequence, args.scaffold, str(params),
        str(args.dangles), *map(str, args.cell)], text=True))
    result.update(sequence=args.sequence, scaffold=args.scaffold, cell=args.cell,
                  temperature_c=37, dangles=args.dangles, pf_scale=1.0,
                  oracle_manifest=manifest, parameter_file=str(params),
                  parameter_sha256=digest(params), binary_sha256=digest(binary),
                  probe_source_sha256=digest(source),
                  object_sha256={str(p.relative_to(oracle)): digest(p) for p in objects},
                  compiler=subprocess.check_output([compiler, "--version"], text=True).splitlines()[0])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
