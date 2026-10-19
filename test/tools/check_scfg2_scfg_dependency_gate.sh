#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"

cd "$repo_root"
if [[ -e src/scfg ]]; then
  echo "scfg2/scfg dependency gate failed: retired src/scfg tree exists" >&2
  exit 1
fi

for retired_path in \
    src/CParty.cc \
    src/W_final.cc \
    src/W_final.hh \
    src/sparse_tree.cc \
    src/sparse_tree.hh \
    src/part_func.hh \
    src/part_func_init.cc \
    src/part_func_energy.cc \
    src/part_func_rescale.cc \
    src/h_globals.cc \
    src/h_globals.hh \
    src/h_externs.hh \
    src/pf_globals.hh \
    src/pf_penalty_globals.hh \
    src/base_types.hh \
    src/part_func_driver.cc \
    src/part_func_driver \
    src/pseudo_loop.cc \
    src/pseudo_loop.hh \
    src/pseudo_loop_backtrack; do
  if [[ -e "$retired_path" ]]; then
    echo "scfg2/scfg dependency gate failed: retired runtime path exists: $retired_path" >&2
    exit 1
  fi
done

if rg -n '#include "scfg/' CMakeLists.txt src; then
  echo "scfg2/scfg dependency gate failed: active runtime references retired scfg" >&2
  exit 1
fi

if rg -n '#include ["<](part_func|sparse_tree|h_globals|h_externs|pf_globals|pf_penalty_globals|base_types)\.hh|\bW_final_pf\b' src/scfg2 -g '*.hh' -g '*.cc'; then
  echo "scfg2/scfg dependency gate failed: active kernel references retired owner support" >&2
  exit 1
fi

echo "scfg2_scfg_dependency_gate=ok"
