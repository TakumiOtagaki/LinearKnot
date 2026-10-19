#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: pf_exact_legacy_alignment_random_constrained_scfg2.sh <legacy_cparty_bin> <scfg2_pf_tool> <param_file> <seed> <count> <lengths>

Generates random sequence / pk-free-constraint pairs and compares constrained PF
energies between:
  1. legacy CParty `-d2 -r G` scaffold PF route
  2. `scfg2_restricted_partition_function_tool --legacy-compatible`

This is the direct legacy-vs-scfg2 oracle helper for constrained exact PF.

lengths: comma-separated list of lengths (e.g. 30,50,80)

Env:
  PF_ABS_TOL=1e-3
  PF_REL_TOL=1e-5
  ALIGNMENT_MIN=10
  MAX_SKIP_PERCENT=5
  DIAG_CSV=/tmp/pf_energy_exact_legacy_scfg2_random_diag.csv
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 6 ]]; then
  usage
  exit 2
fi

legacy_bin="$1"
scfg2_tool="$2"
param_file="$3"
seed="$4"
count="$5"
lengths_csv="$6"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GENERATOR="$ROOT_DIR/test/tools/gen_random_seq_struct.py"
PARSER="$ROOT_DIR/test/tools/parse_cparty_stdout.sh"
STATS="$ROOT_DIR/test/tools/subtools/pf_alignment_stats.py"

PF_ABS_TOL="${PF_ABS_TOL:-1e-3}"
PF_REL_TOL="${PF_REL_TOL:-1e-5}"
ALIGNMENT_MIN="${ALIGNMENT_MIN:-10}"
MAX_SKIP_PERCENT="${MAX_SKIP_PERCENT:-5}"
DIAG_CSV="${DIAG_CSV:-}"

require_exec() { [[ -x "$1" ]] || { echo "error: missing executable: $1" >&2; exit 2; }; }
require_file() { [[ -f "$1" ]] || { echo "error: missing file: $1" >&2; exit 2; }; }

require_exec "$legacy_bin"
require_exec "$scfg2_tool"
require_file "$param_file"
require_exec "$GENERATOR"
require_exec "$PARSER"
require_file "$STATS"

IFS=',' read -r -a lengths <<< "$lengths_csv"

if [[ -n "$DIAG_CSV" ]]; then
  printf 'case_id,length,seq,G,pf_legacy,pf_scfg2,abs_diff,rel_diff,matched\n' > "$DIAG_CSV"
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

abs_diffs="$tmpdir/abs_diffs.txt"
rel_diffs="$tmpdir/rel_diffs.txt"
: > "$abs_diffs"
: > "$rel_diffs"

total=0
compared=0
matched=0
failed=0
skipped=0
legacy_exec_skipped=0
scfg2_exec_skipped=0
legacy_parse_skipped=0
identity_skipped=0
nonfinite_skipped=0

calc_abs() { awk -v a="$1" -v b="$2" 'BEGIN{d=a-b; if(d<0)d=-d; print d}'; }
calc_rel() { awk -v a="$1" -v b="$2" 'BEGIN{d=a-b; if(d<0)d=-d; denom=(a<0?-a:a); if(denom==0) denom=(b<0?-b:b); if(denom==0) denom=1; print d/denom}'; }

for length in "${lengths[@]}"; do
  cases_file="$tmpdir/cases_${length}.tsv"
  if ! python3 "$GENERATOR" --seed "$seed" --count "$count" --length "$length" --prefix "rand${length}" > "$cases_file"; then
    echo "error: failed to generate cases for length $length" >&2
    exit 1
  fi

  while IFS=$'\t' read -r case_id seq g; do
    if [[ "$case_id" == "case_id" ]]; then
      continue
    fi
    total=$((total + 1))
    if [[ -z "$seq" || -z "$g" ]]; then
      skipped=$((skipped + 1))
      identity_skipped=$((identity_skipped + 1))
      continue
    fi
    if [[ ${#seq} -ne ${#g} ]]; then
      skipped=$((skipped + 1))
      identity_skipped=$((identity_skipped + 1))
      continue
    fi

    legacy_out="$tmpdir/legacy_${case_id}.txt"
    scfg2_out="$tmpdir/scfg2_${case_id}.txt"
    if ! env -u SCFG_DENSITY2_VIEW "$legacy_bin" -d2 -P "$param_file" -r "$g" "$seq" > "$legacy_out" 2>/dev/null; then
      skipped=$((skipped + 1))
      legacy_exec_skipped=$((legacy_exec_skipped + 1))
      continue
    fi
    if ! "$scfg2_tool" --legacy-compatible "$seq" "$g" "$param_file" > "$scfg2_out" 2>/dev/null; then
      skipped=$((skipped + 1))
      scfg2_exec_skipped=$((scfg2_exec_skipped + 1))
      continue
    fi

    if ! legacy_parsed=$("$PARSER" "$legacy_out" 2>/dev/null); then
      skipped=$((skipped + 1))
      legacy_parse_skipped=$((legacy_parse_skipped + 1))
      continue
    fi
    IFS=$'\t' read -r parsed_seq parsed_g _ _ _ legacy_pf <<< "$legacy_parsed"
    scfg2_pf=$(tr -d '[:space:]' < "$scfg2_out")

    if [[ "$parsed_seq" != "$seq" || "$parsed_g" != "$g" ]]; then
      skipped=$((skipped + 1))
      identity_skipped=$((identity_skipped + 1))
      continue
    fi

    finite_number_re='^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
    if [[ ! "$legacy_pf" =~ $finite_number_re || ! "$scfg2_pf" =~ $finite_number_re ]]; then
      skipped=$((skipped + 1))
      nonfinite_skipped=$((nonfinite_skipped + 1))
      continue
    fi

    compared=$((compared + 1))
    abs=$(calc_abs "$legacy_pf" "$scfg2_pf")
    rel=$(calc_rel "$legacy_pf" "$scfg2_pf")
    printf '%s\n' "$abs" >> "$abs_diffs"
    printf '%s\n' "$rel" >> "$rel_diffs"

    ok=$(awk -v a="$abs" -v r="$rel" -v at="$PF_ABS_TOL" -v rt="$PF_REL_TOL" 'BEGIN{print (a<=at || r<=rt)?1:0}')
    if [[ "$ok" == "1" ]]; then
      matched=$((matched + 1))
      status=1
    else
      failed=$((failed + 1))
      status=0
    fi

    if [[ -n "$DIAG_CSV" ]]; then
      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$case_id" "$length" "$seq" "$g" "$legacy_pf" "$scfg2_pf" "$abs" "$rel" "$status" >> "$DIAG_CSV"
    fi
  done < "$cases_file"
done

printf 'skip_reasons legacy_exec=%d scfg2_exec=%d legacy_parse=%d identity=%d nonfinite=%d\n' \
  "$legacy_exec_skipped" "$scfg2_exec_skipped" "$legacy_parse_skipped" \
  "$identity_skipped" "$nonfinite_skipped"

if (( total > 0 && skipped * 100 > total * MAX_SKIP_PERCENT )); then
  echo "error: skipped $skipped/$total cases exceeds MAX_SKIP_PERCENT=$MAX_SKIP_PERCENT" >&2
  exit 1
fi

python3 "$STATS" \
  "$abs_diffs" "$rel_diffs" "$compared" "$ALIGNMENT_MIN" "$PF_ABS_TOL" "$PF_REL_TOL" \
  "$matched" "$failed" "$skipped" "$total"
