#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: compare_pf_energy_random.sh <legacy_cparty_bin> <current_cparty_bin> <param_file> <seed> <count> <lengths>

Generates random sequence / pk-free-constraint pairs and compares PF energies
between legacy and current binaries under -d2 with `-r G`.

This helper is oracle-facing on the inside-first refactor branch, so it clears
any inherited `SCFG_DENSITY2_VIEW` setting before invoking either binary.

lengths: comma-separated list of lengths (e.g. 30,50,80)

Env:
  PF_ABS_TOL=1e-3
  PF_REL_TOL=1e-5
  ALIGNMENT_MIN=10
  DIAG_CSV=/tmp/pf_energy_diag.csv
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
current_bin="$2"
param_file="$3"
seed="$4"
count="$5"
lengths_csv="$6"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GENERATOR="$ROOT_DIR/test/tools/gen_random_seq_struct.py"
PARSER="$ROOT_DIR/test/tools/parse_cparty_stdout.sh"

PF_ABS_TOL="${PF_ABS_TOL:-1e-3}"
PF_REL_TOL="${PF_REL_TOL:-1e-5}"
ALIGNMENT_MIN="${ALIGNMENT_MIN:-10}"
DIAG_CSV="${DIAG_CSV:-}"

require_exec() { [[ -x "$1" ]] || { echo "error: missing executable: $1" >&2; exit 2; }; }
require_file() { [[ -f "$1" ]] || { echo "error: missing file: $1" >&2; exit 2; }; }

require_exec "$legacy_bin"
require_exec "$current_bin"
require_file "$param_file"
require_exec "$GENERATOR"
require_exec "$PARSER"

IFS=',' read -r -a lengths <<< "$lengths_csv"

if [[ -n "$DIAG_CSV" ]]; then
  printf 'case_id,length,seq,G,pf_legacy,pf_current,abs_diff,rel_diff,matched\n' > "$DIAG_CSV"
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
      continue
    fi

    legacy_out="$tmpdir/legacy_${case_id}.txt"
    current_out="$tmpdir/current_${case_id}.txt"
    if ! env -u SCFG_DENSITY2_VIEW "$legacy_bin" -d2 -P "$param_file" -r "$g" "$seq" > "$legacy_out" 2>/dev/null; then
      skipped=$((skipped + 1))
      continue
    fi
    if ! env -u SCFG_DENSITY2_VIEW "$current_bin" -d2 -P "$param_file" -r "$g" "$seq" > "$current_out" 2>/dev/null; then
      skipped=$((skipped + 1))
      continue
    fi

    if ! legacy_parsed=$("$PARSER" "$legacy_out" 2>/dev/null); then
      skipped=$((skipped + 1))
      continue
    fi
    if ! current_parsed=$("$PARSER" "$current_out" 2>/dev/null); then
      skipped=$((skipped + 1))
      continue
    fi

    legacy_pf=$(printf '%s' "$legacy_parsed" | awk -F $'\t' '{print $6}')
    current_pf=$(printf '%s' "$current_parsed" | awk -F $'\t' '{print $6}')

    if [[ -z "$legacy_pf" || -z "$current_pf" ]]; then
      skipped=$((skipped + 1))
      continue
    fi

    compared=$((compared + 1))
    abs=$(calc_abs "$legacy_pf" "$current_pf")
    rel=$(calc_rel "$legacy_pf" "$current_pf")
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
        "$case_id" "$length" "$seq" "$g" "$legacy_pf" "$current_pf" "$abs" "$rel" "$status" >> "$DIAG_CSV"
    fi
  done < "$cases_file"

done

python3 "$ROOT_DIR/test/tools/subtools/pf_alignment_stats.py" \
  "$abs_diffs" "$rel_diffs" "$compared" "$ALIGNMENT_MIN" "$PF_ABS_TOL" "$PF_REL_TOL" \
  "$matched" "$failed" "$skipped" "$total"
