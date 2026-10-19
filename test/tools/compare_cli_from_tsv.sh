#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: compare_cli_from_tsv.sh <legacy_cparty_bin> <current_cparty_bin> <cases_tsv>

Compares CLI stdout (seq/restricted/MFE/energy) between legacy and current binaries
for each row in cases_tsv (header: case_id seq G).
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 3 ]]; then
  usage
  exit 2
fi

legacy_bin="$1"
current_bin="$2"
cases_tsv="$3"

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PARSER="$ROOT_DIR/test/tools/parse_cparty_stdout.sh"
PARAM_FILE="$ROOT_DIR/params/rna_DirksPierce09.par"

die() {
  echo "error: $*" >&2
  exit 1
}

require_exec() {
  local path="$1"
  [[ -x "$path" ]] || die "missing executable: $path"
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || die "missing file: $path"
}

require_exec "$legacy_bin"
require_exec "$current_bin"
require_file "$cases_tsv"
require_exec "$PARSER"
require_file "$PARAM_FILE"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

stdout_legacy="$tmpdir/legacy.txt"
stdout_current="$tmpdir/current.txt"
abs_diffs="$tmpdir/abs_diffs.txt"
rel_diffs="$tmpdir/rel_diffs.txt"
: > "$abs_diffs"
: > "$rel_diffs"

total=0
matched=0
skipped=0
failed=0
mfe_failed=0
pf_failed=0
mfe_energy_failed=0
pf_structure_failed=0
pf_energy_failed=0
legacy_pf_inf=0
current_pf_inf=0
legacy_pf_nan=0
current_pf_nan=0
legacy_pf_has_square=0
current_pf_has_square=0
legacy_pf_has_round=0
current_pf_has_round=0
legacy_pf_finite=0
current_pf_finite=0
sample_limit=5
sampled=0
skip_legacy_run=0
skip_current_run=0
skip_legacy_parse=0
skip_current_parse=0
legacy_backtrack_fail=0
current_backtrack_fail=0
backtrack_mismatch=0
legacy_parse_error=""
current_parse_error=""
pf_abs_tol="${PF_ABS_TOL:-1e-3}"
pf_rel_tol="${PF_REL_TOL:-1e-5}"
mfe_abs_tol="${MFE_ABS_TOL:-1e-6}"
mfe_rel_tol="${MFE_REL_TOL:-1e-6}"

is_finite_number() {
  local v="$1"
  case "$v" in
    inf|+inf|nan|+nan|-nan) return 1 ;;
  esac
  return 0
}

calc_abs_diff() {
  awk -v a="$1" -v b="$2" 'BEGIN{d=a-b; if(d<0)d=-d; print d}'
}

calc_rel_diff() {
  awk -v a="$1" -v b="$2" 'BEGIN{d=a-b; if(d<0)d=-d; denom=(a<0?-a:a); if(denom==0) denom=(b<0?-b:b); if(denom==0) denom=1; print d/denom}'
}

compare_float_with_tol() {
  local abs="$1"
  local rel="$2"
  local abs_tol="$3"
  local rel_tol="$4"
  awk -v a="$abs" -v r="$rel" -v at="$abs_tol" -v rt="$rel_tol" 'BEGIN{print (a<=at || r<=rt) ? 1 : 0}'
}

run_case() {
  local bin="$1"
  local g="$2"
  local seq="$3"
  local out="$4"
  local cli_args="${CPARTY_CLI_ARGS:--d2}"
  local -a cli_args_array
  # shellcheck disable=SC2206
  cli_args_array=($cli_args)
  "$bin" "${cli_args_array[@]}" -P "$PARAM_FILE" -r "$g" "$seq" > "$out" 2>/dev/null
}

has_backtracking_failed() {
  grep -qi "backtracking failed" "$1"
}

parse_stdout() {
  "$PARSER" "$1" 2>&1
}

field_at() {
  printf '%s' "$1" | awk -F $'\t' -v idx="$2" '{print $idx}'
}

record_pf_structure_stats() {
  local structure="$1"
  if [[ "$structure" == *"["* || "$structure" == *"]"* ]]; then
    legacy_pf_has_square=$((legacy_pf_has_square + 1))
  fi
  if [[ "$structure" == *"("* || "$structure" == *")"* ]]; then
    legacy_pf_has_round=$((legacy_pf_has_round + 1))
  fi
  if [[ "$structure" == "current" ]]; then
    : # placeholder to avoid accidental call
  fi
}

record_pf_energy_stats() {
  local energy="$1"
  if [[ "$energy" == "inf" || "$energy" == "+inf" ]]; then
    legacy_pf_inf=$((legacy_pf_inf + 1))
  elif [[ "$energy" == "nan" || "$energy" == "+nan" || "$energy" == "-nan" ]]; then
    legacy_pf_nan=$((legacy_pf_nan + 1))
  else
    legacy_pf_finite=$((legacy_pf_finite + 1))
  fi
}

record_pf_structure_stats_current() {
  local structure="$1"
  if [[ "$structure" == *"["* || "$structure" == *"]"* ]]; then
    current_pf_has_square=$((current_pf_has_square + 1))
  fi
  if [[ "$structure" == *"("* || "$structure" == *")"* ]]; then
    current_pf_has_round=$((current_pf_has_round + 1))
  fi
}

record_pf_energy_stats_current() {
  local energy="$1"
  if [[ "$energy" == "inf" || "$energy" == "+inf" ]]; then
    current_pf_inf=$((current_pf_inf + 1))
  elif [[ "$energy" == "nan" || "$energy" == "+nan" || "$energy" == "-nan" ]]; then
    current_pf_nan=$((current_pf_nan + 1))
  else
    current_pf_finite=$((current_pf_finite + 1))
  fi
}

compare_case() {
  local case_id="$1"
  local seq="$2"
  local g="$3"

  if ! run_case "$legacy_bin" "$g" "$seq" "$stdout_legacy"; then
    skipped=$((skipped + 1))
    skip_legacy_run=$((skip_legacy_run + 1))
    return 0
  fi
  if ! run_case "$current_bin" "$g" "$seq" "$stdout_current"; then
    skipped=$((skipped + 1))
    skip_current_run=$((skip_current_run + 1))
    return 0
  fi

  local legacy_bt=0
  local current_bt=0
  if has_backtracking_failed "$stdout_legacy"; then
    legacy_bt=1
    legacy_backtrack_fail=$((legacy_backtrack_fail + 1))
  fi
  if has_backtracking_failed "$stdout_current"; then
    current_bt=1
    current_backtrack_fail=$((current_backtrack_fail + 1))
  fi
  if (( legacy_bt == 1 || current_bt == 1 )); then
    if (( legacy_bt == 1 && current_bt == 1 )); then
      matched=$((matched + 1))
    else
      failed=$((failed + 1))
      backtrack_mismatch=$((backtrack_mismatch + 1))
      if (( sampled < sample_limit )); then
        echo "mismatch case_id=$case_id (backtracking)" >&2
        if (( legacy_bt == 1 )); then
          echo "legacy:  Backtracking failed" >&2
        else
          echo "legacy:  ok" >&2
        fi
        if (( current_bt == 1 )); then
          echo "current: Backtracking failed" >&2
        else
          echo "current: ok" >&2
        fi
        sampled=$((sampled + 1))
      fi
    fi
    return 0
  fi

  local legacy_parsed
  local current_parsed
  if ! legacy_parsed=$(parse_stdout "$stdout_legacy"); then
    skipped=$((skipped + 1))
    skip_legacy_parse=$((skip_legacy_parse + 1))
    if [[ -z "$legacy_parse_error" ]]; then
      legacy_parse_error="$legacy_parsed"
    fi
    return 0
  fi
  if ! current_parsed=$(parse_stdout "$stdout_current"); then
    skipped=$((skipped + 1))
    skip_current_parse=$((skip_current_parse + 1))
    if [[ -z "$current_parse_error" ]]; then
      current_parse_error="$current_parsed"
    fi
    return 0
  fi

  local legacy_seq legacy_restricted legacy_mfe_structure legacy_mfe_energy legacy_pf_structure legacy_pf_energy
  local current_seq current_restricted current_mfe_structure current_mfe_energy current_pf_structure current_pf_energy
  legacy_seq=$(field_at "$legacy_parsed" 1)
  legacy_restricted=$(field_at "$legacy_parsed" 2)
  legacy_mfe_structure=$(field_at "$legacy_parsed" 3)
  legacy_mfe_energy=$(field_at "$legacy_parsed" 4)
  legacy_pf_structure=$(field_at "$legacy_parsed" 5)
  legacy_pf_energy=$(field_at "$legacy_parsed" 6)

  current_seq=$(field_at "$current_parsed" 1)
  current_restricted=$(field_at "$current_parsed" 2)
  current_mfe_structure=$(field_at "$current_parsed" 3)
  current_mfe_energy=$(field_at "$current_parsed" 4)
  current_pf_structure=$(field_at "$current_parsed" 5)
  current_pf_energy=$(field_at "$current_parsed" 6)

  local mfe_ok=1
  local pf_ok=1
  local pf_energy_ok=1
  local mfe_energy_ok=1

  if [[ "$legacy_seq" != "$current_seq" || "$legacy_restricted" != "$current_restricted" ]]; then
    mfe_ok=0
    pf_ok=0
  else
    if [[ "$legacy_mfe_structure" != "$current_mfe_structure" ]]; then
      mfe_ok=0
    else
      mfe_energy_ok=0
      if is_finite_number "$legacy_mfe_energy" && is_finite_number "$current_mfe_energy"; then
        local mfe_abs_diff mfe_rel_diff
        mfe_abs_diff=$(calc_abs_diff "$legacy_mfe_energy" "$current_mfe_energy")
        mfe_rel_diff=$(calc_rel_diff "$legacy_mfe_energy" "$current_mfe_energy")
        mfe_energy_ok=$(compare_float_with_tol "$mfe_abs_diff" "$mfe_rel_diff" "$mfe_abs_tol" "$mfe_rel_tol")
      else
        if [[ "$legacy_mfe_energy" == "$current_mfe_energy" ]]; then
          mfe_energy_ok=1
        else
          mfe_energy_ok=0
        fi
      fi
      if (( mfe_energy_ok == 0 )); then
        mfe_ok=0
      fi
    fi
    if [[ "$legacy_pf_structure" != "$current_pf_structure" ]]; then
      pf_ok=0
    else
      pf_energy_ok=0
      if is_finite_number "$legacy_pf_energy" && is_finite_number "$current_pf_energy"; then
        local abs_diff rel_diff
        abs_diff=$(calc_abs_diff "$legacy_pf_energy" "$current_pf_energy")
        rel_diff=$(calc_rel_diff "$legacy_pf_energy" "$current_pf_energy")
        echo "$abs_diff" >> "$abs_diffs"
        echo "$rel_diff" >> "$rel_diffs"
        pf_energy_ok=$(compare_float_with_tol "$abs_diff" "$rel_diff" "$pf_abs_tol" "$pf_rel_tol")
      else
        pf_energy_ok=0
      fi
      if (( pf_energy_ok == 0 )); then
        pf_ok=0
      fi
    fi
  fi

  record_pf_structure_stats "$legacy_pf_structure"
  record_pf_structure_stats_current "$current_pf_structure"
  record_pf_energy_stats "$legacy_pf_energy"
  record_pf_energy_stats_current "$current_pf_energy"

  if (( mfe_ok == 1 && pf_ok == 1 )); then
    matched=$((matched + 1))
  else
    failed=$((failed + 1))
    if (( mfe_ok == 0 )); then
      mfe_failed=$((mfe_failed + 1))
      if (( mfe_energy_ok == 0 )); then
        mfe_energy_failed=$((mfe_energy_failed + 1))
      fi
    fi
    if (( pf_ok == 0 )); then
      pf_failed=$((pf_failed + 1))
      if [[ "$legacy_pf_structure" != "$current_pf_structure" ]]; then
        pf_structure_failed=$((pf_structure_failed + 1))
      fi
      if (( pf_energy_ok == 0 )); then
        pf_energy_failed=$((pf_energy_failed + 1))
      fi
    fi
    if (( sampled < sample_limit )); then
      echo "mismatch case_id=$case_id" >&2
      echo "legacy:  $legacy_parsed" >&2
      echo "current: $current_parsed" >&2
      sampled=$((sampled + 1))
    fi
  fi
}

while IFS=$'\t' read -r case_id seq g; do
  if [[ "$case_id" == "case_id" ]]; then
    continue
  fi
  total=$((total + 1))
  compare_case "$case_id" "$seq" "$g"
done < "$cases_tsv"

echo "compare_total=$total"
echo "compare_matched=$matched"
echo "compare_failed=$failed"
echo "compare_skipped=$skipped"
echo "compare_mfe_failed=$mfe_failed"
echo "compare_mfe_energy_failed=$mfe_energy_failed"
echo "compare_pf_failed=$pf_failed"
echo "compare_pf_structure_failed=$pf_structure_failed"
echo "compare_pf_energy_failed=$pf_energy_failed"
if [[ -s "$abs_diffs" ]]; then
  abs_top3=$(sort -nr "$abs_diffs" | head -n 3 | paste -sd, -)
  abs_median=$(sort -n "$abs_diffs" | awk '{a[NR]=$1} END{if(NR==0){print ""} else {m=int((NR+1)/2); print a[m]}}')
  rel_top3=$(sort -nr "$rel_diffs" | head -n 3 | paste -sd, -)
  rel_median=$(sort -n "$rel_diffs" | awk '{a[NR]=$1} END{if(NR==0){print ""} else {m=int((NR+1)/2); print a[m]}}')
  echo "pf_abs_top3_max=${abs_top3}"
  echo "pf_rel_top3_max=${rel_top3}"
  echo "pf_abs_median=${abs_median}"
  echo "pf_rel_median=${rel_median}"
  echo "pf_abs_tol=${pf_abs_tol}"
  echo "pf_rel_tol=${pf_rel_tol}"
fi
echo "skip_legacy_run=$skip_legacy_run"
echo "skip_current_run=$skip_current_run"
echo "skip_legacy_parse=$skip_legacy_parse"
echo "skip_current_parse=$skip_current_parse"
if [[ -n "$legacy_parse_error" ]]; then
  echo "legacy_parse_error_sample=$legacy_parse_error"
fi
if [[ -n "$current_parse_error" ]]; then
  echo "current_parse_error_sample=$current_parse_error"
fi
echo "legacy_pf_inf=$legacy_pf_inf"
echo "legacy_pf_nan=$legacy_pf_nan"
echo "legacy_pf_finite=$legacy_pf_finite"
echo "current_pf_inf=$current_pf_inf"
echo "current_pf_nan=$current_pf_nan"
echo "current_pf_finite=$current_pf_finite"
echo "legacy_pf_has_square=$legacy_pf_has_square"
echo "legacy_pf_has_round=$legacy_pf_has_round"
echo "current_pf_has_square=$current_pf_has_square"
echo "current_pf_has_round=$current_pf_has_round"
echo "legacy_backtracking_failed=$legacy_backtrack_fail"
echo "current_backtracking_failed=$current_backtrack_fail"
echo "backtracking_mismatch=$backtrack_mismatch"
compare_min="${COMPARE_MIN_MATCHED:-100}"
compared=$((matched + failed))
echo "compare_compared=$compared"
echo "compare_min=$compare_min"

if (( total < 1 )); then
  echo "compare gate failed: no cases executed" >&2
  exit 1
fi
if (( compared < compare_min )); then
  echo "compare gate failed: compared=$compared (<$compare_min)" >&2
  exit 1
fi
if (( failed != 0 )); then
  echo "compare gate failed: mismatched=$failed" >&2
  exit 1
fi

exit 0
