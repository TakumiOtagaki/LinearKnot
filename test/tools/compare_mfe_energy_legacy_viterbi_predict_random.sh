#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: compare_mfe_energy_legacy_viterbi_predict_random.sh <legacy_cparty_bin> <viterbi_predict_tool> <param_file> <seed> <count> <lengths>

Random expansion + alignment report between legacy MFE (with -r constraints) and
SCFG viterbi prediction (energy + structure). Input G is a random pk-free structure.
lengths: comma-separated list of lengths (e.g. 30,50,80)
Environment:
  ALIGNMENT_PK_FREE_ONLY=1  Compare pk-free structures only.
  ALIGNMENT_PK_ONLY=1       Compare pk ([]) structures only.
USAGE
}

topology_family_of() {
  local st="$1"
  if [[ "$st" == *"("* || "$st" == *")"* ]]; then
    if [[ "$st" == *"["* || "$st" == *"]"* ]]; then
      printf 'density2_pk'
      return
    fi
    printf 'pk_free'
    return
  fi
  if [[ "$st" == *"["* || "$st" == *"]"* ]]; then
    printf 'other_pk'
    return
  fi
  printf 'none'
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
tool_bin="$2"
param_file="$3"
seed="$4"
count="$5"
lengths_csv="$6"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GENERATOR="$ROOT_DIR/test/tools/gen_random_seq_struct.py"
MATH_TOOL="$ROOT_DIR/test/tools/subtools/mfe_alignment_math.py"
DIAG_TOOL="$ROOT_DIR/test/tools/subtools/mfe_alignment_diag.py"
EXPORTER="$ROOT_DIR/test/tools/export_legacy_replay_artifact.py"

abs_tol=1e-7
rel_tol=1e-7
alignment_min="${ALIGNMENT_MIN:-10}"
diag_csv="${ALIGNMENT_DIAG_CSV:-}"
alignment_abs_median_max="${ALIGNMENT_ABS_MEDIAN_MAX:-}"
alignment_abs_p90_max="${ALIGNMENT_ABS_P90_MAX:-}"
alignment_abs_p95_max="${ALIGNMENT_ABS_P95_MAX:-}"
alignment_rel_median_max="${ALIGNMENT_REL_MEDIAN_MAX:-}"
alignment_rel_p90_max="${ALIGNMENT_REL_P90_MAX:-}"
alignment_rel_p95_max="${ALIGNMENT_REL_P95_MAX:-}"
alignment_rel_top1_max="${ALIGNMENT_REL_TOP1_MAX:-}"
alignment_pk_free_only="${ALIGNMENT_PK_FREE_ONLY:-0}"
alignment_pk_only="${ALIGNMENT_PK_ONLY:-0}"
alignment_require_zero_mismatch="${ALIGNMENT_REQUIRE_ZERO_MISMATCH:-0}"
skip_log="${ALIGNMENT_SKIP_LOG:-}"

log_skip() {
  if [[ -n "$skip_log" ]]; then
    printf '%s\t%s\n' "$case_id" "$1" >> "$skip_log"
  fi
}

if [[ ! -x "$legacy_bin" ]]; then
  echo "error: legacy bin not executable: $legacy_bin" >&2
  exit 2
fi
if [[ ! -x "$tool_bin" ]]; then
  echo "error: viterbi tool not executable: $tool_bin" >&2
  exit 2
fi
if [[ ! -f "$param_file" ]]; then
  echo "error: param file not found: $param_file" >&2
  exit 2
fi
if [[ ! -f "$GENERATOR" ]]; then
  echo "error: missing generator script: $GENERATOR" >&2
  exit 2
fi
if [[ ! -x "$EXPORTER" ]]; then
  echo "error: missing replay exporter: $EXPORTER" >&2
  exit 2
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

abs_file="$tmpdir/abs_diffs.txt"
rel_file="$tmpdir/rel_diffs.txt"
touch "$abs_file" "$rel_file"
if [[ -n "$diag_csv" ]]; then
  mkdir -p "$(dirname "$diag_csv")"
  printf 'case_id,length,topology_family,seq,G,legacy_mfe_structure,legacy_mfe_energy,scfg_viterbi_energy,scfg_viterbi_structure,structure_mismatch,abs_diff,rel_diff,mismatch\n' > "$diag_csv"
fi
if [[ -n "$skip_log" ]]; then
  mkdir -p "$(dirname "$skip_log")"
  printf 'case_id\treason\n' > "$skip_log"
fi

alignment_compared=0
alignment_mismatched=0
structure_compared=0
structure_mismatched=0
skipped=0
finite_count=0
nonzero_energy_count=0

IFS=',' read -r -a lengths <<<"$lengths_csv"
for length in "${lengths[@]}"; do
  if [[ -z "$length" ]]; then
    continue
  fi
  case_id_prefix="rand${length}"
  cases_file="$tmpdir/cases_${length}.tsv"
  if ! python3 "$GENERATOR" --seed "$seed" --count "$count" --length "$length" --prefix "$case_id_prefix" > "$cases_file"; then
    echo "error: failed to generate cases for length $length" >&2
    exit 1
  fi

  while IFS=$'\t' read -r case_id seq g; do
    if [[ -z "$seq" || -z "$g" ]]; then
      skipped=$((skipped + 1))
      log_skip "empty_seq_or_g"
      continue
    fi

    artifact_file="$tmpdir/legacy_artifact_${case_id}.kv"
    if ! "$EXPORTER" "$legacy_bin" "$param_file" "$seq" "$g" --case-id "$case_id" --format kv >"$artifact_file" 2>/dev/null; then
      skipped=$((skipped + 1))
      log_skip "legacy_failed"
      continue
    fi

    legacy_seq=$(sed -n 's/^seq=//p' "$artifact_file" | head -n1)
    legacy_mfe_structure=$(sed -n 's/^legacy_reported_mfe_r=//p' "$artifact_file" | head -n1)
    legacy_mfe_energy=$(sed -n 's/^legacy_reported_mfe_energy=//p' "$artifact_file" | head -n1)
    if [[ -z "$legacy_seq" || -z "$legacy_mfe_structure" || -z "$legacy_mfe_energy" ]]; then
      skipped=$((skipped + 1))
      log_skip "legacy_parse_failed"
      continue
    fi

    if [[ "$legacy_seq" != "$seq" ]]; then
      alignment_mismatched=$((alignment_mismatched + 1))
      log_skip "legacy_seq_mismatch"
      continue
    fi

    invalid_chars=$(printf '%s' "$legacy_mfe_structure" | tr -d '.()[]')
    if [[ -n "$invalid_chars" ]]; then
      skipped=$((skipped + 1))
      log_skip "legacy_mfe_invalid_chars"
      continue
    fi

    topology_family=$(topology_family_of "$legacy_mfe_structure")
    if [[ "$alignment_pk_free_only" == "1" && "$topology_family" != "pk_free" ]]; then
      skipped=$((skipped + 1))
      log_skip "pk_free_only"
      continue
    fi
    if [[ "$alignment_pk_only" == "1" && "$topology_family" == "pk_free" ]]; then
      skipped=$((skipped + 1))
      log_skip "pk_only"
      continue
    fi
    if [[ "$alignment_pk_only" == "1" && "$topology_family" == "none" ]]; then
      skipped=$((skipped + 1))
      log_skip "pk_only_none"
      continue
    fi

    tool_stdout="$tmpdir/tool_${case_id}.out"
    tool_stderr="$tmpdir/tool_${case_id}.err"
    if ! "$tool_bin" "$seq" "$g" "$param_file" >"$tool_stdout" 2>"$tool_stderr"; then
      skipped=$((skipped + 1))
      log_skip "viterbi_failed"
      continue
    fi
    tool_output=$(tr -d '\r' < "$tool_stdout")
    if [[ "$tool_output" == *$'\t'* ]]; then
      viterbi_structure=$(printf '%s' "$tool_output" | awk -F $'\t' 'NR==1 {print $1}')
      viterbi_energy=$(printf '%s' "$tool_output" | awk -F $'\t' 'NR==1 {print $2}')
    else
      viterbi_energy=$(printf '%s' "$tool_output" | tr -d '\n')
      viterbi_structure=$(grep -m1 -E '^\[scfg-viterbi-structure\]' "$tool_stderr" | sed -nE 's/.* structure=([^ ]+).*/\1/p' || true)
    fi

    if ! python3 "$MATH_TOOL" check_finite "$legacy_mfe_energy" "$viterbi_energy"; then
      skipped=$((skipped + 1))
      log_skip "non_finite_energy"
      continue
    else
      finite_count=$((finite_count + 1))
    fi

    if python3 "$MATH_TOOL" check_nonzero "$legacy_mfe_energy"; then
      nonzero_energy_count=$((nonzero_energy_count + 1))
    fi

    abs_diff=$(python3 "$MATH_TOOL" abs_diff "$legacy_mfe_energy" "$viterbi_energy")
    scale=$(python3 "$MATH_TOOL" scale "$legacy_mfe_energy" "$viterbi_energy")
    rel_diff=$(python3 "$MATH_TOOL" rel_diff "$abs_diff" "$scale")

    printf '%s\n' "$abs_diff" >> "$abs_file"
    printf '%s\n' "$rel_diff" >> "$rel_file"

    alignment_compared=$((alignment_compared + 1))

    mismatch=0
    if ! python3 "$MATH_TOOL" check_within_tol "$abs_diff" "$scale" "$abs_tol" "$rel_tol"; then
      alignment_mismatched=$((alignment_mismatched + 1))
      mismatch=1
    fi

    structure_mismatch=""
    if [[ -n "$viterbi_structure" ]]; then
      structure_compared=$((structure_compared + 1))
      if [[ "$viterbi_structure" == "$legacy_mfe_structure" ]]; then
        structure_mismatch=0
      else
        structure_mismatch=1
        structure_mismatched=$((structure_mismatched + 1))
      fi
    fi

    if [[ -n "$diag_csv" ]]; then
      printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$case_id" "$length" "$topology_family" "$seq" "$g" "$legacy_mfe_structure" "$legacy_mfe_energy" \
        "$viterbi_energy" "${viterbi_structure:-}" "${structure_mismatch:-}" \
        "$abs_diff" "$rel_diff" "$mismatch" >> "$diag_csv"
    fi
  done < <(tail -n +2 "$cases_file")

done

read -r median_abs p90_abs p95_abs <<< "$(python3 "$MATH_TOOL" stats "$abs_file")"
read -r median_rel p90_rel p95_rel <<< "$(python3 "$MATH_TOOL" stats "$rel_file")"
top3_abs=$(python3 "$MATH_TOOL" top3 "$abs_file")
top3_rel=$(python3 "$MATH_TOOL" top3 "$rel_file")
top1_rel=$(printf '%s\n' "$top3_rel" | awk -F',' 'NR==1 {print $1}')

finite_rate=0.0
if [[ "$alignment_compared" -gt 0 ]]; then
  finite_rate=$(python3 "$MATH_TOOL" finite_rate "$finite_count" "$alignment_compared")
fi

echo "alignment_compared=$alignment_compared"
echo "alignment_mismatched=$alignment_mismatched"
echo "structure_compared=$structure_compared"
echo "structure_mismatched=$structure_mismatched"
echo "skipped=$skipped"
echo "finite_rate=${finite_rate}%"
echo "nonzero_energy_count=$nonzero_energy_count"
echo "abs_top3_max=$top3_abs"
echo "rel_top3_max=$top3_rel"
echo "abs_median=$median_abs"
echo "abs_p90=$p90_abs"
echo "abs_p95=$p95_abs"
echo "rel_median=$median_rel"
echo "rel_p90=$p90_rel"
echo "rel_p95=$p95_rel"
echo "alignment_contract=report"
echo "alignment_min=$alignment_min"
if [[ -n "$diag_csv" ]]; then
  echo "alignment_diag_csv=$diag_csv"
  python3 "$DIAG_TOOL" "$diag_csv"
fi

if (( alignment_compared < alignment_min )); then
  echo "alignment report gate failed: alignment_compared=$alignment_compared (<$alignment_min)" >&2
  exit 1
fi

if [[ "$alignment_require_zero_mismatch" == "1" ]] && (( alignment_mismatched != 0 )); then
  echo "alignment report gate failed: alignment_mismatched=$alignment_mismatched (expected 0)" >&2
  exit 1
fi

if ! python3 "$MATH_TOOL" check_rate_ge "$finite_rate" 100.0; then
  echo "alignment report gate failed: finite_rate=${finite_rate}% (expected 100%)" >&2
  exit 1
fi

if [[ -n "$alignment_rel_top1_max" ]] && ! python3 "$MATH_TOOL" check_le "$top1_rel" "$alignment_rel_top1_max"; then
  echo "alignment report gate failed: rel_top3_max_first=$top1_rel (expected <= $alignment_rel_top1_max)" >&2
  exit 1
fi

if ! python3 "$MATH_TOOL" check_thresholds \
  "$median_abs" "$p90_abs" "$p95_abs" "$median_rel" "$p90_rel" "$p95_rel" \
  "$alignment_abs_median_max" "$alignment_abs_p90_max" "$alignment_abs_p95_max" \
  "$alignment_rel_median_max" "$alignment_rel_p90_max" "$alignment_rel_p95_max"; then
  echo "alignment report gate failed: distribution threshold exceeded" >&2
  exit 1
fi
