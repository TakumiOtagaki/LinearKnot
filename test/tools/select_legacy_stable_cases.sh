#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: select_legacy_stable_cases.sh <legacy_cparty_bin> <seed> <count> <lengths> <out_tsv>

Generates (seq, G) cases from pk_free templates in multi.secstruct and keeps
only those that legacy CParty executes successfully and whose stdout parses.

Lengths is a comma-separated list (e.g. "30,50,80").
Outputs TSV with header: case_id  seq  G
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 5 ]]; then
  usage
  exit 2
fi

legacy_bin="$1"
seed="$2"
count="$3"
lengths_csv="$4"
out_tsv="$5"

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PARSER="$ROOT_DIR/test/tools/parse_cparty_stdout.sh"
SEED_DATASET="${CPARTY_MULTI_SECSTRUCT:-$ROOT_DIR/test/multi.secstruct}"
PARAM_FILE="$ROOT_DIR/params/rna_DirksPierce09.par"

if [[ ! -x "$legacy_bin" ]]; then
  echo "error: legacy binary not executable: $legacy_bin" >&2
  exit 1
fi
if [[ ! -x "$PARSER" ]]; then
  echo "error: missing parser script: $PARSER" >&2
  exit 1
fi
if [[ ! -f "$SEED_DATASET" ]]; then
  echo "error: seed dataset not found: $SEED_DATASET" >&2
  exit 1
fi
if [[ ! -f "$PARAM_FILE" ]]; then
  echo "error: missing parameter file: $PARAM_FILE" >&2
  exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

stdout_legacy="$tmpdir/legacy.txt"
cases_tsv="$tmpdir/cases.tsv"

IFS=',' read -r -a lengths <<< "$lengths_csv"
if [[ ${#lengths[@]} -eq 0 ]]; then
  echo "error: no lengths provided" >&2
  exit 1
fi

mkdir -p "$(dirname "$out_tsv")"
printf 'case_id\tseq\tG\n' > "$out_tsv"

total=0
kept=0
skipped=0

python3 - "$SEED_DATASET" "$seed" "$count" "$lengths_csv" > "$cases_tsv" <<'PY'
import csv
import random
import re
import sys

seed_dataset = sys.argv[1]
rng_seed = int(sys.argv[2])
count = int(sys.argv[3])
lengths_csv = sys.argv[4]

lengths = [int(x) for x in lengths_csv.split(",") if x]
if not lengths:
    raise SystemExit("error: no lengths provided")

allowed_pairs = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")]
alphabet = ("A", "U", "G", "C")

templates = []
with open(seed_dataset, encoding="utf-8") as fh:
    chunks = fh.read().strip().split("\n\n")
for chunk in chunks:
    lines = [line.strip() for line in chunk.splitlines() if line.strip()]
    if len(lines) < 3 or not lines[0].startswith(">"):
        continue
    header = lines[0][1:]
    seq = lines[1]
    structure = lines[2]

    metadata = {}
    for token in header.split("|")[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        metadata[key.strip()] = value.strip()

    if metadata.get("expected") != "valid":
        continue
    if metadata.get("family") != "pk_free":
        continue
    if re.search(r"[^().]", structure):
        continue

    case_id = header.split("|", 1)[0].strip()
    templates.append((case_id, structure))

if not templates:
    raise SystemExit("error: no pk_free templates found in seed dataset")

min_len = min(len(st) for _, st in templates)
for L in lengths:
    if L < min_len:
        raise SystemExit(f"error: requested length {L} < minimum template length {min_len}")

rng = random.Random(rng_seed)

def randomize_sequence(structure: str) -> str:
    seq = ["?"] * len(structure)
    stack = []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("unbalanced structure")
            j = stack.pop()
            b1, b2 = rng.choice(allowed_pairs)
            seq[j] = b1
            seq[i] = b2
        elif ch == ".":
            seq[i] = rng.choice(alphabet)
        else:
            raise ValueError(f"unsupported char: {ch}")
    if stack:
        raise ValueError("unbalanced structure")
    return "".join(seq)

TURN = 3

def enforce_turn(structure: str, turn: int) -> str:
    """Replace too-short pairs with dots so all pairs satisfy j - i - 1 >= turn."""
    if turn <= 0:
        return structure
    chars = list(structure)
    stack = []
    for i, ch in enumerate(chars):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("unbalanced structure")
            j = stack.pop()
            if i - j - 1 < turn:
                chars[j] = "."
                chars[i] = "."
        elif ch != ".":
            raise ValueError(f"unsupported char: {ch}")
    if stack:
        raise ValueError("unbalanced structure")
    return "".join(chars)

def build_structure(target_len: int) -> str:
    structure_parts = []
    total = 0
    while total < target_len:
        remaining = target_len - total
        if remaining < min_len:
            structure_parts.append("." * remaining)
            total = target_len
            break

        _, tmpl = rng.choice(templates)
        if len(tmpl) > remaining:
            continue

        structure_parts.append(tmpl)
        total += len(tmpl)

        if total < target_len and rng.random() < 0.3:
            max_gap = min(3, target_len - total)
            gap = rng.randint(0, max_gap)
            if gap:
                structure_parts.append("." * gap)
                total += gap

    structure = "".join(structure_parts)
    structure = enforce_turn(structure, TURN)
    return structure

writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
writer.writerow(["case_id", "seq", "G"])
for L in lengths:
    for i in range(1, count + 1):
        structure = build_structure(L)
        if "(" not in structure:
            raise SystemExit("error: generated structure without base pairs")
        seq = randomize_sequence(structure)
        writer.writerow([f"len{L}_rand_{i:04d}", seq, structure])
PY

while IFS=$'\t' read -r case_id seq g; do
  if [[ "$case_id" == "case_id" ]]; then
    continue
  fi
  total=$((total + 1))

  if ! "$legacy_bin" -d2 -P "$PARAM_FILE" -r "$g" "$seq" > "$stdout_legacy" 2>/dev/null; then
    skipped=$((skipped + 1))
    continue
  fi
  if ! "$PARSER" "$stdout_legacy" >/dev/null 2>&1; then
    skipped=$((skipped + 1))
    continue
  fi

  printf '%s\t%s\t%s\n' "$case_id" "$seq" "$g" >> "$out_tsv"
  kept=$((kept + 1))
done < "$cases_tsv"

echo "select_total=$total"
echo "select_kept=$kept"
echo "select_skipped=$skipped"
echo "select_output=$out_tsv"

if (( kept < 1 )); then
  echo "error: no stable cases found" >&2
  exit 1
fi
