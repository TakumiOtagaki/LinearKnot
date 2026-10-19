#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: extract_random_case_tsv.sh <input_tsv> <case_id> <output_tsv>

Extracts exactly one case row from a random-case TSV of the form:
  case_id<TAB>seq<TAB>G

and writes a new TSV with the original header plus the selected row.
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

input_tsv="$1"
case_id="$2"
output_tsv="$3"

[[ -f "$input_tsv" ]] || { echo "error: missing input_tsv: $input_tsv" >&2; exit 2; }

python3 - "$input_tsv" "$case_id" "$output_tsv" <<'PY'
import csv
import sys
from pathlib import Path

input_tsv = Path(sys.argv[1])
case_id = sys.argv[2]
output_tsv = Path(sys.argv[3])

with input_tsv.open() as f:
    rows = list(csv.reader(f, delimiter="\t"))

if not rows:
    raise SystemExit(f"error: empty input_tsv: {input_tsv}")

header = rows[0]
matches = [row for row in rows[1:] if row and row[0] == case_id]

if len(matches) != 1:
    raise SystemExit(f"error: expected exactly 1 match for {case_id}, found {len(matches)}")

with output_tsv.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerow(matches[0])

print(f"extracted_case_id={case_id}")
print(f"output_tsv={output_tsv}")
PY
