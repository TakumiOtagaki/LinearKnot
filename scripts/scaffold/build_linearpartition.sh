#!/usr/bin/env bash
# Build LinearPartition (Zhang et al. 2020) for scaffold construction.
#   scripts/scaffold/build_linearpartition.sh [out_dir]
# Output: <out_dir>/LinearPartition/linearpartition (default out_dir:
# build/third_party). LINEARPARTITION_URL overrides the clone source.
set -euo pipefail

url="${LINEARPARTITION_URL:-https://github.com/LinearFold/LinearPartition.git}"
commit="${LINEARPARTITION_COMMIT:-b450fb3}"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
out="${1:-$root/build/third_party}"
src="$out/LinearPartition"

if [[ ! -d "$src/.git" ]]; then
  mkdir -p "$out"
  git clone -q "$url" "$src"
fi
git -C "$src" checkout -q --detach "$commit"
make -C "$src" -s linearpartition
echo "linearpartition: $src/linearpartition"
