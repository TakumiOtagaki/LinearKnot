#!/usr/bin/env bash
# Build the legacy CParty oracle: upstream CParty at a pinned commit plus the
# patches in patches/ (see README.md). Output:
#   <out>/src            patched upstream checkout
#   <out>/cparty-oracle  wrapper used by the parity scripts
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
out="${1:-$root/build/legacy_cparty_oracle}"
src="$out/src"
mapfile -t config < <(python3 "$here/oracle_config.py" config)
[[ ${#config[@]} -ge 3 ]] || { echo "invalid oracle configuration" >&2; exit 2; }
UPSTREAM_URL="${config[0]}"
UPSTREAM_COMMIT="${config[1]}"
patches=()
for name in "${config[@]:2}"; do patches+=("$here/patches/$name"); done

if [[ ! -d "$src/.git" ]]; then
  mkdir -p "$out"
  git clone -q "$UPSTREAM_URL" "$src"
fi
if ! git -C "$src" cat-file -e "$UPSTREAM_COMMIT^{commit}" 2>/dev/null; then
  git -C "$src" fetch -q origin
fi
# src is the disposable generated checkout; discard the previous patch set
# before switching commits (a plain checkout rejects those staged patches).
git -C "$src" checkout -q --detach --force "$UPSTREAM_COMMIT"
git -C "$src" clean -qfdx
git -C "$src" apply --index "${patches[@]}"

cmake -S "$src" -B "$src/build" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$src/build" -j8
python3 "$here/oracle_config.py" record "$src" "$out/manifest.json"

# Sampling (used only for BPP / MEA / centroid output) can abort with
# "backtracking impossible ..." and exit before the energies are printed.
# The parity scripts read only the energies, so the wrapper disables it.
cat > "$out/cparty-oracle" <<WRAP
#!/usr/bin/env bash
exec "$src/build/CParty" -s 0 "\$@"
WRAP
chmod +x "$out/cparty-oracle"
echo "legacy oracle: $out/cparty-oracle"
