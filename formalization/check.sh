#!/usr/bin/env bash
# Read-only verification with an existing, matching Lean/mathlib installation.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
: "${LEAN_BIN:?Set LEAN_BIN to the Lean 4.35.0 executable}"
: "${MATHLIB_ROOT:?Set MATHLIB_ROOT to the existing mathlib checkout}"
EXPECTED=56a005855ea17e0b0958f73ec801ae140b0a1679
LEAN_VERSION="$("$LEAN_BIN" --version)"
printf '%s\n' "$LEAN_VERSION" | grep -F 'version 4.35.0,' >/dev/null
printf '%s\n' "$LEAN_VERSION" | grep -F 'commit eb6cafb7290e8bd5f2937522d8d804f58f02bafc,' >/dev/null
[[ "$(git -C "$MATHLIB_ROOT" rev-parse HEAD)" == "$EXPECTED" ]] || { echo 'mathlib revision mismatch' >&2; exit 2; }
export LEAN_PATH="$MATHLIB_ROOT/.lake/build/lib/lean:$(find "$MATHLIB_ROOT/.lake/packages" -path '*/lib/lean' -type d | paste -sd: -)"
"$LEAN_BIN" "$ROOT/Q024EdgeCovectors.lean"
