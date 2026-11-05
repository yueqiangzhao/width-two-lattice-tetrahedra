#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")"
for f in evidence/results/diagonal_comparison.json evidence/results/enum_22k.json; do
  gzip -dc -- "${f}.gz" > "$f"
done
sha256sum -c ORIGINAL_PAYLOAD_SHA256SUMS
