#!/usr/bin/env bash
# Exact finite replay using unchanged mathematical cores; Python standard library only.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
OUT="${1:-replay}"
mkdir -p -- "$OUT"
PYTHON="${PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1
"$PYTHON" outputs/Q024/enumerate.py --ks 2 3 4 5 6 --output "$OUT/enum_22k.json"
"$PYTHON" outputs/Q024/check_finite.py --finite "$OUT/enum_22k.json" --dataset work/q024/sources/width_2_tetrahedra_v2.txt --output "$OUT/finite_check.json"
"$PYTHON" outputs/Q024/lift_circuits.py --output "$OUT/circuit_certificate.json"
"$PYTHON" outputs/Q024/check_certificates.py --circuit "$OUT/circuit_certificate.json" --finite "$OUT/enum_22k.json" --output "$OUT/computation_certificate.json"
"$PYTHON" work/q024/check_strict_constructor.py
cp work/q024/check_strict_constructor.json "$OUT/strict_comparison.json"
"$PYTHON" work/q024/check_diagonal_constructor.py --output "$OUT/diagonal_comparison.json"
"$PYTHON" outputs/Q024/diagonal_audit_replay.py
cp outputs/Q024/diagonal_audit_replay.json "$OUT/diagonal_lattice_check.json"
printf '%s\n' 'Completed exact finite replay. These checks do not establish the unbounded theorem or a Lean formalization.'
