# Computational verification materials for width-two lattice tetrahedra

Supplement to *Width-two lattice tetrahedra and their multi-widths*, by **Yueqiang Zhao**, prepared for submission to *Mathematika* (5 November 2026).

This public repository contains exact Python enumeration/checker scripts, recorded finite certificates, and a **limited** Lean 4 formalization. The mathematical proof for arbitrary multi-width parameters is in the manuscript; finite computations and six formalized covector facts do **not** constitute a Lean proof of the main classification.

## Files

- `outputs/Q024/` — the original enumeration and classification scripts (Python standard library).
- `work/q024/` — independent strict/diagonal checkers.
- `evidence/results/` — precomputed exact finite checks and certificates. Two large JSON files are stored as losslessly compressed `.json.gz` files to make the public snapshot lightweight.
- `formalization/` — Lean code, historical compiler outputs, and an exact scope statement (`SCOPE.md`).
- `evidence/COMPUTATION_SCOPE.json` — actual tested parameter ranges and limits.
- `replay.sh` — replay script; **requires the third-party dataset described below** for complete execution.
- `SHA256SUMS` — checksums for the committed repository files (before restoring compressed JSON).
- `ORIGINAL_PAYLOAD_SHA256SUMS` — checksums of the original uncompressed code and data (excluding the README, which was adapted for GitHub).

## Verify and restore published certificates

From the repository root:

```bash
sha256sum -c SHA256SUMS
bash restore_large_json.sh
```

The second command expands `evidence/results/diagonal_comparison.json.gz` and `evidence/results/enum_22k.json.gz` **without altering the JSON bytes** and verifies the original payload against `ORIGINAL_PAYLOAD_SHA256SUMS`.

## Run the complete finite replay

The independent comparison uses G. Hamm's version-2 tetrahedron dataset, **not redistributed here**:

https://doi.org/10.5281/zenodo.10050853

After obtaining the corresponding original dataset, copy it to:

`work/q024/sources/width_2_tetrahedra_v2.txt`

The recorded expected SHA-256 of the source dataset is:

`3a12fb3f19a3ecdc2bb8500dbf67ab84a615af6f54ce2bd906eb19285023fd0a`

Then run:

```bash
python3 --version
bash replay.sh replay
```

The replay creates/updates finite JSON results. The source dataset needs to be present for complete reproduction; a successful finite run supports the **specific tested parameters**, not the general symbolic proof.

## Formalization boundary

The `formalization/` directory covers six algebraic statements on integral covectors / adjugate reconstruction. The historical Lean compile record reports success with the specified Lean/mathlib revisions, but the formalization was not rerun for this GitHub publication. It does not cover saturation of short directions, multi-width classification, normal form uniqueness, all finite enumeration implementations, or the nine-branch general formula.

## Public snapshot and licensing

Prepared as supplementary source material for a *Mathematika* submission. Refer to a particular Git commit hash when citing the exact submitted snapshot. No repository-wide license has yet been selected by the author; do not infer open-source reuse rights solely from public GitHub availability. Original third-party data remain hosted by their own publisher.
