# Q024 scoped Lean interface and applicability

## Actual assessment and result

The fixed existing environment was inspected on 28 October 2026: Lean 4.35.0 (compiler commit eb6cafb7290e8bd5f2937522d8d804f58f02bafc), mathlib commit 56a005855ea17e0b0958f73ec801ae140b0a1679. This was not a conclusion based merely on the absence of a manuscript Lean project.

The relevant library has general linear groups (`LinearAlgebra/Matrix/GeneralLinearGroup/Defs.lean`), integer lattices (`Algebra/Module/ZLattice/Basic.lean`), linear independence and ranks, and the exact matrix identities `Matrix.mul_adjugate` and `Matrix.adjugate_mul` (`LinearAlgebra/Matrix/Adjugate.lean`, lines 263 and 268 in this revision). Searches for Hermite normal form, successive minima, lattice width, and multiwidth did not locate a directly usable manuscript-level interface. The integer-lattice and matrix APIs do not alone supply Q024's width-specific saturation, finite box, or HNF orbit argument.

One high-value small interface did have a faithful direct translation. `Q024EdgeCovectors.lean` now proves six statements by Lean kernel checking: the adjugate numerator identity, an integer-solvability divisibility criterion, uniqueness, exact coordinate reconstruction by integer division, the combined unique-existence criterion, and correctness of the reconstructed vector. Actual compiler exit was 0. Every printed axiom list contains only `propext`, `Classical.choice`, and `Quot.sound`; no sorry, admitted theorem, or additional axiom is used. The checker does not use native_decide or an external oracle.

## Exact source mapping

The mathematical interface is

    M : Matrix n n Int, det(M) != 0,
    alpha : n -> Int,
    (exists unique integral row u, u M = alpha)
        iff det(M) divides every coordinate of alpha adj(M).

It applies to arbitrary finite square dimension, in particular the manuscript's 3-dimensional tetrahedra and 2-dimensional planar configurations. It does NOT require determinant +1, -1, or positive determinant; negative and nonunit determinants are covered. The reconstruction theorem uses exact integer division of divisible numerators.

This matches the manuscript's `main.tex` lines 141–149 and `sections/boundary.tex` lines 535–545 in the unchanged user-source baseline. In the unchanged `outputs/Q024/enumerate.py`, `edge_matrix` lines 35–38 stores edges as COLUMNS. `short_covectors` lines 107–124 explicitly uses the ROW convention `u E = a`, computes numerator coordinate `sum(values[i] * adj[i][j])`, filters divisibility by `det`, then returns the integer quotient. Thus there is no transposition or column-covector substitution in the Lean statement. The Lean file proves this algebraic interface, not semantic equivalence of the Python determinant/adjugate implementations, `%`, `//`, or program execution with Lean. On divisible numerators the intended quotient is exact regardless of the determinant sign.

## What remains outside this formal result

* The all-integral-direction multi-width definition, affine invariance, width norm, independent-direction rank criterion, and the saturated short-plane reduction are not formalized here. Real/integer lattice APIs alone would not justify replacing independent realizing directions by a lattice basis.
* The geometric reduction to the 13 planar projections and the finite height boxes for k = 2, 3, 4 is not formalized. The algebraic recovery lemma is one necessary interface, not the full bounded-to-unbounded geometric bridge.
* The finite edge-value enumeration, oscillation bounds, sign quotient, rank tests, exact counts 17/45/47, and the Python-to-Lean data/enumeration certificate bridge are not formalized.
* GL3(Z) and integer matrices exist in mathlib, but no ready HNF normal-form/classification bridge was located for the 24-vertex-order canonical key used here. That algorithm, the complete affine-orbit criterion, and the representative classification are not formalized.
* The strict and diagonal parameter classifications and the nine unrestricted counting formulas are not formalized.

No manuscript mathematical text or original computational script was changed by this supplement. This is a partial algebraic formalization only, not a formal verification of Q024, a new universal proof verdict, or human review. Kernel compilation alone does not establish source-to-manuscript correspondence; that correspondence is a separate check.

## Reproduce

With an existing compiled mathlib checkout at the exact revision and the matching Lean executable, set `MATHLIB_ROOT` and `LEAN_BIN`, then run `bash check.sh`. The script checks those versions, adds only the existing mathlib/package compiled-library paths, and invokes Lean on the source. It does not download dependencies, install software, or modify the shared mathlib checkout. `compile.log` is the actual successful compiler output including the axiom prints; `compile.exit` records 0. These are scoped compiler evidence, not a complete-manuscript gate.
