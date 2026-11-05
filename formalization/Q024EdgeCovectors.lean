import Mathlib.LinearAlgebra.Matrix.Adjugate
import Mathlib.Data.Int.Basic

/-!
A narrowly scoped algebraic interface for Q024's exact direction enumerator.
The edge vectors are COLUMNS of M, the covector u is a ROW, and alpha = u M.
No positivity of det M is assumed. No polytope, width, HNF, or Python-code
correctness theorem is claimed by this file.
-/
namespace Q024
open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Exact integer reconstruction numerator for a row covector. -/
theorem edge_numerators (M : Matrix n n ℤ) (u : n → ℤ) :
    (u ᵥ* M) ᵥ* M.adjugate = M.det • u := by
  rw [vecMul_vecMul, M.mul_adjugate, vecMul_smul, vecMul_one]

/-- All integer solutions to u M = alpha, including the non-unimodular case,
are detected by coordinatewise integer divisibility of alpha adj(M). -/
theorem integral_covector_exists_iff (M : Matrix n n ℤ) (hM : M.det ≠ 0)
    (alpha : n → ℤ) :
    (∃ u : n → ℤ, u ᵥ* M = alpha) ↔
      ∀ j, M.det ∣ (alpha ᵥ* M.adjugate) j := by
  constructor
  · rintro ⟨u, rfl⟩ j
    rw [edge_numerators]
    exact ⟨u j, rfl⟩
  · intro h
    choose u hu using h
    have hn : alpha ᵥ* M.adjugate = M.det • u := by
      funext j
      exact hu j
    refine ⟨u, ?_⟩
    have he : M.det • (u ᵥ* M) = M.det • alpha := by
      rw [← smul_vecMul, ← hn, vecMul_vecMul, M.adjugate_mul,
        vecMul_smul, vecMul_one]
    funext j
    exact mul_left_cancel₀ hM (congrFun he j)

/-- Nonzero determinant gives uniqueness without requiring determinant ±1. -/
theorem integral_covector_unique (M : Matrix n n ℤ) (hM : M.det ≠ 0)
    (u v : n → ℤ) (h : u ᵥ* M = v ᵥ* M) : u = v := by
  have he := congrArg (fun a => a ᵥ* M.adjugate) h
  rw [edge_numerators, edge_numerators] at he
  funext j
  exact mul_left_cancel₀ hM (congrFun he j)

/-- The coordinate formula is integer Euclidean division; on divisible
numerators it is exact for positive or negative nonzero determinant. -/
theorem integral_covector_reconstruction (M : Matrix n n ℤ) (hM : M.det ≠ 0)
    (alpha u : n → ℤ) (h : u ᵥ* M = alpha) (j : n) :
    u j = (alpha ᵥ* M.adjugate) j / M.det := by
  rw [← h, edge_numerators]
  simp only [Pi.smul_apply, smul_eq_mul]
  exact (Int.mul_ediv_cancel_left (u j) hM).symm

/-- Combined existence and uniqueness criterion at the enumerator interface. -/
theorem integral_covector_existsUnique_iff (M : Matrix n n ℤ) (hM : M.det ≠ 0)
    (alpha : n → ℤ) :
    (∃! u : n → ℤ, u ᵥ* M = alpha) ↔
      ∀ j, M.det ∣ (alpha ᵥ* M.adjugate) j := by
  constructor
  · intro h
    exact (integral_covector_exists_iff M hM alpha).mp h.exists
  · intro h
    obtain ⟨u, hu⟩ := (integral_covector_exists_iff M hM alpha).mpr h
    refine ⟨u, hu, ?_⟩
    intro v hv
    exact integral_covector_unique M hM v u (hv.trans hu.symm)

/-- The exact vector returned after the divisibility filter solves the edge
value equation. This is algebraic correctness, not a proof of Python execution. -/
theorem integral_reconstruction_solves (M : Matrix n n ℤ) (hM : M.det ≠ 0)
    (alpha : n → ℤ) (h : ∀ j, M.det ∣ (alpha ᵥ* M.adjugate) j) :
    (fun j => (alpha ᵥ* M.adjugate) j / M.det) ᵥ* M = alpha := by
  obtain ⟨u, hu⟩ := (integral_covector_exists_iff M hM alpha).mpr h
  have hr : (fun j => (alpha ᵥ* M.adjugate) j / M.det) = u := by
    funext j
    exact (integral_covector_reconstruction M hM alpha u hu j).symm
  rw [hr, hu]

end Q024

#print axioms Q024.edge_numerators
#print axioms Q024.integral_covector_exists_iff
#print axioms Q024.integral_covector_unique
#print axioms Q024.integral_covector_reconstruction

#print axioms Q024.integral_covector_existsUnique_iff
#print axioms Q024.integral_reconstruction_solves
