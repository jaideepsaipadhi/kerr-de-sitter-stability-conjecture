# Bridge Part 2, Steps C–D, redone on Giorgi's Teukolsky system (A, F, B, X)

Sources: Giorgi, arXiv:2002.07228: Thm 6.1 (the system T_s, M_s) and Prop 5.7 (the relations (75)–(78)).
This replaces the DGS-based version, whose principal part is non-diagonal (see `correction_principal_symbol.md`).

## Setting (fixed mode)

Fix a mode e^{-iωt+imφ} with ω real and nonzero. The horizontal tensors A, F, B, X become spin-weighted functions of (r, θ), and the operators reduce as follows:

- T_s: elliptic in (r, θ) on the exterior, with **diagonal** principal part Δ∂_r² + (angular Laplacian).
- M_s: of order ≤ 1 in (A, F, B, X).
- (c)D̂⊗ acting on 1-tensors: a first-order ODE operator in θ, depending on r. At a = 0 it reduces to ð from spin weight 1 to spin weight 2.

## Eliminating X (Lemma X)

Relation (76) reads

  (c)D̂⊗X + ½(3H + H̄)⊗X = −2[(c)∇₄F + (3/2 trX + 1/2 tr X̄)F] − 2(F)P·A.

Write 𝒟_r for the operator on the left, at fixed r and fixed (ω, m).

- At a = 0, 𝒟_r = ð : spin 1 → spin 2. Its kernel is the ℓ = 1 spin-1 harmonics, one-dimensional when |m| ≤ 1 and zero otherwise. Its cokernel is the kernel of ð̄ on spin 2, which is zero.
- For a ≠ 0, 𝒟_r is a compact (lower-order) deformation, so the Fredholm index stays 1 for |m| ≤ 1 and 0 otherwise. **Assume (H_X): the cokernel of 𝒟_r is zero for every r > r₊.** Then the kernel has dimension exactly equal to the index, and

  X = 𝒢_r[F, ∂F, A] + k(r)·e_r(θ),

  where 𝒢_r is a right inverse of order ≤ 0 in θ. Acting on ∂F it gives an operator of order ≤ 1 in total. The second term, k(r)·e_r(θ), lies in the kernel and appears only when |m| ≤ 1.
- The scalar k(r) is fixed by projecting (78), ∇₃X + ½ tr X̄ X = −D·F − H̄·F − 2B, onto e_r. This gives a first-order ODE in r for k, with a source that is first order in (F, B).

## The reduced system

(A, F, B) satisfy an elliptic system with diagonal principal part, coupled by operators of order ≤ 1. These operators are pseudodifferential in θ, but each acts on a single r-slice. The system is also coupled to the ODE for k(r).

## Unique continuation inward (Step D)

Assume (A, F, B, X) = 0 for r > R. Use the Carleman weight e^{τφ(r)}, with φ chosen so that the standard weighted estimate holds for Δ∂_r² + (angular Laplacian):

  τ³‖e^{τφ}u‖² + τ‖e^{τφ}∇_{r,θ}u‖² ≲ ‖e^{τφ}Pu‖².

- The coupling terms of order ≤ 1 act on each r-slice, so they commute with the weight. Their contribution is bounded by ‖e^{τφ}(u, ∇u)‖ times a constant uniform in r on compact r-intervals. That constant requires the bound on 𝒢_r to be uniform in r, which follows from (H_X) together with continuity in r. These terms are therefore absorbed once τ is large.
- The ODE for k: a one-dimensional Carleman estimate for a first-order ODE gains τ provided φ′ ≠ 0. Its source is first order in (F, B), which is controlled by the τ‖e^{τφ}∇u‖ term.

Hence the vanishing propagates across every sphere r = const, all the way to r₊. This proves Step D.

## Rellich and Carleman at infinity (Step C)

At large r the horizontal structure becomes integrable and 𝒟_r → ð, so (H_X) holds automatically for r ≫ 1. The polynomial-weight argument goes through for the diagonal system plus order-≤1 coupling, as planned.

## Open items

- **(H_X), the real check.** No cokernel of 𝒟_r for any r > r₊ and any mode (ω, m). This is a θ-ODE question for each (r, ω, m). It can fail only at isolated values, and it can be checked numerically, then proved: for example, by a positivity argument for 𝒟_r𝒟_r^*, or by perturbation from a = 0 combined with an argument excluding eigenvalue crossings.
- **The dictionary.** Translate the bridge's hypotheses into Giorgi's variables: the DGS pair vanishing should correspond to (A̲, F̲, B̲) = 0. Giorgi's footnote 8 and Remark 5.5 give the Newman–Penrose forms of f, b, x.
- **Bridge Part 1 in Giorgi's variables.** Either redo it there, or show via the dictionary that the hypotheses are equivalent.
