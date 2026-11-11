# REFEREE4 (8 Nov 2026, late): independent agent referee + follow-up check — FATAL GAP FOUND

## Agent referee (one subagent, unanchored): verdict MAJOR REVISION, 17 issues
Serious: (1) λ conventions inconsistent; (2) angular inequality H ≥ H0 + E not derived; (3) adjoint identity only
sampled; (4) Ψ4 bookkeeping in the reconstruction lemma wrong as written; (5) divergent Euler integrals (need
Pochhammer/Hankel contours); (11) Theorem B uses T with r2 instead of HPV's critical radius; (12) δ(Ψ2φ³) ≡ 0
(wrong invariant in the gluing step); (13) acceleration/axis steps unverifiable; (15) generalized-mode correction
term wrong.
Checked correct by the agent: window algebra, currents, R1/R2, R3 saturate/solve (spot checks of N and of the seam),
the upper-half-plane architecture, Goldberg–Sachs/reality/ξ/ζ, Iyer–Wald.

## Follow-up on issue (2): the angular bound λ ≥ 3/2 is FALSE (fatal for Theorem A as written)
* step128_angular_identity.py (exact): with the Suzuki–Takasugi–Umetsu angular operator (same λ_s as the radial
  equation used in step39), H − H0 − K_α = α(m+s)²/(2(1−x)) + α(m−s)²/(2(1+x)) + E + αs
  + (1+α)m(αm − 4(1+α)ξ)/(1+αx²). Relative to the E certified in step47, the net mξ coefficient is −2, not +2.
  step47 certified the wrong operator, and its box ξ ∈ [0, |m|] makes ξ look favourable.
* angular_check/ang.py (Galerkin in spin-weighted harmonics; floating point). Kerr check: λ(ξ=0) = 2 ✓, and λ drops
  like 2 − 6.7ξ for (m, s) = (2, ±2), matching first-order spheroidal perturbation theory.
* angular_check/win.py, chk.py: near-extremal KdS (roots r0 = 0.9–0.99, r1 = 1, r2 = 10), top of the window,
  m = 2: λ ≈ −2.4 … −4.1. Then B = λC0 + c1 lies far below the floor 5.5t − 1.25, and even L(t) < 0, so the
  cosmological-side current (f₊ = 1) also fails.
* CONSEQUENCE: Theorem A for real σ inside the window is NOT proved for near-extremal parameters. The event lemma's
  hypothesis B ≥ 5.5t − 1.25 does not hold there, and Prop. 3.7's L(t) > 0 fails. The upper-half-plane argument and
  Theorem B inherit this gap. The manuscript's "unconditional" claims are withdrawn until this is repaired.
* Earlier "λ ≥ 2 numerically" (LEMMAS items 27, 45) used a different λ convention or the same sign error; void.

## What a repair needs
A version of the dual two-current argument valid when L(t) ≤ 0 and B is small or negative. This is the near-extremal,
strongly superradiant regime. It is a genuine new step, not bookkeeping.
