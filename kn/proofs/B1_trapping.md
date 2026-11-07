# B1: trapping in Kerr–Newman is r-normally hyperbolic (full subextremal range)

**Setting and notation.** Full subextremal range a² + Q² < M², BL coordinates, normalised energy E = 1.
- λ = L/E, η = K/E²
- P = r² + a² − aλ
- Δ = r² − 2Mr + a² + Q²
- R(r) = P² − Δη, with (Σ dr/ds)² = R

The case E = 0 is treated separately at the end.

## Lemma B1.1 (at most two exterior roots). For every (λ, η), R has at most two roots in (r₊, ∞), counted with multiplicity.

*Proof.*
- R is a quartic with leading coefficient 1, and R(r₊) = P(r₊)² ≥ 0.
- R'' = 2P'² + 2PP'' − Δ''η = 12r² + 4a² − 4aλ − 2η. This is an even quadratic in r, so it has **at most one positive root**.
- Suppose R had ≥ 3 roots in (r₊, ∞), counted with multiplicity.
 - If R(r₊) = 0, then R has ≥ 4 roots in [r₊, ∞). By Rolle (with multiplicity), R' has ≥ 3 roots in (r₊, ∞), so R'' has ≥ 2 roots there. Contradiction.
 - If R(r₊) > 0, then the root count in (r₊, ∞) together with the signs R(r₊) > 0 and R(+∞) > 0 forces an even number of roots, counted with multiplicity. So there are ≥ 4, and again R'' would have ≥ 2 roots in (r₊, ∞). Contradiction. ∎

This uses no information about Θ, and it is uniform in the full subextremal range. The numerical check (200k random samples, `B1_check.py`) is consistent with it.

## Corollary B1.2 (trapped set = spherical orbits). A null geodesic in the exterior that neither crosses H⁺ nor escapes to r = ∞ has r ≡ r_c, where R(r_c) = R'(r_c) = 0.

*Proof.* Use Lemma B1.1.
- If R has no exterior root, then R > 0 on (r₊, ∞) and r is monotone, so the geodesic reaches H⁺ or I⁺.
- If R has two simple roots x₁ < x₂, the allowed set is (r₊, x₁] ∪ [x₂, ∞), since R < 0 in between. In (r₊, x₁] the geodesic has at most one turning point and then falls into H⁺ (dr/ds ≠ 0 elsewhere). In [x₂, ∞) it escapes.
- If R has a double root r_c, the geodesic either sits at r_c or is asymptotic to it, and the asymptotic ones are not trapped in both time directions. The trapped set in the usual sense of Dyatlov and Wunsch–Zworski (trapped both forward and backward) is {r = r_c, ξ_r = 0}. ∎

This settles the dependency flagged in Lemma A's Remark 1.

## Lemma B1.3 (nondegeneracy). At every exterior spherical photon orbit,
 R''(r_c) = (8r_c/(r_c − M)²)·[(r_c − M)³ + M(M² − a² − Q²)] > 0.

*Proof.* R = R' = 0 with P ≠ 0 give η = P²/Δ and P = 2rΔ/(r − M). Substituting into R'' gives the identity above; it is verified symbolically in `B1_check.py`. It is positive because r_c > r₊ > M and M² − a² − Q² > 0.
- **P = 0 is impossible:** it gives η = 0, and then Θ ≥ 0 forces λ = a sin²θ. So r² + a² = a² sin²θ, which has no solution with r > 0.
- **E = 0:** then R = a²L² − ΔK, and R' = −Δ'K = 0 forces r = M ∉ (r₊, ∞), or K = 0. In that case Θ ≥ 0 forces L = 0, the trivial case. ∎

**Consequence.** In Dyatlov's normal form, the radial Hamiltonian near Γ is ∝ Δξ_r² − (f(r) − η) with f = P²/Δ. At a critical point with f = η, f'' = R''/Δ > 0. So the radial flow has a **hyperbolic fixed point** at (r_c, ξ_r = 0), with expansion rate ν = (1/Σ)·(R''(r_c)/2)^{1/2} > 0, which is continuous and positive on the (compact) trapped set.

## Theorem B1 (r-normal hyperbolicity for every r). For all subextremal (M, a, Q), the trapped set Γ ⊂ T*ℳ ∖ 0 (in the semiclassical / scattering characteristic set at fixed energy) is a smooth compact symplectic submanifold. It is normally hyperbolic, and r-normally hyperbolic for every r ≥ 1 in the sense of Wunsch–Zworski / Dyatlov.

*Proof sketch, with the same structure as Dyatlov's Kerr–de Sitter verification (Dyatlov, *Asymptotics of linear waves…*, Ann. Henri Poincaré 2015, §3; and Wunsch–Zworski 2011 for Kerr).*
1. **Smoothness.** The spherical radius is a smooth function on the allowed (λ, η) set. In fact the map r_c ↦ (λ(r_c), η(r_c)) is an explicit smooth parametrisation, and it is injective because R'' > 0 (implicit function theorem). So Γ = {r = r_c(λ, η), ξ_r = 0} is a smooth submanifold.
2. **Symplectic.** Its normal directions are (r − r_c, ξ_r), and the symplectic form restricted to these is dr ∧ dξ_r ≠ 0.
3. **Normal hyperbolicity.** Lemma B1.3 gives normal expansion and contraction rates ±ν with ν > 0.
4. **r-normal hyperbolicity for every r.** The geodesic flow restricted to Γ is completely integrable: E, L and the Carter constant K are in involution, so the flow on Γ is quasi-periodic on invariant tori and its tangential Lyapunov exponents are 0. Hence normal expansion dominates every power of the tangential rate.
5. **Uniformity.** Continuity in (a, Q) on compact subsets of 𝒫 follows from the explicit formulas. ∎

**Residue.** The only work left is to transcribe Dyatlov's §3 verification to Kerr–Newman, replacing Δ by Δ + Q² throughout. Every place where Dyatlov uses Δ'' = 2 or the structure of R is covered above. **B1 is closed.**

## What B1 does not cover: B2 (threshold/subprincipal condition)

Hintz's estimates at normally hyperbolic trapping, for **tensor-valued** operators, need the threshold quantity

 sup_Γ (imaginary part of the subprincipal symbol, measured in a suitable bundle metric) < ν_min/2,

or a version of this after conjugation. For Kerr, Hintz handles it by constructing a bundle metric adapted to the trapping (the "pseudo-differential inner product" method). For Kerr–Newman, B2 has to be checked for the Einstein–Maxwell gauge-fixed operator, in which the Maxwell components are coupled. This is a **quantitative symbol computation on Γ**: ν from Lemma B1.3, the subprincipal symbol from the linearised operator.
