# A7: from real-axis mode stability to the closed upper half-plane (the continuity argument)

**Reading of A7.** "Threshold strip" means the frequencies with Im σ ∈ (0, δ) near the superradiant band 0 ≤ Re σ/m ≤ Ω_H. No energy identity controls them, and in Kerr they are handled through separability (Whiting / Shlapentokh-Rothman). The clean way to dispose of them, and in fact of the **whole** open upper half-plane at once, is the deformation argument. It needs no estimate inside the strip at all. Below it is written as a proof, with its inputs isolated.

## Setting

- Parameter space 𝒫 = {(M, a, Q): M > 0, a² + Q² < M²}. Normalise M = 1, so 𝒫 ≅ the open unit disc in (a, Q). It is connected and contains the Reissner–Nordström (RN) segment {a = 0}.
- For p ∈ 𝒫, let L_p(σ) be the spectral family of Hintz's constraint-damped, gauge-fixed linearised Einstein–Maxwell operator on Kerr–Newman (KN), acting on Vasy-type variable-order (or b-/scattering-) Sobolev spaces 𝒳^s_p → 𝒴^{s−1}_p.
- A **mode** (resonance) is a σ with Im σ ≥ 0, σ ≠ 0, and ker L_p(σ) ≠ 0.
- **Goal (MS):** for every p ∈ 𝒫 there is no mode. This is exactly AHW24 Thm 6.1(1) in Hintz's black-box form.

## Inputs (each is a separate item of the programme)

- **(I1) Fredholm analytic family, continuous in p.** L_p(σ) is Fredholm of index 0 for Im σ > −c(p). It is analytic in σ and continuous in p in operator norm on a fixed space; one can identify the spaces via a p-dependent diffeomorphism fixing the horizon at r = r₊(p), the standard device in Hintz's papers. Fredholmness needs radial-point estimates at the horizon and at infinity, and normally hyperbolic trapping, **uniformly** on compact subsets of 𝒫. ← **B1** (r-normal hyperbolicity of KN trapping), **B2**.
- **(I2) High frequency, uniform.** For every compact K ⊂ 𝒫 there are R_K and C_K with ‖L_p(σ)⁻¹‖ ≤ C_K for p ∈ K, Im σ ≥ 0, |σ| ≥ R_K. ← semiclassical estimates: trapping (B1), together with Lemma A, which separates superradiance from trapping and so prevents a superradiant obstruction at high frequency; plus **A6** for |m| ≫ |σ|.
- **(I3) No real modes.** For every p ∈ 𝒫 and every σ ∈ ℝ∖{0}, ker L_p(σ) = 0. ← **A5** (real-axis reduction), via the bridge. The bridge reduces this to gauge-invariant real modes; C4 / route (α) closes the I⁺ side.
- **(I4) Zero frequency, uniform.** For every compact K there is ε_K > 0 such that L_p(σ) is invertible for p ∈ K and 0 < |σ| ≤ ε_K, Im σ ≥ 0. ← **A11** (zero-energy structure, including the δQ, δP charge modes, with a uniform low-energy resolvent expansion in the style of Hintz Ch. on σ → 0).
- **(I5) Base case.** MS holds for RN, i.e. p = (1, 0, Q) for every |Q| < 1. This is classical: coupled grav-EM mode stability of RN, via the decoupled Moncrief/Chandrasekhar potentials, which are positive (Giorgi's RN papers [17, 18] in her numbering). It must be verified that this is stated for Hintz's gauge-fixed operator, not only for the gauge invariants. Going from one to the other is the bridge again, with a = 0, where it is classical.

## Proof of (MS) from (I1)–(I5)

Fix p₁ ∈ 𝒫. Join it to p₀ = (1, 0, 0) by a path γ: [0,1] → 𝒫, for instance a straight segment in the (a, Q)-disc; its image is compact. Let K = γ([0,1]), and take R = R_K and ε = ε_K from (I2) and (I4). Consider the half-annulus

 U = {σ: Im σ > −η, ε < |σ| < R},

where η > 0 is still to be chosen.

*Step 1: no modes on ∂U ∩ {Im σ ≥ 0}, uniformly in t.*
- On |σ| = R: by (I2).
- On |σ| = ε: by (I4).
- On the real segments ε ≤ |σ| ≤ R: by (I3).

*Step 2: open up to Im σ > −η.* The set {(t, σ): t ∈ [0,1], σ ∈ ∂U, Im σ ≥ 0} is compact, and on it L is invertible. Invertibility is an open condition (I1), so there is η > 0 for which L_{γ(t)}(σ) is invertible on the whole of ∂U, i.e. also on the short arcs with Im σ ∈ (−η, 0], for every t. Shrink η below inf c(γ(t)), so that the Fredholm property holds on all of Ū.

*Step 3: the count is constant.* For an analytic Fredholm family that is invertible on the boundary of a bounded domain, the Gohberg–Sigal (logarithmic residue) count

 N(t) = (1/2πi) tr ∮_{∂U} L_{γ(t)}(σ)⁻¹ ∂_σ L_{γ(t)}(σ) dσ

is a non-negative integer. It equals the number of resonances in U, counted with multiplicity, and it is continuous in t, because L and L⁻¹ are norm-continuous on ∂U, uniformly. So N(t) is constant.

*Step 4: modes in the upper half-plane can enter U only through ∂U.* Let N₊(t) be the number of resonances in U₊ = U ∩ {Im σ > 0}. Its boundary consists of the arcs (Step 1) and the real segments (I3), and none of them carries a resonance for any t. The same argument as Step 3, applied to U₊ instead of U, gives **N₊ constant in t**. (Step 2 is not needed for U₊ itself; it only fixed the Fredholm strip.)

*Step 5: base case.* N₊(0) = 0 by (I5), and the modes with |σ| ≥ R or |σ| ≤ ε are excluded by (I2) and (I4). Hence N₊(1) = 0. Together with (I2), (I3) and (I4), p₁ has no modes with Im σ ≥ 0 and σ ≠ 0. ∎

## What this shows about the threshold strip

The strip needs **no direct analysis.** A growing mode near the superradiant band can only be created by a resonance crossing the real axis, which (I3) forbids, or arriving from infinity or from zero, which (I2) and (I4) forbid. This is the KN analogue of Shlapentokh-Rothman's continuity argument (Kerr, 2015). The only differences are that the base case is RN instead of Schwarzschild, and that Gohberg–Sigal replaces his ODE argument, since KN does not separate.

## Honest assessment: where the difficulty went

A7 is now **soft**. The difficulty has moved entirely into:
- **(I3) = A5: no real modes in the superradiant band 0 < σ/m < Ω_H.** The I⁺ route (α) handles σ = mΩ_H, where the energy flux forces zero radiation. *Inside* the band, the flux identity
 ω²·|rad at I⁺|² + ω(ω − mΩ_H)·|rad at H⁺|² = 0
 does **not** force zero radiation, because the second coefficient is negative. In Kerr this is exactly where Whiting's transformation is used. **This is the core open problem of KN mode stability.** It is where the "new idea" has to live, and the earlier near-extremal and balance-criterion numerics were all probing it.
- (I2): uniform high-frequency estimates. These rest on B1 plus Lemma A; Lemma A is proved, B1 is open but expected.
- (I4) = A11: the zero-energy analysis.

## Residue for A7 itself

- (a) Hintz's spaces must be identified across p uniformly. This is standard and is done in Hintz's Kerr paper across a ∈ (−M, M); the extension to (a, Q) is mechanical.
- (b) The Gohberg–Sigal count for Fredholm (not compact-perturbation) families: the analytic Fredholm theorem plus the logarithmic-residue formula (Gohberg–Sigal 1971; see also Dyatlov–Zworski, Appendix C).
- (c) (I5) in gauge-fixed form, as noted above.
