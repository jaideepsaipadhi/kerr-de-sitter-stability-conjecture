# Bridge Part 1: R1, R3 closed; R2 reduced; C4 re-assessed

## R3: the coefficient of δφ₂⁰ (closed)

Peeling: Φ₁ ~ r⁻², Ψ₂ ~ r⁻³, Ψ₃ ~ r⁻², Φ₂ ~ r⁻¹. Both terms of ψ₋₁ ∝ 2Φ₁Ψ₃ − 3Ψ₂Φ₂ are O(r⁻⁴).

Linearise about the background. In the principal frame the background has Ψ₃ = Φ₂ = 0, so the r⁻⁴ coefficient is

 2φ⁰₁,bg δψ₃⁰ − 3ψ⁰₂,bg δφ₂⁰.

For Kerr–Newman, Ψ₂ = −M/(r − ia cosθ)³ + Q²/((r − ia cosθ)³(r + ia cosθ)), so ψ⁰₂,bg = −M. Once δψ₃⁰ = 0 (step 1), the coefficient of δφ₂⁰ is **3M ≠ 0**.

The DGS normalisation and any frame changes only contribute nonzero overall factors, by the boost-weight remark in R1. ∎

## R1: Kinnersley → Bondi frame change (closed, under hypothesis (c))

Write the background change of frame as boost/spin λe^{iϑ} ∘ null rotation about l (parameter c) ∘ null rotation about n (parameter b). Only **boundedness** of c and b is needed; no decay rate. In both frames l is asymptotically the outgoing null direction and the spin phases converge, so λ → const ≠ 0 and ϑ → const, and c, b = O(1). This holds for any Bondi frame built on the stationary good cuts.

**Boost/spin.** Ψ₄ has boost weight −2. Both terms of ψ₋₁ have boost weight −1 (Ψ₃ and Φ₂ have weight −1; Φ₁ and Ψ₂ have weight 0). So each quantity rescales by a single nonzero factor. Perturbing λ multiplies the background values Ψ₄,bg = 0 and ψ₋₁,bg = 0, so it contributes nothing.

**Null rotation about l** (n → n + c̄m + cm̄ + |c|²l):
- Ψ₃' = Ψ₃ + 3cΨ₂ + O(c²Ψ₁), Φ₂' = Φ₂ + 2cΦ₁ + c²Φ₀, Φ₁' = Φ₁ + cΦ₀, Ψ₂' = Ψ₂ + 2cΨ₁ + c²Ψ₀.
- The terms in ψ₋₁ linear in c cancel **identically**: 2Φ₁·3cΨ₂ − 3Ψ₂·2cΦ₁ = 0. This is the algebraic reason ψ₋₁ is gauge invariant. So the δc contributions are exactly zero.
- The terms quadratic in c, once linearised, contain δΨ₁, δΦ₀ or δΨ₀ (the background values vanish), multiplied by Φ₁,bg or Ψ₂,bg. They are therefore O(r⁻⁶) or smaller, which is ≪ r⁻⁴.
- Ψ₄' = Ψ₄ + 4cΨ₃ + …. With c = O(1) and δΨ₃ = O(r⁻²) (peeling, hypothesis (b)), the r⁻¹ coefficient is unchanged.

**Null rotation about n** (l → l + b̄m + bm̄ + |b|²n):
- Ψ₄ and Φ₂ are invariant. Ψ₃' = Ψ₃ + b̄Ψ₄, Φ₁' = Φ₁ + b̄Φ₂, Ψ₂' = Ψ₂ + 2b̄Ψ₃ + b̄²Ψ₄.
- Linearised, the only new term in ψ₋₁ is 2Φ₁,bg b̄ δΨ₄. It vanishes **identically**, because hypothesis (c) is δΨ₄ ≡ 0 everywhere, not just at leading order.
- Every other term multiplies a background quantity that vanishes (Φ₂,bg, Ψ₃,bg, Ψ₄,bg).

So the leading coefficients of δΨ₄ and δψ₋₁ in the Bondi frame are nonzero constant multiples of the Kinnersley ones. ∎

## R2: Bondi gauge for outgoing modes (reduced to citations)

- Outgoing modes of Hintz's gauge-fixed operator are, by the scattering-calculus resolvent mapping properties of Hintz's framework, of the form e^{iωr*} × (polyhomogeneous in 1/r).
- A linearised gauge transformation into Bondi gauge is found by solving transport ODEs along the outgoing null cones. Polyhomogeneous inputs give polyhomogeneous outputs (Chruściel–MacCallum–Singleton, *polyhomogeneous Bondi*).
- The proposition only uses the leading terms of δΨ₄ (r⁻¹), δΨ₃ (r⁻²) and δΦ₂ (r⁻¹), plus the asymptotic equations (2.52)–(2.56), which hold at leading order in the polyhomogeneous setting. Logarithms first appear in Ψ₀/Ψ₁ at subleading orders.
- Since Part 2 now runs through the horizon, the positive-spin I⁺ steps 4–7 are **no longer needed**, except for C4 (below).

**Residue:** cite the precise resolvent statement (Hintz arXiv:2606.28253 §§ on outgoing resolvent mapping) and the CMS theorem. This is a citation item, not new mathematics.

## C4 re-assessed: an honest risk (the main open item in the bridge)

At ω = mΩ_H (m ≠ 0, ω real):

1. **Radiation at I⁺ vanishes for free.** Take the T-energy flux identity for a real-frequency mode. The flux through H⁺ is weighted by ω(ω − mΩ_H) = 0, so the flux through I⁺ vanishes, i.e. δσ⁰ = δφ₂⁰ = 0. (The canonical-energy version of C1 is used with the roles of H⁺ and I⁺ swapped.)
2. **But neither end closes the argument.**
   - *Horizon side:* at σ = 0 the horizon flux carries no information, so the horizon traces (C2) are not forced to vanish.
   - *Infinity side:* we would need Rellich-type unique continuation from I⁺ for a **non-separable, non-self-adjoint system with long-range first-order terms**. That is exactly the route that failed earlier (bridge_part2_carleman.md).
   - In Kerr this frequency is handled using separability: Shlapentokh-Rothman's argument is ODE-based. Kerr–Newman has no separability and no Teukolsky–Starobinsky identities.
3. **Structural handle.** At ω = mΩ_H the mode satisfies L_K δ = 0, where K = ∂_t + Ω_H∂_φ is the horizon generator. So we need a *linearised rigidity* statement: a smooth, K-invariant, horizon-regular, bounded linearised Einstein–Maxwell solution with no radiation at I⁺ is pure gauge plus a Kerr–Newman parameter change. In the region where K is timelike (inside the velocity-of-light surface) the K-reduced equations are **elliptic**. That is promising, and it is the same structure exploited in rigidity proofs (Alexakis–Ionescu–Klainerman; Wong–Yu for Kerr–Newman). Beyond that surface the reduced problem is hyperbolic. There the zero radiation from step 1 is the input, but some form of unique continuation from infinity is still needed.

**Candidate routes for C4** (to try next, in this order):
- (α) Mode-level unique continuation from I⁺ for **real** ω only. Avoid the polynomial Carleman, which failed because of complex ω, and use a Kato/Hörmander Rellich theorem for systems after conjugating away the long-range first-order terms with e^{iω r* log}-type phases. The first-order coefficients are O(1/r) and **real-ω specific**, which is where the earlier obstruction may not apply. (This needs checking: the earlier failure was the wrong-sign τ³/r⁴ term.)
- (β) A K-invariant ellipticity / linearised-rigidity argument in the ergoregion-complement, modelled on Wong–Yu.
- (γ) Continuity in parameters: exclude real modes at ω = mΩ_H via a Kerr deformation argument. This is weak, because it does not give the full subextremal range.

**Bridge status:** C1–C3 (with citation residue), R1 and R3 are closed. R2 is reduced to citations. **C4 is open and now the largest risk in the bridge.**

---
## C4, route (α): Rellich uniqueness from I⁺ for real ω, run on Giorgi's Regge–Wheeler system

**Why the earlier Carleman failed, and why this is different.** The earlier attempt worked with the Teukolsky-type variables. There the first-order term ~ (2s/r)∂_r is **long-range**: r·|W| does not tend to 0. On top of that, complex ω produced the wrong-sign τ³/r⁴ term. Every Rellich/Carleman theorem for non-self-adjoint operators excludes long-range complex first-order terms. Here instead we use Giorgi's **generalized Regge–Wheeler system** (arXiv:2002.07228, Thm 7.3, eqs. (94)–(95)) for p and q_F, at **real** ω. For a mode e^{−iωt+imφ}, at large r:

| term | size | class |
|---|---|---|
| □̇₁, □̇₂ principal part | metric δ + O(M/r) | long-range **metric**, real, smooth |
| ω²g^{tt} part | ω²(1 + 4M/r + …) | long-range **real** potential (Coulomb-type) |
| −i(2a cosθ/\|q\|²)∇_t p | 2aω cosθ / r² | **short-range** complex, zeroth order |
| V₁, V₂ | real | long-range real potentials |
| 4Q² q̄³/\|q\|⁵ D·q_F, ½ q³/\|q\|⁵ D⊗p | r⁻² · (r⁻¹∂_θ) | **short-range** first-order coupling |
| L_p[B,F], L_{q_F}[B,F] | coefficients O(\|a\|/r²), O(\|a\|/r³) | short-range, **but nonlocal**: involves B, F, X, not p, q_F |
| a m ω cross terms | O(r⁻³) | short-range |

So, apart from the L-terms, the system is exactly in the class of **Koch–Tataru** (CMP 2006, *Carleman estimates and absence of embedded eigenvalues*). That class is −Δ_g + W·∇ + V, with g a long-range perturbation of Euclidean space, V = V_L (real, long-range) + V_S (complex, short-range), and W complex short-range. It is vector-valued with smooth coupling. The Rellich-type conclusion: at energy λ = ω² > 0, if liminf_{R→∞} R⁻¹∫_{R<|x|<2R}|u|² = 0, then u ≡ 0 outside a compact set.

**The decay hypothesis.** We need p and q_F = o(r⁻¹) on average. Part 1's proposition gives δσ⁰ = δψ_n⁰ = δφ_j⁰ = 0 for all leading peeling coefficients. This uses hypothesis (c), which in the bridge is supplied by the negative-spin vanishing. B, F and X are polynomial combinations of these with background coefficients, so all of them gain one order of decay. p = f₁(∇₃B + C₁B) and q_F = f₂(∇₃F + C₂F) are built so that RW solutions behave like r⁻¹ e^{iωr*}. With the leading coefficients of B and F gone, p and q_F are O(r⁻² log r) by polyhomogeneity (R2). This gives the averaged decay condition.
**To verify (D1):** the exact r-weights of B and F under peeling, and that the r⁻¹ coefficients p⁰ and q_F⁰ are linear in the vanished coefficients. This is a mechanical peeling bookkeeping check.

**The crux (D2): the nonlocal L-terms.** B, F and X are recovered from p, q_F through transport equations along e₃ (ingoing), from (90)–(91) and (78). At large r, for a mode, ∇₃ ≈ −(∂_r + iω) + O(1/r). Since B and F decay (D1), the transport equations integrate from infinity:
 B(r) = −∫_r^∞ e^{iω(s−r)} (p/f₁)(s) ds·(1 + O(1/s)),
and similarly for F and X. Two ways to handle this:
- (D2-a) **Enlarge the system.** Take u = (p, q_F, B, F, X) and run the Carleman estimate for a mixed second-order/transport system. With a radial weight e^{τφ(r)}, φ′ > 0, the transport part ∂_r + iω + O(1/r) has a Carleman estimate with a **gain** of τ^{1/2} in the weighted norm. The coupling back into the RW equations comes with coefficients O(|a|/r²), and the Koch–Tataru estimate controls τ^{1/2}·r^{-1/2}-weighted quantities. So the L-terms are absorbed for large r and τ, provided the coefficient decay r⁻² beats the r⁻¹ that short-range absorption requires. It does, by one power. This is the standard way of adding ODE components to a Carleman estimate (cf. Hörmander Vol. III 17.2 with first-order systems).
- (D2-b) **Eliminate B, F, X.** The tail integral is a Volterra operator. Under the zero-radiation decay, ‖r⁻² ∫_r^∞ |g|‖ ≲ ‖r^{-1−ε} g‖ in the weighted spaces, so the L-terms become a nonlocal short-range perturbation. This needs Koch–Tataru's estimate to tolerate a nonlocal but Volterra-type short-range term, which is plausible because their argument is a weighted L² estimate with radial weights.

I expect (D2-a) to be cleaner. In both versions, the point is that **every coefficient that couples to B, F or X is O(r⁻²)**. That follows from Giorgi's explicit fall-offs (Prop. 7.6: Y, W = O(|a|/r²) or better) and from the multiplicative factors f₁ and f₂.

**Payoff.** If (D1) and (D2) hold, route (α) gives unique continuation from I⁺ for **every real ω ≠ 0**, not only at ω = mΩ_H. So for real frequencies the I⁺ route and the horizon route become redundant, which is a safety margin. Complex ω (Im ω > 0) still uses the horizon route (C1–C3), where the flux identity has no degeneracy.

**Then, after vanishing near infinity:** unique continuation into the interior uses the same elliptic/Fuchsian step as before (Step D). It applies because at the mode level, away from the horizon and the ergosurface, the system is elliptic in (r, θ) for |ω| small relative to …. Correction: it is not elliptic in general. The mode-reduced wave operator is of real principal type in (r, θ) where the K-reduction is hyperbolic. Interior propagation of zero must therefore use either analyticity (Holmgren/John: the coefficients are real-analytic, and Holmgren's theorem needs only non-characteristic surfaces) or the original Step D argument. **Holmgren suffices.** The surfaces {r = const} are non-characteristic for the mode-reduced operator wherever Δ ≠ 0, since the coefficient of ∂_r² is Δ/|q|² > 0. Zero data on r = R₀ ⇒ zero on {r ≥ r₊ + ε}, sweeping inward through r = const, by analytic coefficients + Holmgren applied to the PDE in (r, θ) for each fixed mode.

**C4 status:** reduced to (D1), a bookkeeping check, and (D2), a Carleman absorption of nonlocal O(r⁻²)-coupled terms. This is a definite, standard-looking technical step, not a missing idea. The main risk is that Koch–Tataru's precise hypotheses on the long-range part (derivative bounds on g and V_L) must hold for the tensor-valued RW system in Giorgi's frame. The background quantities are polynomial in 1/r, so this should be routine.
