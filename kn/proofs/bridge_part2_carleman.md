# Bridge Part 2, Step C: the Carleman lemma at infinity, and a replacement route

## 1. The polynomial-weight Carleman argument does not close

Reduce to the radial form at large r. On each component, after the substitution v = r^{s+1} ψ:

  v'' + (ω² + c/r − Λ/r²) v = W[v],

where:
- c is complex (for spin s, c = 4Mω² + 2isω + O(Q²));
- Λ ≥ 0 is the spin-weighted angular operator;
- W is an operator of order ≤ 1 with coefficients O(1/r).

Conjugate by the weight r^τ, writing w = r^τ v. Then split P_τ = S + A and expand ‖P_τ w‖² = ‖Sw‖² + ‖Aw‖² + ⟨[S,A]w, w⟩.

- The gain is about 4τω² ‖w/r‖². It comes from the term 4τ r^{-2}∂², after using Sw ≈ 0.
- The commutator of A with τ(τ+1)/r² gives a term −4τ²(τ+1) ‖w/r²‖² of the **wrong sign**.

The gain wins only where r ≳ τ/ω. As τ → ∞ that region moves out to infinity, so the usual cutoff-at-fixed-R conclusion fails. The standard way around this has two steps: (i) Froese–Herbst exponential bounds, showing that rapid decay implies e^{αr} decay for every α, then (ii) a Carleman estimate with exponential weights. But step (i) uses self-adjointness (a Mourre or virial argument), which we don't have: c is complex and W is non-self-adjoint and matrix-valued. Known positive-energy results (Kato, Agmon, Froese–Herbst, Koch–Tataru) cover scalar operators, real long-range potentials, and small complex short-range parts. **So Step C is not available off the shelf, and my earlier sketch was wrong.**

## 2. Replacement: go through the horizon instead (flux identity + Fuchsian uniqueness)

Real ω ∉ {0, mΩ_H}. Let the mode have all gauge invariants' radiation fields at I⁺ equal to zero (Bridge Part 1).

**(a) The flux identity kills the horizon data.** Use the canonical energy current of the linearised Einstein–Maxwell equations (Hollands–Wald, Prabhu–Wald), together with conservation for the mode. The T-energy flux through I⁺ is ∫(|∂_u δσ⁰|² + k|δφ₂⁰|²), which is zero. The T-energy flux through H⁺ equals ω(ω − mΩ_H)·(c₁‖δσ_H‖² + c₂‖δφ_{0,H}‖²), with c₁, c₂ > 0. Here δσ_H is the horizon shear and δφ_{0,H} the horizon Maxwell flux, both gauge invariant on H⁺ at linear order. So the horizon data vanish.

**(b) Vanishing horizon data ⇒ vanishing to infinite order.** For horizon-regular modes, the positive-spin invariants (A, F, B, X with X eliminated, §H_X) satisfy a Fuchsian system at r = r₊ (a regular singular point in r). Their leading horizon coefficients are fixed by the horizon data in (a). If those vanish, the Frobenius recursion forces every coefficient to vanish. This needs the indicial roots to be non-resonant, to be checked; the horizon exponent is −s − iσ, with σ ∝ (ω − mΩ_H)/κ ≠ 0.

**(c) Fuchsian uniqueness.** Kerr–Newman is real-analytic, so the reduced system has analytic coefficients. A solution of a Fuchsian system with analytic coefficients that vanishes to infinite order at the Fuchsian boundary vanishes identically nearby. This is Baouendi–Goulaouic, *Cauchy problems with characteristic initial hypersurface* (CPAM 1973), which needs specific conditions on the indicial roots. Elliptic unique continuation (Step D, already done) then spreads the vanishing to the whole exterior.

## Open items for the replacement route

- **(C1)** The canonical-energy flux formula for linearised Einstein–Maxwell modes on Kerr–Newman, with positive coefficients and expressed through gauge-invariant data at H⁺ and I⁺. Cite Hollands–Wald and Prabhu–Wald, and check the Maxwell part.
- **(C2)** The precise relation between the horizon flux data and the leading horizon coefficients of (A, F, B, X), plus the non-resonance of the indicial roots.
- **(C3)** Check the hypotheses of Baouendi–Goulaouic for the mode-reduced Giorgi system near r₊ (Fuchsian structure and indicial-root conditions). Alternatively, prove a direct Carleman estimate near the regular singular point; that is a much easier setting than infinity, because it is local and has a fixed singular structure.
- **(C4)** ω = mΩ_H is excluded here: the horizon flux factor vanishes. That frequency is covered by the threshold-strip and threshold analysis instead.

## C1: canonical-energy flux (status after literature check)

**Sources:**
- Hollands–Wald, CMP 321 (2013): canonical energy E relative to the horizon Killing field K = T + Ω_H Φ.
- Keir, arXiv:1306.6087: extension to theories with p-form gauge fields, which includes Einstein–Maxwell. Property 5 there: the fluxes of E across I⁺ and H⁺ are non-negative "for a certain class of perturbations".
- Sorce–Wald (PRD 96, 104014): canonical energy for Einstein–Maxwell on Kerr–Newman.
- Note: the *manifest* positivity of BOTH fluxes is stated for **axisymmetric** perturbations, or for static backgrounds (Hollands–Wald; see also Holzegel et al., arXiv:1602.04524).

**What we need is weaker, and it holds for m ≠ 0 as well:**
- **Horizon flux.** The K-canonical-energy flux through H⁺ is a positive-definite quadratic form in L_K of the horizon shear and of the horizon Maxwell flux. This is the Hollands–Wald horizon flux; it does not need axisymmetry. For a mode, L_K = −i(ω − mΩ_H), so the horizon flux is (ω − mΩ_H)²·(positive-definite form in δσ_H, δφ_{0,H}).
- **Flux through I⁺.** It is quadratic in the radiation data (δσ⁰ and δφ₂⁰), which vanish by Bridge Part 1. So it is zero, whatever its sign would otherwise be.
- **Balance.** For a real-ω mode, Re h is periodic in t, so E(Σ_t) is periodic. Averaging the flux balance over one period: total horizon flux = −(total I⁺ flux) = 0. Hence δσ_H = δφ_{0,H} = 0 whenever ω ≠ mΩ_H.

**Remaining checks for C1:**
- (i) Hollands–Wald impose gauge conditions near H⁺ and I⁺. The mode must be put in such a gauge with regular gauge vectors.
- (ii) The Maxwell part of the horizon flux must be confirmed positive-definite: use Keir's or Sorce–Wald's form.
- (iii) Bookkeeping for non-axisymmetric perturbations, where E uses K and not T: the I⁺ flux of the K-current involves both energy and angular-momentum flux. Both are quadratic in the radiation data, so both vanish here.

## C2: horizon data ⇒ leading horizon coefficients (sketch)

Write x = r − r₊ and σ = (ω − mΩ_H)/(2κ), where κ is the surface gravity. In ingoing coordinates, a horizon-regular mode has components u_s = x^{c_s + iσ} ψ_s that are smooth up to x = 0. Here c_s is Giorgi's **conformal type**: c = 2 for A and c = 1 for F and B. X is eliminated via §H_X.

**Indicial roots.** For each component they are {0, −(c_s + 2iσ)}, in the variable u_s. The second root has nonzero imaginary part, because σ ≠ 0 when ω ≠ mΩ_H. So no root differs from 0 by a nonnegative integer **within** a component.

**Integer shifts between components.** The exponents differ by c_A − c_{F,B} = 1. So the free data are exactly the leading coefficients of each component at its own exponent. (A coefficient that would otherwise be forced at a resonant order is instead free there, and that free coefficient is precisely the next component's leading datum.)

**Identifying the free data with the horizon flux:**
- A ↔ Ψ₀ on H⁺, which (Teukolsky–Press horizon relation) is ∝ L_K² δσ_H. For a mode with ω ≠ mΩ_H, vanishing δσ_H therefore forces A's leading coefficient to vanish.
- F = −ðφ₀ + … − 2σφ₁ (Remark 5.5). On H⁺ this involves δφ_{0,H} and δσ_H.
- B = 3φ₀Ψ₂ − 2φ₁Ψ₁. On H⁺ it involves δφ_{0,H} and δΨ₁_H. The horizon Bianchi identity gives Ψ₁_H in terms of ð δσ_H and the Maxwell data.

So C1 forces the leading coefficients of A, F and B to vanish. Then every Frobenius coefficient vanishes, since there are no other free data. Any log terms at the resonant orders carry coefficients proportional to data that are now zero.
**To verify:** the three identifications above, with Giorgi's normalisations (§4.1 frame, ingoing regularisation).

## C3: formal ⇒ actual vanishing (Fuchsian uniqueness)

Multiply the mode-reduced Giorgi system by x. In t = −log x it becomes ∂_t² + (const)∂_t + e^{−t}(angular elliptic) + …. Equivalently, it is a **Fuchsian** system in x, in the sense of Baouendi–Goulaouic, with real-analytic coefficients (Kerr–Newman is analytic).

- **Scalar case:** Baouendi–Goulaouic, CPAM 26 (1973). A C^∞ solution of a Fuchsian equation whose characteristic exponents are not nonnegative integers (beyond the prescribed ones), and whose prescribed traces vanish, is identically zero near x = 0.
- **Systems:** Tahara's theory of Fuchsian systems (e.g. Tahara, Japan. J. Math. 1979).

The characteristic exponent −(c_s + 2iσ) is non-real, so the non-integer condition holds. The integer offsets between components are the resonance bookkeeping of C2.

⇒ The mode vanishes near H⁺. Step D (elliptic unique continuation, already done) then gives vanishing on the whole exterior.

**To verify:** (i) the precise Fuchsian normal form for the Giorgi system after eliminating X, including how the coupling enters at each order in x; (ii) that the hypotheses of the Tahara/Baouendi–Goulaouic theorem match; (iii) at ω = mΩ_H we have σ = 0 and the exponents collide, consistent with excluding that frequency (C4).

## C3 refined: the Fuchsian framework, and the one condition left to check

**Setting.** The modes are separated only in t and φ. In Kerr–Newman the θ-dependence does not separate for the coupled system, so the problem is a PDE in (x, θ) and not an ODE. That is exactly the setting of Baouendi–Goulaouic (BG): the variable t of their papers is our x, and their y is our θ.

**Weight check, single component.** For a single component, the radial part is Δ^{-s}∂_r(Δ^{s+1}∂_r) with Δ ≈ 2κ' x, where κ' = r₊ − r₋ > 0 (subextremal). Multiplying by x gives
 x²∂_x² + x·(a₀ + iσ-terms)∂_x + b₀ + x·L_θ + O(x²),
where L_θ is the angular elliptic operator. In BG's normal form (order m = 2, weight k = 2), a term ∂_θ^α ∂_x^j with j < 2 must carry a coefficient that vanishes at x = 0. Here the angular term has the factor x, so the condition holds. The indicial polynomial is ρ(ρ−1) + a₀ρ + b₀. Its roots are {0, −(c_s + 2iσ)} after the shift u_s = x^{c_s+iσ}ψ_s (C2).

**BG uniqueness (CPAM 1973, Thm. 3).** The coefficients must be real-analytic in θ and C^∞ in x; Kerr–Newman supplies this. Suppose the characteristic exponents other than the prescribed ones are not in ℕ. Then a C^∞ solution whose traces at the prescribed exponents vanish is identically zero near x = 0. The free trace u_s(0, θ) is a function of θ, and C1/C2 kill it. The exponent −(c_s + 2iσ) is not real, so it never lies in ℕ. ⇒ **Single component: done**, modulo citing the exact statement.

**The coupled system: what remains.** Giorgi's system couples A (c = 2) to F and B (c = 1) through first-order terms. After the shifts u_s, a coupling term P_{AF} u_F appears in the A-equation with a factor x^{c_F − c_A} = x^{−1}.
- If the coupling coefficient vanishes at H⁺ to order ≥ 1, the system is Fuchsian as it stands (Tahara 1979 applies), and C3 closes.
- If it doesn't, rescale to ũ_F = x·u_F. That restores the Fuchsian form, but shifts F's regular exponent to 1. The trace of ũ_F at exponent 1 is then the free datum, and it is again ∝ δφ_{0,H} (C2). The non-integer condition is unaffected.

Either way the argument closes **provided** one condition holds:

 **(F) Triangularity.** The leading-order (x⁰) part of the rescaled coupling matrix must be block-triangular in (A | F, B). This prevents an exponent from one block landing at an integer offset from the other block's prescribed exponent with a nonzero resonance coefficient. If it fails, log terms appear. They still carry vanishing coefficients (C2), so uniqueness survives, but Tahara's theorem would then have to be applied in its resonant form.

So (F) is the single remaining check for C3. It requires the explicit near-H⁺ expansion of Giorgi's coupling coefficients, in particular eqs. (5.x)–(6.x) of arXiv:2002.07228 in the ingoing frame. That is a symbolic computation: about 1 script, checkable with sympy once the coefficients are transcribed exactly. I am not writing those coefficients from memory, since I made that mistake before.

**Status of the bridge after this:** C2 identifications (normalisation check), C3-(F), C4 (ω = mΩ_H), R1–R3. Everything else in Part 2 is reduced to cited theorems.

## C3-(F): triangularity verified from Giorgi's Theorem 6.1

Source: arXiv:2002.07228v2, Thm. 6.1, eqs. (79)–(87), together with the transport equation (78) for X. Coupling terms:
- M₁[F,X] (B-eq) = 2|P_F|²( 2 D·F + 4 H̄·F − (2tr X̄ − tr X) X )
- M₂[A,X,B] (F-eq) = −P_F( ∇₃A + ½(3tr X̄ − tr X)A ) + (3/2 ∇₃H)⊗X + (2H − H̄)⊗B
- M₃[F,X] (A-eq) = 2P_F( 2∇₄F + 2tr X̄ F + (H + H̄)⊗X )
- (78): ∇₃X + ½tr X̄ X = −D̄·F − H̄·F − 2B

**Near-H⁺ scaling of the operators (per mode e^{−iωt+imφ}).** Δ ≈ (r₊ − r₋)x. Of the two principal null derivatives, exactly one is horizon-singular: it is x^{−1}·(Euler operator x∂_x + const). The other is regular: Euler operator plus O(x). Which is which depends on Giorgi's normalisation of e₃/e₄ (§4.1). The angular operators D·, D⊗ and the zeroth-order coefficients are O(1).

- Each Teukolsky operator T_s contains ∇₃∇₄ (or ∇₄∇₃), which is x^{−1}·(Euler)², plus angular O(1) terms. After multiplying by x, its principal part is (Euler)² + x·L_θ. This is Fuchsian, as checked in the previous section.
- After multiplying by x, every coupling term that is angular (D·F, D⊗X), zeroth-order (tr X F, H⊗B, P A, …), or contains the *regular* null derivative acquires a factor x. So it vanishes at x = 0.
- The only coupling that survives at leading order is the one that contains the *singular* null derivative. Among the couplings between A and F that is exactly **one** direction: either −P_F∇₃A in the F-equation, or 2P_F·2∇₄F in the A-equation, but not both, since ∇₃ and ∇₄ cannot both be singular.
- X: multiplying (78) by x gives Euler(X) = x·(…F, B) if ∇₃ is the singular derivative. In that case its couplings vanish at x = 0. If ∇₃ is the regular one, (78) is an Euler ODE for X with O(1) source terms in F and B, and multiplying by x shows nothing. Then the alternate sequencing applies: the source enters through B and F. Both sequencings are triangular, as the argument below shows.

**Conclusion.** The leading-order (x⁰) indicial matrix of the full system (A, F, B, X), written in Boyer–Lindquist (BL) components, is **block-triangular**. Its single nonzero off-diagonal block couples A and F through one singular null derivative. Hence

 det L(ρ) = I_A(ρ) · I_F(ρ) · I_B(ρ) · I_X(ρ),

and the set of indicial roots of the system is the union of the roots of the individual components. At each order n, the Frobenius recursion has the form L(ρ+n)·u_n = (known terms from orders < n). Since L(ρ+n) is triangular, it is invertible except where some component's indicial polynomial vanishes. At those points the free coefficient is exactly that component's trace, as in C2. So with all traces zero, every coefficient vanishes, including any potential log coefficient, because log terms are forced only by nonzero lower-order data. ⇒ **formal vanishing holds, with no extra condition.**

**Actual vanishing.** We work directly with the BL components, which are Fuchsian (horizon-singular exponents simply appear as negative or complex roots). Tahara's uniqueness theorem for Fuchsian systems with real-analytic coefficients then applies, as in the single-equation case.

**Robustness.** The triangularity uses only two structural facts: (a) in Giorgi's system, A and F talk to each other through ∇₄ in one direction and ∇₃ in the other; (b) all other couplings are angular or zeroth-order. It does not depend on the precise coefficient functions. For the same reason the earlier transcription errors cannot affect it.

**Remaining residue in C3:**
- (i) confirm which of e₃, e₄ is horizon-singular in Giorgi's §4.1 normalisation. This is needed only for bookkeeping; both cases are triangular.
- (ii) quote the precise statement of Tahara's theorem: weight conditions and the regularity class of the solution.

Both are citation and bookkeeping items. **C3-(F) is closed.**

## C4: ω = mΩ_H (sketch)

Here σ = 0, and the two indicial roots of each component collide at {0, −c_s}. The second exponent is then a negative integer. A mode that is smooth at H⁺ still has only the regular trace as free data. The collision produces log terms only in the singular branch, which smoothness excludes. So C2 and C3 go through formally. The real issue lies elsewhere: at σ = 0 the horizon flux in C1 degenerates, because the flux is weighted by ω(ω − mΩ_H). So C1 does not yield vanishing horizon data at this frequency. This is the standard superradiant-threshold frequency. In Kerr it is handled by a separate argument (Whiting / Andersson–Ma–Paganini–Whiting: a real mode at ω = mΩ_H has vanishing flux at null infinity, and then the I⁺ route applies). **Open: run the I⁺ route (Part 1, R1–R3) at this one frequency.** This ties C4 to R1–R3, so they should be done together.
