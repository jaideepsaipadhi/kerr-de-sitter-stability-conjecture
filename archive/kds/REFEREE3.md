# Third hostile referee pass, end to end, on the final manuscript (8 Nov 2026, night) — in-house, not independent

Object: writeup/kds_mode_stability.tex (17 pp.), checked line by line against LEMMAS.md, the certificate code, and the
cited papers (HPV 2508.06620, Hintz 2608.25868, read directly).

## Errors found in this pass and fixed in the manuscript
1. **Kernel lemma was false as stated** (draft of item 124). For Im σ ≠ 0, a mode annihilated by (Þ')⁴ CAN be regular at
   both horizons when σ lies on a discrete "ladder" {mΩ₂ + iNκ₂♯}: there the two Frobenius exponents at r₂ differ by a
   real number in the weight lattice. Fixed: the lemma now excludes the ladder Z_m. Compactness is proved sector by
   sector in m. On a fixed sector, ladder points in the strip 0 < Im σ ≤ C₀ (energy bound) satisfy |σ| ≤ |m|Ω₂ + C₀,
   so they cannot carry a QNM to infinity.
2. **SdS base case at ladder points.** Handled by choosing a generic SdS base point. On a separated solution, (Þ')⁴ is
   the radial Teukolsky–Starobinsky operator, so it kills a solution only at zeros of 𝒞_ℓ(σ). These zeros avoid
   iNκ₂♯ off a countable union of nowhere dense sets of ΛM².
3. **Superradiant threshold ω = mΩ_e, s = +2.** Here ρ₁ = −2 and the termwise Riemann–Liouville factor has Gamma poles,
   so the old "injective on the closed window" claim did not cover this point (it lies outside the window but is used
   by the Wronskian argument). Fixed: Prop. 3.4 now proves ỹ ≠ 0 with the right branches. The singular part is handled
   by the Pochhammer transform and the analytic part by the convergent integral.
4. **Theorem B (1) for T-modes that are not ∂φ-eigenfunctions.** The per-(σ, m) gauge vectors do not obviously sum to a
   smooth field. Fixed: a single local formula X = Y_T/(iσ) with T = α_T ξ + β_T ζ, and J[h] = 0 from L_T J = −iσJ.
5. **Theorem B (2), the part with m ≠ 0 at T-frequency 0.** Fixed by the mean-zero φ-antiderivative of −Y_φ.
6. σ = 0, m = 0 identities re-verified (new script step125_sigma0_check.py, exact): the TS identity and the
   determinant [λ² + 6(1−α)λ + 8 − 4α + 8α²]². Added the needed λ > 0 and the sign of both roots. Fixed the wording
   ("4 boundary conditions" → coefficient system).
7. Sign of the dual exponent at z = 1 (−e₁, not e₁) corrected.

## Load-bearing inputs from the literature (correctly cited, not re-proved here)
L1. Fredholm theory for the Teukolsky family on the full subextremal range (Vasy; Petersen–Vasy; the remark in HPV §2.3
    that subprincipal terms only shift thresholds). Includes continuity in b and the energy bound Im σ ≤ C₀.
L2. HPV's high-energy estimate for the gauge-fixed operator, uniform for b near b₀ (HPV Rem. 2.12, Thm 2.13;
    HV18 §5.1.2).
L3. Wald's adjoint method with Λ, and the ORG formula Ψ₄ = ½(Þ')⁴Φ̄. The adjoint identity O'† ∘ P = P ∘ O,
    P = Ψ₂^{−4/3}, is now proved in the paper from the Λ-free Bianchi identities.
L4. Euler transform of Heun equations (Kazakov–Slavyanov, Takemura).
L5. Teukolsky–Starobinsky identities on SdS, and the form of 𝒞_ℓ (base case only).
L6. HV18 §7 (UEMS on SdS).

## Points an expert should look at first (no known error; least redundant)
P1. The base-case genericity argument (L5 plus analyticity in ΛM²).
P2. Prop. 7.9, Step 5: smoothness across the horizons of the σ = m = 0 gauge vector. This is argued via regular adapted
    coordinates, not computed.
P3. Prop. 3.5 (the real dual): the computer-algebra identity C1 plus the exponent bookkeeping.
P4. The angular two-level bound (interval certificate C2).

## Verdict
No step is known to be false after the fixes above. Every gap found in this pass was closed in the text. The argument
depends on L1–L6, which are standard but have not been re-derived. VERDICT: PASSES the in-house end-to-end referee
with these dependencies stated. It has not been checked by an independent human expert.

## Addendum (gaps closed for an unconditional statement)
* F1–F3 for Teukolsky: Lemma 6.1 (regular frame) makes Teukolsky a stationary principally scalar operator with
  smooth coefficients on a line bundle over M. HV18 Thm 5.4, whose Fredholm part does not use trapping, together with
  §5.1.2 (uniformity) and Prop. 5.11 (continuity of QNMs) then apply verbatim. The dynamics hypotheses come from
  PV25/PV24.
* HPV uniformity: HV18 §5.1.2 states that the hypotheses, and hence the high-energy estimates (5.7)/(5.8), hold
  uniformly for w near w0. HPV Rem. 2.12 and Thm 3.1 verify the hypotheses for L_b.
* Kernel lemma: the weight exponents are integers, which gives the ladder {mΩ₂ + iNκ₂♯ : N ≥ 1}.
* SdS base case: explicit and exact. (Þ')⁴(r⁴R) = Δ²(D†)⁴Δ²R/(16r⁴) and D⁴Δ²(D†)⁴Δ²R = [(λ+4)(λ+6)]² + 144M²σ² (times R);
  Λ drops out (step126/127). The failure set is countable in ΛM².
* Collars (HPV's M extends past the horizons): Lemma 7.1 proves vanishing of δΨ₀,₄ in the collars by the
  domain-of-dependence energy argument. The σ = m = 0 gauge vector is glued across the horizons by a pointwise
  injectivity argument, giving a smooth X on all of M.
* ORG formula: Λ-free, since only Þ'-derivatives and Þ'ρ' = ρ'² enter and no commutators.
