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
