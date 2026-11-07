# B2: the threshold (subprincipal) condition at the trapped set: reduction

**What is needed.** Hintz's estimates at normally hyperbolic trapping for bundle-valued operators (Hintz, *Normally hyperbolic trapping and tensor-valued waves*, 2017) need the following. On the trapped set Γ, the subprincipal operator S_sub (a first-order operator along the Hamiltonian flow H_p, acting on the fibre) must satisfy

 sup_Γ σ(S_sub; bundle metric B) < ν_min/2,

where σ(·; B) is the growth rate of the S_sub-transport in a suitable positive bundle metric B, which may be chosen pseudo-differentially, and ν_min is the minimal normal expansion rate. B1.3 gives ν explicitly: ν = Σ⁻¹(R''(r_c)/2)^{1/2}.

For the gauge-fixed linearised Einstein–Maxwell operator on u = (h, Ȧ):

 S_sub = ∇^{π*E}_{H_p} (the connection part) + C(ξ) (the h ↔ A first-order coupling) + Γ_damp(ξ) (the constraint-damping terms).

## (a) Connection part: growth 0 (as in Kerr), from the Killing–Yano tensor

Kerr–Newman is Petrov type D with a Maxwell field aligned to the principal null directions, and it admits a **Killing–Yano 2-form** Y. Along any null geodesic γ, the Penrose–Walker construction gives parallel-propagated objects: Y(γ̇) and *Y(γ̇), together with γ̇. With the conserved K = |Y(γ̇)|², these determine a frame along γ in which parallel transport is a rotation in a fixed 2-plane. On Γ the flow is quasi-periodic (B1), so the parallel transport on T𝓜, and therefore on every tensor bundle, is **bounded**. Hence the connection part contributes growth 0, for every r and in the full subextremal range. This is the same mechanism as in Kerr. **Residue:** cite the existence of a KY tensor for Kerr–Newman (standard; e.g. Frolov–Krtouš–Kubizňák, Living Rev. Rel. 2017) and the Penrose–Walker frame.

## (b) Constraint damping: choosable

The damping terms are at our disposal. In Hintz's construction their size is a free parameter, and B3 needs them only to be nonzero with the right sign at the horizon and at infinity. Choose them supported **away from the trapping region** (a cutoff in r), so that Γ_damp = 0 on Γ. **Check:** this is compatible with B3, the constraint-damping mode analysis. Hintz uses exactly this freedom in the Kerr case, so I expect no obstruction.

## (c) The genuinely new term: the Einstein–Maxwell coupling C(ξ)

Linearising gives:
- In Ric(g) = 2F∘F − ½g|F|², the right side contributes 2(F∘dȦ + dȦ∘F − ½ g⟨F, dȦ⟩) + (terms in h·F·F). These are **first order in Ȧ**.
- In Maxwell, d*F = 0 contributes (∇h)·F terms, **first order in h**.

So C(ξ) is an off-diagonal fibre endomorphism, homogeneous of degree 1 in ξ, of size ~ |F||ξ| ~ (Q/r²)|ξ|. It is **comparable to ν|ξ|-scale quantities on Γ**, so it cannot be dismissed as lower order.

**Structural fact.** Einstein–Maxwell is Lagrangian. Its linearised operator is formally self-adjoint for the (indefinite) fibre pairing G = (trace-reversed metric pairing on h) ⊕ (c·g⁻¹ on Ȧ), with an appropriate positive constant c. Hence C(ξ) is **G-Hermitian up to an imaginary factor**: σ₁(C)(ξ) = i·𝒞(ξ), with 𝒞 G-symmetric. Its contribution to transport is therefore G-unitary. A G-unitary flow does not by itself give bounded growth in a **positive** metric when G is indefinite: it can grow if G is indefinite on an invariant subspace.

**Threshold question (B2-c, precise).** On each trapped null covector ξ ∈ Γ, combine the zero-growth Killing–Yano transport with the coupling transport e^{s·i𝒞(ξ)}. Does the combined transport have **exponential growth rate < ν/2 in some positive metric**? Equivalently, are the eigenvalues of i𝒞(ξ), conjugated into the Killing–Yano frame, purely imaginary, or do they have real parts bounded by ν/2?

Notes:
- In Schwarzschild/Kerr (F = 0) the term is absent.
- In **Reissner–Nordström** the gauge-invariant Regge–Wheeler system of Giorgi has a positive conserved energy for the full subextremal range. This strongly suggests that the coupling is "non-growing" on Γ for a = 0.
- For a ≠ 0, Giorgi's generalised RW system (94)–(95) also has the coupling with the **adjoint structure in a positive metric** (conjugate coefficients, D· ↔ D⊗). This is the positive-metric analogue of what B2-c needs, and suggests that **working in the gauge-invariant RW variables at the trapped set** gives growth 0. Hintz's Kerr proof likewise transfers trapping estimates through gauge-invariant quantities in part.

**Plan for B2-c.** Check the threshold in Giorgi's RW variables, not in (h, Ȧ). On Γ the principal symbol is scalar (□), the first-order i(2a cosθ/|q|²)∇_t terms are purely imaginary times a real coefficient (no growth), and the coupling is 4Q² q̄³/|q|⁵ D· versus −½q³/|q|⁵ D⊗. Its symbol on Γ is ξ_θ-linear with conjugate coefficients. In a positive metric |p|² + w|q_F|² the coupling's subprincipal symbol is Hermitian (zero growth) for a suitable w > 0; the computation below gives w = 4Q². Then lift the trapping estimate from (p, q_F) back to (h, Ȧ) using the transport equations (90)–(91) and (78), which are non-trapping, since they are transports along e₃.
**(1) Done (see below).** (2) the lift back to the full system, which is a standard composition argument of the kind in Hintz's Kerr linear-stability paper.

**Status of B2:**
- (a) closed modulo citation;
- (b) design choice;
- (c) reduced to the anti-Hermitian check in Giorgi's weighting, which looks favourable, plus a lift.

### B2-c(1): the Hermitian check (done)

- Subprincipal symbol of the coupling: σ = i[[0, α ξ·],[β ξ⊗, 0]], with α = 4Q² q̄³/|q|⁵ and β = −½ q³/|q|⁵.
- Contraction identity: for a 1-form p and a symmetric traceless 2-tensor q, ⟨ξ⊗̂p, q⟩ = 2⟨p, ξ·q⟩.
- The adjoint of q ↦ iαξ·q, taken from (q-space, weight w) to (p-space, weight 1), is p ↦ −(iᾱ/(2w)) ξ⊗̂p.
- Hermiticity of σ requires this to equal iβ ξ⊗̂p, i.e. w = −ᾱ/(2β).
- Using ᾱ = 4Q² q³/|q|⁵, this gives w = −(4Q² q³/|q|⁵) / (2·(−½) q³/|q|⁵) = **4Q² > 0**.
- The phases cancel exactly because the coefficients are conjugate (Giorgi's structure).

The other first-order term, i(2a cosθ/|q|²)∇_t, has symbol −(2a cosθ/|q|²)ξ_t. That is real, hence Hermitian, so it causes no growth. Zeroth-order terms, including −(3/2)(H−H̄)⊗̂p, are below the threshold order.

⇒ **In the (p, q_F) variables with metric |p|² + 4Q²|q_F|², the subprincipal operator on Γ has growth 0. The threshold condition holds with full margin ν/2, for the full subextremal range.**

The remaining B2 work is the lift (2). The L-terms are first order in B, F and X, and those are given by non-trapping transports (the D2-type argument, now at the trapped set). The tensor-bundle Killing–Yano point (a) is still needed for the frame in which ξ· and ξ⊗ are transported.
