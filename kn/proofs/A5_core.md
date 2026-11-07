# A5: real modes in the superradiant band (the core problem)

## Statement

For p ∈ 𝒫, m ∈ ℤ∖{0}, and real ω with 0 < ω/m < Ω_H(p): there is no outgoing, horizon-regular solution of the linearised Einstein–Maxwell equations of the form e^{−iωt+imφ}(…).

By the bridge, this is equivalent to: no nonzero gauge-invariant solution (A, F, B, X) or (p, q_F) of that form.

## Why every tool used so far fails here

1. **Flux balance.**
 ω²|rad_{I⁺}|² + ω(ω − mΩ_H)|rad_{H⁺}|² = 0.
 In the band the second coefficient is negative, so the identity is consistent with both amplitudes being nonzero. The K-energy (K = ∂_t + Ω_H∂_φ) gives the same identity multiplied by (ω − mΩ_H)/ω, so it carries no new information.
2. **Unique continuation from either end (C1–C3, route α)** needs one end to carry zero radiation. Flux balance does not provide that.
3. **Whiting's transform (Kerr)** is an integral transform in (t, r) built on the separated radial and angular ODEs. Kerr–Newman does not separate. In the gauge-invariant coupled system, the coupling coefficient q³/|q|⁵ mixes r and θ, so the angular operator depends on r.

## A dimension count, and why A5 is a genuine risk rather than an unproven formality

- For fixed m, a real mode is a zero of a holomorphic Fredholm determinant: one complex equation, i.e. **2 real conditions**.
- The unknowns are (a, Q, ω) ∈ 𝒫 × ℝ: **3 real unknowns**.
- So generically the set of real modes is a **curve** in (a, Q, ω)-space, unless some structure forbids it. Compare Kerr: 2 unknowns (a, ω) against 2 conditions gives generically isolated points, and Whiting shows the set is empty.
- If such a curve exists in Kerr–Newman, it is exactly the locus where a resonance crosses into the upper half-plane. Then KN would be **mode unstable** in the region beyond it, and A7's continuity argument would correctly fail.

So A5 is not "certainly true, just unproven". The structure that would forbid the curve has to be found. The known physics evidence points to stability: Dias–Godazgar–Santos (DGS) 2015 found all QNMs with Im ω < 0 across the parameter grid they computed. But that is numerical evidence on a grid, and a curve of crossings could slip between grid points only if it were very short, which is unlikely.

## Candidate ideas, ranked

**(i) Symmetric-hyperbolic / "twisted" energy at fixed real frequency.**
At fixed real ω, the mode problem for Giorgi's (p, q_F) system is a stationary problem on the (r, θ) half-strip. Look for a sesquilinear current J = Re⟨u, 𝒜u′⟩ + … whose divergence is ≥ 0 in the band, such that its boundary fluxes at H⁺ and I⁺ have the **same** sign (not opposite signs, as with T and K). In Kerr, Whiting's transform implicitly produces such a current at the level of separated ODEs: it is the conserved Wronskian of the transformed equation, where the potential has no superradiant sign problem.
Test case: does the separated-Kerr version of such a current exist in **local** form, i.e. as a differential, not integral, operator in r? If yes, it might lift to KN without separation.
*Known obstruction:* Andersson–Ma–Paganini–Whiting (2017) indicate that Whiting's transform is genuinely nonlocal. So expect (i) to need a pseudodifferential current, i.e. a commutator/Carleman-type estimate in a frequency-localised sense.

**(ii) Real-frequency Carleman-type estimate in the ergoregion with the r-weight chosen frequency-dependently.**
Shlapentokh-Rothman's quantitative mode stability for Kerr (2015) gives explicit resolvent bounds on the real axis. Reproving them with estimates that do not use separation, only a microlocal splitting into "superradiant" and "non-superradiant" frequencies, is exactly where Lemma A helps: in the superradiant band there is **no trapping**. So in the band the real-frequency problem is **non-trapping**. One may hope to get a resolvent bound from a Mourre-type or escape-function argument. That would rule out real modes, since a real mode is a resolvent pole. The idea is: in a non-trapping region, a real mode would need its energy to be "stuck", but every null bicharacteristic escapes to H⁺ or I⁺. Then propagation of singularities at the radial sets gives zero radiation, and route α/C3 finishes.
*Gap:* this argument controls **high** frequency only. At bounded real ω there are no bicharacteristics to speak of. Lemma A is semiclassical. For bounded ω the escape-function method gives nothing.

**(iii) Perturbation from the Kerr axis plus "no-crossing" monotonicity.**
At a resonance on the real axis, compute ∂_Q Im σ at Q = 0, starting from Kerr's known (Whiting) absence of real modes. This alone cannot exclude a crossing curve away from Q = 0.

**(iv) A coupled-system Whiting transform.**
Teukolsky–Starobinsky identities are unavailable. But for **real** ω the coupled radial–angular problem is a bounded-from-below elliptic-in-θ family. One could try to separate in a weaker sense: project the angular dependence onto the eigenbasis of the r-dependent angular operator, giving an infinite ODE system with off-diagonal couplings, and apply a matrix Whiting transform. The off-diagonal terms are proportional to Q²a, which is small in no regime of the full range, so this likely fails beyond small Qa.

## Assessment

None of (i)–(iv) is a proof. (ii) is the most promising at high |ω|, (i) at bounded ω, and neither covers both. **A5 at bounded real ω in the band is the genuine new-idea gap in KN mode stability. Everything else in the programme is either proved, reduced, or expected.**

## The cheapest decisive check

Because of the dimension count, a short targeted computation can tell whether we are trying to prove something true:
- use the validated coupled spectral solver (`kn/numerics/qnm_coupled.py`, validated against Schwarzschild and Kerr a = 0.99);
- track the **least-damped** m = 1, 2 coupled modes along a coarse (a, Q) grid near the extremal boundary, where superradiance is strongest;
- look for Im ω changing sign.
DGS covered much of this already. The value of rerunning it is to target the near-extremal corner, where our own near-extremal analysis found coupled transitions at a* ≈ 0.36–0.5.
