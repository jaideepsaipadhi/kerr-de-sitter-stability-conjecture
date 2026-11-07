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

---
## Test of idea (i): local currents, in the simplest model (Kerr, s = 0, separated)

**Model.** u'' + V(r*)u = 0, with real V, V → ω² as r* → ∞ and V → k² as r* → −∞, where k = ω − mΩ_H < 0 in the band.
An outgoing mode has u ~ A_I e^{iωr*} at infinity and u ~ A_H e^{−ikr*} at the horizon.

**General local quadratic current.**
J = 𝒜|u'|² + ℬ Re(ūu') + 𝒞|u|² + 𝒟 Im(ūu').
Using the ODE, J' = (𝒜' + ℬ)|u'|² + (ℬ' − 2𝒜V + 2𝒞)Re(ūu') + (𝒞' − ℬV)|u|² + 𝒟' Im(ūu').
- J' ≥ 0 for all u forces 𝒟' = 0. The 𝒟-part is then just the conserved Wronskian, i.e. the flux identity ω|A_I|² = |k||A_H|², and contributes nothing new.
- The real part 𝒬 = 𝒜|u'|² + ℬRe(ūu') + 𝒞|u|² must have a PSD derivative. Its endpoint values on pure exponentials are c_I|A_I|² and c_H|A_H|², with c_I = 𝒜_∞ω² + 𝒞_∞ and c_H = 𝒜_H k² + 𝒞_H.

**Condition for a contradiction.** 𝒬(∞) ≥ 𝒬(H) together with the flux ratio |A_I|²/|A_H|² = |k|/ω gives a contradiction iff
 **c_I |k| < c_H ω.**
In the simplest admissible family (ℬ = 0, 𝒞 = 𝒜V, needing 𝒜' ≥ 0 and (𝒜V)' ≥ 0), this becomes 𝒜_∞ ω < 𝒜_H |k| with 𝒜 nondecreasing. That is possible only if ω < |k|, i.e. **ω < mΩ_H/2, the lower half of the band**, and only if V stays positive wherever 𝒜 > 0. The Kerr potential is negative in a barrier region, so this fails. With ℬ ≠ 0 the feasibility question becomes a Riccati / Sturm-comparison problem.

**Verdict on (i).**
- A local current proving A5 would, already for Kerr s = 0, give an elementary proof of Kerr mode stability on the real axis.
- The decay literature does not contain one. Dafermos–Rodnianski–Shlapentokh-Rothman use exactly this class of currents (y|u'|², h-currents) and still had to call Whiting for bounded superradiant frequencies.
- Moreover, for *general* real potentials with these asymptotics, real superradiant modes do exist (1-D scattering with a negative-energy horizon end). So any proof must use fine structure of V, not just the pointwise sign conditions available to local currents.

⇒ **Local currents: very likely dead.** A proof of A5 must be **nonlocal** in r, as Whiting's is, or must avoid working at fixed ω altogether.

## Where this leaves A5

The realistic options are:
- **(N) Numerics, to decide truth.** By the dimension count, a crossing curve, if it exists, is visible with a targeted sweep near the extremal corner: least-damped m = 1, 2 coupled modes, looking for a sign change of Im ω. The validated solver already exists.
- **(W) A nonlocal transform for the coupled system.** Find an integral transform in (t, r), at the PDE level in (r, θ), that maps Giorgi's (p, q_F) system to one with no superradiant sign. Whiting's kernel is explicit for Kerr. The first concrete step is to write Whiting's transform in a form that does not separate θ (an operator-valued kernel) and see where the KN coupling breaks it.
- **(P) Accept A5 as a hypothesis**, numerically supported, and finish the rest of the programme conditionally. That would give "KN is nonlinearly stable for all subextremal parameters, provided mode stability holds". This is the same logical shape as Hintz's Kerr theorem, where mode stability is a black-box input, except that for KN the input is unproven.
