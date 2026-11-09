# Second hostile referee pass, end to end (8 Nov 2026, evening) — in-house, not independent

Scope: PROOF.md, LEMMAS.md items 1–121, the write-up, the certificate code (cert2/), REFEREE.md (first pass).
New checks run this pass: cert2/test_bern.py (Bernstein enclosure soundness, 40 random polynomials, 6 random splits each:
0 failures); step122_table_pd.py (GHP tetrad table on the accelerating NUT Λ-PD metric, where τ ≠ −τ̄': max deviation
1.6e-41); ζ = K·ξ on the NUT PD metric (α = 0) at two points: constant −0.99504 ∂σ (Killing).

## A. Status of every first-pass objection
| item | first pass | now |
|---|---|---|
| I.1 B.2 at one geometry | fatal | SUPERSEDED: Theorem B no longer uses jets (items 90–103); argument is uniform in (Λ, a), a ≠ 0 |
| I.2 jet closure / Khavkine order | fatal | SUPERSEDED (no jets) |
| I.3 lines m = 0, σ = 0 | fatal | RESOLVED: item 96 covers (σ, m) ≠ (0,0) incl. both lines; (0,0) by item 98 |
| II.1 continuity path | serious | PATH FIXED (B below); the analytic input remains OPEN (C below) |
| II.2 c1 identity | serious | RESOLVED, item 89 (exact, 7 s) |
| II.3 zf -> 0 corner | serious | RESOLVED: the new R3 construction is certified on the full compact box incl. ζ -> 0 (items 119, 121) |
| II.4 B.4 assembly | serious | RESOLVED: items 99–101 (induction, t*-polynomial gauge, Iyer–Wald charges); gauge 1-form explicit |
| II.5 HPV match | serious | RESOLVED, see D |
| II.6 item 70 lemma | — | MOOT (superseded) |
| II.7 transform | serious | κ1 ≠ 0 on the CLOSED window now PROVED: w1 = (P1Ω1 + P3Ω3)/(P1 + P3) < Ω1 since Ω3 < Ω1 (|r3| > r1), so ω = mΩ_H is outside; Pochhammer write-up of item 17 still to be written out (not a mathematical gap: standard) |
| II.8 integer exponents at ω = mΩ_H, mΩ_c | serious | RESOLVED: both outside the window (item 88); in the real dual the z = 1 exponents differ by e1 = 2(κ1 + κ3), and at κ1 = 0, Im e1 = 2 Im κ3 ≠ 0 (K(r3) = m ρ3(Ω1 − Ω3) ≠ 0), so no logarithms; there f1 > 0 > f2, the Wronskian kills the mode. Same at κ2 = 0 |
| minor: item 72 Λ vs (1 − α) | minor | noted; sign/conclusion unchanged |
| minor: seam equality | minor | MOOT in R3 (new construction has f'(z0) > 0); R1/R2 have f'(z0) > 0 (item 47) |
| minor: α ≤ 7 − 4√3 | minor | still cited, not proved here — a standard bound; to write |

## B. The continuity path (II.1) — fixed
Use the roots as parameters, as Hintz 2608.25868 §6: keep r1, r2 (horizons) fixed and scale r0 -> s r0, s ∈ (0, 1],
r3 = −(r0 + r1 + r2).  Ordering and distinctness are preserved, a² = −(Λ/3) r0r1r2r3 ≥ 0 with Λ > 0 from the Vieta
quadratic (positive root always exists), and s -> 0 gives a = 0 (SdS).  So every subextremal KdS is joined to SdS
inside the subextremal set.

## C. OPEN: the analytic input of A.11 (upper half plane) — the one remaining dependency
A.11 needs, for the s = ±2 Teukolsky operator on full-range KdS: (i) discreteness and continuity of QNMs in the
parameters; (ii) no QNMs escaping to infinity in Im σ ≥ 0, locally uniformly.  (i) follows from Vasy's framework (the
subprincipal part only shifts threshold regularities; Petersen–Vasy 2608.26034 §2 set up the Fredholm theory for the
Teukolsky operator, including KdS).  Invertibility for Im σ ≫ 1 is standard.  But (ii) near the real axis at high
frequency needs trapping estimates, which require a subprincipal-symbol bound at the trapped set; HPV 2508.06620
Thm 3.1 verify such a bound for the gauge-fixed linearised Einstein operator, NOT for Teukolsky.  Two ways to close:
  (a) verify the subprincipal condition for Teukolsky at the trapped set of full-range KdS (as HPV do for L_b), or
  (b) run the continuity argument for HPV's gauge-fixed operator L_b (their Thm 2.2: discrete QNMs, finitely many in
      Im σ ≥ 0; Thm 3.1: trapping bound) and use their Thm 2.7 mechanism (constraint damping) to transfer.
Until one of these is written, Theorem A is proved for REAL σ (all σ ≠ 0 and σ = 0), and Im σ > 0 rests on (C).
Consequence for Theorem B: item 96 is proved for every mode with Im σ ≥ 0 whose Teukolsky scalars vanish; the
vanishing for Im σ > 0 uses Theorem A there, hence (C).

## D. Theorem B vs HPV Assumption 1.2 (II.5) — RESOLVED
HPV25 is public (arXiv 2508.06620).  Assumption 1.2: (1) σ ≠ 0, Im σ ≥ 0, (DRic − Λ)u = 0, (L_T + iσ)u = 0 ⇒ u = δ*ω;
(2) (L_T)^k u = 0 ⇒ u = g'(b') + δ*ω; T = ∂t* + (a/(r0² + a²))∂φ*.  Note u need NOT be a ∂φ-eigenfunction.
Repair to item 96 for (1), with no Fourier decomposition: J[u] is a background Killing field commuting with T, and
L_T J[u] = J[L_T u] = −iσ J[u], so σ ≠ 0 ⇒ J[u] = 0.  Then ξ, ζ are first-order Killing fields; T = αξ_R + βζ with
constants (ξ_R = −Ξ∂t, ζ = −Ξ(a²∂t + a∂φ), a ≠ 0), so −iσ u = L_T u = −L_{αδξ + βδζ} g and u = δ*(ω), ω the dual of
(αδξ + βδζ)/(iσ) — smooth, local in u.  For (2): the φ-average u0 is stationary and axisymmetric up to the t*-ladder
(items 98–101); u − u0 has (L_T)^k(u − u0) = 0 with every φ-Fourier mode of ∂t-frequency Ω0 m ≠ 0, so the same
first-order-Killing argument applied to (L_T)^{k−1}(u − u0) and the t*-ladder of item 99 make it pure gauge.
(This is exactly HPV's own reduction in their proof of Thm 2.7.)

## E. New parts (items 88–121): checked
* 88, 89: re-derived; exact.
* 90–97 (Theorem B, σ, m): the GHP table now verified on a NON-KdS type D metric (τ ≠ −τ̄'), removing the τ/−τ̄
  ambiguity flagged in item 93; the item-91 relations on the same metric (item 97 R5); ζ Killing on NUT PD.
  Linearisation: every identity used is a consequence of the linearised Ricci/Bianchi identities with
  δκ = δσ = δκ' = δσ' = δΨ0,1,3,4 = 0 and Φ = 0 — valid for any linearised Λ-vacuum solution in the tetrad gauge.
* 98–103 (σ = m = 0): Carter-ansatz linearity and solution space exact (item 103).
* 117–121 (R3): (C) on phase 2 by the Wronskian — sign chain re-derived; (P1) Descartes after Möbius with the
  certified sign pattern; (P2) bootstrap re-derived line by line (Φ(ẑ) = f̂Ĥ/ŝ − 2ŝL̂f̂; Ψ' = 2s1L'(z1) < 0;
  Ψ(z0) = Ĥ/ŝ − 2∫s|L'|; the u − log(1+u) ≤ u²/2 bound; K ≤ 0 ⇒ N < 0 via monotonicity in A and B).
  Requirements of the two-current argument (item 12/20): f(1) = 0 (F_C = μ/s² vanishes linearly), f > 0, (C), (D),
  seam f'(z0) ≥ 0 — all met.
* Certificates: exact integer arithmetic; enclosure soundness tested; the only floating point is in chart
  bookkeeping (box labels), never in a sign decision.  Denominator signs are read off fixed-sign factors and confirmed
  at random interior points (a lemma, not a sampling claim, because the factors are products of fixed-sign terms).

## F. Verdict
Theorem A: PROVED for all real σ (every regime of the event lemma, σ = 0, endpoints, outside the window).
  Im σ > 0: by continuity, CONDITIONAL on the analytic input C.
Theorem B (= HPV Assumption 1.2): PROVED given Theorem A, for every subextremal KdS (a ≠ 0; a = 0 is HV18).
So the final result is complete modulo ONE literature-level input: high-frequency control of s = ±2 Teukolsky QNMs in
Im σ ≥ 0 on full-range KdS (or the gauged-Einstein version of the continuity argument).  Plus independent checking.

## G. Update (item 123): the §C gap reduced
The counting argument along the path is PROVED given compactness (K); σ = 0 handled by openness of invertibility.
(K) and the SdS base case are reduced, via Teukolsky–Starobinsky + Wald reconstruction + gauge fixing, to HPV's
high-frequency bound for L_b, plus three concrete sub-lemmas: (C2) horizon regularity of S†Φ, (C3) non-vanishing
of the TS constant in Im σ ≥ 0 at large |σ| and the finitely many algebraically special frequencies, (C4) local
uniformity of HPV Thm 2.2.  No Teukolsky trapping estimate is needed any more.
