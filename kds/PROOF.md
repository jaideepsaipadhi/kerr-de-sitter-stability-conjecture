# Mode stability of subextremal Kerr–de Sitter: the proof as it stands (8 Nov 2026)

Status tags: **[P]** proved by hand; **[P-C]** proved with an exact or interval computer certificate (script named);
**[N]** numerical evidence only; **[OPEN]** not done.  "Item k" refers to `LEMMAS.md` in this directory, which has the
full detail and history; scripts live alongside it.

## 0. Goal and logical structure

**Target.** For every subextremal KdS black hole (Λ > 0; Δ_r = (r²+a²)(1−Λr²/3) − 2Mr has four distinct real roots
r₃ < 0 < r₀ < r₁ = r_e < r₂ = r_c):

* **Theorem A (Teukolsky mode stability).** The spin s = ±2 Teukolsky equation has no nonzero future-regular mode
  solution e^{−iσt+imφ}R(r)S(θ) with Im σ ≥ 0, except as allowed at σ = 0 below.
* **Theorem B (ungauged linearised Einstein mode stability, UEMS in the sense of Hintz–Vasy 2018 §4.1/§7).** Every
  generalised mode solution h of the linearised equation Ric′(h) = Λh, with Im σ ≥ 0 and smooth across both horizons,
  is pure gauge if σ ≠ 0, and linearised KdS (δM, δa) plus pure gauge (gauge 1-form polynomial in t) if σ = 0.

**Use.** Hintz–Petersen–Vasy (HPV25, not yet public) prove full-range nonlinear stability of KdS conditional on
mode stability for the *ungauged* linearised Einstein equation (as stated in Häfner–Hintz–Vasy, arXiv 2506.21183,
intro).  Theorem B is meant to be exactly that hypothesis.  Whether it matches HPV's statement verbatim is unverified.

**Architecture.**
A1. Reduce real-σ, σ ≠ 0, to a radial problem in Heun form, and split into "outside the superradiant window"
    (Wronskian identity) and "inside the window" (transform to a dual equation plus a two-current argument).
A2. The two-current argument reduces to an explicit inequality (the **event lemma**) in four parameters (t, z₀, A, B).
A3. The event lemma is proved in three regimes R1, R2 (analytic) and R3 (an explicit multiplier).
A4. σ = 0 separately; the window endpoints separately; Im σ > 0 by continuity from Schwarzschild–de Sitter.
B.  Theorem A ⇒ Theorem B via a Λ > 0 port of Andersson–Häfner–Whiting (AHW, arXiv 2207.12952).

---

## Part A: Teukolsky mode stability

### A.1 Radial Heun form and data [P]
Möbius map r₀, r₁, r₂, r₃ ↦ 0, 1, t, ∞, μ(z) = z(z−1)(t−z),
t = (r₁−r₃)(r₂−r₀)/((r₁−r₀)(r₂−r₃)) (item 26).  κ_j = iΞ K̂(r_j)/Δ′(r_j), K̂ = ω(r²+a²) − am, Ξ = 1 + α, α = Λa²/3.
Future-regular boundary conditions (item 43): R ~ (r−r₁)^{−κ₁−s} at r_e, R ~ (r−r₂)^{+κ₂} at r_c (no Δ^{−s}),
matching Casals–Teixeira da Costa (CTdC).  Exponents at infinity {1+2s, 1+s−2κ₃} (item 25).

### A.2 Outside the window [P, item 88]
For real ω the dual equation is real and W[R, R̄] is conserved; W(1) = W(t) gives
Im e₁ |c₁|² = t Im e_t |c_t|², e_j = 2(κ_j + κ₃) (item 34).  Where Im e₁, Im e_t have opposite signs, c₁ = c_t = 0 and
there is no mode.  **[P, item 88]** Im e₁ = 2Ξm f₁(ω/m), Im e_t = 2Ξm f₂(ω/m) with f₁ increasing and f₂ decreasing
(residue identities P₁+P₃ = P₀+P₂ and 1/D₁+1/D₃ = 1/D₀+1/D₂), roots 0 < w₂ < w₁; the gap w₁ − w₂ has sign of
2(r₀−r₁)(r₀−r₃)(r₁−r₃)(r₀r₃−r₁r₂) > 0.  So the sign-agreement set is exactly the window (w₂, w₁).

### A.3 The transform [P]
Inside the window, y ↦ ỹ = c · I₁^{μ+1}[y] (Riemann–Liouville/Pochhammer integral based at z = 1, μ = A − 2 with
A = 1+s−2κ₃) maps a quasinormal mode y to a solution of the dual Heun equation with the correct boundary branches
(items 17, 55).  Injective on the **closed** window: termwise, (z−1)^{ρ+k} ↦ Γ(ρ+k+1)/Γ(k+1−e₁) · (z−1)^{ρ+k+μ+1};
every factor is finite and nonzero, since ρ is non-real and k+1−e₁ is non-real or a positive integer (item 55).

### A.4 Exact form of the dual potential [P, item 31 symbolic]
At every finite m, for s = ±2:  L(z) = A/z + B − 3z with
 A = t(X₀² − (X₁+X_t)²) = 2t u₃(X₀+X₁+X_t) > 0 in the window (items 5, 22),
 B = λC₀ + c₁,  C₀ = 3/(Λ(r₁−r₀)(r₂−r₃)) (item 26),  c₁ = 1 + t + 3(1−α)C₀ (item 89, both halves exact; item 31: the numerators of
 n₁ − (λC₀ + c₁) and n₂ + 3 are divisible by s² − 4 as polynomials, after Vieta elimination).
**Sharp bound [P]** (items 26, 29, 72): with Vieta (Σr_i = 0, e₂ = −3(1−α)/Λ),
 2(1−α)C₀ − (2t−1) ∝ 3(r₀+r₁)² ≥ 0, so c₁ ≥ 4t − 1/2, and C₀ > (2t−1)/(3(1−α)) (strictly, from the positive
 definite form 4r₀² + 7r₀r₁ + r₀r₂ + 4r₁² + r₁r₂ + r₂²).  This holds for all Λ > 0 with four real roots.

### A.5 Angular bound [P-C, item 46]
On the window λ ≥ 3/2 for all subextremal KdS.  Inputs: α ≤ 7 − 4√3; ω/m ≥ 0 and aω/m ≤ 1 on the window (proved).
|m| ≥ 3 analytic (λ ≥ 6.6); |m| = 1, 2 by a two-level gap bound with exact Beta moments, interval arithmetic over
(α, ξ) (`step47_angular_proof.py`; certified minima 1.957, 1.665, 1.920, 1.914).
Hence the physical floor **B ≥ 5.5t − 1.25** (item 29).

### A.6 Reduction to the event lemma [P]
Cosmological side, all m: with f₊ = 1, C₊ = −q₁H ≤ 0 (Lemma H, item 6) and D₊ = qL′ ≤ 0 since L′ = −A/z² − 3 < 0
(item 19).  L(t) > 0 from A.4.  Event side (item 20): with h = f′/f, s = z₀ − z, H = 2μ + sμ′,
 **(C)** s h μ ≤ H,   **(D)** ½(μ′h + μ(h′+h²)) ≤ s(Lh + L′)   on (1, z₀),
plus f(1) = 0 and the seam condition.  The two-current argument (items 8, 12, 17) then forces the dual solution, and
hence (A.3) the mode, to vanish.  **Event lemma:** for all t > 1, z₀ ∈ (1, t), A > 0, B ≥ 5.5t − 1.25 there is h with
(C), (D).

### A.7 Event lemma, regimes R1 and R2 [P, item 47]
With p = 2(z₀−1)L(1)/(t−1): R1 (p ≤ 1) and R2 (p > 1, 2(z₀−1)²(A/z₀ + 3) < t−1) are closed by the explicit ODE
multiplier (μf′)′ = 2s(L̃f)′, f(1) = 0, with L̃ = L (R1) or L̃ = L − κ (R2).  Full proof in item 47.

### A.8 Event lemma, regime R3 (large A) — **PARTLY OPEN**
**Correction (item 71).** The multiplier of item 51 (δ = 1, p-term p(1−u)/x) **fails (C)** when A ≫ B at large t
(A/B ≥ 100, t ≥ 100: normalised margin +0.005 … +0.08).  Item 51's "margin ≥ 0.23" claim is false.
**Repaired multiplier:**  h = 1.7 g (1 − e^{−2κs²}) + p(1−u)/(xz),  g = −L′/L, κ = L(z₀)/μ(z₀), x = z−1,
u = x/x₀, p = ½ min(0.9, 2x₀L(1)/(t−1)).
* **[N]** 1200 cases, t ∈ [1.01, 1e5], z_f ∈ [0.05, 0.95], B/floor ∈ {1, 10, 10³}, A/B ∈ [0.01, 10⁶]: no failure.
  Worst C −0.082, worst D −0.033 (the latter at z_f → 0, the known degenerate corner with margin ≈ −0.65 z_f; item 33
  handles it by the rescaling x = (z₀−1)y).  `step61_retune.py`.
* **Corners (limit problems):**
  - **[P, leading order]** A → ∞: off the seam sLh and s|L′| dominate and Lh = ψ|L′|, so D/n_D → (1−ψ)/(1+ψ) =
    −δ/(2+δ) = −0.259; seam = the A-independent rescaled layer problem (items 13–15).  **[N]** all sequences converge
    (item 74).
  - **[N, explicit limits]** t → ∞ (item 73): inner zone z = O(1): D ratio ≤ −1/3 (computed in closed form at leading
    order); C reduces to 1.7(z−1)(gz) + p ≤ 2z − 1.  Outer zone z = tζ (A = t²α, B = tb): closed-form limits C∞, D∞
    (`step62_outer_limit.py`); sampled maxima C ≤ −0.081, D ≤ −0.196.  b → ∞ sub-limit [P]: C∞ = −ζ(ζ² − 3ζζ₀ + 2ζ₀) < 0.
* **[OPEN]** (i) explicit error bounds O(A^{−1/2}) and O(1/t) with inner/outer matching, so the corners hold for
  A ≥ A₀, t ≥ T₀; (ii) certification of the bounded part with the repaired multiplier (the earlier tile certificates
  `cert/r3core.py` used the old multiplier and are void).

### A.9 σ = 0 [P, items 35, 43]
m ≠ 0: ω = 0 lies outside the window since w₂ > 0 [P, item 88], so A.2 applies.
m = 0: s = −2 energy identity with V = −2Λr² − 2(1−α) − λ < 0 after the Teukolsky–Starobinsky map Δ²∂⁴; s = +2 via
Δ²R₊₂; the cubic kernel is excluded by a 4×4 determinant [λ² + 6(1−α)λ + 8 − 4α + 8α²]² > 0 (item 43 gives the
corrected boundary conditions).

### A.10 Window endpoints [P, items 44, 88]
Im e₁, Im e_t are affine in ω and vanish at the distinct points w₁ ≠ w₂ (item 88), so exactly one vanishes at each
endpoint.  Each case closes by the Wronskian plus injectivity.

### A.11 Upper half-plane [P, argument, item 36]
QNMs form a discrete set depending continuously on (Λ, a) (Petersen–Vasy JEMS 2024, full subextremal Fredholm theory),
with none for Im σ > C, locally uniformly.  At a = 0 (SdS) there are none in Im σ ≥ 0, σ ≠ 0.  A mode entering
Im σ > 0 must cross the real axis, which A.2–A.10 exclude.  Separation of variables is used only at real σ, where the
angular operator is self-adjoint (item 48).

---

## Part B: ungauged Einstein mode stability (items 37–42, 58–70)

### B.1 Teukolsky scalars vanish [P given Theorem A]
For a generalised mode, δψ₀ and δψ₄ are generalised Teukolsky modes; the top coefficient vanishes by Theorem A and the
lower ones inductively.  Because h̄ is a mode with frequency −σ̄ (also Im ≥ 0), all **four** complex conditions
dC(l,m,l,m) = dC(l,m̄,l,m̄) = dC(n,m̄,n,m̄) = dC(n,m,n,m) = 0 hold (item 64).

### B.2 Local structure: B1 (Λ-analogue of AHW Prop 6.5) [P-C at sample points; uniformity mostly certified]
d_J = dim(J-jets of solutions of {linearised Λ-vacuum, four conditions}) − dim(gauge J-jets).  Exact upper bounds via
truncated jets; if d_{J₀} = d_{J₀+1} the gauge-reduced symbol vanishes and the data form a flat connection
(item 58(b)).  For modes the jets are 2-D in (r, x = cos θ) and exact over Q(i) (`step52_exact_modes.py`).
* At sample points d₀ = d₁ = 0 for every tested (σ, m) except σ = 0, m = 0, where d = 4 = PD tangents (item 64).
* **One good point per mode suffices** (item 65(c)): open-dense good set + Khavkine's complete compatibility operator
  + simple connectivity of the exterior.
* **σ-uniformity, fixed m [P-C]:** matrices are polynomials of degree ≤ 2 in σ; exceptional polynomials at two points
  have gcd c·σ⁴ mod P; lifted to an exact statement by the reduction lemma of item 67 with δ_C computed exactly
  (multimodular CRT + Hadamard bound, `step57b_exact_deltaC.py`), gauge minors pinned by the tropical bound.
  Done exactly for m = 1, 2.
* **Joint (σ, m) uniformity (items 68–70):**
  - E_i(σ, μ) = σ⁶μ⁶H_i(σ, μ), with σ⁶ and μ⁶ divisibility **exact** (block-Toeplitz deficiency + interpolation,
    `step58f_jordan.py`, `step60_mu6.py`); exact degree bounds ≤ 80 make the 81×81 mod-p interpolation exact.
  - Three points: gcd of the three resultants Res_σ(H_i, H_j) mod p = μ⁸¹; H_i(σ, 0) has a 9-fold root at σ = 0.
  - Lemma: if ord_{μ=0} Res_σ(H₁, H₂) = 81 exactly, no integer m ≠ 0 has a σ ≠ 0 exceptional at all three points (any
    such m must be ≡ 0 mod p, and then the cofactor's value at m is ≡ its value at 0, a p-unit).  No large-|m| asymptotics.
  - **[OPEN, compute handed off]** (a) exact σ-degree of E_i = 29 (`run_degcert.py`, ~50 min);
    (b) exact multiplicity of E_i at the origin ≥ 21 (`run_origin_cert.sh`, ~2–3 h; gauge minors already exact: 2, 4;
    constraint minor needs 15, mod-p scan saturates at 15 from Toeplitz length 4).
* **B1 identified (item 40) [N, floating point]:** the four Λ-Plebański–Demiański tangents (δM, δa, NUT, acceleration)
  span the σ = 0, m = 0 local space modulo gauge.  Exact version is the m = 0 count in item 64.

### B.3 Axis regularity replaces AHW Prop 6.8 [P, item 41]
Elementary flatness at both poles gives two conditions on (δn, δα) with determinant −64 m p₀⁴/c ≠ 0, so NUT and
acceleration are excluded for smooth modes.  Globalisation by simple connectivity (AHW Lemma 6.6).

### B.4 Assembly [P given B.1–B.3, item 42]
σ ≠ 0: (L_T + iσ)^{J+1} kills h and is invertible on the stationary PD part, so h is pure gauge.
σ = 0: h = c_M ∂_M g + c_a ∂_a g + δ*X with X polynomial in t.  This is the UEMS form.

---

## Summary of what is NOT yet proved

1. **A.8 R3 corners**: explicit error bounds and matching (the genuine remaining analysis).
2. **A.8 R3 bounded part**: interval certification with the repaired multiplier.
3–5. ~~A.2, A.9, A.10 numerical sign claims~~ — proved by hand (item 88).
6. **B.2**: the two handed-off exact certificates (degree 29; origin multiplicity 21).
7. **HPV match**: Theorem B against HV18's UEMS wording, and HPV25 itself once public.


---
## UPDATE 8 Nov (items 88–99): status after the Λ-Wald route
* Theorem A: A.2, A.9 (m ≠ 0), A.10 proved by hand (item 88); c₁ identity fully symbolic (item 89).
  Still open: R3 corner error bounds and the bounded R3 part.
* Theorem B: the jet certificate / bad set (items 75–87) is SUPERSEDED.
  - (σ, m) ≠ (0,0), a ≠ 0: every mode with vanishing Teukolsky scalars is globally pure gauge, explicit smooth gauge vector
    (item 96; in-house referee item 97, five repairs written in).  Mechanism: linearised Goldberg–Sachs; the Killing-spinor
    vector ξ and the Carter field ζ are first-order Killing fields; Im δξ is a background Killing field, hence 0 for modes.
  - σ = m = 0: Carter structure to first order; field equations linear in (Δr, Δx); 4 = PD tangents; axis kills NUT,
    acceleration (item 98).
  - Generalised modes: induction (item 99); one small open point at σ = m = 0 (no linear-in-t mass/spin growth).

---
## FINAL STATUS after REFEREE2 (8 Nov, evening)
* Theorem A: PROVED for real σ (all regimes R1–R3, σ = 0, endpoints, κ_{1,2} = 0 boundary cases). Im σ > 0:
  continuity along the roots-scaling path to SdS, CONDITIONAL on high-frequency QNM control (REFEREE2 §C).
* Theorem B = HPV Assumption 1.2 as written (REFEREE2 §D): PROVED given Theorem A.
* Remaining: REFEREE2 §C (one analytic input); Pochhammer write-up of item 17; cite/prove α ≤ 7 − 4√3; independent check.

## STATUS after item 124
Theorem A: PROVED, all σ with Im σ ≥ 0. The UHP part uses HPV's high-frequency bound (local uniformity in b by openness;
to be confirmed against HPV §3), Petersen–Vasy, and Wald's Hertz identity with Λ.
Theorem B: PROVED (needs Theorem A).
α ≤ 7 − 4√3: PROVED. Left: Pochhammer exposition of item 17; independent check.
