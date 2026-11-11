# Hostile referee report on PROOF.md / LEMMAS.md (8 Nov 2026)

Checks run (scratch scripts in /tmp/ref, no repo files modified): sympy on A.4/item 72; B.3 determinant by hand;
the R3 repaired multiplier (`step32_four_param.check(..., 0.7, p, 2.0, zdecay=True)`) on 1500 random log-uniform
points plus a corner grid; root-geometry sweeps (200k samples) for the window/sign structure; `step39_c1_symbolic.py`
(it timed out at 240 s after verifying only the n2 half).

## (I) Fatal or likely fatal

**I.1 B.2 is checked at ONE black hole only, but Theorem B claims all subextremal KdS.** Items 61, 64–70: every exact
rank, exceptional polynomial, resultant and gcd is computed at Λ = 1/20, a = 3/5, M = 1, at three points (r, x). Nothing
covers other values of (Λ, a): there is no symbolic dependence on the geometry, no exceptional-set argument in (Λ, a), and
no semicontinuity argument (item 59 even says semicontinuity goes the wrong way). So Theorem B is proved, at best, for one
geometry. PROOF.md's "uniformity mostly certified" and the summary list hide this. Fix: carry (Λ, a) as extra symbols, or
use a fixed-geometry argument plus an exceptional-variety analysis in (Λ, a, σ, m). This is a much larger computation.

**I.2 The jet-closure step (item 58(b)) does not follow, and its open status is not reported.** From d₀ = d₁ = 0 (true
dimensions, and only an upper bound at K = 4) the item concludes "the gauge-reduced symbol vanishes ... by prolongation at
every higher order". That is not valid. Equal dimensions at orders J₀ and J₀+1 say nothing about order J₀+2 unless the
gauge-reduced system is formally integrable or involutive. Proving that is exactly item 57, which is still **[OPEN]**
(Cartan–Kuranishi or an explicit Khavkine compatibility order), and it is missing from PROOF.md's "NOT yet proved" list.
Concretely: h is locally pure gauge iff K₁(h) = 0, where Khavkine's K₁ has some order N. To get K₁(h)(q) = 0 you need the
N-jet of h at q to be a gauge jet, i.e. d_N = 0. You do not need d₀ or d₁. N is not known. The uniformity certificates
(items 66–70) only control J ∈ {0, 1}, so even with I.1 fixed they would not prove B1 unless N ≤ 1. Fix: get an explicit
bound on the order of K₁ (or prove involutivity at an explicit order), then certify d_J = 0 up to that order uniformly.

**I.3 m = 0 with σ ≠ 0 is not covered at all, and σ = 0 with |m| ≥ 3 is not covered.** E_i = σ⁶μ⁶H_i (items 70, 71).
So E_i(σ, 0) ≡ 0 for every σ and E_i(0, m) ≡ 0 for every m: on both lines every point is "exceptional", and the
three-point argument gives no good point there. The lemma in item 70 only treats m ≠ 0, σ ≠ 0. What exists elsewhere:
- σ = 0 exact checks only for m = 0, 1, 2 (items 64, 66);
- m = 0, σ ≠ 0 only at two exact sample values (σ = 1/2, i/4; item 61).

Axisymmetric modes, which is where the PD/NUT/acceleration directions live, and stationary modes with |m| ≥ 3 are
therefore unproved. B.4 needs both. Fix: use different minors that do not vanish identically on these lines (work with
H_i, or redo the compression on the slices m = 0 and σ = 0), then a separate fixed-m (or fixed-σ) certificate.

## (II) Serious gaps or mislabelled status

**II.1 The continuity argument in A.11 uses a false geometric premise.** Item 36 says "for fixed Λ the subextremal a-range
is an interval containing 0". But item 45 records that subextremal KdS reaches ΛM² ≈ 0.177 > 1/9. For 1/9 < ΛM² < 0.177,
a = 0 is not subextremal, so there is no path to SdS at fixed ΛM². Fix: use a path in the 2-parameter region (it must be
shown connected and to reach the SdS segment), and get uniform Fredholm/high-frequency bounds along that path. Also
needed and not written:
- the counting step (Gohberg–Sigal/Rouché on 0 ≤ Im σ ≤ C, with |Re σ| large excluded uniformly);
- that Petersen–Vasy's theory applies to the s = ±2 Teukolsky operator in the function spaces where "future-regular"
  separated modes live (threshold regularity depends on s and Im σ);
- that the Fredholm resonant states at real σ are exactly the separated ODE modes of A.1.

**II.2 The A.4 c₁ identity is not reproduced, and the ledger contradicts itself.** PROOF.md labels it [P, item 31
symbolic]. Item 72, written later, says "The c1 identity itself (item 22) is still numerical". I ran `step39_c1_symbolic.py`
for 240 s: it printed "n2 + 3: remainder mod (s²−4) = 0" and a numeric n₁ check, but never reached the symbolic n₁
divisibility. That half is the one that supplies B = λC₀ + c₁, which everything downstream needs. Fix: run it to
completion, archive the log, and fix the contradictory text in item 72.

**II.3 The z_f → 0 corner of the event lemma has no uniform margin and is not proved.** At t ∈ [1+10⁻⁶, 10⁷] with
z_f = 10⁻⁶, the worst normalised D margin is −7.4·10⁻⁷ (my grid), consistent with "≈ −0.65 z_f". The margin goes to 0, so
no interval certificate of the "bounded part" can close near z_f = 0. A.8 says item 33 "handles it by the rescaling", but
item 33 is [Reduction + N] with "TODO (math)". By item 33, z_f → 0 is the window edge, so it is a physical regime. It must
be listed in PROOF.md §Summary as a separate open corner lemma.

**II.4 B.4 assembly has gaps for generalised modes and for the form of the gauge field.**
(a) Item 42(2) says the local structure holds "for ANY linearised solution". It was established (at samples) only for exact
modes e^{−iσt+imφ}H(r, x) (items 60–61). Generalised modes with t^j terms need an induction: write the top coefficient
h_J = δ*X_J; subtract δ*(t^J e^{−iσt} X_J); note that the commutator term from δ* acting on t^J drops the degree. None
of this is written.
(b) At σ = 0, a top coefficient c·t^J ∂_M g (J ≥ 1) has to be excluded, or shown to be of UEMS form. Not addressed.
(c) UEMS in HV18 needs the gauge 1-form to be a generalised mode with the same σ. B.4 only obtains some X with δ*X = h
(defined up to Killing fields, from simple connectivity). Fix: project onto the e^{−iσt} t^j sector, using that
[L_T, δ*] = 0 and that the Killing algebra is finite-dimensional.

**II.5 HPV hypothesis match is unresolved, and the confidence figure is a guess.** Item 63 puts 30% on HPV needing the
gauge-fixed statement, i.e. constraint-damping 1-form modes plus a gauge-fixed zero-energy space. PROOF.md §0 calls
Theorem B "meant to be exactly that hypothesis". Until HPV25 is public, the honest claim is UEMS in the HV18 sense, and
nothing more.

**II.6 Item 70: what the "m ≡ 0 mod p ⇒ contradiction" lemma silently needs.** The logic is correct as stated: common root ⇒
all three resultants vanish at m ⇒ m mod p is a root of μ⁸¹ ⇒ p | m ⇒ T(m) ≡ T(0) ≢ 0. Beyond exactness of ord = 81
(equivalently mult₀ E_i ≥ 21, giving a ≥ 9·9), it also needs all of the following:
- exact deg_σ H_i, i.e. formal degrees honest (run_degcert, handed off);
- R'₁₂ ∈ Z[i][μ] computed with the same formal degree over Z and mod p;
- exact σ⁶ and μ⁶ divisibility (done, items 69(i) and 71(a));
- the bound "ord over Z ≤ ord mod p", which is valid only for the honest-degree reduction.

Also, the inequality a ≥ I_(0,0)(H₁, H₂) uses ord_μ Res = Σ of intersection multiplicities on the line μ = 0. With a
leading-coefficient drop at μ = 0 this is still ≥, so it is fine. The lemma is not wrong, but by I.1–I.3 it only covers
one geometry, J ≤ 1 and m ≠ 0, σ ≠ 0. Both certificates are still owed (PROOF.md says so).

**II.7 A.3: well-definedness of the transform is asserted, not proved.** μ + 1 = A − 1 has Re = s, i.e. −2 for s = −2, so
the Riemann–Liouville integral diverges. It has to be the Pochhammer/analytically continued version. The claim that its
continuation to z = t lands on the correct (analytic) branch of the dual is item 17, still "[P, modulo write-up]". The
termwise injectivity of item 55 also needs κ₁ ≠ 0, i.e. ω < mΩ_H strictly on the closed window. That inequality is only
sampled (item 34, "window subset [0, Ω_H]"); it is not proved.

**II.8 A.1/A.2: integer exponent differences.** At ω = mΩ_H (κ₁ = 0) and ω = mΩ_c (κ₂ = 0), the horizon exponents differ
by the integer s, so log terms can appear. The Frobenius coefficient c₁ and the flux identity in A.2 then need a separate
check. Both frequencies lie outside the window, so A.2 is the step affected, and it is written without this case.

## (III) Minor

- Item 72's "sharp form" is wrong as written. Sympy gives
  C₀ − (2t−1)/(2(1−α)) = 3(r₀+r₁)² / (2(1−α)(r₁−r₀)(r₂−r₃)).
  The item has Λ where (1−α) should be. Sign and conclusion are unchanged.
- A.2, A.9 and A.10 are labelled numerical, but they are elementary and should just be proved:
  - Im e₁ and Im e_t are affine in ω, with slopes (r₁²+a²)/Δ′₁ + (r₃²+a²)/Δ′₃ > 0 and (r₂²+a²)/Δ′₂ + (r₃²+a²)/Δ′₃.
    The second slope was negative on all 200k sampled geometries, including near-extremal ones; it needs a one-line
    proof. Given that, the sign-agreement set is the interval between the two zeros, which is the window.
  - Its lower endpoint equals item 46's "lo" (r₀, r₁ form). I checked this numerically; it is a Σ r_i^k/Δ′(r_i) = 0
    identity.
  - A.9: at ω = 0 the zeros reduce to sign(1/Δ′₁ + 1/Δ′₃) vs sign(1/Δ′₂ + 1/Δ′₃). These are opposite because
    (r₀−r₃)(r₁−r₃) > r₂² > (r₂−r₀)(r₂−r₁) (no failures in 200k samples).
  - A.10 needs only that the two zeros differ (window nondegenerate), plus strict containment in (0, mΩ_H).
- A.5 cites α ≤ 7 − 4√3 as a "known bound" without a reference or proof.
- The seam condition holds only with equality: h(z₀) = 0 (ψ(z₀) = 0, u = 1), so f′(z₀) = 0. Item 47 says f′(z₀) > 0 for
  R1/R2. Confirm that the two-current argument (items 8, 12) accepts equality.
- `step32_four_param.py` docstring: the zdecay branch uses p(1−u)/((z−1)z) while the docstring says p(1−u)/(z−1). PROOF.md
  A.8 has the correct form.
- Notation: μ is used for the Heun weight z(z−1)(t−z), for A − 2 (A.3), and for the symbol replacing m (B.2).
- The B.3 determinant uses the PD parameter "m", not the azimuthal number. Say so.

## (IV) Checked and found correct

- A.4: C₀ − (2t−1)/(3(1−α)) = f/(Λ(r₁−r₀)(r₂−r₃)Q), Q = −e₂ > 0, with f the stated form. The form is positive definite
  (minors 4, 15/4, 7/2). Hence c₁ ≥ 4t − ½ and B ≥ 5.5t − 1.25, given λ ≥ 3/2 and 1 − α ≤ 1.
- B.3: det = −64 m p₀⁴/c, verified by hand from the stated dF1 and dF2.
- A.8 A → ∞ bulk limit: D/n_D → (1−ψ)/(1+ψ) = −δ/(2+δ) = −0.2593 at δ = 0.7. The algebra is right; error bounds are still
  owed, as stated.
- R3 repaired multiplier: no failure on 1500 random points (t up to 10⁶, z_f ∈ [10⁻³, 1), B/floor up to 10⁷, A/B from
  10⁻⁶ to 10⁹) or on the corner grid (t = 1+10⁻⁶ … 10⁷, z_f = 10⁻⁶ … 1−10⁻⁶, A/B up to 10¹²). Worst C −0.0816; worst D
  degenerates at z_f → 0 (see II.3).
- Lemma H endpoint values; item 67's reduction lemma (valid as stated for fixed m); B.1's use of h̄ to get the four
  complex conditions.
