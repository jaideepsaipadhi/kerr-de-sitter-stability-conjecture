# (Q) in a model: vacuum plane waves with a weak null singularity

## The model

Rosen form of a vacuum plane wave on {U < 0}:
 g = −2 dU du + γ_AB(U) dx^A dx^B,  γ = Eᵀ E,  E(U) ∈ GL(2).
- Vacuum (R_UU = 0) means E'' = H E, with H(U) symmetric and traceless (H is the tidal tensor, i.e. the α of this model).
- Shear and expansion are given by χ_E := E' E⁻¹ (its traceless part is χ̂, its trace is trχ).
- **Weak null singularity at U = 0:** choose E such that χ_E ~ U⁻¹ |log(−U)|^{−1−2δ}. Then E ∈ C⁰ ∩ W^{1,1}, E' ∉ L², and the Rosen metric is C⁰ up to U = 0 with Γ ∉ L². Exactly as in Dafermos–Luk, the shear is L¹ but not L². The area factor satisfies Raychaudhuri (trace of E'' = HE), and stays continuous because χ̂ ∈ L¹, with |χ̂|² merely producing an integrable blow-up of trχ.

This is the cleanest setting with the same singular structure as the Kerr Cauchy horizon (one null direction, shear L¹ \ L², metric C⁰).

## Symmetry-adapted charts

Look for an extension chart of the form X = Ẽ(U) x (u and U unchanged up to lower order), with Ẽ ∈ C⁰ invertible. Put M := E Ẽ⁻¹ and K := Ẽ' Ẽ⁻¹. A direct substitution, dx = Ẽ⁻¹(dX − K X dU), gives

 g = −2 dU du + |M (dX − K X dU)|²,

so in the chart (U, u, X):
- g_XX = MᵀM;
- g_UX = −MᵀM K X;
- g_UU = |M K X|².

**Continuity of g̃ at U = 0** (for X ≠ 0 in the chart neighbourhood) forces **M ∈ C⁰ and K ∈ C⁰.**
**Γ̃ ∈ L²_loc** forces:
- ∂_U g_XX ∈ L², where ∂_U(MᵀM) = Mᵀ[(χ_E − K) + (χ_E − K)ᵀ]M, using M' = (χ_E − K) M;
- ∂_U g_UX ∈ L², which contains K' X and so needs **K ∈ H¹**.

Hence χ_E − K ∈ L² with K ∈ H¹ ⊂ L², and therefore **χ_E ∈ L²**. This contradicts χ_E ∉ L².

> **Within symmetry-adapted charts, (Q) holds: the singular shear cannot be gauged away in L².**

## The mechanism, which is the new structural point

The "gauge" quantity K that would cancel the singular shear appears **undifferentiated in the metric** (in g_UX, multiplied by the transverse coordinate). Continuity of g̃ therefore makes it continuous, and Γ̃ ∈ L² makes it H¹. An H¹ (hence L²) gauge can never cancel a non-L² shear. In the general reduction "χ_DL = pure gauge + L²", the pure-gauge part is (∂Φ)⁻¹∂²Φ. **The claim to prove in general is that its shear component is controlled, in H¹, by quantities that appear undifferentiated in g̃.** In the model, that quantity is the mixed U–transverse metric component.

## The general version to prove next

**Conjecture (Q′).** In any extension chart (C⁰ metric, Γ̃ ∈ L²) near a point of ι(CH⁺), let V be the generator field and Y the transverse Jacobi fields, which are bounded by Lemma F1. Then the transverse connection along V,
 𝒦 := (Y⁻¹ ∂_U Y)_{transverse part}, computed in the chart,
is **continuous**, and moreover H¹ along generators. Equivalently, 𝒦 is determined by the chart components g̃(V, Y) and g̃(∂_U Y, Y) up to L² errors.

**Why it is plausible.** In the model, 𝒦 = K is read off from g_UX = g̃(∂_U, ∂_X). In general, g̃(V, ∂_a) is a continuous chart quantity. The DL structure gives g(∂_U, ∂_θ) = 0 (orthogonality of the double null frame), so ∂_U[g̃(V, Y)] = 0 identically. Expanding this in the chart, g̃(∂_U V, Y) + g̃(V, ∂_U Y) + (∂_U g̃)(V, Y) = 0. With ∂_U g̃ ∈ L² along generators (Fubini) and ∂_U V = −Γ̃VV ∈ L², this gives **g̃(V, ∂_U Y) ∈ L²**: the component of ∂_U Y **along V** is L². The danger is the remaining **transverse** component of ∂_U Y, which is χ̂. That is exactly the component the model shows must be a continuous gauge plus L².

**Remaining gap for (Q′).** Show that the transverse part of ∂_U Y (in the chart) cannot have a non-L² piece unless a **metric component** discontinuity appears. The model's argument uses the transverse coordinate X as a *multiplier*. In general the analogue is to use **Lemma F1's bi-Lipschitz bound in the transverse directions**, i.e. variations of the generators across neighbouring θ, as a substitute for linearity in X.

## Status

- (Q) holds in the plane-wave model for symmetry-adapted charts. The mechanism is identified: the gauge must be continuous because it enters the metric undifferentiated.
- Next: remove the symmetry-adapted assumption in the model, i.e. arbitrary charts on the plane wave, using Lemma F1 plus Heisenberg-group translations to average. Then transfer to Dafermos–Luk via (Q′).

---
# Removing the symmetry assumption: what happens for general charts

## 1. Why the symmetric argument does not extend directly

The symmetric proof used linearity, X = Ẽ(U)x, which made the gauge term K = Ẽ'Ẽ⁻¹ appear **undifferentiated** in g_UX, as K·X. For a general chart Φ⁻¹: (U, u, X) ↦ (U, u, x), that link is gone.

**The genuine gauge freedom.** Continuity of g̃ = (∂Φ⁻¹)ᵀ g (∂Φ⁻¹) only fixes the Jacobian J = ∂Φ⁻¹ **up to a pointwise Lorentz transformation**: J = Λ J₀ with Λ(x) ∈ O(g) measurable. Then
- Γ̃ ∈ L² constrains Λ⁻¹∂Λ, up to compensation;
- J must still be a Jacobian, so curl(Λ J₀) = 0.

**So (Q) for general charts is a Lorentzian rigidity problem:** can a Lorentz-valued field Λ, subject to curl(ΛJ₀) = 0 and Λ⁻¹∂Λ ∈ L² modulo compensation, produce the non-L² rotation needed to cancel χ̂?
- The Riemannian analogue is **rigid** (Reshetnyak; Friesecke–James–Müller).
- The Lorentzian one is **not rigid in general**, because O(1,3) is non-compact and the problem is hyperbolic.

This is a real obstruction to a cheap proof, and a real possibility for a counterexample.

## 2. A frame-invariant reduction (proved)

χ̂ is a tensor on the screen bundle V^⊥/V.
- Under the residual frame freedom it changes only by:
 - null rotations about V: χ̂ is invariant;
 - screen rotations: |χ̂| is invariant;
 - boosts: χ̂ is rescaled by a factor that Lemma F1 bounds above and below.
- So |χ̂|² is invariant up to bounded factors.
- In the chart, χ̂(e_A, e_B) = g̃(e_A(V^a)∂_a, e_B) + Γ̃(e_A, V, e_B), with e_A and V bounded (Lemma F1).

Hence

 **(Q) ⇔ the chart components V^a of the generator field have L² derivatives in the screen directions** (modulo L²).

Equivalently, **the optical function u, which is Lipschitz in the extension because ∇u = −V is bounded, must lie in H² in screen directions.** Its Hessian is ∇²u = ∂²u − Γ̃∂u, and its screen part is −χ.

## 3. Candidate route via semiconvexity, and why it is circular as it stands

- In smooth spacetimes, achronal boundaries such as ∂J⁺(S) have null second fundamental form bounded **above**: light cones of points slightly to the past on a generator touch from inside.
- Our singular shear is traceless-dominated, so one eigenvalue of χ → **+∞**. An upper bound would therefore give a contradiction.
- **But** the comparison cone's own χ is governed by the same Riccati equation along the same generator, with forcing α, so the comparison constant degenerates exactly as α blows up.
- Without an a-priori constraint on α from Γ̃ ∈ L², this is circular.
- **What would break the circularity:** a low-regularity version of the comparison. Light cones of a C⁰ metric with Γ̃ ∈ L² should have screen second fundamental forms in L² near their vertex, averaged over cones (Fubini over vertices). If true, comparing the averaged cone bound with the eigenvalue → +∞ of N_u's χ gives (Q). **This is the cleanest candidate new idea, and it is chart-free.**

## 4. Honest status

| claim | status |
|---|---|
| (Q) in the plane-wave model, symmetry-adapted charts | **proved** |
| reduction of (Q) to screen-H² regularity of the optical function in an arbitrary extension | **proved** (sketch) |
| general (Q) | **open**: equivalent to a Lorentzian rigidity / null-hypersurface regularity statement |
| the averaged light-cone L² lemma (§3) | **conjectural**: proposed new idea |

## 5. Next decisive test: try to build a counterexample in the plane wave

Use Brinkmann-type coordinates and a **null rotation about V whose parameter oscillates in the screen directions at a U-dependent frequency.** Null rotations are the only part of O(1,3) that Lemma F1 does not control. Then check whether curl(ΛJ₀) = 0 together with Γ̃ ∈ L² can coexist with screen derivatives of V^a ∉ L².
- **If yes:** this is a counterexample model, and Christodoulou's formulation is in danger.
- **If the curl constraint forbids it:** that computation is the rigidity proof, and it should lift.

---
# The null-rotation test turns into a proof mechanism: the KORN ARGUMENT

## Idea

Let ω be a covector field on the extension chart. **Korn's inequality** is a statement about the **flat** symmetric gradient sym∂ω = ½(∂_a ω_b + ∂_b ω_a) in chart coordinates. It does not involve the metric g̃ at all, so it is **elliptic whatever the Lorentzian signature**. Specifically (Nečas–Lions lemma / second Korn inequality on a bounded Lipschitz domain D ⊂ ℝ⁴): if ω ∈ L²(D) and sym∂ω ∈ L²(D), then ∂ω ∈ L²(D).

**The point:** it gives a way to control the antisymmetric part of a derivative using only the symmetric part. The singular shear of a Cauchy horizon lives in an antisymmetric slot of a suitable field. That is the rigidity statement §1 asked for, obtained from the flat chart structure rather than from the Lorentz group, which is not rigid.

## The field: the Jacobi (Lie-dragged) screen field Y

Take Y = ∂_θ in Dafermos–Luk (DL) coordinates, so [V, Y] = 0 and ∇_V Y = ∇_Y V = χ(Y). Compute ∇Y♭ in the DL null frame (V, e₃, e_A):
- (V, e_A) slot: (∇_V Y)(e_A) = χ(Y, e_A).
- (e_A, V) slot: (∇_{e_A} Y)(V) = −g(Y, ∇_{e_A}V) = −χ(e_A, Y).
- So the **symmetric part vanishes in the (V, e_A) slot**, since χ is symmetric because V is hypersurface-orthogonal. The **antisymmetric part equals χ(Y, ·)**, which is singular, ∉ L².
- The other slots of sym∇Y♭:
 - (V, V): g(∇_V Y, V) = 0.
 - Screen–screen: these are angular derivatives of γ_AB(u, U, θ), which are integrals in U of ∂_θχ ~ U⁻¹log^{−1−2δ}, hence **bounded**.
 - Slots involving e₃: these involve χ̲, η and ζ, which are **regular at CH⁺_A** away from the bifurcation corner, at least in L². **Check against the DL norms.**

⇒ **sym∇Y♭ ∈ L², but the antisymmetric part of ∇Y♭ ∉ L².**

## The contradiction

Suppose an extension chart has g̃ ∈ C⁰ and Γ̃ ∈ L². Let ω := chart components of Y♭. Then
 sym∂ω = sym∇̃ω + Γ̃·ω.
The first term is the tensor sym∇Y♭ ∈ L² transported by the frame change. The second is L²·bounded. Korn gives ∂ω ∈ L², hence ∇̃ω = ∂ω − Γ̃ω ∈ L². Its antisymmetric part is χ(Y, ·) ∈ L², which **contradicts χ̂ ∉ L²**. ∎ (modulo the gap below)

## Check in the plane-wave model

Here Y can be taken to be the Killing field ξ = ∂_x, for which sym∇ξ♭ = 0 **exactly**, while ∇ξ♭ ∝ f′e^{2f}(dU ∧ dx). Korn then forces f′ ∈ L², a contradiction. This **covers all charts, not just symmetry-adapted ones**, subject to the same gap.

## THE GAP (now precise and quantitative)

Korn needs **ω ∈ L²** and the frame change to map L² tensors to L² chart components. Both follow if **Y and V have chart components in L∞** (bounded above, and bounded below where inverses are needed). Lemma F1 gives this only **along almost every generator**, with constants exp(C ∫|Γ̃∘Φ|dU) that are **not uniform** across generators. Γ̃ ∈ L² is supercritical in 4D, critical being L⁴, so naive bootstrapping (Korn + Sobolev) loses integrability at each step.

**Ways to close it:**
1. **Smallness + exceptional set.** On a small enough neighbourhood, ∫|Γ̃∘Φ| dU ≤ ε except on a set E of generators whose measure → 0. Off E the bounds are uniform. Use a **Korn inequality with a removed set**: on John domains with holes, or in a weighted / capacity form. Here the work is to show that E cannot carry the whole non-L² shear. That should follow because the shear lower bound is **ℓ = 2-averaged over the sphere**, so it lives on a set of positive measure that is not concentrated on small sets.
2. **Maximal-function control.** sup over U of ∫|Γ̃∘Φ| is dominated by a strong maximal function of |Γ̃|², whose weak-L¹ bound might suffice together with the exponential only on the exceptional set.
3. **Linear (not exponential) bounds.** Use the null structure: d(1/F)/dU = −(Γ̃VV-components), whose integral is **linear** in Γ̃. So F's lower bound degrades only linearly, which is fine. The upper bound of F is the only exponential place, and it can be handled by choosing the normalising field T adapted to V.

## Status

- **(Q) has a candidate proof via Korn's inequality, which is new to this problem as far as I know.** The rigidity that the Lorentz group lacks is supplied by the flat chart structure.
- It is complete in the model and in the DL setting **modulo**: (a) the uniform-boundedness gap above, and (b) checking that the e₃-slots of sym∇Y♭ are L² in the DL norms.
- If (a) and (b) close, then **Christodoulou's formulation holds near i⁺ wherever Luk–Sbierski's lower bound gives χ̂ ∉ L².**

---
# Gap (a): what closes, what doesn't

## Improvement 1: Korn in L^p, and the shear is not even in L^{1+}

- Korn's inequality holds in L^p for every 1 < p < ∞ (on Lipschitz domains; Nečas, Ornstein's counterexample being p = 1).
- The singular shear χ̂ ~ U⁻¹|log U|^{−1−2δ} is in L¹, but **not in L^p for any p > 1**.
- So we only need sym∂ω ∈ L^p and ω ∈ L^p **for some p > 1**.
- With Γ̃ ∈ L², the term Γ̃·ω lies in L^p as soon as ω ∈ L^q with q > 2.

So the uniform L∞ requirement drops to:

> **(A_q)** the chart components of V, Y and the inverse frame lie in L^q_loc for some q > 2.

## Improvement 2: linear (not exponential) control where the null structure allows it

- **F = −g̃(V, T).** Since d(1/F)/dU = −(Γ̃VV)/F²·F² is linear in Γ̃, we get
 1/F ≤ 1/F₀ + ∫|Γ̃∘Φ|dU.
 The **lower** bound on F is therefore polynomially integrable (q < ∞ for all q, by Fubini). The **upper** bound fails exactly when ∫|Γ̃∘Φ| dU reaches 1/F₀, which happens on a set of generators of measure ≤ ε²|U₀|F₀²·‖Γ̃‖²_{L²}.
- **Y** satisfies a linear ODE with coefficient |χ| + |Γ̃||V|. Here |χ| is uniformly L¹ (this is the DL bound), and |Γ̃||V| is the only problematic part.

## What remains

On the **exceptional set E** of generators, V may concentrate in the chart (F → ∞), and Y may grow exponentially. **(A_q) on E is not implied by Γ̃ ∈ L² alone**, because L² is supercritical in 4D (the critical exponent is L⁴). This is a genuine analytic gap, not bookkeeping.

## What is proved: a conditional theorem that already goes beyond Lipschitz

> **Theorem (conditional on (A_q) and the e₃-slot check).** Assume the Luk–Sbierski lower bound, so that χ̂ ∉ L^p_loc for every p > 1 near a point of CH⁺. Then the DL interior admits **no** extension with g̃ ∈ C⁰, Γ̃ ∈ L²_loc, in which the DL generator and Jacobi fields have chart components in L^{2+ε}_loc.

This is strictly stronger than Lipschitz inextendibility:
- Lipschitz extensions have bounded geodesic tangents, so they satisfy (A_q) automatically. Luk–Sbierski's theorem is the special case **Γ̃ ∈ L∞**.
- Here Γ̃ is only L², and the frame is only L^{2+}-integrable, not bounded.

## Route to removing (A_q): two concrete options

1. **Prove (A_q) from the eikonal structure.** V = −∇u, with u Lipschitz in the extension (u is continuous and its level sets are achronal null hypersurfaces). A **Lipschitz** u gives ∂u ∈ L∞, hence V = −g̃⁻¹∂u ∈ **L∞**, with no ODE needed. So the real question is **whether the DL optical function u extends to a Lipschitz function on the extension**.
 - For achronal null hypersurfaces in a C⁰ spacetime, the level sets are Lipschitz graphs (Chruściel–Grant).
 - The **spacing** between level sets, i.e. |∂u| bounded and bounded below, needs a uniform lower bound on the "null distance" between N_u and N_{u'}. That is continuity of the time-separation structure, which a C⁰ metric does provide (Lorentzian length structure).
 - **This looks like the right route.** It replaces the generator ODE (supercritical) by causal-theoretic Lipschitz bounds, which need only g̃ ∈ C⁰.
 - Y can then be handled the same way through the angular DL coordinates θ, which are Lie-dragged along V: θ is constant along generators, so ∂θ is controlled by the Lipschitz structure of the generator foliation.
2. **Exceptional-set Korn.** Use Korn with a removed set E of small measure and show the ℓ = 2-averaged shear lower bound cannot concentrate on E.

**Next:** route 1, i.e. prove that u (and θ) extend as Lipschitz functions near p in any C⁰ extension. If true, (A_∞) holds, the theorem becomes unconditional (modulo the e₃ check), and **(Q) is settled: Christodoulou's formulation holds near i⁺ under the Luk–Sbierski lower bound.**

---
# Route 1 (Lipschitz extension of u) checked: it collapses back onto gap (a)

**What C⁰ causality gives.** In a chart where g̃ is C⁰-close to Minkowski, every causal curve, parametrised by chart time t, is Lipschitz with a uniform constant. Achronal null hypersurfaces are Lipschitz graphs (Chruściel–Grant). So the level sets N_c of u are **individually** Lipschitz graphs.

**What it does not give.** Lipschitz-ness of u also needs the **spacing** between N_c and N_{c+Δ} to be ≳ Δ in the chart. Equivalently, the DL parameter u along the ingoing null generators (the e₃ curves, which near p run along CH⁺) must be comparable to chart time. Its chart speed F₃ = −g̃(e₃, T) satisfies the same linear-in-Γ̃ transport. Fubini gives control only along **a.e.** e₃-curve, and the supercritical problem returns.

At the pure C⁰ level this cannot be fixed. Continuous extensions across weak null singularities are known to be non-unique; that is the failure of uniqueness cited in Luk–Sbierski. So some C⁰ extensions do reparametrise the DL coordinates wildly. **Any proof of (A_q) must use Γ̃ ∈ L².** Γ̃ ∈ L² is supercritical (L⁴ is critical in 4D), which is exactly why the generator/Jacobi ODEs give only a.e. control.

## The core problem, stated sharply

> **Proposition (A) (open).** Let (M̃, g̃) be a C⁰ extension of the DL interior near p ∈ ι(CH⁺) with Γ̃ ∈ L²_loc. Then the DL coordinate fields (∂_U, ∂_u, ∂_θ) have chart components in L^q_loc for some q > 2 near p.

**(A) + Korn ⇒ (Q) ⇒ Christodoulou's SCC near i⁺** under the Luk–Sbierski lower bound. This is the conditional theorem in the previous section.

**Why (A) is plausible.**
- The DL fields satisfy linear ODEs along generators, with coefficients |Γ̃| (L² in 4D) plus quantities that are **uniformly** L¹ in the DL parameter (χ, ω̲, ζ).
- A wild reparametrisation would need Γ̃ to concentrate along whole null generator segments. That is a 1-dimensional concentration in 4D, which **is** compatible with Γ̃ ∈ L². So (A) is **not** automatic.
- **Candidate counterexample shape:** a chart that "pinches" a thin tube of generators near p. Γ̃ ~ (tube radius)^{−1} on a tube of cross-section ρ³ (3 transverse dimensions × U-length) gives ‖Γ̃‖²_{L²} ~ ρ³·ρ⁻² = ρ, which is small. A countable union of such pinches with Σρ_k < ∞ is in L². So an L²-Christoffel chart **can** distort the DL frame wildly on a null set of generators.
- **But** Korn only needs L^q integrability (q > 2) of the frame, not boundedness. Pinches with frame distortion ~ ρ^{−s} on tubes of volume ρ³ contribute ~ ρ^{3−qs} to ∫|frame|^q. So whether (A) can fail is a **quantitative** competition between the Christoffel cost of a pinch and the frame-distortion gain.

**Next step.** Compute that competition exactly in the plane-wave model. Construct the optimal "pinch" chart (a localised reparametrisation of the generators) and compare its Γ̃-L² cost with its frame-L^q distortion.
- If every pinch with finite Γ̃-cost has frame distortion in L^{2+}, then (A) holds, plausibly via a covering argument.
- If not, the pinch is the mechanism by which Christodoulou's formulation could fail.
