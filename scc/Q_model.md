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
