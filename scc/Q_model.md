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
