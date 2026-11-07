# Christodoulou-type SCC near Kerr: scoping (Nov 2026)

## Correction to my earlier answer

I said the interior blow-up step was wide open. **That is out of date.** Two papers, both from spring 2026:

- **Luk–Sbierski, arXiv:2604.04877 (Apr 2026)** — *The formation of a weak null singularity in the interior of generic rotating black holes.*
  - Setting: a characteristic IVP on H⁺ ∪ H₁, with data settling down to subextremal, strictly rotating Kerr.
  - Hypotheses: (i)–(ii) upper bounds (Dafermos–Luk type); (iii)(a) a weighted L²-averaged **lower bound on the ℓ = 2 modes of the dynamical s = +2 Teukolsky field ψ₀ on H⁺**; (iii)(b) auxiliary upper bounds on the ℓ > 2 modes and on T ψ₀, with slightly improved rates. Rates: q = 11 (conformally regular I⁺ class, modified Price law) or q = 7 (Gajic–Kehrberger class); see their Remark 4.8.
  - Conclusion: a weak null singularity near i⁺, where the metric is **C⁰-extendible but not C^{0,1}_loc (Lipschitz) extendible**.
  - Method: Dafermos–Luk C⁰ stability, plus Sbierski's linear s = +2 instability (Ann. PDE 2023, strengthened in arXiv:2604.06283), plus Sbierski's Lipschitz-inextendibility criterion (arXiv:2409.18838). The nonlinear correction is shown to be perturbative.
- **Gurriaran, arXiv:2603.17911 (Mar 2026)** — the same conclusion (Lipschitz SCC near i⁺), via Klainerman–Szeftel gauges and the Ma–Zhang-type route, under a nonlinear Price-law-type estimate.

## What "closing SCC near Kerr" needs, item by item

| # | Requirement | Literature status | New idea needed? |
|---|---|---|---|
| E1 | Exterior converges to Kerr, quantitatively, full subextremal range | **Hintz (nonlinear Kerr, full range)** | no |
| E2 | Construct the event horizon of Hintz's solution and put the data on H⁺ into Dafermos–Luk / Luk–Sbierski event-horizon gauge with the upper bounds (i), (ii), (iii)(b) at rate q | Hintz, *Horizons of some asymptotically stationary spacetimes* (arXiv:2411.12568, Duke); Dafermos–Luk [22] (near-horizon region). The match of **rates** (q = 7 or 11 vs. Hintz's decay rates) still has to be checked | probably no: bookkeeping between gauges and rates |
| E3 | **Generic nonlinear lower bound (iii)(a)**: for an open dense set of Cauchy data near Kerr, ψ₀^{ℓ=2} on H⁺ satisfies the weighted L² lower bound | Linear: sharp Teukolsky asymptotics with a generically nonzero coefficient (Ma–Zhang, slowly rotating; Millet arXiv:2302.06946, full range). Nonlinear: **open**. Nonlinear tails may change the rate (Luk–Oh [73], Kehrberger), and this is why q = 11 rather than Price | **moderate**: see plan below |
| I1 | Interior near i⁺: weak null singularity, Lipschitz-inextendible | **Luk–Sbierski 2026; Gurriaran 2026** (conditional on E2, E3) | no |
| I2 | **Global interior**: the whole bifurcate Cauchy horizon, including the region near the bifurcation sphere and both ends for two-ended data; every boundary point singular; no other piece of boundary | Dafermos–Luk I covers only the region near i⁺. The global bifurcate version is announced or planned, not available as far as I can tell | **hard but expected**: global extension of the Dafermos–Luk scheme |
| F1 | **Formulation gap.** Christodoulou's SCC asks for inextendibility with **Christoffel symbols in L²_loc**, which is weaker regularity than Lipschitz. Luk–Sbierski exclude only C^{0,1} extensions | Low-regularity inextendibility theory: Sbierski (C⁰ for Schwarzschild, Lipschitz for weak null singularities). Uniqueness of extensions **fails** at C⁰ ([5] in Luk–Sbierski). Nothing is known at the Christoffel-L² level for weak null singularities | **YES: this is the genuine new-idea gap** |
| F2 | Mass-inflation-type blow-up (the spin −2 / Ψ₄ side), the natural mechanism behind F1 | Linear: Gurriaran arXiv:2503.24114 (generic linearised curvature singularity, spin −2). Nonlinear: open | likely needed for F1 |

## Assessment

- After Hintz, Luk–Sbierski and Gurriaran, **Lipschitz SCC near i⁺** is conditional only on E2 and E3, which are exterior-side, perturbative items. Doing them gives: *generic small perturbations of subextremal Kerr data develop a weak null singularity near i⁺ across which the metric is Lipschitz-inextendible.* This is the most **closable** target. It is a real theorem, but it is a gluing of existing pieces plus one genuinely nonlinear estimate (E3).
- **Christodoulou's formulation** needs F1, and I2 for the global statement. F1 is where a new idea is required: an inextendibility argument that works for extensions with only L² Christoffel symbols, at a weak null singularity in vacuum without symmetry. Even in spherical symmetry this level was reached only via mass inflation and Hawking-mass blow-up (Luk–Oh; Gautam). Outside symmetry there is no Hawking mass, so a substitute monotone or blow-up quantity is needed.

## Plan for E3 (generic nonlinear lower bound)

1. **Smooth dependence.** Hintz's Nash–Moser solution map, data ↦ solution, is C¹ (indeed tame-smooth) into weighted spaces. **Check:** that the late-time leading coefficient functional is continuous, or C¹, in those spaces.
2. **Leading asymptotics.** Hintz's low-frequency resolvent and polyhomogeneity framework gives ψ₀|_{H⁺} ~ c(data)·v^{−q/2…} (+ logs) for the ℓ = 2 projection. Here c(data) = ℓ_lin(data) + O(‖data‖²), with ℓ_lin the linear coefficient (Millet / Ma–Zhang), plus a **nonlinear tail** contribution that may dominate the rate. This is the subtle point: the q = 11 vs Price discrepancy.
3. **Genericity.** If the coefficient functional c at the correct (dominant) rate is C¹ with nonzero differential at Kerr, then {c ≠ 0} is open and dense near Kerr data. If the nonlinear tail dominates and its coefficient is quadratic, nonvanishing needs a separate argument: show the quadratic form is not identically zero, which is a computation.
4. **Weighted L² lower bound.** Pointwise asymptotics with c ≠ 0 imply (iii)(a).

## Plan for F1 (the new idea), first thoughts only

- The blow-up actually proved is curvature blow-up in weighted L² along the Cauchy horizon (ψ₀, s = +2). For Christoffel-L² inextendibility one needs a quantity that is (a) **finite** for any extension with Γ ∈ L²_loc and (b) **infinite** for our solution.
- Candidates: (a′) a **renormalised Hawking-mass / quasi-local mass** on the double-null spheres, generalising the spherical mass inflation; (b′) a **geodesic-deviation / focusing integral** along transverse null geodesics, which is ∫ Ric-type and could be controlled by ‖Γ‖_{L²}; (c′) Luk–Oh–Shlapentokh-Rothman's scattering-theoretic blow-up of the **spin −2** field (Ψ₄, which captures mass inflation), combined with an integrated Raychaudhuri argument.
- In vacuum, **mass inflation itself has not been proved without symmetry**. Spin −2 blow-up is only linear (Gurriaran 2025). So F1 likely needs **F2 (nonlinear mass inflation in vacuum)** first. That is the true core.
