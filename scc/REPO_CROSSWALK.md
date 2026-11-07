# Crosswalk: this session's Level-3 work vs the-strong-cosmic-censorship-conjecture repo (read 2026-11-07)

Main sources:
- scc/working/KERR_LEDGER.md
- KERR_K3_DMC.md (§§35–54)
- paper_k4/{shear_obstruction, chart_upgrade, riccati_monotone, null_normal_form, unified_holonomy}.tex
- done/scc_C5.md
- working/BOOST_CLOSE(+REFEREE).md

| This session | Repo counterpart | Effect |
|---|---|---|
| Lemma F1, u-line lemma (bounded e₃/u-side boost; transverse cost ≥ ∫dU/U) | **Lemma A′** (K3 §35, reset counting; "uniform e3/u-chart cap ⇒ bounded tilt, cost ≳ χ²") | Repo version PROVED (sketch), refereed. Supersedes ours. |
| Correction: the Korn singular term is χ·e^{−s} (boost-weighted) | **Chart shear identity** χ̲ = λ[g̃(∂_{e_A}L̂, e_B) + Γ̃] (shear_obstruction.tex), plus Thm: (L-REG_sh) ⇔ finite λ-weighted shear cost | Confirms the correction. Repo: this hypothesis is "as hard as the target". |
| Open boost/θ issue (e₄-side headroom, "A bounded?") | **Lemma B** (§§42–54): e₄ has V-growing headroom (≲ 1/(Vℓ)). Optimal-control numerics give an integrable cost for a tilted sawtooth frame on one slice. A V-consistency argument (§53) suggests the cheap strategy fails across slices. The 2-D test is inconclusive. | Our remaining gap **is** Lemma B. |
| Topology / collapse lemma | **Lemma C** (§39): push-up (causal plainness) ⇒ injectivity along CH generators, which excludes collapse. Corners = bubbling = push-up failure. Cameron–Sbierski Ex. 3.3-type bubbling metrics can lie in W^{1,2}. | Closes the non-bubbling case. |
| Purely geometric Korn route | **Warped escape** (chart_upgrade.tex §, TWO_CHECKS Check 2): a non-vacuum example with shear ∉ L² that **has** a C⁰ ∩ W^{1,2} extension (sketch, with an NEC caveat) | **Kills any vacuum-free proof**, including the Korn argument as I posed it. Vacuum must enter, e.g. via the Riccati identity ∇_L(Jχ̂) = −Jα (riccati_monotone.tex). |
| Item 3 corner (bifurcation sphere) | **Two-sided Hardy rigidity** (§50A): non-L² shear on both CH components near S_B ⇒ cost ≥ ∫dV/(UV) = ∞ | Closes the inextendibility side of the corner (sketch), given both-sided blue-shift input. |
| α ∉ H⁻¹ test (F-a) | **C5** (done/scc_C5.md): linear, gauge-invariant, any gauge, via Gurriaran Ψ̂₊₂ + 1D duality | The linear version is done. The nonlinear version is open. |

## Resulting status of Christodoulou-SCC near Kerr

- **Near i⁺:** PROVED (sketch, repo) for extensions that are **causally plain** at the boundary point, modulo SH-a. SH-a is the nonlinear Price-law lower bound, i.e. item E3 here.
- **Remaining genuine gap:** exclude **bubbling** (non-causally-plain) Christodoulou-class extensions, i.e. Lemma B or a replacement. Per the warped escape, the argument must use the vacuum equations.
- **Whole CH:**
  - the corner closes via two-sided Hardy, given non-L² shear on both components near S_B;
  - this needs a both-ends blue-shift analysis (DL's unpublished bifurcation paper);
  - it also needs global C⁰ stability.

## Suggested next attack: Lemma B with vacuum

1. Feed the Riccati transfer identity ∇_L(Jχ̂) = −Jα (riccati_monotone) into the e₄-side cost.
2. The cheap sawtooth frame in §52 is unconstrained by curvature. In vacuum, its tilt modulation must be compatible with α's ℓ = 2 coherence (null_normal_form: Luk-class singularities are "slow" modes, ρχ̂ → 0).
3. Test: add the vacuum constraint (α determined by χ̂ via Riccati, not free) to the §52 optimal-control problem, and see whether the integrable-cost strategy survives.
