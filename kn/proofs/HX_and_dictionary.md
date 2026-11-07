# (H_X) and the DGS ↔ Giorgi dictionary

Source: Giorgi, arXiv:2002.07228, Remark 5.5, Prop 5.7 and its footnote 8 (the Newman–Penrose forms), §2.5, §4.1.

## 0. Convention fix

In the PDF text extraction, **complex-conjugation bars are lost**. Read literally, footnote 8's relation for x has the operator −δ + (3β + 3τ − π + α). Using the Kinnersley identity α = π − β̄ turns this into 3β − β̄ + 3τ. That combination has a cotθ coefficient (3/r̄ − 1/r̄*)/(2√2), which is complex and singular at the poles, and so cannot occur in a smooth GHP operator.

With the bars restored, (3β + 3τ − π̄ + ᾱ), and using ᾱ = π̄ − β, the operator becomes

  −δ + 2β + 3τ.

This is exactly the GHP form δ − (p − q)β + … for weights (p, q) = (3, 1), i.e. spin 1. Its cotθ coefficient is real.
**To verify:** check against Giorgi's tensorial relation (76) by computing (c)D̂⊗X + ½(3H + H̄)⊗X directly in her frame (§4.1).

## 1. The operator (Kinnersley frame, mode e^{-iωt+imφ})

The ingredients are

  δ = (1/(√2 r̄))(∂_θ − K_θ),  K_θ = m/sinθ − aω sinθ,  r̄ = r + ia cosθ,
  2β = cotθ/(√2 r̄),  3τ = −3ia sinθ/(√2 r̄ r̄*).

Hence

  𝒟_r x := −δx + (2β + 3τ)x = −(1/(√2 r̄)) [ ∂_θ − K_θ − cotθ + 3ia sinθ/r̄* ] x.

In Chandrasekhar's notation this is L†_{-1} plus the smooth term 3ia sinθ/r̄*. It maps spin weight 1 to spin weight 2.

## 2. (H_X): kernel and cokernel, with explicit solutions

**Kernel.** Solve f′ = (K_θ + cotθ − 3ia sinθ/r̄*) f. The solution is

  f = sinθ · tan^m(θ/2) · e^{aω cosθ} · (r − ia cosθ)^{−3}.

This is regular at both poles exactly when |m| ≤ 1, with the same pole exponents as at a = 0. So the kernel is **one-dimensional for |m| ≤ 1 and zero otherwise**, for every r > r₊, every ω, and every a, Q.

**Cokernel.** This is the kernel of the formal adjoint with respect to sinθ dθ, which is −L₂ plus a smooth term. Its solutions behave like θ^{−m−2} at θ = 0 and (π−θ)^{m−2} at θ = π. Both cannot be regular at once for any m, so the **cokernel is zero for all r, ω, m**.

⇒ **(H_X) holds**, conditional on the convention fix in §0. The right inverse 𝒢_r is uniformly bounded on compact r-intervals, because the coefficients depend smoothly on r and the kernel and cokernel have constant dimension. The leftover kernel coefficient k(r) appears only when |m| ≤ 1.

## 3. Dictionary (Remark 5.5 and its underlined versions)

| Giorgi | Newman–Penrose | DGS |
|---|---|---|
| A̲ | Ψ₄ | ψ₋₂ ∝ r̄*⁴ Ψ₄ |
| B̲ | 3φ₂Ψ₂ − 2φ₁Ψ₃ | ψ₋₁ ∝ r̄*³(2Φ₁Ψ₃ − 3Ψ₂Φ₂)/Φ₁ (same quantity, up to a nonvanishing factor) |
| F̲ | −δ̄φ₂ + …  (underlined f) | none |
| X̲ | Δφ₂ + …  (underlined x) | none |

**DGS pair = 0 ⟺ (A̲, B̲) = 0.** The other two then vanish as well:
- **F̲:** the underlined version of (75) reads P̄_F(∇₄A̲ + …) = ½D̂⊗B̲ + 3H̲⊗B̲ − (3P̄ + 2|P_F|²)F̲. With A̲ = B̲ = 0 this forces (3P̄ + 2|P_F|²)F̲ = 0. The coefficient 3Ψ₂ + 2φ₁φ̄₁ is nonzero on the exterior (compare the earlier factor 3r̄(Mr̄ − Q²), whose only zero, r = Q²/M, lies inside the horizon). **To verify:** the nonvanishing of exactly this combination.
- **X̲:** the underlined version of (77) gives (3P − 2|P_F|²)X̲ = 0, so X̲ = 0 provided that coefficient is nonzero. **To verify.**

So the bridge hypothesis "DGS pair vanishes" is equivalent to the whole negative-spin Giorgi set (A̲, F̲, B̲, X̲) vanishing. Bridge Part 2's unique-continuation argument runs on the positive-spin set (A, F, B, X), with X eliminated via §2.
