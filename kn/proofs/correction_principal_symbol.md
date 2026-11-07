# Correction: the DGS system does NOT have a scalar principal part

Source: Dias–Godazgar–Santos, eqs. (4), (7) and (8). Check script: `dgs_symbol.py`.

## The principal symbol

Fix a mode (ω, m) and take the principal symbol in (r, θ). Write X = √Δ ξ_r and Y = ξ_θ, and set

  p = Q² α₋ r̄*,  q = Q² α₊ r̄*.

The symbol is

    M(ξ) = [ (1+p)X² + (1−q)Y²      −(p+q)√Δ XY            ]
           [ −(p+q)XY/√Δ             (1−q)X² + (1+p)Y²      ]

- The coupling operators H₋₂ and H₋₁ are second order, so the off-diagonal entries are principal.
- det M = (1+p)(1−q)(X²+Y²)².
- The eigenvalues of M(ξ) are (1+p)(X²+Y²) and (1−q)(X²+Y²). Each is the scalar symbol times a function of position only. The eigenvectors, however, depend on the direction of ξ.
- 1+p = 3r̄(Mr̄ − Q²)·α₋ and 1−q = 3r̄(Mr̄ − Q²)·α₊. Both vanish only at r̄ = Q²/M, i.e. r = Q²/M on the equator. Since Q²/M < M < r₊, that point lies inside the horizon.
- So the reduced system is **elliptic** throughout the exterior, but its principal part is **not diagonal**. In the spacetime symbol it is also not the scalar wave operator times the identity.
- At large r, p and q are O(Q²/(Mr)) → 0. Near null infinity the principal coupling is therefore a perturbation of relative size O(1/r).

## Consequences for earlier claims

1. **Bridge Part 2, Step D.** Aronszajn–Cordes applies to diagonal principal parts, so it does not apply as stated. Possible fixes:
   - (a) Use a pseudodifferential decoupling. The eigenvalues are distinct (p ≠ −q generically) and the eigenprojectors are smooth and homogeneous of order 0, so the system decouples microlocally into two scalar elliptic operators. Calderón-type unique continuation should then apply. Points where p = −q need separate care.
   - (b) **Preferred:** use Giorgi's generalized Regge–Wheeler variables for Kerr–Newman instead (arXiv:2002.07228, JDG 2024). Check whether that system has a diagonal (scalar wave) principal part with lower-order coupling. If it does, Steps B–D go through as originally written.
2. **R7** is answered **negatively** for the DGS variables.
3. **Steps B and C (Rellich and Carleman at infinity).** The principal coupling is O(1/r) at infinity, so it is a second-order perturbation of relative size 1/r. A polynomial-weight Carleman estimate controls ∇u but not ∇²u directly. Either work in variables with a diagonal principal part (option b), or decouple first (option a).
4. **The real-axis reduction** (Rellich and unique continuation for the coupled system) has the same issue and the same fix.
5. **Unaffected.**
   - The numerics (near-zone weights, transition curve, balance tests) used the full DGS operators exactly.
   - Lemma A concerns null geodesics. The characteristic set is still the null cone, because the determinant is the square of the scalar symbol.
   - The Gap 2 trapping argument and Hintz's framework use the gauge-fixed metric and potential system, whose principal part is the wave operator times the identity.
