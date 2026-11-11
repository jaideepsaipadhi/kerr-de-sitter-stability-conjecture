# Large-m (semiclassical) reduction: lemmas so far

Notation. KdS roots r3 < 0 < r0 (inner) < r1 = r_e < r2 = r_c. Mobius map z: r0, r1, r2, r3 -> 0, 1, t, infinity.
mu(z) = z(z-1)(t-z).  omega = m w, lambda = m^2 ell, m -> infinity.
kappa_j = i Xi Khat(r_j)/Delta'(r_j), Khat = w(r^2+a^2) - a;  u_j = Im kappa_j;  sum_j u_j = 0.

Status tags: [P] proved by hand, [N] verified numerically only.

1. [N] Leading-order dual = original geodesic Heun equation with Euler parameter A1 = -2 kappa_3
   (checked to 1e-13).  [P] once the exponents at infinity, {0, -2 kappa_3} at this order, are
   established; this needs writing out.
   Dual exponents: 2(kappa_j + kappa_3), j = 0, 1, t.  ABt = 4 kappa_3^2.
2. [P] (spin decouples) The spin terms are O(m) and drop out of all leading-order conditions.
3. [P] ell enters only through an accessory shift: L2 = L2^0(z) + ell C0, with
   C0 = mu (dr/dz)^2 / Delta > 0 constant.  [N] checked: constant in z to 1e-14.
4. [N -> P once 5 holds] L2 = n0/z + n1 exactly (n2 = 0; checked to 1e-18).
   Reason: qt^2 - P2 is quartic and vanishes at z = 1, t.  n2 = 0 needs to be shown from the
   leading coefficient of P2 (the exponent at infinity) together with the linear qt.
5. [P] n0 = t (X0^2 - (X1+Xt)^2), X_j = |kappa_j + kappa_3|.  In the window, X0 - X1 - Xt = 2 u3, so
   n0 = 2 t u3 (X0+X1+Xt) > 0, because u3 > 0 for all w > Omega_3 and Omega_3 < Omega_c (|r3| > r_c).
6. [P] Lemma H: H(z) = 2mu - (z-z0) mu' > 0 on [1,t] for every z0 in (1,t).  Linear in z0; at z0 = 1,
   H = (z-1)[(z-1)^2 + (t-1)]; at z0 = t, H = (t-z)[t(t-1) - (t-z)^2].
7. [P] Lemma G: G(z) = H(z) - (z-1)(t-z)(z0-z) > 0 for 1 < z < z0 < t.  Linear in z0; at z0 = z, G = 2mu;
   at z0 = t, G = z(t-z)(t-1).
8. [P] Reduction.  With qt linear, f+ chosen between its two bounds, and f- likewise:
   cosmological side <=> L2'/L2 <= H/(mu (z-z0));  holds when n0 >= 0 (L2' <= 0).
   event side <=> -L2'/L2 <= H/(mu (z0-z));  with L2 = n0/z + n1, n0 > 0, n1 >= 0 this follows from G >= 0.
   Hence: LARGE-m CONDITIONS HOLD whenever n0 >= 0 [proved] and n1(ell) >= 0.
9. [P] n1 = ell C0 EXACTLY.  Expanding the dual at z = infinity gives n1 = -q2 (all kappa_3 terms cancel),
   and AB2 = 0 makes  -q2 = -mu(z)[Qz2(z) + P1(z)^2/4]  for every z.
   Key identity: the residues rho_j = Khat(r_j)/Delta'(r_j) sum to 0 over all four roots, so by partial
   fractions and the Mobius map  Xi Khat r'(z)/Delta = sum_{j=0,1,t} Xi rho_j/(z - z_j) = -i P1/2.
   Hence Xi^2 Khat^2 r'^2/Delta^2 = -P1^2/4: at ell = 0 the geodesic potential cancels the gauge term exactly.
   [N] confirmed: ell_crit = 0 to 1e-15 in six geometries.
   => L2(z) = n0/z + ell mu r'^2/Delta, with n0 = 2t u3 (X0+X1+Xt) > 0 and C0 = mu r'^2/Delta > 0.
   => ALL LEADING-ORDER LARGE-m CONDITIONS HOLD FOR EVERY ell > 0 (no angular bound needed at leading order).
10. [P] AB2 = 0, A1 = -2 kappa_3, n2 = 0.
    By the identity in item 9, Qz2 + P1^2/4 = ell C0/(z(z-1)(z-t)), so F(z) = AB2 z - q2 is the constant
    ell C0, which gives AB2 = 0 and q2 = -ell C0.  Then A1 is in {0, S1} with S1 = -2 kappa_3; the transform
    is the nonzero branch (the branch choice is part of the injectivity item).
    Dual: ABt2 = S1^2 = 4 kappa_3^2 and Pt1 ~ 4 kappa_3/z at infinity, so the 1/z^2 terms of Qt2 cancel,
    deg P2 <= 3, deg N <= 1, n2 = 0.
    => THE LEADING-ORDER LARGE-m STATEMENT IS PROVED, apart from the branch identification.
11. OPEN: next-order (O(1/m)) error control.  [N] step 24: the free-qt LP margin does not shrink with m
    (m = 4, 16, 64; tau ~ 1e-2).  Two cases at m = 64 / m = 4 have tau > 0 but fail the fine check, which is
    consistent with a seam boundary layer of width ~ 1/m^2 that the LP grid does not resolve.
    The crux: near z0, D's leading part q (L f)' vanishes, so the O(m) term (1/2) q1 (mu f')' decides.
    Flattening f there costs q L' f > 0 because L is decreasing, so the layer needs a matched multiplier,
    e.g. f'/f = -L'/L inside the layer.  This needs ell bounded below
    uniformly (ell >= c > 0), which is much weaker than the old ell_min requirement.

## Finite m: explicit multipliers (steps 25-26)
12. [Construction] With qt linear (exact finite-m horizon values) and L = (qt^2 - P4)/mu:
    cosmological side f+ = 1.  C+ = -q1 H <= 0 automatically (Lemma H); D+ = q L' <= 0 needs L' <= 0 on [z0,t].
      At large m, L' = -m^2 n0/z^2 + O(m), and n0 >= 2 t u3 X0 is bounded below on the closed window,
      so this holds for m >= m* (explicit once the O(m) term is bounded).
    event side f'/f = (1+delta) g (1 - exp(-kappa m^4 (z0-z)^2)) + p (1-u)/(z-1),
      g = -L'/L, kappa = L2(z0)/mu(z0), u = (z-1)/(z0-1).
      The Gaussian factor is the exact solution of the layer inequality (1/2) mu psi' >= sigma L2 (1-psi),
      sigma = m^2 (z0-z); it gives f'(z0) = 0, so the seam holds with equality.  The p-term makes f(1) = 0;
      near z = 1 it needs p <~ 2 (z0-1) L(1)/mu'(1).
13. [N] Fine-grid check (geometric refinement to 1e-14 at z0) passes for ALL 280 tested cases:
    m = 1, 2, 4, 16, 64; 7 geometries including a = 0.999 and Lambda = 1e-3; window fractions 0.02-0.98;
    ell in {ell_min, 5}; with (delta, p) from a list of five.  Margins at m = 8...512 are stable
    (C ~ -0.1...-0.33, D ~ -0.16...-0.8 relative).
14. OPEN for a full proof of the real-frequency superradiant case:
    (a) explicit m* with rigorous O(1/m) bounds (all leading margins are now closed-form positive);
    (b) certification for m < m* on the compact parameter set (interval arithmetic / Bernstein);
    (c) a finite-m eigenvalue lower bound (ell bounded below) and L > 0;
    (d) the finite-m dual data are currently a numerical polynomial fit (residual ~1e-11) and should be
        made symbolic.

## Spin +2 and the upper half-plane (steps 28, literature)
15. [N] s = +2 works directly (no Teukolsky-Starobinsky transfer needed).  The real dual is on the other Euler branch
    (gauge (0,0,0) with B instead of A); step3.real_branch now picks the branch automatically.
    The dual window coincides with CTdC (3.37) to 10 digits (3 geometries), and the explicit layer construction
    closes 27/27 cases (m = 1, 4, 32; window fractions 0.05/0.5/0.95; ell = ell_min).
    Leading order: spin drops out (item 2), so the large-m lemmas carry over unchanged.
16. [Reduction, to write] The upper half-plane follows by continuity from Schwarzschild-de Sitter, as in Hintz 2608.25868.
    HPV25 Thm 2.2: for each b, finitely many resonances with Im sigma >= 0, each finite-dimensional; the HV
    framework (high-energy estimates via normally hyperbolic trapping, Dyatlov) holds for all b in B (Remark 2.12);
    solution operators depend continuously on b (Thms 2.11, 2.13).  A mode in Im sigma > 0 can only appear by
    crossing the real axis or coming from infinity; the second is excluded locally uniformly by the high-energy
    estimates.  So Assumption 1.2 reduces to: (real sigma != 0) + (sigma = 0) + (SdS: known, Kodama-Ishibashi /
    Hintz-Vasy), with the continuity argument run for the gauged operator whose non-decaying modes are the
    ungauged ones (HPV constraint damping, Thm 2.7).

## The transform: contour, exponents, injectivity (step 29)
17. [P, modulo write-up] Let y be a QNM in the Heun gauge (0,0,0): exponent rho = 1 - delta at z = 1, analytic at z = t.
    Transform: ytilde = c * I_1^{mu+1}[y], mu = A - 2 (A = the real-dual root; B for s = +2).  This is the
    Riemann-Liouville integral based at z = 1, Pochhammer-regularised; a contour around an analytic point gives 0,
    so the base point must be 1.
    (a) ytilde ~ B(mu+1, rho+1) (z-1)^{A-delta} at z = 1, which is the dual's non-analytic branch 1 - delta_tilde.
        The Beta coefficient is nonzero.
    (b) ytilde is analytic at z = t: write ytilde = (z-1)^{mu+1} int_0^1 (1-tau)^mu y(1+(z-1)tau) dtau.  y is analytic
        near the segment (1, t], and 0 and infinity are far away, so this is analytic in z near t.
    (c) Injective: I^{mu+1} is inverted by I^{-mu-1} (Riemann-Liouville semigroup; Titchmarsh uniqueness), and
        c != 0 because mu and rho are not integers.  Im rho = -2 Im kappa_1 != 0 and Im mu = Im A != 0 on the open window.
    [N] Termwise check: the transformed Frobenius series solves the dual Heun equation (residual ~1e-17), its exponent
        equals 1 - delta_tilde exactly, |B| ~ 0.3-1.5; checked for s = -2 (two geometries) and s = +2.
    => A nonzero real-frequency QNM in the window gives a nonzero dual solution with exactly the boundary behaviour
       used in the two-current argument.  Since the dual argument forces that solution to vanish, no such QNM exists.

## Exact finite-m structure and the 4-parameter reduction (steps 30-32)
18. [N, fit residual 1e-12 to 1e-15; to be derived] At EVERY finite m, for s = -2 and s = +2:
      L(z) = A/z + B - 3z,   A = m^2 n0^(2),   B = m^2 ell C0 + c1,   c1 = 1 + t + 3(1 - alpha) C0,  alpha = Lam a^2/3.
    That is: n0 has no lower-order terms, n2 = -3 exactly, and c1 is geometric (independent of w, ell, s;
    checked in 4 geometries).  C0 = (t-1) r'(1)/Delta'(r1) (evaluate mu r'^2/Delta as z -> 1).
19. [P given 18] COSMOLOGICAL SIDE, ALL m: L' = -A/z^2 - 3 < 0 on [z0, t].  With f+ = 1, C+ = -q1 H <= 0
    (Lemma H) and D+ = q L' <= 0.  No error term.  Still needed: L(t) > 0, i.e. A/t + B - 3t > 0, which holds
    because B >= c1 > 3t - 1 - t + ... (check: c1 - 3t = 1 - 2t + 3(1-alpha)C0 > 0 when C0 > (2t-1)/(3(1-alpha))).
20. [Reduction] EVENT SIDE: q1 cancels, so the conditions depend only on (t, z0, A, B):
      s h <= H/mu,   (1/2)(mu' h + mu(h' + h^2)) <= s (L h + L'),   s = z0 - z,  h = f'/f.
    [N] With the layer multiplier h = (1+delta) g (1 - exp(-beta kappa s^2)) + p(1-u)/(z-1), kappa = L(z0)/mu(z0):
        649/666 grid points close; all 17 failures have B <= 3.5, far below the physical floor B >= c1 ~ 1 + 3.7t.
        Everything at B = 10(1+t) and B = 1e4 closes, for t in [1.05, 265], z0 anywhere, A in [1e-3, 1e5].
    => Error control and certification reduce to ONE lemma in four parameters (t, z0, A, B), plus the
       physical constraint B >= c1(t, C0, alpha), with no reference to KdS beyond these numbers.
21. TODO: derive item 18 analytically (the item-9 identity should extend: the spin terms give the -3z and c1);
    prove C0 > (2t-1)/(3(1-alpha)) (it gives L(t) > 0); then prove or certify the 4-parameter event lemma on
    {B >= c1}.  The m -> infinity corner (A, B -> infinity) is the proved leading-order statement.

## Exact exponents; the event lemma on the physical region (steps 33-35)
22. [N to 1e-11; P once written] Exact finite-m dual data: the real-dual root is A = 1 + s - 2 kappa_3, and the other
    root is B_inf = 1 + 2s (= -3 for s = -2, 5 for s = +2).  Hence the dual exponents at z = 0, 1, t are exactly
    2(kappa_j + kappa_3): spin cancels, so n0 = t(X0^2 - (X1+Xt)^2) exactly at every m.  The exponent difference at
    infinity is B_inf - 1 = +-4, which gives n2 = -3.  c1 = 1 + t + 3(1-alpha) C0 [N to 1e-13, 12 random geometries];
    this is the accessory datum and is still to derive.
23. [N] C0 / ((2t-1)/(3(1-alpha))) >= 1.50 on 12 random geometries (-> 1.5 in the Kerr limit, i.e. C0 ~ t - 1/2).
    Hence L(t) > 0 with margin, and c1 >= 4t - 1/2.
24. [N] Event lemma on B >= 4t - 1/2: 2366/2376 points close (t in [1.01, 300], z0-fraction in [0.001, 0.999],
    A in [1e-4, 1e6]) with p = c (z0-1) L(1)/mu'(1) adaptive.  The 10 failures: z0-fraction 0.999 (next to the
    cosmological horizon), t >= 40, B exactly at the floor.  All close at B >= 1.3 (4t - 1/2).
    Physically B = c1 + lambda C0, so B >= 1.3 (4t - 1/2) follows from lambda >~ 1.2 (C0 ~ t).
    => NEEDED: the eigenvalue bound lambda >= 1.2 for the spin-2 KdS angular operator (a = 0: lambda = l(l+1) - 4 >= 2).

## C0 bound and the angular eigenvalue (steps 36-37)
25. [P] Item 1 (exponents at infinity): R ~ (r-r3)^x with x in {kappa3, -kappa3-s}, and r - r3 ~ 1/z.  The gauge
    factors contribute z^{-kappa3+2s+1}, so A = x - kappa3 + 2s + 1, giving {1+2s, 1+s-2kappa3}.
26. [P] C0 = 3/(Lam (r1-r0)(r2-r3)) and t = (r1-r3)(r2-r0)/((r1-r0)(r2-r3)), from the Mobius constants
    (r - r_j = (r_j - r3)(z - z_j)/(zinf - z), r - r3 = (r0 - r3) zinf/(zinf - z)).
    With 3(1-alpha)/Lam = -e2 (Vieta), C0 > (2t-1)/(3(1-alpha)) <=> F > 0, where
      F = 4 r0^2 + 7 r0 r1 + r0 r2 + 4 r1^2 + r1 r2 + r2^2   (every term positive).  PROVED.
    => L(t) > 0 at every m, and c1 >= 4t - 1/2.
27. [N] Lowest angular eigenvalue (STU eq. 3.1 converted to our convention: lam_ours = lam98 - s(1-alpha); Chebyshev
    collocation with regular endpoints; reproduces l(l+1) - 4 at a = 0): min over the window, m <= 6, s = +-2,
    8 geometries: lam_min >= 2.000 (as low as 2.0002 for Lam = 1e-4, a = 0.99).
    Required (step 35): lam >= lam_req = (1.3(4t - 1/2) - c1)/C0 <= 1.21 (-> 1.2 in the Kerr limit).
    Smallest margin: 0.80.  TODO: prove lam >= 1.2 (an energy estimate for the angular operator over the window).

## c1 derived (step 38)
28. [P up to a final symbolic simplification; N to 1e-16 at 5 points, both spins, frequencies outside the window too]
    Closed-form chain with no fitting:
      n1 = ABt(1+t) - qt + Tt(1 - St/2) + (7/2)(1+t)   (from P4 = p4 z^4 + p3 z^3 + ..., n1 = p3 + p4(1+t)),
      q  = AB(1+t) + T(1 - S/2) - c3,  c3 = 2 zinf a_{-2} - k a_{-1}  (Laurent data of the r-normal form at r3;
           the Mobius Schwarzian vanishes),
      A = 1+s-2kappa3, B_inf = 1+2s, dual parameters as in step 2.
    Result: n1 = lambda C0 + 1 + t + 3(1-alpha) C0 for every omega (not only in the window), every s = +-2.
    n2 = ABt + St/2 - St^2/4 + 3/4 = -3 exactly (confirmed in the same chain).
    Remaining: reduce the identity n1 - (lambda C0 + 1 + t + 3(1-alpha)C0) = 0 symbolically, modulo the Vieta
    relations (e1 = 0, e2 = -(3/Lam)(1-alpha), e4 = -3a^2/Lam).  This is mechanical.
=> Every ingredient of L(z) = A/z + B - 3z is now derived:
     A = m^2 n0^(2) = t(X0^2 - (X1+Xt)^2) > 0 (proved),   B = lambda C0 + 1 + t + 3(1-alpha)C0,
     C0 = 3/(Lam (r1-r0)(r2-r3)) (proved).

## Certification build (cert/)
29. [P] Sharper C0 bound: 2(1-alpha) C0 >= 2t - 1  <=>  G = 3 (r0 + r1)^2 >= 0.  Hence c1 >= 4t - 1/2 and
    lambda C0 >= lambda (t - 1/2).  With lambda >= 1.5 (numerics: >= 2.000), the physical floor is
    B >= 5.5 t - 1.25 (proved, given lambda >= 1.5).
30. [Verifier] cert/verifier.py (points) and cert/boxes2.py (boxes) prove Cc < 0 and x*Dd < 0 on [1, z0]
    with mpmath interval arithmetic.  Ingredients:
    - cancellation-free forms: x h, and x^2(h' + h^2) = x^2(g psi)' - p + (xh)^2, using d/dz[p(1-u)/x] = -p/x^2;
    - dependency-safe forms, using mu' = mu/x - x(2z - t):
        C  = s[(xh - 1) mu/x + x(2z - t)] - 2x mu/x,
        xD = (1/2)(mu/x (xh + X2) - x(2z - t) xh) - s(L xh + x L');
    - a center-fixed multiplier in tau = (z-1)/X0 (x h and x^2 h' are then center quantities, by the chain rule);
    - (A,B) via the 4 box corners (C does not depend on (A,B); xD is linear in (A,B)).
    Endpoint signs are automatic: C(z=1) = X0(t-1)(p-1) and xD(1) = (1/2)(t-1)p^2 - X0 L(1) p.
    [P, computer-assisted] RIGOROUSLY CERTIFIED so far:
      t in [2,3] x zf in [0.2,0.8] x A in [0.5,50] x Bf in [1,3]   (109 boxes, 29 s)
      and 8 tiles with t in [1.2,1.6], zf in [0.05,0.7], A in [0.01,100], Bf in [1,10].
    Bottleneck: zf -> 1 (z0 next to the cosmological horizon; mu(z0) -> 0, so the layer is sharp).  Tiles there are
    mostly certified but hit the box cap.  Next: tau-grading near tau = 1 and a layer-adapted variable; analytic
    treatment of the corners zf -> 0, 1, t -> 1, infinity, and A, B -> infinity.

31. [P] c1 identity, symbolic (step39_c1_symbolic.py). With roots r0,r1,r2 free, r3=-(r0+r1+r2),
    a^2 and alpha eliminated by Vieta, and omega, m, lam, a, s free symbols: the numerators of
    n1 - (lam C0 + 1 + t + 3(1-alpha) C0) and n2 + 3 are divisible by s^2 - 4 as polynomials.
    So the identity is exact for s = +-2 (and false for general s; it is spin-specific).
    Laurent coefficients at r3 taken in closed form; no series, no fitting.

32. [P modulo finite check] Angular bound lam_ours >= 1.5 on the superradiant window (steps 40-41).
    Exact split (verified to 1e-13): U = U_std + Sing_al + E, with Sing_al >= 0 because its numerator is
    >= m^2+s^2+2smx >= (|m|-2)^2, and p = (1+al x^2)(1-x^2) >= 1-x^2.  Hence H >= H0 + E in the form sense,
    where H0 is the a=0 spin-weighted operator (levels l(l+1)-4, l >= max(|m|,2)).
    |m| >= 3: lam_std = m^2+|m|-4 >= 8, inf E >= -C al |m| -> trivially >= 1.5.
    |m| = 2: pointwise, lam >= 2 + inf E >= 1.65.
    |m| = 1: pointwise fails (inf E ~ -1.0 near extremal a), but the two-level gap bound
      lam >= mineig[[2+<E>, sigma],[sigma, 8+inf E]] gives >= 1.87.
    Window check: CTdC window subset of omega/m in [0, Omega_H] on all 4011 sampled geometries.
    Remaining finite piece: <E>, Var(E), inf E are explicit elementary functions of (al, xi) on a compact
    box, so the inequalities are a 2-parameter interval check (deferred with the other finite certifications).

33. [Reduction + N] Event-lemma CORNERS (step42, step42b).  Margins along sequences into each corner converge to
    strictly negative limits: A,B->inf (large m) -0.50; t->1 -0.76; t->inf -0.36; zf->1 -0.50; A->0 -0.445.
    zf->0 (seam at the event horizon = window edge q~(1)->0) closes only with the ADAPTIVE p = c 2(z0-1)L(1)/(t-1),
    and degenerately: margin ~ -K zf.  Corner lemma there = the rescaled problem x = (z0-1) y, divided by z0-1.
    TODO (math): write each limit problem and prove it; they are 1-3 parameter problems with explicit limits.

34. [P modulo convention check] OUTSIDE THE WINDOW (step43).  Real dual + conserved Wronskian W[R, Rbar]; the dual
    mode is the 1-delta_t branch at z=1 and the analytic branch at z=t; W(1) = W(t) gives
    Im e_1 |c_1|^2 = t Im e_t |c_t|^2,  e_j = 2(kappa_j + kappa_3).  So no mode where the signs differ.
    [N] The signs AGREE exactly on the CTdC window and DIFFER everywhere outside it: 249288 samples, no exception.
    So the two-current argument is needed only inside the window, and outside it a single identity suffices.
    Window endpoints (Im e = 0) excluded: measure zero, handled by continuity (item 36).

35. [P] SIGMA = 0.  m != 0: omega = 0 lies outside the window (243/243 geometries), so item 34 applies
    (Im kappa_1, Im kappa_3 != 0 since am != 0, so the transform is non-degenerate).
    m = 0, s = -2: (Delta^{-1}R')' + V Delta^{-2} R = 0, mode ~ Delta^2 at both horizons, V = -2 Lam r^2 - 2(1-al) - lam < 0
      => int Delta^{-1} R'^2 = int V Delta^{-2} R^2 <= 0 => R = 0.
    m = 0, s = +2: [P, symbolic] Delta^2 R_{+2} solves the s=-2 equation (same lam); [P, symbolic] the K=0 TS map
      f -> Delta^2 f'''' preserves the s=-2 equation, and sends the regular branch to the Delta^2 branch, so the
      s=-2 energy identity kills it; kernel = cubic solutions, excluded since the 4x4 determinant is
      [lam^2 + 6(1-al) lam + 8 - 4 al + 8 al^2]^2 > 0 for lam > 0 (= ((lam+2)(lam+4))^2 at al = 0).
    Convention to confirm: future-regular modes are the Delta^{-s} branch at BOTH horizons.

36. [Argument] UPPER HALF-PLANE by continuity.  Teukolsky QNMs on KdS form a discrete set depending continuously on
    (Lam, a) (Fredholm analytic family; Vasy 2013, Petersen-Vasy JEMS 2024 for the full subextremal range), with no
    modes for Im sigma > C uniformly on compact parameter sets (high-frequency/trapping estimates).  At a = 0
    (SdS) there are none in Im sigma >= 0, sigma != 0 (known).  For fixed Lam the subextremal a-range is an
    interval containing 0.  A mode reaching Im sigma > 0 at some a would cross the real axis, but items 23, 34 and 35
    exclude every real sigma.  (The window endpoints need one line: there Im e_j = 0 and a real QNM would give a
    nonzero solution with zero flux; still to close, either directly or because QNMs cannot sit only at these points.)

37. [OPEN, load-bearing] THEOREM B.  HPV25 is conditional on mode stability for the UNGAUGED LINEARIZED EINSTEIN
    equation, not Teukolsky alone (Haefner-Hintz-Vasy 2506.21183, intro).  In Kerr, AHW24 derives it from Teukolsky
    mode stability (Whiting), which is its decisive input.  So the port needs: (B1) a KdS Wald-type lemma:
    psi_0 = psi_4 = 0 for a mode with sigma != 0 => pure gauge; at sigma = 0, => linearized KdS + gauge;
    (B2) the AHW Hertz-potential / TS reconstruction steps with Lam terms; (B3) the sigma = 0 zero-mode classification.

38. [Reduction] THEOREM B, step-by-step against AHW (arXiv 2207.12952), with what changes for Lam > 0:
    (i)   Thm 4.5 (no Teukolsky modes, Im sigma >= 0): = our items 23/34/35/36.  KdS has no infinity; the decay /
          indicial-root step at r = inf is replaced by the cosmological horizon (already symmetric in our setup).
          STILL TO PORT: vanishing beyond the horizons (AHW's energy functional E(r), vanishing to all orders at r+),
          now needed at both r+ and r_c.
    (ii)  Prop 6.5 (psi_0 = psi_4 = 0 => locally Plebanski-Demianski + gauge): rests on completeness of the
          Aksteiner-Backdahl invariants (AABKW, CMP 384, 2021), proved only on KERR.  Khavkine's theorem (1805.03751):
          the KdS Killing operator is of regular finite type, so a complete compatibility complex EXISTS abstractly.
          The gap is identifying its first operator with {psi_0, psi_4, linearized Ricci - Lam h, I_xi, I_zeta}_Lam.
          => THE research item of Theorem B.  Plan: build the Lam-terms into the AB invariants (GHP on KdS) and check
          completeness by Khavkine's jet-rank criterion at a point (a finite linear-algebra computation).
    (iii) Lemma 6.6 (local -> global): uses only simple connectivity of the exterior.  KdS static region
          R_t x (r+, r_c) x S^2 is simply connected: carries over.
    (iv)  Prop 6.8 (decay at infinity forces ndot = cdot = 0, i.e. no NUT and no acceleration): no infinity in KdS.
          Replacement: smoothness at the axis.  A linearized NUT gives a Misner string (g_{t phi} ~ n(cos th +- 1)
          cannot be regular at both poles), and a linearized acceleration gives a conical defect ~ alpha M at one
          pole (Griffiths-Podolsky).  Both are excluded for smooth modes.  To verify by direct linearization of the
          Lam-PD metric.
    (v)   sigma != 0: AHW's argument (g = (i/sigma) delta^* omega from L_xi g = delta^* omega + ...) uses only
          stationarity: carries over verbatim once (ii) holds.
    (vi)  sigma = 0: also uses Thm 4.7 (scalar operator invertible at sigma = 0) and the 1-form Coulomb step, both
          proved at infinity in Kerr.  On KdS: Vasy-framework Fredholm theory plus the two-horizon setup; to check.

39. [N, strong] B1 at the jet level (step44_jet_completeness.py).  At p = (r=3, th=1), a = 0.6, M = 1, using a TPS +
    dual-number linearisation (all checks pass: tetrad null/orthogonal 1e-14, background Lam-vacuum, background psi0 = 0,
    gauge terms lie in the solution set to 1e-16):
      d_J = dim{jets of linearised Lam-vacuum h with dpsi_0 = dpsi_4 = 0, extendable} - dim{gauge jets}
      K=6: {3:4, 4:6 (top-level artifact)};  K=7: {3:4, 4:4, 5:4};  K=8: {3:4, 4:4, 5:4, 6:4}
    IDENTICAL for Lam = 0 (Kerr, where AABKW/AHW give exactly the 4-dim PD tangent m, a, NUT, acceleration) and
    Lam = 0.05.  So locally the solution space of {linearised Einstein, psi_0 = psi_4 = 0} mod gauge has rank 4 on KdS,
    as on Kerr.  The Lam-PD family supplies 4 tangent directions, so B1 (local) holds once those are shown independent
    mod gauge (TODO: plug the 4 PD tangents into the jet map).
    To make rigorous: exact rational arithmetic plus a formal-integrability/stabilisation bound (Khavkine finite type),
    which is a finite computation and goes on the deferred-certification list.  Also repeat at more points.
(i) [P, separated modes] Vanishing beyond the horizons: the radial ODE has regular singular points at r+ and r_c with
    finite exponents, so a smooth solution vanishing on (r+, r_c) vanishes to infinite order at r+ and r_c and hence
    identically on the extended domain.  General (non-separated) modes reduce to this via the Finster-Smoller
    Jordan-block argument used by AHW.

40. [N, strong] B1 IDENTIFIED (step45_pd_tangents.py).  In the PD chart (tau, r, p, sigma) with
    P = k + 2np - eps p^2 + 2 alpha m p^3 - (alpha^2 k + Lam/3) p^4,  Q = k - 2mr + eps r^2 - 2 alpha n r^3 - (alpha^2 k + Lam/3) r^4,
    Omega = 1 - alpha p r, the background (n = alpha = 0) passes Ric = Lam g (2e-13) and psi0 = 0, and the five tangents
    d g/d(m, n, eps, k, alpha) satisfy all constraints (1e-12 to 1e-15).  Their span modulo gauge has rank 4 = d_J
    (J = 3, 4).  So the local space {linearised Lam-vacuum, dpsi_0 = dpsi_4 = 0}/gauge IS the PD tangent space:
    Lam-analogue of AHW Prop 6.5 at the jet level (exact-arithmetic + stabilisation goes on the deferred list).

41. [P] AXIS REGULARITY replaces AHW Prop 6.8 (decay at infinity).  Background n = alpha = 0: P is even, with roots +-p0
    (the axes), P'(p0) = -c < 0, P''(+-p0) = d.  At a root p_j the axis Killing field is eta_j = c_j (p_j^2 d_tau + d_sigma),
    and elementary flatness gives c_j = 2/|P'(p_j)| (the conformal factor cancels).  A smooth metric (no periodic time)
    needs eta_+ = eta_-:  (F1) p_+^2 = p_-^2  and  (F2) |P'(p_+)| = |P'(p_-)|.  These are geometric, so a smooth gauge
    term does not change their linearisations.  Even perturbations (dm, deps, dk) preserve both.  The odd ones
    q_n = 2p dn and q_alpha = 2 m p^3 d alpha shift both roots by the same dp = -q(p0)/P'(p0), and:
      dF1 = 8 p0^2/c dn + 8 m p0^4/c d alpha,
      dF2 = -2(2 d p0/c + 2) dn - 2(2 d m p0^3/c + 6 m p0^2) d alpha,
      det = -64 m p0^4/c != 0   (m != 0, a != 0).
    So a smooth global linearised solution with dpsi_0 = dpsi_4 = 0 has dn = d alpha = 0.  The remaining even directions,
    modulo the scaling gauge, are exactly dM and da.  (a = 0 is Schwarzschild-de Sitter, known.)  Globalisation:
    AHW Lemma 6.6 with a simply connected exterior.
=> THEOREM B, sigma != 0 and the PD/axis part: done, modulo the deferred exact-arithmetic jet certification.
   Remaining for Theorem B: (vi) the sigma = 0 extras (AHW Thm 4.7 analogue and the 1-form Coulomb step on KdS).

42. [P, assembling 34-41] THEOREM B: ungauged Einstein mode stability on subextremal KdS (the Hintz-Vasy UEMS form,
    which is what HPV25 assume).  Let h be a generalized mode, h = e^{-i sigma t} sum_{j<=J} t^j h_j, Im sigma >= 0,
    solving Ric'(h) = Lam h, smooth on the extended domain.
    (1) Teukolsky scalars vanish.  dpsi_0 and dpsi_4 are generalized Teukolsky modes.  The top coefficient is a genuine
        mode, so it is 0 (items 23/34/35/36 and the vanishing beyond the horizons).  The next coefficient then solves the
        homogeneous equation (its source is the top one), and so on down: all vanish.
    (2) Local structure (item 40) holds for ANY linearised solution, not only modes: h = sum c_i dg/dp_i + delta^* X
        locally; the coefficients are constant (local system) and the decomposition is global (simple connectivity).
    (3) Axis regularity (item 41): only the dM and da directions survive.
    (4) sigma != 0: apply (L_T + i sigma)^{J+1}.  It kills h and is invertible on the stationary PD part, giving
        c_M dg/dM + c_a dg/da = delta^*(...).  These are not pure gauge (rank 4 mod gauge, item 40), so c = 0
        and h is pure gauge.
        sigma = 0: h = c_M dg/dM + c_a dg/da + delta^* X, with X allowed polynomial in t: exactly UEMS.
    The gauge-FIXED zero-energy analysis (AHW Thm 4.7, the Coulomb 1-form step) is NOT needed for the ungauged statement.
    Caveat: confirm against HPV25's precise statement once public.

STATUS (end of programme pass): every logical step from "Teukolsky mode stability" through "UEMS on KdS" is in place.
Outstanding, all finite or check-type:
  (a) event-lemma interior tiles [interval cert] + corner limit problems (item 33, math)
  (b) angular bound, 2-parameter interval check (item 32)
  (c) exact-arithmetic jet ranks + Khavkine stabilisation bound, more points (items 39/40)
  (d) window endpoints Im e_j = 0 (item 36)
  (e) conventions: future-regular = Delta^{-s} branch at both horizons (items 34/35)
  (f) angular-mode completeness for complex sigma (non-self-adjoint spheroidal; Finster-Smoller-type, cited as in AHW)

43. [P] (e) CONVENTIONS settled.  Future-regular boundary conditions (regular in ingoing v at r+, outgoing u at r_c):
    R ~ (r-r1)^{-kappa_1 - s}   (Delta^{-s} e^{-ik r*})   at the event horizon,
    R ~ (r-r2)^{+kappa_2}       (e^{+ik r*}, NO Delta^{-s}) at the cosmological horizon  (matches CTdC).
    The transform pipeline (Heun gauge (0,0,0): exponent 1-delta at z=1, analytic at z=t) encodes exactly this, so items
    23/34 stand.
    CORRECTION to item 35 (m = 0, sigma = 0): I had wrongly used Delta^2 at BOTH horizons for s = -2.  Correct statement and
    proof, using only the two exact identities from item 35:
      s = -2: R ~ Delta^2 (analytic) at r+, analytic at r_c.  g := Delta^2 R'''' again solves the s=-2 equation (TS identity),
              and g ~ Delta^2 at BOTH ends, so the energy identity (V < 0) gives g = 0, hence R is cubic, hence R = 0 (det > 0).
      s = +2: f := Delta^2 R_{+2} solves the s=-2 equation, is analytic at r+ and ~ Delta^2 at r_c; the same TS map kills it.
    So sigma = 0 holds with the corrected BCs; the proof is uniform in both spins.

44. [P + N] (d) WINDOW ENDPOINTS.  Im e_1 and Im e_t are affine in omega, so each vanishes at exactly one omega.  [N] At all 994
    sampled endpoints exactly one of them vanishes (ratio <= 9e-14) and the other is nonzero; Im kappa_1, Im kappa_3 != 0 there,
    so the transform is non-degenerate.
    Im e_t = 0, Im e_1 != 0: the Wronskian gives c_1 = 0, so the dual is analytic at z=1; but the transform forces a nonzero
      (Beta-coefficient) non-analytic branch unless y = 0.
    Im e_1 = 0, Im e_t != 0: c_t = 0, so the dual vanishes near t, hence identically (Frobenius uniqueness), hence y = 0 (injectivity).
    Endpoints closed.  Item 36's caveat removed.

45. [BUG FOUND AND FIXED] valid_kds imposed Lam M^2 < 1/9 (the a = 0 bound) and accepted complex roots.  Subextremal KdS
    reaches Lam M^2 ~ 0.177 and alpha -> 7 - 4 sqrt 3 ~ 0.0718.  Fixed (four distinct real roots, r3 < 0 < r0 < r1 < r2).
    step46_resweep.py reran every cheap check on the corrected region (5069 geometries): window subset [0, Omega_H]: OK;
    omega = 0 outside the window: OK; outside-window sign identity: 0 failures; endpoints: one Im e vanishes (ratio 7e-11);
    xi/m <= 0.441; angular margins OK.  The 4-parameter event lemma is uniform in (t, zf, A, B >= floor), and the floor
    is analytic, so it is unaffected.

46. [P, computer-assisted, ~1 min] ANGULAR BOUND lam >= 3/2 on the window, ALL subextremal KdS (replaces item 32).
    Hypotheses, all proved:
      alpha = Lam a^2/3 <= (2 - sqrt3)^2 = 7 - 4 sqrt3 on subextremal KdS (known bound; consistent with the max 0.0715 found).
      omega/m >= 0 on the window: lo = a(|D'0| - |D'1|)/((r1^2+a^2)|D'0| - (r0^2+a^2)|D'1|) with
        |D'0| - |D'1| = (Lam/3)(r1-r0)[f(r0) - f(r1)], f(x) = (x - r3)(r2 - x), a concave parabola with vertex
        -(r0+r1)/2 < 0, so f is decreasing on [r0, r1] and lo > 0.  Hence m xi >= 0.
      xi/m = a omega/m <= 1: hi is a positive-weight average of a/(r2^2+a^2) and a/(r0^2+a^2), each < 1/a.
    |m| >= 3: analytic.  With m xi >= 0 the xi-terms are >= 0, E >= -6 al - 4 al(1+al)|m|, so
      lam >= m^2 + |m| - 4 - 0.43 - 0.31|m| >= 6.6.
    |m| = 1, 2 (step47_angular_proof.py): H >= H0 + P + R, P polynomial, |R| <= al C; two-level (gap) bound with EXACT
      Beta moments <P>, <P^2>; interval arithmetic over al in [0, 0.0718], xi in [0, |m|].  Certified minima:
      1.957 (m=1, s=-2), 1.665 (m=1, s=+2), 1.920 (m=2, s=-2), 1.914 (m=2, s=+2).  All >= 3/2.
    Two-level validity: psi0 is the unique l=2 state of H0 for fixed (m, s), the next level is 8, and
    <H phi, phi> >= (2+p)|c|^2 + (8+Pmin)|phi|^2 - 2 sigma |c||phi|.

47. [P, analytic] EVENT LEMMA via an ODE multiplier, in two regimes (replaces interval certification there).
    With h = f'/f: D <=> (mu f')' <= 2s (Lf)',  C <=> s mu f' <= H f,  s = z0 - z, H = 2mu + s mu'.
    Construction: f solves (mu f')' = 2s (Ltilde f)' with f(1) = 0, i.e.  mu f' = 2s Ltilde f + 2 int_1^z Ltilde f,
    which forces f ~ (z-1)^{ptilde}, ptilde = 2(z0-1) Ltilde(1)/(t-1).
    (R1) p := 2(z0-1)L(1)/(t-1) <= 1: take Ltilde = L.
    (R2) p > 1 and 2(z0-1)^2 (A/z0 + 3) < t-1: take Ltilde = L - kappa, kappa = L(1) - (t-1)/(2(z0-1)) > 0 (ptilde = 1).
      Ltilde(z0) > 0 is exactly the R2 condition, since L(1) - L(z0) = A x0/z0 + 3 x0.
    Proof.  Ltilde > 0, so mu f' > 0: f is positive and increasing, and the seam condition f'(z0) > 0 holds.
      D: (mu f')' = 2s(Ltilde f)' <= 2s(Lf)' because ((L - Ltilde) f)' = kappa f' >= 0.  In R1 this is an equality.
      C: Ltilde is decreasing, so s mu f' <= 2s Ltilde(1) f (s + z - 1) = 2 s x0 Ltilde(1) f <= s(t-1) f.  It remains
         to show H >= s(t-1).  With z = 1+x: H - s(t-1) = x[2(1+x)(t-1-x) + s(2t-4-3x)].  The bracket is linear in s,
         >= 0 at s = 0, and at s = t-1-x it equals (t-1-x)(2t-2-x) >= 0.  QED
    [N] Sanity: the R1 sufficient inequality holds on 3407/3407 random samples.
    REMAINING REGIME: p > 1 and 2(z0-1)^2 (A/z0+3) >= t-1, i.e. L drops by more than (t-1)/(2x0) across [1, z0].
    That is the large-A (large-m, WKB) regime, where the seam boundary layer is essential.  It is covered numerically by
    the layer multiplier (item 20) and asymptotically by the proved leading-order large-m lemma; it needs a uniform
    analytic bridge.

48. [P] (f) NOT NEEDED.  Separation of variables is only used at REAL sigma (items 23, 34, 35, 44, and the Teukolsky input
    to item 42 at real sigma).  There xi = a omega is real, the angular operator is self-adjoint with compact resolvent,
    its eigenfunctions are complete and orthonormal, and every (generalised) mode decomposes exactly.  Im sigma > 0 is
    handled by continuity on the FULL (unseparated) operator (item 36), which needs only the real axis to be mode-free.
    No Keldysh / Finster-Smoller completeness for non-self-adjoint spheroidal operators is required.

49. [Partial] UNIVERSAL ODE MULTIPLIER (step48_universal.py).  Ltilde = c constant gives the EXPLICIT solution
    f' = (1/mu) exp(int 2c(z0-eta)/mu), which depends only on (t, z0, c).  D reduces to ((L-c)f)' >= 0, and A and B
    eliminate exactly: sup over A >= 0 and B >= floor of |L'|/(L-c) is max(1/z, 3/(5.5t - 1.25 - 3z - c)).  So the event lemma
    would follow from ONE 2-parameter inequality in (t, z0), with C checked exactly as 2cs(s + int f/f) <= H.
    [N] It FAILS for t >~ 3.5, mid zf, for a structural reason: the ODE forces (mu f')' = 2csf' >= 0, so D needs (L-c)f
    increasing right up to the seam (the s cancels).  The seam-layer mechanism needs (mu f')' < 0 near z0.
    So: use this explicit multiplier on [1, z1] (it handles the z -> 1 zone automatically), switch to the seam layer on
    [z1, z0], and match at z1.  That is a TWO-zone proof, not three.  Next step.

50. [N, negative] TWO-ZONE HYBRID (step49_hybrid.py): explicit ODE multiplier on [1, z1] plus the A-scaled seam layer on
    [z1, z0] closes only 258/400 random R3-type samples.  Failures concentrate at large t (>~ 15) and large B, and some are
    numerical (gradient/NaN at the seam).  Structural lesson, which is the useful output: for A -> inf, D forces h >= g
    (WKB) everywhere outside a layer of width -> 0 as A grows.  The constant-c ODE profile Psi does not track g = |L'|/L
    when B dominates, so zone I must itself be WKB-like there.  The step-32 layer multiplier (h = (1+delta) g phi + p(1-u)/x)
    remains the right object.
    Plan for an ANALYTIC proof of R3: prove the layer multiplier's two inequalities directly, with
      (i)   bulk: s delta |L'| dominates the O(1) left side once delta, beta are fixed;
      (ii)  seam layer: rescale sigma = sqrt(kappa) s; the limit problem is m-independent (already proved, step 13-15
            leading order), and the corrections are O(kappa^{-1/2}) with explicit constants;
      (iii) z -> 1 zone: the R2 argument applies locally, since L drops little over a short interval.
    A careful multi-page estimate; not closed in this pass.

51. [N, strong] R3 IS CLOSED BY ONE UNIVERSAL LAYER MULTIPLIER (step50, step50b).
    h = 2 g (1 - exp(-2 kappa s^2)) + p (1-u)/x,  g = -L'/L,  kappa = L(z0)/mu(z0),  p = (1/2) min(0.9, 2 x0 L(1)/(t-1))
    (delta = 1, beta = 2, c = 1/2) closes ALL 1260 R3 samples (t in [1.01, 3e3], zf in [0.05, 0.95], A from the R2
    threshold up to 1e4 times it, B from floor to 1e3 times floor), with worst normalised margin 0.23.
    EVENT LEMMA, DECOMPOSED:  R1 (p <= 1) analytic;  R2 analytic;  R3 one fixed explicit multiplier with margin >= 0.23.
    Remaining for R3:
      (a) bounded part: interval certification with this FIXED multiplier (no parameter search; the margin is uniform).
      (b) unbounded directions (A, B -> inf; t -> inf; zf -> 0, 1): the scaling limits of a FIXED multiplier are explicit
          functions, so each corner is a closed-form limit inequality.  Item 33 already showed the limits are negative.
    Layer derivation (for the analytic writeup): with h = G phi, G = (1+delta) g, phi = 1 - E, E = exp(-beta kappa s^2), D is
      mu G beta kappa s E >= s L g (1 - (1+delta)phi) + (1/2)phi (mu'G + mu G') + (1/2) mu G^2 phi^2.
      Leading balance at the seam: (1+delta) beta E >= (1+delta) E - delta (true for beta >= 1, margin delta).
      Bulk: s delta |L'| >= (1/2)(mu'G + mu G' + mu G^2), where the right side is scale-free.
      C: Lemma G with delta capped by min G z/(s mu).

53. [N + partial P] CORNERS with the FIXED R3 multiplier (step51_corners_fixed.py).  Every corner sequence stays
    negative and converges: B->inf -0.379, A,B->inf -0.335, A->inf -0.334, t->inf -0.284, zf->0 (inside R3, which forces
    A ~ x0^-2) -0.379, zf->1 -0.372, t->1 -0.379.
    [P] The recurring -0.3793 is exact: the C condition at z -> 1, where h ~ p/x and p = 0.45, gives (p-1)/(p+1) = -0.55/1.45.
    [P] B -> inf: g -> 0, so h -> p(1-u)/x and D's right side s L h -> +inf off the seam.  At the seam the left side tends to
        -(1/2) mu(z0) p/x0^2 < 0, uniformly in B.  C is B-independent in the limit and is covered by Lemma G (p <= 0.45 < 1).
    TODO (analytic): write the A -> inf / t -> inf limit problems.  For A -> inf this is the proved large-m leading-order
        lemma, plus O(A^{-1/2}) layer corrections.

54. AUDIT NOTES (hostile pass, short).  Dependencies to re-verify before believing the whole:
    - item 23 (two-current window argument) was validated numerically only under the OLD valid_kds region; rerun on the full
      subextremal region (steps 26/28 layer sweeps) to be safe;
    - item 17 transform injectivity: write the Riemann-Liouville / Titchmarsh argument out in full;
    - items 39-40: floating-point ranks only; exact certification plus Khavkine prolongation order still OPEN;
    - item 42: HPV25's exact hypothesis is unverified (paper not public);
    - the cosmological-side lemma (item 19) and the C0 floor (item 36) are analytic: OK.

55. [P] TRANSFORM INJECTIVITY, written out (replaces the sketch in item 17(c)).  Real omega makes every kappa_j purely
    imaginary.  In Heun gauge (0,0,0), y = (z-1)^rho sum_k a_k (z-1)^k near z = 1, with rho = 1 - delta = -2 kappa_1 - s, and
    mu = A - 2 = s - 1 - 2 kappa_3.  The Riemann-Liouville / Pochhammer integral I^{mu+1} based at z = 1 acts termwise:
       (z-1)^{rho+k}  ->  Gamma(rho+k+1)/Gamma(rho+k+mu+2) (z-1)^{rho+k+mu+1},   with rho + k + mu + 2 = k + 1 - e_1,
    e_1 = 2(kappa_1 + kappa_3).  Gamma(rho+k+1) is finite (rho non-real: Im rho = -2 Im kappa_1 != 0 on the window), and
    1/Gamma(k+1-e_1) != 0 (k+1-e_1 is non-real inside the window and equals k+1 at the endpoint e_1 = 0).  So every
    termwise factor is finite and nonzero, and I^{mu+1} y = 0 forces all a_k = 0.  Injective on the CLOSED window.
    The leading coefficient B(mu+1, rho+1) = Gamma(mu+1)Gamma(rho+1)/Gamma(1-e_1) is nonzero for the same reasons.
    No Titchmarsh argument is needed.

56. [Reduction] REMAINING CORNERS.
    t -> inf (zf fixed, A = t a, B = t b, b >= 5.5): with z = t zeta, mu ~ t^3 zeta^2 (1-zeta) and L ~ t(b - 3 zeta) + a/zeta;
      the WKB slope in zeta tends to 3/(b - 3 zeta), and a drops out.  The z -> 1 end tends to the exact (p-1)/(p+1) margin.
      So the corner is a 2-parameter (zf, b) explicit limit inequality, numerically -0.2835.
    A -> inf: the proved large-m leading-order lemma (steps 13-15) plus O(A^{-1/2}) seam-layer corrections; needs the
      constants written out.

57. [OPEN, precise statement] B1 FINITE-TYPE STEP (Khavkine).  Computed: at each point p and each jet order J = 3..6,
    pi_J(solutions) = pi_J(gauge) + span(PD tangents).  Needed: a solution h vanishing modulo gauge to order J0 at p is
    pure gauge near p.  Equivalently, the gauge-invariant content of {linearised Lam-vacuum, dpsi_0 = dpsi_4 = 0} is a
    FLAT CONNECTION on a rank-4 bundle, determined by its J0-jet.  Route: Khavkine 1805.03751 Thm 9 (regular finite type
    Killing operator => complete compatibility complex) gives completeness for SOME K1.  Identifying ker of
    {dpsi_0, dpsi_4} within K1 needs the KdS analogue of AABKW's identification, or a direct proof that the system above
    closes (involutive) at an explicit order N, checked by verifying that the order-(N+1) prolongation adds no new
    conditions.  That verification is again a finite rank computation, but it must be EXACT and paired with Cartan-Kuranishi.

58. [Reduction, P modulo exact ranks] B1 CLOSURE WITHOUT AABKW-TYPE INVARIANTS.  Let S_J(p) be the TRUE space of J-jets at p
    of local solutions of {linearised Lam-vacuum, dpsi_0 = dpsi_4 = 0}, and G_J(p) the J-jets of gauge terms.  Put
    d_J(p) = dim S_J/G_J.
    (a) Upper bounds are computable: S_J ⊂ pi_J(S_K) for every K (truncated jets over-count), so
        dim pi_J(S_K) - dim G_J >= d_J.  Lower bound: the 4 PD tangents are actual solutions, independent mod gauge
        (item 40), so d_J >= 4.
    (b) If d_{J0}(p) = d_{J0+1}(p) = 4, the projection S_{J0+1}/G -> S_{J0}/G is onto between spaces of equal dimension,
        hence injective: the gauge-reduced symbol vanishes at order J0+1, and by prolongation at every higher order.
        So the gauge-invariant data is a finite-type system, a flat connection on a rank-4 bundle of J0-jets, and smooth
        solutions are determined by the J0-jet at a point (ODE uniqueness along paths).
    (c) Then a local solution h has J0-jet in G_{J0} + span(PD) at p, i.e. h - delta^*X - v vanishes mod gauge to order J0
        at p, and so is pure gauge near p (Khavkine finite-type Killing operator: the gauge normal form is smooth).
    Floating point gives d_3 = d_4 = 4 at the test point (items 39/40), so J0 = 3.
    What is still owed: (i) EXACT ranks for dim pi_3(S_7), dim pi_4(S_7) and G_3, G_4 at a rational point (cos th = 3/5,
        r = 3, Lam = 1/20, a = 3/5: all TPS coefficients are rational, so FLINT fmpz_mat ranks are exact);
      (ii) every point, not only a generic one: the ranks are analytic in (r, th), so the bound holds off an analytic
        subset Z; on Z the dimension can only jump UP.  Close it via (b) applied on the complement plus continuity of the
        constant PD coefficients across Z, or by exact ranks over Q(r, cos th) symbolically.

59. [Design note] EXACT-RANK STEP for B1 (item 58(i)): what is needed, and why the shortcuts fail.
    Write the constraint matrix C (rows = unknown jet coefficients), split into low (degree <= J) rows C_L and high rows C_H.
      dim pi_J(S_K) = n_L - rank C + rank C_H.
    So the UPPER bound on d_J that the closure needs requires a lower bound on rank C (mod-p ranks give this) AND an UPPER
    bound on rank C_H.  Upper bounds on rank are never certifiable from mod-p ranks or from verified floating point:
    both only certify LOWER bounds.  An exact null space of C_H over Q is unavoidable.
    Shortcuts considered and rejected:
      - many-prime CRT against the Hadamard bound: thousands of primes, each needing a full rebuild;
      - verified numerics: certifies rank lower bounds only;
      - semicontinuity in Lam from Kerr (AABKW at Lam = 0): rank C_H moves in the wrong direction;
      - explicit gauge elements of S ∩ ker pi_J: they miss the non-extendable top-level jets.
    Feasible route: build C EXACTLY via symbolic coefficients.  dRic and dpsi_{0,4} are second-order operators whose
    coefficients are explicit rational functions of (r, cos th, sin th); take their exact 2-D Taylor series at the rational
    point (r = 3, cos th = 3/5, sin th = 4/5, Lam = 1/20, a = 3/5), assemble C (and the gauge map) over Q, and run FLINT
    exact rank / null space for K = 7, J = 3, 4.  Estimated a few hours of engineering plus compute.  Not done in this pass.

60. [Key simplification] Theorem B only needs B1 for (generalised) MODES, h = e^{-i sigma t + i m phi} H(r, x), x = cos th.
    Jets then live in 2 variables (r, x) instead of 4, and the metric is RATIONAL in (r, x) (no trig).  At K = 7 that is about
    10 x 36 = 360 unknowns instead of 3300: an exact computation is small.  Owed: (a) exact Q(i) linear algebra at fixed
    (sigma, m); (b) uniformity in (sigma, m).  The ranks drop only on an algebraic subset of (sigma, m), so handle the
    generic case by exact symbolic rank over Q(i)(sigma) per m, and the exceptional sigma explicitly.  (c) Large |m|
    needs an m-uniform argument, since m enters polynomially; again symbolic in m.
    Also observed: the "full column rank of C_HH" shortcut fails because the linearised contracted Bianchi identity
    makes the equations dependent, so exact null vectors really are needed.

61. [P at sample points, EXACT] B1 FOR MODES (step52_exact_modes.py).  Exact Q(i) arithmetic (python-flint); mode-reduced 2-D
    jets at r = 3, x = cos th = 3/5, Lam = 1/20, a = 3/5, M = 1.  Exact upper bounds on d_J; truncation artifacts sit only
    at the top two levels J = K-2, K-1 (they move up with K).
      sigma = 1/3 + i/5, m = 1:  K=7 -> d_J = 0 for J <= 4;   K=9 -> d_J = 0 for J <= 6.
      sigma in {1/2, i/4} (m=0), 5/2 + i/9 (m=4), 1/10 (m=-2), 3 + 2i (m=7), 1/1000 and i/1000 (m=1): d_0 = d_1 = 0
      (and d_2 = 0 at K=6).
    With d_{J0} = d_{J0+1} = 0 (J0 = 0), item 58(b) applies on a rank-0 bundle: at these (sigma, m) and this point, every
    mode with dpsi_0 = dpsi_4 = 0 is LOCALLY PURE GAUGE.  This is rigorous: exact ranks, and upper bounds are what (a)
    needs.
    OWED: (i) uniformity over all sigma in the closed upper half-plane, all m, and all points.  The ranks are algebraic in
        (sigma, m, r, x), so this needs symbolic exact rank or explicit exceptional-set analysis.
      (ii) *** sigma = 0 ANOMALY ***: at sigma = 0 the exact counts are stable in K and grow like 2J+1
        (m=1: K=7 and K=9 give 0,3,5,7,9,11,...; m=0, with the 4 extra t- and phi-linear Killing gauges added:
        0,2,5,7,9,11).  Taken at face value that is infinitely many local stationary non-gauge solutions, which
        CONTRADICTS the 4-D count of 4 (which reproduced Kerr's known answer).  sigma = 1/1000 is clean, so sigma = 0
        is an isolated exceptional point of the mode-reduced system.  Either (A) a sigma = 0-specific bug or incomplete
        gauge in the mode code, or (B) the floating 4-D count under-counted (float rank can be inflated by noise).
        MUST be resolved: sigma = 0 is exactly where linearised KdS (dM, da) lives in item 42.
        Diagnostic plan: check that dg/dM lies in S at sigma = 0; compute the exact 4-D count at the rational point;
        list a basis of the excess at J = 1 and test each against Lie derivatives along non-mode vector fields.

62. [Status] R3 core certification (cert/r3core.py, 1 core, 60 s per-tile cap), partial: 69/192 tiles done; 61 fully
    certified, 8 with 4-8 leftover boxes each.  All 8 are at t < 2.5 with large A or zf near 1, consistent with the
    per-tile time cap rather than failures.

63. [Estimate] HPV25 HYPOTHESIS MATCH.  Public descriptions (HHV 2506.21183; the 2606.28253 survey) say HPV assume mode stability
    for the UNGAUGED linearised Einstein equation, which is the Hintz-Vasy UEMS notion of HV18.  UEMS: for Im sigma >= 0, a
    generalised mode solution of Ric'(h) = Lam h is pure gauge (sigma != 0), or linearised KdS plus pure gauge with a
    gauge 1-form allowed to grow polynomially in t (sigma = 0).  Item 42 gives exactly this, given B1 (items 40, 58, 61)
    and Theorem A.  Likely mismatches and how much each would cost:
      - HPV may phrase it on their own Fredholm spaces (variable-order Sobolev at the horizons), not "smooth modes".  Our
        modes are smooth across both horizons (items 43/(i)), so we expect only a regularity bridge (a few lemmas).
      - They may require the result for the GAUGE-FIXED operator (as in HHV for Kerr).  That needs, in addition, the 1-form
        (constraint damping) mode analysis and the gauge-fixed zero-energy space (HHV25 + Hintz 2024b analogues on KdS).
        This is a real extra project, though Hintz-Vasy did it for slowly rotating KdS.
      - The sigma = 0 statement must match their zero-energy description: it depends on resolving item 61(ii).
    Best guess: a ~70% chance the ungauged UEMS form is exactly what is needed, ~30% that a gauge-fixed or 1-form
    supplement is required.
    61(ii) diagnostic, first step: h = dg/dM satisfies ALL 120 exact sigma = 0 constraints (step52b_check_dM.py), so the
    sigma = 0 constraint matrix is consistent.  The anomaly therefore lies in the gauge side or the counting, not in wrong
    equations.  Next: exact 4-D count at the rational point (decides between (A) and (B)), then a basis of the J=1 excess.

64. [P, EXACT] sigma = 0 ANOMALY RESOLVED.  Cause: the mode code imposed only psi_0(h) = dC(h)(l,m,l,m) = 0 and
    psi_4(h) = dC(h)(n,mb,n,mb) = 0 on the COMPLEX mode h.  A real perturbation's psi_0 = 0 is two real conditions, and
    complexified that is FOUR complex conditions: dC(l,m,l,m) = dC(l,mb,l,mb) = dC(n,mb,n,mb) = dC(n,m,n,m) = 0.
    Theorem B is entitled to all four: if h is a mode with Im sigma >= 0, then hbar is a mode with frequency -conj(sigma),
    which is also in the closed upper half-plane, so Theorem A applied to hbar gives the two conjugate-tetrad conditions.
    With all four imposed (step52, CONJ=1), K = 7, exact:
      sigma = 0, m = 1:          d_J = 0 for ALL J <= 5 (even the top-level artifacts are gone);
      sigma = 0, m = 0:          d_J = 0, 2, 4, 4, 4, 4  -> 4 = the PD tangents (dM, da, NUT, acceleration), matching 4-D;
      sigma = 1/3 + i/5, m = 1:  d_J = 0 for all J <= 5.
    Earlier sigma != 0 results (item 61) remain valid: fewer constraints give LARGER upper bounds, so their d = 0 was
    already a valid bound.
    Diagnostics recorded for the referee: the 4-D lift test (step52d) is invalid for the psi blocks (complex psi_0 of the
    complex mode does not vanish on Re h and Im h separately), and its Ricci "violations" were roundoff normalised by a
    ~1e-16 vector (step52k).  Entrywise, the mode operator equals the 4-D one exactly (steps 52e-52j).
    B1 for modes, closure: sigma = 0, m = 0 gives d_2 = d_3 = 4, every other case d_0 = d_1 = 0, so item 58(b) applies.
    STILL OWED: uniformity over all (sigma, m) and all points (exact symbolic ranks / exceptional sets).

65. [P exact at samples + reduction] UNIFORMITY of B1 for modes.
    (a) With all four conditions imposed, K = 4 already gives d_0 = d_1 = 0 (the sigma = 0, m = 0 PD sector needs K >= 5).
    (b) step53_sample.py: 30 random exact samples (Lam M^2 <= 0.09, a in [0.12, 0.98], points across (r+, r_c) x (-1, 1),
        sigma real/imaginary/complex/zero, |m| <= 4), all d_0 = d_1 = 0.
    (c) [P] Points: only ONE good point per mode is needed.  If d_0 = d_1 = 0 at some p for (sigma, m), closure (item 58(b))
        makes h pure gauge near p.  The counts are analytic in the point, so the good set is open and dense.  Khavkine
        1805.03751 (the KdS Killing operator is of regular finite type) gives an abstractly complete compatibility
        operator K1; K1(h) = 0 on a dense set, hence everywhere, so h is locally pure gauge everywhere, and globally pure
        gauge by simple connectivity (AHW Lemma 6.6).
    (d) OWED: for each m, existence of ONE good point for EVERY sigma in the closed upper half-plane.  At a fixed point,
        the constraint and gauge matrices are polynomials in sigma of degree <= 2, so C(sigma) = C0 + C1 sigma + C2 sigma^2
        can be built exactly by interpolating 3 builds.  The bad sigma are the roots of explicit determinant polynomials
        (random compressions), a finite set that can be checked exactly at that point or at a second point.  Then large |m|:
        m enters polynomially, so treat m the same way; an m-uniform bound is needed for |m| -> inf.

66. [P mod-P certificate, m = 1] sigma-UNIFORMITY of B1 for modes (step54, step55_exceptional.py).  K = 4, J in {0, 1}.
    Build C(sigma), G(sigma) EXACTLY as degree-<=2 polynomials (3 exact builds).  At p1 = (r=3, x=3/5) and
    p2 = (r=7/2, x=-1/5), Lam = 1/20, a = 3/5, m = 1:
      rank C = 66 (exceptional polynomial delta_C of degree 11); generic rank C_H over Q(i)(sigma) EXACT (sympy):
      66 (J=0), 60 (J=1); rank G_L = 10, 24.  Generic d bound n_L - r_C + R_H - r_G = 0 for J = 0 and J = 1 at both
      points.  So d_J(sigma) <= 0 at every sigma that is not exceptional at that point.
    Exceptional polynomials E1 = delta_C delta_G0 delta_G1 (at p1) and E2 (at p2): gcd(E1, E2) mod P (P = 998244353 = 1 mod 4,
    i = sqrt(-1) mod P) is EXACTLY c * sigma^4.  So the only possible common exceptional sigma is sigma = 0, and sigma = 0
    is checked directly (exact, K = 7, all four conditions, d = 0; item 64).
    => for m = 1, every sigma in C has a good point; with item 65(c), every m = 1 mode with all four Teukolsky
       conditions is pure gauge.
    Caveat for full rigour: a nonzero common root over Qbar could reduce mod a prime above P to 0 (excluded after deflating
    sigma^4, since the deflated polynomials are coprime mod P) or to infinity (needs exact leading coefficients nonzero
    mod P; corroborate with a second prime, or compute delta exactly by multimodular CRT).
    Other m: same script (python3 step55_exceptional.py <m>, about 10 min per m).  |m| -> inf needs an m-uniform argument
    (m enters polynomially; treat m as a second symbol, or use a large-|m| asymptotic rank argument).
    m = 2: identical outcome (generic d bounds 0 at both points; gcd(E1, E2) = c sigma^4 exactly); sigma = 0, m = 2 is
    exact-checked (K=4: d_0 = d_1 = d_2 = 0).  CERTIFIED (same caveat).  The script's last line prints False only because
    it does not special-case sigma = 0; the exact sigma = 0 check closes it.

## 67. Prime rigor for the σ-uniformity certificates (item 6) — CLOSED for m = 1, 2

Lemma (reduction). E1, E2 ∈ Z[i][σ], 𝔭 | P the prime with i ↦ I_MOD. If (a) 𝔭 ∤ lc(E1), (b) 𝔭 ∤ lowest nonzero
coefficient of E1, (c) gcd(E1 mod 𝔭, E2 mod 𝔭) = u σ^c, then every common complex root of E1, E2 is σ = 0.
Proof: primitive gcd h has lc(h) | lc(E1), so its roots are 𝔭-integral and ≡ 0 mod 𝔭; write E1 = σ^b F with F(0) a
𝔭-unit by (b); a root σ0 ≠ 0 would give F(σ0) = 0 yet F(σ0) ≡ F(0) ≢ 0. ∎
E1 = δ_C δ_G0 δ_G1 (exact minors, exceptional set ⊂ their zeros).
 - δ_G0, δ_G1: degree and order pinned exactly by the tropical (assignment) bound = mod-P values (step57_lc_check.py).
 - δ_C: tropical bound is loose (45 vs 11), so δ_C computed EXACTLY over Z[i] by multimodular CRT with proved
   Hadamard bound (log2 B ≈ 9.1k / 11.9k bits; 147 / 192 primes; ~2 min per point; step57b_exact_deltaC.py).
   Result at both points, m = 1 and m = 2: exact ord_0 = 4, exact degree = 11, both extreme coefficients 𝔭-units.
 - (c) is step55's gcd = c σ^4.
Hence for m = 1, 2: no σ ≠ 0 is exceptional at both test points — EXACT, no probabilistic or mod-P caveat.
The "root at infinity" flag of step55 is a non-issue: it only reflected the degree bound 2r exceeding the true degree.
Method scales to any fixed m at ~5 min per m; it is the finite-m part of the (σ, m) argument (items 4/5).

## 68. Joint (σ, m) resultant test (items 4/5) — fast, but single-prime certificate FAILS; status

step58_resultant_fast.py: E_i(σ, μ) interpolated once as a bivariate polynomial mod p (40×40 grid, spot-checked);
σ-order 6, deg_σ 29, deg_μ 29 at both points. R(μ) = formal Res_σ(E1/σ⁶, E2/σ⁶): degree 1135 (bound 1334), not ≡ 0.
~90 s per prime (was ~1 h).
Result at 4 primes (998244353, 1073741789, 1073741741, 1073741717): R mod p has 51 / 59 / 47 / 59 roots in F_p,
in ± pairs (m ↔ −m symmetry), always including μ = 0. Only integer root with |m| ≤ 2000 at every prime: m = 0.
No common small-height rational root across primes.
Reading: ~50 roots per prime (vs ~1 for a generic polynomial) ⇒ R splits into ~50 irreducible factors over Q(i)
(Chebotarev). So "R mod p has no F_p root" will NEVER hold for any prime — that certificate route is dead.
Heuristic: an integer m ≠ 0 surviving all 4 primes has density ~ Π(#roots/p) ≈ 2⁻⁹⁷. Not a proof.
Rigorous routes still open:
 (A) third test point: g_p = gcd(R12, R13) mod p. If g_p = cμ^k then any integer common m has m ≡ 0 mod p, and
     with the exact μ-order of R12 at 0 (one exact slice computation at m = 0) the lemma of §67 kills m ≠ 0.
 (B) lift caveat for both routes: exact divisibility σ⁶ | E_i (bivariate). Fixed-m exact data (§67) shows order 4 for
     δ_C, so the σ⁶ must be split between δ_C and the gauge minors; needs an exact argument (structural σ-factor
     of gauge rows from ∂_t ↦ −iσ) rather than computation.

## 69. Joint (σ, m) uniformity — reduced to one handed-off certificate

Three test points (r, x) = (3, 3/5), (7/2, −1/5), (4, 1/3) at Λ = 1/20, a = 3/5, M = 1.
(i) σ⁶ | E_i EXACTLY in Q(i)[σ, μ] (step58f): block-Toeplitz deficiency at σ = 0, T₂ of each minor, rank bounded at
    2n+1 integer μ ⇒ all relevant minors vanish identically. δ_C contributes σ⁴, G₂(J=2) σ², G₁ σ⁰. Lemma used:
    ord_σ det A ≥ K r − rank T_K(A₀, A₁, A₂) for every K; transfers to Q(i)(μ) since minors have μ-degree ≤ 2Kr.
(ii) Exact degree bounds deg_σ, deg_μ E_i ≤ 80, 78 (step58g, reversed T₂) ⇒ the 81×81 mod-p interpolation is exact.
(iii) On that exact grid (step58h): deg_σ = deg_μ = 29, R_ij(μ) = formal Res_σ(F_i, F_j), F = E/σ⁶, degree 1135 with
     ord_{μ=0} = 357, and gcd(R12, R13, R23) mod p = μ³⁵⁷ EXACTLY.
Consequence (rigorous once (iv) holds): if an integer m with 0 < |m| < p had a σ ≠ 0 exceptional at all three points,
then m mod p ≠ 0 would be a root of μ³⁵⁷. Impossible. Repeating at k primes: no such m with |m| < ∏ p (≥ 2^120 for 4).
(iv) NEEDED: exact deg_σ E_i = 29 (so formal degrees are honest; excludes exceptional σ with a 𝔭-denominator).
     At μ₀ = 1, point 1: exact reversed T₉ gives deg ≤ 11 + 6 + 12 = 29 ✓. Uniform-in-μ version = run_degcert.py
     (62 + 8 + 11 μ₀ values × 3 points, ~50 min on 2 cores) — HANDED OFF.
Then the (σ, m) question = large-|m| asymptotics with any explicit M₀ < 2^120, plus m = 0 (already exact).

## 70. Item 4 (large |m|) — NO ASYMPTOTICS NEEDED; reduced to an exact local statement at the origin

Structure found (step59, mod p, exact 81-grid):  E_i(σ, μ) = σ⁶ μ⁶ H_i(σ, μ) at all three test points.
 - Res_σ of the F_i = E_i/σ⁶ therefore carries a spurious μ^(6·46) = μ²⁷⁶; for the H_i: R'_ij = Res_σ(H_i, H_j) has degree
   859, ord_{μ=0} = 81, and gcd(R'12, R'13, R'23) mod p = μ⁸¹ EXACTLY.
 - H_i(σ, 0) = σ⁹ · (distinct factors); the common root is σ* = 0, i.e. the ORIGIN (σ, m) = (0, 0), and 81 = 9·9 is the
   intersection multiplicity of two curves with 9-fold points there.
Argument (all integer m ≠ 0, no size bound):  let a = exact ord_{μ=0} R'12 ≤ 81.  If a = 81, write R'12 = μ⁸¹ T with T(0)
a 𝔭-unit.  An integer m ≠ 0 with a σ ≠ 0 exceptional at all three points has m mod p a root of μ⁸¹, so p | m, and then
T(m) ≡ T(0) ≢ 0 mod 𝔭, contradicting T(m) = 0.  So the whole (σ, m) question reduces to a = 81, and
  a ≥ I_(0,0)(H1, H2) ≥ mult_0 H1 · mult_0 H2,   mult_0 H = mult_0 E − 12,
so it suffices that mult_0 E_i ≥ 21 EXACTLY.  Also needed: exact μ⁶ | E_i (same Toeplitz method as σ⁶; mod-p evidence
only so far) and the degree certificate of §69 (now for H, same degrees in σ).
Exact progress (step59b; lower bound on mult_0 via T_K deficiency along EVERY line, interpolated over the slope):
   G1: 2, G2: 4 at all three points (K = 2, saturated).  C: 11 (K=2), 13 (K=3) at point 1; mod-p scan saturates at 15
   from K = 4, so exact K = 4 should give 15 + 2 + 4 = 21. Exact K = 4 for C is ~1 h per point -> run_origin_cert.sh.
STATUS: item 4 = run_degcert.py + run_origin_cert.sh + exact μ⁶ divisibility.  No analytic large-m lemma required.

## 71. CORRECTION to item 51 (R3 multiplier) + μ⁶ exact

(a) μ⁶ | E_i EXACTLY at all three points (step60_mu6.py, forward T₂ in μ, interpolated over σ): all six from δ_C.
    Item 4 now waits only on run_degcert.py + run_origin_cert.sh.
(b) ITEM 51 WAS WRONG.  The fixed multiplier h = 2g(1 − e^{−2κs²}) + p(1−u)/x FAILS condition C for t large and A ≫ B
    (A/B ≥ 100 at t ≥ 100: margin +0.005 … +0.08).  Never sampled before (corner sequences had A ≈ B).
    Mechanism (zone 1 ≪ z ≪ t, ψ = 1+δ, g ≈ 1/z):  C needs z(z−1)h ≤ 2z − 1, i.e. (1+δ)(z−1) + p·z(z−1)w(z) ≤ 2z − 1.
    With w = (1−u)/x this needs p z ≤ 1 (fails for z > 1/p); with δ = 1 the g-term alone saturates up to O(1/z).
    REPAIR:  h = (1+δ) g (1 − e^{−2κs²}) + p (1−u)/(x z),  δ = 0.7  (same p, β = 2).
      Near z = 1 unchanged (exact (p−1)/(p+1) C-margin kept); large-z C condition becomes (1+δ)(z−1) + p ≤ 2z−1, margin ~ (1−δ)z.
    step61_retune.py, 1200 cases (t ∈ [1.01, 1e5], zf ∈ [0.05, 0.95], B/floor ∈ {1, 10, 1e3}, A/B ∈ [0.01, 1e6]): no failures;
      worst C −0.082 (t = 1e5, A/B = 1e6), worst D −0.033 (t = 1.01, zf = 0.05: the known degenerate zf → 0 corner, margin ~ −0.65 zf,
      item 33 rescaling).  δ = 0.5: C −0.14 but D −0.026; δ = 0.85: C −0.040.
    CONSEQUENCES: items 51/53 margins and the parked R3 tile certification must be redone with the repaired multiplier;
    the corner analysis (item 5) now has a new corner A/B → ∞ (with t → ∞), where the inner-zone C limit is the explicit
    inequality (1+δ)(z−1) + p ≤ 2z − 1 (margin (1−δ)z − δ + 1 − p ... > 0 for z ≥ 1 when δ < 1, p < 1).
    Inner-zone D limit (t → ∞, z = O(1), A = ta, B = tb): leading ratio −(a/z² + X)/(3a/z² + X) ≤ −1/3, X = (a/z+b) p w > 0.

## 72. ITEM 23 PROVED (analytic; replaces the 12-geometry numerics; independent of the valid_kds bug)

C0 = 3 / (Λ (r1 − r0)(r2 − r3))  [identified; matches ncoef to 1e-8 on 14 geometries, ω-independent — step: inline check].
With Δ_r = −(Λ/3) ∏(r − r_i):  Σ r_i = 0,  e2 = −3(1−α)/Λ,  t = (r1−r3)(r2−r0)/((r1−r0)(r2−r3)).
Claim C0 > (2t−1)/(3(1−α))  ⟺  (3/2)Σ r_i² − 2(r1−r3)(r2−r0) + (r1−r0)(r2−r3) > 0, which after r3 = −(r0+r1+r2) is
   f = 4r0² + 7r0r1 + r0r2 + 4r1² + r1r2 + r2²,  positive definite (leading minors 4, 15/4, 7/2).   ∎
Sharp form:  C0 − 1.5 (2t−1)/(3(1−α)) = 3 (r0 + r1)² / (2 Λ (r1 − r0)(r2 − r3)) ≥ 0,
   i.e. C0 / bound ≥ 3/2 EXACTLY, equality iff r0 = −r1 — explaining the numerical 1.50 → 1.5 (Kerr limit).
Holds for every Λ > 0 with four real ordered roots: no subextremality restriction enters.
Corollary: c1 = 1 + t + 3(1−α) C0 ≥ 4t − 1/2.  The c1 identity itself (item 22) is still numerical (1e-8 here,
1e-13 before): closed form 3(1−α)C0 = 9(1−α)/(Λ(r1−r0)(r2−r3)); deriving it from the Heun accessory parameter is a
bookkeeping exercise left open.

## 73. t -> inf corner for the REPAIRED multiplier (δ = 0.7): inner + outer limit problems

Scaling A = t²α (α ≥ 0 covers both A ~ tB and A/B -> inf), B = t b (b ≥ 5.5), zf -> ζ0.
OUTER zone z = tζ (step62_outer_limit.py): every term of C is Θ(t³), every term of D is Θ(t); the p-term is O(t⁻²)
  and drops.  With Λ = α/ζ + b − 3ζ, g∞ = (α/ζ² + 3)/Λ, ψ = 1.7(1 − e^{−2K}), K = Λ(ζ0)(ζ0−ζ)²/(ζ0²(1−ζ0)), G = ψ g∞,
  m = ζ²(1−ζ):
     C∞ = (ζ0−ζ) G m − 2m − (ζ0−ζ) m',
     D∞ = ½(m' G + m(G' + G²)) + (ζ0−ζ)(α/ζ² + 3)(1 − ψ).
  Normalised maxima over ζ0 ∈ [0.02, 0.999], b ∈ [5.5, 1e5], α ∈ [0, 1e5], all ζ:  C ≤ −0.081, D ≤ −0.196.
  Worst C at ζ -> 0, α large = the inner/outer overlap, consistent with the inner analysis.
  (A first run gave D = +0.316: sign error in g∞' — caught because the finite-t check passed there; fixed.)
INNER zone z = O(1) (item 71): D ratio ≤ −1/3 at leading order; C ⇔ 1.7(z−1)·(gz) + p ≤ 2z − 1, margin ≥ 0.3(z−1) + (1−p).
Finite-t cross-check: check(...) at t = 1e2..1e4, zf = 0.05, α ∈ {0.01, 0.1, 1}: all pass, D ≤ −0.156.
REMAINING for this corner: (i) sub-limits b -> inf, α -> inf of the outer problem (closed forms); (ii) the
  overlap/matching estimate (both zone limits agree on 1 ≪ z ≪ t); (iii) explicit O(1/t) error bounds so the corner
  is a theorem for t ≥ T0, with the bounded part t ≤ T0 going to the R3 certification.

## 74. A -> inf corner, repaired multiplier: limits are explicit

Sequences A = 10^2..10^10 (B at floor, and A,B -> inf together), t ∈ {1.05, 2, 10, 100}, zf ∈ {0.05, 0.5, 0.95}: all pass,
all converge.  D-margin -> −0.2590 for EVERY (t, zf);  C-margin -> −0.379 (= (p−1)/(p+1), z -> 1 end) or a t-dependent
value ≥ −0.119 at t = 100, zf = 0.95 (joins the t -> inf corner, item 73, where C -> −0.081).
[P, leading order] Bulk D limit = −δ/(2+δ): as A -> inf off the seam, sLh ~ A and s|L'| ~ A dominate the μ-terms (O(1)·h);
  with Lh = ψ|L'| and ψ -> 1+δ, D/nD -> (1−ψ)/(1+ψ) = −δ/(2+δ)  (= −0.2593 at δ = 0.7; = −1/3 at δ = 1, matching item 53).
  Seam (ψ -> 0) is the A-independent rescaled layer problem, items 13–15 (leading order proved); numerically never the max.
REMAINING: explicit O(A^{-1/2}) bounds for (i) the μ-term corrections in the bulk and (ii) the layer/bulk overlap.

## 75. THEOREM B REDESIGN after the referee (REFEREE.md, items I.1–I.3)

Literature check: GHZ 1908.09095 Thm 1 / Wald 1973 reduce "ψ0 = ψ4 = 0 ⇒ PD + gauge" to Wald's linearised type-D
classification, cited without derivation; for Λ = 0 AHW needed AABKW completeness to make it rigorous.  No Λ > 0
analogue exists.  So B1 stays on our route, repaired as follows.

(R1) FORMAL INTEGRABILITY IS NOT NEEDED (answers I.2).  Khavkine 1805.03751: the KdS Killing operator is of regular
     finite type, so it has a complete compatibility operator K1 of some finite order N: h is locally pure gauge iff
     K1(h) = 0.  K1(h)(q) depends only on the N-jet of h at q modulo gauge N-jets.  So if d_N(q) = 0, i.e. every N-jet of
     a true solution at q is a gauge N-jet (an upper bound, certified by truncation at K > N exactly as now), then
     K1(h)(q) = K1(gauge)(q) = 0.  If d_N = 0 on an open dense set of q, then K1(h) ≡ 0 and h is locally pure gauge,
     with no prolongation or propagation argument.  Item 58(b) is replaced by this; the closure claim there is withdrawn.
     NEEDED: an explicit bound on N.  (For Kerr the AABKW invariants have order ≤ 3; Khavkine's construction gives N for
     KdS from the Killing-operator prolongation; expected N ≤ 3, TO DETERMINE.)  Then the certificates must use J = N
     (we have J ∈ {0, 1}; pointwise exact data already show d_J = 0 up to J = 4–6, items 61, 64).

(R2) GEOMETRY UNIFORMITY (answers I.1).  d_N = n_L − rank C + rank C_H − rank G.  The UPPER bound on rank C_H comes from
     explicit gauge rows (step56b), which are solutions for EVERY (Λ, a), so it holds identically.  The lower bounds on
     rank C and rank G hold off the zero set of exact minors, which are polynomials in (Λ, a, r, x, σ, μ).  Requirement:
     for each physical (Λ, a, σ, m) the minors do not vanish identically in (r, x).  That is a statement about a
     subvariety of (Λ, a, σ, μ)-space, to be attacked with the same exact/mod-p machinery with (Λ, a) as extra variables
     (degree bounds by Toeplitz deficiency, as in items 69–71).

(R3) THE LINES σ = 0 AND m = 0 (answers I.3).  On each line the generic ranks differ (σ⁶μ⁶ factor), so each needs its own
     minors and its own 1-parameter exceptional-polynomial certificate (in μ on σ = 0, in σ on m = 0), plus the PD sector
     at σ = m = 0 (d = 4, item 64).

(R4) ASSEMBLY GAPS (B.4, referee II): generalised modes with t^j terms, a top term c t^J ∂_M g at σ = 0, and the form of the
     gauge 1-form.  To write; expected routine, following HV18 §7 and AHW §6.4.

Order of work: (R1) determine N; redo the exact pointwise certificates at J = N; then (R3); then (R2); then (R4).

## 76. Killing stabilisation depth on KdS: N' = 1 (EXACT) — input to the order of Khavkine's K1

step63_killing_stabilisation.py (exact sympy rationals, coordinates (t, r, x = cos th, phi), Lam = 1/20, a = 3/5, M = 1):
unknowns = Killing 1-jet (xi^a, w_ab = nabla_a xi_b) at p, 10-dimensional; conditions L_xi(nabla^j Riem)(p) = 0.
   p = (3, 3/5):     j = 0: rank 6, kernel 4;   j ≤ 1: rank 8, kernel 2.
   p = (7/2, -1/5):  same (step63b_point2.py).
Kernel 2 = the true Killing algebra (d_t, d_phi); both give EXACT zero residual (sanity check).  Since the kernel already
equals the dimension of the Killing algebra, further conditions cannot cut it, so the integrability conditions
stabilise at depth N' = 1 (the kernel-4 stage = Killing fields + the 2-dim type-D isotropy boost/rotation of the Weyl tensor).
(Bug caught on the way: sympy.nsimplify turned exact rationals into algebraic numbers; replaced by exact Rational.)
CONSEQUENCE [inference from Khavkine Thm 9 / Prop 13 structure, TO CONFIRM against the paper]: the restricted flat
connection lives on the rank-2 subbundle cut out by L_xi R = L_xi nabla R = 0; the deformed constraints for a
perturbation h involve nabla(linearised Riem), i.e. 3 derivatives of h, and the curvature of the deformed connection 2–3.
So ord K1 = N' + 2 = 3 (consistent with Kerr, where the AABKW complete invariants have order 3).
=> Theorem B certificates must be run at jet order J = 3 (not {0, 1}).  Pointwise exact data (items 61, 64) already give
   d_J = 0 for J ≤ 4–6 at the sampled (sigma, m), so the J = 3 certificate is expected to pass.

## 77. Pointwise exact certificate at jet order J = 3 = ord K1 (step64_sample_J3.py)

K = 6 truncation (honest through J = 3; truncation only over-counts, so every reported d_J is a valid UPPER bound), all four
Teukolsky conditions, exact Q(i): 21 random rational samples, Lam M^2 ∈ [0.042, 0.163] (max subextremal ≈ 0.177),
a/M ∈ [0.33, 1.07], points across (r_e, r_c) × (−1, 1), σ real / imaginary / complex / ZERO, m ∈ [−4, 4]
(including σ = 0 with |m| = 3, 4 — the sector referee I.3 flagged as unchecked):  d_0 = … = d_3 = 0 (and d_4 = 0) in ALL.
With item 75 (R1) + item 76 (ord K1 = 3): at each sampled (geometry, σ, m) and point, every solution's 3-jet is a gauge
3-jet, so K1(h) vanishes there.  STILL OWED for a proof: uniformity (the good set must be open-dense in the point for
EVERY geometry and EVERY (σ, m), including the lines σ = 0 and m = 0) — items 75 (R2), (R3); and confirming ord K1 = 3.

## 78. The lines m = 0 and σ = 0 at jet order J = 3 (answers referee I.3, modulo exact lift) — step65_line_cert.py

K = 6, all four Teukolsky conditions, Lam = 1/20, a = 3/5, points p1 = (3, 3/5), p2 = (7/2, −1/5).  Every rank bound is
backed by a fixed minor in the exceptional polynomial E, INCLUDING the gauge-row minor behind the R_H upper bound (rank
of the gauge rows can drop at special parameter values; previously treated as automatic).  Only J = 3 minors enter E.
 - Line m = 0 (variable σ; interpolation nodes σ = 1, 2, 3 to avoid the σ = 0 shape change): generic d_3 ≤ 0 at both
   points; deg E = 70, ord_0 = 32; gcd(E1, E2) mod P = c σ³² EXACTLY a monomial.
 - Line σ = 0 (variable m; nodes m = 1, 2, 3): generic d_3 ≤ 0 at both points; deg E = 69, ord_0 = 35;
   gcd(E1, E2) mod P = c m³⁵.
So on each line the only parameter value that is exceptional at both points is the origin (σ, m) = (0, 0), the PD sector
(exact d = 4, item 64).  (Generic J = 0 bound is 2: the known looseness of the gauge-row R_H bound at low J; irrelevant,
since only J = ord K1 = 3 is needed.)
OWED for exactness (as in item 67): exact degree and exact lowest order of E1 on each line (Toeplitz / CRT, as in items
67–71), and then the same for every geometry (item 75 R2).

## 79. Generic (σ, m) certificate at jet order J = 3 (mod P, fixed geometry) — step66_biv_J3.py, biv3.log

K = 6, J = 3 minors only (constraint, low-gauge, gauge-row R_H minor), Lam = 1/20, a = 3/5, three points
(3, 3/5), (7/2, −1/5), (4, 1/3).  Matrices from 6 exact builds at nodes (σ, m) ∈ {(1,1),(2,1),(3,1),(1,2),(1,3),(2,2)}.
  generic d_3 bound = 0 at all three points;  E_i: deg_σ 65, deg_μ 46 (80×80 grid, top rows/cols zero), E_i = σ¹⁹μ⁷H_i;
  R_ij = Res_σ(H_i, H_j): degree 3486, ord_{μ=0} = 64 for all three pairs;
  gcd(R12, R13, R23) mod P = c μ⁶⁴ (monomial).   [613 s total]
With the item-70 lemma: no integer m ≠ 0 has σ ≠ 0 exceptional at all three points, PROVIDED the exact lifts:
  (i) exact σ¹⁹ μ⁷ divisibility (Toeplitz deficiency at σ = 0 / μ = 0, as items 71, 60);
  (ii) exact formal degrees (deg_σ E_i = 65; reversed Toeplitz at fixed μ0 + interpolation, as items 69, step58i);
  (iii) exact ord_{μ=0} R_12 = 64 (origin multiplicity of H_i, as item 70/step59b).
Lines m = 0 and σ = 0: item 78.  Together: at THIS geometry, B1 holds for every (σ, m) except the origin (PD sector),
modulo (i)–(iii) and the line lifts, and modulo ord K1 = 3 (item 76).
REMAINING for Theorem B: those exact lifts; geometry uniformity (item 75 R2); assembly (R4); confirming ord K1 = 3.

## 80. KEY REDUCTION: (σ, m)-uniformity follows from ONE 4D jet statement; geometry uniformity becomes (Λ, a, r, x) only

Float 4D jets (step44, all perturbations, δRic − Λh = 0, δψ0 = δψ4 = 0, real h), point r = 3, θ = 1:
   K = 5: d = {0,0,0, 6*};   K = 6: d = {0,0,0, 4, 6*}   (* = top truncation level, artifact), for Λ = 0 and Λ = 0.05.
   So honest d_3 = 4 = number of Λ-PD tangents (δM, δa, NUT, acceleration).
LEMMA (Killing-field descent) [P, to write in full].  Suppose at q:
   (4D_4)  every 4-jet at q of a local solution = (gauge 4-jet) + Σ c_i (4-jet of PD tangent i).
 Let h be a mode, L_T h = −iσ h, L_Φ h = i m h.  T, Φ are Killing: L_T commutes with the linearised operators, maps
 gauge to gauge (L_T δ*X = δ*[T, X]) and kills the PD tangents (stationary, axisymmetric).  The 3-jet of L_T h at q is
 determined by the 4-jet of h at q.  Hence −iσ (3-jet of h) = 3-jet of L_T h ∈ gauge 3-jets; likewise i m (3-jet of h).
 So for (σ, m) ≠ (0, 0) the 3-jet of every mode solution at q is a gauge 3-jet: d_3^mode = 0 for ALL (σ, m) ≠ 0 at once.
 With ord K1 = 3 (item 76) and item 75 (R1): if (4D_4) holds on an open dense set of q, every non-origin mode with
 vanishing Teukolsky scalars is locally pure gauge.
CONSEQUENCE: items 65–70, 78, 79 (σ/m exceptional polynomials, resultants, lines) are SUPERSEDED for (σ, m) ≠ (0, 0).
 Theorem B's local step reduces to: (4D_4) on an open dense set of points, for every subextremal geometry — a problem in
 (Λ, a, r, x) only — plus the origin (PD sector, item 64) and the assembly (R4).
NEXT: exact / mod-p 4D jet builder in coordinates (t, r, x = cos θ, φ) (metric rational); certify (4D_4) as
 d_4 ≤ 4 (upper bound by truncation, K ≥ 7) + the four PD tangents' 4-jets independent mod gauge (explicit Λ-PD family);
 then geometry uniformity via exceptional polynomials in (Λ, a, r, x).

## 81. 4D jet system assembled from the exact mode matrices (step67_4d_from_modes.py)

Mode coefficient matrices (exact builds at 6 nodes, reduced mod P) -> full 4D jet system via σ <-> i∂_t, m <-> −i∂_φ on
t^j φ^k (block-triangular), total-degree truncation (consistent since all operators have order ≤ 2).
K = 6, Lam = 1/20, a = 3/5, point (3, 3/5): unknowns 2100, equations 980, gauge generators 1320.
  Sanity: gauge images are solutions EXACTLY mod P (G·C = 0).
  Direct mod-P counts (nullspace; heuristic only — P could in principle be unlucky): d = {0, 0, 0, 4, 6*}, matching the
  independent float 4D code (step44) exactly.  So honest d_3 = 4 = PD.
  Rigorous upper bound n_L − rank C + (n_H − rank X_high) − rank G_L is LOOSE: J = 3 gives 6, not 4 (the gauge-row
  bound on rank C_H only sees generators of degree ≥ J+2 and misses other solutions supported in high degree).
FIX (next): lower-bound dim(S ∩ ker π_J) by EXACT elements — compute the mod-P nullspace of C_H, rationally
  reconstruct (CRT over a few primes), and verify v·C = 0 EXACTLY over Q(i) (a cheap matrix–vector product with the exact
  coefficient matrices).  That makes the upper bound tight and rigorous.  Then K = 7 for J = 4, then the descent check
  (item 80, as an exact containment), then geometry uniformity.

## 82. RIGOROUS d_3 ≤ 4 on full 4D jets (one point) — step68_4d_exact.py, ex4d.log

d_J ≤ n_L − rank C + rank C_H − rank G_L, with rank C a lower bound (any) and rank C_H EXACT over Q(i) (upper bound is
what matters), rank G_L a mod-P lower bound.  K = 6, Lam = 1/20, a = 3/5, point (3, 3/5), exact 4D assembly 2100 × 980:
   n_L = 350, rank_Q C = 776, rank_Q C_H = 708, rank_P G_L = 278   =>   d_3 ≤ 4   RIGOROUSLY.   (388 s)
Combined with the PD tangents (4 independent classes, float item 40 / step44) this is d_3 = 4.
K = 7 (needed for honest J = 4, the descent lemma's input): exact 4D assembly 3300 × 1764 done (475 s), but the exact
realified rank (6600 × 3528 dense fmpq) died — memory.  NEEDED: a leaner exact certificate for rank C_H ≤ r
(only an UPPER bound is needed): exhibit n_H − r exact left-null vectors of C_H — mod-P nullspace, CRT/rational
reconstruction, verification v·C_H = 0 by sparse exact products (the assembly is a sparse dict, cheap).  rank C may stay
mod P (lower bound).

## 83. RIGOROUS d_4 ≤ 4 on full 4D jets (one point), K = 7 — step69_K7_lean.py, k7lean.log

Exact rank of C_H via a row-scaled INTEGER realified fmpz_mat (5200 × 3528; peak 5.7 GB, 1135 s):
   J = 4: n_L = 700, rank_P C = 1318, rank_Q C_H = 1124, rank_P G_L = 502   =>   d_4 ≤ 4  RIGOROUSLY.
(Lam = 1/20, a = 3/5, point (3, 3/5).)  With item 82: d_3 ≤ 4, d_4 ≤ 4.
DESCENT INPUT STILL OWED: V_4 (solution 4-jets mod gauge, dim ≤ 4) must be spanned by classes of STATIONARY, AXISYMMETRIC
solutions, so that L_T = ∂_t and L_Φ = ∂_φ kill it (item 80).  Since dim V_4 ≤ 4, it suffices to exhibit 4 such solutions
whose 3-jets are independent modulo gauge 3-jets.  Candidates: the four Λ-PD tangents.  δM, δa are explicit derivatives of
the KdS metric (rational).  NUT and acceleration need the Λ-Plebański–Demiański metric (Griffiths–Podolský form),
linearised at KdS.  Independence check is SMALL and exact: gauge 3-jets come from X-jets of degree ≤ 4 (280 × 350),
so rank_Q[PD 3-jets; G_3] − rank_Q G_3 = 4 is a tiny exact computation.  (Alternative avoiding PD: exact containment
∂_t π_4(S) ⊂ π_3(G) — needs an exact nullspace, heavier.)

## 84. DESCENT HYPOTHESIS (4D_4) PROVED at one point — step70_pd_exact.py, pd70.log

Exact Carter/PD ↔ BL map: tau = (t − aφ)/Ξ, sig = −φ/(aΞ), p = a x, with m = M, n = 0, eps = 1 − Λa²/3, k = a², alpha = 0.
 (0) the pulled-back PD background equals the BL KdS metric EXACTLY (sympy cancel).
 The 5 tangents ∂g/∂(m, n, eps, k, alpha), pulled back (exact rational), 3-jets at (3, 3/5), Λ = 1/20, a = 3/5:
 against the FULL 4D gauge 3-jets (280 generators × 350 jet unknowns, exact over Q(i)):
     rank G_3 = 278 (matches the mod-P rank in item 82),  rank [G_3; PD] = 282  =>  PD rank mod gauge = 4  EXACTLY.
 (A first attempt projecting gauge onto t,φ-independent jets gave rank 100 = everything: projection discards the t·X
  terms' t-dependent parts and is far too generous; replaced by the full 4D gauge system.)
CONSEQUENCE at this point:  V_3 has 4 independent stationary-axisymmetric classes and dim V_3 ≤ 4 (item 82), so
 V_3 = span(PD); PD independent at J = 3 ⇒ independent at J = 4, and dim V_4 ≤ 4 (item 83) ⇒ V_4 = span(PD).
 Hence (4D_4) holds, and by the descent lemma (item 80) every mode solution with vanishing Teukolsky scalars and
 (σ, m) ≠ (0, 0) has 3-jet = gauge 3-jet at this point — for ALL (σ, m) at once.
REMAINING for Theorem B's local step: (U) the same on an open dense set of points for EVERY subextremal geometry
 (exceptional polynomial in (Λ, a, r, x)); (O) ord K1 = 3 confirmed against Khavkine; (P) origin (PD sector);
 (R4) assembly of generalised modes.

## 85. UNIFORMITY MECHANISM: at K = 7 the certificate uses only OPEN conditions (step67 K=7, step71)

step71 (K = 6, J = 3): exact nullity of C_H = 1042 vs gauge_high rank 1040: the 2 extra high solutions are supported in
degree ≥ J+1 = 4 (support analysis s = 4: excess 2; s = 5, 6: excess 0) — truncation artifacts of K = 6 at J = 4, not
top-degree symbol kernel (the 480-dim top kernel lies inside gauge_high).
step67 at K = 7, structural bound d_J ≤ n_L − rank C + (n_H − rank X_high) − rank G_L, ALL three ranks mod-P lower bounds:
   J = 3: 4,  J = 4: 4,  J = 5: 4   (J ≤ 2 loose, irrelevant).
So with K = 7 every ingredient of the descent certificate (item 80) is an OPEN (Zariski) condition in (Λ, a, r, x):
   d_4 ≤ 4 (three rank LOWER bounds),  PD independence mod gauge (rank lower bound, item 84),
 while "gauge images are solutions" and "PD tangents are solutions" are identities in the parameters.
 ⇒ the certificate holds on a Zariski-open dense U ⊂ (Λ, a, r, x) containing (1/20, 3/5, 3, 3/5).
ANALYTICITY (mode equations elliptic in (r, x) in harmonic gauge between the horizons ⇒ K1(h) real-analytic on the
 connected exterior): for a geometry (Λ, a), K1(h) vanishes on the open set U_(Λ,a) ∩ exterior ⇒ everywhere.
 Note U_(Λ,a) ≠ ∅ as soon as the minor polynomials are not identically zero in (r, x) for that (Λ, a); unphysical
 (r, x) also witness this (Zariski density), so physical points are not needed for the witness.
REMAINING CRUX (U'): exclude the algebraic set B = {(Λ, a): the certificate minors vanish identically in (r, x)}
 (dim ≤ 1) from the physical region.  Direct bivariate interpolation in (Λ, a) is too expensive (one K = 7 evaluation ≈ 7 min).

## 86. Uniformity (U'): the analytic route needs a Λ-AABKW theorem — status of Theorem B

AHW 2207.12952 App. A (Kerr): assuming AABKW completeness of {linearised Ricci, ϑΨ±2, I_ξ, I_ζ}, the gradients of I_ξ, I_ζ
are fixed by a gauge-invariant Killing field Im Ā = Aξ + Bζ, giving ODEs (A.25) whose solutions (A.26) are exactly the PD
tangents.  For Λ > 0 every step needs Λ source terms (linearised Bianchi identities, Killing-spinor identities), and the
completeness input is a Λ-version of AABKW (complete local gauge invariants on KdS, for ALL parameters) — a paper-scale
symbolic computation in spinor/GHP form, not available.  Our K = 7 jet certificate is a POINTWISE substitute for it.
STATUS of Theorem B's local step (modulo ord K1 = 3 [item 76, to confirm], origin/PD sector, assembly R4):
  proved at every (Λ, a) outside an algebraic set B of dimension ≤ 1 in the (Λ, a)-plane (Zariski-open certificate +
  analyticity, item 85), and at a = 0 / small |a| by Hintz–Vasy 2018 §7 (SdS UEMS + perturbation).
  B is not known to miss the physical region.  Routes: (i) Λ-AABKW (analytic, all parameters); (ii) a much smaller
  pointwise certificate whose minors can be computed symbolically in (Λ, a).

## 87. Fast mod-P certificate and the bad set B (step52p_modp_modes.py, step72_fast_cert.py, step73_probe.py)

mod-P port of the mode builder (nmod_mat): K = 7 descent certificate in 40 s (was ~450 s); reproduces d_4 ≤ 4 at the
reference point.  Probe of the structural (open-condition) certificate d_4 ≤ 4:
   a = 0 (SdS, Λ = 0.05):   d_4 ≤ 8  — FAILS, as predicted (extra spherical Killing fields): a = 0 ∈ B.
   a = 0.001:               d_4 ≤ 4.
   12 random subextremal geometries (Λ ∈ [0.008, 0.153], a ∈ [0.05, 1.02], points across the exterior): all d_4 ≤ 4.
EVIDENCE: B ∩ physical region = {a = 0}, which is covered by Hintz–Vasy 2018 §7 (SdS UEMS) — and small |a| perturbatively.
NOT A PROOF: B could contain other curves.  A proof needs either a Λ-AABKW theorem or a certificate small enough for
symbolic minors in (Λ, a) (the current minors have size ~1300; degree in a likely ~10^3–10^4, too many evaluations).

## 88. A.2, A.9 (m != 0) and A.10 PROVED BY HAND — the window is an explicit interval
Real omega, w = omega/m (m != 0), a > 0 wlog.  rho_j = r_j^2 + a^2, Omega_j = a/rho_j, D_j = |Delta'(r_j)|.
Delta_r = -(Lam/3) prod (r - r_k) gives sign Delta'(r_j) = (-, +, -, +) for j = 0,1,2,3.  With P_j = rho_j/D_j:
  Im e_1 = 2 Xi m f1(w),  f1 = P1(w - Omega1) + P3(w - Omega3)   (increasing, root w1 = (a/D1 + a/D3)/(P1 + P3)),
  Im e_t = 2 Xi m f2(w),  f2 = -P2(w - Omega2) + P3(w - Omega3).
Residue identities ((r^2+a^2)/Delta and 1/Delta have zero residue at infinity):
  P1 + P3 = P0 + P2,   1/D1 + 1/D3 = 1/D0 + 1/D2.
Hence slope(f2) = P3 - P2 = P0 - P1 and root w2 = a(1/D1 - 1/D0)/(P1 - P0) = item 46's "lo".
(i) P1 > P0: D0/D1 = f(r0)/f(r1) > 1 with f(x) = (x - r3)(r2 - x) decreasing on [r0, r1] (item 46), and rho1 > rho0.
    So f2 is strictly decreasing, and w2 > 0.
(ii) w2 < w1: cross-multiplying, with x_j = 1/D_j, the difference is
    x0x1(rho1-rho0) + x1x3(rho1-rho3) + x0x3(rho3-rho0)   (a^2 cancels: rho_i - rho_j = r_i^2 - r_j^2),
    and D0 D1 D3 times it is exactly (sympy factor, no Vieta needed)
       2 (r0 - r1)(r0 - r3)(r1 - r3)(r0 r3 - r1 r2)  > 0,
    signs (-)(+)(+)(-).  [With r3 = -(r0+r1+r2): -2(r0-r1)(r0+r1)(r0+r2)(r0+2r1+r2)(2r0+r1+r2).]
Consequences (all [P]):
  * Sign agreement {Im e_1 . Im e_t > 0} is EXACTLY the open interval (w2, w1) (the factor m^2 > 0 drops out).
    Outside it the Wronskian identity kills the mode.  This is A.2, with no numerics.
  * w2 > 0, so omega = 0 is outside the window for every m != 0: A.9's numerical claim is proved.
  * At w = w2 only Im e_t vanishes (f1(w2) < 0); at w = w1 only Im e_1 vanishes (f2(w1) < 0).  A.10's "exactly one"
    is proved; item 44's two endpoint arguments then close the endpoints.
  * m = 0, omega real != 0: Im e_1 ∝ omega (P1 + P3), Im e_t ∝ omega (P0 - P1): opposite signs, Wronskian kills it.
  * a = 0: w1 = w2 = 0, empty window.

## 89. c1 IDENTITY (n1 half) PROVED SYMBOLICALLY in 7 s — step89_c1_kappa.py
Item 31's script never finished the n1 half (sympy blow-up).  Rewritten in the kappa variables:
  * scale to Lam = 3 (r -> c r, a -> c a, omega -> omega/c, Lam -> Lam/c^2 leaves kappa_j, t, Lam C0 r^2, lam invariant);
  * K(r3) = -i D1 kappa3 / Xi gives N0/D1^2 = -kappa3^2 - s kappa3 and N1 in closed form; the only other omega is
    W = i Xi omega;
  * the residue identities sum rho_j/Delta'_j = sum 1/Delta'_j = sum r_j^2 rho_j/Delta'_j = sum r_j^2/Delta'_j = 0 and
    sum r_j rho_j/Delta'_j = -3/Lam give the exact linear relations
       sum kappa_j = 0,   sum r_j^2 kappa_j = 0,   sum r_j kappa_j = -W,
    so (kappa_0, kappa_1) are free, and a, m, omega disappear.
Result, as an identity of rational functions in (r0, r1, r2, kappa0, kappa1, lam, s), with r3 = -(r0+r1+r2) and
alpha = 1 + e2 (Vieta at Lam = 3):
   n1 - (lam C0 + 1 + t + 3(1-alpha) C0) = -2 (s^2 - 4) (r0^2 + r0r1 + r0r2 - r1r2) / ((r0 - r1)(r0 + r1 + 2 r2)),
   n2 + 3 = 4 - s^2.
Both vanish at s = +-2.  The c1 identity is now fully proved.

## 90. NEW ROUTE for the bad set B (replaces the jet certificate, items 80–87): Λ-Wald via linearised type D integration
Step 1 [P]: LINEARISED GOLDBERG–SACHS WITH Λ.
  Background: KdS, principal tetrad (l, n, m, m̄), GHP, κ = σ = κ' = σ' = 0, Ψ0 = Ψ1 = Ψ3 = Ψ4 = 0, Φ_ab = 0,
  Ψ2 ≠ 0 everywhere in the exterior.  Perturbation h with δΨ0 = δΨ4 = 0 (gauge invariant: Stewart–Walker, since the
  background values vanish; holds for modes by Theorem A, item B.1).
  (a) Tetrad gauge.  To first order an infinitesimal null rotation about n with parameter d sends δΨ1 -> δΨ1 + 3dΨ2,
      δΨ0 -> δΨ0 (the change is 4dΨ1 = 0 on the background); about l with parameter c, δΨ3 -> δΨ3 + 3c̄Ψ2, δΨ4 unchanged.
      δΨ2 is unchanged by both.  Since Ψ2 ≠ 0, choose d = −δΨ1/(3Ψ2), c̄ = −δΨ3/(3Ψ2): now δΨ0 = δΨ1 = δΨ3 = δΨ4 = 0.
      (Pointwise algebraic, smooth in the data; c, d inherit the mode dependence e^{−iσt+imφ}.)
  (b) The vacuum-Λ Bianchi identities in GHP are those of vacuum (∇R = 0, Φ = 0: Λ drops out).  Linearising
         Þψ1 − ð'ψ0 = −τ'ψ0 + 4ρψ1 − 3κψ2      gives  0 = −3Ψ2 δκ,
         Þ'ψ0 − ðψ1 = ρ'ψ0 − 4τψ1 + 3σψ2       gives  0 =  3Ψ2 δσ,
      and their primes give δκ' = δσ' = 0.  Hence in this tetrad gauge the perturbed spacetime is type D to first
      order with both principal null directions geodesic and shear-free: h solves the LINEARISED type D Λ-vacuum
      GHP system for (ρ, ρ', τ, τ', Ψ2) together with the linearised tetrad/commutator relations.
Step 2 [OPEN, the real work]: FINITE TYPE.  In the type D Λ-vacuum GHP system all four directional derivatives of
  (ρ, ρ', τ, τ', Ψ2) are determined algebraically by those quantities (Þρ = ρ², ðρ = (ρ−ρ̄)τ, Þτ = (τ − τ̄')ρ,
  Bianchi Þψ2 = 3ρψ2 ... ; the remaining derivatives such as ð'ρ, Þ'ρ come from the Ricci identity with the Λ-term
  and the commutators applied to Ψ2 — Debever–Kamran–McLenaghan, J. Math. Phys. 25 (1984), exhaustive integration
  with Λ).  Claim to prove: on the GENERIC branch (τ ≢ 0, ρ ≢ ρ̄), the linearisation of this closed system has
  solution space, modulo boost-spin and diffeomorphism gauge, of dimension 4 = tangents of PD (5 constants − scaling,
  with Λ fixed), and these are δM, δa, NUT, acceleration.
Step 3 [P once step 2 holds]: for KdS, τ and Im ρ are ∝ a·cos θ (ρ = −1/(r − i a cos θ) in Kinnersley normalisation,
  τ ∝ i a sin θ/ρ²-type), so every exterior point with a ≠ 0, cos θ ≠ 0, sin θ ≠ 0 is on the generic branch.
  ⇒ B ∩ {a ≠ 0} = ∅ for EVERY geometry, by hand.  Also retires the jet certificate, the descent lemma and ord K1 = 3.
Risk: step 2 requires the integrability conditions of the closed system to hold identically (formal integrability) at
  the background, so the linearised count equals the nonlinear one; check = a few GHP identities with Λ.

## 91. Step 2 of item 90, first order: the type D Λ-vacuum GHP system is NOT closed at first order
Unknowns U = (ρ, ρ', τ, τ', Ψ2) (+ conjugates), κ = σ = κ' = σ' = 0, Φ = 0, GHP Λ = R/24.
Given directly (Ricci identities + Bianchi):
  Þρ = ρ²,  ðρ = (ρ − ρ̄)τ,  Þτ = (τ − τ̄')ρ,  ðτ = τ²   (and primes: Þ'ρ' = ρ'², ð'ρ' = (ρ' − ρ̄')τ', Þ'τ' = (τ' − τ̄)ρ', ð'τ' = τ'²),
  (f)  Þ'ρ − ð'τ = ρρ̄' − ττ̄ − Ψ2 − 2Λ,   (f') Þρ' − ðτ' = ρ'ρ̄ − τ'τ̄' − Ψ2 − 2Λ,
  ÞΨ2 = 3ρΨ2, ðΨ2 = 3τΨ2, Þ'Ψ2 = 3ρ'Ψ2, ð'Ψ2 = 3τ'Ψ2.
Commutators applied to Ψ2 (weight (0,0)) [derived by hand; to be machine-checked on the KdS tetrad]:
  [Þ, Þ']Ψ2:  Þρ' − Þ'ρ = ττ̄ − τ'τ̄'
  [ð, ð']Ψ2:  ðτ' − ð'τ = ρρ̄' − ρ̄ρ'
  [Þ, ð]Ψ2:   identity (consistent with the given Þτ, ðρ)
  [Þ, ð']Ψ2:  Þτ' = ð'ρ;    prime:  Þ'τ = ðρ'.
  (f) − (f') follows from the two commutator relations, so only 3 independent equations constrain the 4 quantities
  (Þ'ρ, ð'τ, Þρ', ðτ'), and ð'ρ = Þτ', Þ'τ = ðρ' are unconstrained at this order.
=> At first order 3 complex first-derivative quantities are free.  Closure (finite type) must come from the second-order
   prolongation: the commutators applied to ρ, τ, ρ', τ' themselves.  That is the next computation; it is where
   Debever–Kamran–McLenaghan / Kinnersley obtain the relations Ψ2 ∝ ρ³-type (Kerr: Ψ2 = Mρ³ in the Kinnersley gauge).

## 92. Item 91 MACHINE-CHECKED on KdS (step92_ghp_check.py, ~20 s)
Exact KdS metric (Λ = 1/20, a = 3/5, M = 1), Kinnersley-type tetrad (l·n = 1, m·m̄ = −1, signature +−−−), spin
coefficients from Christoffels, GHP derivatives by 40-digit numerical differentiation at (r, θ) = (3, 0.7).
κ, σ, κ', σ' vanish (~1e-125).  All nine first-order relations of item 91 hold to ≤ 5e-42:
  Þρ = ρ², ðρ = (ρ−ρ̄)τ, Þτ = (τ−τ̄')ρ, ðτ = τ², Þ'ρ' = ρ'², Þρ' − Þ'ρ = ττ̄ − τ'τ̄', ðτ' − ð'τ = ρρ̄' − ρ̄ρ',
  Þτ' = ð'ρ, Þ'τ = ðρ'.
Convention fixed by the check (the first attempt failed exactly the primed relations):
  Þ = D − pε − qε̄,  ð = δ − pβ − qᾱ,  Þ' = Δ − pγ − qγ̄,  ð' = δ̄ − pα − qβ̄,  ρ' = −μ, τ' = −π.
The script is a reusable numeric oracle for every GHP identity needed in the second-order prolongation.

## 93. Step 2 continued: Ψ2, tetrad derivative table, the Killing-spinor vector (all machine-checked on KdS)
(a) (f) fixes Ψ2 = −M/(r − i a cos θ)³ (no Ξ factor) and Λ_GHP = Λ/6 (check: F = 1/c³ − Λ/3 to 15 digits).
    So φ := Ψ2^{−1/3} ∝ r − i a cos θ.
(b) GHP derivatives of the tetrad on type D Λ-vacuum (step93_tetrad_table.py; every entry matched numerically):
      Þl = 0,            Þn = −τ'm − τ̄'m̄,   Þm = −τ̄'l,          Þm̄ = −τ'l,
      Þ'l = −τ̄m − τm̄,   Þ'n = 0,            Þ'm = −τn,           Þ'm̄ = −τ̄n,
      ðl = −ρ̄m,          ðn = −ρ'm,          ðm = 0,              ðm̄ = −ρ'l − ρ̄n,
      ð'l = −ρm̄,         ð'n = −ρ̄'m̄,         ð'm = −ρ̄'l − ρn,     ð'm̄ = 0.
(c) Killing-spinor vector: ξ with ξ·l = φρ, ξ·n = −φρ', ξ·m = φτ, ξ·m̄ = −φτ' equals −Ξ ∂_t on KdS
    (step93_killing.py, the other 7 sign patterns fail).  All 10 Killing equations, written with (b), reduce to
    identities or to the already-known relations Þτ' = ð'ρ, Þ'τ = ðρ', C1, C2.  So ξ is Killing on EVERY type D
    Λ-vacuum, by algebra alone — but the Killing equations give NO new first-order relation.
(d) Where closure comes from: ξ is complex in general; on KdS it is a complex constant k times a REAL Killing field.
    Reality, ξ̄ = (k̄/k) ξ, is the algebraic closure:
       k̄ φρ = k φ̄ρ̄,   k̄ φρ' = k φ̄ρ̄',   k̄ φτ = −k φ̄τ̄'   (fixes ρ̄, ρ̄', τ̄' in terms of ρ, ρ', τ, and the phase).
    LINEARISED: the reality of ξ is NOT automatic for a perturbation.  Im δ(ξ/k) solves the linearised Killing
    equation with zero background field, so it is a Killing field of the background: Im δ(ξ/k) = A ∂_t + B ∂_φ with
    constants A, B.  This is the Λ-analogue of AHW 2207.12952 App. A, "Im Ā = Aξ + Bζ", and it is reached here
    WITHOUT the AABKW completeness theorem: steps 1 (item 90) + (b)–(c) give it directly.
NEXT: (i) integrate the linearised reality relations with the constants (A, B, δk) — AHW's ODEs (A.25) with Λ terms;
   (ii) the second Killing field ζ = ∂_φ from the Killing tensor (needed for B); (iii) count: expected 4 = (δM, δa,
   NUT, acceleration) after gauge, giving B ∩ {a ≠ 0} = ∅.

## 94. Step 2, structure of the count (route fixed; three lemmas to prove)
Normalise φ so that ξ = ∇·K is real on the background (ξ = −Ξ∂_t; on KdS with the Kinnersley tetrad φρ = −1).
Given a perturbation h with δΨ0 = δΨ4 = 0, in the gauge of item 90:
 L1 [P, item 93(d)]  Im δ(ξ/k) = A ∂_t + B ∂_φ, constants A, B (any δk).
 L2 [to prove]  Subtracting a suitable combination of PD tangents kills (A, B, δk-phase).  Expected: the acceleration
    tangent has Im δξ ≠ 0 (the Killing-spinor vector of the accelerating PD metric is not real up to a constant —
    classical: reality of ξ characterises the Kerr–NUT–(A)dS / Carter [Ã] subclass, Hughston–Sommers 1973,
    Debever–Kamran–McLenaghan 1984), and one of (A, B) is fixed by the normalisation freedom of φ.  So after L2 the
    perturbed ξ is real to first order.
 L3 [to prove]  With ξ real, the first-order Killing tensor K_ab = (from K_AB K̄_A'B') gives the second Killing field
    ζ = K·ξ; real (ξ, ζ) and the Killing tensor's eigenfunctions (r, cos θ) give Carter canonical coordinates for the
    perturbed metric to first order (each step is a canonical algebraic/ODE construction, so it linearises).
    In Carter form the Λ-vacuum equations force Δ_r, Δ_x to be quartics with linked coefficients; with Λ fixed and the
    coordinate scalings removed, the variations are (δM, δa, δn): the three non-accelerating PD tangents.
 COUNT: 3 (L3) + 1 (acceleration, L2) = 4.  Smooth at every exterior point with a ≠ 0 and sin θ cos θ ≠ 0, because
    every step divides only by Ψ2, τ (∝ a sin θ) and Im ρ (∝ a cos θ).  ⇒ B ∩ {a ≠ 0} = ∅.
CHEAP CHECK for L2 (next): compute Im δξ for the four exact PD tangents (step70) at one point; prediction: zero for
    δM, δa, NUT; nonzero for acceleration.  If NUT also gives Im δξ ≠ 0, L2's bookkeeping changes but not the count.

## 95. L2 VERIFIED: only acceleration makes ξ complex (step95_pd_xi.py, exact nonlinear Λ-PD, 9 s)
Λ-PD in Griffiths–Podolský form (Λ = 1/20, a = 3/5, ε = 1 − Λa²/3, k = a², m = 1), PND-aligned tetrad (κ ~ 1e-42),
Ψ2 from (f), φ = Ψ2^{-1/3}, ξ from item 93(c), at (r, p) = (3, 0.36):
   n = 0,   α = 0:     ξ = (−0.5 − 0.866i) ∂τ                     (real up to a constant)
   n = 1/10, α = 0:     ξ ∝ ∂τ                                      (real up to a constant)
   n = 0,   α = 1/20:   ξ ∝ ∂τ + i(1/20)·(phase) ∂σ, |ξ^σ/ξ^τ| = α, relative phase −π/2   (complex)
   n = 1/10, α = 1/20:  same, |ξ^σ/ξ^τ| = α
So ξ/k = ∂τ ∓ iα ∂σ: Im δ(ξ/k) = ∓δα ∂σ for the acceleration tangent, 0 for δm, δn (and δa, δε via the same family).
L2 [P modulo write-up]: in item 94's L1, B (the ∂φ ∝ ∂σ coefficient) is killed by subtracting δα = ±B times the
  acceleration tangent; A (the ∂t ∝ ∂τ coefficient) is the phase of δk, since ξ/k ∝ ∂τ on the background.  After this
  the perturbed ξ is real to first order.  Matches Hughston–Sommers: reality of ξ ⇔ no acceleration.
Remaining for Theorem B's B-set: L3 (Carter canonical form for the real-ξ perturbation, to first order; count 3).

## 96. THEOREM B LOCAL STEP FOR ALL (σ, m) ≠ (0, 0), a ≠ 0 — BY HAND, NO CERTIFICATE, NO BAD SET
Replaces items 75–87 (jets, descent, Khavkine order, bad set B) for (σ, m) ≠ (0, 0).
Setting: h a (complex) mode of the linearised Λ-vacuum equations on KdS with δΨ0 = δΨ4 = 0 (Theorem A).
(1) Tetrad gauge (item 90): δΨ1 = δΨ3 = 0 by null rotations; then δκ = δσ = δκ' = δσ' = 0.  The perturbed principal
    tetrad is fixed up to boost-spin, so every quantity below is a LINEAR, LOCAL function of h that COMMUTES with
    L_∂t, L_∂φ (the background tetrad is t, φ-independent) and is boost-spin invariant.
(2) ξ := B l + A n − D m − C m̄, A = φρ, B = −φρ', C = φτ, D = −φτ', φ = Ψ2^{−1/3}.  Its 10 Killing equations are
    GHP identities on ANY type D Λ-vacuum (item 93(c), with the tetrad table 93(b)), so they linearise:
       L_ξ h + L_{δξ} g = 0,   δξ = δξ[h].
    On KdS ξ = k ξ_R, ξ_R = −Ξ∂t real (k a fixed constant; φ = r − i a cos θ normalisation).
(3) Reality for free.  J[h] := (δξ[h] − conj(δξ[h̄]))/(2ik) is complex-linear in h; for real h it is Im(δξ/k), and
    taking the imaginary part of (2)/k gives L_{J[h]} g = 0: J[h] is a background Killing field, A'∂t + B'∂φ.
    J commutes with L_∂t, L_∂φ, so for a mode J[h] is a mode of type (σ, m); a constant combination of ∂t, ∂φ is
    annihilated by both, so (σ, m) ≠ (0, 0) forces J[h] = 0.  Hence ξ + δξ is k × (real) to first order, and the
    reality relations  φρ = φ̄ρ̄,  φρ' = φ̄ρ̄',  φτ = −φ̄τ̄'  (normalised φ) hold to first order.
    [No subtraction of PD tangents needed for modes: the acceleration/phase constants of items 94–95 are (0,0)-modes.]
(4) Second Killing field, BY HAND.  Put R = Re φ, I = Im φ.  The reality relations give
       ÞR = −A, Þ'R = B, ðR = ð'R = 0;   ÞI = Þ'I = 0, ðI = iC, ð'I = −iD
    (R, I are the Carter coordinates r, a cos θ).  Define ζ := I²(B l + A n) + R²(D m + C m̄)
    (= K·ξ with K = |φ|²(l⊗n + n⊗l) − R² g).  Its Killing equations, using those of ξ:
       ll, nn, mm, m̄m̄: immediate;
       ln:  sum = −|φ|² φ(ττ̄ − τ'τ̄') = 0, since ττ̄ = τ'τ̄' = −φττ'/φ̄ by reality;
       mm̄:  sum = |φ|² φ(ρρ̄' − ρ'ρ̄) = 0, since φρρ̄' = φ̄ρ̄ρ̄' = φρ'ρ̄ by reality;
       lm:  sum = 2ACφ + |φ|²(ðA + τ̄'A) = 2ACφ − 2ACφ = 0  (ðA = −φρ̄τ, τ̄' = −φτ/φ̄, φρ = φ̄ρ̄ = A);
       lm̄, nm, nm̄: the conjugate/primed computations, identical.
    All algebraic in GHP ⇒ linearises:  L_ζ h + L_{δζ} g = 0, δζ a mode of type (σ, m).
    On KdS (step96_zeta.py, two points): ζ = −0.36216 ∂t − 0.6036 ∂φ = c2 ∂t − Ξa ∂φ  (constant; the R² g term is what
    makes it constant — K = |φ|²(ln+nl) alone is not).
(5) Conclusion.  ξ_R = c0 ∂t with c0 = −Ξ ≠ 0, ζ = c2∂t + c1∂φ with c1 = −Ξa ≠ 0 for a ≠ 0.
      σ ≠ 0:          −iσ c0 k h = L_ξ h = −L_{δξ} g          ⇒ h = L_X g, X = δξ/(iσ c0 k).
      σ = 0, m ≠ 0:   i m c1 h = L_ζ h = −L_{δζ} g            ⇒ h = L_X g, X = −δζ/(i m c1).
    So every mode with (σ, m) ≠ (0, 0) and vanishing Teukolsky scalars is GLOBALLY pure gauge on the exterior, with
    an explicit gauge vector (linear differential expression in h; smooth where Ψ2 ≠ 0, i.e. everywhere outside r = 0).
    No genericity, no open-dense set, no analyticity, no bad set.  HPV's T-frequency-zero modes with m ≠ 0 have
    ∂t-frequency σ = Ω m ≠ 0, so they fall in the first case.
REMAINING for Theorem B: (i) σ = m = 0 (PD sector, items 41, 64; acceleration/NUT excluded by axis regularity);
  (ii) generalised modes (L_T + iσ)^k h = 0 by induction on k; (iii) symbolic check of c1 = −Ξa and of the three
  "identical" Killing components (lm̄, nm, nm̄) — cheap; (iv) referee pass on (1)–(5).

## 97. HOSTILE REFEREE PASS ON ITEM 96 (in-house) — verdict: PASSES, with five repairs now written in
R1 [gap -> fixed] Components nm, nm̄ of ζ's Killing equations are NOT "identical by symmetry" (prime does not map ξ to
   itself).  Done by hand: sum = −2BCφ̄ + |φ|²(ðB + τB), and ðB needs ðρ' (a first derivative NOT fixed by the
   type D equations).  Differentiating the reality relation φρ' = φ̄ρ̄' (legitimate: it holds to first order on an
   open set) with ðφ̄ = −τ̄'φ̄, ðρ̄' = (ρ̄' − ρ')τ̄' gives ð(φρ') = φτρ', so ðB + τB = 2τB and the sum is 0.  lm̄, nm̄ are
   the complex conjugates of lm, nm (ζ real).  All 10 components now verified by hand.
R2 [gap -> fixed] Step (3) needs the LOCAL Killing algebra of KdS (a ≠ 0) to be span{∂t, ∂φ} on every open set.
   Proof: a Killing field preserves Ψ2, hence r and a cos θ (Re/Im of the normalised φ); dr, dx independent, so it is
   tangent to the 2-dim (t, φ)-orbits.  Isotropy at p fixes ∇r, ∇x and must preserve the principal null directions
   (or swap them, impossible in the identity component); l = (timelike part in the orbit plane) + (∂r part), and a
   boost of the orbit plane fixing ∇r, ∇x moves that timelike part, so isotropy is trivial.  dim = 2.
R3 [gap -> fixed] Smoothness of the gauge vector ACROSS the horizons (HPV need ω smooth there).  The Kinnersley
   tetrad is singular at Δ_r = 0, but ξ, ζ (and δξ, δζ) are boost-spin invariant natural constructions from the
   (perturbed) Weyl spinor: K_AB ∝ Ψ2^{−1/3} o_(A ι_B) is fixed covariantly by Ψ_ABCD on type D, ξ ∝ ∇·K.  Computed in
   a horizon-regular principal tetrad (a boost of Kinnersley) they are smooth wherever h is and Ψ2 ≠ 0.
R4 [check done] c1 exactly: ζ = −Ξ(a² ∂t + a ∂φ) on KdS (two geometries: (Λ, a) = (1/20, 3/5) and (3/100, 4/5),
   two points each; ζ_t = −Ξa², ζ_φ = −Ξa to all printed digits) — Carter's Killing-tensor field.  c1 = −Ξa ≠ 0 ⇔ a ≠ 0.
R5 [stated] The commutators used in item 91 ([Þ,Þ'], [ð,ð'], [Þ,ð], [Þ,ð'] on a (0,0) scalar, with κ = σ = Ψ1 = Φ = 0)
   are the standard GHP ones; their consequences were machine-checked on KdS (item 92).  Since the identities are then
   used on PERTURBED metrics, the formulas must be the universal ones, not merely true on KdS: to make this airtight,
   re-check items 91/93 on the accelerating NUT PD metric (item 95's code), a different type D family.
Other objections considered and dismissed:
 * "δξ might not commute with ∂t, ∂φ": every choice (null rotations fixed by δΨ1 = δΨ3 = 0; branch of Ψ2^{−1/3};
   normalisation k) is pointwise and t, φ-equivariant.
 * "J[h] might not be a mode": J is complex-linear (conjugate applied twice) and commutes with real Lie derivatives.
 * "Theorem A gives only δψ0 = δψ4 = 0 for the Teukolsky BCs": item B.1 already supplies this for smooth modes.
 * HPV's T-generalised modes: ∂t-frequency σ = Ωm ≠ 0 for m ≠ 0, covered by the σ ≠ 0 case.
OPEN after the referee: R5 cross-check (cheap); σ = m = 0 sector; generalised-mode induction.
R5 DONE (step97_r5_pd_check.py, 4 s): all nine item-91 relations hold to 2e-42 on the accelerating NUT Λ-PD metric
   (m = 1, n = 1/10, α = 1/20), a type D family different from KdS (ξ genuinely complex there).  The commutator
   consequences are universal, not KdS coincidences.

## 98. σ = m = 0 SECTOR BY HAND: stationary axisymmetric h is (δM, δa, δn, δα) + gauge, for every a ≠ 0
h with L_∂t h = L_∂φ h = 0, δΨ0 = δΨ4 = 0, tetrad gauge of item 90.
(a) J[h] = A'∂t + B'∂φ need not vanish.  Remove A' by a constant phase change of φ̃ (shifts Im δ(ξ/k) by θ·ξ_R ∝ ∂t)
    and B' by subtracting B'/(coefficient) × the acceleration tangent (item 95: its Im δ(ξ/k) is a nonzero multiple
    of ∂σ ∝ ∂φ).  Now reality holds to first order, so ξ̃, ζ̃ are first-order Killing fields (items 93, 96).
(b) Carter structure to first order.  With R = Re φ̃, I = Im φ̃ (normalised):
      ∇R = B l − A n,  ∇I = iD m − iC m̄,  g(∇R, ∇I) = 0,  g(∇R,∇R) = −2AB,  g(∇I,∇I) = −2CD.
    ð(AB|φ|²) = φ³τρ'(φ̄ρ̄ − φρ) = 0 by reality (using ðA = −φρ̄τ, ðB = −φτρ', ð|φ|² = φτ(φ − φ̄)); with ξ- and
    stationarity, d(AB|φ|²) ∝ dR, so  Δr := −2AB|φ|² = Δr(R);  likewise Δx := −2CD|φ|² = Δx(I).
    ξ(R) = ξ(I) = ζ(R) = ζ(I) = 0, ∇R, ∇I are gradients (orthogonal distribution integrable), and
      g(ξ,ξ) = 2AB − 2CD,  g(ξ,ζ) = 2I²AB + 2R²CD,  g(ζ,ζ) = 2I⁴AB − 2R⁴CD
    are all functions of (R, I) through Δr(R), Δx(I), R, I: in coordinates (t, ψ, R, I) adapted to (ξ̃, ζ̃) the
    perturbed metric is, to first order, EXACTLY Carter's separable ansatz with arbitrary Δr(R), Δx(I).
(c) On Carter's ansatz the Λ-vacuum equations are LINEAR in (Δr, Δx) (Carter 1968; Plebański 1975):
      Δr = −(Λ/3)R⁴ + εR² − 2mR + γ,   Δx = −(Λ/3)I⁴ − εI² + 2nI + γ.
    So the linearisation has the same solution space: (δε, δm, δn, δγ); one combination is the joint scaling of
    (R, I) left by the real normalisation of φ̃ (gauge).  3 = (δM, δa, δn).  With (a): + δα.  Total 4 = PD tangents,
    plus a gauge vector built from the coordinate change and the rescalings of ξ, ζ (δξ, δζ are constant combinations
    of ∂t, ∂φ, giving X linear in t, φ: e.g. t ↦ (1 − εc)t).
(d) Axis regularity (item 41/B.3) kills δn and δα.  Result: h = c_M ∂_M g + c_a ∂_a g + L_X g, X polynomial (degree ≤ 1)
    in t — the HV18/HPV σ = 0 form, at every subextremal geometry with a ≠ 0 (a = 0: HV18).
Caveat [to check]: (b) uses reality of FIRST DERIVATIVES (differentiated reality relations); legitimate since the
    relations hold to first order on the whole exterior.

## 99. Generalised modes
σ ≠ 0 or m ≠ 0: induction on k.  If (L_T + iσ)^k h = 0, then h_top = (L_T + iσ)^{k−1}h is a mode with vanishing
Teukolsky scalars (B.1), so h_top = L_{X1} g by item 96, X1 a mode.  Y = t^{k−1}X1/(k−1)! (t-polynomial times mode)
satisfies (L_T + iσ)^{k−1} Y = X1 (L_T = ∂t on the t-polynomial part, ∂t, ∂φ Killing), so h − L_Y g has order k − 1.
⇒ every generalised mode with (σ, m) ≠ (0, 0) is pure gauge, gauge vector polynomial in t.
σ = m = 0, k ≥ 2: the top coefficient is (δM, δa) + gauge by item 98.  OPEN (small): exclude L_T h = c_M ∂_M g + ...
with (c_M, c_a) ≠ 0 (no linear-in-t growth of the mass/spin), as in HV18 for SdS; expected via the conserved
charges (Komar integrals of ξ, ζ), which are t-independent.

## 100. Item 99's open point: no linear-in-t growth of (δM, δa)
Suppose L_T h = v with v = c_M ∂_M g + c_a ∂_a g + L_X g (stationary up to gauge), h = t v + w.
Iyer–Wald: for any solution h and background Killing field K ∈ {ξ_R, ζ}, δH_K[h] = ∮_S (δQ_K − K·Θ(g; h)) is
independent of the 2-sphere S in the exterior (homology between S_0 and S_t = flow of S_0 by ∂t).  On S_t the
integrand is t·(integrand for v) + (t-independent terms, from w and from dt ⊗ v), so constancy in t forces
δH_K[v] = 0 for K = ξ_R, ζ.  The pure-gauge part contributes 0 [standard for X stationary; for X linear in t the
contribution is ∮ Q_{[K,X]}, a Komar integral of a stationary field — CHECK].  For KdS, δH = (δE, δJ) with E = M/Ξ²,
J = aM/Ξ²; ∂(E, J)/∂(M, a) has determinant M/Ξ⁴ ≠ 0, so c_M = c_a = 0.  Then v = L_X g, and h − L_{tX} g is stationary,
and the ladder closes.  [P modulo the one CHECK.]

## 101. THE CHECK (item 100) — RESOLVED: δH_K[L_X g] = 0 for EVERY smooth X, on any closed 2-surface
Iyer–Wald, Λ-vacuum Lagrangian L, background Killing K (L_K g = 0), on-shell background.  With Θ(L_X g) = X·L + dQ_X,
L_X Q_K = X·dQ_K + d(X·Q_K), dQ_K = Θ(L_K g) − K·L = −K·L, and K·dQ_X = L_K Q_X − d(K·Q_X):
   δH_K[L_X g] = ∮ (δ_X Q_K − K·Θ(L_X g)),   δ_X Q_K = L_X Q_K − Q_{[X,K]}
             = ∮ (X·dQ_K − K·X·L) − ∮ (Q_{[X,K]} + L_K Q_X)            (exact terms dropped, S closed)
             = ∮ (K·X·L − K·X·L) − ∮ (Q_{[X,K]} + Q_{[K,X]})             (L_K Q_X = Q_{[K,X]}: Q is natural, L_K g = 0)
             = 0.
No decay or stationarity of X needed.  So item 100 is [P]: no linear-in-t growth of (δM, δa).

## 102. HOSTILE REFEREE PASS ON ITEMS 98–100 (in-house) — verdict: PASSES; three repairs written in
F1 [gap -> fixed] Item 98(b) said "with ξ- and stationarity, d(AB|φ|²) ∝ dR".  Precisely: ∇f = (Þ'f) l + (Þf) n once
   ðf = ð'f = 0, and ∇f ∥ ∇R = Bl − An ⇔ AÞ'f + BÞf = 0 ⇔ ξ̃(f) = 0 (the m, m̄ terms drop).  ξ̃(f̃) = 0 to first order
   because f̃ is a natural invariant of the perturbed geometry and ξ̃ is a first-order Killing field.  So
   df̃ ∧ dR̃ = 0 to first order: f̃ = F(R̃).
F2 [gap -> fixed] Adapted coordinates need [ξ̃, ζ̃] = 0 to first order.  [ξ̃, ζ̃] = ε([δξ, ζ] + [ξ, δζ]).  For stationary
   axisymmetric h, L_ξ h = L_ζ h = 0, so L_{δξ} g = L_{δζ} g = 0: δξ, δζ are background Killing fields, hence commute
   with ξ, ζ.  ✓.  Orthogonal transitivity: ∇R̃, ∇Ĩ are gradients orthogonal to ξ̃, ζ̃ (ξ̃(R̃) = 0 etc.). ✓
F3 [gap -> fixed] Item 99 used BL time t in Y = t^{k−1}X1/(k−1)!; BL t is singular at the horizons.  Use the regular
   time t* (T t* = 1, smooth across both horizons, as in HPV): Y = t*^{k−1}X1/(k−1)! is smooth, and (L_T + iσ)(t*^j X1)
   = j t*^{j−1} X1 as before.  Same in item 100 (S_t := flow of S_0 by T; t* labels it).
Checked and accepted:
 * 98(a): the acceleration tangent has δΨ0 = δΨ4 = 0 (type D) and J ≠ 0 with a nonzero ∂φ part (PD ↔ BL:
   ∂σ = −aΞ∂φ − Ξa²·(1/Ξ)∂t-type combination; ∂φ coefficient −aΞ ≠ 0); its ∂t part is absorbed by the phase of φ̃.
 * 98(b): metric reconstruction g_RR = |φ|²/Δr, g_II = |φ|²/Δx, and the (ξ, ζ)-block from g(ξ,ξ) = (Δx − Δr)/|φ|²,
   g(ξ,ζ) = −(I²Δr + R²Δx)/|φ|², g(ζ,ζ) = (R⁴Δx − I⁴Δr)/|φ|² (from AB = −Δr/(2|φ|²), CD = −Δx/(2|φ|²)): Carter's form.
 * 98(c): linearity of the Λ-vacuum equations in (Δr, Δx) on Carter's ansatz (Plebański 1975) — cited, standard; it is
   the one external input of the σ = m = 0 sector.  [Could be machine-checked in seconds: Ricci of the ansatz with
   symbolic Δr(R), Δx(I).]
 * 98(d): axis determinant −64 m p0⁴/c ≠ 0 (item 41) kills (δn, δα).
 * a ≠ 0 is used (I = a cos θ must be a coordinate); a = 0 is HV18.
 * 101: the identity is the standard covariance of Noether charges; no boundary terms since S is closed.
OPEN: machine-check 98(c) (seconds).

## 103. 98(c) MACHINE-CHECKED AND PROVED (step102_carter_linear.py, 3 s)
Carter ansatz ds² = S(dR²/Δr + dI²/Δx) + (Δx(dt + R²dψ)² − Δr(dt − I²dψ)²)/S, S = R² + I², Δr(R), Δx(I) ARBITRARY.
* Mixed Ricci R^a_b is affine-linear in (Δr, Δx) (second c-derivative under Δ -> cΔ vanishes identically).
* The RR component gives  −S Δr'' + 2RΔr' − 2Δr − 2IΔx' + 2Δx − 2ΛS² = 0;  ∂_I² separates it:
  −2Δr'' − 8ΛR² = const, so Δr = −(Λ/3)R⁴ + εR² + bR + c; the II component gives Δx = −(Λ/3)I⁴ + ε'I² + dI + e;
  substituting back forces ε' = −ε, e = c.  Solution space = {(ε, m, n, γ)}: Δr = −(Λ/3)R⁴ + εR² − 2mR + γ,
  Δx = −(Λ/3)I⁴ − εI² + 2nI + γ, and these solve ALL components (checked).
Since the equations are linear, the linearised solutions are exactly (δε, δm, δn, δγ).  Item 98 has no external input left.

## 104. R3 made ALGEBRAIC: rational cutoff replaces e^{−βκs²} (step104_rational.py, step104_grid.py, ~10 s per run)
Multiplier h = (1+δ) g · w/(1+w) + p(1−u)/(x z),  w = β κ s²,  δ = 0.7  (instead of (1 − e^{−βκs²})).
Same 1200-case grid as item 71 (t ∈ [1.01, 1e5], zf ∈ [0.05, 0.95], B/floor ∈ {1, 10, 1e3}, A/B ∈ [0.01, 1e6]):
   β = 1: 66 fails;  β = 2: 26;  β = 4: 0 (worst D −0.0090);  β = 8: 0 (worst C −0.0823, worst D −0.0178);
   β = 16: 0 (D −0.0040);  β = 40: 33 fails.      Choose β = 8.
Worst cases unchanged in location: C at (t, A/B) -> ∞ (item 73 corner), D at t large, zf = 0.05 (the zf -> 0 corner).
CONSEQUENCE: h, h' are rational in (z; t, z0, A, B), so (C) and (D) are POLYNOMIAL inequalities after clearing the
positive denominators.  The bounded part of R3 becomes an exact polynomial-positivity problem (Bernstein / interval
on polynomials, exact rational arithmetic), and the corners become leading-coefficient sign questions.
NEXT: (i) clear denominators, record the polynomials P_C, P_D and their degrees; (ii) corners: t -> ∞, A -> ∞, zf -> 0
as Newton-polygon leading terms with explicit signs; (iii) bounded part by exact subdivision.

## 105. R3 condition (C) PROVED BY HAND for the rational multiplier (all t, z0, A; B ≥ 5.5t − 1.25)
Polynomial form (step105_polys.py): C has numerator degree 10 (212 terms), D degree 18 (1803 terms); denominators of
known sign (−zL < 0 and a factor −[z0(…) + 8s²z0L(z0)] < 0).
In R3, 2(z0−1)L(1)/(t−1) > 1, so p = ½ min(0.9, ·) = 0.45 exactly.
Proof of (C): s h μ ≤ H.  All three factors of h are ≥ 0 and w/(1+w) < 1, 1 − u ≤ 1, so
   s h μ ≤ s μ (1.7 g + 0.45/(x z)),   x = z − 1,  g = (A/z² + 3)/L.
Divide by μ = z x (t − z) > 0, use μ'/μ = 1/z + 1/x − 1/(t − z):  it suffices that
   1.7 s g + 0.45 s/(xz) ≤ 2 + s/z + s/x − s/(t − z).                                  (*)
 (i)  g − 1/z = (6 − B/z)/L, so g ≤ 1/z when B ≥ 6z; otherwise (only for z > B/6 ≥ (5.5t − 1.25)/6) the excess is
      e ≤ (6 − B/z)/(B − 3z) ≤ (0.5 + 1.25/t)/(2.5t − 1.25).
 (ii) 0.7 s/z + 0.45 s/(xz) − s/x = (s/z)[0.7 − (z − 0.45)/(z − 1)] ≤ −0.3 s/z < 0.
 (iii) s/(t − z) < 1 since z0 < t; and where e ≠ 0, s < t − z < (0.5t + 1.25)/6, so
      1.7 s e ≤ 1.7 (0.5t + 1.25)(0.5 + 1.25/t) / (6(2.5t − 1.25)) ≤ 0.69 (max at t = 1, decreasing).
 Hence LHS(*) − (2 + s/z + s/x − s/(t−z)) ≤ −0.3 s/z + 1 + 0.69 − 2 < 0.   ∎
 [Numerical sanity: the crude inequality has normalised margin ≤ −0.083 on 4e5 random R3 points, step105_crudeC.py.]
 Also holds in R1/R2 wherever p ≤ 0.45 is used.  The t -> ∞ / A/B -> ∞ C-corner of items 71/73 is thereby closed.
REMAINING in R3: condition (D) (degree 18).

## 106. (D) diagnosis: the cutoff-derivative term Gψ' is load-bearing (step106_D_diag.py, 3e5 R3 points, ~5 s)
Exact (D): max normalised margin −0.028 (t ≈ 2e4, zf ≈ 0.12, A/B ≈ 0.014, interior y ≈ 0.7).
Dropping both negative terms ½μGψ' and ½μq': fails (NaN/0 normalisation at z = z0, positive elsewhere).
Dropping only ½μGψ': fails, margin +0.9998 at the seam side (y -> 1, A/B large).
So unlike (C), (D) cannot be proved by discarding signed terms: the seam layer (s² κ = O(1)) needs the exact
balance −sL' ≈ s(A/z² + 3) against ½μGψ' ≈ −1.7 μ g · 8κs/(1+w)², i.e. the rescaled layer problem of items 13–15.
PLAN for (D):  split by w = 8κs²:
  seam zone w ≤ W0: rescale s = σ/√κ; to leading order D becomes the A-independent layer inequality (items 13–15),
     plus explicit remainders bounded by powers of s/x0 and 1/(κ x0²);
  bulk zone w ≥ W0: ψ ≥ W0/(1+W0), and the crude-bounding strategy of item 105 should work (test next, choosing W0).
Bulk test (same points): with w ≥ W0 and BOTH ½μGψ', ½μq' dropped:
   W0 = 1: +0.39 (fails);  W0 = 3: +0.21 (fails);  W0 = 10: −0.058 (holds, 84% of points);  W0 = 30: −0.131.
So take W0 = 10: in the bulk ψ ≥ 10/11 and the crude form
   ½μ'h + ½μ(1.7ψ g' + h²) ≤ sL(1.7ψ − 1)g + sLq,   with g' ≤ g² (L'' ≥ 0),
is the target for a by-hand proof; the bulk lever is w ≥ 10 ⇔ s² L(z0) ≥ 1.25 μ(z0), with L(z) ≥ L(z0) (L' < 0).
The seam zone w ≤ 10 (σ := s√κ ≤ 1.12) carries the remaining difficulty.
Bulk attempt 1 (FAILED, recorded): bounding ½μ'h ≤ ½μh(1/z + 1/x) (dropping −1/(t−z)), ψ ≥ 10/11 in the sL-term,
  ψ < 1 elsewhere, q ≤ 0.45/(xz), g' ≤ g²: normalised margin +0.43 at zf ≈ 0.99, y ≈ 0.92, A/B ≈ 0.01, t ≈ 2e4.
  The −μ/(t−z) part of μ' (z0 near t) is essential; the bulk proof must keep μ' exact and split on its sign.
Bulk attempt 2 (HOLDS numerically, margin −0.092 on all w ≥ 10 points): keep μ' exact and split on its sign:
  (B2)  ½μ'·h± + 0.85 μ g² + ½μ(1.7g + q)²  ≤  0.5455·sLg + sLq     on {w = 8κs² ≥ 10},
        h+ = 1.7g + q (μ' > 0),  h− = (17/11) g + q (μ' < 0),  0.5455 = 1.7·10/11 − 1.
  (B2) implies (D) in the bulk (uses only ψ ∈ [10/11, 1), g' ≤ g², and the two dropped terms being ≤ 0).
  (B2) has no cutoff function left: a rational inequality in (z, t, z0, A, B) plus the side condition
  s² L(z0) ≥ 1.25 μ(z0).  Worst point: t ≈ 1.016, zf ≈ 0.49, y ≈ 0.07 (near z = 1, small t) — not a corner.
B2 split test: (a) g-part  ½μ'c± + 2.295 μg ≤ 0.5455 sL  FAILS (+0.29 near z = 1, t ≈ 1.016: there R3 only gives
  x0 L(1) > (t−1)/2, short of the 1.56(t−1) the g-part alone needs);  (b) q-part  ½μ' + 1.7μg + ½μq ≤ sL  holds (−0.093).
  So near z = 1 the q-term's surplus must pay for the g-part: the proof of B2 has to combine them with a z-dependent
  weight (q ~ 0.45/x dominates as x -> 0).  Not yet found.

## 107. Bulk (D) refined: keep q' (it is cutoff-free) — B3 holds with margin −0.111; first sub-lemma proved
Why B2 was marginal near z = 1: at x -> 0 the 1/x terms are ½μ'q + ½μq² ≈ (0.225 + 0.10125)(t−1)/x, which needs
x0L(1) ≥ 0.725(t−1), more than R3's 0.5(t−1).  The dropped ½μq' ≈ −0.225(t−1)/x cancels the first: with it the
requirement is x0L(1) ≥ 0.225(t−1), which R3 gives with room.
B3 := B2 + ½μq' on the LHS:  numerically max normalised −0.111 on {w ≥ 10} (worst t ≈ 620, zf ≈ 0.71, y ≈ 0.66: interior).
Sub-lemma Q [P]: ½μ(q² + q') ≤ −0.12375 v(t−z)/(xz) − 0.225(t−z)/x0  (v = 1 − u ∈ [0, 1]).
   Proof: q = 0.45v/(xz), q' = 0.45[−1/(x0 xz) − v(2z−1)/(xz)²]; so ½μ(q² + q') = μv(0.2025v − 0.45(2z−1))/(2(xz)²)
   − 0.225μ/(x0 xz), and 0.2025v − 0.45(2z−1) ≤ 0.2025 − 0.45 = −0.2475 (v ≤ 1, z ≥ 1); μ = xz(t−z).  ∎
Also ½μ'q = 0.225 v(t−z)(1/x + 1/z) − 0.225v exactly (μ'/μ = 1/x + 1/z − 1/(t−z), μq = 0.45v(t−z)).
Remaining for B3: the g-terms  ½μ'·c±g + 2.295μg² + 1.7μgq  against  0.5455 sLg + (sLq − q-terms), using the bulk
   lower bound s ≥ √(1.25/κ) together with L(z) ≥ L(z0) (so sL ≥ max(sL0, √(1.25 μ0 L0))).  Case split on z vs z0.
B3 SPLIT (both hold numerically on {w ≥ 10}; with Sub-lemma Q they imply B3, hence bulk (D)):
  (Iq)  (t−z)(0.225(2z−1) − 0.12375) + 0.765 μg  ≤  0.45 sL                         margin −0.036
  (Ig)  ½μ'c± + 2.295 μg  ≤  0.5455 sL + (0.225v + 0.225(t−z)/x0)/g                 margin −0.106
  (Iq) is divided by v/(xz); (Ig) by g.  Both are now inequalities between explicit rational functions, no cutoff.
  Tightest regime of (Iq): zf small (z0 = ζt, ζ ≲ 0.12), A small, z near z0: needs s ≳ 0.49 z0 against the bulk bound
  s ≥ √(1.25/κ) ≈ 0.48 z0 there — a 2% margin; the hand proof must use g ≈ 3/B (not g ≤ 1/z) in that regime.

## 108. (Ig): which terms are load-bearing, and a monotonicity reduction
(Ig) with g·L = N := A/z² + 3:   F := ½c±μ' + 2.295 μN/L − 0.5455 sL − 0.225(t−z)L/(x0N) ≤ 0.
 * Exact μ' (both signs, c+ = 1.7, c− = 17/11): holds, margin −0.103.
 * Dropping the μ' < 0 part: FAILS (+0.17, zf ≈ 0.99: z0 near t, where μ' < 0 pays).
 * Also dropping the slack 0.225(t−z)L/(x0N): FAILS (+0.29, z near 1, t ≈ 1.016).
 * μ' ≤ (t−z)(2z−1): FAILS (+0.43).  So the proof must keep μ' exact and the slack term.
REDUCTION [P]: at fixed (t, z0, A, z), F is strictly DECREASING in B (L increases with B; μN/L decreases; sL and
 (t−z)L/(x0N) increase), and the constraints (B ≥ 5.5t − 1.25; R3: 2x0L(1) > t−1; bulk w = 8L(z0)s²/μ(z0) ≥ 10) are
 all monotone (preserved as B increases).  Hence it suffices to prove F ≤ 0 on the minimal-B boundary:
 B = max(5.5t − 1.25, B_R3, B_w), i.e. either B at the floor, or the point lies on the R3 boundary (handled by R2's
 multiplier side too), or on w = 10 (shared with the seam zone).  One variable eliminated; same holds for (Iq).
A-monotonicity test for F on the floor slice B = 5.5t − 1.25 (step108_A_monotone.py, 2e4 (t, z0, z) lines, 200 A-values
 each): decreasing in A on 93%, interior maximum on 7%, never increasing.  So no clean A-reduction; the worst A is
 interior on a small set.  (d/dA of N/L has the sign of B − 6z, which changes only for z > B/6 ≈ 0.92t.)
(Ig) attempt G1 (FAILED): on z ≤ B/6 use g ≤ 1/z in 2.295μg and in the slack, and ½c±μ' ≤ 0.85(t−z)(2z−1) − 0.7727zx:
  normalised +0.21 at t ≈ 35, zf ≈ 0.73, A/B ≈ 0.09.  For small A/B, g ≈ (A/z² + 3)/B ≪ 1/z, so g ≤ 1/z throws away
  the main gain; the proof must work in the two regimes A ≪ Bz (g ≈ N/B, L ≈ B) and A ≫ Bz (g ≈ 1/z, L ≈ A/z)
  separately, with the crossover A ~ Bz handled explicitly.

## 109. (Ig) on z < B/6 reduced to a 3-variable problem by an EXACT identity + concavity [P]
Identity: N = A/z² + 3 = (L − B + 6z)/z, so  g = N/L = 1/z − d/(zL),  d := B − 6z  (> 0 for z < B/6).
At fixed (t, z0, z) (B at the floor by item 108), A ↔ L is affine, and
   F(L) = K − a/L − bL − c·d/(L − d),   with  K = ½c±μ' + 2.295 μ/z − c,  a = 2.295 μ d/z,  b = 0.5455 s,
                                              c = 0.225 z(t − z)/x0      (all ≥ 0 except K).
Each of −a/L, −bL, −cd/(L − d) is concave on L > d, so F is CONCAVE in L, and F -> −∞ as L -> ∞.  The constraints
(A ≥ 0, both R3 inequalities, bulk w ≥ 10) are each L ≥ (explicit bound) since L(1), L0 are affine increasing in A.
Hence  sup F = max(F(L_min), F(L_crit)) with L_crit the unique stationary point if L_crit > L_min.
Numerics (step109_Ig_concave.py, 3e5 (t, z0, z)):
   boundary F(L_min) ≤ 0: holds, worst −0.057 (2.3e5 points);
   interior with −a/L dropped (L* = d + √(cd/b), bound K − 2√(bcd) − bd): FAILS +0.18 at zf ≈ 0.02, t ≈ 5e4 —
     the −a/L term is needed at the interior max; next: keep it (F'(L) = a/L² − b − cd/(L − d)² = 0, a cubic-type
     equation; or bound a/L from below using L_crit ≤ d + √((cd + a)/b)).
So (Ig) on 84% of the bulk is now a statement about explicit functions of (t, z0, z) only.
UPDATE: at the stationary point F'(L) = 0 forces b < a/L², i.e. L_crit < √(a/b).  Numerically √(a/b) < L_min at ALL
  227,437 sampled (t, z0, z) with z < B/6: the stationary point is never feasible, so F is decreasing on the feasible
  set and  sup F = F(L_min).  Hence (Ig) on z < B/6 ⇐ (J1) √(a/b) < L_min  and  (J2) F(L_min) ≤ 0 (margin −0.057).
  (J1) reads 4.207 x(t − z)(B − 6z)/s < L_min², with L_min ≥ max(B − 3z, L from w = 10, L from R3).
  Both are explicit inequalities in (t, z0, z).  To prove next.
(J1) needs only A ≥ 0 and the bulk condition (R3 inactive): margin −0.090 on 4e5 points (step109_J1.py); worst at
  zf ≈ 0.009 (z0 ≪ t), y ≈ 0.52.  The crude form s·L_min > 4.207 x(t − z) (dropping d/L ≤ 1) FAILS (+0.23 near
  z0 ≈ t, z ≈ 0.91 t, where d = B − 6z ≈ 0.1t ≪ L): the factor d must be kept.  Useful exact facts:
  L ≥ B − 3z ≥ d;  L − L0 = A(1/z − 1/z0) + 3s ≥ 3s;  on the bulk boundary L0 = 1.25 μ0/s².
(J1) sharpened: since L ≥ B − 3z, L² ≥ (B − 3z)·L, so (J1) ⇐ (J1')
     4.207 x(t − z)·d/(B − 3z)  ≤  s·L,   with  s·L ≥ max(s(B − 3z), 1.25 μ0/s)   (bulk; L ≥ L0).
  (J1') with that max: holds, margin −0.072 (4e5 points; worst zf ≈ 0.32, y ≈ 0.59, t ≈ 6e4).
  Replacing the max by the geometric mean √(1.25 μ0 (B − 3z)): FAILS (z0 -> t makes μ0 -> 0).  So the proof is a
  two-case split in s: use s(B − 3z) when s ≥ 4.207 x(t−z)d/(B−3z)², and 1.25μ0/s otherwise; need: the two
  thresholds never leave a gap (an inequality in (t, z0, z) after z = z0 − s, x = x0 − s).

## 110. Certificate route for (D): feasibility findings (8 Nov, evening)
* Naive 5-D interval branch-and-bound on the full (D) (cert2/iv.py, cert2/dcert.py; outward-rounded numpy intervals):
  on the small box t ∈ [2,3], zf ∈ [0.2,0.8], y ∈ [0.02,0.98], A/B ∈ [0.5,50], B/floor ∈ [1,3]: 5.7e5 boxes closed,
  3.3e6 open after 12 s.  Dependency blow-up; NOT viable under the 1-minute rule.  Symbolic expansion of P_D in compact
  variables also times out (> 60 s).
* Reduced pieces (no cutoff, B at floor, one L-variable):
  - (Ig) boundary F(L_min) ≤ 0 with L_min from A ≥ 0 and the bulk only (superset of R3): holds, −0.072, BOTH strips.
    For z ≥ B/6 (d < 0) F is CONVEX in L and -> −∞, so sup = F(L_min) there too [P].
  - (Iq) boundary: FAILS (+0.067 at t ≈ 7, zf ≈ 0.9993, y ≈ 0.978, A = 0) even with the R3 constraints — but the full
    (D) holds there (−0.17).  So the (Iq)/(Ig) split loses too much in that corner; the certificate must target B3
    (the unsplit bulk inequality), which is a function of the single variable L at fixed (t, z0, z).
  - (J1q) stationary point of (Iq) infeasible: −0.26.
PLAN: certify B3 over (T = 1/t, zf, y, λ) with λ compactifying L ∈ [L_min, ∞) — 4-D, short expression; then the seam.

## 111. B3 certificate attempts (cert2/b3cert.py naive intervals; cert2/di.py + b3mv.py mean-value form)
Box t ∈ [2,3], zf ∈ [0.02,0.98], y ∈ [0.02,0.98], B/floor ∈ [1,3], (B−3z)/L ∈ [1e-3,1]: naive intervals 2.1e6 closed /
3.7e6 open in 22 s.  Mean-value form (forward AD over intervals), smaller box t ∈ [2,3], zf ∈ [0.2,0.8], y ∈ [0.1,0.9],
B/floor ∈ [1,3], λ ∈ [0.1,1]: 5.4e5 closed / 6.1e5 open in 56 s.  Diagnostic (cert2/diag2.py): after 18 bisections,
1.3e5 open boxes, 8.9e4 with centres in the domain, ZERO with F > 0 at the centre: no counterexample, the boxes are just
too wide for the 5–10 % margin.  Open boxes are spread over all five variables (no single bad corner).
Conclusion: a 5-variable branch-and-bound at this margin needs ~1e6–1e7 boxes per unit cell — tens of minutes in total.
Not compatible with the 1-minute rule without first removing variables analytically (B, then L) or a tighter
enclosure (Taylor models).

## 112. Option 1: removing B from B3 — convexity in g [P] and where the worst cases live
(a) B3 is NOT monotone in B at fixed (t, z0, z, A) (20k lines: 17298 decreasing, 58 increasing, 439 mixed/261 interior
    max) — the μ' < 0 term and −0.5455 sLg pull the other way.  So item 108's reduction does not transfer to B3.
(b) KEY [P]: at fixed (t, z0, z), B3 depends on (A, B) ONLY through L and g.  At fixed L, moving along A = z(L − B + 3z)
    changes g = 1/z − (B − 6z)/(zL) AFFINELY in B, and B3 = 2.295 μ g² + (affine in g) + const (μ' sign fixed by z):
    B3 is CONVEX in B along each L-level set.  The constraints along the level set are: bulk (L0 increasing in B) gives
    B ≥ B_w; floor B ≥ 5.5t − 1.25; A ≥ 0 gives B ≤ L + 3z; R3 (L(1), A/z0 decreasing in B) give B ≤ B_R3.  So the
    feasible B-set is an interval and the max of B3 is at an END: on the floor, on the bulk boundary w = 10, on A = 0,
    or on an R3 boundary.  (This is the analogue of item 108, now valid for B3.)
(c) Numerics on those surfaces (step112_surfaces.py):
    floor B = 5.5t − 1.25, vary A: max at A = A_min in 19932/19934 lines (2 interior);
    bulk surface w = 10, vary B: max at an END in 4945/4945 (127 at floor, 4818 at A = 0);
    A = 0 surface, vary B: max at the floor in 10227, interior in 316, top in 1.
    => the certificates needed are:  C1 the corner (floor, A_min): 3 variables (t, z0, z);
       C2 the A = 0 surface: 4 variables (t, z0, z, B) but a SHORT expression (g = 3/L, L = B − 3z);
       C3 the R3 boundary surfaces (not yet examined);
       plus a proof that on the floor slice the max is at A_min (2 exceptional lines to look at).
(d) FIRST CELL CERTIFIED: C2 (A = 0 surface) on t ∈ [2,3], zf ∈ [0.05,0.95], y ∈ [0.05,0.95], B/floor ∈ [1,10]:
    certified, 3.05e5 boxes, 16 s (cert2/b3mv.py with λ ≡ 1; mean-value form, outward-rounded intervals).
    The reduction makes the certificate feasible: one unit cell per run, well under a minute.
REMAINING to finish (D)-bulk: tile C2 over t ∈ (1, T∞], zf, y ∈ (0,1), B/floor ∈ [1, ∞) (edges by limit problems);
    C1 (corner, 3-D); C3 (R3 boundary surfaces); the floor-slice A_min claim.  Then the seam zone.
(e) SOUNDNESS FIX in cert2/iv.py: inf·0 inside a product used to be replaced by 0 (unsound when an inverse straddles 0);
    now any NaN in the four products gives (−∞, ∞).  The (d) cell re-certifies after the fix (3.05e5 boxes, 20 s).
    Probes (A = 0 surface, B/floor ∈ [1,10], zf, y ∈ [0.05,0.95] unless stated):
      t ∈ [1.01,1.1]: certified, 2.1e3 boxes, 0.4 s;     zf ∈ [0.005,0.05]: certified 3.3e3, 0.3 s;
      y ∈ [0.005,0.05]: certified 5.7e3, 0.4 s;          y ∈ [0.95,0.995]: certified 5.7e3, 0.6 s;
      t ∈ [5,10]: NOT certified in 45 s (5.0e5 closed, 2.3e5 open);  t ∈ [100,150]: not certified in 45 s.
    So the small-t cells and the zf, y edges are cheap; LARGE t is the cost — the enclosure must use t-scaled variables
    (z = tζ etc., the item 73 outer variables) and the t -> ∞ limit has to be certified as its own (compact) problem.

## 113. A = 0 surface on a COMPACT domain: exact Bernstein branch-and-bound (cert2/a0_scaled.py, a0_bern.py)
* t·B3 on A = 0 is, after T = 1/t, bt = floor/B, an explicit rational function of (T, ζ, y, bt) ∈ [0,1]⁴ whose
  denominator has a fixed sign (−) — so B3 ≤ 0 ⇔ num ≥ 0, num a polynomial of degree (4,5,5,3), 168 terms, for each
  sign branch of μ'.  The constraints (R3 c2, R3 pp, bulk w ≥ 10, and μ' for the branch) are polynomials too.
  This covers t ∈ (1, ∞) and B ∈ [floor, ∞) at once.  Exact integer Bernstein coefficients, de Casteljau splits.
* Run (48 s): 4.7e4 boxes closed; does not finish.  All stuck boxes sit at one corner: T -> 0, bt -> 0, y -> 1
  (t -> ∞, B -> ∞, z -> z0) — num vanishes on the edges {T = bt = 0} and {bt = 0, y = 1}.
* Blow-up at that corner (corner.py): T = aε, bt = bε, 1 − y = cε gives num = ε²·(−264 b ζ)(1650a(ζ − 1)
  − 10200 bζ³ + 6800 bζ² − 12001 cζ²) + O(ε³).  The b² term has the WRONG sign for ζ < 2/3, so the polynomial alone is
  not ≥ 0 there; it is the bulk constraint that saves it: at the corner w ≈ 44c²ε/(b(1 − ζ)), so w ≥ 10 forces
  b ≤ 4.4 c²ε/(1 − ζ) = O(ε), the b² term drops, and the leading form is −264bζ(1650a(ζ−1) − 12001cζ²) > 0.
  Bernstein cannot see this because the bulk polynomial also vanishes at the corner.
NEXT: re-parametrise the corner using the bulk variable (replace bt by w, or by b/(c²ε)) so the constraint becomes a
  box side; then Bernstein closes it.  Same treatment will be needed at the other degenerate corners.
(113, cont.) Re-parametrised with ω = 10/w ∈ (0,1] (bulk constraint = box side; cert2/a0w.py, a0w_run.py, bern_core.py):
  the (T, bt -> 0, y -> 1) corner is gone.  New stuck corner: ζ -> 1, y -> 1 (z0 -> t, z -> z0), all T, ω.
  Blow-up ζ = 1 − aε, y = 1 − cε: num = 150a·[T a²(150ω + 375) + acω(1091 − 791T) + 1360c²ω²(1 − T)]·ε³ + …  (> 0), but
  on the face ζ = 1 itself num ≈ −832320 c⁵ ω³ (1 − T)² ε⁵ < 0, so num < 0 where 1 − ζ ≲ (1 − y)³.
  Checked: those points have B ≈ 3z0 + 10μ0/(8s²) < floor (e.g. t = 7/6: B = 3.50 vs floor 5.17), i.e. the FLOOR
  constraint excludes them — again a constraint whose polynomial vanishes at the corner, invisible to Bernstein.
  Pattern: every stuck corner so far is a point where a constraint and num vanish together; the fix each time is to
  blow up so the binding constraint becomes a box side.  Next: blow up (ζ, y) -> (1, 1) with 1 − ζ = σ(1 − y)².

## 114. C2 CERTIFIED: B3 ≤ 0 on the whole A = 0 surface of the bulk, all t ∈ (1, ∞), B ≥ floor (exact, 5 s)
cert2/run_C2_all.sh — exact integer Bernstein branch-and-bound, three charts covering the compact parameter box:
   chart bt = floor/B  on T = 1/t ∈ [0.05, 1], all (ζ, y, bt):      certified, 349 boxes, 0.1 s
   chart ω = 10/w      on T ∈ [0, 0.05], ζ ∈ [0, 0.9], all (y, ω):  certified, 6504 boxes, 4.7 s
   chart bt            on T ∈ [0, 0.05], ζ ∈ [0.9, 1]:              certified, 8 boxes, 0.0 s
 Rule per box: some constraint (R3 c2, R3 pp, bulk, floor) is certainly violated, or for each sign branch of μ' that
 can occur, the B3 numerator has the certified sign.  Denominator signs fixed (factors of fixed sign / squares).
 Each chart avoids the other chart's degenerate corner (bt-chart: T, bt -> 0, y -> 1; ω-chart: ζ, y -> 1).
Floor surface (B = floor, A free) — in progress (cert2/surf_build.py, surf_run.py):
   chart A = floor·p/(1−p): stuck at (y, p) -> (1, 1) (bulk degenerate) and at T -> 1, ζ -> 0 (the known zf -> 0
   degenerate corner, margin ~ −Kζ: needs the ζ factor divided out);
   chart A from w = 10/ω: stuck at (ζ, y, ω) -> 1 (floor degenerate).
Remaining surfaces: floor (finish), bulk w = 10, R3 boundaries; then the seam zone.

## 115. Floor surface CERTIFIED; generic composition builder; status of the other surfaces
* FLOOR SURFACE (B = floor, all A) CERTIFIED: cert2/surf_driver.py with three charts (A = floor·p/(1−p); A from
  w = 10/p; A from the R3 c2 parameter), each a superset of the feasible set.  46 dyadic (T, ζ, y)-cells, each with one
  certifying chart, logged in cert2/cert_floor.log; cert2/check_tiling.py: the 46 cells are disjoint and have total
  volume exactly 1 (exact rationals) — they tile the cube.  Total run time ≈ 3.5 min in four resumable < 1 min pieces.
* gen_build.py / surf2.py: B3 numerator as a generic polynomial in (A, B, T, ζ, y) (2 s), surfaces by exact rational
  substitution + homogenisation with sp.Poly arithmetic (< 1 s per surface) — replaces symbolic `together` (> 60 s).
  Generic denominator −4000·T y ζ (T − 1)(…)²(…)² > 0, so B3 ≤ 0 ⇔ numerator ≤ 0.
* Bulk surface w = 10 (chart B = floor/p): numerator degree (20, 11, 14, 6); first 45 s: no cell closed — the hardest
  surface (B3's margin is smallest on the bulk boundary, and degree is high).  Item 112(c) numerics say its max over B is
  always at an END (floor or A = 0, both certified); a PROOF of that would remove this surface entirely.
* R3 boundary surfaces: pp = 1 (chart B = floor/p) and c2 = 0 — running (c2: 20 cells closed in 41 s, 4 cells hit the
  depth cap and need another chart).

## 116. Bulk surface w = 10: reduced to a ONE-variable no-interior-maximum statement [P reduction + N]
On w = 10 at fixed (t, z0, z): A = z0(L0* − B + 3z0), L0* = 1.25 μ0/s², so L = L1 − κB (κ = s/z, L1 = (z0/z)(L0* + 3z0) − 3z)
and g = g0 − C u with u = 1/L, g0 = (1 + 1/κ)/z, C = (L1/κ − 6z)/z.  Writing B3 = 2.295μg² + G1 g + G0 with
G1 = ½μ'c± + 1.7μq − 0.5455 sL, G0 = ½μ'q + ½μq² + ½μq' − sLq, one gets EXACTLY
   B3 = f(u) := a2 u² + a1 u + a0 − |E|/u,
   a2 = 2.295μC²,  a1 = −C(4.59μg0 + G1c),  a0 = 2.295μg0² + G1c g0 + 0.5455 sC + G0c,  |E| = s(0.5455 g0 + q),
with G1c = ½μ'c± + 1.7μq, G0c = ½μ'q + ½μq² + ½μq' (u-free).  The feasible B-set on the surface is an interval
[floor, min(B_A=0, B_pp, B_c2)] (pp, c2, A decrease in B along the surface), i.e. an interval of u.
Numerics (cert2/bulk_end_check.py, 6978 lines): max over the interval ALWAYS at an end (0 interior), f convex on 99.3%.
The ends lie on the floor, A = 0 (both certified, items 114–115), or on an R3 surface (in progress).
REMAINING for this surface: prove f has no interior maximum, i.e. f'(u) ≠ 0 wherever f''(u) = 2a2 − 2|E|/u³ ≤ 0 —
  a polynomial statement 2a2u³ + a1u² + |E| > 0 on {u ≤ (|E|/a2)^{1/3}} ∩ interval, low degree in u.
(116, cont.) THREE-VARIABLE FORM [P]: f'' = 2a2 − 2|E|/u³ changes sign once, at u_c = (|E|/a2)^{1/3}; f' is decreasing
  on the concave part u ≤ u_c, so min f' there is ≥ f'(u_c) = 3a2u_c + a1.  Hence f has no interior maximum on the
  feasible u-interval [u_lo, u_hi] (u_lo = 1/L(floor)) provided
     (C*)  whenever |E| > a2 u_lo³ (concave part meets the interval):  a1 ≥ 0  or  27 a2² |E| + a1³ > 0,
  a statement in (t, z0, z) only, independent of which of A = 0 / pp / c2 bounds the interval above.
  Numerics (cert2/bulk_fu.py): f' > 0 on every concave part, 13913 lines, worst normalised 0.015; (C*) itself: only
  779 of 1e5 sampled points have |E| > a2u_lo³; all satisfy (C*), worst normalised 0.0073 — thin margin.
  Next: build (C*) symbolically in (T, ζ, y) and certify by exact Bernstein.
(116, cont.) Building (C*) symbolically FAILED within 60 s: sympy cancel/together times out; a no-cancellation
  rational-polynomial builder (cert2/ratpoly.py, cstar_build.py) gives numerator degree (145, 57, 57) — far too large
  for Bernstein.  The common factors (powers of T, ζ, y, 1 − y, 1 − T, z, z0, t − z0) must be cancelled by hand before
  composing; or (C*) certified by intervals in (t, z0, z) with an explicit large-t limit.  Not done.

## 117. ROUTE AROUND (D): saturate (C), then solve (D) with equality (Hintz 2608.25868 §4.3.3 structure)
In f-language (h = f'/f), (D) ⇔ D̃f := (μf')' − 2s(Lf)' ≤ 0, and (C) ⇔ h ≤ h_C := H/(sμ).
Construction (no ad hoc multiplier, no certificate of (D)):
  Phase 1: f = F_C := exp∫h_C on (1, ẑ] — (C) holds with EQUALITY; (D) ⇔ d(z) := D[h_C] ≤ 0.  F_C ~ (z−1) at z = 1,
           so f(1) = 0 automatically.
  Phase 2: on [ẑ, z0], f2 solves D̃f2 = 0 — (D) holds with EQUALITY — with (f2, f2') = (F_C, F_C') at ẑ; since
           d(ẑ) = 0 the second derivatives match (f is C²).
  (C) on phase 2 [P, Wronskian]: with p = μ exp(−∫2sL/μ) > 0, D̃f = (μ/p)(pf')' − 2sL'f, so W = p(F_C f2' − F_C' f2)
           has W' = −(p/μ) f2 D̃F_C ≤ 0 wherever d ≥ 0 and f2 > 0; W(ẑ) = 0 ⇒ f2'/f2 ≤ h_C on [ẑ, z0).
  What must be PROVED instead of (D):
    (P1) d has exactly one sign change on (1, z0), from − to +  (an explicit rational function of z);
    (P2) f2 > 0 on [ẑ, z0] and the seam condition f2'(z0) > 0  (comparison with a second explicit function).
Numerics (step117_hintz_route.py; 365 random R3 parameter points, t up to 1e4, A/B up to 1e5, B up to 100×floor):
  d < 0 near z = 1 and exactly one sign change: 365/365;  phase 2 positive, (C) holds, seam holds: 365/365.
  (A first run with a finite-difference h_C' mislocated ẑ and showed spurious (C)-violations; with the exact h_C'
   and bisection for ẑ, 0 violations, exactly as the Wronskian argument predicts.)

## 118. (P1) — the sign change of d — reduced to finitely many TWO-variable (t, z0) inequalities
(a) [P] h_C = H/(sμ) = 2/s + μ'/μ = (log(μ/s²))', so the phase-1 function is simply F_C = μ/s².
(b) [P, step118_P1.py] D̃F_C = N(z)/s⁴ with N of degree 6 in z.  In S = z0 − z ∈ (0, x0):
     N = 6μ0² − 6μ0μ0'·S + c2 S² + c3 S³ + 0·S⁴ + 2(B + 4t − 15z0 + 4) S⁵ + 15 S⁶,
     c2 = −4μ0(A/z0 + B)·… (exactly 4(z0−1)(z0−t)(A + Bz0) + R2(t, z0)),  c3 = 2[A(t+1−2z0) + Bμ0' + R3(t, z0)].
     N(0) = 6μ0² > 0 (d > 0 at the seam).  (P1) ⇔ N has exactly one root in (0, x0).
(c) [N] Möbius S = x0 v/(1+v) and Descartes: the 7 transformed coefficients T0..T6 show EXACTLY ONE sign change on
     all 92,444 sampled R3 points (step118_P1_numerics.py; plain root-finding also gives exactly one root, 0 exceptions).
     Pattern: T0, T1 > 0 always; T4, T5, T6 < 0 always; T2, T3 mostly < 0 (T2 > 0 on 3.8%, T3 > 0 on 0.06%).
(d) [P] T0 = 6μ0² > 0;  T1 = 6z0(z0−t)x0²(3z0² + 2z0 − 4tz0 − t) > 0 (both factors negative: the last is convex in z0
     and negative at z0 = 1, t).  For j = 2..6, T_j is AFFINE in (A, B) with both coefficients < 0 (by hand: e.g.
     T4: 6x0³(−3t + 2z0 + 1), 6x0³(z0² + 2z0 − 2tz0 − t); T5: 2x0³(−5t + 2z0 + 3), −2x0³(2tz0 + 3t − 4z0 − 1); T6:
     −2(t−1)x0³ twice; T3: 2x0³(−7t + 6z0 + 1), 2x0³(5z0² + 2z0 − 6tz0 − t); T2: 4(z0−t)x0³, 4z0(z0−t)x0³).
     So each T_j attains its sup over the feasible (A, B)-region at a lower-left vertex of that region (floor, R3 c2, R3 pp
     lines): two-variable functions of (t, z0).
(e) Remaining: (i) T4, T5, T6 < 0 at those vertices — A = 0 alone is NOT enough (sign changes seen), the R3 lower bounds
     on A are needed; (ii) "T3 ≥ 0 ⇒ T2 > 0" on the feasible region (the λ-combination with λ = 2(t−z0)/(7t−6z0−1)
     kills A but is DECREASING in B, so it must be evaluated where T3 = 0 meets the region).  All are 2-D statements:
     exact Bernstein in (T = 1/t, ζ) will be cheap.

## 119. (P1) PROVED (hand reductions + exact 2-D certificate, cert2/run_P1_all.sh, 8/8 certified, < 1 s)
Feasible (A, B) in R3 at fixed (t, z0): B ≥ Bf = 5.5t − 1.25, A ≥ a0 := max(0, Ac2), Ac2 = z0((t−1)/(2x0²) − 3), and the
pp line A + B ≥ P := 3 + (t−1)/(2x0).  Certified: a0 ≥ P − Bf (G0), so pp is REDUNDANT and the feasible set is the
quadrant {A ≥ a0, B ≥ Bf}.  Since T2..T6 are affine with negative A- and B-coefficients (item 118(d)):
  * T4, T5, T6 < 0 on the quadrant ⇔ at its corner (a0, Bf): certified (al4, al5, al6 < 0), both cases a0 = 0 / a0 = Ac2;
  * {T3 ≥ 0} ∩ quadrant is a triangle with vertices (a0,Bf), (a0 + α3/|β3|, Bf), (a0, Bf + α3/|γ3|); T2 > 0 there ⇔
    α2 > 0, α2|β3| − |β2|α3 > 0, α2|γ3| − |γ2|α3 > 0 whenever α3 ≥ 0: certified (al2, v2, v3).
  Hence the Möbius-transformed coefficients of N have exactly one sign change: (+, +, ±, ·, −, −, −) with
  T3 ≥ 0 ⇒ T2 > 0.  By Descartes, N has exactly one root on (0, x0); with N(0) = 6μ0² > 0 and the sign at S = x0
  forced by the single root, d < 0 near z = 1 and d > 0 near z0 with exactly one sign change.  ∎ (P1)
Certificate: step119_build.py builds the 17 two-variable polynomials in (T = 1/t, ζ) (degrees ≤ 11); p1chart.py runs
exact integer Bernstein on four charts — two rectangles and the blow-up of the corner (T, ζ) -> (0, 0) by ζ = Tw and
T = ζw (the corner is where the case gate Ac2 and the targets vanish together); factors u, w, 1 − u, 1 − w that are
positive on the open box are divided out first.  8 runs (2 cases × 4 charts), all certified.
Remaining for the route of item 117: (P2) f2 > 0 on [ẑ, z0] and the seam f2'(z0) > 0.

## 120. (P2): reduced to the seam sign alone [P]; first two proof attempts fail [N]
(a) [P] On phase 2, D̃f2 = 0 reads μf'' = −(μ' − 2sL)f' + 2sL'f.  At any zero of f' with f > 0, μf'' = 2sL'f < 0, so f'
    can only cross 0 downward and never return.  f2'(ẑ) = F_C'(ẑ) = H(ẑ)/ŝ³ > 0 (Lemma H).  Hence
       f2'(z0) > 0  ⇒  f2' > 0 on [ẑ, z0]  ⇒  f2 > 0 and increasing:  (P2) ⇔ the seam sign f2'(z0) > 0.
(b) [P] Φ := μf' − 2sLf satisfies Φ' = 2Lf (from D̃f = 0), and Φ(z0) = μ0 f2'(z0).  So Φ(ẑ) ≥ 0 would suffice;
    Φ(ẑ) = f̂ (Ĥ/ŝ − 2ŝL̂), i.e. J := H − 2s²L > 0 at ẑ.
    [N] FALSE in general: J(ẑ) > 0 on only 3816 of 37019 sampled R3 points (step120_P2_phi.py).  The growth term
    2∫Lf must be used.
(c) [N] Barriers for the phase-2 Riccati h = f2'/f2 (a function k with D[k] ≤ 0, k ≤ h(ẑ) gives h ≥ k):
    numerically h ≥ 5.8 g on phase 2 (g = −L'/L) in all 543 cases (step120_P2_explore.py), but none of k = g, 2g, 3g,
    0.3h_C, 0.5h_C, 0.7h_C is a subsolution on all of phase 2 (each fails in ≥ 60 % of cases, near the seam).
NEXT: a barrier adapted to the seam (k with the right behaviour as s -> 0 — e.g. matching the Frobenius data of the
  phase-2 ODE at z0, where it is regular), or an integral argument using Φ' = 2Lf with a lower bound f ≥ f̂ e^{∫k}.

## 121. (P2) PROVED — the seam sign; hence the event lemma in R3 is PROVED (route of item 117)
(a) [P] Bootstrap.  Suppose f2'(z0) ≤ 0; let z1 ∈ (ẑ, z0] be the first zero of f2' (f2'(ẑ) > 0).  On [ẑ, z1] f2 is
    increasing, so f2 ≥ f̂ := f2(ẑ) > 0, and Φ = μf' − 2sLf has Φ' = 2Lf ≥ 2Lf̂.  Hence
       −2s1L1f2(z1) = Φ(z1) ≥ Φ(ẑ) + 2f̂∫_ẑ^{z1} L = f̂(Ĥ/ŝ − 2ŝL̂ + 2∫_ẑ^{z1} L),
    and since f2(z1) ≥ f̂, L1 > 0:  Ψ(z1) := Ĥ/ŝ − 2ŝL̂ + 2∫_ẑ^{z1} L + 2s1L1 ≤ 0.  But Ψ'(z1) = 2s1L'(z1) < 0, so
    Ψ(z1) ≥ Ψ(z0) = Ĥ/ŝ − 2∫_ẑ^{z0} s|L'| dz  (by parts).  So it suffices that  Q(ẑ) := Ĥ/ŝ − 2∫_ẑ^{z0} s|L'| > 0.
(b) [P] ∫_z^{z0} s|L'| = A(z0/z − 1 − log(z0/z)) + (3/2)s² ≤ s²(A/(2z²) + 3/2) (u − log(1+u) ≤ u²/2), so
    Q(z) ≥ K(z)/s with K := H − s³|L'| = H − s³(A/z² + 3).
(c) [P + cert] K > 0 at ẑ.  By (P1), N ≥ 0 exactly on S ∈ [0, Ŝ], so it suffices that K(S) ≤ 0 ⇒ N(S) < 0.
    K is affine decreasing in A, so K ≤ 0 ⇔ A ≥ A_K(S) := z²(H − 3S³)/S³.  N is affine in A with A-coefficient
    2S²[2x0(z0 − t) + S(t + 1 − 2z0)] < 0 on (0, x0) (affine in S, negative at both ends), and affine in B with
    B-coefficient BC := S²[4z0x0(z0 − t) + 2μ0'S + 2S³] < 0 (certified).  Hence for K ≤ 0, A ≥ 0, B ≥ Bf:
       N(S; A, B) ≤ N(S; max(A_K, 0), Bf) < 0     (certified, both cases A_K ≥ 0 and A_K ≤ 0).
    NO R3 constraint is used.  Certificate: step121_build.py (4 polynomials in (T, ζ, σ = S/x0), degrees ≤ 7),
    cert2/p2cert.py, cert2/run_P2_all.sh: 19 charts covering [0,1]³ — two T-slabs, the corner (T, ζ) -> (0, 0)
    (ζ = Tw / T = ζw), and blow-ups of the degenerate edges (ζ, σ) -> (1, 0) and (T, σ) -> (0, 1); all 19 certified, ~1 min.
    [N] margins: Q(ẑ) normalised ≥ 0.72 on 55k R3 points; the certified inequality itself is tight only on those edges.
(d) CONCLUSION.  For every R3 parameter point the piecewise f (phase 1 F_C = μ/s² on (1, ẑ], phase 2 the solution of
    D̃f = 0 on [ẑ, z0]) is C², positive, has f(1) = 0, satisfies (C) and (D) on (1, z0) and the seam condition
    f'(z0) > 0.  This replaces the rational multiplier of items 104–116 and every certificate built for it.
    THEOREM A's event lemma: R1, R2 by item 47; R3 by items 117–121.

## Item 123 — Attack on REFEREE2 §C (Im σ > 0 in Theorem A): reduction of the high-frequency input
Goal: N(b) := #(s = +2 Teukolsky QNMs with Im σ > 0, with multiplicity) is 0 on every subextremal KdS b.
Path: roots scaling b(s), s ∈ (0,1], to SdS (REFEREE2 §B).

Step 1 (counting, given compactness). Petersen–Vasy give an analytic Fredholm family P_b(σ) on Im σ > −c, continuous in b.
If (K) the QNMs in Im σ ≥ 0 lie in a compact set K_U for b in a neighbourhood U, uniformly, then by Gohberg–Sigal N(b)
is locally constant unless a QNM crosses Im σ = 0.  Crossing at σ ≠ 0 real: excluded by Theorem A (real σ, PROVED).
Crossing at σ = 0: Theorem A at σ = 0 gives P_b(0) invertible; invertibility is open in (b, σ), so a disc around 0 is
QNM-free for b near b0.  Hence N is constant along the path. PROVED given (K) and the base case.

Step 2 (compactness (K) from HPV, the new idea).  Let Ψ be an s = +2 QNM, Im σ ≥ 0, σ ≠ 0.
  (a) Teukolsky–Starobinsky: the 4th-order TS map sends it to an s = −2 QNM at the same σ, nonzero unless the TS
      constant D_b(σ, m, λ) vanishes (algebraically special frequency).
  (b) Wald's adjoint method (valid on any type-D Λ-vacuum: S E = O T holds with Λ) builds h = S†Φ from the −2 mode Φ:
      a linearised Einstein mode at frequency σ with Ψ0[h] = (TS map)(Φ) ≠ 0 off the algebraically special set.
  (c) Gauge-fix h: h' = h + δ*ω, ω solving HPV's gauge-propagation equation at σ (solvable off that operator's QNMs,
      which lie in a compact set: principally scalar wave operator on 1-forms, standard estimates, used in HPV Thm 4.1).
      h' is a QNM of HPV's L_b with nonzero (gauge-invariant) Ψ0.
  (d) HPV Thm 2.2/3.1: QNMs of L_b in Im σ ≥ 0 lie in a compact set; their trapping hypothesis is an open condition,
      so the bound is locally uniform in b.
  So every UHP Teukolsky QNM with large |σ| is an L_b QNM, contradiction for |σ| large: (K) holds,
  CONDITIONAL on the sub-lemmas below.

Step 3 (base case).  SdS: a UHP Teukolsky QNM gives, by 2(a)–(b), a linearised Einstein mode with Ψ0 ≠ 0, contradicting
HV18 (every such mode is pure gauge + linearised KdS, all of which have Ψ0 = 0).  Same sub-lemmas.

Remaining sub-lemmas (the gap is now these, not "uncontrolled high-frequency analysis"):
  (C2) Regularity: S†Φ is smooth across both future horizons when Φ is a future-regular mode (check the
       Hartle–Hawking/regular-tetrad boost weights of S†; in Kerr this is the DHR computation; KdS has two horizons
       of the same type, so the check is the same at each).
  (C3) Algebraically special frequencies: D_b(σ, m, λ) ≠ 0 for Im σ ≥ 0, |σ| large (locally uniformly in b), and the
       finitely many remaining special frequencies treated directly (explicit algebraically special solutions;
       show they are not regular at both horizons).
  (C4) Local uniformity of HPV's Thm 2.2 bound — state the openness argument precisely (their proof uses
       a threshold/subprincipal condition at the trapped set and at the horizons, all strict inequalities).

## Item 124 — Closing REFEREE2 §C sub-lemmas (C2, C3) and the α bound
* α bound PROVED: α = Λa²/3 < 7 − 4√3. Proof by scaling plus AM–GM/QM–AM on the three positive roots: |e4|/p2² ≤ 1/48, with
  equality only when the roots are equal (extremal). Written in the paper as Lemma (rotation bound).
* C3 BYPASSED: no Teukolsky–Starobinsky map is needed. The s = +2 mode Ψ is used directly as an ORG Hertz potential, and
  h̃ = conj(S'†Ψ) has frequency −σ̄ (still in the UHP) with Ψ4[h̃] = ½ Þ'^4 Ψ̄.
* Kernel lemma PROVED: for Im σ ≠ 0, a mode in ker Þ'^4 is a polynomial (degree ≤ 3, along n) times an n-invariant mode,
  which goes as e^{−iσv}. It is therefore regular at r1 but on the irregular branch at r0. The two Frobenius exponents at r_j
  differ by s + 2iK_j/Δ'_j, and Im K_j ≠ 0. So Ψ4[h̃] ≠ 0. This uses both horizons, so it is specific to Λ > 0.
* C2 PROVED: S'† is GHP-covariant, so evaluate it in a tetrad that is smooth across each future horizon; a regular mode is
  smooth there, so h̃ is smooth.
* Gauge-fixing gap closed: a resonance of Υδ* with large |σ| gives a nonzero pure-gauge resonant state δ*ω of L_b
  (L_b δ*ω = δ̃*Υδ*ω = 0), so it is excluded by HPV too.
* Literature inputs: Wald's adjoint/ORG identity with Λ (used in Kerr–AdS by Dias–Santos; the identity is algebraic in Λ);
  Petersen–Vasy Fredholm theory; HPV Thm 2.2/3.1.
* C4: an openness argument is written; it still needs to be checked against HPV §3 (their proof is pointwise in b).
* Theorem A in the UHP: PROVED modulo these literature inputs and C4. Paper §"The upper half plane".

## Item 133 — Full-window L(z₀) > 0: exact certificate (nonnegative-Bernstein)
* New idea: no domination or charts are needed. On the open cube (0,1)^3, p > 0 if all Bernstein coefficients are ≥ 0 and some
  coefficient is > 0, since every Bernstein basis function is strictly positive on the open cube. Zero coefficients are allowed,
  so this handles polynomials that vanish on cube edges. `cert3/nncert.py` also has the subdivision version with the face
  condition, but it was not needed.
* `cert3/fullcert.py`: factor F₁ and P (b eliminated) for each of λQ₂, λH₂, λH₁, map each factor to (x, v = 1/r₂, τ), and
  certify the sign of every factor, including the even-multiplicity ones, which are therefore nonzero.
* Result (`cert3/fullcert_out.txt`): F₁ > 0 and P > 0 for Q₂, H₂ and H₁ on the whole window. Every factor closes with a single
  box, with no subdivision.
* Consequence: L(z₀) > 0 on the full window for all m ≥ 1, and every subextremal KdS, via the concavity-in-1/m reduction
  (item 127).
* Still to re-audit when writing up: that the window is exactly the open cube in (r₀, 1/r₂, τ); that the denominators removed by
  the gcd-reduced class F have known sign; and the reduction "Q > 0 ⇔ F₁ > 0 ∧ P > 0" (item 131).
* The domination machinery (domcert.py, autocert.py, remcert.py, edge*.py) is superseded but kept for the record.

## Item 134 — Audit of item 133 closed; event lemma (P1) on the new range
* Audit (item 133), all three points closed:
  - Denominators: win_red3 reduces e.n·e.d, which has the same sign as e. So nothing with unknown sign was removed.
  - Root choice: e4 = r0 r1 r2 r3 < 0, so b² − e2 b + e4 has exactly one positive root, b+ = (e2 + √D)/2, and b = a² is that
    root. G(b+) = (−Y + F1√D)/2 > 0 ⇐ F1 > 0 ∧ F1²D − Y² > 0. Correct in both directions of Y's sign.
  - Domain: r0 ∈ (0,1), r2 > 1, τ ∈ (0,1) is the open cube in (r0, 1/r2, τ).
  So L(z0) > 0 on the window is PROVED.
* (P1) without the false floor. The old Descartes proof (item 119) needs B ≥ 5.5t − 1.25 and fails with B ≥ 4t − 1/2.
  New facts:
  - [P] N at z = 1 equals −(t−1)² x0² (p − 1), with p = 2x0 L(1)/(t−1). So N(z=1) < 0 exactly on p > 1, i.e. R3. At the
    corner (A = Ac2, L(z0) = 0) the pp line, the R2/R3 gate and L(z0) = 0 all meet, at p = 1.
  - [P] With B = 3z0 − A/z0 + β (β = L(z0) > 0), the Möbius coefficients T2..T6 are affine in (A, β) with
    coefficients ≤ 0 (step134_P1new.py).
  - [N] On R3 ∩ {L(z0) > 0}, N has exactly ONE root in (0, x0): 20000/20000 random points, plus all 91 near-corner points
    checked with mpmath. P1 holds without any B floor.
  - [N] But the single-transform Descartes pattern has more than one sign change on about 0.05% of points
    (T_j < 0 < T_{j+1} for some j ≥ 2), so item 119's proof does not carry over.
* NEXT: prove P1 on the new range by Descartes with one bisection (separate Möbius maps on S ∈ (0, x0/2) and (x0/2, x0)),
  or by "N = 0 ⇒ N' < 0". Then redo P2 (item 121 used B ≥ Bf only in the last certificate: replace it with
  B = 3z0 − A/z0 + β), and R1/R2 (item 47).

## Item 135 — (P1) PROVED on R3 ∩ {L(z0) > 0}, with no B floor
Write B = 3z0 − A/z0 + β, where β = L(z0) > 0 (item 133). The R3 constraints p > 1 and A ≥ Ac2 reduce to the quadrant
{A ≥ a0 := max(0, Ac2), β ≥ 0}, with the corner on p = 1, outside R3. Möbius S = x0 v/(1+v): coefficients T0..T6
(step134_P1new.py, T134.pkl).
* T0 = 6z0³(z0−t)²x0² > 0, and T1 = 6z0²(z0−t)x0²(3z0² + 2z0 − 4tz0 − t) > 0. Neither depends on (A, β).
* T6 = z0 · N(z = 1) = −z0 (t−1)² x0² (p − 1) < 0 on R3. [P, exact factorisation]
* For k = 3, 4, 5, T_k is affine with STRICTLY negative A- and β-coefficients. So {T_k ≥ 0} ∩ quadrant is a closed
  triangle at the corner, or empty. T_{k−1} is affine with nonpositive coefficients, so T_{k−1} > 0 on that triangle iff
  g1 := a|qA| − |pA| b > 0 and g2 := a|qB| − |pB| b > 0 whenever b ≥ 0. Here a and b are T_{k−1} and T_k at the corner,
  and p•, q• are their coefficients.
  CERTIFIED: g1 > 0 and g2 > 0 UNCONDITIONALLY, for k − 1 = 2, 3, 4 in both cases (I: a0 = Ac2, II: a0 = 0). The proof
  factors g and certifies the sign of every factor by nonneg-Bernstein (nnd.py) on (u, R) ∈ (0,1)², with
  case I: t − 1 = 6x0²/u, x0 = u/6 + R/(1−R); case II: t − 1 = x0 + (6x0² − x0)u, x0 = 1/6 + R/(1−R). Script:
  step135g.py, all 12 signs +1. (Floating-point samples had shown about 0.4% "violations"; exact rational re-evaluation
  showed all of them were roundoff.)
* Hence T_k ≤ 0 ⇒ T_{k+1} < 0 for k = 2..5. The sequence (T0, ..., T6) has exactly one sign change, so by Descartes N has
  exactly one root in (0, x0). ∎ (P1)

## Item 136 — (P2) on the new range: item 121's sufficient condition FAILS; P2 itself holds numerically
* Item 121(c) used B ≥ Bf. With B = 3z0 − A/z0 + β, the A-coefficient of N is 2S³(−S² − x0(t − z0)) < 0, but
  "K(Ŝ) > 0" fails on about 67% of (t, z0, S) at the corner, and Q(ẑ) > 0 fails on 168 of 3000 random feasible points
  (step136b.py).
* [N] The seam condition f2'(z0) > 0 itself holds on all 80 points where Q(ẑ) ≤ 0 (ODE integration, step136c.py).
  The bootstrap of item 121(a) is too lossy near L(z0) → 0: it replaces f2 by f̂.
* NEXT: a sharper bootstrap. Keep the growth of f2 (f2 ≥ f̂ e^{∫k} with k a subsolution near ẑ), or compare Φ with
  the explicit phase-1 data on a short interval past ẑ.

## Item 137 — Quantitative margin L(z0) ≥ A/(2z0) PROVED on the window; R1/R2 survive; P2 status
* PROVED (cert3/margin_red.py + fullcert_m.py, same exact pipeline as item 133): on the whole window, for all m ≥ 1,
  L(z0) − A/(2z0) > 0. The proof replaces the A-term coefficient 4 in LQ2 by 2. All six sign certificates (F1 and P for
  Q2, H2, H1) pass with no subdivision. Concavity in 1/m is unaffected, since only the m² coefficient changes.
  At coefficient 1 (margin 3A/(4z0)) the H2-P certificate does not close.
* R1/R2 (item 47) need only L̃ > 0 on [1, z0] (L decreasing ⇒ ⇔ L(z0) > 0 in R1; in R2, L̃(z0) > 0 IS the R2 condition)
  and H ≥ s(t−1). Neither uses the B floor. So the event lemma in R1 ∪ R2 is PROVED unconditionally.
* P2 in R3 with β = L(z0) ≥ A/(2z0):
  - [N] Item 121's Q(ẑ) > 0 still fails on 20 of 3000 random points (step137c.py 0.5). The only lossy step is f2 ≥ f̂
    in the bootstrap.
  - [N] The one-step improvement (freeze Φ ≥ Φ̂, solve the linear lower ODE f_l' = (Φ̂ + 2sL f_l)/μ, then use
    Φ(z0) ≥ Φ̂ + 2∫L f_l) passes 3000/3000 with normalised margin ≥ 0.71. Each Picard iterate is a rigorous lower bound
    (L > 0, s > 0), but it is an integral, not a closed form.
  - NEXT: an explicit lower bound for f2 on [ẑ, z1]. Integrand (Ĥ/ŝ − 2ŝL̂ + 2sL)/μ, with sL decreasing; certify the
    resulting closed-form inequality in (t, z0, S, A) on the margin set.

## Item 138 — Margin L(z0) ≥ A/(2z0) + z0 PROVED; P2 reduced to an explicit comparison criterion
* PROVED (cert3/margin_red.py, LQ2 coefficient 2 and LQ0 term λ̃(1+t−4z0); fullcert_m.py: all six signs +1).
  The expression is (Λ/3)·[L(z0) − A/(2z0) − z0] in the m-quadratic form; λ̃ = −b/e4 = Λ/3 > 0.
  So on the whole window, for all m ≥ 1:  β := L(z0) ≥ A/(2z0) + z0.
* P2 via the Riccati variable w = μ f2'/f2. It satisfies μw' = −w² + 2sLw − 2μs|L'|, with w(ẑ) = ŵ := Ĥ/ŝ > 0, and
  f2'(z0) > 0 ⇔ w(z0) > 0. Because w crosses 0 only downward, it suffices that w > 0 on [ẑ, z0].
  Rigorous comparison chain:
  (i)   Upper bound. μw' ≤ −w² + 2sLw, so w ≤ U_B (Bernoulli: 1/U_B solves μy' = 1 − 2sLy, y(ẑ) = 1/ŵ). Also
        w ≤ max(ŵ, 2ŝL̂), because sL is decreasing and w decreases wherever w > 2sL. Set M := min(U_B, max(ŵ, 2ŝL̂)).
  (ii)  Chord. For 0 ≤ w ≤ M, −w² ≥ −Mw, so while w ≥ 0: μw' ≥ (2sL − M)w − 2μs|L'|.
  (iii) Linear comparison. w ≥ v, where v' = ((2sL − M)/μ)v − 2s|L'| and v(ẑ) = ŵ. This has the explicit solution
        v = e^{I}(ŵ − ∫ 2s|L'| e^{−I}), with I = ∫_ẑ (2sL − M)/μ.
  CRITERION (P2*): min over [ẑ, z0] of v > 0.
* [N] (P2*) holds on 3000/3000 random R3 points with β ≥ A/(2z0) + z0 (step138i.py 0.5 1). Without the +z0 term it
  fails on 4 of 3000, and with only the constant M it fails on 82.
* REMAINING for P2: certify (P2*) over the 4-parameter family (t, z0, A, β). It is a nested-integral inequality, so it
  needs validated quadrature with exact bounds on the integrands (rational functions plus the logs from ∫1/μ).

## Item 139 — Margin L(z0) ≥ A/(2z0) + 2z0 PROVED; P2 by an explicit polynomial barrier, except one corner
* PROVED (cert3/margin_red.py with LQ0 term λ̃(1+t−5z0); fullcert_m.py: all six signs +1): L(z0) ≥ A/(2z0) + 2z0 on the
  window. Margins 3z0, 4z0 and 6z0 do not certify (H2-P and H1-P return None).
* P2, barrier version (replaces the nested-integral criterion of item 138 wherever it works). With
  w = μf2'/f2 and σ = s/ŝ, take k := ŵ σ² (1 − σ + σ²). If
     (B)  μk' + k² − 2sLk + 2μs|L'| ≤ 0   on (ẑ, z0),
  then k is a subsolution of the Riccati equation. With k(ẑ) = ŵ, this gives w ≥ k > 0 on [ẑ, z0), and strict comparison
  gives w(z0) > 0. Condition (B), divided by s, is POLYNOMIAL in (t, z0, A, β, Ŝ, σ).
  - Eliminations (step140_build.py, step140e.py): N is affine in A and in β, so the root condition N(Ŝ) = 0 can be
    solved for A at fixed (t, z0, Ŝ, γ), with β = A/(2z0) + c z0 + γ, γ ≥ 0. No R3 gate is needed numerically.
  - [N] (B) holds on 9000/9000 generic random R3 points with margin ≥ A/(2z0) + z0 (one failure at seed 23, in the corner
    below). The plain k = ŵσ² fails on 2 of about 16000, and k = ŵσ^p with p > 2.1 fails badly.
  - FAILS in the CORNER t ≳ 100, zf = x0/(t−1) ≲ 0.02, ẑ → 1 (Ŝ/x0 → 1): about 7% failures for c = 1 and about 4% for
    c = 1.5 (step141b.py). This is the known degenerate corner (cf. item 33's rescaling x = (z0 − 1)y). It needs a
    barrier adapted to μ ≈ t·x·z near z = 1.
* NEXT: (a) corner barrier in rescaled variables; (b) exact certificate of (B) on the complement, which is a 4- or
  5-variable polynomial with σ and the γ ≥ 0 gate.

## Item 140 — The P2 corner (t ≫ 1, zf small, ẑ → 1): every comparison scheme tried fails on a few percent
Test set: 3000 R3 points with t ∈ [1e2, 3e3], zf ∈ [1e-3, 0.05], margin β ≥ A/(2z0) + 2z0 (step141b.py).
* Polynomial barriers ŵσ^p(1 − cσ(1−σ)) and ŵσ²e^{−aσ(1−σ)} (step142.py): best 2914/3000, union 2981/3000.
* The equilibrium-type barrier θ·2sL (step142b.py): at most 2556/3000. Max-combinations with the polynomial barrier,
  and with the Riccati root k₊ = sL + √(s²L² − 2μs|L'|) (step142c/e/f): at most 2936/3000.
* The nested-integral criterion P2* of item 138 (step142g.py): 2944/3000.
* The TRUE solution is comfortable: w(z0)/ŵ ≈ 0.27–0.33. Early on it tracks 2sL closely (w ≈ 2sL + O(ŵ/20)), then
  flattens (step142d.py). The schemes fail because w hugs the moving root 2sL − μ|L'|/L early on: any barrier below it
  must decrease at least as fast as (2sL)', and polynomial-in-σ shapes overshoot at σ ≈ 0.9.
* Structure: under ω = w/t, Â = A/t, B̂ = B/t, μ ≈ t·x·z, the corner has a t-independent LIMIT PROBLEM
  x z ω' = −ω² + 2(sL/t)ω − 2xz s|L'|/t on [ẑ, z0], with parameters (x0, Â, β̂).
* NEXT: prove P2 in the corner by (i) the t → ∞ limit problem plus an explicit O(1/t) perturbation bound, or (ii) a
  two-zone barrier that follows the root k₊ − Cμ|L'|/L until s²L² = 4μs|L'|, then switches to a polynomial tail
  matched at the switch point.

## Item 141 — Two-zone barriers for P2: zone 1 is free, the tail is the obstruction
* [P] ZONE 1 LEMMA. Wherever D := s²L² − 2μs|L'| > 0, the Riccati right side R(k) = −k² + 2sLk − 2μs|L'| is positive
  strictly between its roots k± = sL ± √D. So ANY nonincreasing k in (k−, k+) is a subsolution: μk' ≤ 0 < R(k).
  Examples: sL; sL + θ√D for θ ∈ (0,1); min(ŵ, ·) wherever ŵ < k+. Min and max of subsolutions are subsolutions in the
  Dini sense, which is enough for ODE comparison.
* So all the difficulty sits in the tail zone near z0, where D ≤ 0 and R < 0 everywhere.
  Tested tails (step143*.py): linear in s (k = cs), linear in z with positive end value (k = k* − ρ(z − z*), ρ by
  fixed point), with switch points scanned. Best: corner 1995/3000, generic 1490/1500. The tail needs
  |k'| ≥ (k² − 2sLk)/μ + 2s|L'|, and 2s|L'| is large at the switch point in the corner (A large).
* Conclusion: the corner needs the actual tail dynamics. Next route: the t → ∞ limit problem (item 140) on the tail
  zone, where D ≤ 0 and the equation is μw' = −(w − sL)² − |D|: a monotone decreasing Riccati equation with an
  explicit lower solution once |D| is bounded above by a quadratic in s.

## Item 142 — Explicit tail bound: rigorous but too lossy in the corner
* [P] TAIL LEMMA. Fix a switch point z_s. Suppose zone 1 (item 141) gives w(z_s) ≥ v_s := min(ŵ, s_sL_s + θ√D_s). Then
  v(z) := v_s − ∫_{z_s}^{z} (v_s² + max(−D,0))/μ is a lower solution for w. Reason: for 0 ≤ v ≤ v_s and 0 ≤ sL ≤ v_s
  we have (v − sL)² ≤ v_s², and the true equation is μw' = −(w − sL)² + D. So
     w(z0) ≥ v_s − v_s² ∫_{z_s}^{z0} dz/μ − ∫_{z_s}^{z0} max(−D,0)/μ dz,
  which is explicit (∫dz/μ by partial fractions; the D-integral is rational plus logs).
* [N] With the switch point scanned (step144*.py): corner 1749/3000, generic 1485/1500. The loss is replacing
  (w − sL)² by v_s²: in the tail w stays near 0.3ŵ while sL → 0, so the true drop is much smaller than the bound.
* Still to try: a tail bound that keeps the Riccati structure. (1/v)' ≤ (1 + max(−D,0)/v_min²)/μ gives an explicit
  lower bound on v once a floor v_min is assumed, which can be closed by bootstrap.

## Item 143 — P2 via the Φ-integral: no tail barrier needed. Generic region 100%, corner about 97%
* [P] Φ-CRITERION. Normalise f2(ẑ) = 1. Then μ0 f2'(z0) = Φ(z0) = Φ̂ + 2∫_ẑ^{z0} L f2, with Φ̂ = Ĥ/ŝ − 2ŝL̂. Here
  L > 0 and f2 > 0 (f2 = exp∫w/μ is finite and positive whatever the sign of w). So ANY lower bound f2 ≥ f_lo ≥ 0 on any
  subinterval gives f2'(z0) > 0 as soon as Φ̂ + 2∫ L f_lo > 0. A tail where w turns negative only costs its (nonnegative)
  contribution, which is dropped.
* Lower bounds used:
  - zone 1 (item 141): f2 ≥ exp∫k/μ with k = min(ŵ, sL + θ√D), on the initial run where D > 0 and k is nonincreasing;
  - tail: w ≥ v, where v' = −(v_s·min(v_s, v + sL) + max(−D,0))/μ, using (v − sL)² ≤ v_s·min(v_s, v + sL) for
    0 ≤ v, sL ≤ v_s. Only the part with v > 0 is used.
* [N] (step145b/c/d.py): generic R3 1500/1500 for every θ tried; corner (t ∈ [1e2, 3e3], zf ≤ 0.05) 2902/3000 with
  θ = 0.9. Zone 1 alone: 2527/3000.
* Status: generic P2 now rests on an explicit nested-integral criterion, which still has to be certified (validated
  quadrature). The corner still loses about 3% and needs the limit-problem analysis or a sharper tail.

## Item 144 — P2: the θ·sL barrier family covers EVERY test point
* [P] θ-BARRIER LEMMA. Let k = θ(z)·sL with θ piecewise C¹ (Dini sense at kinks) and θ(ẑ)·ŝL̂ ≤ ŵ. If
     g := θ(2−θ)s²L² + θμL − (2−θ)μs|L'| − μθ' sL > 0   on (ẑ, z0),
  then k is a strict subsolution of μw' = −w² + 2sLw − 2μs|L'|. Hence w > k ≥ 0 on [ẑ, z0). Near z0, g → θμ0L0 > 0
  forces w(z0) > 0. So f2'(z0) > 0, which is P2.
  (Derivation: (sL)' = −L − s|L'|; R(θsL) = θ(2−θ)s²L² − 2μs|L'|.)
  For constant θ, g is a POLYNOMIAL in (t, z0, A, β, Ŝ, u) after clearing positive denominators. The start condition
  is Ĥ ≥ θŝ²L̂.
* [N] (step146b/c/f.py):
  - constant θ ∈ [0.2, 1.95] (one θ per point): generic 2000/2000, corner 2952/3000;
  - piecewise-linear θ (rise θ0 = ŵ/(ŝL̂) → θ1 on [0, u1], then linear down to θ2): ALL 55 remaining corner points pass
    (55/55).
  So, with the proved margin β ≥ A/(2z0) + 2z0, P2 holds on all 5000 test points via an explicit polynomial barrier.
  The hard points have ŵ < ŝL̂: there θ must start below 1 and rise to about 1.1 by u ≈ 0.3.
* REMAINING to make P2 a theorem: partition the parameter space into cells, fix θ per cell (θ constant on most cells),
  and certify g > 0 by Bernstein, as with items 133–137. The variables (t, z0, A, β, Ŝ, u) are linked by N(Ŝ) = 0:
  eliminate A or β via the affine dependence (item 139).

## Item 145 — P2 reduced to root-free polynomial implications; true on the physical set, false on the abstract 4-parameter set
* Root-free form. By P1, N ≥ 0 exactly on z ∈ [ẑ, z0]. So the θ-barrier (item 144, constant θ) proves P2 provided
     (a) N(z) ≥ 0 ⇒ g_θ(z) > 0   and   (b) N(z) ≥ 0 ⇒ H(z) − θs²L(z) ≥ 0,   for all z ∈ (1, z0)
  ((b) at z = ẑ is the start condition). Both are POLYNOMIAL implications in (t, z0, A, B, z). No root ẑ is needed.
* [N] On sampled R3 points (step147d.py): 29952/30000 with a per-point θ from a grid in [0.3, 1.95].
* g_θ is AFFINE in t (μ is linear in t), so on each t-interval only the endpoints matter.
* [N, NEGATIVE] On the abstract R3 ∩ margin set, (a) FAILS for every θ at extreme points: x0 ≳ 10, A ≳ 10³, t near the
  R3 upper limit t − 1 = min(2x0L(1), 2x0²(A/z0 + 3)), z near 1. Rechecked with exact N via mpmath (step147g.py).
  So the abstract 4-parameter event lemma cannot be closed this way.
* [N, POSITIVE] On PHYSICAL points (window map (r0, r2, frac, m), m = 1..29, Galerkin λ; nearext/phys147.py):
  1851/1851 R3 points pass (a) and (b) with constant θ.
* Consequence: certify P2 on the physical parametrisation (cube variables x, v, τ, u = 1/m, plus the position z),
  reusing the machinery of items 133–139 (b eliminated through b² − e2b + e4 = 0). Large m is the u → 0 face.

## Item 146 — What separates the abstract failures from physics; monotonicity in B
* The abstract (a)-failures of item 145 all have B ≪ t: (B − t)/t ≈ −1, and t/x0 between 6·10³ and 2·10⁸. R3
  (t − 1 < 2x0L(1)) then forces A huge, with B of order −A/(2z0).
* [N] PHYSICAL points with x0 ≥ 5 have (B − t)/t ≥ 0.57 (1279 samples). Restricting the abstract set to B ≥ t
  removes every failure (step148.py). But B ≥ t is FALSE physically in general: for near-extremal, m = 2, top of the
  window, (B − t)/t reaches −1.26 with λ ≈ −4.2 (nearext/btest.py). So the right physical input is a joint (x0, t, B)
  constraint, not a single inequality.
* [P] MONOTONICITY. g_θ is increasing in B wherever L > 0. N is affine in B with negative coefficient (item 121), so
  raising B shrinks {N ≥ 0}. Hence implication (a) is MONOTONE in B: it suffices to check it at the crude angular
  lower bound B_crude = c1 + λ_crude·C0. Implication (b), the start condition, is not monotone and needs separate
  handling.
* Plan for the physical certificate of P2: cube variables (x, v, τ), u = 1/m, and z, with A, t, z0 exact and
  B := B_crude. The crude λ bound is polynomial in (m, α, ξ). b is eliminated as in item 131. (a) is an implication
  with gate N ≥ 0; handle it by λ-multiplier certificates (−g − c·N ≥ 0 type, item 119 style) or by explicit
  parametrisation of {N ≥ 0}.

## Item 147 — P2 (a) reduced BY HAND to one physical inequality: ẑ ≥ (3/5) z0
* [N] On 117819 physical R3 points with B = B_crude (and on 119395 more with B anywhere above B_crude), every constant
  θ ∈ [0.3, 1.1] satisfies both (a) and (b) (nearext/phys149*.py). So ONE θ works universally; take θ = 1.
* [P] HAND LEMMA for θ = 1. Write L = ℓ + sP with ℓ = L(z0) and P = A/(z z0) + 3. Then |L'| = P + As/(z²z0), so
  L − s|L'| = ℓ − As²/(z²z0). With the certified margin ℓ ≥ A/(2z0) + 2z0 (item 139):
     L − s|L'| ≥ A(z² − 2s²)/(2z²z0) + 2z0 > 0   whenever z ≥ √2 s, i.e. z ≥ (2 − √2) z0.
  So g_1 = s²L² + μ(L − s|L'|) > 0 on [(2−√2)z0, z0). Hence (a) holds as soon as ẑ ≥ (2−√2) z0 ≈ 0.586 z0.
* [N] Physically ẑ/z0 ≥ 0.620 (93016 R3 points with B_crude; nearext/phys149d.py). We use the rational threshold
  3/5 > 2 − √2. By P1, ẑ ≥ (3/5)z0 ⇔ N(z = 3z0/5) < 0 (needed only when z0 > 5/3). Exactly:
     N(3z0/5) = N0 + A·NA + B·NB, with
     NA = (16z0²/125)(t(5 − 4z0) + 3z0² − 4z0) < 0 and NB = (16z0³/3125)(t(100 − 75z0) + 54z0² − 75z0) < 0 for z0 > 5/3,
     N0 = (2z0²/3125)(3375t²z0² − 9000t²z0 + 5875t² − 3672tz0³ + 11700tz0² − 9000tz0 + 891z0⁴ − 3672z0³ + 3375z0²).
  N is decreasing in A and in B, so B_crude is the worst case.
* Certificate in progress (cert3/zhat_build.py, zhat_red.py, zhat_cert.py): λ̃N in the window variables, quadratic in m.
  The m-coefficients are q2, q1, q0, with q0 < 0 numerically (concave in u = 1/m). Targets are q2 > 0, m = 2 (Hhalf)
  and m = 1 (H1). First obstacle: after b-reduction, q2's F1 core (degree (62, 62, 12)) is NOT sign-definite on the
  whole cube (40+/20− at random points). So the sign test must use the full rule (F1 < 0 ∧ Y < 0 ∧ P < 0 also
  certifies), and possibly the gate z0 > 5/3.
* Still needed for P2: the start condition (b) at θ = 1, H(ẑ) ≥ ŝ²L̂. It holds on all physical points; no proof yet.

## Item 148 — Gate-free threshold z* = 1 + (31/50)x0²/z0; fast Bernstein sign test; leading m-coefficient PROVED
* Threshold. The gate z0 > 5/3 of item 147 is removed by z* := 1 + (31/50) x0²/z0 ∈ (1, z0).
  [P] z* ≥ (2−√2)z0 for every x0 > 0. Writing z0 = 1 + x0:
      z0 + 0.62x0² − (2−√2)z0² = 0.4142 − 0.1716x0 + 0.0342x0² (coefficients rounded),
  with negative discriminant, so it is positive. Hence the hand lemma (item 147) makes g_1 > 0 on [z*, z0).
  [N] Physically ẑ ≥ z*: 0 violations in 99461 R3 points (nearext/phys149d.py). The simpler 1 + 0.6x0 fails on 159.
* FAST BERNSTEIN TEST (cert3/fastbern.py). The Bernstein coefficients on [0,1]³ have the signs of the coefficients
  of rev(shift(rev(p))): rev(p) = xⁿp(1/x) per variable, and shift is x ↦ x + 1 (one flint compose). This replaces
  the Fraction arrays; the (162, 164, 38) polynomial is tested in 35 s. subdiv_pos does dyadic subdivision exactly
  (p(x/2)·2ⁿ, p((x+1)/2)·2ⁿ).
* Build (cert3/zhat2_*.py). λ̃·(−N(z*)) = m²q2 + mq1 + q0 at B = B_crude (m ≥ 2), plus the m = 1 variant H1.
  - [P] q2 > 0 on the whole window cube: F1 > 0 by factor signs (core (93, 93, 19) passes one Bernstein box);
    P = g²·Pr with Pr > 0 by one Bernstein box (degree (162, 164, 38)) and g nonvanishing (all factors definite).
    So N(z*) < 0 for all m ≥ m0(geometry): P2(a) holds for large m everywhere.
  - Hhalf (m = 2): its F1 core (110, 110, 24) is NEGATIVE on the face τ = 1 (window top) and positive at τ ≤ 0.99.
    [N] Physically N(z*) ≥ 0 does occur near the window top, but ONLY at points outside R3 (66/180 test points, all
    with R3 false). So the m = 1, 2 certificates need the R3 gate (or a split near τ → 1, where R1/R2 take over).
* Still open for P2: the R3-gated certificates for m = 1, 2 (and q0 ≤ 0 for concavity); start condition (b).

## Item 149 — P2 condition (a) PROVED for all m, by hand, from a new certified margin L(z0) ≥ A/(2z0) + 2t
* [P, certificate] L(z0) ≥ A/(2z0) + 2t on the whole window, for all m ≥ 1 (cert3/margin_red.py with LQ2 coefficient 2
  and LQ0 = −K(1+7α) + λ̃(1 − t − 3z0); fullcert_m.py: all six signs +1). Physically min (ℓ − A/(2z0))/t ≈ 2.63
  (nearext/ell153.py), so the certified constant is close to sharp.
* [N] Physically g_1 > 0 on ALL of (1, z0), not just on [ẑ, z0): 206179 points, normalised margin ≥ 0.745
  (nearext/g152.py). Abstractly, ℓ ≥ A/(2z0) + 2t ⇒ g_1 > 0: 0 failures in 3·10⁶ samples (step153.py).
* [P] HAND PROOF of (a), θ = 1, on the whole interval. Let ℓ = L(z0), a := A/(2z²z0), s = z0 − z, z ∈ (1, z0).
  Exactly: L = ℓ + s(A/(zz0) + 3) and L − s|L'| = ℓ − As²/(z²z0).
  With ℓ ≥ A/(2z0) + 2t:  L − s|L'| ≥ a(z² − 2s²) + 2t  and  L ≥ az(z + 2s) + 2t =: X + 2t.
  AM–GM: (X + 2t)² ≥ 8tX. Hence
     g_1 = s²L² + μ(L − s|L'|) ≥ 8t·az(z+2s)s² − 2μas² + 2tμ ≥ 2tμ > 0,
  because 8t·z(z+2s) ≥ 2μ = 2z(z−1)(t−z) ⇔ 4t(z + 2s) ≥ (z − 1)(t − z), which is true termwise. ∎
  So the barrier k = sL is a strict subsolution on all of (1, z0), for every m and every subextremal KdS. No
  root ẑ, no R3 gate, no computer beyond the margin certificate.
* REMAINING for P2: only the START condition (b), ŵ ≥ ŝL̂, i.e. H(ẑ) ≥ ŝ²L(ẑ) at the phase-1/phase-2 switch.
  [N] Physically true (phys150c.py, all m = 1..3 points, window up to frac 0.99999). Multiplier route (step154*.py):
  T := H − s²L and N are affine in (A, γ); certify T − κN ≥ 0 on the feasible polyhedron. Infeasible only for
  x0 ≲ 0.05 with t ≈ z0 (about 2% of abstract samples): that corner needs either physics or a smaller θ there.

## Item 150 — Start condition (b): status
* [N] Physically the start ratio ŵ/(ŝL̂) = Ĥ/(ŝ²L̂) is ≥ 1.13 on 84434 R3 points (median 1.35; nearext/start155.py).
  So θ = 1 works with room to spare, including at x0 ≈ 0.006.
* [P] Facts used: at ẑ the phase-1 Riccati residual vanishes, i.e. w_C = H/s satisfies μw' = −w² + 2sLw − 2μs|L'|
  exactly there (d = 0 ⇔ D̃F_C = 0). In the multiplier form (b) ⇐ T − κN ≥ 0, with T := H − s²L.
  The choice κ = κ_min = T_B/N_B removes B entirely: Φ = T·N_B − T_B·N is B-FREE and AFFINE in A, with
  A-coefficient 2(t−z)s³(z−1)/z > 0 (step156.py). So κ_min would need an UPPER bound on A; larger κ needs a LOWER
  bound on B (the margin) and gains from N_A < 0.
* [N, NEGATIVE] Abstractly (margin ℓ ≥ A/(2z0) + 2t + γ, A ≥ a0, p > 1), the implication N ≥ 0 ⇒ H − θs²L ≥ 0 fails
  on about 1.7% of samples for EVERY θ ∈ [0.3, 1]. All of them sit in the NEAR-NARIAI corner t − 1 ≈ x0 ≪ 1
  (step155*.py, step157.py). A larger margin constant does not help. g_θ > 0 itself holds abstractly for θ ≥ 0.3,
  so only the start is at issue.
* So (b) needs one physical input that the abstract parameters do not see in the near-Nariai corner, e.g. a
  certified relation between A and (t − 1)/x0 there, or a direct physical certificate of Φ-type
  inequalities in the window variables (as in item 148, with the fast Bernstein test).

## Item 151 — (b) in root form: T* := T|_{N=0} is B-free; the condition is "B*(z) ≥ B_crude ⇒ T*(z) ≥ 0"
* N is affine in B, so N(z) = 0 ⇔ B = B*(z) (explicit rational in t, z0, A, z). At the true root ẑ, B*(ẑ) = B_true ≥ B_crude.
  So (b) ⇐ [for all z ∈ (1, z0): B*(z) ≥ B_crude ⇒ T*(z) ≥ 0], where T* = (H − s²L)|_{B=B*} is B-FREE.
* [N] A constant-κ multiplier over the whole interval (T − κN ≥ 0) fails on almost all physical points (b158.py):
  it is far too strong past ẑ, where N > 0 grows.
* [N] Root form: 37997/38090 physical R3 points pass (b158b.py). The 93 failures have not yet been separated from
  floating-point blow-up: the B* denominator vanishes inside (1, z0), which caused the divide warnings. Next: recheck
  them exactly. If they are genuine, restrict to the sub-interval where the denominator keeps its sign (P1 gives a
  single root, so only one branch is relevant).

## Item 152 — (b): the 93 "failures" are endpoint roundoff; exact structure is a pair of single crossings
* [P] As z → z0: s²B*(z) → −(3/2) z0(z0−t)(z0−1) and T* → z0(t−z0)(z0−1)/2 > 0 (symbolic limits). So T* > 0 near z0.
  Every one of the 93 flagged points (and the new ones) sits at u = 1 − 10⁻⁶ or 1 − 10⁻⁴, i.e. at z0, where the
  float evaluation cancels. Genuine failures among ~70k physical R3 points: 0 (nearext/b159.py).
* [N] Shape (nearext/b161.py and the sign strips): B*(z) is monotone increasing on (1, z0) (36191/36191 points). T*
  has a single − → + crossing at z_T, and z_T < ẑ. So
     (b) ⇔ T*(ẑ) ≥ 0 ⇔ z_T ≤ ẑ ⇔ [ T*(z) < 0 ⇒ N(z; B_true) < 0 ] ⇐ [ T*(z) < 0 ⇒ N(z; B_crude) < 0 ]
  (N is decreasing in B). Neither T* nor N(·; B_crude) involves the Galerkin eigenvalue.
* A single multiplier N − κT* < 0 with κ ∝ s² is infeasible on most points (b160.py), so the gate must be handled by
  covering: on each box of (x, v, τ, u = 1/m, w), with z = 1 + x0w, show either T* > 0 or N(·; B_crude) < 0 by the
  fast Bernstein test. The two zero sets are separated (physical start ratio ≥ 1.13), so the subdivision terminates.
  Complication: b = (e2 + √D)/2 enters both. Treat √D as an extra variable ς with ς² = D and box bounds, or reduce
  b on each box via the (F1, P) sign rule.

## Item 153 — T* is affine in A, so (b) is an UPPER bound on A at ẑ
* [P] T* = (n0 + A·nA)/den with den = 2z·(bracket of N_B) < 0 and nA = 2(t−z)s³(z−1) > 0. So T* is DECREASING in A,
  and T*(z) ≥ 0 ⇔ A ≤ A_T(z) := z·Q(t,z0,z)/(2(t−z)s³(z−1)) (Q explicit, step: see the print in this item's
  session). Hence (b) ⇔ A ≤ A_T(ẑ).
  As z → z0, A_T → ∞ (s³ in the denominator), consistent with T* > 0 near z0.
* [N] Physically: z_T := min{z : A_T(z) ≥ A} exists, and N(z_T; A, B_crude) < 0 (so z_T < ẑ) on 96887/96887 R3 points
  (nearext/bmin163.py).
* [N, NEGATIVE] What the abstract corner is missing is not a margin. With ℓ ≥ 0.6A/z0 + 2t the abstract failures stay at
  about 1.6%. Physically ℓ ≥ 0.63A/z0 + 2t holds (min excess 0.13·A/z0; nariai164.py), and near Nariai ℓ/t ≥ 11. But
  no simple bound of the form ℓ ≳ t/x0 holds globally. The separating information is a joint (A, ℓ, t, x0)
  constraint, so (b) needs a certificate in the window variables, either the covering of item 152 or z_T-parametrised.

## Item 154 — (b) by z_T-parametrisation: works for m ≥ 3; m = 1, 2 need direct certificates
* Lower bound on B without the linear m-term: since ℓ0 ≥ m, λ ≥ m² + m − 4 − 4α − 4αΞm − 2Ξ²aw·m(m+2). So
  B ≥ min(Blow1, Blow2) with Blow1 = m²(b2 + b1) + b0 (valid when b1 ≤ 0) and Blow2 = m²b2 + b0 (valid when b1 ≥ 0);
  b2 = C0(1 − 2Ξ²aw), b1 = C0(1 − 4αΞ − 4Ξ²aw), b0 = −C0(4 + 4α) + c1. N is decreasing in B, so it suffices to
  prove N < 0 for BOTH bounds, with no gate on the sign of b1.
  [N] With integer m and this weaker B: 0 violations on 101339 physical points (nearext/bmin165.py).
* z_T-parametrisation: set A = A_T(z_T) (so m² = A_T(z_T)/Â). Then everything is explicit in (geometry, z_T).
  Gates: z_T is a "record" point of A_T (first crossing), R3, and m² ≥ m0².
  [N] (nearext/zt166/167.py): with m² ≥ 9, about 0.3% flagged points remain, all at m² ≥ 3·10³ with z_T → z0, i.e.
  the float-cancellation regime of item 152. With m² ∈ [1, 9) (non-integer m included) the relaxation fails often. For
  m = 1, 2 the true λ floor (ℓ0 = 2) is much stronger than m² + m and must be used directly.
* Plan: (i) m ≥ 3: certify N(z_T; A_T(z_T), Blow_i) < 0 on the record set in window variables × z_T (exact; the
  near-z0 end is handled by the limits of item 152); (ii) m = 1, 2: fixed-m covering certificate (T* > 0 or N < 0 per
  box) in (x, v, τ, w).

## Item 155 — (b) REPLACED by a weaker start condition (b') that needs no physics beyond the certified margin
* [P] BARRIER k = min(ŵ, sL). sL is a strict subsolution on all of (1, z0) (item 149). The constant ŵ is a
  subsolution wherever R(ŵ; z) = −ŵ² + 2sLŵ − 2μs|L'| > 0. A min of subsolutions is a subsolution (Dini). And
  k(ẑ) = min(ŵ, ŝL̂) ≤ ŵ. So P2 holds provided
     (b')  R(ŵ; z) > 0 for z ∈ [ẑ, z_s], where z_s is the point at which sL drops to ŵ.
  If ŵ ≥ ŝL̂ (the old (b)), the segment is empty and (b') holds vacuously.
* [P] Sufficient form. On the segment ŵ ≤ sL ≤ k+, so R(ŵ) > 0 ⇔ ŵ > k− = 2μs|L'|/k+. Since k+ ≥ sL, this follows from
     ŵ·L(z) > 2μ(z)|L'(z)|   on [ẑ, z_s]     (which also forces D > 0 there).
* [N] ABSTRACT, with only the certified margin ℓ ≥ A/(2z0) + 2t (no R3, no physics): (b') holds on every sample. That
  includes all the old (b) failures: near-Nariai 536/536 (step168.py), global 109/109 (step168b.py), margin-only
  105/105 (step168d.py). The sufficient form ŵL > 2μ|L'| also holds on all of them (step168f.py). Without the margin
  (R3 only, or nothing) (b') fails, so the margin is the active input.
* Status: (b') is a statement about abstract (t, z0, A, ℓ) with the root ẑ, with no b and no Galerkin λ. To close it:
  parametrise by ŝ (B = B*(ŝ) from N(ŝ) = 0), impose the gate γ ≥ 0, take z on the segment, and certify
  ŵL − 2μ|L'| > 0 by covering with the fast Bernstein test. Or find the hand proof: the bound
  |L'|/L ≤ 2z0/(z(z+2s)) + 3/(2t) from the margin reduces the leading term to ẑ ≥ √2·ŝ, so the μ̂' term and the
  segment bound s ≥ ŵ/L must close the rest.

## Item 156 — (b') reduced to a purely algebraic 4/5-variable statement, numerically clean
* [N] The sufficient form ŵL(z) > 2μ(z)|L'(z)| holds on ALL of [ẑ, z0), not just on the segment, on every abstract
  sample with the margin (3349/3349; step169b.py). Replacing ŵ by min_z H/s (root-free) FAILS (30%; step170.py),
  so the value at the root matters.
* [P] Reduction to the margin-minimal B. Q := ŵL − 2μ|L'| is affine and increasing in B (∂Q/∂B = ŵ > 0). At the
  true root B = B*(ŝ), and the margin says B* ≥ Bm(A) := 3z0 − A/(2z0) + 2t. Since N is decreasing in B,
  B* ≥ Bm ⇔ N(ŝ; A, Bm) ≥ 0. Writing g := N(ŝ; A, Bm):
     Q(B*) = Q(Bm) + ŵ·g/|N_B|, so (b') ⇐ [ g ≥ 0 ⇒ Q(Bm) > 0 ].
  Statement S': for t > z0 > 1, A ≥ 0, ŝ ∈ (0, x0), z ∈ [z0 − ŝ, z0):
     N(ŝ; A, Bm(A)) ≥ 0  ⇒  (Ĥ/ŝ)·L_m(z) − 2μ(z)(A/z² + 3) > 0,   L_m = A/z + Bm − 3z.
  No physics, no b, no λ, no root.
* [N] S' holds on 192969/192969 gated samples, t up to 10⁷·x0 and A up to 10⁹ (step171.py). Without the gate it fails
  on about 3.5%, so the gate is genuinely used.
* [P] STRUCTURE for the certificate. Both g and Q(Bm) are AFFINE in A. For fixed (t, z0, ŝ, z) the gate set
  {A ≥ 0 : g ≥ 0} is a ray or a segment, so S' reduces to endpoint checks:
  (i) A = 0 when g(0) ≥ 0; (ii) A = −g0/g1 when it lies in range; (iii) the recession condition q1 ≥ 0 when the ray is
  unbounded. Each is a 4-variable rational inequality in (t, z0, ŝ, z) with no square roots. Next: certify these
  three by the fast Bernstein test after compactification.

## Item 157 — (b') PROVED. Hence P2, and the EVENT LEMMA IN EVERY REGIME
* Exact form (step172.py): with g = 2z0·N(ŝ; A, Bm(A)) = g0 + g1A and 2ŝzz0·Q(Bm) = q0 + q1A:
  g1 = 2ŝ²(−ŝ³ + ŝ(t − z0²) − 2z0(z0−1)(t−z0)). Every piece is polynomial in (t, z0, ŝ, z); no roots, no square roots.
* Compactification (step173_build.py): z0 = 1 + X/(1−X), t = z0 + x0·T/(1−T), ŝ = x0σ, z = z0 − ŝρ, with
  (X, T, σ, ρ) ∈ (0,1)⁴. Denominator signs: g0: (T−1)²(X−1)⁷ < 0; g1: (T−1)(X−1)⁵ > 0; q0: (T−1)²(X−1)⁶ > 0;
  W := q1g0 − q0g1: (T−1)³(X−1)¹¹ > 0. Degrees ≤ (10, 3, 10, 4).
* CERTIFIED by exact 4-D Bernstein covering (cover4.py, fast rev-shift-rev coefficients; per box: all coefficients
  ≥ 0, and every face not on the cube boundary carries a positive coefficient):
  - g1 < 0 everywhere (1 box). So for g0 ≥ 0 the admissible A form the segment [0, A*], A* = −g0/g1; for g0 < 0
    there is none.
  - (I) g0 ≥ 0 ⇒ q0 > 0: every box has either g0 < 0 or q0 > 0 (51 boxes, 0.4 s). This gives Q(Bm) > 0 at A = 0.
  - (II) g0 ≥ 0 ⇒ W > 0 (25 boxes, 0.6 s). This gives Q(A*) = −W/g1 > 0.
  - The face ρ = 1 (z = ẑ), which the open cube excludes: q0 > 0 and W > 0 there unconditionally (1 box each).
  Q is affine in A, so Q(Bm) > 0 on the whole admissible segment. ∎ S'
* Chain: S' ⇒ Q(B*) = Q(Bm) + ŵg/|N_B| > 0 (margin: B* = B_true ≥ Bm; N_B < 0 by item 121). Then ŵL > 2μ|L'| on
  [ẑ, z0). On the segment where ŵ ≤ sL, this gives D > 0 and k− < ŵ < k+, i.e. R(ŵ) > 0, which is (b'). So
  k = min(ŵ, sL) is a subsolution from ẑ to z0, hence w > 0 there and f2'(z0) > 0: P2.
* EVENT LEMMA, all regimes, for every subextremal KdS and every m ≥ 1 on the window:
  R1/R2 (item 47, needs only L(z0) > 0: item 133); R3: P1 (item 135), P2(a) (item 149), P2(b') (this item).
  Computer inputs are all exact: the margin certificate L(z0) ≥ A/(2z0) + 2t (items 137/149), plus the covering
  certificates above.
* BUG NOTE: fastbern.half had the shift/scale order reversed for right halves (it computed p(x/2 + 1) instead of
  p((x+1)/2)). Fixed now. No earlier certified claim used subdivision; items 133–149 all passed in a single box.
  Only the diagnostic zhat2_sub.py "False" result of item 148 is void.

## Item 158 — Cosmological lemma for the bent q: structure, and a physical verification with the crude B
* Setting. f+ = 1, q̃ = (z − z0)Q with Q = q1 + κ(t − z), κ ≥ 0. Changing q changes L, since q²/μ − L is fixed:
  L̃ = L + κ(z−z0)²(2q1 + κ(t−z))/(z(z−1)). Conditions on (z0, t):
     C̃ = μ'q̃ − 2μq̃' ≤ 0,   D̃ = q̃L̃' + ½(μq̃'')' = q̃L̃' − κμ' ≤ 0,   and L̃(t) > 0 at the endpoint.
  (For κ = 0 these reduce to item 19: C+ = −q1H ≤ 0 by Lemma H, and D+ = qL' ≤ 0.)
* [P] MONOTONICITY. C̃ and D̃ do not involve B (L̃' does not). L̃(t) = L(t) + 2κq1(t−z0)²/(t(t−1)) is increasing in B.
  So the crude angular bound B_crude suffices, exactly as on the event side.
* [P] C̃ = −QH + 2κμ(z−z0) = −q1H − κE(z), with E = (t−z)H − 2μ(z−z0) = (t−z)·F and
  F = z0(2tz − t − z²) − tz − z³ + 2z², which is affine in z0. F(z0 = z) = 2z(z−1)(t−z) > 0 and
  F(z0 = 1) = (z−1)(t − z²). So E ≥ 0, hence C̃ ≤ 0 for every κ ≥ 0, unless z0 is close to 1 with z² > t. In that
  case C̃ ≤ 0 needs κ/q1 ≤ H/|E|, which is a bound on the bend relative to q1.
* [N] PHYSICAL check (nearext/cos174.py, 244834 points, m = 1..39, crude B): q1 > 0 always. L(t) < 0 occurs at
  1923 points. With the bend κ = 1.05 × (minimal κ making L̃(t) > 0), C̃ ≤ 0 and D̃ ≤ 0 hold on ALL of them. With κ = 0
  everywhere else, the old item 19 applies.
* REMAINING. A proof needs two physical relations: (i) κ/q1, i.e. (−L(t))·t(t−1)/(q1²(t−z0)²), against H/|E| in the
  region z0 ≈ 1, z² > t; (ii) D̃ ≤ 0, where the correction term in L̃' competes with L' = −A/z² − 3 and with κμ'
  (μ' < 0 near t). Both are certificate targets in the window variables (q1 is explicit there:
  q1 = −(t−1)Im(e1)/(2(z0−1))).

## Item 159 — Cosmological bend: abstract reductions fail; physics is needed in D̃
* [P] At z → z0, q̃ → 0, so D̃(z0) = −κμ'(z0). Hence the bend needs μ'(z0) ≥ 0: z0 must lie left of the maximum of μ.
  [N] Physically μ'(z0)/(z0t) ≥ 0.687 at every point needing a bend (these have zf ≤ 0.35; nearext/cos176.py).
* [N] Physical ratios at the bend-needing points (κ = minimal bend):
  κ(t−z0)/q1 ≤ 0.60;  κt(t−1)/(q1(t−z0)(A/t²+3)) ≤ 0.18;  At/(q1z0)² ≥ 10⁻³;  −L(t)/t up to 500.
* [N, NEGATIVE] On the abstract set with ALL these constraints imposed, C̃ ≤ 0 always holds but D̃ ≤ 0 fails on
  8–15% of samples (step175/177/178.py). The failures have large q1 against small A and t ≫ z0. So D̃ needs a
  physical relation that these scalar ratios miss.
* [P] Useful identity for the next step: A = t(X0² − Y²) and q1 = tY/z0 with Y = X1 + Xt (the residue of q² − P4 at
  z = 0), so q1 and A are tied through X0. A certificate in the window variables (with κ = c·κ_min built in) is the
  route. The bend is needed only for small zf, m ≥ 2, near-extremal.

## Item 160 — Cosmological side: sign-changing f+ instead of bending q. NEGATIVE
* Idea: with the currents of Prop. "Reduction", J^{f+}(t) = L(t)f(t)q(t)|B_t|². If L(t) < 0, take f+ with f(z0) > 0,
  f'(z0) ≤ 0 and f(t) < 0. Then J(t) > 0 when B_t ≠ 0, monotonicity forces Re(u'ū)(z0) < 0, and that contradicts the
  event side. This would need C^f ≤ 0 and D^f ≤ 0 on (z0, t).
* [N] On 300 physical bend-needing points, f = 1 − c(z−z0)^p (p = 1, 2, 3, 4, 6; f(t) = −δ) satisfies C^f, D^f ≤ 0 at
  NONE of them (nearext/cosf180.py). Near t both terms of D^f = q(f'L + fL') + ½q1(μf')' turn positive when L, f < 0.
  The bent-q construction (items 158–159) remains the route.

## Item 161 — Cosmological D̃ as a LOWER BOUND on A: physical margin 1.44
* [P] D̃ is affine in A with coefficient −q̃/z² < 0 on (z0, t). So D̃ ≤ 0 ⇔ A ≥ A_req(t, z0, q1, κ; z), and the
  cosmological lemma needs A ≥ max_z A_req.
* [N] Physically, with κ = 1.05·κ_min (κ_min from L(t) at B_crude), A/max_z A_req ≥ 1.44 on all 1742 bend-needing
  points (5th percentile 2.55, median 31); at 8% of them A_req ≤ 0, so D̃ ≤ 0 holds for every A ≥ 0
  (nearext/cos181.py).
* m-structure: q1 ∝ m, A ∝ m², B_crude = m²b2 + mb1 + b0, κ ∝ (−L(t))/q1 = O(m). So D̃ is a low-degree polynomial in m
  with coefficients in (x, v, τ, z), and the certificate has the same shape as items 148/157 (leading coefficient,
  plus endpoint values at m = 2, 3, plus concavity), with the gate L(t) < 0 and the √D reduction.
* This is the last mathematical gap in Theorem A.

## Item 162 — Cosmological certificate: exact algebraic data assembled
* [P, numerically confirmed to 1e-15] q1 = tY/z0 with Y = X1 + Xt. In the window variables
  Y = −Ξ·a·m·(f1 + f2)/λ̃ (f1 = U1 + U3 < 0, f2 = U2 + U3 < 0, λ̃ = Λ/3; checked against the direct κ_j formulas,
  cert3/cos_dbg.py). Only even powers of a appear: q1² = t²Ξ²b·m²(f1+f2)²/(λ̃²z0²), polynomial in b.
* [P] Write K := κq1 = 1.05·(−L(t))·t(t−1)/(2(t−z0)²), which contains no q1, and Q2 := q1². Then
  q1C̃ = μ'(z−z0)(Q2 + K(t−z)) − 2μ(Q2 + K(t−z) − K(z−z0)), and
  q1³D̃ = E/(z²(z−1)²), where E is a polynomial whose (A, B, Q2, K)-monomials are only
  K³, Q2K, Q2K², Q2², Q2²K, AQ2K, AQ2² (step182.py, P182.pkl). B enters only through K.
* m-structure: A, Q2 ∝ m²; K is quadratic in m (through L(t) at B_crude). E has degree ≤ 6 in m. The bend gate
  L(t) < 0 is a quadratic in m; it holds only on a bounded m-range for each geometry, because the m² coefficient
  Â/t + b2 of L(t) is positive.
* Certificate to build: in the variables (x, v, τ) plus w (z = z0 + (t−z0)w), for m in the gate range, show E < 0 and
  q1C̃ ≤ 0. Bounded m-range plus the polynomial structure in m means the m = 2..m_max(geometry) cases can be handled
  through u = 1/m ∈ [0, 1/2] with the gate. Expected size: the largest certificate so far (5 variables, √D reduction).

## Item 163 — Cosmological certificate: the gate L(t) < 0 is unavoidable; m-range unbounded
* [N] With the minimal-type bend K = 1.05(−L(t))t(t−1)/(2(t−z0)²) evaluated without the gate, E > 0 or q1C̃ > 0 on
  most points with L(t) > 0 (a negative "bend" is bad). On the gate L(t) < 0: 0 failures in 597 points
  (nearext/cos183.py).
* [N] No constant-coefficient multiplier −E + c·L(t)·W ≥ 0 exists for the normalisations W ∈ {q1⁴t², q1⁴t(t−z0)z0,
  q1⁴t³/z0}: the required c ranges are disjoint by about 30 orders of magnitude (cos184.py).
* [N] A universal bend K = θq1²/(t−z0) (no gate) fails broadly for every θ ∈ [0.01, 1] (cos185.py).
* [N] The bend is needed for m up to at least 103 (r0 → 1; nearext sampling with m ≤ 200). So the m-range is bounded for
  each geometry but not uniformly.
* So the certificate must handle the gate by COVERING: on each box of (x, v, τ, u = 1/m) either L(t) > 0 is certified
  (no bend; item 19 applies) or E < 0 and q1C̃ ≤ 0 are certified for all w ∈ (0, 1). Every sign is a b-expression
  (F1·b + F0), decided per box by the (F1, Y, P) rule.

## Item 164 — Crude angular bound PROVED by hand
λ ≥ λ_cr = ℓ0(ℓ0+1) − 4 − 4α − 4αm(1+α) − 2(1+α)²ξ(m+2), for s = ±2, m ≥ 1, α ≥ 0, ξ ≥ 0.
H = H0 + K_α + V with K_α ≥ 0 and H0 ≥ ℓ0(ℓ0+1) − s². Exactly (sympy, 0 residual),
V = α(m+s)²/(2(1−x)) + α(m−s)²/(2(1+x)) + 𝒩/(1+αx²). And 𝒩 − T(1+αx²), with T = λ_cr − ℓ0(ℓ0+1) + 4, is a sum of
manifestly nonnegative terms on |x| ≤ 1, |s| = 2 (regrouping verified exactly for s = ±2). Written into the paper as
Prop. "Crude angular bound". [N] Galerkin λ − λ_cr ≥ 2·10⁻⁴ on 1561 points.

## Item 165 — Referee issue (1), λ conventions: RESOLVED
* λ_ours (lowest eigenvalue of H = H0 + K_α + V; equal to ℓ(ℓ+1) − 4 at a = 0) is the λ in B = λC0 + c1. The c1
  identity uses the radial potential "+ s(1−α) − λ", which for s = +2 and a = 0 gives Teukolsky's A = (ℓ−2)(ℓ+3).
* step126's SdS radial uses "−λ_script" with no s(1−α) term, so λ_script = λ_ours − 2. Its constant
  [(λ_script+4)(λ_script+6)]² + 144M²σ² is therefore [(λ_ours+2)(λ_ours+4)]² + 144M²σ² = [(ℓ−1)ℓ(ℓ+1)(ℓ+2)]² + 144M²σ².
  This is the standard TS constant and matches the referee's formula.
* Paper fixed: C_ℓ and σ_ℓ are now written in λ_ℓ = ℓ(ℓ+1) − 4 (our convention).

## Item 166 — Referee issue (3), adjoint identity: PROVED by hand (paper updated)
O'^† = (Þ − ρ̄)(Þ' + 3ρ') − (ð − τ̄')(ð' + 3τ') + V0, from (Þ − ρ)^† = −(Þ − ρ̄) and (Þ' − 4ρ' − ρ̄')^† = −(Þ' + 3ρ')
(and primes/ð analogues). With ÞP = −4ρP etc. (P = Ψ2^{−4/3}, Bianchi ÞΨ2 = 3ρΨ2 …) this gives O'^†∘P = P∘O exactly.
No sampling is needed.

## Item 167 — Referee issue (11), HPV's T: addressed (paper Remark "Choice of T")
HPV25 (arXiv 2508.06620, eq. (9)) use T' = ∂t* + a/(r0²+a²)∂φ*, with r0 a critical point of μ/(r²+a²).
T' − T = c∂φ*, which commutes with everything. To transfer the statement:
* Gauge-fix the T'-mode (HV18 §10).
* The generalized resonant space of the gauge-fixed operator is finite-dimensional and ∂φ*-invariant, so the mode is a
  finite sum of m-components.
* Each component is a T-mode at σ' − cm, with the same imaginary part. For generalized σ ≠ 0 use HPV Thm 2.7's
  reduction. When σ' − cm = 0 with m ≠ 0, the parameter part vanishes (axisymmetry).
To double-check against HV18 §10 when re-refereeing: whether the gauge equation is solvable on T'-modes.

## Item 168 — Referee issues (15) and (12): FIXED in the paper
* (15) Generalized-mode induction, part (2). The gauge vector of the stationary step is X = X0 + t*K, with X0 stationary
  and K Killing (from Step 5). The correct correction is Z = t*^{k−1}X0/(k−1)! + t*^k K/k!, so that
  (ad ∂t)^{k−1}Z = X. The charge argument now also justifies c_M = c_a = 0 at the top order:
  δH_K[(L∂t)^{k−1}h0] = d/dτ of a conserved charge = 0.
* (12) Gluing across horizons. δ(Ψ2φ³) ≡ 0 because φ := Ψ2^{−1/3}, so it cannot fix the constants. They are now matched
  with the Iyer–Wald charges: δQ_K − K·Θ is closed and smooth on M, so δH_K[h] is the same on spheres in all three
  regions, and it vanishes on L_Xg. The nonzero Jacobian of (E, J) in (M, a) then gives (c_M, c_a) equal on all regions.

## Item 169 — Referee issue (4), Ψ4 bookkeeping: rewritten in the paper
* Convention h_xy := h_ab x^a y^b. S'^† produces h_Ψ = X n⊗n + Y n⊙m̄ + Z m̄⊗m̄, whose nonzero components are h_ll,
  h_lm, h_mm. Linearized Ψ4 uses only h_nn, h_nm̄, h_m̄m̄, so Ψ4[h_Ψ] = 0. For h̃ = conj(h_Ψ), Ψ4[h̃] comes from
  h̃_m̄m̄ alone. The result is Ψ4[h̃] = c0·Þ'^4(P̄Ψ̄) with a universal c0 ≠ 0 (½ in KSFK), and only c0 ≠ 0 is used.
  Wald's step is also spelled out: E S'^† = T'^† O'^†.
* To recheck at the re-referee: the claim that the Λ-terms do not enter (only Þ' and the Ricci identities Þ'ρ' = ρ'²
  are used). The value of c0 itself is irrelevant.

## 170. Euler transform (referee issue 5) — rewritten
The real RL integral is always divergent at one endpoint (Re ν = s−1, Re ρ₁ = −s). Replaced it with the Pochhammer loop γ_z around {1, z}. On powers the normalisation c(β) ∝ 1/(Γ(−β)Γ(−ν)Γ(β+ν+2)) is entire and nonzero for β = ρ₁+k in the window (ν ∉ Z, since κ₃ ∉ R). Threshold κ₁ = 0: the analytic part is killed by Pochhammer, so use a single loop around w = z (s = −2, or s = +2 with a₀ = a₁ = 0), and Pochhammer on the polar part otherwise. Boundary terms vanish because the contour is closed. Status: fixed in the paper.

## 171. Axis/acceleration step (referee issue 13) — rewritten
Replaced the unverifiable determinant "−64Mp₀⁴/c" with a triangular system.
(i) Misner string: the axis generators are ∂ψ + I±²∂t̃. Under δΔ_I = 2δn·I, both roots shift by δn/(c+λ̃a²) (sympy-checked), so the generators differ by 4a·δI·∂t̃. Since t̃ is not periodic, δn = 0 (a ≠ 0; at a = 0 it is the classical NUT case).
(ii) C-metric conical condition at first order gives 4M·δα = 0.
Diagonal entries 4a/(c+λ̃a²) and 4M. Residual soft spot: the claim that the acceleration tangent doesn't move the axis generators at first order is asserted from Griffiths–Podolský, not recomputed.

## 172. Cosmological side, C̃ ≤ 0: the w-variable ELIMINATED (proved by hand + positive-coefficient identity)
Setting: q̃ = (z−z0)(q1 + κ(t−z)), with H = 2μ − μ'(z−z0) and F as in item 158, so q1C̃ = −Q2·H − K(t−z)F (K = κq1, Q2 = q1²).
* [P] Substitute z0 = 1+a, z = 1+a+b, t = 1+a+b+c (a, b, c > 0). Then
  H = a²b + 2a²c + 2ab² + 2abc + ab + 2ac + b³ + b² + bc > 0, and
  G := (t−z0)H + (t−z)F = a²b² + 2a²bc + 4a²c² + 2ab³ + 2ab²c + ab² + 4abc² + 2abc + 4ac² + b⁴ + b³ + b²c + 2bc² ≥ 0,
  all coefficients positive (sympy). So (t−z)(−F) ≤ (t−z0)H on 1 < z0 < z < t. Numerically sup = 0.895.
* [P] Hence, wherever F < 0, q1C̃ ≤ −H(Q2 − K(t−z0)); wherever F ≥ 0, q1C̃ ≤ −Q2H < 0.
  LEMMA: C̃ ≤ 0 on (z0, t) as soon as κ(t−z0) ≤ q1.
* So the C̃ half of hyp:cos reduces to ONE scalar inequality with no w:
      (S_C)   1.05·(−L(t))·t(t−1) ≤ 2 q1² (t−z0)   on the gate L(t) < 0.
  Physically κ(t−z0)/q1 ≤ 0.60 (item 159). This is a 4-variable (x, v, τ, 1/m) certificate, one dimension down.
* [N] Bend-shape variants φ = ((z−z0)/(t−z0))^p (nearext/cos190.py, 1257 gate points): C̃ ≤ 0 for every p, but
  D̃ ≤ 0 fails at 437/767/975 points for p = 1/2/4. p = 0 (the linear bend) has 0 failures. Keep p = 0.
* Reformulation used for intuition: C̃ = −2μ^{3/2}(q̃/√μ)', so C̃ ≤ 0 ⇔ q̃²/μ = V + L̃ is nondecreasing.
Remaining for hyp:cos: (S_C) [4D, scalar] and D̃ ≤ 0 [5D, still needs w].

## 173. Cosmological side, D̃: exact structure after the (a, b, c) substitution
z0 = 1+a, z = 1+a+b, t = 1+a+b+c; E = q1³D̃·z²(z−1)² (P182), r := K/Q2 = κ/q1. [P, sympy]
  E =  −A·b(a+b)²·Q2(Q2 + cK)                                   [A-term, < 0]
       −3b·z²(z−1)²·Q2²                                         [< 0]
       + z²(z−1)²·g·Q2K,  g = z(z−1) − (t−z)(5z − 3z0 − 1)      [indefinite]
       + b²·{ 2(2a²+2ab+2a+b)Q2²K  −  h₂·Q2K²  −  c·h₃·K³ }
  h₂ = a²b − 6a²c + 2ab² − 6abc + ab − 6ac + b³ + b² − 3bc,  h₃ = same pattern with 6→2, 3→1.
* Every non-A term except −3bz²(z−1)²Q2² carries a factor r. So D̃ ≤ 0 ⇔ A ≥ r·[α(a,b,c,r) + Q2·β(a,b,c,r)], where
  α comes from the g-term (with the −3b term subtracted) and β from the b²-bracket. Both are explicit.
* b → 0 (z → z0): E → z0²(z0−1)²·Q2K·[z0(z0−1) − (t−z0)(2z0−1)] = −z0(z0−1)·μ'(z0)·Q2K, which matches item 159
  (the bend needs μ'(z0) ≥ 0).
* Certificate target, now explicit: on the gate, with r ≤ 1/(t−z0) (from S_C), show
  sup_b [α + Q2β]·r ≤ A, with b ∈ (0, t−z0). This is the 5D statement, but in polynomial form linear in A and Q2, with
  the w-dependence confined to α and β.

## 174. D̃: split sup over b is nearly lossless; the bend ratio r must stay a variable
Divide E by Q2² (r = K/Q2, T = t−z0). Then D̃ ≤ 0 ⇔ A ≥ P1(b) + Q2·P2(b) for all b ∈ (0, T), where
  P1 = z²(g·r − 3b)/(b(1+cr)),   P2 = b·r(2p − h2·r − c·h3·r²)/(z1²(1+cr))   (z1 = z−1, c = T−b, p = 2a²+2ab+2a+b).
[N] (nearext/cos191.py, 1229 gate points):
* exact A/max_b(P1 + Q2P2) has min 1.468 (reproduces item 161's 1.44, which validates the decomposition);
* SPLIT bound A ≥ α + Q2·β⁺ with α = sup_b P1 and β = sup_b P2 has min margin 1.371, so it is nearly lossless;
* replacing r by the S_C cap 0.7/T DESTROYS it (min 0.001). The actual r is usually far below the cap, so r must
  remain a certificate variable.
REDUCED TARGET for hyp:cos (D̃ half): A ≥ α(a, T, r) + Q2·β(a, T, r), where α and β are 3-variable functions defined by a
one-dimensional sup. Next step: replace α and β by explicit polynomial upper bounds in (a, T, r), via Bernstein on b
with a margin slack of about 1.37, then certify the physical inequality in the window variables (4D, no w).

## 175. hyp:cos reduced to three scalar conditions plus one 1-variable sup (exact identities, sympy-checked)
Notation: a = z0−1, T = t−z0, b = z−z0, c = t−z, r = κ/q1 = K/Q2, p = 2az + b, z1 = z−1.
[P] Identities:  g = b(2+4a+b−5c) − M0, where M0 := T(1+2a) − a(1+a) ∝ μ'(z0);
                 2p − h2 r − c h3 r² = (1+rc)·[p(2+rc) − r b z z1].
Hence  P1 = z²(gr − 3b)/(b(1+cr)) ≤ z²[r(2+4a+T) − 3]/(1+cr) whenever M0 ≥ 0, and
       P2 = b r [p(2+rc) − r b z z1]/z1²   EXACTLY (the (1+cr) cancels).
THEOREM (reduction). On the gate L(t) < 0, hyp:cos (C̃ ≤ 0, D̃ ≤ 0, L̃(t) > 0 with κ = 1.05 κ_min) follows from:
  (S_C)  rT ≤ 1                          [item 172]
  (S_μ)  M0 = T(1+2a) − a(1+a) ≥ 0        [μ'(z0) ≥ 0, item 159]
  (S_α)  r(2+4a+T) ≤ 3                   [kills the P1 part: α ≤ 0]
  (S_β)  A ≥ Q2 · sup_{0<b<T} b r[(2az+b)(2+r(T−b)) − r b z(z−1)]/(z−1)²
[N] (cos192/193): (S_μ) has 0 violations in 1224 gate points; with (S_α), α ≤ 0 holds at every point; (S_β) has min
margin 1.40. Dropping the −r b z z1 term in (S_β) loses a factor of up to 1.3e5. It is essential, since it is the
self-limiting part of the bend.
Remaining: certify (S_C), (S_μ), (S_α), (S_β) in the window variables (x, v, τ, 1/m); (S_β) carries one sup over b.

## 176. (S_μ) reduced to "gate ⇒ zf ≤ 1/2"
[P] The max of μ = z(z−1)(t−z) is at z* = [(1+t) + √(t²−t+1)]/3. Moreover z* − 1 ≥ (t−1)/2 ⇔ √(t²−t+1) ≥ (t+1)/2
⇔ 3(t−1)² ≥ 0. So μ'(z0) ≥ 0 whenever zf := (z0−1)/(t−1) ≤ 1/2 (and zf* runs from 1/2 at t→1 to 2/3 at t→∞).
In window terms, zf ≤ 1/2 ⇔ q1v ≤ −qtv ⇔ Im e1 ≤ t·Im(−et)  (e1 = 2(κ1+κ3), et = 2(κ2+κ3)).
[N] cos194.py: 268158 physical points, 648 on the gate, max zf on the gate = 0.356 (margin 0.14 below 1/2).
REMAINING for S_μ: certify "L(t) < 0 ⇒ zf ≤ 1/2". Equivalently, zf > 1/2 ⇒ L(t) ≥ 0, which is a statement with no b
and no bend.

## 177. (S_μ) target "zf > 1/2 ⇒ L(t) ≥ 0": genuinely gate-linked; healthy margin
[N] zf is m-independent (κ ∝ m cancels in z0), and zf > 1/2 on 32% of the window (max 0.99998). So S_μ cannot be
proved without the gate.
[N] On zf > 1/2, normalized L(t)/(A/t + |B| + 3t) ≥ 0.217 (77135 points, m ≤ 80). The worst case is m = 2. For
zf > 0.45 the bound is 0.172, and for zf > 0.4 it is 0.074.
PLAN: the same machinery as the margin certificate (cert3/margin_red.py: variables r0, r2, b = a², τ; b reduced mod
b² − e2b + e4; L(t) as LQ2·m² + LQ1·m + LQ0 with the crude λ). Add the constraint 2(z0 − 1) ≥ (t − 1) as a polynomial
sign on the box (clear the denominators of z0 = t/R), and certify LQ(t) ≥ 0 on boxes where the constraint can hold.
Margin 0.2 means the subdivision depth should be similar to fullcert_m.

## 178. S_μ certificate: build data and a NEGATIVE on multipliers
* [P] Window expression for λ̃·L(t) (same normalization as margin_red; m ≥ 2, crude λ):
    LQ2' = 4Ξ²b·U0U3/λ̃ + K(1 − 2Ξ²b·y),   LQ1' = LQ1 (unchanged),   LQ0' = −K(1+7α) + λ̃(1 − 2t).
  Derived from margin_red via λ̃A/(t m²) = (2z0/t)·[λ̃A/(2z0 m²)] and λ̃c1 = λ̃(1+t) + 3(1−α)K.
* [N, NEGATIVE] S-procedure L(t) + β(1 − 2zf)·S ≥ 0 on the whole window, with S ∈ {A/t, t, A/t + 3t}: the required β
  ranges are disjoint (S = A/t+3t: need β ≥ 6.9, allowed β ≤ 0.67). So the constraint zf ≥ 1/2 has to enter through
  box COVERING (on each (x, v, τ) box, either 2(z0−1) < t−1 is certified, or Q2'/H2'/H1' ≥ 0), not through a
  multiplier.

## 179. S_μ covering certificate: built, validated, 99% of the window closed
* cert3/smu_build.py: exact window polynomials for λ̃L(t) (Q2' = m² coeff, H2' = value at m = 2, H1' = m = 1) and
  G = 2(z0−1) − (t−1), b-reduced to (F1, Y, P). smu_chk.py: H2' and z0 match the direct pdq numerics to 1e-13 at 3
  points. Concavity in 1/m is automatic (LQ0' = −K(1+7α) + λ̃(1−2t) < 0).
* Factorisation: every odd-multiplicity factor except the cores is linear and sign-definite on the cube. Cores: G has
  (2,2,1)², (4,4,2) in Y, (5,5,2)²; L has (24,24,2) in F1, (26,26,2) in Y, (24,26,4) in P.
* Sign rule at b₊: value > 0 iff (P > 0, F1 > 0) or (P < 0, Y < 0); value < 0 iff (P > 0, F1 < 0) or (P < 0, Y > 0).
  Sampling shows G < 0 is realised through BOTH patterns, so both are needed.
* cert3/smu_cover4.py (resumable via margin/smu_open.json): Bernstein covering in (x, v = 1/r2, τ). A box closes if
  G < 0 is certified or Q2', H2', H1' > 0 are certified (F1/P rule). Result: open volume 0.0110 (122 boxes at total
  level ≤ 15). 90 boxes are interior, 24 touch x = 1 (r0 → 1, near-extremal), and a few touch τ = 1.
* Adding the (P < 0, Y < 0) branch for L (smu_cover5.py) makes each box far slower (the Y cores have degree 26), so it
  closed nothing in 45 s.
* STATUS: S_μ certificate ~99% done. The remaining interior boxes need depth; the x = 1 boxes need an asymptotic
  (near-extremal) argument, since a Bernstein covering cannot close on a face where G and L degenerate.

## 180. S_μ near-extremal face x = 1: first look
* [P] Each target has the linear factor (r0 − 1) (sign-definite on the open cube, removed). The cores do NOT vanish
  identically on x = 1.
* [N] On the face, the H2 cores of F1 and Y have constant sign (−). The P core of H2 and all three G cores change sign
  on the face. So the face boxes straddle the curves {G-core = 0} and {H2-P core = 0} on x = 1. Whether these cross
  inside {zf ≥ 1/2} is the question.
* [N, unreliable] A float check of the b₊-values near x = 1 is dominated by cancellation (|terms| ~ 1e50). Exact or
  mpmath evaluation is needed, and the straightforward mpmath loop is too slow (> 60 s for 625 points). Next: restrict the
  cores to x = 1 symbolically (2-variable polynomials in v, τ) and decide the face question exactly (resultant of the
  G-core and the H2-P core on the face, plus sign checks on each arc).

## 181. S_μ covering: everything certified except a corner neighbourhood
* [P, exact grid] On the face x = 1 − 10⁻² and x = 1 − 10⁻⁴ (29×29 rational grid in (v, τ), exact fmpq evaluation with
  the b₊ sign rule): zf ≥ 1/2 ⇒ Q2', H2', H1' > 0 at every point. There is no genuine failure on the near-extremal face.
* New box rules (cheap, no P needed): F1 > 0 and Y < 0 ⇒ value > 0; F1 < 0 and Y > 0 ⇒ value < 0. This works since
  value = (F1√D − Y)/2.
* cert3/smu_cover7.py (resumable, polys persisted as dicts, largest-first order, atomic state writes):
  open volume 1.1e-2 → 2.7e-9 over ~11 runs of < 1 min. The 23 remaining boxes all lie in
  x ≥ 0.996, v ≤ 0.0015, τ ≥ 0.9985, i.e. a neighbourhood of the CORNER (r0, 1/r2, τ) = (1, 0, 1):
  extremal, Λ → 0 (Kerr), top of the window.
* REMAINING for S_μ: a corner lemma (blow up at (1, 0, 1), e.g. x = 1 − s·p, v = s·q, τ = 1 − s·r with
  p + q + r = 1, and certify at s → 0 plus a Bernstein check on the blown-up chart). Everything else in the window is
  certified.

## 182. (S_μ) PROVED (certificate): on the gate L(t) < 0, zf ≤ 1/2, hence μ'(z0) ≥ 0
* Corner lemma, cert3/smu_corner.py: the corner cube [1−1/256, 1] × [0, 1/256] × [1−1/256, 1] in (x, v, τ) is covered
  by three projective charts. A: x = 1−σ/N, v = σβ/N, τ = 1−σγ/N. B: v leading. C: 1−τ leading. N = 256 and
  (σ, β, γ) ∈ [0,1]³. The cores are substituted exactly, powers of σ are removed, and each chart is covered with the
  same box rules. Results: B closed (8 boxes), C closed (1 box), A closed after resumable runs (~80 boxes, max level 8).
* On the σ = 0 face of chart A, the L-cores give F1 > 0, Y > 0 with the H2-P core changing sign (so L(t) < 0 does occur
  at the extremal-Kerr corner). On that set G < 0 at every sampled face point. That is consistent: the gate lives at
  zf < 1/2.
* Main covering (smu_cover7): every box outside the corner cube is closed. The last 23 open boxes lie in
  x ≥ 0.996 > 1−1/256, v ≤ 0.0015 < 1/256, τ ≥ 0.9985, i.e. inside the corner cube.
* Ingredients: exact window polynomials validated against direct numerics (item 179); b₊ sign rules
  (F1, Y, P); concavity in 1/m (LQ0' < 0) ⇒ checking Q2', H2' (m = 2) and H1' (m = 1) covers all m ≥ 1.
=> S_μ holds on the whole window. Remaining for hyp:cos: S_C, S_α, S_β.

## 183. S_C and S_α: gate-free reformulation, polynomials built, S_C main cube certified, ultracold edge open
* [P] The S_C and S_α conditions hold automatically where L(t) ≥ 0 (their left sides are then ≤ 0). So they are
  equivalent to GLOBAL inequalities with no gate:
    W_C := 2λ̃Q2·T − 1.05·t(t−1)·(−λ̃L(t)) ≥ 0,   W_α := 6λ̃Q2·T²/(2+4a+T) − 1.05·t(t−1)·(−λ̃L(t)) ≥ 0,
  with λ̃Q2/m² = t²Ξ²b(f1+f2)²/(λ̃z0²). Both are quadratics in m with the concave L-part, so checking Q2 (m²), H2 (m = 2)
  and H1 (m = 1) suffices. cert3/sca_build.py; validated against direct numerics to 1e-12 at 3 points (sca_chk.py).
* cert3/gcover.py: generic resumable covering (box index mode, lazy descend, largest-first, depth-6 cap per subtree,
  cached factor cores). W_C on the main cube: ALL CLOSED outside the edge region x, v ∈ [7/8, 1].
* Edge x, v → 1 (r0, r2 → r1: ULTRACOLD triple-root limit). Chart EA (x = 1−σ/8, v = 1−σβ/8) with sub-charts at
  (σ, τ) → 0: EA1 closed. EA2 (τ = s/8, σ = sg/8) is stuck, because the leading forms of the F1/Y/P cores CHANGE SIGN
  in g there. Most likely the discriminant D = e2² − 4e4 → 0 at the triple root, so the (F1, Y, P) sign rule loses
  resolution.
* [N] Physics is fine there: 224546 points with r0 = 1−ε0, r2 = 1+ε2 (ε up to 0.3), normalized W_C ≥ 0.89 (797 gate
  points). So the obstruction is in the b-reduction, not in the inequality.
* NEXT: in the ultracold edge, either (i) avoid the b-reduction by working with b as a variable (b₊ is then a smooth
  function near the triple root?), or (ii) factor the leading D-power out of F1, Y, P in the EA2 chart. Then S_α
  (same machinery, polynomials ready) and S_β.

## 184. Ultracold edge for W_C: diagnosis
* D = e2² − 4e4 does NOT vanish at the triple root (r0 = r2 = 1, r3 = −3: e2 = −6, e4 = −3, D = 48). So the b₊
  sign rule is not the problem; item 183's guess was wrong.
* Chart EA (x = 1−σ/8, v = 1−σβ/8, τ free): the σ = 0 face is the single physical edge x = v = 1, but the chart
  polynomials depend on the approach direction β. Open boxes cluster (a) along the whole σ = 0 face and (b) in a band
  β ∈ [0.5, 0.75], τ ∈ [0.25, 0.5] extending into σ > 0. This suggests that in the limit, value(W_C) vanishes to leading
  order along a curve β(τ), i.e. a genuine second-order degeneracy (both sides of S_C → 0 at the same rate as
  r0, r2 → r1).
* Comparing leading forms needs a consistent σ-order for F1 and Y (the per-core power stripping in gcover.chart is
  sign-correct for box tests, but it hides the relative orders). Not done yet.
* [N] Physically W_C is fine near the edge (normalized ≥ 0.89, item 183).
* NEXT: a direct asymptotic expansion of W_C at the ultracold edge (r0 = 1−ε·p, r2 = 1+ε·q), using the Nariai/ultracold
  scaling of t − 1, T, q1, L(t) in ε, to get the leading coefficient as an explicit function of (p/q, τ). Then a
  2-variable certificate for that leading coefficient, plus a remainder bound.

## 185. Ultracold asymptotics of W_C: the q1² term dominates by a full power of ε
[N] r0 = 1−ε·cos(θπ/2), r2 = 1+ε·sin(θπ/2), m = 2, several τ: t−1 stays O(1) (it depends on the direction θ only),
2λ̃q1²T ~ ε⁻², while 1.05·t(t−1)·(−L(t)) ~ ε⁻¹, so their ratio → 0 linearly in ε (≤ 0.048 at ε = 1e-2, ≤ 5e-4 at
ε = 1e-4). L(t) < 0 does occur there (θ = 0.8, τ = 0.95), but W_C ≫ 0. So S_C is NOT degenerate at the ultracold
edge; the covering is stuck only on the chart's sub-faces:
  β = 0 (v = 1, the Nariai face r2 = r1, where t → 1) meets σ = 0, and τ = 0 (bottom of the window, q1 → 0) meets σ = 0.
EA2 (τ-blowup) is slow (degree inflation) but closing; EA main stalls at the (σ, β) = (0, 0) corner, which needs one
more blowup or a hand bound.
HAND ROUTE (preferred): near the edge prove λ̃q1²T ≥ c·ε⁻²·(τ-factor) and −λ̃L(t) ≤ C·ε⁻¹, with explicit window
formulas; then W_C ≥ 0 for ε ≤ ε0 with explicit ε0 ≥ 1/8 needed to match the chart.

## 186. W_C ultracold region: covering does not converge efficiently
* gcover.py now supports nested blow-ups (CH = 'EA:ijl,...') and a variable edge-chart scale. With the exclusion
  square shrunk to [31/32, 1]², the main-cube covering stalls in thin strips along x ≈ 0.98, v ∈ [0.92, 0.96] and
  mirror (box count 97 → 411 over 4 runs; open vol 4.0e-3 → 2.3e-3). The degeneracy is not confined to a corner square:
  it runs along the edges x = 1 and v = 1 near the corner.
* The EA chart (scale 1/8) is slow (degree inflation; open vol 0.4375 → 0.421).
* CONCLUSION: brute covering near the ultracold corner costs too much for < 1 min runs. Use the hand route of item 185:
  near r0, r2 → r1, 2λ̃q1²T ~ ε⁻² beats 1.05·t(t−1)(−λ̃L) ~ ε⁻¹. Write explicit bounds in a neighbourhood
  {x, v ≥ 7/8}, which is where the main covering already closes everything else (item 183, scale-8 run).

## 187. S_C: exact z0-FREE form (identity, checked to 1e-14)
Throughout the window Im e1 < 0 and Im et < 0 (sampled), and zf = |Im e1|/(|Im e1| + t|Im et|),
q1 = −(t−1)Im e1/(2(z0−1)). Substituting gives
    rT ≤ 1   ⇔   2.1·(−L(t)) ≤ |Im et|·(|Im e1| + t|Im et|)          (S_C exact)
and, in the window, Im e1 = 2Ξ a m f1/λ̃, Im et = 2Ξ a m f2/λ̃ (checked to 1e-15).
So λ̃·[Im et(Im e1 + t Im et) + 2.1 L(t)] = m²·4Ξ²b f2(f1 + t f2)/λ̃ + 2.1(m²LQ2' + mLQ1 + LQ0'):
no z0, no R. cert3/scn_build.py, validated (scn_chk.py) to 1e-14.
[N] Margins on the gate (816 pts): exact 1.54; dropping the |Im e1| term 1.54; also dropping A ≥ 0, 1.43.
So on the gate, S_C ⇐ t·(Im et)² ≥ 2.1(3t − B). But that weaker form is NOT global (fails where L ≥ 0 only
because of A), so it would need the gate inside the covering.
Covering of the exact form (gcover, scN): it again stalls only at the ultracold corner (x, v → 1), and the
edge charts do not converge at depth 6 (Bernstein pessimism at degree ~56).
Scaling there: Im e1, Im et ~ 1/D'(r) blow up as roots merge, while L ~ ε⁻¹. The right side dominates by
powers of ε, so a hand corner lemma is the right tool. The z0-free form makes it tractable.

## 188. Ultracold corner for S_C: margin grows like 1/ε
[N] r0 = 1−ε cos θ, r2 = 1+ε sin θ, all τ, m ≤ 60, the sufficient form t·(Im et)² ≥ 2.1(3t − B) on the gate:
  ε ∈ [0.03, 0.1]: min ratio 96.6;  ε ∈ [0.01, 0.03]: 255;  ε ∈ [0.001, 0.01]: 1229.
The gate near the corner occurs only at the top of the window (τ ≈ 0.96) with θ ≈ 1.5 and m = 2.
So a corner lemma with very crude constants is enough (a factor of ~100 to spare at ε = 0.1).
Ingredients still to write: (i) a lower bound for |Im et| on the gate near the corner (it blows up like 1/D'(r2)
unless the numerator ω(r2²+a²) − am is small; S_μ, i.e. |Im e1| ≤ t|Im et| on the gate, controls that case);
(ii) an upper bound 3t − B ≤ C/ε from the crude angular bound and C0 = 3/(Λ(r1−r0)(r2−r3)).

## 189. Ultracold corner: scaling structure (toward a hand lemma)
Im κ_j = Ξ m (r_j²+a²)(w − Ω_j)/D'(r_j), Ω_j = a/(r_j²+a²), w ∈ [w2, w1],
w1 = (P1Ω1 + P3Ω3)/(P1+P3), P_j = (r_j²+a²)/|D'(r_j)|.
Corner limit r0, r2 → 1 (r3 → −3): 3λ̃² + 6λ̃ − 1 = 0, so λ̃ = 2/√3 − 1, a² = 3λ̃.
D'(r0), D'(r1), D'(r2) ~ ε²; D'(r3) = O(1); Ω1 − Ω2 ~ ε; w1 = Ω1 + O(ε²).
Hence at the top of the window |Im et| ~ (Ω1−Ω2)/D'(r2) ~ 1/ε, so t(Im et)² ~ ε⁻², while 3t − B ~ C0 ~ ε⁻¹.
That is the observed 1/ε margin growth (item 188).
OBSTRUCTION for a clean hand lemma: if w2 < Ω2 < w1, the numerator of Im κ2 vanishes inside the window, so |Im et|
can be small. Then t(Im et)² gives no lower bound, and one must show the gate cannot occur there (L(t) ≥ 0, through
the A term or B). This is a genuine multi-scale analysis (the three horizons merge at rates set by θ, and the window
endpoints move with ε).

## 190. S_C at the ultracold corner: GATE-FREE sufficient form via S_μ
[P] Put F := |Im e1| + t|Im et| = −(Im e1 + t·Im et), which is linear in w; in the window F = −2Ξam(f1 + t f2)/λ̃.
On the gate S_μ (PROVED, item 182) gives zf ≤ 1/2 ⇔ t|Im et| ≥ |Im e1| ⇒ t|Im et| ≥ F/2. Hence
|Im et|·F ≥ F²/(2t). With −L(t) ≤ 3t − B (since A ≥ 0), S_C (item 187) follows from
      (G_C)   F² ≥ 4.2·t·(3t − B),
which needs NO gate. This also closes the item-189 obstruction: Im e1 and Im et cannot both be small, because
their numerators vanish at w = Ω1 and w = Ω2 respectively.
[N] Near the corner (r0 = 1−ε cos θ, r2 = 1+ε sin θ, θ ∈ [0, π/2], all τ, m ≤ 60), the ratio F²/(4.2t(3t−B)) has
min 6.9 for ε ∈ [0.13, 0.25], 16 for ε ∈ [0.03, 0.13], and 69 for ε ≤ 0.03 (from 8.4-factor runs, ×2).
Globally (G_C) is FALSE (min 0.71 at 0.13% of points, away from the corner), so it is used only on the corner
square, where the W_C covering stalls; outside the square, W_C is certified (item 183).
cert3/scg_build.py builds the (G_C) window polynomials (degree ~56); edge-chart covering again slow.
NEXT: prove (G_C) on the corner square by hand. F and B are explicit and linear/quadratic in w, so the minimum over
w sits at a kink (w = Ω1 or Ω2) or at a window endpoint; this reduces to finitely many 2-parameter (ε, θ)
inequalities.

## 191. (G_C) on the corner chart: face definite, covering too slow
* In chart EA (x = 1−σ/8, v = 1−σβ/8), the σ = 0 face of every (G_C) core is sign-DEFINITE on the whole closed face,
  edges and corners included (F1 < 0, Y < 0, P < 0 ⇒ value > 0). So, unlike W_C, there is no degenerate sub-face,
  and a covering must terminate.
* Throughput is the obstacle: ~3 boxes closed per 42 s run (exact fmpz Bernstein, 9 cores of degree ≥ 56 after the
  chart substitution). Open vol 0.375 after one run.
* To finish this cheaply: (a) a hand argument using that F is linear and 3t − B is linear in w (so (G_C) is a convex
  quadratic in w, whose sign follows from a discriminant inequality in (ε, θ, m)); or (b) a faster engine
  (float Bernstein with outward rounding / interval arithmetic) for the chart covering.

## 192. (G_C) corner: τ eliminated exactly (quadratic in τ), 2D covering, still throughput-limited
* [P] Each (G_C) target (Q2, H2, H1), after clearing denominators, is a QUADRATIC q(τ) = c0 + c1τ + c2τ² with
  coefficients polynomial in (r0, r2, b) (cert3/scg_tau.py). So q ≥ 0 on [0, 1] ⇔ q(0) > 0, q(1) > 0, and
  [c2 ≤ 0 or c1 ≥ 0 or c1 + 2c2 ≤ 0 or 4c0c2 − c1² > 0]. Each condition is b-reduced to (F1, Y, P).
  This turns the corner problem into a 2D one, in (x, v) or the chart (σ, β).
* [N] On the chart face σ = 0 (the ultracold limit), q(0), q(1), the discriminant and c1 are all strictly positive
  for every β (sampled at 21 points). So the 2D covering must terminate.
* cover2d.py, chart EA: open vol 0.20 after one 40 s run (14 boxes closed). With only {q0, q1, v1} it is worse
  (0.48): the discriminant branch is needed in the interior. Throughput is again the limit (the disc core has degree
  ~110 before the chart substitution).
* The faster-engine route (interval float Bernstein) is now the practical bottleneck for finishing S_C. The
  mathematics is in place: S_C on the window = [W_C covering outside the corner square (item 183)] + [(G_C) on the
  corner square via S_μ (item 190)], and (G_C) on the square reduces to this 2D sign problem.

## 193. Fast float-Bernstein engine; (G_C) corner chart closes except two isolated spurious points
* cert3/fbern2.py: exact rational Bernstein coefficients (via the rev–shift–rev trick), normalised and converted to
  float64, then de Casteljau splitting. All operations are convex combinations, so the coefficient error after k
  splits is ≤ (k+4)·1e-13 relative to max|b|, and a sign is accepted only beyond that tolerance (rigorous up to
  float rounding, which is bounded). About 100× faster than the flint compose route.
* BUG FOUND AND FIXED: the τ-coefficients must be taken from the numerator n with the sign of d's constant
  (d = c·(square), and c < 0 for H1). Using n·d instead introduced spurious zeros from d. Now scg_tau.py asserts
  that d is a square and uses sign(c)·n.
* Extra closure branch: q(τ) ≥ c0 − |c1| − |c2| (the four sign combinations).
* cover2d_f.py, chart EA: ~20k boxes closed. The only remaining open boxes converge to two isolated points,
  (σ, β) ≈ (0.17658, 0.82755) and (0.29920, 0.05645). There, for the Q2 target (m² coefficient) at τ = 0, the
  b-reduced q0 has F1 = Y = 0, i.e. the numerator n vanishes at BOTH roots b±. Physically the m² coefficient there is
  ≈ 2.0e4 > 0 (and 6.9e3 at the second point), so these are removable zeros of n (n and d share a zero at b₊ at
  isolated points), not failures.
* NEXT: locate the shared factor (probably the y2 denominator ρ1Pa0 − ρ0Pa1 or a U-denominator) and cancel it
  symbolically before the τ-split, or handle the two points with a local chart. Then run chart EB, which mirrors EA.

## 194. (S_C) PROVED (certificate). Corner charts EA and EB both close.
* Root cause of the "stuck points" (item 193): the τ-condition polynomials (q0, q1, disc, c2neg, v1, v2) carry
  SQUARE factors in (r0, r2, b), e.g. (1,1,1,0)², (1,2,1,0)², (2,2,2,0)², which vanish at b₊ on curves. After
  b-reduction these become spurious near-zeros of F1, Y, P. Fix (scg_tau.py): factor each condition polynomial
  and drop the even-multiplicity factors before the b-reduction (sign unchanged).
* cover2d_f.py (float Bernstein, fbern2.py): chart EA closed (~21k boxes, 3 runs); chart EB closed (~24.5k boxes,
  3 runs). EA ∪ EB = the corner square {r0 ≥ 7/8, 1/r2 ≥ 7/8}.
* ASSEMBLY of S_C (rT ≤ 1 on the gate):
  - outside the corner square: the W_C covering, exact fmpz Bernstein (item 183);
  - inside it: S_μ (item 182) + A ≥ 0 + exact identity (item 187) ⇒ S_C ⇐ (G_C), and (G_C) ⇔ for each target
    (Q2, H2, H1) the τ-quadratic is ≥ 0 on [0, 1], certified via the closure rules above.
* Caveats for the re-referee: (i) the corner part uses the float engine (tolerance (k+4)·1e-13 relative; should be
  re-run in exact arithmetic or with interval rounding for the final version); (ii) A ≥ 0 on the window is used as in
  items 147ff.
Remaining for hyp:cos: S_α and S_β.

## 195. S_α: exact z0-free form, polynomials, 3D float covering; obstacle = catastrophic cancellation in the b₊ rule
* [P, identity checked to 1e-14] On the gate, with F = |Im e1| + t|Im et|:
    S_α ⇔ 2.1·(−L(t))·[2F + (t−1)(F + 3|Im e1|)] ≤ 3t(t−1)·(Im et)²·F.
  It is gate-free as W_α := 3t(t−1)(Im et)²F + 2.1·L(t)·G ≥ 0 (G > 0). After dividing by m·a it is a quadratic in m
  with the concave L-part, so Q2/H2/H1 suffice. Physical margin 4.77 (770 gate points).
  cert3/sal_build.py (odd parts of n and d only; squares dropped); validated (sal_chk.py) to 1e-12.
* cert3/cover3d_f.py: 3D float-Bernstein covering (fbern3.py), resumable, depth cap per subtree. It excludes the
  ultracold square and the extremal-Kerr strip {x ≥ 7/8, v ≤ 1/8}, to be handled by charts EA/EB/KA/KB (implemented).
  Weakly definite cores are now detected on the whole cube and treated as constants.
* Main cube: open vol ~8e-6, but stuck at level 36 on a surface near τ ≈ 1 and near (x, v) ≈ (0.985, 0.865).
  Cause (mpmath at (0.9858, 0.8643, 0.9763)): for Q2, F1 ≈ 1.6e11, Y ≈ 1.2e12, and F1√D − Y ≈ 2.2e7, a relative
  cancellation of ~1e-5. The value at b₊ is tiny compared with its two components (n(b₋) is huge), so the
  (F1, Y, P) rule needs boxes of size ~1e-5 per axis (level ~50). This is a defect of the b-reduction, not of S_α
  (the value is positive).
* NEXT: avoid the b-reduction near this surface. Options: (i) a parametrisation in which a² is polynomial (e.g. use
  (r0, a) or (λ̃, a) as coordinates and solve for r2 rationally, if possible); (ii) scale out the factor making n(b₊)
  small (find the polynomial relation responsible); (iii) interval evaluation of n at b₊ with a tight √D enclosure.

## 196. Rational parametrisation attempt for b = a² (to kill the b₊ cancellation): NEGATIVE so far
* For fixed r0, the relation b² − e2·b + e4 = 0 is a CUBIC curve in (b, r2). The rational point (b, r2) = (0, 0) and
  the lines r2 = s·b leave a residual quadratic in b with discriminant (s−1)(r0s−1)(r0s² + 3r0s + 3s + 1) (sympy),
  so this pencil does not rationalise it.
* With coordinates (r0, r1 = 1, a), λ̃ and M come out rational (Δ(r0) = Δ(1) = 0 is linear in them), but r2 then
  solves a quadratic. Some square root is unavoidable in any parametrisation by two horizons plus one physical
  parameter.
* Open options: a singular point of the cubic (if it is singular for every r0, the curve is rational and a
  parametrisation exists); or work at b₊ with a 4th variable z = √D and an interval enclosure. The cancellation
  needs an enclosure of z of relative width ~1e-5 near the stuck surface (τ ≈ 1).
* [P] UPDATE: for each fixed r0 the cubic b² − e2·b + e4 = 0 is NONSINGULAR, both affine and at infinity (gradients at
  (0:1:0) and (1:0:0) are nonzero). So it is a smooth plane cubic, of genus 1: no rational parametrisation exists.
  The b₊ square root cannot be removed by a change of coordinates; the cancellation must be handled with a √D
  enclosure.

## 197. S_α covering: per-factor b₊ signs remove the cancellation; edge charts; one degenerate point left
* FIX for item 195: instead of b-reducing the PRODUCT n·d (which creates the 1e-5 cancellation), factor n and d,
  keep the odd-multiplicity factors, and decide the sign of EACH b-dependent factor at b₊ separately with the
  (F1, Y, P) rule; b-free factors are signed directly. cert3/sal_build2.py, cover3d_g.py.
* Main cube: CLOSED (float engine). Excluded regions are each handled by blow-up charts:
  VT (v → 0, τ → 1; Kerr limit at the top of the window): both charts closed in 1 box.
  XT (x → 1, τ → 1; extremal at the top of the window): closes except near the edge point
     (x, τ) = (1, 1), v ≈ 0.86, where on the chart face σ = 0 the main numerator core (14, 16, 9, 3) has F1, Y, P all
     changing sign in v. This suggests the value of the Q2 (m²) coefficient of W_α vanishes in the limit at an isolated
     edge point v = v*.
  EA/EB (ultracold) and KA/KB (extremal-Kerr strip): not yet run for S_α.
* OPEN QUESTION before more covering: is the Q2 coefficient of W_α genuinely → 0 at (x, v, τ) → (1, v*, 1)? It is
  3t(t−1)b·ê²F̂λ̃ + 2.1·LQ2'·G', with both terms ≥ 0, so this requires ê → 0 AND LQ2' → 0 there. If genuine, near
  that point S_α relies on the lower-order m-terms, and a point blow-up (with algebraic v*) is needed.

## 198. S_α near the edge x = τ = 1: the degenerate point is an artifact; it persists in an A-free variant
* [N] Physically W_α/m³ at m = 1000 near (x, τ) → (1, 1), v ∈ [0.84, 0.87] is ~1e8 / 1e11 / 1e14 at ε = 1e-2 / 1e-3 /
  1e-4, dominated by 3t(t−1)(Im et)²F. There is NO genuine zero; the chart-face zero is an ORDER-MIXING artifact
  (the leading σ-orders of F1 and Y differ, and the dominant one's leading coefficient vanishes along a curve v*(γ),
  v* ≈ 0.8175 at γ = 0.5).
* A-free variant (salB): W_α ≥ 3t(t−1)(Im et)²F − 2.1(3t−B)·G. Physical margin ≥ 4.03 on the whole XT region
  (118k points). Built (salB_build*.py), but its XT-chart covering stalls in the same band v ≈ 0.83–0.87, which
  borders the ultracold square (v ≥ 7/8).
* So the XT charts need either the ultracold exclusion pushed down to v ≥ 0.8 (with EA/EB covering more), or a
  curve blow-up along v*(γ). Not done.
STATUS of hyp:cos: S_μ ✔, S_C ✔ (float caveat), S_α: main cube ✔, VT ✔, XT/EA/EB/KA/KB open; S_β not started.

## 199. S_α: ultracold region pushed to [3/4, 1]²; everything closes except the ultracold charts
* XT charts (x ≥ 7/8, τ ≥ 7/8) with v ∈ (1/8, 3/4): CLOSED (the v ≤ 1/8 part is covered by VT, the v ≥ 3/4 part by EA/EB).
* KA/KB (extremal-Kerr strip x ≥ 7/8, v ≤ 1/8): closed in 1 box each. VT: closed.
* EA/EB at scale 1/4 (the square [3/4, 1]²): open vol 0.030 / 0.044, stuck at small σ (the ultracold face), β ∈ [0.88, 1].
* Sufficient forms that would avoid the face degeneracy FAIL on [3/4, 1]²:
   - via S_μ and |u| ≤ F: 3(t−1)F² ≥ 8.4t(3t−B)(4t−2) has min 1.7e-3 (Nariai side, v → 1);
   - A-free W_α has min ~1e-10 (τ → 0).
  So the corner needs the exact W_α plus sub-blow-ups of the ultracold face, as for S_C (items 183–194).
STATUS hyp:cos: S_μ ✔; S_C ✔; S_α ✔ outside the ultracold square [3/4, 1]² (all float-engine); S_β not started.

## 200. S_α ultracold square: progress and the last stuck pocket
* [P, cert] The m² coefficient of W_α is positive on the WHOLE window, for free: f1 < 0, f2 < 0 and LQ2' > 0 are each
  certified on the full cube in ONE box (cert3/aux_build.py, cover3d_g NOEX). Then Q2 = 3t(t−1)b·ê²F̂λ̃ +
  2.1·LQ2'·G' > 0 (F̂ > 0, G' > 0 follow from f1, f2 < 0 and t > 1). So only H2 (m = 2) and H1 (m = 1) remain.
* OR-closure (cover3d_h.py, env ALT): a box closes if W_α's H2, H1 > 0 OR an alternative sufficient set holds:
  (a) L(t) > 0 at m = 1, 2 (smu2_res; then W_α ≥ 0 trivially), or (b) the A-free form (salB2_res).
  Chart EA (x − 1 dominant, scale 1/4): CLOSED in 1 box via (a).
  Chart EB: open vol 0.026 → 4.7e-4. The remaining pocket is EB σ ∈ [0.047, 0.30], β ∈ [0, 0.46], τ ≥ 0.68, i.e. r2 ≈
  r1 + (1–8)%, r0 much closer to 1 than r2, upper window. There the main numerator core of H1/H2 is ill-conditioned
  at b₊ (F1, Y, P all undetermined at level 27) and neither (a) nor (b) closes.
STATUS hyp:cos: S_μ ✔, S_C ✔, S_α ✔ except this pocket (float engine), S_β not started.

## 201. S_α pocket: a well-margined sufficient form exists, but certificates stay ill-conditioned there
* [N] In the pocket (r2 = 1/(1−σ/4), r0 = 1−σβ/4, σ ∈ [0.04, 0.31], β ∈ [0, 0.47], τ ≥ 0.67, m = 1, 2):
  G_α: 3(t−1)F² ≥ 8.4·t(3t−B)(4t−2) holds with min ratio 7.2 (125k points; 25% of them on the gate).
  G_α ⇒ S_α on the gate, using S_μ (t|Im et| ≥ F/2), |Im e1| ≤ F and A ≥ 0. Polys: gal_build.py, validated (gal_chk) to 1e-13.
* As an OR-alternative in chart EB it closes almost nothing (open vol 4.63e-4 → 4.58e-4).
* Diagnosis: at a pocket point the H2 main numerator core is small at b₊ (rel 5e-3), while the factors (r0 − 1) and
  (r0 − r2)-type are small (rel 4e-3, 2e-2). The core carries a HIDDEN power of the distance to r0 = 1 (the β = 0 face
  of chart EB, i.e. extremal r0 → r1) that is not stripped by the σ-blow-up. So the β = 0 face of EB is degenerate
  for τ ≳ 0.68. That needs a second blow-up (β, 1−τ) → 0 inside EB, or stripping the β-power by a weighted chart
  (x = 1 − σβ/4 with β = ρ^k).
* [N] Stripping monomial powers of β (and of τ) after the EB chart substitution changes nothing (open vol
  4.56e-4 → 4.54e-4). The degeneracy is not a monomial factor; most likely the face polynomial on {β = 0} (extremal
  r0 = r1) vanishes along a curve in (σ, τ) near τ ≈ 0.68–1. Next: compute that face curve explicitly (the restriction
  of the H2/H1 cores to β = 0) and blow up along it, or prove the pocket by hand via G_α (ratio ≥ 7.2) with the
  extremal-limit expansion.

## 202. The S_α pocket face curve located
On the EB-chart face β = 0 (extremal r0 = r1), the H2 main core's F1, Y, P all vanish TOGETHER along a curve
τ*(σ): none for τ ≥ 0.6 at σ = 0.05; τ* ≈ 0.83 at σ = 0.15; τ* ≈ 0.99 at σ = 0.3 (face11.py).
Interpretation: as r0 → r1, D'(r1) → 0, so Im e1 ~ (w − Ω1)/D'(r1) dominates unless w ≈ Ω1, and the curve is the
locus where the leading (w − Ω1)/β term cancels. So S_α at the pocket is a two-scale problem
(β → 0 with w − Ω1 = O(β)). Plan: blow up along w = Ω1 (i.e. τ = τ*(σ) + β·s) in the EB chart, or do the hand
expansion of G_α with the scaling w − Ω1 = β·s.
* [N] Physically W_α is smooth and hugely positive straight across the face curve: at r2 = 1/(1−0.0375), τ ∈ [0.78, 0.88],
  r0 = 1 − ε: W ~ 6–8e10 (ε = 1e-3, m = 1) and 6–8e16 (ε = 1e-5), with no dip near τ* ≈ 0.83. So the face-curve
  sign change of the core is a representation artifact of the order-mixing kind: the true leading behaviour of n near
  r0 = r1 is not captured by "strip the minimal β-power per monomial" (W ~ ε⁻², d ~ ε³). The right fix is to expand n
  in β exactly at b₊ (the b-reduction mixes orders), i.e. b-reduce AFTER the chart substitution and β-expansion,
  not before.
* CORRECTION to the line above: stripping positive monomials never changes the SIGNS of F1, Y, P, so "reduce after
  substituting" would not help by itself. The real picture on the face β = 0: across τ*(σ), (F1, Y, P) go from
  (+, +, +) to (−, −, −). The b₊-value is positive on BOTH sides (P > 0, F1 > 0 ⇒ +; P < 0, Y < 0 ⇒ +), and it
  vanishes to leading face order ON the curve. So boxes straddling the curve can never be closed by the sign rule. What
  is needed is the next-order term in β along the curve, i.e. a genuine blow-up τ = τ*(σ) + β·s, with τ*(σ) algebraic
  (implicit in the face polynomial). Not attempted.

## 203. S_α pocket: a gate-only shortcut via S_C (math OK), covering still stalls
* [P] S_α ⇐ S_C ∧ (2 + 4a + T ≤ 3T) ⇔ S_C ∧ (t−1)(1 − 3zf) ≥ 1. In window terms (u, et < 0):
  K := (t−1)(2f1 − t·f2) + f1 + t·f2 > 0 (positive factor 2Ξ/λ̃ removed; m- and a-free). kal_build.py.
* [N] On the pocket's gate points (44k, m = 2): (t−1)(1 − 3zf) ≥ 2.1. Off the gate K can be very negative (fine:
  S_α is only needed on the gate).
* Per-m OR-closure (cover3d_h.py, PERM=1): a box closes if K > 0, OR [for m = 1, 2: (W_α(m) > 0 or L(t, m) > 0)];
  Q2 > 0 holds globally (item 200).
* Result: EB open vol 4.48e-4 → 4.40e-4; ~1k boxes per run. The pocket boxes still touch the β = 0 face (r0 = r1,
  where t → ∞); presumably K's cleared form degenerates there too.
REALISTIC STATUS: the remaining S_α pocket needs either the curve blow-up (item 202) or a direct near-extremal
hand argument (t → ∞ as r0 → r1, so (t−1)(1 − 3zf) → ∞ unless zf → 1/3; the gate forces zf ≤ 1/2 only).
* [N, NEGATIVE] The K-route (item 203) FAILS near extremality: on gate points with r0 = 1−ε, v ∈ [0.75, 0.995]:
  max zf = 0.319 / 0.339 / 0.346 for ε ∈ [1e-2, 1e-1] / [1e-3, 1e-2] / [1e-5, 1e-3], and (t−1)(1−3zf) reaches −925.
  So the gate does NOT keep zf below 1/3 near r0 = r1; K is useless there. The pocket genuinely needs W_α (or G_α) at
  the two-scale limit (w − Ω1 = O(β)), i.e. the curve blow-up.

## 204. S_β: first look
S_β (item 175): A ≥ Q2·sup_{0<s<T} s·r[(2az+s)(2+r(T−s)) − r s z(z−1)]/(z−1)², with z = z0+s, a = z0−1, T = t−z0,
r = κ/q1. Since Q2·r = K := 1.05(−L(t))t(t−1)/(2T²): S_β ⇔ A ≥ K·Φ, Φ := sup_s φ(s),
φ(s) = s[(2az+s)(2+r(T−s)) − r s z c]/c², c = z−1 = a+s.
[N] 1398 gate points: min A/(KΦ) = 1.338. The argmax over s is interior and widely spread (median s/T ≈ 0.19,
10% at 0.018, 10% at s = T).
[N, NEGATIVE] Closed-form bounds of Φ all fail badly: using 2az + s ≤ 2zc gives min ratio 5.8e-4; t(4+2rT) gives
4.3e-4; 6t gives 3.5e-4. The cancellation between (2az+s)(2+r(T−s)) and r s z c is essential (the loss is ~t when
a ≪ s).
PLAN: S_β ⇔ the QUARTIC in s, Ψ(s) := A·c² − K·s[(2az+s)(2+r(T−s)) − r s z c] ≥ 0 on [0, T] (on the gate). Certify
it in 4D (x, v, τ, s/T) with OR-closure "L(t) > 0" for off-gate boxes. This needs a 5th variable in the F-class /
b-reduction pipeline (currently fixed at (r0, r2, b, τ)). Alternatively, eliminate s through quartic-nonnegativity
conditions on an interval (Sturm/discriminant), which is heavy.
* m-structure of S_β (worked out): with q := Q2/m², α := A/m² (both m-free) and κ̂ := K/m² = c_K(k2 + k1u + k0u²)
  (u = 1/m, k0 < 0, so concave in u),
    Q2·Ψ/m⁴ = g(κ̂) := q·α·c² − 2q·κ̂·s(2az+s) − κ̂²·s·E,   E := (2az+s)(T−s) − s z c.
  S_β is g(κ̂) ≥ 0 for all s ∈ [0, T] and all κ̂ ∈ (0, κ̂_max(geometry)] that occur (gate ⇔ κ̂ > 0). Since g is a
  concave/convex quadratic in κ̂ depending on the sign of E, the natural certificate is in (x, v, τ, s/T, κ̂/κ̂_max): 5D,
  with a b-reduction. This is the largest certificate of the four; not started.

## 205. S_β reduced from 5D to two 4D statements (m = 1, 2)
[P] Put ρ := rT = κT/q1 (so ρ ≤ 1 is S_C, proved) and Γ := A·T/q1², which is m-FREE (A, q1² ∝ m²). Writing
s = Tσ, φ(s) = φ0(s) + ρ·φ1(s) is affine in ρ, with
  φ0 = 2s(2az+s)/c²,  φ1 = s[(2az+s)(1−σ) − σ z c]/c²   (c = z−1, z = z0+s, a = z0−1).
S_β ⇔ Γ ≥ ρ·sup_s φ = sup_s (ρφ0 + ρ²φ1) =: H(ρ), and H is nondecreasing in ρ ≥ 0 [N: 0 violations on 1406 points].
By the S_C identity, ρ_m = 2.1·(−λ̃L/m²)/(λ̃·ê(û + tê)) with an m-free denominator, while
λ̃L/m² = LQ2' + LQ1·u + LQ0'·u² (u = 1/m) is CONCAVE in u (LQ0' < 0). So ρ_m is convex in u ∈ [0, 1/2], and
ρ_m ≤ max(ρ_2, ρ_∞ < 0) for all m ≥ 2. Hence
   S_β for ALL m  ⇐  Γ ≥ H(ρ̄),  ρ̄ := max(ρ_1, ρ_2)   ⇐   [Γ ≥ H(ρ_1)] ∧ [Γ ≥ H(ρ_2)].
Each bracket is a 4D statement in (x, v, τ, σ): Γ·c² ≥ ρφ0c² + ρ²φ1c² for σ ∈ [0, 1] (quartic in s), with an
OR-closure "ρ_m ≤ 0" (that m is off the gate).
[N] Exact margin min 1.287 (1406 gate points). The crude ρ := 1 substitute FAILS (2.3e-4), so ρ_m must be kept.
NEXT: extend the polynomial pipeline to a 5th symbol s (σ), then run the 4D float covering for m = 1 and m = 2.
