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
