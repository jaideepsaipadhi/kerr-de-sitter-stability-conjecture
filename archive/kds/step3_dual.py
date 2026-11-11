"""
Step 3: KdS radial Teukolsky -> Heun form -> Euler-transform dual -> normal-form potential.

Question: for real omega and real lambda, is there a gauge choice + Euler transform
(A or B) whose dual normal-form potential Qt(z) is REAL on (1, t), i.e. on the exterior
r_e < r < r_c?  (CTdC / HIZ say the m2<->m3 mass symmetry achieves this.)

Conventions (CTdC): roots r3 < r0 < r1 < r2 of Delta; r1 = event horizon, r2 = cosmological.
Mobius map z = zinf (r - r0)/(r - r3):  r0->0, r1->1, r2->t, r3->infinity.
"""
import itertools
import mpmath as mp
import sympy as sp

mp.mp.dps = 30


def kds_roots(Lam, a, M):
    # Delta = -(Lam/3) r^4 + (1 - Lam a^2/3) r^2 - 2 M r + a^2
    co = [-Lam/3, 0, 1 - Lam*a**2/3, -2*M, a**2]
    rts = sorted([mp.re(x) for x in mp.polyroots([mp.mpf(c) for c in co], maxsteps=200, extraprec=60)])
    return rts  # r3 < r0 < r1 < r2


def teukolsky_heun(Lam, a, M, s, omega, m, lam):
    """Return dict with t, Heun params (g,d,e,A,B,q) for each of the 8 gauge choices,
    plus exponents.  Done symbolically in z with numeric coefficients."""
    Lam, a, M = mp.mpf(Lam), mp.mpf(a), mp.mpf(M)
    r3, r0, r1, r2 = kds_roots(Lam, a, M)
    zinf = (r1 - r3)/(r1 - r0)
    t = zinf*(r2 - r0)/(r2 - r3)
    alpha = Lam*a**2/3
    Xi = 1 + alpha

    z = sp.symbols('z')
    rr = sp.symbols('r')
    F = lambda x: sp.Float(str(x), 35)
    D = (rr**2 + F(a)**2)*(1 - F(Lam)*rr**2/3) - 2*F(M)*rr
    Kq = F(omega.real)*(rr**2 + F(a)**2) - F(a)*m
    if omega.imag != 0:
        Kq = (F(omega.real) + sp.I*F(omega.imag))*(rr**2 + F(a)**2) - F(a)*m
    om = F(omega.real) + sp.I*F(omega.imag)
    XiF = F(Xi)
    Q = (XiF**2*Kq**2 - sp.I*s*XiF*Kq*sp.diff(D, rr))/D + 4*sp.I*s*XiF*om*rr \
        - sp.Rational(2, 3)*F(Lam)*(s + 1)*(2*s + 1)*rr**2 + s*(1 - F(alpha)) - F(lam)
    # ODE: R'' + p R' + qq R = 0 in r
    p_r = (s + 1)*sp.diff(D, rr)/D
    q_r = Q/D

    # change variables r -> z
    r_of_z = (F(r0)*F(zinf) - F(r3)*z)/(F(zinf) - z)
    drdz = sp.diff(r_of_z, z)
    d2rdz2 = sp.diff(drdz, z)
    # R_r = R_z / r'(z);  R_rr = (R_zz - R_z r''/r') / r'^2
    # R_zz + [p r' - r''/r'] R_z + q r'^2 R = 0
    Pz = sp.together(p_r.subs(rr, r_of_z)*drdz - d2rdz2/drdz)
    Qz = sp.together(q_r.subs(rr, r_of_z)*drdz**2)

    # exponents at r_j (kappa_j = i Xi K_j / Delta'_j):  {kappa, -kappa - s}
    Dp = sp.diff(D, rr)
    exps = {}
    for name, rj in [('0', r0), ('1', r1), ('t', r2), ('inf', r3)]:
        Kj = complex(sp.N(Kq.subs(rr, F(rj)), 30))
        Dpj = complex(sp.N(Dp.subs(rr, F(rj)), 30))
        kap = 1j*complex(Xi)*Kj/Dpj
        exps[name] = (kap, -kap - s)
    # The point z = zinf (r = infinity) is removable: exponents {2s+1, 2s+2} in (z - zinf).
    return dict(t=t, zinf=zinf, Pz=Pz, Qz=Qz, z=z, exps=exps, roots=(r3, r0, r1, r2))


def gauge_to_heun(data, choice):
    """choice = (i0,i1,it) selects exponent index at 0,1,t to factor out.
    R = z^{e0} (z-1)^{e1} (z-t)^{et} (z-zinf)^{2s+1} f.  Return Heun params of f."""
    z, t, zinf = data['z'], data['t'], data['zinf']
    e0 = data['exps']['0'][choice[0]]
    e1 = data['exps']['1'][choice[1]]
    et = data['exps']['t'][choice[2]]
    return e0, e1, et


def heun_params_from_operator(Pz, Qz, z, t, e0, e1, et, ez, s):
    """Substitute R = phi f, phi = z^e0 (z-1)^e1 (z-t)^et (z-zinf)^ez; read off Heun params
    numerically and check the result really is Heun."""
    tt = sp.Float(str(t), 35)
    E0, E1, Et = [sp.Float(str(x.real), 30) + sp.I*sp.Float(str(x.imag), 30) for x in (e0, e1, et)]
    # log-derivative of phi
    L1 = E0/z + E1/(z - 1) + Et/(z - tt) + ez['e']/(z - ez['zinf'])
    L2 = sp.diff(L1, z) + L1**2  # phi''/phi
    # f'' + (P + 2 L1) f' + (Q + P L1 + L2) f = 0
    Pf = sp.together(Pz + 2*L1)
    Qf = sp.together(Qz + Pz*L1 + L2)
    g = 1 + 2*E0 + complex(0)
    # Heun: P = g/z + d/(z-1) + e/(z-t);  Q = (AB z - q)/(z(z-1)(z-t))
    # Read g,d,e as residues of Pf; check no residue at zinf
    def res(expr, pt):
        return complex(sp.N(sp.limit((z - pt)*expr, z, pt) if False else ((z - pt)*expr).subs(z, pt + sp.Float('1e-25')), 25))
    def numres(expr, pt, h=1e-12):
        f_ = sp.lambdify(z, expr*(z - pt), 'mpmath')
        return complex(f_(mp.mpf(pt) + mp.mpf(h)))
    gres = numres(Pf, 0)
    dres = numres(Pf, 1)
    eres = numres(Pf, t)
    zres = numres(Pf, ez['zinf'])
    # Q check: (Q * z(z-1)(z-t)) should be linear in z
    Qf_l = sp.lambdify(z, Qf*z*(z - 1)*(z - tt), 'mpmath')
    pts = [mp.mpf('0.37'), mp.mpf('1.9'), mp.mpf('2.6'), mp.mpf('-3.1')]
    vals = [complex(Qf_l(p)) for p in pts]
    slope = (vals[1] - vals[0])/complex(pts[1] - pts[0])
    icpt = vals[0] - slope*complex(pts[0])
    lin_err = max(abs(vals[k] - (slope*complex(pts[k]) + icpt)) for k in range(4))
    return dict(g=gres, d=dres, e=eres, AB=slope, q=-icpt, zres=zres, lin_err=lin_err)


def euler_dual(h, t, which):
    g, d, e, AB, q = h['g'], h['d'], h['e'], h['AB'], h['q']
    # A, B from A+B = g+d+e-1, A B = AB
    S = g + d + e - 1
    disc = mp.sqrt(mp.mpc(S)**2 - 4*mp.mpc(AB))
    A_ = (S + disc)/2 if which == 0 else (S - disc)/2
    A_ = complex(A_)
    t = complex(t)
    gt = -A_ + g + 1
    dt = -A_ + d + 1
    et = -A_ + e + 1
    ABt = (A_ - 2)*(2*A_ - d - e - g)
    qt = A_**2*t + A_**2 - A_*d*t - A_*e - A_*g*t - A_*g - A_*t - A_ + d*t + e + g*t + g + q
    return dict(g=gt, d=dt, e=et, AB=ABt, q=qt, A=A_)


def normal_form_Q(h, t, zval):
    g, d, e, AB, q = h['g'], h['d'], h['e'], h['AB'], h['q']
    P = g/zval + d/(zval - 1) + e/(zval - t)
    Pp = -g/zval**2 - d/(zval - 1)**2 - e/(zval - t)**2
    Qh = (AB*zval - q)/(zval*(zval - 1)*(zval - t))
    return Qh - Pp/2 - P**2/4


def real_branch(h, t):
    """Return the Euler branch (0 = A, 1 = B) whose dual normal-form potential is real on (1,t)."""
    best = None
    for which in (0, 1):
        ht = euler_dual(h, complex(t), which)
        zs = [1 + (t.real - 1)*k/9 for k in range(1, 9)]
        im = max(abs(normal_form_Q(ht, complex(t), complex(zz)).imag) for zz in zs)
        sc = max(abs(normal_form_Q(ht, complex(t), complex(zz)).real) for zz in zs)
        if best is None or im/sc < best[0]:
            best = (im/sc, which, ht)
    assert best[0] < 1e-8, f'no real branch (rel Im {best[0]:.1e})'
    return best[1], best[2]


def scan(Lam, a, M, s, omega, m, lam):
    data = teukolsky_heun(Lam, a, M, s, omega, m, lam)
    t = float(data['t'])
    zinf = data['zinf']
    out = []
    for choice in itertools.product([0, 1], repeat=3):
        e0, e1, et = gauge_to_heun(data, choice)
        h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'], e0, e1, et,
                                      {'e': 2*s + 1, 'zinf': sp.Float(str(zinf), 35)}, s)
        # original normal form on (1,t)
        zs = [1 + (t - 1)*k/10 for k in range(1, 10)]
        imQ_orig = max(abs(normal_form_Q(h, complex(t), complex(zz)).imag) for zz in zs)
        for which in (0, 1):
            ht = euler_dual(h, complex(t), which)
            Qs = [normal_form_Q(ht, complex(t), complex(zz)) for zz in zs]
            imQ = max(abs(x.imag) for x in Qs)
            reQ = max(abs(x.real) for x in Qs)
            out.append((choice, which, imQ, reQ, imQ_orig, h['lin_err'], abs(h['zres'])))
    return out, data


if __name__ == '__main__':
    import sys
    # A generic subextremal KdS black hole, superradiant real frequency
    Lam, a, M = 0.05, 0.6, 1.0
    s, m, lam = -2, 2, 7.3
    r3, r0, r1, r2 = kds_roots(Lam, a, M)
    print('roots r3,r0,r1,r2 =', [float(x) for x in (r3, r0, r1, r2)])
    Xi = 1 + Lam*a**2/3
    Om1 = a/(r1**2 + a**2); Om2 = a/(r2**2 + a**2)
    print('horizon angular velocities (event, cosmo):', float(Om1), float(Om2))
    omega = complex(m*float((Om1 + Om2)/2))  # inside superradiant band m*Om2 < omega < m*Om1
    print('omega =', omega)
    res, data = scan(Lam, a, M, s, omega, m, lam)
    print(f"{'gauge':>10} {'A/B':>4} {'max|Im Qdual|':>15} {'max|Re Qdual|':>15} {'max|Im Qorig|':>15} {'Heun lin err':>13} {'res@zinf':>10}")
    for choice, which, imQ, reQ, imQo, le, zr in res:
        print(f'{str(choice):>10} {"AB"[which]:>4} {imQ:15.3e} {reQ:15.3e} {imQo:15.3e} {le:13.2e} {zr:10.2e}')
