"""
Step 38: derive n1 (hence c1) in closed form, with no fitting.

Chain:
 1. Teukolsky normal form in r:  Q_r = q_r - p_r'/2 - p_r^2/4, with p_r = (s+1)Delta'/Delta, q_r = Q/Delta.
    Laurent coefficients at r3:  Q_r = a_{-2}/(r-r3)^2 + a_{-1}/(r-r3) + ...
 2. Mobius map: r - r3 = k/(zinf - z), k = (r0 - r3) zinf, and the Schwarzian vanishes, so
    Q_T(z) = Q_r(r(z)) r'(z)^2 and Q_T = c2/z^2 + c3/z^3 + ... with c2 = a_{-2}, c3 = 2 zinf a_{-2} - k a_{-1}.
 3. The Heun normal form equals Q_T, so the 1/z^3 coefficients give  q = AB(1+t) + T(1 - S/2) - c3,
    where S = g + d + e and T = d + t e.
 4. Dual parameters (step 2): gt = g - A + 1, dt = d - A + 1, et = e - A + 1, ABt = (A-2)(2A - d - e - g),
    qt = A^2 t + A^2 - A d t - A e - A g t - A g - A t - A + d t + e + g t + g + q.
 5. n1 = ABt(1+t) - qt + Tt(1 - St/2) + (7/2)(1+t),  n2 = ABt + St/2 - St^2/4 + 3/4  (should be -3).
Compare n1 with m^2 ell C0 + 1 + t + 3(1-alpha) C0 (ell = lambda/m^2) at 40-digit precision.
"""
import mpmath as mp
import sympy as sp
from step3_dual import kds_roots

mp.mp.dps = 50
r = sp.symbols('r')


def chain(Lam, a, M, s, omega, m, lam):
    r3, r0, r1, r2 = [mp.mpf(x) for x in kds_roots(Lam, a, M)]
    Lam_, a_, M_ = [mp.mpf(str(x)) for x in (Lam, a, M)]
    om = mp.mpc(omega); lam_ = mp.mpf(str(lam))
    al = Lam_*a_**2/3; Xi = 1 + al
    zinf = (r1 - r3)/(r1 - r0); t = zinf*(r2 - r0)/(r2 - r3)
    # exact exponents
    D = lambda x: (x**2 + a_**2)*(1 - Lam_*x**2/3) - 2*M_*x
    Dp = lambda x: 2*x*(1 - Lam_*x**2/3) - (x**2 + a_**2)*2*Lam_*x/3 - 2*M_
    K = lambda x: om*(x**2 + a_**2) - a_*m
    kap = [1j*Xi*K(x)/Dp(x) for x in (r0, r1, r2, r3)]
    g, d, e = [1 + 2*kap[j] + s for j in range(3)]
    A = 1 + s - 2*kap[3]; Binf = 1 + 2*s
    AB = A*Binf
    # Laurent coefficients of Q_r at r3, computed with sympy at high precision
    F = lambda v: sp.Float(mp.nstr(v, 45), 45) if not isinstance(v, mp.mpc) else (sp.Float(mp.nstr(v.real, 45), 45) + sp.I*sp.Float(mp.nstr(v.imag, 45), 45))
    eps = sp.symbols('eps')
    R3 = F(r3)
    Ds = -(F(Lam_)/3)*(r - R3)*(r - F(r0))*(r - F(r1))*(r - F(r2))
    Dps = sp.diff(Ds, r)
    Ks = F(om)*(r**2 + F(a_)**2) - F(a_)*m
    Qs = (F(Xi)**2*Ks**2 - sp.I*s*F(Xi)*Ks*Dps)/Ds + 4*sp.I*s*F(Xi)*F(om)*r - sp.Rational(2, 3)*F(Lam_)*(s + 1)*(2*s + 1)*r**2 + s*(1 - F(al)) - F(lam_)
    pr = (s + 1)*Dps/Ds; qr = Qs/Ds
    Qr = qr - sp.diff(pr, r)/2 - pr**2/4
    ser = sp.series(sp.together(Qr.subs(r, R3 + eps)), eps, 0, 0).removeO()
    am2 = complex(sp.N(ser.coeff(eps, -2), 40)); am1 = complex(sp.N(ser.coeff(eps, -1), 40))
    am2, am1 = mp.mpc(am2), mp.mpc(am1)
    k = (r0 - r3)*zinf
    c3 = 2*zinf*am2 - k*am1
    S = g + d + e; T = d + t*e
    q = AB*(1 + t) + T*(1 - S/2) - c3
    # dual
    gt, dt, et = g - A + 1, d - A + 1, e - A + 1
    ABt = (A - 2)*(2*A - d - e - g)
    qt = A**2*t + A**2 - A*d*t - A*e - A*g*t - A*g - A*t - A + d*t + e + g*t + g + q
    St = gt + dt + et; Tt = dt + t*et
    n1 = ABt*(1 + t) - qt + Tt*(1 - St/2) + mp.mpf(7)/2*(1 + t)
    n2 = ABt + St/2 - St**2/4 + mp.mpf(3)/4
    C0 = 3/(Lam_*(r1 - r0)*(r2 - r3))
    pred = lam_*C0 + 1 + t + 3*(1 - al)*C0
    return n1, n2, pred, am2


if __name__ == '__main__':
    import random
    random.seed(3)
    for (Lam, a, s, omega, m, lam) in [(0.05, 0.6, -2, 0.1, 2, 7.3), (0.001, 0.95, -2, 0.31, 3, 20.0), (0.05, 0.6, 2, 0.1, 2, 7.3),
                                       (0.02, 0.9, -2, 0.07, 1, 2.5), (0.0003, 0.4, 2, 0.02, 5, 40.0)]:
        n1, n2, pred, am2 = chain(Lam, a, 1.0, s, omega, m, lam)
        print(f'Lam={Lam} a={a} s={s:+d}: |n1 - (lam C0 + 1 + t + 3(1-alpha) C0)| / |n1| = {float(abs(n1 - pred)/abs(n1)):.1e},'
              f'  n2 = {complex(n2):.12f},  a_-2 = {complex(am2):.6f}')
