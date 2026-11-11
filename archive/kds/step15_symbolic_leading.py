"""
Step 15: closed-form leading-order (m -> infinity) dual potential P2(z).

Leading order of the original equation (spin drops out):
    Q_r2 = Xi^2 Khat^2/Delta^2 - ell/Delta,   Khat = w (r^2 + a^2) - a.
In z (Mobius map), Q_z2 = Q_r2 (dr/dz)^2 at leading order.
Original Heun data at leading order (gauge (0,0,0), i.e. exponent kappa_j at each finite point):
    g1 = 2 k0/m, d1 = 2 k1/m, e1 = 2 kt/m,  k_j = i Xi Khat(r_j) m / Delta'(r_j)
    AB2 z - q2 = z (z-1)(z-t) [Q_z2 + P1^2/4],  P1 = g1/z + d1/(z-1) + e1/(z-t)
(the right-hand side is checked to be linear in z).
Dual, leading order:
    A1 solves A^2 - S1 A + AB2 = 0 with S1 = g1 + d1 + e1 (branch fixed by matching numerics)
    gt1 = g1 - A1, dt1 = d1 - A1, et1 = e1 - A1, ABt2 = A1 (2 A1 - S1)
    qt2 = q2 + A1^2 (t+1) - A1 (d1 t + e1 + g1 t + g1)
    Qt2 = (ABt2 z - qt2)/(z (z-1)(z-t)) - Pt1^2/4,   P2 = mu^2 Qt2,   mu = z (z-1)(t-z)
Validation: compare with the m-extrapolated p4 from step 13.
"""
import numpy as np
import sympy as sp
from step3_dual import kds_roots
from step5_window import ctdc_window
from step13_largem import leading

z, r = sp.symbols('z r')


def P2_closed_form(Lam, a, M, w, ell, branch=0):
    import mpmath as mp
    mp.mp.dps = 40
    r3, r0, r1, r2 = [mp.mpf(x) for x in kds_roots(Lam, a, M)]
    Lam, a, M, w, ell = [mp.mpf(x) for x in (Lam, a, M, w, ell)]
    zinf = (r1 - r3)/(r1 - r0); t = zinf*(r2 - r0)/(r2 - r3)
    Xi = 1 + Lam*a**2/3
    D = lambda r: (r**2 + a**2)*(1 - Lam*r**2/3) - 2*M*r
    Dp = lambda r: 2*r*(1 - Lam*r**2/3) - (r**2 + a**2)*2*Lam*r/3 - 2*M
    Kh = lambda r: w*(r**2 + a**2) - a
    Qr2 = lambda r: Xi**2*Kh(r)**2/D(r)**2 - ell/D(r)
    rz = lambda zz: (r0*zinf - r3*zz)/(zinf - zz)
    drdz = lambda zz: (r0*zinf - r3*zinf)/(zinf - zz)**2
    Qz2 = lambda zz: Qr2(rz(zz))*drdz(zz)**2
    kap = [1j*Xi*Kh(rj)/Dp(rj) for rj in (r0, r1, r2)]
    g1, d1, e1 = [2*k for k in kap]
    P1 = lambda zz: g1/zz + d1/(zz - 1) + e1/(zz - t)
    F = lambda zz: zz*(zz - 1)*(zz - t)*(Qz2(zz) + P1(zz)**2/4)
    pts = [mp.mpf(x) for x in ('0.31', '1.7', '-2.3', '4.9', '0.77')]
    vals = [F(p) for p in pts]
    AB2 = (vals[1] - vals[0])/(pts[1] - pts[0]); c0 = vals[0] - AB2*pts[0]
    lin_err = max(abs(vals[k] - (AB2*pts[k] + c0)) for k in range(5))/max(abs(v) for v in vals)
    q2 = -c0
    S1 = g1 + d1 + e1
    disc = mp.sqrt(S1**2 - 4*AB2)
    A1 = (S1 + disc)/2 if branch == 0 else (S1 - disc)/2
    gt1, dt1, et1 = g1 - A1, d1 - A1, e1 - A1
    ABt2 = A1*(2*A1 - S1)
    qt2 = q2 + A1**2*(t + 1) - A1*(d1*t + e1 + g1*t + g1)
    Pt1 = lambda zz: gt1/zz + dt1/(zz - 1) + et1/(zz - t)
    Qt2 = lambda zz: (ABt2*zz - qt2)/(zz*(zz - 1)*(zz - t)) - Pt1(zz)**2/4
    mu = lambda zz: zz*(zz - 1)*(t - zz)
    zs = [1 + (t - 1)*mp.mpf(k)/8 for k in range(1, 8)]
    P2v = [mu(zz)**2*Qt2(zz) for zz in zs]
    coef = mp.matrix([[zz**j for j in range(5)] for zz in zs[:5]])**-1*mp.matrix(P2v[:5])
    quart_err = max(abs(sum(coef[j]*zs[k]**j for j in range(5)) - P2v[k]) for k in (5, 6))/max(abs(v) for v in P2v)
    return float(t), [complex(coef[j]) for j in range(5)], float(max(lin_err, quart_err)), dict(A1=complex(A1))


if __name__ == '__main__':
    for (Lam, a, wf, ell) in [(0.05, 0.6, 0.5, 3.0), (0.001, 0.95, 0.9, 1.0), (0.05, 0.9877, 0.2, 10.0), (0.02, 0.9, 0.98, 1.5)]:
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
        w = lo + wf*(hi - lo)
        ref = leading(Lam, a, wf, ell)
        for br in (0, 1):
            t, P2, rem, info = P2_closed_form(Lam, a, 1.0, w, ell, branch=br)
            P2r = np.array(P2)
            err = np.max(np.abs(P2r[:len(ref['p4'])] - ref['p4']))/np.max(np.abs(ref['p4']))
            print(f'Lam={Lam} a={a} wfrac={wf} ell={ell} branch={br}: max|Im P2|={np.max(np.abs(P2r.imag)):.1e}  '
                  f'rel. difference from extrapolated numerics={err:.2e}   (division remainder {rem:.1e})')


def P2_explicit(Lam, a, M, w, ell):
    """Explicit leading-order dual with A1 = -2 kappa_3:
         gt1 = 2(k0+k3), dt1 = 2(k1+k3), et1 = 2(kt+k3),  ABt2 = 4 k3^2,
         qt2 = q2 + 4 k3^2 (t+1) + 2 k3 (d1 t + e1 + g1 t + g1)."""
    import mpmath as mp
    t, P2ref, err, info = P2_closed_form(Lam, a, M, w, ell, branch=1)
    r3, r0, r1, r2 = [float(x) for x in kds_roots(Lam, a, M)]
    Xi = 1 + Lam*a**2/3
    Dp = lambda r: 2*r*(1 - Lam*r**2/3) - (r**2 + a**2)*2*Lam*r/3 - 2*M
    Kh = lambda r: w*(r**2 + a**2) - a
    k0, k1, kt, k3 = [1j*Xi*Kh(rj)/Dp(rj) for rj in (r0, r1, r2, r3)]
    return abs(info['A1'] + 2*k3)/abs(k3), (k0, k1, kt, k3), t


if __name__ == '__main__':
    print('\nCheck A1 = -2 kappa_3, and the kappa sum rule:')
    for (Lam, a, wf, ell) in [(0.05, 0.6, 0.5, 3.0), (0.001, 0.95, 0.9, 1.0), (0.05, 0.9877, 0.2, 10.0), (0.01, 0.99, 0.99, 100.0)]:
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
        rel, ks, t = P2_explicit(Lam, a, 1.0, lo + wf*(hi - lo), ell)
        print(f'  Lam={Lam} a={a}: |A1 + 2k3|/|k3| = {rel:.1e},  |sum k_j| = {abs(sum(ks)):.1e}')
