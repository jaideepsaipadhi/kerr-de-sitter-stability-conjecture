"""
Step 5: is the dual still superradiant, and does that match the CTdC window?

In the dual (0,0,0)A, the image of a QNM behaves like (z-1)^{A-delta} at z = 1 and is
analytic at z = t (in the Heun gauge).  In Hintz normalisation R = mu^{-1/2} u this gives
  qt(1) = i mu'(1) (A - delta)/2,   qt(t) = i |mu'(t)| (epsilon - A)/2,
with Hintz's outgoing convention (exponent -i q/mu' at z_e and -i q/|mu'| at z_c).

Flux (Wronskian) identity for a real potential: the dual is "superradiant" iff
qt(1) qt(t) < 0.

Check against CTdC (3.37), using
kappa_j = |Delta'(r_j)| / (2 Xi (r_j^2+a^2)) and omega_j = a/(r_j^2+a^2), with
index 0 = inner horizon, 1 = event horizon, 2 = cosmological horizon:
  (w1/k1 - w0/k0)/(1/k1 - 1/k0) < omega/m < (w2/k2 + w0/k0)/(1/k2 + 1/k0).
"""
import mpmath as mp
import sympy as sp
from step3_dual import kds_roots, teukolsky_heun, heun_params_from_operator, euler_dual, real_branch

mp.mp.dps = 30


def dual_q_signs(Lam, a, M, s, omega, m, lam):
    data = teukolsky_heun(Lam, a, M, s, complex(omega), m, lam)
    t = complex(data['t']).real
    e = data['exps']
    h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'],
                                  e['0'][0], e['1'][0], e['t'][0],
                                  {'e': 2*s + 1, 'zinf': sp.Float(str(data['zinf']), 35)}, s)
    _, ht = real_branch(h, complex(t))
    A = ht['A']
    mup1 = -3 + 2*(1 + t) - t          # mu'(1) = t - 1 > 0
    mupt = -3*t**2 + 2*(1 + t)*t - t   # mu'(t) = t - t^2 < 0
    q1 = 1j*mup1*(A - h['d'])/2
    qt = 1j*abs(mupt)*(h['e'] - A)/2
    return q1, qt, t


def ctdc_window(Lam, a, M):
    r3, r0, r1, r2 = kds_roots(Lam, a, M)
    Xi = 1 + Lam*a**2/3
    Dp = lambda r: -(4*Lam/3)*r**3 + 2*(1 - Lam*a**2/3)*r - 2*M
    k = lambda r: abs(Dp(r))/(2*Xi*(r**2 + a**2))
    om = lambda r: a/(r**2 + a**2)
    lo = (om(r1)/k(r1) - om(r0)/k(r0))/(1/k(r1) - 1/k(r0))
    hi = (om(r2)/k(r2) + om(r0)/k(r0))/(1/k(r2) + 1/k(r0))
    return float(lo), float(hi), float(om(r2)), float(om(r1))


if __name__ == '__main__':
    import numpy as np
    for (Lam, a, M) in [(0.05, 0.6, 1.0), (0.02, 0.9, 1.0), (0.1, 0.3, 1.0)]:
        s, m, lam = -2, 2, 7.3
        lo, hi, Om2, Om1 = ctdc_window(Lam, a, M)
        print(f'\nLam={Lam}, a={a}, M={M}: original superradiant band omega/m in ({Om2:.5f}, {Om1:.5f})')
        print(f'   CTdC window (3.37): ({lo:.5f}, {hi:.5f})')
        grid = np.linspace(Om2*0.5, Om1*1.3, 41)
        flips = []
        prev = None
        for x in grid:
            q1, qt, t = dual_q_signs(Lam, a, M, s, m*x, m, lam)
            assert abs(q1.imag) < 1e-8 and abs(qt.imag) < 1e-8, (q1, qt)
            sr = q1.real*qt.real < 0
            if prev is not None and sr != prev:
                flips.append(x)
            prev = sr
            print(f'   omega/m={x:.5f}  qt(1)={q1.real:+.4f}  qt(t)={qt.real:+.4f}  dual superradiant: {sr}')
        print('   sign changes of dual-superradiance near omega/m =', [round(f, 5) for f in flips])
