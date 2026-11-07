"""
Step 4: write each real dual in Hintz form and inspect the structure that matters.

Dual normal form: u'' + Qt(z) u = 0.  Put R = mu^{-1/2} u with mu(z) = z(z-1)(t-z) > 0 on (1,t).
Then (mu R')' + W R = 0, where W = mu Qt - (mu'^2 - 2 mu mu'')/(4 mu).
P4 := mu W is a real polynomial of degree <= 4 (checked numerically).

Hintz form:  -(mu R')' - (qt^2/mu) R + Lt(z) R = 0, i.e. P4 = qt^2 - Lt mu,
with qt a real quadratic and Lt(z) linear.  The interpolation conditions are
qt(z_j)^2 = P4(z_j) at z_j = 0, 1, t, so P4(z_j) >= 0 is needed (oscillatory exponents).

What the script reports for each dual:
  - P4(0), P4(1), P4(t)
  - for each sign choice of qt(z_j) = +-sqrt(P4(z_j)): whether qt changes sign on (1,t)
    (that is superradiance in the dual), and the coefficients of Lt
"""
import itertools
import mpmath as mp
import sympy as sp
from step3_dual import kds_roots, teukolsky_heun, gauge_to_heun, heun_params_from_operator, euler_dual, normal_form_Q

mp.mp.dps = 30


def real_duals(Lam, a, M, s, omega, m, lam, tol=1e-6):
    data = teukolsky_heun(Lam, a, M, s, omega, m, lam)
    t = complex(data['t'])
    duals = {}
    for choice in itertools.product([0, 1], repeat=3):
        e0, e1, et = gauge_to_heun(data, choice)
        h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'], e0, e1, et,
                                      {'e': 2*s + 1, 'zinf': sp.Float(str(data['zinf']), 35)}, s)
        for which in (0, 1):
            ht = euler_dual(h, t, which)
            tr = t.real
            zs = [1 + (tr - 1)*k/13 for k in range(1, 13)]
            imQ = max(abs(normal_form_Q(ht, t, complex(zz)).imag) for zz in zs)
            if imQ < tol:
                key = round(normal_form_Q(ht, t, complex(1 + (tr - 1)*0.37)).real, 6)
                duals.setdefault(key, (choice, which, ht))
    return data, list(duals.values())


def P4_of(ht, t):
    t = t.real
    mu = lambda z: z*(z - 1)*(t - z)
    mup = lambda z: -3*z**2 + 2*(1 + t)*z - t
    mupp = lambda z: -6*z + 2*(1 + t)
    def P4(z):
        Q = normal_form_Q(ht, complex(t), complex(z)).real
        return mu(z)**2*Q - (mup(z)**2 - 2*mu(z)*mupp(z))/4
    # fit degree-4 polynomial through 9 points, check residual
    import numpy as np
    zs = np.linspace(1.05, t - 0.05, 9)
    vals = np.array([P4(z) for z in zs])
    coef = np.polyfit(zs, vals, 6)
    return coef, P4


def analyse(Lam, a, M, s, omega, m, lam, verbose=True):
    import numpy as np
    data, duals = real_duals(Lam, a, M, s, omega, m, lam)
    t = complex(data['t']).real
    out = []
    for choice, which, ht in duals:
        coef, P4 = P4_of(ht, complex(t))
        deg_hi = max(abs(coef[0]), abs(coef[1]))  # z^6, z^5 coefficients should vanish
        p4 = np.poly1d(coef[2:])
        Pj = [p4(0.0), p4(1.0), p4(t)]
        rec = dict(choice=choice, which='AB'[which], deg_hi=deg_hi, P4_at=Pj, options=[])
        if min(Pj) >= -1e-8:
            sq = [np.sqrt(max(x, 0)) for x in Pj]
            for signs in itertools.product([1, -1], repeat=2):  # overall sign of qt irrelevant: fix sign at z=0
                vals = [sq[0], signs[0]*sq[1], signs[1]*sq[2]]
                qt = np.poly1d(np.polyfit([0.0, 1.0, t], vals, 2))
                rem = p4 - qt*qt
                mu = np.poly1d([-1, 1 + t, -t, 0])
                Lneg, remainder = np.polydiv(rem, mu)  # rem = -Lt mu
                Lt = -Lneg
                zz = np.linspace(1, t, 400)
                sr = bool(np.any(qt(zz[:-1])*qt(zz[1:]) < 0))
                rec['options'].append(dict(qt=qt.coeffs, Lt=np.atleast_1d(Lt.coeffs if hasattr(Lt, 'coeffs') else Lt),
                                            superradiant=sr, rem_err=float(np.max(np.abs(remainder)))))
        out.append(rec)
    if verbose:
        print(f't = {t:.6f}')
        for rec in out:
            print(f"dual gauge={rec['choice']} {rec['which']}: |z^6,z^5 coeff| = {rec['deg_hi']:.1e}; "
                  f"P4(0),P4(1),P4(t) = {[round(float(x), 5) for x in rec['P4_at']]}")
            for o in rec['options']:
                print(f"     qt coeffs {np.round(o['qt'], 5)}, Lt coeffs {np.round(o['Lt'], 5)}, "
                      f"qt changes sign on (1,t): {o['superradiant']}, division remainder {o['rem_err']:.1e}")
    return out


if __name__ == '__main__':
    Lam, a, M = 0.05, 0.6, 1.0
    s, m, lam = -2, 2, 7.3
    r3, r0, r1, r2 = kds_roots(Lam, a, M)
    Om1 = a/(r1**2 + a**2); Om2 = a/(r2**2 + a**2)
    for frac in (0.1, 0.5, 0.9):
        omega = complex(m*float(Om2 + frac*(Om1 - Om2)))
        print(f'\n=== omega = {omega.real:.5f} (fraction {frac} through original superradiant band)')
        analyse(Lam, a, M, s, omega, m, lam)
