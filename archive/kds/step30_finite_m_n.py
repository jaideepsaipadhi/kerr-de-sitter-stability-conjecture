"""Step 30: exact finite-m coefficients of L = n0/z + n1 + n2 z (no extrapolation).  Fit each n_i as a polynomial in m
(degree <= 2) at fixed w = omega/m and ell = lambda/m^2, and compare the m^2 coefficients with the closed forms
n0^(2) = 2 t u3 (X0+X1+Xt), n1^(2) = ell C0, n2^(2) = 0."""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
P = np.polynomial.polynomial

def ncoef(Lam, a, w, ell, m, s=-2):
    b = dual_hintz_form(Lam, a, 1.0, s, m*w, m, ell*m*m, sign0=-1)
    t, p4 = b['t'], b['p4']
    qa, qb = P.polyval(1.0, b['q']), P.polyval(t, b['q'])
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1])
    N, rem = P.polydiv(P.polysub(P.polymul(q, q), p4), P.polymul([-1.0, 1.0], [t, -1.0]))
    N = list(N) + [0, 0, 0]
    return np.array(N[:3]), float(np.max(np.abs(rem))), t

for (Lam, a, wf, ell, s) in [(0.05, 0.6, 0.5, 1.2, -2), (0.001, 0.95, 0.2, 1.0, -2), (0.01, 0.99, 0.5, 0.9, -2), (0.05, 0.6, 0.5, 1.2, +2)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0); w = lo + wf*(hi - lo)
    ms = np.array([3, 5, 8, 13, 21, 34])
    Ns = np.array([ncoef(Lam, a, w, ell, int(m), s)[0] for m in ms])
    print(f'Lam={Lam} a={a} wf={wf} ell={ell} s={s}:')
    for i in range(3):
        cf = np.polyfit(ms, Ns[:, i], 2)
        resid = np.max(np.abs(np.polyval(cf, ms) - Ns[:, i]))/np.max(np.abs(Ns[:, i]))
        print(f'    n{i}(m) = {cf[0]:+.6f} m^2 {cf[1]:+.6f} m {cf[2]:+.6f}   (fit residual {resid:.1e})')
