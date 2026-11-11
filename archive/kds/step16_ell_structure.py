"""
Step 16: check that L2(z; ell) = L2^0(z) + ell*C0 with C0 a positive constant (independent of z,
and of w), and that H(z) = 2 mu - (z - z0) mu' > 0 (the ell -> infinity condition).
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step13_largem import leading
from step9_freeq import valid_kds
P = np.polynomial.polynomial

def L2_on_grid(Lam, a, wf, ell, zz):
    d = leading(Lam, a, wf, ell)
    t, qa, qb, p4 = d['t'], d['qa'], d['qb'], d['p4']
    q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1])
    bump = P.polymul([-1.0, 1.0], [t, -1.0])
    Nq, _ = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
    return P.polyval(zz, Nq)/zz, t, q

for (Lam, a) in [(0.05, 0.6), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99), (0.1, 0.3)]:
    if not valid_kds(Lam, a): continue
    C0s = []
    for wf in (0.1, 0.5, 0.9):
        _, t, _ = L2_on_grid(Lam, a, wf, 1.0, np.array([1.5]))
        zz = np.linspace(1.001, t - 0.001, 50)
        A, _, q = L2_on_grid(Lam, a, wf, 1.0, zz)
        B, _, _ = L2_on_grid(Lam, a, wf, 11.0, zz)
        dif = (B - A)/10.0
        C0s.append((dif.mean(), dif.std()/abs(dif.mean())))
        z0 = -q[0]/q[1]
        mu = lambda z: z*(z - 1)*(t - z); mup = lambda z: -3*z**2 + 2*(1 + t)*z - t
        H = 2*mu(zz) - (zz - z0)*mup(zz)
    print(f'Lam={Lam} a={a}: C0 per w = {[f"{c:.6f} (z-spread {s:.1e})" for c, s in C0s]}; min H = {H.min():+.3e}')
