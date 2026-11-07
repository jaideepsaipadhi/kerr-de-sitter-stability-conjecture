"""
Step 41: |m| = 1 angular bound via a two-level (gap) lower bound.
H >= H0 + E (form sense; dropped terms are >= 0).  H0 = a=0 spin-weighted operator: levels 2 (l=2), 8 (l=3), ...
psi0 (l=2, s, m=+-1): |psi0|^2 ∝ (1-x)^{|m+s|} (1+x)^{|m-s|}.  With p = <E>, sigma^2 = <E^2> - p^2:
  lam >= smaller eigenvalue of [[2 + p, sigma], [sigma, 8 + inf E]]  (Cauchy-Schwarz on the off-diagonal block).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from step40_angular_E import E
from step3_dual import kds_roots
from step9_freeq import valid_kds

xg, wg = leggauss(200)
def bound(al, xi, m, s):
    w = wg*(1 - xg)**abs(m + s)*(1 + xg)**abs(m - s); w = w/w.sum()
    e = E(xg, al, xi, m, s); p = (w*e).sum(); sig = np.sqrt(max((w*e*e).sum() - p*p, 0))
    lo = 8 + E(np.linspace(-1, 1, 4001), al, xi, m, s).min()
    a11, a22 = 2 + p, lo
    return (a11 + a22)/2 - np.sqrt(((a11 - a22)/2)**2 + sig**2)

worst = (9, None)
for Lam in np.geomspace(1e-6, 0.11, 50):
    for a in np.linspace(0.01, 1.2, 150):
        if not valid_kds(Lam, a): continue
        r3, r0, r1, r2 = [float(q) for q in kds_roots(Lam, a, 1.0)]
        al = Lam*a*a/3; OmH = a/(r1**2 + a**2)
        for u in np.linspace(0, OmH, 30):
            for s in (-2, 2):
                for m in (1, -1):
                    v = bound(al, a*m*u, m, s)
                    if v < worst[0]: worst = (v, Lam, a, u/OmH, s, m)
print('|m| = 1: worst two-level lower bound on lam =', worst)
