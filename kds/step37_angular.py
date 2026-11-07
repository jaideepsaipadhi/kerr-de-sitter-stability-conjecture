"""
Step 37: lowest spin-weighted KdS angular eigenvalue vs the value needed by the event lemma.

Angular equation (Suzuki-Takasugi-Umetsu, gr-qc/9805064, eq. 3.1), x = cos(theta), xi = a omega:
  d/dx[(1+al x^2)(1-x^2) S'] + { lam98 - s(1-al) + ((1+al)^2/al) xi^2 - 2 al x^2
     + (1+al)/(1+al x^2) [2s(al m - (1+al) xi) x - ((1+al)^2/al) xi^2 - 2 m (1+al) xi + s^2]
     - (1+al)^2 m^2/((1+al x^2)(1-x^2)) - (1+al)(s^2 + 2 s m x)/(1-x^2) } S = 0.
Our radial separation constant: lam_ours = lam98 - s(1-al)  (matches l(l+1) - s^2 at a = 0).
Regular solutions: S = (1-x)^{|m+s|/2} (1+x)^{|m-s|/2} y(x), with y smooth.  Chebyshev collocation in y.
Requirement (step 35): B = c1 + lam C0 >= 1.3 (4t - 1/2), i.e. lam >= lam_req = (1.3(4t - 1/2) - c1)/C0.
"""
import numpy as np
from scipy.linalg import eig
from step3_dual import kds_roots
from step5_window import ctdc_window
from step9_freeq import valid_kds


def cheb(N):
    x = np.cos(np.pi*np.arange(N + 1)/N)
    c = np.hstack([2, np.ones(N - 1), 2])*(-1)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def lowest_lambda(al, xi, m, s, N=60):
    """Collocation on all Chebyshev nodes.  The operator is multiplied by (1-x^2), so the endpoint rows impose the
    regularity (Frobenius) condition automatically.  Coefficients are evaluated at nodes pulled in by 1e-11, which
    resolves the cancelling singular terms numerically."""
    D, x = cheb(N); D2 = D @ D
    xe = x*(1 - 1e-11)
    k1, k2 = abs(m + s)/2, abs(m - s)/2
    p = (1 + al*xe**2)*(1 - xe**2); pp = 2*al*xe*(1 - xe**2) - 2*xe*(1 + al*xe**2)
    l1 = -k1/(1 - xe) + k2/(1 + xe)
    l2 = l1**2 + (-k1/(1 - xe)**2 - k2/(1 + xe)**2)
    V = (-s*(1 - al) + ((1 + al)**2/al)*xi**2 - 2*al*xe**2
         + (1 + al)/(1 + al*xe**2)*(2*s*(al*m - (1 + al)*xi)*xe - ((1 + al)**2/al)*xi**2 - 2*m*(1 + al)*xi + s**2)
         - (1 + al)**2*m**2/((1 + al*xe**2)*(1 - xe**2)) - (1 + al)*(s**2 + 2*s*m*xe)/(1 - xe**2))
    w = (1 - xe**2)
    c2 = w*p; c1 = w*(pp + 2*p*l1); c0 = w*(p*l2 + pp*l1 + V)
    Aop = np.diag(c2) @ D2 + np.diag(c1) @ D + np.diag(c0)
    ev = eig(-Aop, np.diag(w), right=False)
    ev = ev[np.isfinite(ev)]
    ev = np.sort(ev[np.abs(ev.imag) < 1e-6*np.maximum(1, np.abs(ev))].real)
    return ev[0] - s*(1 - al)


if __name__ == '__main__':
    # sanity check: a -> 0 (al -> 0 with xi -> 0) should give lam_ours = l(l+1) - s^2 = 2 for l = 2
    print('check a ~ 0, s = -2, m = 1: lowest lam_ours =', round(lowest_lambda(1e-8, 1e-6, 1, -2), 5), '(expect 2)')
    print('check a ~ 0, s = -2, m = 2: lowest lam_ours =', round(lowest_lambda(1e-8, 1e-6, 2, -2), 5), '(expect 2)')
    print('check a ~ 0, s = -2, m = 3: lowest lam_ours =', round(lowest_lambda(1e-8, 1e-6, 3, -2), 5), '(expect 8)')
    worst = (np.inf, None)
    for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99), (0.001, 0.999), (0.0001, 0.99)]:
        if not valid_kds(Lam, a): continue
        r3, r0, r1, r2 = [float(q) for q in kds_roots(Lam, a, 1.0)]
        al = Lam*a*a/3
        C0 = 3/(Lam*(r1 - r0)*(r2 - r3)); t = (r1 - r3)*(r2 - r0)/((r1 - r0)*(r2 - r3))
        c1 = 1 + t + 3*(1 - al)*C0
        lam_req = (1.3*(4*t - 0.5) - c1)/C0
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
        mins = []
        for m in (1, 2, 3, 4, 6):
            for wf in (0.0, 0.5, 1.0):
                for s in (-2, 2):
                    lamv = lowest_lambda(al, a*m*(lo + wf*(hi - lo)), m, s)
                    mins.append(lamv)
        lmin = min(mins)
        if lmin - lam_req < worst[0]: worst = (lmin - lam_req, (Lam, a))
        print(f'Lam={Lam} a={a}: t={t:.3f} C0={C0:.3f}  lam_req={lam_req:+.3f}   min lam over window, m <= 6, s = +-2: {lmin:.4f}')
    print('smallest margin lam_min - lam_req:', worst)
