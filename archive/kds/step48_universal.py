"""Step 48: universal ODE multiplier.  f' = w = (1/mu) exp(int^z 2c (z0-eta)/mu(eta) d eta), f = int_1^z w.
Event lemma <= [exists c in (0, min(c0, 2.5t - 1.25)), c0 = (t-1)/(2 x0)] with
   Psi_c(z) := f'/f  >=  max(1/z, 3/(5.5t - 1.25 - 3z - c))   on (1, z0].
Scan (t, zf); report the best c and the worst normalised margin."""
import numpy as np
from scipy.integrate import cumulative_trapezoid as ctz

def psi_margin(t, zf, c, n=6000):
    x0 = zf*(t - 1); z0 = 1 + x0
    x = np.concatenate([np.geomspace(1e-12*x0, 1e-3*x0, n//2), np.linspace(1e-3*x0, x0, n//2)[1:]])
    z = 1 + x
    mu = z*x*(t - z)
    k = 2*c*(z0 - z)/mu
    # exponent: int k dz; near x->0, k ~ p/x with p = 2c x0/(t-1); subtract p/x analytically
    p = 2*c*x0/(t - 1)
    reg = k - p/x
    I = ctz(reg, z, initial=0) + p*np.log(x/x[0])
    w = np.exp(I - I.max())/mu
    # f = int_1^z w, with w ~ C x^{p-1} near 0: add the analytic piece below x[0]
    f = ctz(w, z, initial=0) + w[0]*x[0]/p
    Psi = w/f
    rhs = np.maximum(1/z, 3/(5.5*t - 1.25 - 3*z - c))
    mD = np.min((Psi - rhs)/(Psi + rhs))
    F = ctz(f, z, initial=0) + f[0]*x[0]/(p + 1)
    s = z0 - z; mup = -3*z**2 + 2*(1 + t)*z - t; H = 2*mu + s*mup
    lhsC = 2*c*s*(s*f + F)/f
    mC = np.min((H - lhsC)/(H + lhsC))
    return min(mD, mC)

worst = (9, None); fails = []
for t in np.geomspace(1.001, 1e4, 60):
    for zf in np.concatenate([np.geomspace(1e-6, 0.1, 15), np.linspace(0.15, 0.999999, 25)]):
        x0 = zf*(t - 1); cmax = 2.5*t - 1.25
        best = max(((psi_margin(t, zf, c), c) for c in cmax*np.geomspace(1e-4, 0.999, 30)), key=lambda v: v[0])
        if best[0] < worst[0]: worst = (best[0], (t, zf, best[1], cmax))
        if best[0] < 0: fails.append((round(t, 4), zf, best))
print('worst normalised margin over (t, zf):', worst)
print('failures:', len(fails), fails[:10])
