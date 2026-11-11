"""
Step 32: the event-side conditions as a 4-parameter problem.

Exact finite-m structure (step 30-31, fitted to 1e-15): L(z) = A/z + B - 3z, with
  A = m^2 n0^(2) > 0,   B = m^2 ell C0 + 1 + t + 3(1-alpha) C0 > 0.
q1 cancels from both event-side conditions (C and D are each linear in q1), so with s = z0 - z:
  C <= 0  <=>  s h <= H/mu,              H = 2 mu + s mu'
  D <= 0  <=>  (1/2)(mu' h + mu (h' + h^2)) <= s (L h + L'),   where h = f'/f
Also needed: h(z0) >= 0 (seam), f(1) = 0 (h ~ p/(z-1) near 1), and L > 0 on [1, z0].
Multiplier: h = (1+delta) g (1 - exp(-beta kappa s^2)) + p (1-u)/(z-1),  g = -L'/L,  kappa = L(z0)/mu(z0).
We scan (t, z0, A, B) over wide ranges, including values far below anything physical, and report where no
(delta, p, beta) works.
"""
import numpy as np


def check(t, z0, A, B, delta, p, beta, zdecay=False, rational=False):
    mu = lambda z: z*(z - 1)*(t - z); mup = lambda z: -3*z**2 + 2*(1 + t)*z - t
    L = lambda z: A/z + B - 3*z; Lp = lambda z: -A/z**2 - 3; Lpp = lambda z: 2*A/z**3
    x0 = z0 - 1
    s_ = np.concatenate([np.geomspace(1e-13*x0, 1e-2*x0, 3000), np.linspace(1e-2*x0, x0*(1 - 1e-9), 3000)])
    z = np.sort(np.concatenate([z0 - s_, 1 + np.geomspace(1e-12*x0, 1e-2*x0, 2000)]))
    z = z[(z > 1) & (z < z0)]
    s = z0 - z
    if np.min(L(z)) <= 0:
        return False, 'L<=0'
    kap = L(z0)/mu(z0)
    g = -Lp(z)/L(z); gp = -(Lpp(z)*L(z) - Lp(z)**2)/L(z)**2
    if rational:
        w = beta*kap*s**2; psi = (1 + delta)*w/(1 + w); psip = -(1 + delta)*2*beta*kap*s/(1 + w)**2
    else:
        E = np.exp(-beta*kap*s**2); psi = (1 + delta)*(1 - E); psip = -(1 + delta)*E*2*beta*kap*s  # d/dz
    u = (z - 1)/x0
    if zdecay:
        q = (z - 1)*z
        h = g*psi + p*(1 - u)/q; hp = gp*psi + g*psip + p*(-(1/x0)/q - (1 - u)*(2*z - 1)/q**2)
    else:
        h = g*psi + p*(1 - u)/(z - 1); hp = gp*psi + g*psip - p/(z - 1)**2
    muv, mupv = mu(z), mup(z)
    H = 2*muv + s*mupv
    Cc = s*h*muv - H                 # need <= 0
    nC = np.abs(s*h*muv) + np.abs(H)
    Dd = 0.5*(mupv*h + muv*(hp + h**2)) - s*(L(z)*h + Lp(z))   # need <= 0
    nD = 0.5*(np.abs(mupv*h) + np.abs(muv*(hp + h**2))) + s*(np.abs(L(z)*h) + np.abs(Lp(z)))
    ok = np.max(Cc/nC) <= 0 and np.max(Dd/nD) <= 0
    return ok, (float(np.max(Cc/nC)), float(np.max(Dd/nD)))


