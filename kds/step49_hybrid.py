"""Step 49: two-zone hybrid multiplier.  Zone I [1, z1]: explicit ODE multiplier (constant c).  Zone II [z1, z0]:
h = Psi(z1) phi(s)/phi(s1), phi = 1 - exp(-beta kappa s^2), kappa = L(z0)/mu(z0)  (continuous at z1, h(z0) = 0).
Conditions: C: s h mu <= H;  D: (1/2)(mu' h + mu (h' + h^2)) <= s (L h + L').  Random (t, zf, A, B >= floor)."""
import numpy as np
from scipy.integrate import cumulative_trapezoid as ctz
rng = np.random.default_rng(1)

def zoneI(t, x0, c, x):
    z = 1 + x; mu = z*x*(t - z); z0 = 1 + x0
    p = 2*c*x0/(t - 1); k = 2*c*(z0 - z)/mu
    I = ctz(k - p/x, z, initial=0) + p*np.log(x/x[0])
    w = np.exp(I - I.max())/mu
    f = ctz(w, z, initial=0) + w[0]*x[0]/p
    return w/f

def margins(t, zf, A, B, c, frac, beta, n=4000):
    x0 = zf*(t - 1); z0 = 1 + x0; s1 = frac*x0
    x = np.concatenate([np.geomspace(1e-12*x0, 1e-3*x0, n//2), np.linspace(1e-3*x0, x0*(1 - 1e-12), n)[1:]])
    z = 1 + x; s = z0 - z
    mu = z*x*(t - z); mup = -3*z**2 + 2*(1 + t)*z - t; H = 2*mu + s*mup
    L = A/z + B - 3*z; Lp = -A/z**2 - 3
    Psi = zoneI(t, x0, c, x)
    I1 = s >= s1
    z1i = np.argmin(np.abs(s - s1))
    kap = (A/z0 + B - 3*z0)/(z0*x0*(t - z0))
    phi = 1 - np.exp(-beta*kap*s**2); phi1 = 1 - np.exp(-beta*kap*s1**2)
    h = np.where(I1, Psi, Psi[z1i]*phi/phi1)
    hp = np.gradient(h, z)
    Cc = s*h*mu - H
    Dd = 0.5*(mup*h + mu*(hp + h*h)) - s*(L*h + Lp)
    nD = 0.5*(np.abs(mup*h) + np.abs(mu*hp) + mu*h*h) + s*(np.abs(L*h) + np.abs(Lp))
    sel = slice(5, -5)       # drop gradient end effects
    return max(np.max((Cc/(np.abs(s*h*mu) + H))[sel]), np.max((Dd/nD)[sel]))

tot = ok = 0; fails = []
for _ in range(400):
    t = 1 + 10**rng.uniform(-2, 3); zf = rng.uniform(0.02, 0.98)
    A = 10**rng.uniform(-2, 5); B = (5.5*t - 1.25)*10**rng.uniform(0, 3)
    best = 9
    for c in (2.5*t - 1.25)*np.array([0.02, 0.1, 0.3, 0.6]):
        for frac in (0.05, 0.2, 0.5):
            for beta in (0.25, 1, 4, 25):
                v = margins(t, zf, A, B, c, frac, beta)
                best = min(best, v)
                if best < 0: break
            if best < 0: break
        if best < 0: break
    tot += 1; ok += best < 0
    if best >= 0: fails.append((round(t, 3), round(zf, 3), f'{A:.2g}', f'{B:.3g}', round(best, 4)))
print(f'hybrid closes {ok}/{tot}');
for f in fails[:15]: print('  fail', f)
