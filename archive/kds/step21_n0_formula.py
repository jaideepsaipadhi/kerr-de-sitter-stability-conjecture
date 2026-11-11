"""Step 21: check n0 = t (X0^2 - (X1 + Xt)^2) with X_j = |kappa_j + kappa_3|, and look at the signs of
u_j = Im kappa_j inside the window, to find the hypothesis under which X0 >= X1 + Xt."""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from step20_exact_n import n_coeffs
    from step17_ellmin import ell_min
from step5_window import ctdc_window
from step3_dual import kds_roots

def kappas(Lam, a, w):
    r3, r0, r1, r2 = [float(x) for x in kds_roots(Lam, a, 1.0)]
    Xi = 1 + Lam*a**2/3
    Dp = lambda r: 2*r*(1 - Lam*r**2/3) - (r**2 + a**2)*2*Lam*r/3 - 2
    return [Xi*(w*(r**2 + a**2) - a)/Dp(r) for r in (r0, r1, r2, r3)]   # u_j = Im kappa_j

worst = 0; signs = set()
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.001, 0.95), (0.05, 0.9877), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999), (0.1, 0.05)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    for wf in np.linspace(0.001, 0.999, 9):
        w = lo + wf*(hi - lo)
        n0, n1, n2, rem, t, _, _ = n_coeffs(Lam, a, w, ell_min(Lam, a, w))
        u0, u1, ut, u3 = kappas(Lam, a, w)
        X0, X1, Xt = abs(u0 + u3), abs(u1 + u3), abs(ut + u3)
        pred = t*(X0**2 - (X1 + Xt)**2)
        worst = max(worst, abs(pred - n0)/(abs(n0) + 1e-300))
        signs.add((np.sign(u0), np.sign(u1), np.sign(ut), np.sign(u3), np.sign(u0 + ut), np.sign(u0 + u1)))
print(f'max relative error of the n0 formula: {worst:.2e}')
print('sign patterns (u0, u1, ut, u3, u0+ut, u0+u1) seen in the window:', signs)
