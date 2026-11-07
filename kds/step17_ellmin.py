"""
Step 17: the leading-order conditions at the physical lower end ell = ell_min(w).

  ell_min(w) = min over theta of  Xi^2 (1 - a w sin^2 theta)^2 / (Delta_theta sin^2 theta),
  Delta_theta = 1 + (Lam a^2/3) cos^2 theta.

This is the classical turning-point bound for the angular equation at large m.  It has to be
proved for the spin-weighted KdS angular operator; a 1/m correction would not affect the
leading-order statement.  Phi and Psi are monotone in ell (step 16), so ell_min and
ell = infinity (proved by hand) are the only cases to check.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import minimize_scalar
from step13_largem import check
from step5_window import ctdc_window
from step14_cosmo_general import psi_margin
from step9_freeq import valid_kds

def ell_min(Lam, a, w):
    al = Lam*a**2/3; Xi = 1 + al
    f = lambda th: Xi**2*(1 - a*w*np.sin(th)**2)**2/((1 + al*np.cos(th)**2)*np.sin(th)**2)
    return minimize_scalar(f, bounds=(1e-3, np.pi/2), method='bounded').fun

worst = (np.inf, None)
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999)]:
    if not valid_kds(Lam, a): print('invalid', Lam, a); continue
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    out = []
    for wf in (0.001, 0.02, 0.2, 0.5, 0.8, 0.98, 0.999):
        w = lo + wf*(hi - lo)
        lm = ell_min(Lam, a, w)
        r = check(Lam, a, wf, lm)
        ps, L2min, _, _, _, _ = psi_margin(Lam, a, wf, lm)
        ev = r.get('event_Phi_rel', np.nan); evL = r.get('event_L2min', np.nan)
        out.append(f'{wf}: ell_min={lm:.3f} Phi={ev:+.3f} L2ev={evL:+.2e} Psi={ps:+.3f} L2c={L2min:+.2e}')
        m_ = min(ev, ps)
        if m_ < worst[0]: worst = (m_, (Lam, a, wf, lm))
        if not (r['signs'] and ev >= 0 and evL > 0 and ps >= 0 and L2min > 0):
            out[-1] += '   <-- FAIL'
    print(f'Lam={Lam} a={a}'); [print('   ', o) for o in out]
print('worst relative margin:', worst)
