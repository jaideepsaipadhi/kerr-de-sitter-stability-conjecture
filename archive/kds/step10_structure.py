"""Step 10: look for structure.  Does the simplest choice of qt already work?
  family A: g = 0 (qt is the straight line through the two horizon values)
  family B: g a constant, scanned
For each case: does score > 0, and does it pass the cutting-plane fine check?"""
import sys, numpy as np
from step6_hintz_lp import dual_hintz_form
from step5_window import ctdc_window
from step9_freeq import make_form, zero_in, score, valid_kds
import step8_piecewise as S8
P = np.polynomial.polynomial

def fine_ok(info):
    for side, key in (('cosmo', 'solc'), ('event', 'sole')):
        sol = info[key]; extra = {}
        for it in range(8):
            fc, bad = S8.fine_check(info['form'], info['z0'], side, sol, return_points=True)
            if not bad: break
            for i_, pts in bad.items():
                extra[i_] = np.concatenate([extra.get(i_, np.array([])), pts])
            tau, sol2 = S8.piecewise_lp(info['form'], info['z0'], side, 8, 4, grading=1.0, extra=extra)
            if tau is None or tau < 0: return False
            sol = sol2
        if bad: return False
    return True

cases = [(0.05, 0.6, 1, 2.0, 0.05), (0.05, 0.6, 1, 2.0, 0.95), (0.02, 0.9, 1, 50.0, 0.95), (0.02, 0.9, 1, 2.0, 0.95),
         (0.001, 0.95, 2, 10.0, 0.35), (0.001, 0.95, 2, 10.0, 0.65), (0.001, 0.95, 1, 2.0, 0.95),
         (0.05, 0.9877, 1, 2.0, 0.05), (0.05, 0.9877, 4, 50.0, 0.95), (0.1, 0.3, 4, 2.0, 0.95), (0.02, 0.9, 8, 120.0, 0.5)]
for (Lam, a, m, lam, fr) in cases:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    base = dual_hintz_form(Lam, a, 1.0, -2, m*(lo + fr*(hi - lo)), m, lam, sign0=-1)
    t = base['t']
    sA, iA = score(base, [0.0])
    okA = sA > 0 and fine_ok(iA)
    # family B: constant g scan
    q1, qt_ = P.polyval(1.0, base['q']), P.polyval(t, base['q'])
    scale = abs(qt_ - q1)/(t - 1)**2
    best = None
    for gc in np.linspace(-4, 4, 17)*scale:
        sB, iB = score(base, [gc])
        if best is None or sB > best[0]: best = (sB, gc, iB)
    okB = best[0] > 0 and fine_ok(best[2])
    zA = iA['z0'] if iA else float('nan')
    zlin = 1 + (t - 1)*q1/(q1 - qt_)
    print(f"{(Lam,a,m,lam,fr)}: t={t:.3f} q(1)={q1:+.3f} q(t)={qt_:+.3f} | linear qt: tau={sA:+.1e} fine={okA} z0={zA:.4f} "
          f"| best const g/scale={best[1]/scale:+.2f}: tau={best[0]:+.1e} fine={okB}", flush=True)
