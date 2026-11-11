"""Step 43: outside the window.  Real dual, Wronskian conservation: Im e_1 (t-1)|c1|^2 + Im e_t t(t-1)|c_t|^2 = 0,
e_j = 2(kappa_j + kappa_3).  Mode stability follows wherever sign Im e_1 = sign Im e_t.  Check that the set where the
signs DIFFER is exactly the CTdC window (lo, hi) in omega/m, in our conventions, for s = -2 and +2 (spin-free)."""
import numpy as np
from step3_dual import kds_roots
from step5_window import ctdc_window
from step9_freeq import valid_kds
def ime(Lam, a, om, m, R):
    r3, r0, r1, r2 = R
    Xi = 1 + Lam*a*a/3
    Dp = lambda r: -(4*Lam/3)*r**3 + 2*(1 - Lam*a*a/3)*r - 2*1.0
    kap = lambda r: Xi*(om*(r*r + a*a) - a*m)/Dp(r)          # Im kappa
    return 2*(kap(r1) + kap(r3)), 2*(kap(r2) + kap(r3))
bad = n = 0; worst = 0
for Lam in np.geomspace(1e-5, 0.11, 15):
    for a in np.linspace(0.02, 1.2, 30):
        if not valid_kds(Lam, a): continue
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0); R = [float(x) for x in kds_roots(Lam, a, 1.0)]
        for m in (1, 2, 5):
            for u in np.linspace(-0.5, 0.6, 221):          # u = omega/m
                e1, et = ime(Lam, a, u*m, m, R)
                inside = lo < u < hi; differ = e1*et < 0
                n += 1
                if inside != differ and min(abs(u - lo), abs(u - hi)) > 1e-9:
                    bad += 1
print(f'{n} (geometry, m, omega) samples; mismatches between "signs differ" and CTdC window: {bad}')
