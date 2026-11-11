"""Step 46: redo the cheap numerical checks on the CORRECT subextremal region (four real roots; Lam M^2 up to ~0.18)."""
import numpy as np
from step5_window import ctdc_window
from step9_freeq import valid_kds
from step40_angular_E import E
from step41_angular_m1 import bound
def roots(L, a): return np.sort(np.roots([-L/3, 0, 1 - L*a*a/3, -2, a*a]).real)
G = [(L, a) for L in np.linspace(1e-4, 0.2, 80) for a in np.linspace(0.01, 1.6, 160) if valid_kds(L, a)]
print('geometries:', len(G), ' max Lam:', max(g[0] for g in G), ' max alpha:', max(L*a*a/3 for L, a in G))
xs = np.linspace(-1, 1, 1001)
bad_inc = bad_sign = bad_om0 = 0; endp = 0; ximax = 0; ang1 = 9; ang2 = 9; ang3 = 9
for (L, a) in G:
    lo, hi, _, _ = ctdc_window(L, a, 1.0); r3, r0, r1, r2 = roots(L, a)
    al = L*a*a/3; Xi = 1 + al; OmH = a/(r1**2 + a**2)
    Dp = lambda r: -(4*L/3)*r**3 + 2*(1 - al)*r - 2
    k = lambda r, u: Xi*(u*(r*r + a*a) - a)/Dp(r)
    bad_inc += not (0 <= lo < hi <= OmH + 1e-14)
    bad_om0 += lo < 0 < hi
    for u in np.linspace(-0.3, 0.5, 161):
        e1 = 2*(k(r1, u) + k(r3, u)); et = 2*(k(r2, u) + k(r3, u))
        if min(abs(u - lo), abs(u - hi)) > 1e-9 and ((lo < u < hi) != (e1*et > 0)): bad_sign += 1
    for u in (lo, hi):
        e1 = 2*(k(r1, u) + k(r3, u)); et = 2*(k(r2, u) + k(r3, u))
        endp = max(endp, min(abs(e1), abs(et))/max(abs(e1), abs(et)))
    ximax = max(ximax, a*hi)
    for u in np.linspace(lo, hi, 9):
        for s in (-2, 2):
            ang1 = min(ang1, bound(al, a*u, 1, s))
            ang2 = min(ang2, 2 + E(xs, al, 2*a*u, 2, s).min())
            ang3 = min(ang3, 8 + E(xs, al, 3*a*u, 3, s).min())
print('window inside [0, Omega_H] failures:', bad_inc, '| omega=0 inside window:', bad_om0)
print('outside-window sign test failures (window <=> signs agree):', bad_sign)
print('endpoints: max min/max |Im e| ratio:', endp)
print('max xi/m over window:', ximax)
print('angular lower bounds: |m|=1 two-level', ang1, '| |m|=2 pointwise', ang2, '| |m|=3 pointwise', ang3, '(need >= 1.5)')
