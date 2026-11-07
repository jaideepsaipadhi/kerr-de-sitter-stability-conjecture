"""Step 33: (a) the real-dual root is A = 1 + s - 2 kappa_3 and the other root is B = 1 + 2s, exactly at finite m;
(b) c1 = 1 + t + 3(1-alpha) C0 on random geometries;  (c) the bound C0 > (2t-1)/(3(1-alpha)), which gives L(t) > 0."""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, sympy as sp
from step3_dual import teukolsky_heun, heun_params_from_operator, real_branch, kds_roots
from step5_window import ctdc_window
from step30_finite_m_n import ncoef
from step9_freeq import valid_kds

def kap3(Lam, a, om, m):
    r3 = float(kds_roots(Lam, a, 1.0)[0]); Xi = 1 + Lam*a**2/3
    Dp = 2*r3*(1 - Lam*r3**2/3) - (r3**2 + a**2)*2*Lam*r3/3 - 2
    return 1j*Xi*(om*(r3**2 + a**2) - a*m)/Dp

rng = np.random.default_rng(1)
print('(a) exponents:')
for (Lam, a, s, m) in [(0.05, 0.6, -2, 3), (0.001, 0.95, -2, 7), (0.05, 0.6, 2, 3), (0.02, 0.9, 2, 5)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0); om = m*(lo + 0.4*(hi - lo))
    data = teukolsky_heun(Lam, a, 1.0, s, complex(om), m, 9.0)
    e = data['exps']
    h = heun_params_from_operator(data['Pz'], data['Qz'], data['z'], data['t'], e['0'][0], e['1'][0], e['t'][0],
                                  {'e': 2*s + 1, 'zinf': sp.Float(str(data['zinf']), 35)}, s)
    which, ht = real_branch(h, complex(data['t']))
    A = ht['A']; S = h['g'] + h['d'] + h['e'] - 1; Bv = S - A
    k3 = kap3(Lam, a, om, m)
    print(f'   s={s:+d}: |A - (1+s-2k3)| = {abs(A - (1 + s - 2*k3)):.1e},  |B - (1+2s)| = {abs(Bv - (1 + 2*s)):.1e}')

print('(b) c1 = 1 + t + 3(1-alpha) C0 and (c) C0 bound, random geometries:')
worst_c1 = 0; min_ratio = np.inf; cnt = 0
while cnt < 12:
    Lam = 10**rng.uniform(-4, np.log10(0.11)); a = rng.uniform(0.01, 1.0)
    if not valid_kds(Lam, a): continue
    try:
        lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    except Exception: continue
    if not (hi > lo): continue
    cnt += 1
    w = lo + 0.5*(hi - lo)
    n4, _, t = ncoef(Lam, a, w, 1.0, 4); n7, _, _ = ncoef(Lam, a, w, 1.0, 7)
    C0 = (n7[1] - n4[1])/(49 - 16); c1 = n4[1] - 16*C0
    al = Lam*a*a/3
    pred = 1 + t + 3*(1 - al)*C0
    worst_c1 = max(worst_c1, abs(c1 - pred)/abs(pred))
    ratio = C0/((2*t - 1)/(3*(1 - al)))
    min_ratio = min(min_ratio, ratio)
    print(f'   Lam={Lam:.2e} a={a:.3f}: t={t:9.3f} C0={C0:9.4f}  c1 rel.err={abs(c1 - pred)/abs(pred):.1e}  C0 / bound = {ratio:.4f}')
print(f'max c1 relative error {worst_c1:.1e};  min C0/bound {min_ratio:.4f} (need > 1)')
