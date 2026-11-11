"""Step 19: L2 = n0/z + n1 + n2 z.  Report signs and check the sufficient conditions for L2' <= 0 on [1,t]
(n2 <= 0 and n0 >= 0, or n2 t^2 <= n0), together with L2(t) > 0."""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from step13_largem import leading
    from step17_ellmin import ell_min
from step5_window import ctdc_window
P = np.polynomial.polynomial
rows = []
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    for wf in (0.001, 0.2, 0.5, 0.8, 0.999):
        for ellmode in ('min', 'big'):
            w = lo + wf*(hi - lo); lm = ell_min(Lam, a, w) if ellmode == 'min' else 50.0
            d = leading(Lam, a, wf, lm)
            t, qa, qb, p4 = d['t'], d['qa'], d['qb'], d['p4']
            q1 = (qb - qa)/(t - 1); q = np.array([qa - q1, q1])
            bump = P.polymul([-1.0, 1.0], [t, -1.0])
            N, rem = P.polydiv(P.polysub(P.polymul(q, q), p4), bump)
            n0, n1, n2 = (list(N) + [0, 0, 0])[:3]
            cond = (n2 <= 0 and n0 >= 0) or (n2*t*t <= n0)
            Lt = n0/t + n1 + n2*t
            rows.append((Lam, a, wf, ellmode, n0, n1, n2, t, cond, Lt))
print(f'{"Lam":>7} {"a":>6} {"wfrac":>6} {"ell":>4} {"n0":>11} {"n1":>11} {"n2":>11} {"n2 t^2 <= n0?":>14} {"L2(t)":>10}')
for r in rows:
    print(f'{r[0]:7g} {r[1]:6g} {r[2]:6g} {r[3]:>4} {r[4]:+11.4e} {r[5]:+11.4e} {r[6]:+11.4e} {str(r[8]):>14} {r[9]:+10.3e}')
print('all monotone:', all(r[8] for r in rows), ' all L2(t)>0:', all(r[9] > 0 for r in rows),
      ' sign(n0)>0 all:', all(r[4] > 0 for r in rows), ' sign(n2):', set(np.sign([r[6] for r in rows])))
