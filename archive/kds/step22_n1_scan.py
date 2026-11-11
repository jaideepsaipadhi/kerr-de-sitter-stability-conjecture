"""Step 22: where is n1(ell_min) smallest?  Report n1, n0, and n1 relative to n0/t (the L2(t) > 0 scale) across the window."""
import warnings; warnings.filterwarnings('ignore')
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from step20_exact_n import n_coeffs
    from step17_ellmin import ell_min
from step5_window import ctdc_window
res = []
for (Lam, a) in [(0.05, 0.6), (0.02, 0.9), (0.1, 0.3), (0.001, 0.95), (0.05, 0.9877), (0.08, 0.5), (0.01, 0.99), (0.0001, 0.99), (0.001, 0.999), (0.1, 0.05), (0.11, 0.01), (0.111, 0.001)]:
    try: lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    except Exception: continue
    for wf in (0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999):
        w = lo + wf*(hi - lo); lm = ell_min(Lam, a, w)
        n0, n1, n2, rem, t, _, _ = n_coeffs(Lam, a, w, lm)
        res.append((n1/(abs(n0) + abs(n1)), Lam, a, wf, n0, n1, t, lm))
res.sort()
for r in res[:12]:
    print(f'rel n1={r[0]:+.3e}  Lam={r[1]} a={r[2]} wfrac={r[3]}  n0={r[4]:+.3e} n1={r[5]:+.3e} t={r[6]:.3f} ell_min={r[7]:.4f}')
