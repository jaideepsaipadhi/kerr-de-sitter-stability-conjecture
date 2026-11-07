"""Step 31: the constant c1 in n1 = m^2 ell C0 + c1.  How does it depend on w, ell, s?  Compare with t, C0 and the roots."""
import warnings; warnings.filterwarnings('ignore')
import numpy as np
from step30_finite_m_n import ncoef
from step5_window import ctdc_window
from step3_dual import kds_roots
from step16_ell_structure import L2_on_grid
for (Lam, a) in [(0.05, 0.6), (0.001, 0.95), (0.01, 0.99), (0.08, 0.5)]:
    lo, hi, _, _ = ctdc_window(Lam, a, 1.0)
    vals = []
    for wf in (0.2, 0.7):
        for ell in (1.0, 3.0):
            w = lo + wf*(hi - lo)
            n5, _, t = ncoef(Lam, a, w, ell, 5); n9, _, _ = ncoef(Lam, a, w, ell, 9)
            A2 = (n9[1] - n5[1])/(81 - 25); c1 = n5[1] - 25*A2
            vals.append((wf, ell, c1, A2/ell))
    r3, r0, r1, r2 = [float(x) for x in kds_roots(Lam, a, 1.0)]
    C0 = vals[0][3]
    print(f'Lam={Lam} a={a}: t={t:.4f}  C0={C0:.6f}  c1 values (wf, ell, c1): {[(v[0], v[1], round(v[2], 6)) for v in vals]}')
    print(f'      c1 - 3t = {vals[0][2] - 3*t:.6f};  (c1 - 3t)/C0 = {(vals[0][2] - 3*t)/C0:.6f};  c1/C0 = {vals[0][2]/C0:.6f};  '
          f'c1/t = {vals[0][2]/t:.6f}; (c1-2t)/C0={(vals[0][2]-2*t)/C0:.6f}; (c1-(1+t))/C0 = {(vals[0][2]-(1+t))/C0:.6f}')
