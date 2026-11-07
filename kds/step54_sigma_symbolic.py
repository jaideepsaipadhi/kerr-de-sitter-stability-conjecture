"""B1 sigma-uniformity at a fixed point and fixed m (K = 4): build C(sigma), G(sigma) EXACTLY as polynomials of degree <= 2 by
interpolating three exact builds, then compute GENERIC ranks over Q(i)(sigma) by fraction-free elimination (sympy DomainMatrix),
and the exceptional polynomials where rank C or rank G_L can drop."""
import sys, time
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.matrices import DomainMatrix
import step52_exact_modes as S
m = int(sys.argv[1]) if len(sys.argv) > 1 else 1
K = 4; t0 = time.time()
pts = [Fr(0), Fr(1), Fr(2)]
builds = [S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), (sv, Fr(0)), m, K) for sv in pts]
print('built 3 exact samples', f'{time.time()-t0:.0f}s', flush=True)
sig = sp.Symbol('s')
dom = sp.QQ_I[sig]
def poly_mat(key_r, key_i, rows=None, cols=None):
    A0 = builds[0][key_r]; nr, nc = A0.nrows(), A0.ncols()
    rows = range(nr) if rows is None else rows; cols = range(nc) if cols is None else cols
    out = []
    for i in rows:
        row = []
        for j in cols:
            v = [complex(0)]*3
            vals = [sp.Rational(int(b[key_r][i, j].p), int(b[key_r][i, j].q)) + sp.I*sp.Rational(int(b[key_i][i, j].p), int(b[key_i][i, j].q)) for b in builds]
            # quadratic through (0,v0),(1,v1),(2,v2)
            c0 = vals[0]; c2 = (vals[2] - 2*vals[1] + vals[0])/2; c1 = vals[1] - vals[0] - c2
            row.append(dom.convert(sp.expand(c0 + c1*sig + c2*sig**2)))
        out.append(row)
    return DomainMatrix(out, (len(out), len(out[0]) if out else 0), dom)
res0 = builds[0]; nk, comps, monsK, mons = res0['nk'], res0['comps'], res0['monsK'], res0['mons']
J = 1
low = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J]
high = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J]
C = poly_mat('CR', 'CI'); CH = poly_mat('CR', 'CI', rows=high); GL = poly_mat('GR', 'GI', cols=low)
print('matrices', C.shape, CH.shape, GL.shape, f'{time.time()-t0:.0f}s', flush=True)
Kf = dom.get_field()
for name, A in (('C', C), ('C_H', CH), ('G_L', GL)):
    r = A.convert_to(Kf).rank()
    print(f'generic rank {name} over Q(i)(sigma) = {r}', f'{time.time()-t0:.0f}s', flush=True)
