"""K = 7, J = 4 rigorous bound: d_4 <= n_L - rank_P C + rank_Q C_H - rank_P G_L, with C_H's exact rank from an
integer (row-scaled) realified fmpz_mat (lighter than fmpq).  Mode coefficients exact; 4D assembly exact (step68)."""
import sys, time, resource, numpy as np, flint
from fractions import Fraction as Fr
from math import lcm
import step68_4d_exact as X, step67_4d_from_modes as F
K = 7; J = int(sys.argv[1]) if len(sys.argv) > 1 else 4
t0 = time.time(); Lam, a, r0, x0 = Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)
D = X.mode_coeffs_exact(Lam, a, r0, x0, K); Dm = F.mode_coeffs(Lam, a, r0, x0, K)
Cm, Gm, tu, tx = F.assemble(Dm, K); ent, nU, nE, tot = X.assemble_exact(D, K, None)
print('assembled', nU, nE, round(time.time() - t0), 's', flush=True)
rC = F.frank(Cm); low = [n for n, d in enumerate(tot) if d <= J]; high = [n for n, d in enumerate(tot) if d > J]
rGL = F.frank(Gm[:, low]); print('mod-P ranks: C', rC, 'G_L', rGL, round(time.time() - t0), 's', flush=True)
hi = {u: n for n, u in enumerate(high)}; nr = len(high)
rows = [dict() for _ in range(2*nr)]
for (u, e), (x, y) in ent.items():
    if u not in hi: continue
    k = hi[u]
    if x: rows[k][e] = x; rows[nr + k][nE + e] = x
    if y: rows[k][nE + e] = y; rows[nr + k][e] = -y
M = flint.fmpz_mat(2*nr, 2*nE)
for i, rw in enumerate(rows):
    if not rw: continue
    L = 1
    for v in rw.values(): L = lcm(L, v.denominator)
    for j, v in rw.items(): M[i, j] = int(v*L)
del rows
print('integer C_H built', 2*nr, 'x', 2*nE, 'maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024, round(time.time() - t0), 's', flush=True)
rCH = M.rank()//2
print(f'J={J}: n_L={len(low)} rank_P C={rC} rank_Q C_H={rCH} rank_P G_L={rGL} => d_J <= {len(low) - rC + rCH - rGL} (RIGOROUS)  '
      f'maxrss MB {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024}  [{round(time.time() - t0)} s]', flush=True)
