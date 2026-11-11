"""Does gauge_high + (top-degree principal-symbol kernel) account for all of S ∩ ker pi_J?
K = 6, J = 3, point (3, 3/5), Lam = 1/20, a = 3/5.  Exact nullity of C_H known: n_H - 708 (item 82).
 - top kernel: unknowns of total degree exactly K; rank of their rows of C computed EXACTLY (small); the kernel vectors are
   truncated solutions supported in degree K > J.
 - combined rank of [gauge rows X with deg >= J+2 restricted to high cols ; top-kernel basis] mod P (a LOWER bound)."""
import time, numpy as np, flint
from fractions import Fraction as Fr
import step67_4d_from_modes as F, step68_4d_exact as X
K, J = 6, 3; t0 = time.time(); Lam, a, r0, x0 = Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)
Dm = F.mode_coeffs(Lam, a, r0, x0, K); Cm, Gm, tu, tx = F.assemble(Dm, K)
D = X.mode_coeffs_exact(Lam, a, r0, x0, K); ent, nU, nE, tot = X.assemble_exact(D, K, None)
P = F.P
high = [n for n, d in enumerate(tu) if d > J]; top = [n for n, d in enumerate(tu) if d == K]
rtopQ = X.realified(ent, top, nE).rank()//2
print(f'n_H = {len(high)}, exact nullity of C_H = {len(high) - 708}; top-degree unknowns {len(top)}, exact rank C_top = {rtopQ}, '
      f'top kernel dim = {len(top) - rtopQ}  [{time.time()-t0:.0f}s]', flush=True)
# mod-P top kernel basis (left null vectors of C_top), embedded in the high coordinates
A = flint.nmod_mat(Cm.shape[1], len(top), [int(v) for v in Cm[top, :].T.flat], P)
nsp, nul = A.nullspace(); print('mod-P top kernel dim', nul, flush=True)
hpos = {u: k for k, u in enumerate(high)}
TK = np.zeros((nul, len(high)), dtype=np.int64)
for c in range(nul):
    for k, u in enumerate(top): TK[c, hpos[u]] = int(nsp[k, c])
xh = [n for n, d in enumerate(tx) if d >= J + 2]
GH = Gm[np.ix_(xh, high)]
r_g = F.frank(GH); r_comb = F.frank(np.vstack([GH, TK]))
print(f'rank gauge_high = {r_g}; rank [gauge_high ; top kernel] = {r_comb}; target nullity = {len(high) - 708}  [{time.time()-t0:.0f}s]')
print('support degree s: nullity of C restricted to unknowns of degree >= s (exact)  vs  gauge supported there (mod-P rank)')
for s in range(J + 1, K + 1):
    rows = [n for n, d in enumerate(tu) if d >= s]
    nullQ = len(rows) - X.realified(ent, rows, nE).rank()//2
    xs = [n for n, d in enumerate(tx) if d >= s + 1]
    rg = F.frank(Gm[np.ix_(xs, rows)])
    print(f'  s={s}: unknowns {len(rows)}, exact nullity {nullQ}, gauge rank {rg}, excess {nullQ - rg}  [{time.time()-t0:.0f}s]', flush=True)
