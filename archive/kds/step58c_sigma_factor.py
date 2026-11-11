"""Caveat (B): is sigma^6 | E_i structural?  For each minor, count rows/cols whose sigma^0 part (keys '00','01','02')
vanishes identically EXACTLY (then sigma divides that row/col as a polynomial), factor them out, and redo the exact
tropical lower bound on the sigma-order of the remaining determinant."""
import random, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linear_sum_assignment
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
pt = (Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5))
bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
nz = lambda c: c[0] != 0 or c[1] != 0
def sord(Q, i, j):   # exact sigma-order of entry (i, j) as polynomial in (sigma, mu)
    for k, keys in ((0, ('00', '01', '02')), (1, ('10', '11')), (2, ('20',))):
        if any(nz(Q[key][i][j]) for key in keys): return k
    return None
tot = 0
for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
    r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
    V = np.full((r, r), 99.0)
    for a, i in enumerate(pr):
        for b, j in enumerate(cc):
            o = sord(Q, i, j)
            if o is not None: V[a, b] = o
    rowmin = V.min(axis=1); V2 = V - rowmin[:, None]; colmin = V2.min(axis=0); V3 = V2 - colmin[None, :]
    ri, ci = linear_sum_assignment(V3); rest = V3[ri, ci].sum()
    lb = rowmin.sum() + colmin.sum() + rest; tot += lb
    print(f'{name}: size {r}; sigma pulled from rows {rowmin.sum():.0f}, from cols {colmin.sum():.0f}, tropical rest {rest:.0f} '
          f'=> exact sigma-order >= {lb:.0f}')
print(f'exact sigma-order of E >= {tot:.0f}  (mod-p order is 6; need >= 6)')
