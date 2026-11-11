"""Upper bound for rank C_H without symbolic algebra: rows of the gauge map for X of degree >= J+2 lie in the left kernel of C_H.
R_H <= n_H - rank(those rows on the high columns).  Compare with the exact R_H (66 for J=0, 60 for J=1) at m=1."""
from fractions import Fraction as Fr
import flint
import step52_exact_modes as S
res = S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), (Fr(1, 3), Fr(1, 5)), 1, 4)
mons = res['mons']; nk, comps, monsK = res['nk'], res['comps'], res['monsK']
M = len(mons)
for J in (1, 2):
    high = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J]
    xrows = [f*M + i for f in range(4) for i in range(M) if sum(mons[i]) >= J + 2]
    GR = S.sub_cols(S.sub_rows(res['GR'], xrows), high); GI = S.sub_cols(S.sub_rows(res['GI'], xrows), high)
    rk = S.crank(GR, GI)
    CRH = S.sub_rows(res['CR'], high); CIH = S.sub_rows(res['CI'], high)
    print(f'J={J}: n_H = {len(high)}, rank(high-degree gauge rows) = {rk}, bound R_H <= {len(high) - rk}; exact rank C_H = {S.crank(CRH, CIH)}')
