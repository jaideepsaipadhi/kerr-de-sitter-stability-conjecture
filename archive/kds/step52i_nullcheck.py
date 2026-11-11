import flint
from fractions import Fraction as Fr
import step52_exact_modes as S
res = S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), (Fr(0), Fr(0)), 1, 4)
R = S.realify(res['CR'], res['CI']); RT = R.transpose()
Rr, rk = RT.rref()
piv = []; row = 0
for j in range(RT.ncols()):
    if row < rk and Rr[row, j] != 0: piv.append(j); row += 1
print('rank', rk, ' pivot values (first 8):', [Rr[i, piv[i]] for i in range(min(8, len(piv)))])
free = [j for j in range(RT.ncols()) if j not in piv]
f = free[0]; v = flint.fmpq_mat(RT.ncols(), 1); v[f, 0] = 1
for i, pj in enumerate(piv): v[pj, 0] = -Rr[i, f]/Rr[i, pj]
w = RT*v
print('residual of first null vector (exact):', sum(1 for i in range(w.nrows()) if w[i, 0] != 0), 'nonzero entries')
