import sys
from fractions import Fraction as Fr
import flint, step52_exact_modes as S
K = 5
for sig, m in (((Fr(0), Fr(0)), 1), ((Fr(1, 3), Fr(1, 5)), 1)):
    res = S.run(Fr(1, 20), Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), sig, m, K)
    GR, GI, CR, CI = res['GR'], res['GI'], res['CR'], res['CI']
    PR = GR*CR - GI*CI; PI = GR*CI + GI*CR
    nz = sum(1 for i in range(PR.nrows()) for j in range(PR.ncols()) if PR[i, j] != 0 or PI[i, j] != 0)
    print(f'sigma={sig} m={m}: gauge rows violating exact constraints: {nz} nonzero entries of {PR.nrows()*PR.ncols()}')
