import random, flint
from fractions import Fraction as Fr
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
pt = (Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5))
bs = build(*pt); Cq = coeff_mats(bs, 'CR', 'CI'); Cm = to_mod(Cq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
r, pr, pc = rank_piv(eval_mod(Cm, s0, u0))
al = 7
A = {d: [[0]*r for _ in range(r)] for d in (0, 1, 2)}
for a, i in enumerate(pr):
    for b, j in enumerate(pc):
        g = lambda k: int(Cm[k][i, j])
        A[0][a][b] = g('00'); A[1][a][b] = (al*g('10') + g('01')) % P; A[2][a][b] = (al*al*g('20') + al*g('11') + g('02')) % P
for K in range(2, 13):
    n = K*r
    M = flint.nmod_mat(n, n, [A[(x//r) - (y//r)][x % r][y % r] if 0 <= (x//r) - (y//r) <= 2 else 0 for x in range(n) for y in range(n)], P)
    print(f'K={K}: mod-p deficiency {n - M.rank()}', flush=True)
