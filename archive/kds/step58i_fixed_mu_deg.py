"""Exact sigma-degree of each selected minor at FIXED integer mu0:  deg det A(., mu0) = 2r - ord_tau det(tau^2 A(1/tau)),
and ord_tau >= K r - rank T_K (equality once K >= max partial multiplicity).  Upper bound is exact (Q-rank).
argv: point index, K, list of mu0."""
import sys, random, time, flint
from fractions import Fraction as Fr
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
ip, KT = int(sys.argv[1]), int(sys.argv[2]); mus = [int(x) for x in sys.argv[3].split(',')]
which = sys.argv[4].split(',') if len(sys.argv) > 4 else ['C', 'G1', 'G2']
pt = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))][ip]
t0 = time.time()
bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
X = {0: [('00', 0), ('01', 1), ('02', 2)], 1: [('10', 0), ('11', 1)], 2: [('20', 0)]}
for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
    if name not in which: continue
    r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
    n = KT*r
    for mu0 in mus:
        y = Fr(mu0)
        Bl = [[[tuple(sum((Q[key][i][j][p]*y**e for key, e in X[d]), Fr(0)) for p in (0, 1)) for j in cc] for i in pr] for d in (2, 1, 0)]
        ent = []
        for a in range(2*n):
            for b in range(2*n):
                aa, bb = a % n, b % n; bi, bj = aa // r, bb // r
                if bj > bi or bi - bj > 2: x_, y_ = 0, 0
                else: x_, y_ = Bl[bi - bj][aa % r][bb % r]
                v = Fr(x_ if (a < n) == (b < n) else (-y_ if a < n else y_)); ent.append(flint.fmpq(v.numerator, v.denominator))
        rk = flint.fmpq_mat(2*n, 2*n, ent).rank()//2
        print(f'pt{ip} {name} mu0={mu0} K={KT}: deficiency {n-rk} => deg_sigma <= {2*r-(n-rk)}  [{time.time()-t0:.0f}s]', flush=True)
