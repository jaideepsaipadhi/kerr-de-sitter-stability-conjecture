"""Exact lower bound on mult_(0,0) det A for each selected minor: along the line (sigma, mu) = (alpha t, t),
ord_t det >= K r - rank T_K(A_0, A_1(alpha), A_2(alpha)); minors of T_K are polynomials in alpha of degree <= 2Kr, so the rank
bound at 2Kr+1 integer alpha holds for every line (homogeneity covers beta = 0).  mult_0 = min over lines of ord_t."""
import sys, random, time, flint
from fractions import Fraction as Fr
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
ip, KT = int(sys.argv[1]), int(sys.argv[2]); which = sys.argv[3].split(',') if len(sys.argv) > 3 else ['C', 'G1', 'G2']
pt = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))][ip]
t0 = time.time()
bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
tot = 0
for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
    if name not in which: continue
    r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
    n = KT*r; kmax = -1
    for t in range(2*n + 1):
        al = Fr(t + 1)
        X = {0: [('00', 1)], 1: [('10', al), ('01', 1)], 2: [('20', al*al), ('11', al), ('02', 1)]}
        Bl = [[[tuple(sum((Q[key][i][j][p]*w for key, w in X[d]), Fr(0)) for p in (0, 1)) for j in cc] for i in pr] for d in (0, 1, 2)]
        ent = []
        for a in range(2*n):
            for b in range(2*n):
                aa, bb = a % n, b % n; bi, bj = aa // r, bb // r
                if bj > bi or bi - bj > 2: x_, y_ = 0, 0
                else: x_, y_ = Bl[bi - bj][aa % r][bb % r]
                v = Fr(x_ if (a < n) == (b < n) else (-y_ if a < n else y_)); ent.append(flint.fmpq(v.numerator, v.denominator))
        kmax = max(kmax, flint.fmpq_mat(2*n, 2*n, ent).rank()//2)
    tot += n - kmax
    print(f'pt{ip} {name}: T_{KT} deficiency on every line = {n - kmax} => mult_0 >= {n - kmax}  [{time.time()-t0:.0f}s]', flush=True)
print(f'pt{ip}: mult_0 E >= {tot} (need 21)', flush=True)
