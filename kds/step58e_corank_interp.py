"""Caveat (B), exact and fast.  For the selected minor A (size r) of C, G_L(J=1), G_L(J=2):
rank A(0, mu0) <= k at 2r+1 distinct rationals mu0  ==>  every (k+1)-minor of A(0, mu) (deg_mu <= 2r) vanishes
identically  ==>  rank over Q(i)(mu) <= k  ==>  sigma^(r-k) | det A(sigma, mu) in Q(i)[sigma, mu].
Exact ranks over Q via the realification [[Re, -Im], [Im, Re]] (rank = 2 * complex rank), python-flint fmpq_mat."""
import sys, random, time, flint
from fractions import Fraction as Fr
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
pts = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))]
for pt in pts:
    t0 = time.time()
    bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
    low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
    Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
    random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
    tot = 0
    for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
        r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
        kmax = -1
        for t in range(2*r + 1):
            m0 = Fr(t + 1, 1)
            vals = [[(Q['00'][i][j][0] + Q['01'][i][j][0]*m0 + Q['02'][i][j][0]*m0*m0,
                      Q['00'][i][j][1] + Q['01'][i][j][1]*m0 + Q['02'][i][j][1]*m0*m0) for j in cc] for i in pr]
            ent = []
            for a in range(2*r):
                for b in range(2*r):
                    x, y = vals[a % r][b % r]
                    v = x if (a < r) == (b < r) else (-y if a < r else y)
                    ent.append(flint.fmpq(v.numerator, v.denominator))
            rk2 = flint.fmpq_mat(2*r, 2*r, ent).rank()
            assert rk2 % 2 == 0
            kmax = max(kmax, rk2//2)
        tot += r - kmax
        print(f'  {name}: size {r}, max rank at sigma=0 over {2*r+1} integer mu = {kmax} => sigma^{r-kmax} | det exactly', flush=True)
    print(f'pt {pt}: sigma^{tot} | E exactly; need 6: {tot >= 6}  [{time.time()-t0:.0f}s]', flush=True)
