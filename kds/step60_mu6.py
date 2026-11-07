"""Exact DEGREE upper bounds for each selected minor of A(sigma, mu) (entries of total degree <= 2), via the
Toeplitz deficiency of the reversed matrix polynomial:  deg_x det A <= 2r - (K r - rank T_K(reversed in x)),
rank bound transferred to all values of the other variable by vanishing of minors at 2*K*r+1 integer points.
mode 'sig': degree in sigma (uniform in mu);  mode 'mu': degree in mu (uniform in sigma)."""
import sys, random, time, flint
from fractions import Fraction as Fr
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
mode, KT, ip = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); FWD = len(sys.argv) > 4
pt = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))][ip]
t0 = time.time()
bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
# coefficient of x^d as polynomial in the other variable y: list of (key, power of y)
if mode == 'sig': X = {0: [('00', 0), ('01', 1), ('02', 2)], 1: [('10', 0), ('11', 1)], 2: [('20', 0)]}
else:             X = {0: [('00', 0), ('10', 1), ('20', 2)], 1: [('01', 0), ('11', 1)], 2: [('02', 0)]}
tot = 0
for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
    r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
    n = KT*r; kmax = -1
    for t in range(2*n + 1):
        y = Fr(t + 1)
        Bl = []
        for d in ((0, 1, 2) if FWD else (2, 1, 0)):
            Bl.append([[tuple(sum((Q[key][i][j][p]*y**e for key, e in X[d]), Fr(0)) for p in (0, 1)) for j in cc] for i in pr])
        ent = []
        for a in range(2*n):
            for b in range(2*n):
                aa, bb = a % n, b % n; bi, bj = aa // r, bb // r
                if bj > bi or bi - bj > 2: x_, y_ = 0, 0
                else: x_, y_ = Bl[bi - bj][aa % r][bb % r]
                v = Fr(x_ if (a < n) == (b < n) else (-y_ if a < n else y_)); ent.append(flint.fmpq(v.numerator, v.denominator))
        rk2 = flint.fmpq_mat(2*n, 2*n, ent).rank(); kmax = max(kmax, rk2//2)
    ub = (n - kmax) if FWD else 2*r - (n - kmax); tot += ub
    print(f'  {name}: size {r}, reversed T_{KT} deficiency {n-kmax} => deg_{mode} det <= {ub}  [{time.time()-t0:.0f}s]', flush=True)
print(f'pt {ip}: deg_{mode} E <= {tot}  (mod-p value 29; grid 40 needs <= 39)', flush=True)
