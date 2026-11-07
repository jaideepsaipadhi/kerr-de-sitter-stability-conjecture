"""Caveat (B) via coranks: ord_sigma det A(sigma, mu) >= corank of A(0, mu) over Q(i)(mu) (sigma monic => divisibility in
Q(i)[sigma, mu]).  Compute exact ranks of the three selected minors at sigma = 0 over QQ_I(mu) (sympy DomainMatrix)."""
import sys, random, time, sympy as sp
from fractions import Fraction as Fr
from sympy.polys.matrices import DomainMatrix
from step56_m_uniform import build, coeff_mats, to_mod, eval_mod, EV
from step55_exceptional import P, rank_piv
pt = [(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)), (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5)), (Fr(1, 20), Fr(3, 5), Fr(4), Fr(1, 3))][int(sys.argv[1])]
t0 = time.time()
bs = build(*pt); b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI'); Cm, Gm = to_mod(Cq), to_mod(Gq)
random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
u = sp.Symbol('u'); dom = sp.QQ_I[u]
q = lambda c: sp.Rational(c[0].numerator, c[0].denominator) + sp.I*sp.Rational(c[1].numerator, c[1].denominator)
tot = 0
for name, Q, Qm, cols in (('C', Cq, Cm, None), ('G1', Gq, Gm, low[1]), ('G2', Gq, Gm, low[2])):
    r, pr, pc = rank_piv(eval_mod(Qm, s0, u0, cols=cols)); cc = [cols[c] for c in pc] if cols else pc
    data = [[dom.convert(sp.expand(q(Q['00'][i][j]) + q(Q['01'][i][j])*u + q(Q['02'][i][j])*u**2)) for j in cc] for i in pr]
    rk = DomainMatrix(data, (r, r), dom).convert_to(dom.get_field()).rank()
    tot += r - rk
    print(f'{name}: size {r}, exact rank at sigma=0 over Q(i)(mu) = {rk}, corank {r - rk}  [{time.time()-t0:.0f}s]', flush=True)
print(f'pt {pt}: sigma^{tot} | E exactly (need 6): {tot >= 6}')
