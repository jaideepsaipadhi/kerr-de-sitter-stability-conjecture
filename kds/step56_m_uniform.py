"""
B1 uniformity in (sigma, m) jointly.  mu = m as a complex variable; C, G are polynomials of total degree <= 2 in (sigma, mu),
built EXACTLY from 6 exact evaluations.  Two points p1, p2.
 (1) generic count bound over Q(i)(sigma, mu): n_L - r_C + R_H - r_G <= 0  (R_H exact, sympy over QQ_I(s,u)).
 (2) E_k(sigma, mu) = product of fixed maximal minors (C, G_L for J=0,1) at p_k, mod P.
 (3) Deflate by sigma^v; R(mu) = Res_sigma(E1', E2') mod P; leading coefficients lc_k(mu).
     If R(mu) != 0 and R and lc1*lc2 have NO roots in F_P, then no integer m admits a common exceptional sigma != 0.
 Output also the sigma = 0 line data for a separate check.
"""
import sys, time, random
from fractions import Fraction as Fr
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
import step52_exact_modes as S
from step55_exceptional import P, I_MOD, fr_mod, gauss_mod, rank_piv, det_mod, interp, pgcd, pmul

EV = [(1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (2, 2)]
def build(Lam, a, r0, x0, K=4):
    bs = {ev: S.run(Lam, a, Fr(1), r0, x0, (Fr(ev[0]), Fr(0)), Fr(ev[1]), K) for ev in EV}
    return bs
MONS = [lambda s, u: 1, lambda s, u: s, lambda s, u: s*s, lambda s, u: u, lambda s, u: u*u, lambda s, u: s*u]
KEYS = ['00', '10', '20', '01', '02', '11']
VM = sp.Matrix([[sp.Rational(f(*ev)) for f in MONS] for ev in EV]); VMI = VM.inv()
VMI = [[Fr(int(sp.fraction(VMI[i, j])[0]), int(sp.fraction(VMI[i, j])[1])) for j in range(6)] for i in range(6)]

def coeff_mats(bs, kr, ki, rows=None, cols=None):
    """Exact coefficients of p = c00 + c10 s + c20 s^2 + c01 u + c02 u^2 + c11 s u for each entry (Gaussian rationals)."""
    A = bs[EV[0]][kr]; nr, nc = A.nrows(), A.ncols()
    rows = list(range(nr)) if rows is None else rows; cols = list(range(nc)) if cols is None else cols
    out = {k: [[None]*len(cols) for _ in rows] for k in ('00', '10', '20', '01', '02', '11')}
    for ii, i in enumerate(rows):
        for jj, j in enumerate(cols):
            v = {ev: (Fr(int(bs[ev][kr][i, j].p), int(bs[ev][kr][i, j].q)), Fr(int(bs[ev][ki][i, j].p), int(bs[ev][ki][i, j].q))) for ev in EV}
            vals = [v[ev] for ev in EV]
            for kk, k in enumerate(KEYS):
                c = (sum(VMI[kk][t]*vals[t][0] for t in range(6)), sum(VMI[kk][t]*vals[t][1] for t in range(6)))
                out[k][ii][jj] = c
    return out

def to_mod(Cq):
    return {k: np.array([[gauss_mod(*c) for c in row] for row in M], dtype=np.int64) for k, M in Cq.items()}

def eval_mod(Cm, s, u, rows=None, cols=None):
    M = (Cm['00'] + Cm['10']*s + Cm['20']*(s*s % P) + Cm['01']*u + Cm['02']*(u*u % P) + Cm['11']*(s*u % P)) % P
    if rows is not None: M = M[rows]
    if cols is not None: M = M[:, cols]
    return M

def sym_rank2(Cq, rows):
    s_, u_ = sp.symbols('s u'); dom = sp.QQ_I[s_, u_]
    mon = {'00': 1, '10': s_, '20': s_**2, '01': u_, '02': u_**2, '11': s_*u_}
    data = []
    for i in rows:
        row = []
        for j in range(len(Cq['00'][0])):
            e = 0
            for k, mk in mon.items():
                c = Cq[k][i][j]
                if c[0] or c[1]: e += (sp.Rational(c[0].numerator, c[0].denominator) + sp.I*sp.Rational(c[1].numerator, c[1].denominator))*mk
            row.append(dom.convert(sp.expand(e)))
        data.append(row)
    return DomainMatrix(data, (len(data), len(data[0])), dom).convert_to(dom.get_field()).rank()

def point_data(Lam, a, r0, x0, t0):
    bs = build(Lam, a, r0, x0)
    b0 = bs[EV[0]]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
    low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (1, 2)}
    high = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J] for J in (1, 2)}
    Mb = len(mons)
    Cq = coeff_mats(bs, 'CR', 'CI'); Gq = coeff_mats(bs, 'GR', 'GI')
    Cm, Gm = to_mod(Cq), to_mod(Gq)
    random.seed(3); s0, u0 = random.randrange(P), random.randrange(P)
    rC, prC, pcC = rank_piv(eval_mod(Cm, s0, u0))
    minors = [(Cm, prC, pcC)]; info = {'rC': rC}
    for J in (1, 2):
        rG, prG, pcG = rank_piv(eval_mod(Gm, s0, u0, cols=low[J]))
        minors.append((Gm, prG, [low[J][c] for c in pcG]))
        xrows = [f*Mb + i for f in range(4) for i in range(Mb) if sum(mons[i]) >= J + 2]
        RH = len(high[J]) - rank_piv(eval_mod(Gm, s0, u0, rows=xrows, cols=high[J]))[0]   # valid UPPER bound for rank C_H
        info[J] = (rG, RH, len(low[J]) - rC + RH - rG)
    print(f'  point r={r0} x={x0}: rank C = {rC};  J -> (rank G_L, R_H upper bound, d bound): '
          f'{info[1]}, {info[2]}   [{time.time()-t0:.0f}s]', flush=True)
    return minors, info

def E_eval(minors, s, u):
    v = 1
    for (Mm, rr, cc) in minors: v = v*det_mod(eval_mod(Mm, s, u, rows=rr, cols=cc)) % P
    return v

def univ_in_sigma(minors, u, D):
    xs = list(range(1, D + 2)); return interp(xs, [E_eval(minors, x, u) for x in xs])

def deflate(p):
    k = 0
    while k < len(p) and p[k] == 0: k += 1
    return p[k:], k

def resultant(a, b):
    """Resultant mod P via Sylvester determinant (a, b coefficient lists low->high)."""
    da, db = len(a) - 1, len(b) - 1
    n = da + db; M = np.zeros((n, n), dtype=np.int64)
    for i in range(db):
        for j, c in enumerate(reversed(a)): M[i, i + j] = c
    for i in range(da):
        for j, c in enumerate(reversed(b)): M[db + i, i + j] = c
    return det_mod(M)

def has_root_in_Fp(poly):
    """gcd(poly, x^P - x) != 1 ?  poly low->high mod P."""
    poly = [c % P for c in poly]
    while len(poly) > 1 and poly[-1] == 0: poly.pop()
    if len(poly) <= 1: return poly[0] == 0
    def mulmod(a, b):
        r = pmul(a, b)
        # reduce mod poly
        inv = pow(poly[-1], P - 2, P); d = len(poly) - 1
        while len(r) - 1 >= d:
            f = r[-1]*inv % P; sh = len(r) - 1 - d
            for k in range(len(poly)): r[sh + k] = (r[sh + k] - f*poly[k]) % P
            while len(r) > 1 and r[-1] == 0: r.pop()
            if len(r) - 1 < d: break
        return r
    # x^P mod poly by square-and-multiply
    result = [1]; base = [0, 1]; e = P
    while e:
        if e & 1: result = mulmod(result, base)
        base = mulmod(base, base); e >>= 1
    xp_minus_x = result + [0]*max(0, 2 - len(result)); xp_minus_x[1] = (xp_minus_x[1] - 1) % P
    g = pgcd(poly[:], xp_minus_x)
    return len(g) - 1 >= 1

if __name__ == '__main__':
    t0 = time.time()
    p1 = (Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)); p2 = (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5))
    M1, I1 = point_data(*p1, t0); M2, I2 = point_data(*p2, t0)
    okd = all(I[J][2] <= 0 for I in (I1, I2) for J in (1, 2))
    print('generic (sigma, mu) d bounds <= 0 at both points:', okd, flush=True)
    # degrees in sigma at a random mu, and in mu at a random sigma
    Ds = 2*(I1['rC'] + I1[1][0] + I1[2][0])            # safe bound for sigma-degree of E (entries degree <= 2)
    random.seed(9); ur = random.randrange(P)
    e1 = univ_in_sigma(M1, ur, Ds); e2 = univ_in_sigma(M2, ur, Ds)
    e1d, v1 = deflate(e1); e2d, v2 = deflate(e2)
    print(f'deg_sigma E1 = {len(e1)-1} (sigma^{v1} factor), deg_sigma E2 = {len(e2)-1} (sigma^{v2}) at random mu', flush=True)
    d1, d2 = len(e1d) - 1, len(e2d) - 1
    ds1, ds2 = len(e1) - 1, len(e2) - 1
    # mu-degrees at a random sigma
    sr = random.randrange(P); Du = Ds
    f1 = interp(list(range(1, Du + 2)), [E_eval(M1, sr, x) for x in range(1, Du + 2)])
    f2 = interp(list(range(1, Du + 2)), [E_eval(M2, sr, x) for x in range(1, Du + 2)])
    du1, du2 = len(f1) - 1, len(f2) - 1
    print(f'deg_mu E1 = {du1}, deg_mu E2 = {du2}', flush=True)
    DR = d1*du2 + d2*du1 + 2
    def Rat(u):
        a = univ_in_sigma(M1, u, ds1)[v1:v1 + d1 + 1]; b = univ_in_sigma(M2, u, ds2)[v2:v2 + d2 + 1]
        a = (a + [0]*(d1 + 1))[:d1 + 1]; b = (b + [0]*(d2 + 1))[:d2 + 1]
        return resultant(a, b), a[-1], b[-1]
    us = list(range(1, DR + 2)); vals = [Rat(u) for u in us]
    R = interp(us, [v[0] for v in vals]); L1 = interp(us, [v[1] for v in vals]); L2 = interp(us, [v[2] for v in vals])
    def peval(p, x): return sum(c*pow(x, k, P) for k, c in enumerate(p)) % P
    chk = all(peval(R, x) == Rat(x)[0] for x in (DR + 7, DR + 19))
    print('R interpolant verified at 2 extra points:', chk, flush=True)
    print(f'R(mu) degree {len(R)-1} (bound {DR}); R identically zero: {R == [0]}  [{time.time()-t0:.0f}s]', flush=True)
    print('R has a root in F_P:', has_root_in_Fp(R), ' lc1 root in F_P:', has_root_in_Fp(L1), ' lc2 root in F_P:', has_root_in_Fp(L2))
