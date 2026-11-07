"""
B1 sigma-uniformity for fixed m (K = 4, J in {0, 1}), two points p1, p2.
For each point: build C(sigma), G(sigma) exactly as degree-<=2 polynomials (3 exact builds), take generic ranks
  R_C, R_G (lower bounds from a random evaluation mod P; exact R_H from sympy over Q(i)(sigma)),
check the generic count d_J = n_L - R_C + R_H - R_G <= 0, and compute exceptional polynomials mod P:
  delta_C(sigma) = det of a fixed maximal nonsingular minor of C(sigma),  delta_G likewise for G_L(sigma).
Exceptional set at p (where rank C or rank G_L can drop) is contained in roots(delta_C delta_G).
If, as homogeneous forms of their degree bounds, E1 = delta_C(p1) delta_G0(p1) delta_G1(p1) and E2 (same at p2) have no common
projective root mod P, then over Qbar(i) no sigma is exceptional at both points: every sigma has a good point.
(Reduction mod a prime above P maps a common root over Qbar to a common projective root mod P.)
"""
import sys, time, random
from fractions import Fraction as Fr
import numpy as np, sympy as sp
from sympy.polys.matrices import DomainMatrix
import step52_exact_modes as S
P = 1000000009 if (1000000009 % 4 == 1) else None
# choose a prime p = 1 mod 4 with p < 2^30 (products fit int64 after mod)
P = 998244353   # = 119*2^23+1, p % 4 == 1
I_MOD = pow(3, (P - 1)//4, P)    # 3 is a primitive root mod 998244353; this is sqrt(-1)
assert I_MOD*I_MOD % P == P - 1
def fr_mod(q): return (q.numerator % P)*pow(q.denominator % P, P - 2, P) % P
def gauss_mod(re, im): return (fr_mod(re) + I_MOD*fr_mod(im)) % P

def build_point(Lam, a, r0, x0, m, K=4):
    pts = [Fr(0), Fr(1), Fr(2)]
    bs = [S.run(Lam, a, Fr(1), r0, x0, (sv, Fr(0)), m, K) for sv in pts]
    def coeffs(kr, ki, rows=None, cols=None):
        A = bs[0][kr]; nr, nc = A.nrows(), A.ncols()
        rows = list(range(nr)) if rows is None else rows; cols = list(range(nc)) if cols is None else cols
        Cq = [[[None]*len(cols) for _ in rows] for _ in range(3)]   # exact Gaussian-rational coefficients c0,c1,c2
        for ii, i in enumerate(rows):
            for jj, j in enumerate(cols):
                v = [(Fr(int(b[kr][i, j].p), int(b[kr][i, j].q)), Fr(int(b[ki][i, j].p), int(b[ki][i, j].q))) for b in bs]
                c0 = v[0]; c2 = ((v[2][0] - 2*v[1][0] + v[0][0])/2, (v[2][1] - 2*v[1][1] + v[0][1])/2)
                c1 = (v[1][0] - v[0][0] - c2[0], v[1][1] - v[0][1] - c2[1])
                Cq[0][ii][jj], Cq[1][ii][jj], Cq[2][ii][jj] = c0, c1, c2
        return Cq
    r0d = bs[0]; nk, comps, monsK, mons = r0d['nk'], r0d['comps'], r0d['monsK'], r0d['mons']
    lowJ = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in (0, 1)}
    highJ = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J] for J in (0, 1)}
    return dict(C=coeffs('CR', 'CI'), G=coeffs('GR', 'GI'), low=lowJ, high=highJ)

def mod_mat(Cq, sig, rows=None, cols=None):
    nr, nc = len(Cq[0]), len(Cq[0][0])
    rows = range(nr) if rows is None else rows; cols = range(nc) if cols is None else cols
    s2 = sig*sig % P
    M = np.zeros((len(rows), len(cols)), dtype=np.int64)
    for ii, i in enumerate(rows):
        for jj, j in enumerate(cols):
            v = 0
            for k, sk in ((0, 1), (1, sig), (2, s2)):
                c = Cq[k][i][j]
                if c[0] or c[1]: v = (v + gauss_mod(*c)*sk) % P
            M[ii, jj] = v
    return M

def rank_piv(M):
    A = M.copy() % P; nr, nc = A.shape; r = 0; prow = []; pcol = []
    rowsidx = list(range(nr))
    for c in range(nc):
        piv = next((i for i in range(r, nr) if A[i, c] % P), None)
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]; rowsidx[r], rowsidx[piv] = rowsidx[piv], rowsidx[r]
        inv = pow(int(A[r, c]), P - 2, P); A[r] = A[r]*inv % P
        for i in range(nr):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c]*A[r]) % P
        prow.append(rowsidx[r]); pcol.append(c); r += 1
        if r == nr: break
    return r, prow, pcol

def det_mod(M):
    A = M.copy() % P; n = A.shape[0]; d = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if A[i, c] % P), None)
        if piv is None: return 0
        if piv != c: A[[c, piv]] = A[[piv, c]]; d = -d
        d = d*int(A[c, c]) % P; inv = pow(int(A[c, c]), P - 2, P)
        for i in range(c + 1, n):
            if A[i, c]: A[i] = (A[i] - A[i, c]*inv % P*A[c]) % P
    return d % P

def interp(xs, ys):
    """Coefficients (low -> high) of the polynomial through (xs, ys) mod P (Newton form)."""
    n = len(xs); coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1])*pow((xs[i] - xs[i - j]) % P, P - 2, P) % P
    poly = [0]*n
    for i in range(n - 1, -1, -1):           # Horner expansion
        poly = [(-xs[i]*poly[0] + coef[i]) % P] + [(poly[k - 1] - xs[i]*poly[k]) % P for k in range(1, n)] if False else poly
    # straightforward expansion
    res = [0]
    for i in range(n - 1, -1, -1):
        res = [((res[k - 1] if k >= 1 else 0) - (xs[i]*res[k] if k < len(res) else 0)) % P for k in range(len(res) + 1)]
        res[0] = (res[0] + coef[i]) % P
    while len(res) > 1 and res[-1] == 0: res.pop()
    return res

def exceptional_poly(Cq, rows=None, cols=None):
    random.seed(1); s0 = random.randrange(P)
    r, prow, pcol = rank_piv(mod_mat(Cq, s0, rows, cols))
    rr = [rows[i] for i in prow] if rows is not None else prow
    cc = [cols[j] for j in pcol] if cols is not None else pcol
    D = 2*r
    xs = list(range(1, D + 2)); ys = [det_mod(mod_mat(Cq, x, rr, cc)) for x in xs]
    return r, interp(xs, ys), D

def pgcd(a, b):
    a = a[:]; b = b[:]
    def trim(p):
        while len(p) > 1 and p[-1] == 0: p.pop()
        return p
    a, b = trim(a), trim(b)
    while len(b) > 1 or b[0] != 0:
        if len(b) == 1 and b[0] == 0: break
        # a mod b
        inv = pow(b[-1], P - 2, P)
        while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
            f = a[-1]*inv % P; sh = len(a) - len(b)
            for k in range(len(b)): a[sh + k] = (a[sh + k] - f*b[k]) % P
            trim(a)
            if len(a) < len(b) or (len(a) == 1 and a[0] == 0): break
        a, b = b, a
        if len(b) == 1 and b[0] == 0: break
    return trim(a)

def pmul(a, b):
    out = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] = (out[i + j] + x*y) % P
    return out

def sym_rank(Cq, rows):
    sig = sp.Symbol('s'); dom = sp.QQ_I[sig]
    data = []
    for i in rows:
        row = []
        for j in range(len(Cq[0][0])):
            e = 0
            for k in range(3):
                c = Cq[k][i][j]
                if c[0] or c[1]: e += (sp.Rational(c[0].numerator, c[0].denominator) + sp.I*sp.Rational(c[1].numerator, c[1].denominator))*sig**k
            row.append(dom.convert(sp.expand(e)))
        data.append(row)
    return DomainMatrix(data, (len(data), len(data[0])), dom).convert_to(dom.get_field()).rank()

def analyse(Lam, a, r0, x0, m):
    t0 = time.time(); B = build_point(Lam, a, r0, x0, m)
    out = {}
    rC, dC, DC = exceptional_poly(B['C'])
    E = dC; Dtot = DC
    for J in (0, 1):
        rG, dG, DG = exceptional_poly(B['G'], cols=B['low'][J])
        E = pmul(E, dG); Dtot += DG
        RH = sym_rank(B['C'], B['high'][J]); nL = len(B['low'][J])
        out[J] = (rG, RH, nL - rC + RH - rG)
    print(f'  point (Lam={Lam}, a={a}, r={r0}, x={x0}) m={m}: rank C = {rC} (deg delta_C {len(dC)-1}/{DC}); '
          f'(rank G_L, generic rank C_H, d bound) J=0: {out[0]}, J=1: {out[1]}; {time.time()-t0:.0f}s', flush=True)
    return rC, out, E, Dtot

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    p1 = (Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5)); p2 = (Fr(1, 20), Fr(3, 5), Fr(7, 2), Fr(-1, 5))
    rC1, rG1, E1, D1 = analyse(*p1, m)
    rC2, rG2, E2, D2 = analyse(*p2, m)
    g = pgcd(E1, E2)
    both_inf = (len(E1) - 1 < D1) and (len(E2) - 1 < D2)
    print('gcd coefficients (low->high):', g[:6], '... len', len(g)); 
    print(f'm={m}: gcd(E1, E2) mod P has degree {len(g)-1};  common root at infinity possible: {both_inf}')
    okd = all(v[2] <= 0 for v in list(rG1.values()) + list(rG2.values()))
    print('generic d bounds <= 0 at both points:', okd)
    print('NO sigma is exceptional at both points (mod-P certificate):', len(g) - 1 == 0 and not both_inf)
    print('=> every sigma has a good point for this m:', okd and len(g) - 1 == 0 and not both_inf)
