"""Generic (sigma, m) certificate at jet order J = 3 (= ord K1), K = 6, three test points (mod P).
Matrices are polynomials of total degree <= 2 in (sigma, mu = m): coefficients from 6 exact builds at nodes EV
(never sigma = 0 or m = 0, where shapes change).  E_i = product of the J = 3 minors (constraint, low-gauge, gauge-row
R_H minor).  E_i(sigma, mu) interpolated on a G x G grid (top rows/cols checked zero).  Then strip sigma^a mu^b and
compute gcd of the pairwise resultants Res_sigma(H_i, H_j)."""
import time, random, numpy as np, flint
from fractions import Fraction as Fr
import step52_exact_modes as S
from step55_exceptional import P, gauss_mod, rank_piv, interp

EV = [(1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (2, 2)]
MON = lambda s, u: [1, s, u, s*s, s*u, u*u]
V = [[Fr(c) for c in MON(s, u)] for s, u in EV]
def inv6(A):
    n = 6; M = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x/pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f*y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]
VI = inv6(V)
def fdet(M):
    n = M.shape[0]; return int(flint.nmod_mat(n, n, [int(v) for v in M.flat], P).det())

def build(Lam, a, r0, x0, K=6):
    bs = [S.run(Lam, a, Fr(1), r0, x0, (Fr(s), Fr(0)), u, K) for s, u in EV]
    def coeffs(kr, ki):
        A = bs[0][kr]; nr, nc = A.nrows(), A.ncols()
        Cm = [np.zeros((nr, nc), dtype=np.int64) for _ in range(6)]
        for i in range(nr):
            for j in range(nc):
                vr = [Fr(int(b[kr][i, j].p), int(b[kr][i, j].q)) for b in bs]
                vi = [Fr(int(b[ki][i, j].p), int(b[ki][i, j].q)) for b in bs]
                if not any(vr) and not any(vi): continue
                for k in range(6):
                    cr = sum(VI[k][q]*vr[q] for q in range(6)); ci = sum(VI[k][q]*vi[q] for q in range(6))
                    if cr or ci: Cm[k][i, j] = gauss_mod(cr, ci)
        return Cm
    b0 = bs[0]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
    low3 = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= 3]
    high3 = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > 3]
    Mb = len(mons); xr = [f*Mb + i for f in range(4) for i in range(Mb) if sum(mons[i]) >= 5]
    return dict(C=coeffs('CR', 'CI'), G=coeffs('GR', 'GI'), low=low3, high=high3, xr=xr)

def ev(Cm, s, u, rows=None, cols=None):
    mon = [m_ % P for m_ in MON(s, u)]
    M = sum(Cm[k]*mon[k] % P for k in range(6)) % P
    if rows is not None: M = M[rows]
    if cols is not None: M = M[:, cols]
    return M

def point(Lam, a, r0, x0, G, t0):
    B = build(Lam, a, r0, x0); print(f'  build ({r0}, {x0}) [{time.time()-t0:.0f}s]', flush=True)
    random.seed(9); s0, u0 = random.randrange(P), random.randrange(P)
    minors = []
    for Cm, rows, cols in ((B['C'], None, None), (B['G'], None, B['low']), (B['G'], B['xr'], B['high'])):
        r, pr, pc = rank_piv(ev(Cm, s0, u0, rows, cols))
        rr = [rows[i] for i in pr] if rows is not None else pr; cc = [cols[j] for j in pc] if cols is not None else pc
        minors.append((Cm, rr, cc))
    nL, nH = len(B['low']), len(B['high'])
    rC, rG, rX = (len(mm[1]) for mm in minors); d3 = nL - rC + (nH - rX) - rG
    E = lambda s, u: int(np.prod([fdet(ev(Cm, s, u, rr, cc)) for Cm, rr, cc in minors], dtype=object) % P)
    xs = list(range(1, G + 1))
    rows_ = [(interp(xs, [E(s, u) for u in xs]) + [0]*G)[:G] for s in xs]
    cols_ = [(interp(xs, [rows_[i][l] for i in range(G)]) + [0]*G)[:G] for l in range(G)]
    C = [[cols_[l][k] for l in range(G)] for k in range(G)]
    ks = [k for k in range(G) if any(C[k])]; ls = [l for l in range(G) if any(C[k][l] for k in range(G))]
    os_ = min(ks); om = min(ls)
    print(f'  pt ({r0},{x0}): generic d_3 bound {d3}; deg_sigma {max(ks)}, deg_mu {max(ls)} (grid {G}); '
          f'sigma-order {os_}, mu-order {om}  [{time.time()-t0:.0f}s]', flush=True)
    H = [[C[k][l] for l in range(om, max(ls) + 1)] for k in range(os_, max(ks) + 1)]
    return d3, H, max(ks) < G - 3 and max(ls) < G - 3

def resultant(F1, F2):
    n1, n2 = len(F1) - 1, len(F2) - 1
    l1 = max(l for r in F1 for l, c in enumerate(r) if c); l2 = max(l for r in F2 for l, c in enumerate(r) if c); DR = n2*l1 + n1*l2
    def sylv(u):
        a = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F1]; b = [sum(c*pow(u, l, P) for l, c in enumerate(r)) % P for r in F2]
        n = n1 + n2; M = np.zeros((n, n), dtype=np.int64)
        for k in range(n2):
            for j, c in enumerate(reversed(a)): M[k, k + j] = c
        for k in range(n1):
            for j, c in enumerate(reversed(b)): M[n2 + k, k + j] = c
        return fdet(M)
    us = list(range(1, DR + 2)); return flint.nmod_poly(interp(us, [sylv(u) for u in us]), P)

if __name__ == '__main__':
    t0 = time.time(); G = 80; Lam, a = Fr(1, 20), Fr(3, 5)
    out = [point(Lam, a, r0, x0, G, t0) for r0, x0 in ((Fr(3), Fr(3, 5)), (Fr(7, 2), Fr(-1, 5)), (Fr(4), Fr(1, 3)))]
    print('generic d_3 <= 0 at all points:', all(o[0] <= 0 for o in out), ' grids honest:', all(o[2] for o in out), flush=True)
    Hs = [o[1] for o in out]
    Rs = [resultant(Hs[i], Hs[j]) for i, j in ((0, 1), (0, 2), (1, 2))]
    for f in Rs:
        c = f.coeffs(); print('  R: degree', f.degree(), 'order at 0', next(k for k, x in enumerate(c) if int(x)), flush=True)
    g = Rs[0].gcd(Rs[1]).gcd(Rs[2])
    print(f'gcd of the three resultants mod P: degree {g.degree()}, monomial: {all(int(c) == 0 for c in g.coeffs()[:-1])}  [{time.time()-t0:.0f}s]')
