"""Line certificates for Theorem B at jet order J <= 3 (= ord K1), K = 6.
Mode m fixed, sigma the variable (covers the line m = 0, and any fixed m).  At a point p:
  d_J(sigma) <= n_L(J) - rank C + R_H(J) - rank G_L(J),
  R_H(J) <= n_H(J) - rank(gauge rows of degree >= J+2 restricted to high columns)    [gauge rows are solutions: valid]
Each rank lower bound holds off the zeros of one fixed maximal minor (chosen at a random sigma mod P).  E(sigma) = product
of all these minors; it is computed mod P by interpolation (entries have degree <= 2 in sigma).
Output: the generic d bound for J = 0..3, deg/order of E, and gcd(E at p1, E at p2) mod P."""
import sys, time, random
import numpy as np, flint
from fractions import Fraction as Fr
import step52_exact_modes as S
from step55_exceptional import P, gauss_mod, rank_piv, interp, pgcd, pmul

def fdet(M):
    n = M.shape[0]; return int(flint.nmod_mat(n, n, [int(v) for v in M.flat], P).det())

def build(Lam, a, r0, x0, m, K):
    pts = [Fr(1), Fr(2), Fr(3)]   # avoid sigma = 0 (extra stationary gauge rows at m = 0)
    if LINE == 'sigma0':          # variable = m (mu), sigma = 0; nodes m = 1, 2, 3
        bs = [S.run(Lam, a, Fr(1), r0, x0, (Fr(0), Fr(0)), int(mv), K) for mv in pts]
    else:
        bs = [S.run(Lam, a, Fr(1), r0, x0, (sv, Fr(0)), m, K) for sv in pts]
    def coeffs(kr, ki):
        A = bs[0][kr]; nr, nc = A.nrows(), A.ncols()
        Cm = [np.zeros((nr, nc), dtype=np.int64) for _ in range(3)]
        for i in range(nr):
            for j in range(nc):
                v = [(Fr(int(b[kr][i, j].p), int(b[kr][i, j].q)), Fr(int(b[ki][i, j].p), int(b[ki][i, j].q))) for b in bs]
                c2 = ((v[2][0] - 2*v[1][0] + v[0][0])/2, (v[2][1] - 2*v[1][1] + v[0][1])/2)
                c1 = (v[1][0] - v[0][0] - 3*c2[0], v[1][1] - v[0][1] - 3*c2[1])
                c0 = (v[0][0] - c1[0] - c2[0], v[0][1] - c1[1] - c2[1])
                for k, c in enumerate((c0, c1, c2)):
                    if c[0] or c[1]: Cm[k][i, j] = gauss_mod(*c)
        return Cm
    b0 = bs[0]; nk, comps, monsK, mons = b0['nk'], b0['comps'], b0['monsK'], b0['mons']
    low = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J] for J in range(4)}
    high = {J: [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J] for J in range(4)}
    Mb = len(mons)
    xrows = {J: [f*Mb + i for f in range(4) for i in range(Mb) if sum(mons[i]) >= J + 2] for J in range(4)}
    return dict(C=coeffs('CR', 'CI'), G=coeffs('GR', 'GI'), low=low, high=high, xrows=xrows)

def ev(Cm, s, rows=None, cols=None):
    M = (Cm[0] + Cm[1]*s + Cm[2]*(s*s % P)) % P
    if rows is not None: M = M[rows]
    if cols is not None: M = M[:, cols]
    return M

def minor_poly(Cm, s0, rows=None, cols=None):
    r, pr, pc = rank_piv(ev(Cm, s0, rows, cols))
    rr = [rows[i] for i in pr] if rows is not None else pr
    cc = [cols[j] for j in pc] if cols is not None else pc
    D = 2*r; xs = list(range(1, D + 2))
    return r, interp(xs, [fdet(ev(Cm, x, rr, cc)) for x in xs])

def analyse(Lam, a, r0, x0, m, K=6):
    t0 = time.time(); B = build(Lam, a, r0, x0, m, K)
    random.seed(5); s0 = random.randrange(P)
    rC, E = minor_poly(B['C'], s0)
    info = {}
    for J in range(4):
        rG, eG = minor_poly(B['G'], s0, cols=B['low'][J])
        rX, eX = minor_poly(B['G'], s0, rows=B['xrows'][J], cols=B['high'][J])
        if J == 3: E = pmul(pmul(E, eG), eX)            # only J = N = 3 is needed (item 75 R1, item 76)
        RH = len(B['high'][J]) - rX
        info[J] = len(B['low'][J]) - rC + RH - rG
    while len(E) > 1 and E[-1] == 0: E.pop()
    o = next(k for k, c in enumerate(E) if c)
    print(f'  pt (r={r0}, x={x0}) m={m}: generic d bounds J=0..3: {[info[J] for J in range(4)]}; '
          f'deg E {len(E)-1}, ord_0 E {o}  [{time.time()-t0:.0f}s]', flush=True)
    return info, E

LINE = sys.argv[2] if len(sys.argv) > 2 else 'fixed_m'
if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    Lam, a = Fr(1, 20), Fr(3, 5)
    i1, E1 = analyse(Lam, a, Fr(3), Fr(3, 5), m)
    i2, E2 = analyse(Lam, a, Fr(7, 2), Fr(-1, 5), m)
    g = pgcd(E1, E2)
    mono = all(c == 0 for c in g[:-1])
    print(f'm={m}: gcd(E1, E2) mod P degree {len(g)-1}; is c*sigma^k: {mono}')
    print('generic d_3 bound <= 0 at both points:', i1[3] <= 0 and i2[3] <= 0)
    print('generic d bounds <= 0 at both points:', all(v <= 0 for v in list(i1.values()) + list(i2.values())))
