"""Exact-input 4D jet system assembled from the mode coefficient matrices (mod P), and a rigorous UPPER bound on d_J.
Mode operator C(sigma, m) = sum c_ab sigma^a m^b (exact builds at 6 nodes, step66).  With sigma <-> i d_t, m <-> -i d_phi:
   L(t^j phi^k H) = sum_ab c_ab i^a (-i)^b [j!/(j-a)!][k!/(k-b)!] t^(j-a) phi^(k-b) (C_ab H).
4D unknowns: t^j phi^k x (comp, (r,x)-monomial mu), j+k+deg mu <= K.  Equations: t^j' phi^k' x (block, nu), j'+k'+deg nu <= K-2.
Gauge generators: t^j phi^k x (X_f, mu), j+k+deg mu <= K+1, same rule with the gauge coefficient matrices.
Bound (all ranks mod P; each inequality is rigorous):
   d_J <= n_L - rank C + (n_H - rank Xhigh) - rank G_L,
   Xhigh = gauge generators of total degree >= J+2 (their images have no jets of degree <= J, so they lie in S ∩ ker pi_J)."""
import sys, time, itertools, numpy as np, flint
from fractions import Fraction as Fr
from math import factorial
import step52_exact_modes as S
from step55_exceptional import P, gauss_mod, rank_piv
from step66_biv_J3 import EV, VI

IU = next(pow(g, (P - 1)//4, P) for g in range(2, 100) if pow(g, (P - 1)//2, P) == P - 1)
EXP = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]          # monomial order [1, s, u, s^2, s u, u^2]

def mode_coeffs(Lam, a, r0, x0, K):
    bs = [S.run(Lam, a, Fr(1), r0, x0, (Fr(s), Fr(0)), u, K) for s, u in EV]
    def co(kr, ki):
        A = bs[0][kr]; nr, nc = A.nrows(), A.ncols(); out = [dict() for _ in range(6)]
        for i in range(nr):
            for j in range(nc):
                vr = [Fr(int(b[kr][i, j].p), int(b[kr][i, j].q)) for b in bs]
                vi = [Fr(int(b[ki][i, j].p), int(b[ki][i, j].q)) for b in bs]
                if not any(vr) and not any(vi): continue
                for k in range(6):
                    cr = sum(VI[k][q]*vr[q] for q in range(6)); ci = sum(VI[k][q]*vi[q] for q in range(6))
                    if cr or ci: out[k][(i, j)] = gauss_mod(cr, ci)
        return out, nr, nc
    C, nU, nE = co('CR', 'CI'); G, nX, _ = co('GR', 'GI')
    b0 = bs[0]
    return dict(C=C, G=G, nU=nU, nE=nE, nX=nX, nk=b0['nk'], comps=b0['comps'], monsK=b0['monsK'], mons=b0['mons'])

def assemble(D, K):
    mons = D['mons']; deg = lambda i: sum(mons[i])
    monsK = D['monsK']; nk = D['nk']; ncomp = len(D['comps'])
    keep = [i for i in range(len(mons)) if deg(i) <= K - 2]; nblk = D['nE']//len(keep)
    Xmons = [i for i in range(len(mons)) if deg(i) <= K + 1]; nXf = D['nX']//len(Xmons)
    # mode index -> (degree)
    udeg = [deg(monsK[j]) for c in range(ncomp) for j in range(nk)]
    edeg = [deg(keep[j]) for b in range(nblk) for j in range(len(keep))]
    xdeg = [deg(Xmons[j]) for f in range(nXf) for j in range(len(Xmons))]
    tp = [(j, k) for j in range(K + 2) for k in range(K + 2 - j)]
    U = [(j, k, i) for (j, k) in tp for i in range(D['nU']) if j + k + udeg[i] <= K]
    Eq = [(j, k, e) for (j, k) in tp for e in range(D['nE']) if j + k + edeg[e] <= K - 2]
    X = [(j, k, i) for (j, k) in tp for i in range(D['nX']) if j + k + xdeg[i] <= K + 1]
    Ui = {u: n for n, u in enumerate(U)}; Ei = {e: n for n, e in enumerate(Eq)}
    def fac(a, b, j, k):
        if a > j or b > k: return 0
        return pow(IU, a, P)*pow(P - IU, b, P) % P*(factorial(j)//factorial(j - a)) % P*(factorial(k)//factorial(k - b)) % P
    Cm = np.zeros((len(U), len(Eq)), dtype=np.int64)
    for ku, (a, b) in enumerate(EXP):
        for (i, e), c in D['C'][ku].items():
            for (j, k) in tp:
                f = fac(a, b, j, k)
                if not f: continue
                u = (j, k, i); q = (j - a, k - b, e)
                if u in Ui and q in Ei: Cm[Ui[u], Ei[q]] = (Cm[Ui[u], Ei[q]] + c*f) % P
    Xi = {x: n for n, x in enumerate(X)}
    Gm = np.zeros((len(X), len(U)), dtype=np.int64)
    for ku, (a, b) in enumerate(EXP):
        for (x, i), c in D['G'][ku].items():
            for (j, k) in tp:
                f = fac(a, b, j, k)
                if not f: continue
                xv = (j, k, x); u = (j - a, k - b, i)
                if xv in Xi and u in Ui: Gm[Xi[xv], Ui[u]] = (Gm[Xi[xv], Ui[u]] + c*f) % P
    tot_u = [j + k + udeg[i] for (j, k, i) in U]; tot_x = [j + k + xdeg[i] for (j, k, i) in X]
    return Cm, Gm, tot_u, tot_x

def frank(M):
    if M.size == 0: return 0
    return flint.nmod_mat(M.shape[0], M.shape[1], [int(v) for v in M.flat], P).rank()

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    D = mode_coeffs(Fr(1, 20), Fr(3, 5), Fr(3), Fr(3, 5), K); print('mode coefficient matrices', round(time.time() - t0), 's', flush=True)
    Cm, Gm, tu, tx = assemble(D, K); print('4D assembled: unknowns', Cm.shape[0], 'equations', Cm.shape[1], 'gauge', Gm.shape[0], round(time.time() - t0), 's', flush=True)
    # sanity: gauge images are solutions (G C = 0 mod P)
    GC = (Gm.astype(object) @ Cm.astype(object)) % P
    print('gauge in solutions (G*C == 0 mod P):', not GC.any(), flush=True)
    rC = frank(Cm)
    for J in range(K - 1):
        low = [n for n, d in enumerate(tu) if d <= J]; high = [n for n, d in enumerate(tu) if d > J]
        xh = [n for n, d in enumerate(tx) if d >= J + 2]
        rX = frank(Gm[np.ix_(xh, high)]); rGL = frank(Gm[:, low])
        print(f'J={J}: d_J <= {len(low) - rC + len(high) - rX - rGL}   [{round(time.time() - t0)} s]', flush=True)

def direct_counts(Cm, Gm, tu, K):
    """Heuristic (mod-P nullspace) counts dim pi_J(S) - rank pi_J(G); equals the exact count unless P is unlucky."""
    n = Cm.shape[0]
    A = flint.nmod_mat(Cm.shape[1], n, [int(v) for v in Cm.T.flat], P)   # A v = 0  <=>  v C = 0
    nsp, nul = A.nullspace()
    Sb = np.array([[int(nsp[i, j]) for j in range(nul)] for i in range(n)], dtype=np.int64)   # n x nul
    out = {}
    for J in range(K - 1):
        low = [q for q, d in enumerate(tu) if d <= J]
        out[J] = frank(Sb[low, :]) - frank(Gm[:, low])
    return out
