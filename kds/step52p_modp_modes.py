"""
Step 52: EXACT jet counts for B1, mode-reduced (LEMMAS 58-60).
h = E * H(r, x),  E = exp(-i sigma t + i m phi),  x = cos(theta),  coordinates (t, r, x, phi).
Jets of H in the 2 variables (y_r, y_x) = (r - r0, x - x0) at a rational point; all background quantities are rational in (r, x).
Exact arithmetic: background TPS coefficients in Fraction; linear maps as pairs (Re, Im) of flint.fmpq_mat; complex ranks by
realification (rank_C M = rank_Q [[Re, Im], [-Im, Re]] / 2).
Outputs d_J = dim pi_J(S_K) - dim pi_J(G)  (an UPPER bound for the true d_J, by item 58(a)), exactly.
Expected: sigma = 0, m = 0 -> 4 (the PD tangents are stationary and axisymmetric); sigma != 0 -> 0.
"""
import itertools, sys, time
from fractions import Fraction as Fr
import flint, os
PMOD = int(os.environ.get('PRIME', 998244353))
def QM(n, d): return int(n) % PMOD * pow(int(d) % PMOD, PMOD - 2, PMOD) % PMOD
CONJ = os.environ.get('CONJ', '1') == '1'   # impose also the conjugate-tetrad conditions (Theorem A applied to hbar)

def setup(N):
    mons = [(i, d - i) for d in range(N + 1) for i in range(d, -1, -1)]
    idx = {m: k for k, m in enumerate(mons)}
    pairs = [(i, j, idx[(mi[0] + mj[0], mi[1] + mj[1])]) for i, mi in enumerate(mons) for j, mj in enumerate(mons) if sum(mi) + sum(mj) <= N]
    return mons, idx, pairs

class Ctx:
    def __init__(s, N):
        s.N = N; s.mons, s.idx, s.pairs = setup(N); s.M = len(s.mons)
        s.D = []
        for v in range(2):
            Dm = flint.nmod_mat(s.M, s.M, PMOD)
            for k, m in enumerate(s.mons):
                if m[v] > 0:
                    mm = list(m); mm[v] -= 1; Dm[k, s.idx[tuple(mm)]] = m[v]
            s.D.append(Dm)
    # ---- background (complex Fraction lists) ----
    def zero(s): return [Fr(0)]*s.M
    def const(s, c): v = s.zero(); v[0] = Fr(c); return v
    def var(s, v, c0): w = s.const(c0); w[s.idx[(1, 0) if v == 0 else (0, 1)]] = Fr(1); return w
    def mul(s, f, g):
        out = s.zero()
        for i, j, k in s.pairs:
            if f[i] and g[j]: out[k] += f[i]*g[j]
        return out
    def add(s, *fs): return [sum(c) for c in zip(*fs)]
    def scal(s, c, f): return [c*a for a in f]
    def recip(s, f):
        c = f[0]; u = list(f); u[0] = Fr(0); out = s.const(1/c); term = s.const(1/c)
        for _ in range(s.N):
            term = s.scal(Fr(-1)/c, s.mul(term, u)); out = s.add(out, term)
        return out
    def bdiff(s, f, mu):  # background derivative; mu: 0=t,1=r,2=x,3=phi
        if mu in (0, 3): return s.zero()
        v = mu - 1; out = s.zero()
        for k, m in enumerate(s.mons):
            if m[v] > 0:
                mm = list(m); mm[v] -= 1; out[s.idx[tuple(mm)]] += m[v]*f[k]
        return out
    def Mat(s, f):        # multiplication by real background f acting on row vectors
        A = flint.nmod_mat(s.M, s.M, PMOD)
        for i, j, k in s.pairs:
            if f[j]: A[i, k] += QM(f[j].numerator, f[j].denominator)
        return A

# complex numbers as (re, im) pairs of Fractions; complex bg as (re_list, im_list)
class CB:
    def __init__(s, re, im=None): s.re = re; s.im = im
class CP:   # complex pert: (Re, Im) fmpq_mat, rows = unknowns
    def __init__(s, re, im): s.re, s.im = re, im
    def __add__(s, o): return CP(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return CP(s.re - o.re, s.im - o.im)
    def scale(s, c):
        cr, ci = [QM(v.numerator, v.denominator) for v in c]
        return CP(s.re*cr - s.im*ci, s.re*ci + s.im*cr)

def run(Lam, a, Mm, r0, x0, sigma, m, K):
    N = K + 1; C = Ctx(N); M = C.M
    r = C.var(0, r0); x = C.var(1, x0)
    al = Fr(Lam)*a*a/3; Xi = 1 + al
    rr = C.mul(r, r); xx = C.mul(x, x); one = C.const(1)
    s2 = C.add(one, C.scal(-1, xx))
    rho2 = C.add(rr, C.scal(a*a, xx))
    rpa = C.add(rr, C.const(a*a))
    Dr = C.add(C.mul(rpa, C.add(one, C.scal(-Fr(Lam)/3, rr))), C.scal(-2*Mm, r))
    Dth = C.add(one, C.scal(al, xx))
    irho2 = C.recip(rho2)
    A1 = C.scal(-1/Xi**2, C.mul(Dr, irho2)); A2 = C.scal(1/Xi**2, C.mul(C.mul(Dth, s2), irho2))
    g = [[C.zero() for _ in range(4)] for _ in range(4)]
    g[0][0] = C.add(A1, C.scal(a*a, A2))
    g[0][3] = g[3][0] = C.add(C.scal(-a, C.mul(A1, s2)), C.scal(-a, C.mul(A2, rpa)))
    g[3][3] = C.add(C.scal(a*a, C.mul(A1, C.mul(s2, s2))), C.mul(A2, C.mul(rpa, rpa)))
    g[1][1] = C.mul(rho2, C.recip(Dr)); g[2][2] = C.mul(rho2, C.recip(C.mul(Dth, s2)))
    det = C.add(C.mul(g[0][0], g[3][3]), C.scal(-1, C.mul(g[0][3], g[0][3]))); idet = C.recip(det)
    gi = [[C.zero() for _ in range(4)] for _ in range(4)]
    gi[0][0] = C.mul(g[3][3], idet); gi[3][3] = C.mul(g[0][0], idet); gi[0][3] = gi[3][0] = C.scal(-1, C.mul(g[0][3], idet))
    gi[1][1] = C.recip(g[1][1]); gi[2][2] = C.recip(g[2][2])
    # background Christoffels (real)
    dg = [[[C.bdiff(g[p][q], mu) for mu in range(4)] for q in range(4)] for p in range(4)]
    Gb = [[[None]*4 for _ in range(4)] for _ in range(4)]
    for A in range(4):
        for b in range(4):
            for c in range(b, 4):
                acc = C.zero()
                for d in range(4):
                    if any(gi[A][d]):
                        acc = C.add(acc, C.mul(gi[A][d], C.add(dg[d][c][b], dg[d][b][c], C.scal(-1, dg[b][c][d]))))
                Gb[A][b][c] = Gb[A][c][b] = C.scal(Fr(1, 2), acc)
    # mode wave vector: d_mu(E F) = E (d_mu F + k_mu F), k = (-i sigma, 0, 0, i m) as complex pairs
    sr, si = sigma
    kvec = [(si, -sr), (Fr(0), Fr(0)), (Fr(0), Fr(0)), (Fr(0), Fr(m))]
    # unknowns
    monsK = [k for k, mo in enumerate(C.mons) if sum(mo) <= K]; nk = len(monsK)
    comps = [(p, q) for p in range(4) for q in range(p, 4)]; ncol = len(comps)*nk
    def hcomp(p, q):
        if p > q: p, q = q, p
        k0 = comps.index((p, q)); R = flint.nmod_mat(ncol, M, PMOD)
        for j, i in enumerate(monsK): R[k0*nk + j, i] = 1
        return CP(R, flint.nmod_mat(ncol, M, PMOD))
    Z = lambda: CP(flint.nmod_mat(ncol, M, PMOD), flint.nmod_mat(ncol, M, PMOD))
    Mats = {}
    def MatC(f):
        key = id(f)
        if key not in Mats: Mats[key] = (C.Mat(f), f)
        return Mats[key][0]
    def pmul(P, f):   # pert * real background
        A = C.Mat(f); return CP(P.re*A, P.im*A)
    def pdiff(P, mu):
        if mu in (1, 2):
            D = C.D[mu - 1]; return CP(P.re*D, P.im*D)
        return P.scale(kvec[mu])
    # dg pert = h ; d(g^-1) = -g^-1 h g^-1
    H = [[hcomp(p, q) for q in range(4)] for p in range(4)]
    GIp = [[None]*4 for _ in range(4)]
    for p in range(4):
        for q in range(4):
            acc = Z()
            for u in range(4):
                for v in range(4):
                    w = C.mul(gi[p][u], gi[v][q])
                    if any(w): acc = acc - pmul(H[u][v], w)
            GIp[p][q] = acc
    dH = [[[pdiff(H[p][q], mu) for mu in range(4)] for q in range(4)] for p in range(4)]
    # Gamma pert: dGam^a_bc = 1/2 [ dg^{ad} (..bg..) + g^{ad} (dh terms) ]
    Gp = [[[None]*4 for _ in range(4)] for _ in range(4)]
    for A in range(4):
        for b in range(4):
            for c in range(b, 4):
                acc = Z()
                for d in range(4):
                    bgt = C.add(dg[d][c][b], dg[d][b][c], C.scal(-1, dg[b][c][d]))
                    if any(bgt): acc = acc + pmul(GIp[A][d], bgt)
                    if any(gi[A][d]):
                        acc = acc + pmul(dH[d][c][b] + dH[d][b][c] - dH[b][c][d], gi[A][d])
                half = (Fr(1, 2), Fr(0))
                Gp[A][b][c] = Gp[A][c][b] = acc.scale(half)
    # tetrad
    z = C.zero()
    lv = [C.scal(Xi, rpa), Dr, z, C.const(Xi*a)]
    nv = [C.scal(Xi, rpa), C.scal(-1, Dr), z, C.const(Xi*a)]
    mre = [z, z, C.scal(-1, C.mul(s2, Dth)), z]; mim = [C.scal(a*Xi, s2), z, z, C.const(Xi)]
    low = lambda v: [C.add(*[C.mul(g[p][q], v[q]) for q in range(4)]) for p in range(4)]
    ll, nl = low(lv), low(nv)
    # Riemann pert, Ricci pert and dpsi0, dpsi4
    keep = [k for k, mo in enumerate(C.mons) if sum(mo) <= K - 2]
    Ric = [[Z() for _ in range(4)] for _ in range(4)]
    psi0 = Z(); psi4 = Z(); psi0b = Z(); psi4b = Z()
    # complex vector components: m = mre + i mim; mbar = mre - i mim
    def cvec(sign): return [(mre[k], C.scal(sign, mim[k])) for k in range(4)]
    mC, mbC = cvec(1), cvec(-1)
    def cmulb(a_, b_):   # complex bg * complex bg
        return (C.add(C.mul(a_[0], b_[0]), C.scal(-1, C.mul(a_[1], b_[1]))), C.add(C.mul(a_[0], b_[1]), C.mul(a_[1], b_[0])))
    def pmulc(P, w):     # pert * complex bg
        A = C.Mat(w[0]); B = C.Mat(w[1])
        return CP(P.re*A - P.im*B, P.re*B + P.im*A)
    for e in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(c + 1, 4):
                    # pert of R^e_{bcd} = d_c G^e_db - d_d G^e_cb + G^e_cf G^f_db - G^e_df G^f_cb (linearised)
                    R = pdiff(Gp[e][d][b], c) - pdiff(Gp[e][c][b], d)
                    Rb = C.add(C.bdiff(Gb[e][d][b], c), C.scal(-1, C.bdiff(Gb[e][c][b], d)))
                    for f in range(4):
                        R = R + pmul(Gp[e][c][f], Gb[f][d][b]) + pmul(Gp[f][d][b], Gb[e][c][f])
                        R = R - pmul(Gp[e][d][f], Gb[f][c][b]) - pmul(Gp[f][c][b], Gb[e][d][f])
                        Rb = C.add(Rb, C.mul(Gb[e][c][f], Gb[f][d][b]), C.scal(-1, C.mul(Gb[e][d][f], Gb[f][c][b])))
                    for (cc, dd, sg) in ((c, d, 1), (d, c, -1)):
                        Rs = R.scale((Fr(sg), Fr(0)))
                        if e == cc: Ric[b][dd] = Ric[b][dd] + Rs
                        for (vl, vu, A2v, out) in ((ll, lv, mC, 0), (nl, nv, mbC, 1), (ll, lv, mbC, 2), (nl, nv, mC, 3)):
                            # w = vl_e * A2v^b * vu^cc * A2v^dd   (complex)
                            w = cmulb(cmulb((vl[e], z), A2v[b]), cmulb((vu[cc], z), A2v[dd]))
                            if not (any(w[0]) or any(w[1])): continue
                            contrib = pmulc(Rs, w)
                            # + h_{ae} vu^a R^e_{b cc dd} * A2v^b vu^cc A2v^dd
                            for a_ in range(4):
                                ww = cmulb(cmulb((C.mul(vu[a_], Rb), z), A2v[b]), cmulb((vu[cc], z), A2v[dd]))
                                if sg < 0: ww = (C.scal(-1, ww[0]), C.scal(-1, ww[1]))
                                if any(ww[0]) or any(ww[1]): contrib = contrib + pmulc(H[a_][e], ww)
                            if out == 0: psi0 = psi0 + contrib
                            elif out == 1: psi4 = psi4 + contrib
                            elif out == 2: psi0b = psi0b + contrib
                            else: psi4b = psi4b + contrib
    # constraint matrix over C: columns = equations (keep degrees), rows = unknowns
    def cols(P):
        Re = flint.nmod_mat(ncol, len(keep), PMOD); Im = flint.nmod_mat(ncol, len(keep), PMOD)
        for i in range(ncol):
            for jj, j in enumerate(keep):
                Re[i, jj] = P.re[i, j]; Im[i, jj] = P.im[i, j]
        return Re, Im
    blocks = []
    for (p, q) in comps:
        blocks.append(cols(Ric[p][q] - hcomp(p, q).scale((Fr(Lam), Fr(0)))))
    blocks.append(cols(psi0)); blocks.append(cols(psi4))
    if CONJ: blocks.append(cols(psi0b)); blocks.append(cols(psi4b))
    neq = sum(B[0].ncols() for B in blocks)
    CR = flint.nmod_mat(ncol, neq, PMOD); CI = flint.nmod_mat(ncol, neq, PMOD); off = 0
    for (Re, Im) in blocks:
        for i in range(ncol):
            for j in range(Re.ncols()):
                CR[i, off + j] = Re[i, j]; CI[i, off + j] = Im[i, j]
        off += Re.ncols()
    # gauge map: X_f mode, degree <= N: (delta^* X)_pq = d_p X_q + d_q X_p - 2 Gb^f_pq X_f
    nX = 4*M
    GR = flint.nmod_mat(nX, ncol, PMOD); GI = flint.nmod_mat(nX, ncol, PMOD)
    for f in range(4):
        for i in range(M):
            Xv = C.zero(); Xv[i] = Fr(1)
            for k0, (p, q) in enumerate(comps):
                vr = C.zero(); vi = C.zero()
                for (aa, bb) in ((p, q), (q, p)):
                    if bb == f:
                        if aa in (1, 2): vr = C.add(vr, C.bdiff(Xv, aa))
                        else:
                            kr, ki = kvec[aa]; vr = C.add(vr, C.scal(kr, Xv)); vi = C.add(vi, C.scal(ki, Xv))
                vr = C.add(vr, C.scal(-2, C.mul(Gb[f][p][q], Xv)))
                for j, mi in enumerate(monsK):
                    GR[f*M + i, k0*nk + j] = QM(vr[mi].numerator, vr[mi].denominator); GI[f*M + i, k0*nk + j] = QM(vi[mi].numerator, vi[mi].denominator)
    if sigma == (Fr(0), Fr(0)) and m == 0:
        # extra stationary gauge: X = t xi, phi xi (xi = d_t, d_phi Killing): h_pq = delta_p^A xi_q + delta_q^A xi_p
        ext = flint.nmod_mat(nX + 4, ncol, PMOD); exi = flint.nmod_mat(nX + 4, ncol, PMOD)
        for i in range(nX):
            for j in range(ncol): ext[i, j] = GR[i, j]; exi[i, j] = GI[i, j]
        row = nX
        for Aa in (0, 3):
            for bb in (0, 3):
                for k0, (p, q) in enumerate(comps):
                    v = C.zero()
                    if p == Aa: v = C.add(v, g[q][bb])
                    if q == Aa: v = C.add(v, g[p][bb])
                    for j, mi in enumerate(monsK): ext[row, k0*nk + j] = QM(v[mi].numerator, v[mi].denominator)
                row += 1
        GR, GI = ext, exi
    return dict(CR=CR, CI=CI, GR=GR, GI=GI, monsK=monsK, mons=C.mons, comps=comps, nk=nk, ncol=ncol)

def realify(Re, Im):   # complex (rows x cols) -> real 2rows x 2cols with rank doubled
    r, c = Re.nrows(), Re.ncols(); R = flint.nmod_mat(2*r, 2*c, PMOD)
    for i in range(r):
        for j in range(c):
            a_, b_ = Re[i, j], Im[i, j]
            if a_: R[i, j] = a_; R[r + i, c + j] = a_
            if b_: R[i, c + j] = b_; R[r + i, j] = -b_
    return R

def crank(Re, Im): return realify(Re, Im).rank()//2

def sub_rows(A, rows):
    B = flint.nmod_mat(len(rows), A.ncols(), PMOD)
    for ii, i in enumerate(rows):
        for j in range(A.ncols()): B[ii, j] = A[i, j]
    return B
def sub_cols(A, cols_):
    B = flint.nmod_mat(A.nrows(), len(cols_), PMOD)
    for i in range(A.nrows()):
        for jj, j in enumerate(cols_): B[i, jj] = A[i, j]
    return B

def counts(res, K):
    CR, CI = res['CR'], res['CI']; nk, comps, monsK, mons = res['nk'], res['comps'], res['monsK'], res['mons']
    rkC = crank(CR, CI); out = {}
    for J in range(0, K - 1):
        low = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) <= J]
        high = [k*nk + j for k in range(len(comps)) for j, i in enumerate(monsK) if sum(mons[i]) > J]
        rkH = crank(sub_rows(CR, high), sub_rows(CI, high)) if high else 0
        dimS = len(low) - rkC + rkH
        dimG = crank(sub_cols(res['GR'], low), sub_cols(res['GI'], low))
        out[J] = dimS - dimG
    return out

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    LAM = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 20)
    cases = eval(sys.argv[2]) if len(sys.argv) > 2 else [((Fr(0), Fr(0)), 0), ((Fr(1, 3), Fr(1, 5)), 1), ((Fr(2, 7), Fr(0)), 2), ((Fr(0), Fr(0)), 1)]
    for sigma, m in cases:
        t0 = time.time()
        res = run(LAM, Fr(3, 5), Fr(1), Fr(3), Fr(3, 5), sigma, m, K)
        print(f'K={K} sigma={sigma[0]}+{sigma[1]}i m={m}: exact d_J upper bounds {counts(res, K)}  ({time.time()-t0:.0f}s)', flush=True)
