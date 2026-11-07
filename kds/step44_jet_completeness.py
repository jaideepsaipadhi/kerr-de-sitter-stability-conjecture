"""
Step 44: jet-level test of the Lambda-analogue of AHW Prop 6.5 (B1):
  linearized Lambda-vacuum h with dpsi_0 = dpsi_4 = 0  =>  h = delta^* X + (finite-dim family)   locally.
At a point p of the KdS static region, for polynomial h of degree <= K (10 components):
  S_K = { h : jets of (dRic - Lam h), dpsi_0, dpsi_4 vanish to degree K-2 at p },
  G_K = { delta^* X truncated, X of degree <= K+1 }  (subset of S_K: a built-in consistency check).
Report  d_J = rank pi_J(S_K) - rank pi_J(G_K)  (jets that extend, modulo gauge), for J <= K-2.
Expected for Kerr (Lam = 0): stabilises at 4 (PD tangent: m, a, NUT, acceleration).
Machinery: truncated power series (TPS) in y = x - p, degree <= N = K+1; dual numbers for the linearisation,
batched over all basis columns of h.
"""
import itertools, sys
import numpy as np

def setup(N):
    mons = [m for d in range(N + 1) for m in itertools.product(range(d + 1), repeat=4) if sum(m) == d]
    idx = {m: i for i, m in enumerate(mons)}
    M = len(mons)
    pairs = [(i, j, idx[tuple(a + b for a, b in zip(mi, mj))]) for i, mi in enumerate(mons) for j, mj in enumerate(mons)
             if sum(mi) + sum(mj) <= N]
    P = np.array(pairs)
    D = []
    for mu in range(4):
        Dm = np.zeros((M, M))
        for i, m in enumerate(mons):
            if m[mu] > 0:
                mm = list(m); mm[mu] -= 1; Dm[i, idx[tuple(mm)]] = m[mu]
        D.append(Dm)        # coefficient vector c (row) -> c @ Dm is the derivative
    return mons, idx, M, P, D

class T:
    """TPS: bg (M,) complex/real; pert (ncol, M) or None."""
    def __init__(s, bg, pert=None): s.bg, s.pert = bg, pert

def build(N):
    mons, idx, M, P, D = setup(N)
    def mulb(f, g):                       # background x background
        out = np.zeros(M, dtype=np.result_type(f, g)); np.add.at(out, P[:, 2], f[P[:, 0]]*g[P[:, 1]]); return out
    def Mat(f):                           # matrix of multiplication by f acting on row vectors
        Mf = np.zeros((M, M), dtype=f.dtype); np.add.at(Mf, (P[:, 0], P[:, 2]), f[P[:, 1]]); return Mf
    def mul(a, b):
        bg = mulb(a.bg, b.bg); pert = None
        if a.pert is not None: pert = a.pert @ Mat(b.bg)
        if b.pert is not None:
            t = b.pert @ Mat(a.bg); pert = t if pert is None else pert + t
        return T(bg, pert)
    def add(a, b, cb=1.0):
        pert = a.pert if b.pert is None else (b.pert*cb if a.pert is None else a.pert + cb*b.pert)
        return T(a.bg + cb*b.bg, pert)
    def dif(a, mu):
        return T(a.bg @ D[mu], None if a.pert is None else a.pert @ D[mu])
    def const(c):
        v = np.zeros(M); v[0] = c; return v
    def var(mu, c0):
        v = np.zeros(M); v[0] = c0; v[idx[tuple(int(k == mu) for k in range(4))]] = 1.0; return v
    def recip(f):
        c = f[0]; u = f.copy(); u[0] = 0; out = const(1/c); term = const(1/c)
        for _ in range(N):
            term = -mulb(term, u)/c; out = out + term
        return out
    def sincos(f):
        c = f[0]; u = f.copy(); u[0] = 0
        s_u = np.zeros(M); c_u = const(1.0); term = const(1.0)
        for k in range(1, N + 1):
            term = mulb(term, u)/k
            if k % 2 == 1: s_u = s_u + (-1)**((k - 1)//2)*term
            else: c_u = c_u + (-1)**(k//2)*term
        return np.sin(c)*c_u + np.cos(c)*s_u, np.cos(c)*c_u - np.sin(c)*s_u
    return dict(mons=mons, idx=idx, M=M, D=D, mulb=mulb, Mat=Mat, mul=mul, add=add, dif=dif, const=const,
                var=var, recip=recip, sincos=sincos)

def kds_metric(E, mulb, recip, Lam, a, Mm, r0, th0):
    M = E['M']
    r = E['var'](1, r0); th = E['var'](2, th0)
    s, c = E['sincos'](th)
    al = Lam*a*a/3; Xi = 1 + al
    rr = mulb(r, r); cc = mulb(c, c); ss = mulb(s, s)
    rho2 = rr + a*a*cc
    Dr = mulb(rr + a*a*E['const'](1), E['const'](1) - Lam/3*rr) - 2*Mm*r
    Dth = E['const'](1) + al*cc
    irho2 = recip(rho2)
    A1 = -mulb(Dr, irho2)/Xi**2; A2 = mulb(mulb(Dth, ss), irho2)/Xi**2
    rpa = rr + a*a*E['const'](1)
    g = [[np.zeros(M) for _ in range(4)] for _ in range(4)]
    g[0][0] = A1 + a*a*A2
    g[0][3] = g[3][0] = -a*mulb(A1, ss) - a*mulb(A2, rpa)
    g[3][3] = a*a*mulb(A1, mulb(ss, ss)) + mulb(A2, mulb(rpa, rpa))
    g[1][1] = mulb(rho2, recip(Dr)); g[2][2] = mulb(rho2, recip(Dth))
    z = np.zeros(M)
    lv = [Xi*rpa, Dr, z, Xi*a*E['const'](1)]
    nv = [Xi*rpa, -Dr, z, Xi*a*E['const'](1)]
    mv = [1j*a*Xi*s, z.astype(complex), Dth.astype(complex), 1j*Xi*recip(s)]
    return g, lv, nv, mv


def run(Lam, a, Mm, K, r0, th0, metric_fn=None, extra=None):
    N = K + 1; E = build(N); M = E['M']; mulb, recip = E['mulb'], E['recip']
    if metric_fn is not None:
        g, lv, nv, mv = metric_fn(E, mulb, recip)
    else:
        g, lv, nv, mv = kds_metric(E, mulb, recip, Lam, a, Mm, r0, th0)
    # inverse metric: block (t,phi) 2x2 and diagonal r, theta
    det = mulb(g[0][0], g[3][3]) - mulb(g[0][3], g[0][3]); idet = recip(det)
    gi = [[np.zeros(M) for _ in range(4)] for _ in range(4)]
    gi[0][0] = mulb(g[3][3], idet); gi[3][3] = mulb(g[0][0], idet); gi[0][3] = gi[3][0] = -mulb(g[0][3], idet)
    gi[1][1] = recip(g[1][1]); gi[2][2] = recip(g[2][2])
    # unknowns: h_{ab}, a<=b, monomials of degree <= K
    monsK = [i for i, m in enumerate(E['mons']) if sum(m) <= K]
    comps = [(p, q) for p in range(4) for q in range(p, 4)]
    ncol = len(comps)*len(monsK)
    def hcomp(p, q):
        if p > q: p, q = q, p
        k = comps.index((p, q)); P = np.zeros((ncol, M))
        for j, i in enumerate(monsK): P[k*len(monsK) + j, i] = 1.0
        return P
    mul, add, dif = E['mul'], E['add'], E['dif']
    G = [[T(g[p][q], hcomp(p, q)) for q in range(4)] for p in range(4)]
    # delta g^{-1} = - g^{-1} h g^{-1}
    GI = [[None]*4 for _ in range(4)]
    for p in range(4):
        for q in range(4):
            pert = np.zeros((ncol, M))
            for u in range(4):
                for v in range(4):
                    if not (gi[p][u].any() and gi[v][q].any()): continue
                    pert -= hcomp(u, v) @ E['Mat'](mulb(gi[p][u], gi[v][q]))
            GI[p][q] = T(gi[p][q], pert)
    dG = [[[dif(G[p][q], mu) for mu in range(4)] for q in range(4)] for p in range(4)]
    Gam = [[[None]*4 for _ in range(4)] for _ in range(4)]
    for A_ in range(4):
        for b in range(4):
            for c_ in range(b, 4):
                acc = None
                for d in range(4):
                    if not GI[A_][d].bg.any() and not GI[A_][d].pert.any(): continue
                    t = add(add(dG[d][c_][b], dG[d][b][c_]), dG[b][c_][d], -1.0)
                    t = mul(GI[A_][d], t); acc = t if acc is None else add(acc, t)
                acc = T(0.5*acc.bg, 0.5*acc.pert)
                Gam[A_][b][c_] = Gam[A_][c_][b] = acc
    z = np.zeros(M)
    mb = [np.conj(v) for v in mv]
    low = lambda v: [sum(mulb(g[p][q], v[q]) for q in range(4)) for p in range(4)]
    ll, nl = low(lv), low(nv)
    chk = max(np.abs(sum(mulb(ll[p], lv[p]) for p in range(4))).max(), np.abs(sum(mulb(ll[p], mv[p]) for p in range(4))).max(),
              np.abs(sum(mulb(low(mv)[p], mv[p]) for p in range(4))).max())
    # Riemann R^e_{bcd} = d_c Gam^e_{db} - d_d Gam^e_{cb} + Gam^e_{cf} Gam^f_{db} - Gam^e_{df} Gam^f_{cb}
    Ric = [[None]*4 for _ in range(4)]
    psi0 = np.zeros((ncol, M), dtype=complex); psi4 = np.zeros((ncol, M), dtype=complex)
    bpsi0 = np.zeros(M, dtype=complex)
    for e in range(4):
        for b in range(4):
            for c_ in range(4):
                for d in range(c_ + 1, 4):
                    R = add(dif(Gam[e][d][b], c_), dif(Gam[e][c_][b], d), -1.0)
                    for f in range(4):
                        R = add(R, mul(Gam[e][c_][f], Gam[f][d][b]))
                        R = add(R, mul(Gam[e][d][f], Gam[f][c_][b]), -1.0)
                    for (cc_, dd_, sg) in ((c_, d, 1.0), (d, c_, -1.0)):
                        if e == cc_:
                            Ric[b][dd_] = T(sg*R.bg, sg*R.pert) if Ric[b][dd_] is None else add(Ric[b][dd_], R, sg)
                        # delta(R_{abcd}) l^a m^b l^c m^d with R_{abcd} = g_{ae} R^e_{bcd}
                        for (vec_low, A2v, B2v, out) in ((ll, mv, lv, 'p0'), (nl, mb, nv, 'p4')):
                            w = mulb(mulb(vec_low[e], A2v[b]), mulb(B2v[cc_], A2v[dd_]))*sg
                            if not w.any(): continue
                            contrib = R.pert @ E['Mat'](w)
                            # + h_{ae} R^e_{bcd} vec^a ... : h_{ae} vec^a
                            vecu = lv if out == 'p0' else nv
                            for a_ in range(4):
                                ww = mulb(mulb(vecu[a_], A2v[b]), mulb(B2v[cc_], A2v[dd_]))*sg
                                if ww.any():
                                    contrib = contrib + hcomp(a_, e) @ E['Mat'](mulb(R.bg.astype(complex), ww))
                            if out == 'p0': psi0 += contrib; bpsi0 += mulb(R.bg, w)
                            else: psi4 += contrib
    # constraints: jets of degree <= K-2
    keep = [i for i, m in enumerate(E['mons']) if sum(m) <= K - 2]
    rows = []
    for (p, q) in comps:
        Ep = Ric[p][q].pert - Lam*hcomp(p, q)
        rows.append(Ep[:, keep])
    for X in (psi0, psi4):
        rows.append(X[:, keep].real); rows.append(X[:, keep].imag)
    Cmat = np.hstack(rows)          # (ncol, nconstraints): h-coefficient vector c satisfies c @ Cmat = 0
    # background checks
    bric = max(np.abs(Ric[p][q].bg - Lam*g[p][q])[keep].max() for (p, q) in comps)
    # gauge: (delta^* X)_{pq} = d_p X_q + d_q X_p - 2 Gam^f_{pq} X_f, X lower, degree <= K+1 = N
    nX = 4*M
    Gh = np.zeros((nX, ncol))
    for f in range(4):
        for i in range(M):
            Xv = np.zeros(M); Xv[i] = 1.0
            for k, (p, q) in enumerate(comps):
                v = np.zeros(M)
                if p == f: v += Xv @ E['D'][q]
                if q == f: v += Xv @ E['D'][p]
                for ff in range(4):
                    if ff == f: v -= 2*mulb(Gam[ff][p][q].bg, Xv)
                Gh[f*M + i, k*len(monsK):(k + 1)*len(monsK)] = v[monsK]
    gauge_resid = np.abs(Gh @ Cmat).max()/max(1, np.abs(Gh).max()*np.abs(Cmat).max())
    extras = []
    if extra is not None:
        for hx in extra(E, mulb, recip):
            v = np.zeros(ncol)
            for k, (p, q) in enumerate(comps):
                v[k*len(monsK):(k + 1)*len(monsK)] = hx[p][q][monsK]
            extras.append(v)
    return dict(extras=np.array(extras) if extras else None, E=E, Cmat=Cmat, Gh=Gh, monsK=monsK, comps=comps, chk=chk, bric=bric, bpsi0=np.abs(bpsi0[keep]).max(),
                gauge_resid=gauge_resid)

def counts(res, K, tol=1e-9):
    Cmat, Gh, monsK, mons = res['Cmat'], res['Gh'], res['monsK'], res['E']['mons']
    U, sv, Vt = np.linalg.svd(Cmat.T, full_matrices=True)      # constraints x unknowns
    rank = (sv > tol*sv[0]).sum()
    S = Vt[rank:].T                                            # null space basis (unknowns x dimS)
    out = {}
    nk = len(monsK)
    for J in range(0, K - 1):
        sel = [k*nk + j for k in range(len(res['comps'])) for j, i in enumerate(monsK) if sum(mons[i]) <= J]
        def rk(Mx):
            if Mx.size == 0: return 0
            s_ = np.linalg.svd(Mx, compute_uv=False); return int((s_ > tol*max(s_[0], 1e-300)).sum())
        out[J] = rk(S[sel, :]) - rk(Gh[:, sel].T)
    return out, sv

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for (Lam, a) in [(0.0, 0.6), (0.05, 0.6)]:
        res = run(Lam, a, 1.0, K, 3.0, 1.0)
        print(f'Lam={Lam} a={a} K={K}: tetrad null/orth check {res["chk"]:.1e}, background Ric-Lam g {res["bric"]:.1e}, '
              f'background psi0 {res["bpsi0"]:.1e}, gauge in S_K residual {res["gauge_resid"]:.1e}', flush=True)
        out, sv = counts(res, K)
        print('   d_J (solution jets mod gauge, degree <= J):', out, flush=True)
