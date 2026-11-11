"""Stabilisation depth N' of the Killing integrability conditions on KdS (input to the order of Khavkine's K1).
At a point p, a Killing field is determined by its 1-jet (xi^a, w_ab = nabla_a xi_b, antisymmetric): 10 unknowns.
Conditions L_xi (nabla^j Riem)(p) = 0, j = 0, 1, 2, are linear in (xi, w), using
   (L_xi T)_{a1..an} = xi^e nabla_e T_{a1..an} + sum_i T_{a1..e..an} nabla_{a_i} xi^e,   nabla_a xi^e = w_a^e.
EXACT rational arithmetic: coordinates (t, r, x = cos th, phi), metric rational in (r, x), point r = 3, x = 3/5,
Lam = 1/20, a = 3/5, M = 1.  Prints the kernel dimension after imposing j = 0, ..., k."""
import itertools, time, sympy as sp
t0 = time.time()
t, r, x, ph = sp.symbols('t r x phi'); X = [t, r, x, ph]
Lam, a, M = sp.Rational(1, 20), sp.Rational(3, 5), 1
r0, x0 = sp.Rational(7, 2), sp.Rational(-1, 5)
Dr = (r**2 + a**2)*(1 - Lam*r**2/3) - 2*M*r; Dx = 1 + Lam*a**2*x**2/3; Xi = 1 + Lam*a**2/3
rho2 = r**2 + a**2*x**2; s2 = 1 - x**2
g = sp.zeros(4)
g[0, 0] = -(Dr - Dx*a**2*s2)/(rho2*Xi**2)
g[0, 3] = g[3, 0] = a*s2*(Dr - Dx*(r**2 + a**2))/(rho2*Xi**2)
g[3, 3] = s2*(Dx*(r**2 + a**2)**2 - Dr*a**2*s2)/(rho2*Xi**2)
g[1, 1] = rho2/Dr
g[2, 2] = rho2/(Dx*s2)          # dth^2 = dx^2/(1-x^2)
pt = {r: r0, x: x0}
# Taylor data of g at p up to order 4 in (r, x): we need Gamma to order 3, Riem to order 2.
def tay(e, order=4):
    u, v = sp.symbols('u v')
    ser = sp.series(sp.series(e.subs({r: r0 + u, x: x0 + v}), u, 0, order + 1).removeO(), v, 0, order + 1).removeO()
    return sp.expand(ser), u, v
u, v = sp.symbols('u v')
G = sp.zeros(4)
for i in range(4):
    for j in range(i, 4):
        if g[i, j] != 0:
            G[i, j] = G[j, i] = tay(g[i, j])[0]
def trunc(e, n=4):
    p = sp.Poly(sp.expand(e), u, v); return sum(c*u**k*v**l for (k, l), c in p.terms() if k + l <= n)
var = {0: None, 1: u, 2: v, 3: None}
d = lambda e, i: 0 if var[i] is None else sp.diff(e, var[i])
print('metric jets done', round(time.time() - t0), 's', flush=True)
# inverse metric as truncated series: g = g0 + E, ginv = sum (-g0^{-1} E)^k g0^{-1}
g0 = G.subs({u: 0, v: 0}); g0i = g0.inv(); E = G - g0
Gi = g0i; term = g0i
for k in range(1, 5):
    term = (-g0i*E*term).applyfunc(lambda e: trunc(e)); Gi = Gi + term
Gi = Gi.applyfunc(lambda e: trunc(e))
print('inverse done', round(time.time() - t0), 's', flush=True)
Gam = [[[trunc(sum(Gi[a_, l]*(d(G[l, b], c) + d(G[l, c], b) - d(G[b, c], l)) for l in range(4))/2, 3)
         for c in range(4)] for b in range(4)] for a_ in range(4)]
print('Christoffels done', round(time.time() - t0), 's', flush=True)
# Riemann R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}  (order 2)
Rup = {}
for a_, b, c, dd in itertools.product(range(4), repeat=4):
    if c < dd:
        e = d(Gam[a_][dd][b], c) - d(Gam[a_][c][b], dd) + sum(Gam[a_][c][e_]*Gam[e_][dd][b] - Gam[a_][dd][e_]*Gam[e_][c][b] for e_ in range(4))
        e = trunc(e, 2); Rup[(a_, b, c, dd)] = e; Rup[(a_, b, dd, c)] = -e
    elif c == dd:
        Rup[(a_, b, c, dd)] = 0
Rl = {(a_, b, c, dd): trunc(sum(G[a_, e_]*Rup[(e_, b, c, dd)] for e_ in range(4)), 2) for a_, b, c, dd in itertools.product(range(4), repeat=4)}
print('Riemann done', round(time.time() - t0), 's', flush=True)
def cov(T, n, order):   # covariant derivative of a covariant n-tensor (dict of truncated series), new index first
    out = {}
    for idx in itertools.product(range(4), repeat=n + 1):
        e_, rest = idx[0], idx[1:]
        val = d(T[rest], e_)
        for i in range(n):
            for f in range(4):
                val -= Gam[f][e_][rest[i]]*T[rest[:i] + (f,) + rest[i + 1:]]
        out[idx] = trunc(val, order)
    return out
DR = cov(Rl, 4, 1); print('nabla R done', round(time.time() - t0), 's', flush=True)
DDR = cov(DR, 5, 0); print('nabla^2 R done', round(time.time() - t0), 's', flush=True)
at0 = lambda e: sp.Rational(e.subs({u: 0, v: 0})) if hasattr(e, 'subs') else sp.Rational(e)
gp = g0; gpi = g0i
# unknowns: xi^0..3 (index 0..3), w_ab for a<b (index 4..9); w_a^e = w_af g^{fe}
pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
def wlow(a_, b):
    if a_ == b: return [0]*10
    vec = [0]*10; i, j = min(a_, b), max(a_, b); vec[4 + pairs.index((i, j))] = 1 if a_ < b else -1; return vec
def wup(a_, e_):   # coefficient vector of w_a^e
    vec = [0]*10
    for f in range(4):
        cf = gpi[f, e_]
        if cf != 0:
            wl = wlow(a_, f); vec = [vec[k] + cf*wl[k] for k in range(10)]
    return vec
def lie_rows(T0, DT0, n):
    rows = []
    for idx in itertools.product(range(4), repeat=n):
        vec = [0]*10
        for e_ in range(4): vec[e_] += DT0[(e_,) + idx]
        for i in range(n):
            for e_ in range(4):
                c = T0[idx[:i] + (e_,) + idx[i + 1:]]
                if c != 0:
                    wv = wup(idx[i], e_); vec = [vec[k] + c*wv[k] for k in range(10)]
        rows.append(vec)
    return rows
R0 = {k: at0(e) for k, e in Rl.items()}; DR0 = {k: at0(e) for k, e in DR.items()}; DDR0 = {k: at0(e) for k, e in DDR.items()}
# nabla^3 R at p is not available at this truncation; j = 0, 1 use (R, nabla R) and (nabla R, nabla^2 R)
rows = []
for j, (T0, DT0, n) in enumerate([(R0, DR0, 4), (DR0, DDR0, 5)]):
    rows += lie_rows(T0, DT0, n)
    rk = sp.Matrix(rows).rank()
    print(f'after L_xi nabla^{j} Riem = 0: rank {rk}, kernel dim {10 - rk}  [{round(time.time() - t0)} s]', flush=True)
EOF_MARK = None
# --- sanity: d_t and d_phi must lie in the kernel ---
def kill_vec(k):
    xi = [0]*4; xi[k] = 1
    vec = [sp.Integer(0)]*10; vec[k] = 1
    for (i, j) in pairs:   # w_ij = nabla_i xi_j = (d_i g_jk - d_j g_ik)/2 for a coordinate Killing field
        vec[4 + pairs.index((i, j))] = at0((d(G[j, k], i) - d(G[i, k], j))/2)
    return vec
M0 = sp.Matrix(rows)
for k, nm in ((0, 'd_t'), (3, 'd_phi')):
    res = M0*sp.Matrix(kill_vec(k)); print(nm, 'max |residual| =', max(abs(e) for e in res))
