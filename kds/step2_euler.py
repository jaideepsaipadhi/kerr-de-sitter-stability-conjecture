"""
Step 2: Euler integral transforms of Heun's equation, derived directly.

Heun operator (singular points 0, 1, t, infinity):
  H_w = w(w-1)(w-t) d^2 + [g(w-1)(w-t) + d w(w-t) + e w(w-1)] d + (AB w - q),
with g + d + e = A + B + 1 (Fuchs relation).

We look for a Heun operator Ht in z, with the same singular points, satisfying
  Ht_z (z-w)^mu = H_w^dagger (z-w)^mu        (formal adjoint in w).
Then ytilde(z) = Int_C (z-w)^mu y(w) dw solves Ht ytilde = 0 whenever H y = 0,
provided C is a closed (Pochhammer) contour, or the boundary terms vanish.

The script solves for every admissible (mu, gt, dt, et, ABt, qt).
"""
import sympy as sp

z, w, t = sp.symbols('z w t')
g, d, e, A, B, q = sp.symbols('gamma delta epsilon A B q')
mu, gt, dt, et, ABt, qt = sp.symbols('mu gt dt et ABt qt')

P2 = lambda x: x*(x - 1)*(x - t)
P1 = lambda x, g_, d_, e_: g_*(x - 1)*(x - t) + d_*x*(x - t) + e_*x*(x - 1)

# Represent c(z,w) (z-w)^(mu-j) as {j: c}.  Differentiation in w and z by hand.
def dw(f):
    out = {}
    for j, c in f.items():
        out[j] = out.get(j, 0) + sp.diff(c, w)
        out[j + 1] = out.get(j + 1, 0) + c*(mu - j)   # d/dw (z-w)^(mu-j) = -(mu-j)(z-w)^(mu-j-1)... sign below
    return out

def dw(f):
    out = {}
    for j, c in f.items():
        out[j] = out.get(j, 0) + sp.diff(c, w)
        out[j + 1] = out.get(j + 1, 0) - c*(mu - j)
    return out

def dz(f):
    out = {}
    for j, c in f.items():
        out[j] = out.get(j, 0) + sp.diff(c, z)
        out[j + 1] = out.get(j + 1, 0) + c*(mu - j)
    return out

def mul(c0, f):
    return {j: c0*c for j, c in f.items()}

def add(*fs):
    out = {}
    for f in fs:
        for j, c in f.items():
            out[j] = out.get(j, 0) + c
    return out

K = {0: sp.Integer(1)}
adj = add(dw(dw(mul(P2(w), K))), mul(-1, dw(mul(P1(w, g, d, e), K))), mul(A*B*w - q, K))
lhs = add(mul(P2(z), dz(dz(K))), mul(P1(z, gt, dt, et), dz(K)), mul(ABt*z - qt, K))
diff_ = add(lhs, mul(-1, adj))
# multiply through by (z-w)^(2-mu): term j -> (z-w)^(2-j)
expr = sp.expand(sum(c*(z - w)**(2 - j) for j, c in diff_.items()))
poly = sp.Poly(expr, z, w)
eqs = [c for c in poly.coeffs()]
eqs.append(g + d + e - A - B - 1)  # Fuchs relation for the original equation

sols = sp.solve(eqs, [mu, gt, dt, et, ABt, qt, B], dict=True)
for k, sol in enumerate(sols):
    print(f'--- solution {k}')
    for kk, v in sol.items():
        print(f'   {kk} = {sp.factor(v)}')
