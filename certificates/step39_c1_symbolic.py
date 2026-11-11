"""
Step 39: SYMBOLIC proof of the c1 identity
    n1 = lam*C0 + 1 + t + 3(1-alpha)*C0,   n2 = -3,
following the exact chain of step38 but with every quantity a symbol.
Parametrisation: roots r0, r1, r2 free, r3 = -(r0+r1+r2) (no cubic term in Delta),
Lam free; Vieta gives  e2 = -(3/Lam)(1-alpha),  e4 = r0 r1 r2 r3 = -3 a^2/Lam,  e3 = -6M/Lam.
a is kept as a symbol and a^2 is eliminated via Vieta at the end. omega, m, lam, s are free symbols.
"""
import sympy as sp

r, eps = sp.symbols('r eps')
r0, r1, r2, Lam, a, om, m, lam, s = sp.symbols('r0 r1 r2 Lam a omega m lam s')
r3 = -(r0 + r1 + r2)
roots = [r0, r1, r2, r3]
e2 = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1, 4))
e4 = r0*r1*r2*r3
a2 = -Lam*e4/3                       # Vieta
al = 1 + Lam*e2/3                    # from e2 = -(3/Lam)(1-alpha)
Xi = 1 + al

Ds = -(Lam/3)*sp.Mul(*[(r - x) for x in roots])
Dps = sp.diff(Ds, r)
Dp = lambda x: Dps.subs(r, x)
K = lambda x: om*(x**2 + a2) - a*m
kap = [sp.I*Xi*K(x)/Dp(x) for x in roots]
g, d, e = [1 + 2*kap[j] + s for j in range(3)]
A = 1 + s - 2*kap[3]; AB = A*(1 + 2*s)
zinf = (r1 - r3)/(r1 - r0); t = zinf*(r2 - r0)/(r2 - r3)

Ks = om*(r**2 + a2) - a*m
Qs = (Xi**2*Ks**2 - sp.I*s*Xi*Ks*Dps)/Ds + 4*sp.I*s*Xi*om*r \
     - sp.Rational(2, 3)*Lam*(s + 1)*(2*s + 1)*r**2 + s*(1 - al) - lam
pr = (s + 1)*Dps/Ds; qr = Qs/Ds
Qr = qr - sp.diff(pr, r)/2 - pr**2/4

# Laurent coefficients at r3, in closed form (no series): eps = r - r3,
# Delta = D1 eps (1 + h eps + ...), h = sum_{j!=3} 1/(r3 - rj) = Delta''/(2 Delta') at r3.
D1 = Dp(r3)
h = sum(1/(r3 - x) for x in roots[:3])
Nf = Xi**2*Ks**2 - sp.I*s*Xi*Ks*Dps              # q_r = Nf/Delta^2 + Rf/Delta
Rf = 4*sp.I*s*Xi*om*r - sp.Rational(2, 3)*Lam*(s + 1)*(2*s + 1)*r**2 + s*(1 - al) - lam
N0 = Nf.subs(r, r3); N1 = sp.diff(Nf, r).subs(r, r3)
am2 = N0/D1**2 + (s + 1)/2 - (s + 1)**2/4
am1 = (N1 - 2*h*N0)/D1**2 + Rf.subs(r, r3)/D1 - (s + 1)**2*h/2
k = (r0 - r3)*zinf
c3 = 2*zinf*am2 - k*am1
S = g + d + e; T = d + t*e
q = AB*(1 + t) + T*(1 - S/2) - c3
gt, dt, et = g - A + 1, d - A + 1, e - A + 1
ABt = (A - 2)*(2*A - d - e - g)
qt = A**2*t + A**2 - A*d*t - A*e - A*g*t - A*g - A*t - A + d*t + e + g*t + g + q
St = gt + dt + et; Tt = dt + t*et
n1 = ABt*(1 + t) - qt + Tt*(1 - St/2) + sp.Rational(7, 2)*(1 + t)
n2 = ABt + St/2 - St**2/4 + sp.Rational(3, 4)
C0 = 3/(Lam*(r1 - r0)*(r2 - r3))
pred = lam*C0 + 1 + t + 3*(1 - al)*C0


def zero(expr, name):
    ex = sp.together(expr)
    nu, _ = sp.fraction(ex)
    nu = sp.expand(nu).subs(a**2, a2)   # only even powers of a survive after this
    nu = sp.expand(sp.expand(nu).subs(a**2, a2))
    q_, rem = sp.div(sp.Poly(nu, s), sp.Poly(s**2 - 4, s)); print(name, ': remainder mod (s^2-4) =', sp.expand(rem.as_expr()), flush=True)


# numeric sanity before the heavy symbolic step
import random
random.seed(1)
vals = {r0: sp.Rational(13, 10), r1: sp.Rational(29, 10), r2: sp.Rational(71, 10), Lam: sp.Rational(1, 20),
        om: sp.Rational(1, 7), m: 2, lam: sp.Rational(37, 5), s: -2}
vals[a] = sp.sqrt(a2.subs(vals))
print(flush=True); print('numeric check n1-pred:', sp.N((n1 - pred).subs(vals), 30), ' n2+3:', sp.N((n2 + 3).subs(vals), 30))
zero(n2 + 3, 'n2 + 3')
zero(n1 - pred, 'n1 - pred')
