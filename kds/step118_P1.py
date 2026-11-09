import sympy as sp
z,t,z0,A,B=sp.symbols('z t z0 A B',positive=True)
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z
F=mu/s**2
Dt=sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z)
nu,dn=sp.fraction(sp.factor(sp.together(Dt)))
print('denominator',dn)
P=sp.Poly(sp.expand(nu),z)
print('degree in z',P.degree())
print('factored numerator:',sp.factor(nu))
S=sp.Symbol('S',positive=True)
Ns=sp.Poly(sp.expand(nu.subs(z,z0-S)),S)
for k,c in sorted(zip([m[0] for m in Ns.monoms()],Ns.coeffs())):
    print('S^%d :'%k, sp.factor(c))
