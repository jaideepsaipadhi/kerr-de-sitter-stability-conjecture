import sympy as sp
z,t,z0,A,B,S,v,be=sp.symbols('z t z0 A B S v beta')
s=z0-z; mu=z*(z-1)*(t-z); L=A/z+B-3*z
F=mu/s**2
Dt=sp.diff(mu*sp.diff(F,z),z)-2*s*sp.diff(L*F,z)
nu,dn=sp.fraction(sp.factor(sp.together(Dt)))
N=sp.expand(nu.subs(z,z0-S)); x0=z0-1
M=sp.Poly(sp.expand(sp.expand(N.subs(S,x0*v/(1+v)))*(1+v)**6 ),v) if False else None
num=sp.together(N.subs(S,x0*v/(1+v))); num=sp.Poly(sp.expand(sp.fraction(num)[0]),v)
print('dn',dn,'deg',num.degree())
Bs=3*z0-A/z0+be
for k in range(num.degree()+1):
    c=sp.expand(num.coeff_monomial(v**k).subs(B,Bs)*z0)
    cA=sp.factor(c.coeff(A)); cb=sp.factor(c.coeff(be)); c0=sp.factor(c.subs({A:0,be:0}))
    print(k,'| A:',cA,'| beta:',cb,'| at 0:',c0)
