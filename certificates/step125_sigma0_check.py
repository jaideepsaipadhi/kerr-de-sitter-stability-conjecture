# C5: sigma = 0, m = 0 identities (exact).  (i) if (R'/D)' + V R/D^2 = 0 then g = D^2 R'''' solves the same equation;
# (ii) the linear system for cubic solutions has determinant [lam^2 + 6(1-al)lam + 8 - 4al + 8al^2]^2.
import sympy as sp, itertools
r,L,a,M,lam,al=sp.symbols('r Lambda a M lam alpha')
D=(r**2+a**2)*(1-L*r**2/3)-2*M*r
V=-2*L*r**2-2*(1-L*a**2/3)-lam
R0,R1=sp.symbols('R0 R1'); R2=(sp.diff(D,r)*R1-V*R0)/D
d=lambda e: sp.diff(e,r)+sp.diff(e,R0)*R1+sp.diff(e,R1)*R2
e=R0
for k in range(4): e=sp.together(d(e))
g=sp.together(D**2*e); gp=d(g); gpp=d(gp)
print('(i) residual zero:', sp.simplify(sp.numer(sp.together(D*gpp-sp.diff(D,r)*gp+V*g)))==0)
Da=sp.expand(D.subs(L,3*al/a**2)); Va=V.subs(L,3*al/a**2)
c=sp.symbols('c0:4'); Rc=sum(c[i]*r**i for i in range(4))
eqs=sp.Poly(sp.expand(Da*sp.diff(Rc,r,2)-sp.diff(Da,r)*sp.diff(Rc,r)+Va*Rc),r).all_coeffs()
Mx=sp.Matrix([[sp.diff(q,ci) for ci in c] for q in eqs])
print('(ii) det:', sp.factor(Mx.det()))
