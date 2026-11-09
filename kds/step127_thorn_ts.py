# a=0 (SdS), Kinnersley: check (Th')^4 (r^4 R) = c(r) * Delta^2 (D^+)^4 (Delta^2 R) as operators on modes e^{-i sigma t}.
import sympy as sp
r,L,sg=sp.symbols('r Lambda sigma'); M=1
D=r**2-2*M*r-L*r**4/3; R=sp.Function('R')(r)
n_d=lambda f: -(D/(2*r**2))*(sp.diff(f,r)+sp.I*sg*r**2/D*f)
gam=-D/(2*r**3)+sp.diff(D,r)/(4*r**2)
F=r**4*R
e=n_d(F)-4*gam*F; e=n_d(e)-2*gam*e; e=n_d(e); e=n_d(e)+2*gam*e
Dp=lambda f: sp.diff(f,r)+sp.I*sg*r**2/D*f
T=D**2*R
for k in range(4): T=Dp(T)
T=D**2*T
ratio=sp.simplify(sp.expand(e)/sp.expand(T))
print(ratio)
