# SdS Teukolsky-Starobinsky composition constant for s=+2 (a=0): check D^4 Delta^2 (D^+)^4 Delta^2 R = C R on solutions.
import sympy as sp, time
t0=time.time()
r,L,sg,lam=sp.symbols('r Lambda sigma lam'); M=1; s=2
D=r**2-2*M*r-L*r**4/3; Dp=sp.diff(D,r); K=r**2*sg
V=(K**2-sp.I*s*K*Dp)/D+4*sp.I*s*sg*r-sp.Rational(2,3)*L*(s+1)*(2*s+1)*r**2-lam
R0,R1=sp.symbols('R0 R1')
# Delta^{-s}(Delta^{s+1} R')' + V R = 0  =>  D R'' + (s+1) Dp R' + V R = 0
R2=-((s+1)*Dp*R1+V*R0)/D
def d(e): return sp.together(sp.diff(e,r)+sp.diff(e,R0)*R1+sp.diff(e,R1)*R2)
def lin(e):  # express as a*R0+b*R1, simplify coefficients
    e=sp.expand(sp.together(e)); n,dn=sp.fraction(sp.together(e))
    a=sp.cancel(sp.diff(n,R0)/dn); b=sp.cancel(sp.diff(n,R1)/dn); return a*R0+b*R1
Dd=lambda e: lin(d(e)+sp.I*K/D*e)    # D^+ = d/dr + iK/Delta
Dm=lambda e: lin(d(e)-sp.I*K/D*e)    # D   = d/dr - iK/Delta
e=lin(D**2*R0)
for k in range(4): e=Dd(e)
e=lin(D**2*e)
for k in range(4): e=Dm(e); print(k, round(time.time()-t0,1), flush=True)
a=sp.factor(sp.cancel(sp.diff(e,R0))); b=sp.factor(sp.cancel(sp.diff(e,R1)))
print('R1 coeff:', b); print('C =', a)
