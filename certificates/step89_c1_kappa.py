# Item 89: c1 identity in kappa variables (Lam=3 by scaling). Runs in ~7 s.
import sympy as sp, time
t0=time.time()
r0,r1,r2,lam,s,k0,k1,k2=sp.symbols('r0 r1 r2 lam s k0 k1 k2')
r3=-(r0+r1+r2); R=[r0,r1,r2,r3]
e2=sum(R[i]*R[j] for i in range(4) for j in range(i+1,4))
Lam=3; al=1+e2
k3=-(k0+k1+k2); kap=[k0,k1,k2,k3]
W=-sum(R[j]*kap[j] for j in range(4))       # W = i Xi omega
r=sp.Symbol('r')
D1=sp.diff(-sp.Mul(*[(r-x) for x in R]),r).subs(r,r3)
h=sum(1/(r3-x) for x in R[:3])
g,d,e=[1+2*kap[j]+s for j in range(3)]
A=1+s-2*k3; AB=A*(1+2*s)
zinf=(r1-r3)/(r1-r0); t=zinf*(r2-r0)/(r2-r3)
am2=-k3**2-s*k3+(s+1)/2-(s+1)**2/4
Rf3=4*s*r3*W-sp.Rational(2,3)*Lam*(s+1)*(2*s+1)*r3**2+s*(1-al)-lam
am1=(-4*k3*r3*W-2*s*r3*W)/D1-2*s*h*k3+2*h*(k3**2+s*k3)+Rf3/D1-(s+1)**2*h/2
k=(r0-r3)*zinf
c3=2*zinf*am2-k*am1
S=g+d+e; T=d+t*e
q=AB*(1+t)+T*(1-S/2)-c3
ABt=(A-2)*(2*A-d-e-g)
qt=A**2*t+A**2-A*d*t-A*e-A*g*t-A*g-A*t-A+d*t+e+g*t+g+q
St=(g-A+1)+(d-A+1)+(e-A+1); Tt=(d-A+1)+t*(e-A+1)
n1=ABt*(1+t)-qt+Tt*(1-St/2)+sp.Rational(7,2)*(1+t)
n2=ABt+St/2-St**2/4+sp.Rational(3,4)
C0=3/(Lam*(r1-r0)*(r2-r3))
pred=lam*C0+1+t+3*(1-al)*C0

E=sp.expand(n1-pred)
P=sp.Poly(E,k0,k1,k2,lam,s)
expr=sum(sp.cancel(c)*sp.Mul(*[v**e for v,e in zip((k0,k1,k2,lam,s),mon)]) for mon,c in zip(P.monoms(),P.coeffs()))
k2v=-((r0**2-r3**2)*k0+(r1**2-r3**2)*k1)/(r2**2-r3**2)
F=sp.Poly(sp.expand(expr.subs(k2,k2v)),k0,k1,lam,s)
for mon,c in zip(F.monoms(),F.coeffs()):
    print(mon, sp.factor(sp.cancel(c)))
print(round(time.time()-t0,1))
