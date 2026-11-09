import sympy as sp, time, pickle
t0=time.time()
r0,r2,lam,tau=sp.symbols('r0 r2 lam tau'); V=(r0,r2,lam,tau)
P_=lambda e: sp.Poly(e,*V,domain='QQ')
class Q:
    def __init__(s,n,d=None): s.n=n if isinstance(n,sp.Poly) else P_(n); s.d=d if d is not None else P_(1)
    def __add__(s,o): o=o if isinstance(o,Q) else Q(o); return Q(s.n*o.d+o.n*s.d, s.d*o.d) if s.d!=o.d else Q(s.n+o.n,s.d)
    __radd__=__add__
    def __neg__(s): return Q(-s.n,s.d)
    def __sub__(s,o): return s+(-(o if isinstance(o,Q) else Q(o)))
    def __rsub__(s,o): return Q(o)-s
    def __mul__(s,o): o=o if isinstance(o,Q) else Q(o); return Q(s.n*o.n,s.d*o.d)
    __rmul__=__mul__
    def inv(s): return Q(s.d,s.n)
    def __truediv__(s,o): o=o if isinstance(o,Q) else Q(o); return s*o.inv()
    def __rtruediv__(s,o): return Q(o)*s.inv()
r1=1; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
e2=sp.expand(sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4))); e4=sp.expand(r0*r1*r2*r3)
a2=Q(-lam*e4); al=Q(lam)*a2; Xi=1+al
Pi=[sp.expand(sp.Mul(*[rs[j]-rs[k] for k in range(4) if k!=j])) for j in range(4)]
Dabs=[Q(lam*Pi[0]),Q(-lam*Pi[1]),Q(lam*Pi[2]),Q(-lam*Pi[3])]
rho=[Q(sp.expand(x**2))+a2 for x in rs]; P=[rho[j]/Dabs[j] for j in range(4)]
y2=(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); y1=(P[1]/rho[1]+P[3]/rho[3])/(P[1]+P[3])
y=y2+Q(tau)*(y1-y2)
u=[-P[0]*(y-1/rho[0]), P[1]*(y-1/rho[1]), -P[2]*(y-1/rho[2]), P[3]*(y-1/rho[3])]
f1=u[1]+u[3]; f2=u[2]+u[3]
t=Q(sp.expand((r1-r3)*(r2-r0)),sp.Poly(sp.expand((r1-r0)*(r2-r3)),*V)); C0=1/Q(sp.expand(lam*(r1-r0)*(r2-r3)))
R=(f1+t*f2)/(f1+f2); z0=t/R; aw=a2*y
Q2=4*Xi*Xi*u[0]*u[3]*a2*R + C0*(1-2*Xi*Xi*aw)
print('Q2 built',round(time.time()-t0,1),'deg n',Q2.n.degree(r0),Q2.n.degree(r2),Q2.n.degree(lam),Q2.n.degree(tau),'terms',len(Q2.n.terms()),flush=True)
