# Endpoint w = w2 inequalities for L(z0) > 0 (item 129). Variables r0 in (0,1), r1 = 1, r2 > 1, lam = Lambda/3.
# Constraint g(lam) = lam^2 e4 + lam e2 + 1 = 0 (unique positive root).  a^2 = -lam e4, alpha = lam a^2.
import sympy as sp, time, pickle
t0=time.time()
r0,r2,lam=sp.symbols('r0 r2 lam',positive=True); r1=sp.Integer(1); r3=-(r0+r1+r2)
rs=[r0,r1,r2,r3]
e2=sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); e4=r0*r1*r2*r3
a2=-lam*e4; al=lam*a2; Xi=1+al
Pi=[sp.Mul(*[rs[j]-rs[k] for k in range(4) if k!=j]) for j in range(4)]
Dabs=[lam*sp.Abs(1)*(-Pi[0]), lam*Pi[1], lam*(-Pi[2]), lam*Pi[3]]   # signs of Delta'(r_j) = -lam*Pi_j: (-,+,-,+) -> |.|
# check: Delta' (r_j) = -lam Pi_j ; sign(-Pi_j)=(-,+,-,+) so |Delta'_j| = (+Pi0? ) handled: D_j = lam*|Pi_j|
Dabs=[lam*Pi[0], -lam*Pi[1], lam*Pi[2], -lam*Pi[3]]
rho=[x**2+a2 for x in rs]; P=[rho[j]/Dabs[j] for j in range(4)]
# work with y = w/a ; Omega_j = a/rho_j  -> w - Omega_j = a (y - 1/rho_j)
y2=(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0])
u0=-P[0]*(y2-1/rho[0]); u3=P[3]*(y2-1/rho[3])    # times a
t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); C0=1/(lam*(r1-r0)*(r2-r3))
# at w2: z0 = t, R = 1.  A/z0 = 4 Xi^2 m^2 u0 u3 a^2 ;  a w = a^2 y2
aw=a2*y2
Q2=4*Xi**2*u0*u3*a2 + C0*(1-2*Xi**2*aw)
Q1=C0*(1-4*al*Xi-4*Xi**2*aw)
Q0=-C0*(1+7*al)+1+t-3*t
exprs={'Q2':Q2,'H2':Q2+Q1/2+Q0/4,'H1':Q2+Q1+Q0+4*C0}
out={}
for k,e in exprs.items():
    n,d=sp.fraction(sp.together(e)); n=sp.expand(n); d=sp.factor(d)
    print(k,'deg num', sp.Poly(n,r0,r2,lam).total_degree(),'den',d, round(time.time()-t0,1),flush=True)
    out[k]=(n,d)
pickle.dump(out,open('w2_raw.pkl','wb'))
