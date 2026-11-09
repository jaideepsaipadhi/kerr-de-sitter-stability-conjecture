# Full window: w = a*y, y = y2 + tau (y1 - y2), tau in [0,1].  Same three inequalities Q2>0, H2>0, H1>0.
import sympy as sp, time, pickle
t0=time.time()
r0,r2,lam,tau=sp.symbols('r0 r2 lam tau',positive=True); r1=sp.Integer(1); r3=-(r0+r1+r2)
rs=[r0,r1,r2,r3]
e2=sp.expand(sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4))); e4=sp.expand(r0*r1*r2*r3)
a2=-lam*e4; al=lam*a2; Xi=1+al
Pi=[sp.Mul(*[rs[j]-rs[k] for k in range(4) if k!=j]) for j in range(4)]
Dabs=[lam*Pi[0], -lam*Pi[1], lam*Pi[2], -lam*Pi[3]]
rho=[x**2+a2 for x in rs]; P=[rho[j]/Dabs[j] for j in range(4)]
y2=(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); y1=(P[1]/rho[1]+P[3]/rho[3])/(P[1]+P[3])
y=y2+tau*(y1-y2)
u=[-P[0]*(y-1/rho[0]), P[1]*(y-1/rho[1]), -P[2]*(y-1/rho[2]), P[3]*(y-1/rho[3])]   # times a
f1=u[1]+u[3]; f2=u[2]+u[3]
t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); C0=1/(lam*(r1-r0)*(r2-r3))
R=(f1+t*f2)/(f1+f2); z0=t/R
aw=a2*y
Q2=4*Xi**2*u[0]*u[3]*a2*R + C0*(1-2*Xi**2*aw)
Q1=C0*(1-4*al*Xi-4*Xi**2*aw)
Q0=-C0*(1+7*al)+1+t-3*z0
g=sp.Poly(lam**2*e4+lam*e2+1,lam); D=sp.expand(e2**2-4*e4)
out={}
for k,e in {'Q2':Q2,'H2':Q2+Q1/2+Q0/4,'H1':Q2+Q1+Q0+4*C0}.items():
    n,d=sp.fraction(sp.together(e)); n=sp.expand(n)
    dl=sp.Poly(n,lam).degree()
    pr=sp.prem(sp.Poly(n,lam),g); cs=sp.Poly(pr,lam).all_coeffs(); cs=([0]*(2-len(cs)))+cs
    F1,F0=[sp.expand(c) for c in cs]
    out[k]=(sp.factor(d),dl,F1,F0)
    print(k,'deg_lam',dl,'F1 tdeg',sp.Poly(F1,r0,r2,tau).total_degree(),round(time.time()-t0,1),flush=True)
pickle.dump((out,e2,e4),open('win_red.pkl','wb'))
