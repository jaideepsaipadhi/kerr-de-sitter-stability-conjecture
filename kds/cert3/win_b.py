# Full-window build in the variable b = a^2 (lambda = -b/e4 eliminated; monic constraint b^2 - e2 b + e4 = 0).
import sympy as sp, time, pickle
t0=time.time()
r0,r2,b,tau=sp.symbols('r0 r2 b tau')
r1=1; r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
e2=sp.expand(sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4))); e4=sp.expand(r0*r1*r2*r3)
Pi=[sp.expand(sp.Mul(*[rs[j]-rs[k] for k in range(4) if k!=j])) for j in range(4)]
Pa=[Pi[0],-Pi[1],Pi[2],-Pi[3]]          # |Pi_j|, D_j = lam*Pa_j
rho=[x**2+b for x in rs]
y2=sp.cancel((Pa[0]-Pa[1])/(rho[1]*Pa[0]-rho[0]*Pa[1]))
y1=sp.cancel((Pa[3]+Pa[1])/(rho[1]*Pa[3]+rho[3]*Pa[1]))
y=y2+tau*(y1-y2)
U=[(1-rho[0]*y)/Pa[0],(rho[1]*y-1)/Pa[1],(1-rho[2]*y)/Pa[2],(rho[3]*y-1)/Pa[3]]   # = lam*u_j * (1/a factor removed)
F1_=U[1]+U[3]; F2_=U[2]+U[3]
t=sp.cancel((r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)))
R=sp.cancel((F1_+t*F2_)/(F1_+F2_)); z0=sp.cancel(t/R)
print('R,z0',round(time.time()-t0,1),flush=True)
Xi=1-b**2/e4; al=-b**2/e4; lamv=-b/e4; K=1/((r1-r0)*(r2-r3))
# lam*Q2: 4 Xi^2 a^2 u0 u3 R * lam = 4 Xi^2 b (U0/lam)(U3/lam) R lam = 4 Xi^2 b U0 U3 R / lam
LQ2=4*Xi**2*b*U[0]*U[3]*R/lamv + K*(1-2*Xi**2*b*y)
LQ1=K*(1-4*al*Xi-4*Xi**2*b*y)
LQ0=-K*(1+7*al)+lamv*(1+t-3*z0)
ex={'Q2':LQ2,'H2':LQ2+LQ1/2+LQ0/4,'H1':LQ2+LQ1+LQ0+4*K}
g=sp.Poly(b**2-e2*b+e4,b); out={}
for k,e in ex.items():
    n,d=sp.fraction(sp.cancel(sp.together(e)))
    print(k,'together',round(time.time()-t0,1),flush=True)
    rem=sp.Poly(sp.expand(n),b).rem(g); cs=rem.all_coeffs(); cs=[0]*(2-len(cs))+cs
    out[k]=(sp.factor(d),[sp.expand(c) for c in cs])
    print(k,'reduced',round(time.time()-t0,1),flush=True)
pickle.dump((out,e2,e4),open('win_b.pkl','wb'))
