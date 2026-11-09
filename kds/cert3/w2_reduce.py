import sympy as sp, pickle, time, numpy as np
t0=time.time()
r0,r2,lam=sp.symbols('r0 r2 lam',positive=True); r1=1; r3=-(r0+r1+r2)
rs=[r0,r1,r2,r3]
e2=sp.expand(sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4))); e4=sp.expand(r0*r1*r2*r3)
g=sp.Poly(lam**2*e4+lam*e2+1,lam)
raw=pickle.load(open('w2_raw.pkl','rb')); red={}
for k,(n,d) in raw.items():
    rem=sp.Poly(n,lam).rem(g)   # careful: g not monic; use pseudo-division
    pr=sp.prem(sp.Poly(n,lam),g)  # e4^k * n = q g + pr
    cs=sp.Poly(pr,lam).all_coeffs()
    if len(cs)==1: cs=[0]+cs
    F1,F0=[sp.expand(c) for c in cs]
    red[k]=(F1,F0); print(k,'reduced',round(time.time()-t0,1),flush=True)
pickle.dump((red,e2,e4),open('w2_red.pkl','wb'))
# numeric sign survey
f=lambda e: sp.lambdify((r0,r2),e,'numpy')
D=e2**2-4*e4
for k,(F1,F0) in red.items():
    X=-F1*e2+2*e4*F0
    fF1,fX,fD=f(F1),f(X),f(D)
    R0,R2=np.meshgrid(np.linspace(0.02,0.98,60),1/np.linspace(0.02,0.98,60))
    s=np.sqrt(fD(R0,R2)); val=fX(R0,R2)-fF1(R0,R2)*s
    print(k,'F1 sign range',np.sign(fF1(R0,R2)).min(),np.sign(fF1(R0,R2)).max(),' X sign',np.sign(fX(R0,R2)).min(),np.sign(fX(R0,R2)).max(),' X-F1 s >0 all:',bool((val>0).all()))
