exec(open('step141.py').read().split("P=pts(")[0])
import numpy as np
rng=np.random.default_rng(99); out=[]
while len(out)<3000:
    tt=10**rng.uniform(2,3.5); zz=1+(tt-1)*10**rng.uniform(-3,-1.3); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=AA/(2*zz)+2*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    c=np.array([g(tt,zz,AA,BB) for g in cf],float); r=np.roots(c)
    r=[x.real for x in r if abs(x.imag)<1e-9*abs(x) and 0<x.real<x0]
    if len(r)!=1: continue
    out.append((tt,zz,AA,BB,r[0]))
for c in (0.5,1.0,1.5):
    kf=lambda x,c=c: x**2*(1-c*x*(1-x)); dkf=lambda x,c=c: 2*x*(1-c*x*(1-x))+x**2*(-c*(1-2*x))
    print(c,sum(cond(Q,kf,dkf) for Q in out),'/',len(out),flush=True)
