import sympy as sp, numpy as np, pickle
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
f=sp.lambdify((t,z0,A,be),T)
rng=np.random.default_rng(8); viol={j:0 for j in range(1,6)}; n=0
for _ in range(200000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1); need=(tt-1)/(2*x0)-3*x0-AA*x0/zz
    bb=max(need,0)+10**rng.uniform(-3,3)*(1+abs(need))
    c=np.array(f(tt,zz,AA,bb),float); sc=np.abs(c).max(); n+=1
    for j in range(1,6):
        if c[j]<-1e-9*sc and c[j+1]>1e-9*sc: viol[j]+=1
print(n,viol)
import mpmath as mp
mp.mp.dps=40
rng=np.random.default_rng(8); out={}
for _ in range(200000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1); need=(tt-1)/(2*x0)-3*x0-AA*x0/zz
    bb=max(need,0)+10**rng.uniform(-3,3)*(1+abs(need))
    c=np.array(f(tt,zz,AA,bb),float); sc=np.abs(c).max()
    if any(c[j]<-1e-9*sc and c[j+1]>1e-9*sc for j in range(1,6)):
        cs=[mp.mpf(x) for x in f(mp.mpf(tt),mp.mpf(zz),mp.mpf(AA),mp.mpf(bb))][::-1]
        r=[x for x in mp.polyroots(cs,maxsteps=300,extraprec=300) if abs(mp.im(x))<1e-25 and mp.re(x)>0]
        out[len(r)]=out.get(len(r),0)+1
print(out)
