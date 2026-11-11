import numpy as np
rng=np.random.default_rng(3); tot=bad=0; bz=[]
for i in range(200000):
    tt=1+10**rng.uniform(-2,3.5); zz=1+(tt-1)*10**rng.uniform(-3,-0.0001); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3)); AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=AA/(2*zz)+2*zz+10**rng.uniform(-4,2); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    q=1+x0*np.linspace(1e-6,1-1e-6,600); s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    g=s*s*L*L+mu*L-mu*s*Lp; tot+=1
    if g.min()<=0: bad+=1; bz.append(((q[g.argmin()]-1)/x0, x0/(tt-1)))
print(tot,bad, np.percentile([b[0] for b in bz],[0,50,100]) if bz else '')
