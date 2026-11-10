import numpy as np, sys
CC,DD=0.5,1.0
rng=np.random.default_rng(5); tot=bad=0; badroot=[]
for i in range(200000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3))
    AA=a0+10**rng.uniform(-3,4)*rng.integers(0,2)*(a0+1)
    bb=CC*AA/zz+DD*zz+10**rng.uniform(-4,1); BB=3*zz-AA/zz+bb
    if 2*x0*(AA+BB-3)<=tt-1: continue
    Sh=x0*rng.uniform(0.001,0.999); zh=zz-Sh; mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1)
    w0=(2*mu_+Sh*mup)/Sh
    sg=np.linspace(1e-6,1,400); q=zz-Sh*sg; s=Sh*sg; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    g=-2*mu*w0/Sh**2+w0**2*s**3/Sh**4-2*L*w0*s**2/Sh**2+2*mu*Lp
    tot+=1; bad+= g.max()>0
print(tot,bad)
