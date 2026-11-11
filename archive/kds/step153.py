import numpy as np
rng=np.random.default_rng(14); n=3000000
x0=10**rng.uniform(-5,5,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,9,n)
AA=10**rng.uniform(-6,10,n)*rng.integers(0,2,n); gam=10**rng.uniform(-6,9,n)*rng.integers(0,2,n)
ell=AA/(2*zz)+2*tt+gam
u=rng.uniform(1e-8,1-1e-8,n); q=1+x0*u; s=zz-q
L=ell+s*(AA/(q*zz)+3); Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
g=s*s*L*L+mu*(L-s*Lp)
print('bad',(g<=0).sum(),'of',n)
b=g<=0
if b.any(): print(np.percentile(x0[b],[0,50,100]),np.percentile((tt/zz)[b],[0,50,100]),np.percentile(u[b],[0,50,100]))
