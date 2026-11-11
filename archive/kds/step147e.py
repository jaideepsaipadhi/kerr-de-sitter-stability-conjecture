import numpy as np
rng=np.random.default_rng(9); bad=0; n=400000
tt=1+10**rng.uniform(-3,4,n); zz=1+(tt-1)*10**rng.uniform(-4,0,n)*0.99999; x0=zz-1
AA=10**rng.uniform(-4,7,n); gam=10**rng.uniform(-5,4,n)*rng.integers(0,2,n)
u=rng.uniform(1e-6,1-1e-6,n); q=1+x0*u; s=zz-q
L=AA*(1/q-1/(2*zz))+5*zz-3*q+gam; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
for th in (1.3,1.5,1.7):
    g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
    print(th,(g<=0).sum(),'of',n)
