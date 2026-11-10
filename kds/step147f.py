import numpy as np
rng=np.random.default_rng(10); n=1000000
x0=10**rng.uniform(-4,4,n); zz=1+x0
AA=10**rng.uniform(-5,8,n); gam=10**rng.uniform(-6,5,n)*rng.integers(0,2,n)
u=rng.uniform(1e-7,1-1e-7,n); q=1+x0*u; s=zz-q
L=AA*(1/q-1/(2*zz))+5*zz-3*q+gam; Lp=AA/q**2+3
L1=AA*(1-1/(2*zz))+5*zz-3+gam
for th in (1.3,1.5,1.7):
  for nm,tt in (('t=z0',zz),('tmax',1+np.minimum(2*x0*L1,2*x0**2*(AA/zz+3)))):
    mu=q*(q-1)*(tt-q)
    g=th*(2-th)*s*s*L*L+th*mu*L-(2-th)*mu*s*Lp
    print(th,nm,(g<=0).sum())
