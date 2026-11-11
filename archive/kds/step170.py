import numpy as np
rng=np.random.default_rng(15); n=300000; bad=0; badz0=0
x0=10**rng.uniform(-4,4,n); zz=1+x0; tt=zz+x0*10**rng.uniform(-6,7,n)
AA=10**rng.uniform(-6,9,n)*rng.integers(0,2,n); gam=10**rng.uniform(-6,8,n)*rng.integers(0,2,n)
ell=AA/(2*zz)+2*tt+gam
u=np.linspace(1e-6,1-1e-6,400)
B=0
for i in range(n):
    q=1+x0[i]*u; s=zz[i]-q; mu=q*(q-1)*(tt[i]-q); mup=(q-1)*(tt[i]-q)+q*(tt[i]-q)-q*(q-1)
    wC=(2*mu+s*mup)/s; wmin=wC.min()
    L=ell[i]+s*(AA[i]/(q*zz[i])+3); Lp=AA[i]/q**2+3
    Q=wmin*L-2*mu*Lp
    if Q.min()<=0: bad+=1
    if i>40000: break
print('tested',i,'bad',bad)
