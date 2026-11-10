import numpy as np, time
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(9); t0=time.time(); n=0; bad=0; worst=9
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,80))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1
    q=1+x0*np.linspace(1e-5,1-1e-5,300); s=zz-q; L=AA/q+BB-3*q; Lp=AA/q**2+3; mu=q*(q-1)*(tt-q)
    g=s*s*L*L+mu*(L-s*Lp); n+=1
    r=(g/(s*s*L*L+mu*L+mu*s*Lp)).min(); worst=min(worst,r); bad+=r<=0
print(n,'bad',bad,'worst normalised',worst)
