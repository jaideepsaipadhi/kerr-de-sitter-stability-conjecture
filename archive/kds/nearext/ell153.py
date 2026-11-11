import numpy as np, time
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(9); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,80))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    ell=AA/zz+BB-3*zz; R.append((ell/tt, (ell-AA/(2*zz))/tt, (ell-AA/(2*zz)-2*zz)/(tt-zz+1e-300), BB/tt))
R=np.array(R); print(len(R)); print('min ell/t',R[:,0].min(),' min (ell-A/2z0)/t',R[:,1].min(),' min (ell-A/2z0-2z0)/(t-z0)',R[:,2].min(),'min B/t',R[:,3].min())
