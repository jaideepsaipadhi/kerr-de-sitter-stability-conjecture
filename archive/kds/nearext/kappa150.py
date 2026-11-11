import numpy as np, time
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(7); t0=time.time(); pos=[]; neg=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.99999); m=int(rng.integers(1,4))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; zs=1+0.62*x0*x0/zz
    Ns=Nf(tt,zz,AA,BB,zz-zs); N1=Nf(tt,zz,AA,BB,x0)
    r=Ns/N1
    (pos if N1>0 else neg).append((r,x0,frac))
pos=np.array(pos); neg=np.array(neg)
print('N1>0:',len(pos),'max ratio',pos[:,0].max() if len(pos) else None,' N1<0:',len(neg),'min ratio',neg[:,0].min())
