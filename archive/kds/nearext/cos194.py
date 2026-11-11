import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(21); t0=time.time(); G=[]; n=0
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(1e-5,0.99999); m=int(rng.integers(1,120))
    try: t,z0,A,B,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    n+=1; Lt=A/t+B-3*t
    if Lt<0: G.append(((z0-1)/(t-1), m, t))
G=np.array(G); print("all"); print(n,len(G),'max zf on gate',G[:,0].max())
