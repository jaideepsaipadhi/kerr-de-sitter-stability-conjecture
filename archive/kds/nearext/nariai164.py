import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(4); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,0); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; ell=AA/zz+BB-3*zz
    R.append(((BB)/(AA/zz), (ell-AA/(2*zz))/tt, (tt-1)/x0, x0, (ell-AA/(2*zz)-2*tt)/(AA/zz)))
R=np.array(R); print(len(R),'t-1 range',np.percentile((R[:,2]*R[:,3]),[0,50,100]))
print('B/(A/z0) min',R[:,0].min(),' (ell-A/2z0)/t min',R[:,1].min(),' (ell-A/2z0-2t)/(A/z0) min',R[:,4].min())
sel=R[:,3]<0.05
print('phys x0<.05: A*x0 range',np.percentile((R[sel,0]*0+1),[0]) )
