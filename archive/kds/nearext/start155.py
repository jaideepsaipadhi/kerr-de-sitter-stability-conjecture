import numpy as np, time
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(11); t0=time.time(); R=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=10**rng.uniform(0.02,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; p=2*x0*(AA+BB-3)/(tt-1)
    if not (p>1 and 2*x0**2*(AA/zz+3)>=tt-1): continue
    u=np.linspace(1e-6,1-1e-6,2000); q=1+x0*u; s=zz-q; N=Nf(tt,zz,AA,BB,s)
    i=np.where(N>=0)[0].min(); zh=q[i]; Sh=zz-zh
    mu_=zh*(zh-1)*(tt-zh); mup=(zh-1)*(tt-zh)+zh*(tt-zh)-zh*(zh-1); H=2*mu_+Sh*mup; Lh=AA/zh+BB-3*zh
    R.append((H/(Sh*Sh*Lh),x0,tt,m,frac))
R=np.array(R); print(len(R),'min ratio',R[:,0].min(),'pct',np.percentile(R[:,0],[0,1,50]))
i=R[:,0].argmin(); print('worst',R[i])
print('min x0 in R3',R[:,1].min(), 'R3 pts with x0<0.05:',(R[:,1]<0.05).sum())
j=R[:,1].argmin(); print(R[j])
rr=(R[:,2]-(1+R[:,1]))/R[:,1]; print('min (t-z0)/x0 in R3 phys',rr.min())
