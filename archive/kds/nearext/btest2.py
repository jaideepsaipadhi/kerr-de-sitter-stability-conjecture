import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
rng=np.random.default_rng(3); t0=time.time(); R=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.999); m=int(rng.integers(1,30))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=d['z0']-1
    if x0<5: continue
    R.append((d['B']/(d['A']/d['z0']),(d['B']-d['t'])/d['t'],x0,d['t']/x0))
R=np.array(R); print(len(R),'phys x0>=5: B/(A/z0)',np.percentile(R[:,0],[0,50,100]),'(B-t)/t',np.percentile(R[:,1],[0,50,100]),'t/x0',np.percentile(R[:,3],[0,50,100]))
