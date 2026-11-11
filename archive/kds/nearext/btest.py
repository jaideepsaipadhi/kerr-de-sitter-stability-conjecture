import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
rng=np.random.default_rng(2); t0=time.time(); worst=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.999); m=int(rng.integers(1,8))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    worst.append(((d['B']-d['t'])/d['t'],r0,r2,frac,m,d['lam']))
worst.sort(); print(len(worst)); [print([round(float(x),4) for x in w]) for w in worst[:6]]
