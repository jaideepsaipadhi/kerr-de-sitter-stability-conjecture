import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
rng=np.random.default_rng(12); t0=time.time(); rows=[]
while time.time()-t0<40:
    r0=1-10**rng.uniform(-2.5,0); r2=10**rng.uniform(0.1,2.5); frac=rng.uniform(0.0,1.0); m=int(rng.integers(1,6))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    t,z0,A,B=d['t'],d['z0'],d['A'],d['B']; Lz0=A/z0+B-3*z0
    rows.append((Lz0/(abs(A)/z0+abs(B)+3*z0),r0,r2,frac,m,d['lam'],(z0-1)/(t-1)))
rows.sort(); [print([round(float(v),4) for v in r]) for r in rows[:8]]; print(len(rows))
