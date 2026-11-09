import numpy as np, time
exec(open('/tmp/chk.py').read().split("for r0 in (0.9")[0])
def cos_check(d,kap,n=4000):
    t,A,B,z0,q1=d['t'],d['A'],d['B'],d['z0'],abs(d['q1'])
    z=z0+(t-z0)*np.linspace(1e-6,1,n)
    mu=z*(z-1)*(t-z); mup=(z-1)*(t-z)+z*(t-z)-z*(z-1)
    g=q1+kap*(t-z); qt=(z-z0)*g; qtp=g-kap*(z-z0)
    Lt=A/z+B-3*z+kap*(z-z0)**2*(2*q1+kap*(t-z))/(z*(z-1))
    Ltp=np.gradient(Lt,z)
    C=mup*qt-2*mu*qtp; D=qt*Ltp-kap*mup
    return bool(C.max()<=0), bool(D[2:-2].max()<=0), bool(Lt[-1]>0)
rng=np.random.default_rng(3); t0=time.time(); N=0; bad=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-2.5,0); r2=10**rng.uniform(0.1,2.5); frac=rng.uniform(0.02,0.999); m=int(rng.integers(1,6))
    try: d=data(r0,1.0,r2,m,frac)
    except Exception: continue
    t,z0,q1=d['t'],d['z0'],abs(d['q1'])
    Lt=d['Lt']; k0=max(0.0,-Lt*t*(t-1)/(2*q1*(t-z0)**2))
    res=None
    for c in [0,1.2,2,4,8]:
        kap=k0*c if k0>0 else 0.0
        if k0>0 and c==0: continue
        r=cos_check(d,kap)
        if all(r): res=(c,r); break
    N+=1
    if res is None: bad.append((round(r0,4),round(r2,1),round(frac,3),m,round(Lt,2),cos_check(d,1.2*k0)))
print('N',N,'bad',len(bad)); [print(b) for b in bad[:10]]
