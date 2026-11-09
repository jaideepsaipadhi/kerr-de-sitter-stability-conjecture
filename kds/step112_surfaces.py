import numpy as np
exec(open('/tmp/mb.py').read().split("stats=")[0])
rng=np.random.default_rng(12)
res={'bulk':{'lo':0,'hi':0,'int':0},'A0':{'lo':0,'hi':0,'int':0}}
for _ in range(15000):
    t=1+10**rng.uniform(-2,5); z0=1+rng.uniform(0,1)*(t-1); z=1+rng.uniform(0,1)*(z0-1)
    x0=z0-1; s=z0-z; mu0=z0*x0*(t-z0); fl=5.5*t-1.25
    # (iii) A=0, vary B >= floor
    Bs=fl*10**np.linspace(0,4,200); f,ok=F(t,z0,z,np.zeros_like(Bs),Bs)
    if ok.sum()>=3:
        fo=np.where(ok,f,-np.inf); i=np.argmax(fo); idx=np.where(ok)[0]
        res['A0']['lo' if i==idx[0] else 'hi' if i==idx[-1] else 'int']+=1
    # (ii) bulk boundary: L0 = 1.25 mu0/s^2 ; vary B, A = z0*(L0 - B + 3 z0)
    L0=1.25*mu0/s**2; Bs=fl*10**np.linspace(0,4,200); As=z0*(L0-Bs+3*z0)
    f,ok=F(t,z0,z,As*(1+1e-12),Bs); ok=ok&(As>=0)
    if ok.sum()>=3:
        fo=np.where(ok,f,-np.inf); i=np.argmax(fo); idx=np.where(ok)[0]
        res['bulk']['lo' if i==idx[0] else 'hi' if i==idx[-1] else 'int']+=1
print(res)
