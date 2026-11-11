import numpy as np, time
from scipy.optimize import linprog
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
rng=np.random.default_rng(7); t0=time.time(); rows=[]
while time.time()-t0<30:
    r0=1-10**rng.uniform(-3,0)*0.999; r2=10**rng.uniform(0.05,3); frac=rng.uniform(0.001,0.99999); m=int(rng.integers(1,4))
    try: tt,zz,AA,BB=pdata(r0,1.0,r2,m,frac)
    except Exception: continue
    x0=zz-1; zs=1+0.62*x0*x0/zz
    Ns=Nf(tt,zz,AA,BB,zz-zs); N1=Nf(tt,zz,AA,BB,x0); G2=(tt-1)-2*x0**2*(AA/zz+3)
    sc=abs(Ns)+abs(N1)
    rows.append((Ns/sc,N1/sc,G2*(tt-1)*x0**2*zz**2/sc, G2*tt*x0*zz**3/sc, N1*x0/sc, N1/zz/sc))
R=np.array(rows); print(len(R))
# variables: k1,k2a,k2b,k1b,k1c >=0 ; constraint Ns - k1*N1 - k1b*N1x0 - k1c*N1/z0 - k2a*G2a - k2b*G2b <= -1e-3
A_ub=np.column_stack([-R[:,1],-R[:,2],-R[:,3],-R[:,4],-R[:,5]]); b_ub=-1e-4-R[:,0]
res=linprog(np.ones(5),A_ub=A_ub,b_ub=b_ub,bounds=[(0,None)]*5,method='highs')
print(res.status,res.message,res.x)
