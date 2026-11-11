import sympy as sp, numpy as np, pickle
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
f=[sp.lambdify((t,z0,A,be),Tk) for Tk in T]
rng=np.random.default_rng(4); need=[]
for _ in range(20000):
    tt=1+10**rng.uniform(-4,5); zz=1+(tt-1)*rng.uniform(1e-4,1-1e-4)
    a0=max(0.0,zz*((tt-1)/(2*(zz-1)**2)-3)); sc=a0/zz+3*zz
    lo,hi=0.0,1e3
    ok=lambda c: all(f[k](tt,zz,a0,c*sc)<=0 for k in (4,5,6))
    if ok(0): need.append(0); continue
    if not ok(hi): need.append(np.inf); continue
    for _ in range(50):
        m=(lo+hi)/2; lo,hi=(lo,m) if ok(m) else (m,hi)
    need.append(hi)
need=np.array(need); print(np.percentile(need,[50,90,99,100]))
