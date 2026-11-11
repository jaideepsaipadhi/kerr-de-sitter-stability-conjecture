import sympy as sp, numpy as np, pickle,sys
t,z0,A,be=sp.symbols('t z0 A beta'); T=pickle.load(open('T134.pkl','rb'))
f=[sp.lambdify((t,z0,A,be),Tk) for Tk in T]
k0=float(sys.argv[1]); c0=float(sys.argv[2])
rng=np.random.default_rng(5); bad=0;n=0
for _ in range(100000):
    tt=1+10**rng.uniform(-2,3); zz=1+(tt-1)*rng.uniform(1e-2,1-1e-2); x0=zz-1
    Bf=k0*tt+c0; P=3+(tt-1)/(2*x0)
    a0=max(0.0,zz*((tt-1)/(2*x0**2)-3), P-Bf)
    Bc=max(Bf,3*zz-a0/zz)
    if any(f[k](tt,zz,a0,Bc-(3*zz-a0/zz))>0 for k in (4,5,6)): bad+=1
print(k0,c0,'bad',bad)
