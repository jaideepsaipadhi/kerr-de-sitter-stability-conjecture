import numpy as np, time, warnings, pickle, sympy as sp; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
n,d,C=pickle.load(open('/home/claude/kds/P182.pkl','rb'))
z,t,z0,A,B,Q2,K=sp.symbols('z t z0 A B Q2 Khat',real=True)
fE=sp.lambdify((z,t,z0,A,Q2,K),n); fC=sp.lambdify((z,t,z0,Q2,K),C)
TH=(0.01,0.03,0.1,0.3,1.0); OK={th:0 for th in TH}; rng=np.random.default_rng(8); t0=time.time(); st={'n':0,'Lt<0':0,'E_bad_gate':0,'E_bad_nogate':0,'C_bad_gate':0,'C_bad_nogate':0}
while time.time()-t0<40:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=AA/tt+BB-3*tt
    zg=zz+(tt-zz)*np.linspace(1e-6,1-1e-6,800)
    st['n']+=1
    for th in TH:
        Kv=th*q1*q1/(tt-zz); Ltil=Lt+2*Kv*(tt-zz)**2/(tt*(tt-1))
        Ev=fE(zg,tt,zz,AA,q1*q1,Kv); Cv=fC(zg,tt,zz,q1*q1,Kv)
        OK[th]+= (Ev.max()<=0) and (Cv.max()<=0) and Ltil>0
print(st['n'],OK)
