import numpy as np, time, warnings, pickle, sympy as sp; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
n,d,C=pickle.load(open('/home/claude/kds/P182.pkl','rb'))
z,t,z0,A,B,Q2,K=sp.symbols('z t z0 A B Q2 Khat',real=True)
fE=sp.lambdify((z,t,z0,A,Q2,K),n); fC=sp.lambdify((z,t,z0,Q2,K),C)
LO=[np.inf]*3; HI=[-np.inf]*3; rng=np.random.default_rng(8); t0=time.time(); st={'n':0,'Lt<0':0,'E_bad_gate':0,'E_bad_nogate':0,'C_bad_gate':0,'C_bad_nogate':0}
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    Lt=AA/tt+BB-3*tt; Kv=1.05*(-Lt)*tt*(tt-1)/(2*(tt-zz)**2)
    zg=zz+(tt-zz)*np.linspace(1e-6,1-1e-6,800)
    Ev=fE(zg,tt,zz,AA,q1*q1,Kv); Cv=fC(zg,tt,zz,q1*q1,Kv)
    st['n']+=1; g=Lt<0; st['Lt<0']+=g
    for k,W in enumerate((q1**4*tt**2, q1**4*tt*(tt-zz)*zz, q1**4*tt**3/zz)):
        r=Ev/(Lt*W)
        if g: LO[k]=min(LO[k],(-Ev/(-Lt*W)).min())
        else:
            if (Ev>0).any(): HI[k]=max(HI[k],r[Ev>0].max())
print('need c >= ',HI,' and c <= ',LO)
