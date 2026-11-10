import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('cos176.py').read().split("rng=np.random.default_rng(2)")[0])
rng=np.random.default_rng(5); t0=time.time(); pts=[]
while time.time()-t0<25 and len(pts)<300:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    if AA/tt+BB-3*tt<0: pts.append((tt,zz,AA,BB,q1))
print('bend pts',len(pts))
def test(P,f,fp,fpp):
    tt,zz,AA,BB,q1=P; z=zz+(tt-zz)*np.linspace(1e-7,1,3000)
    mu=z*(z-1)*(tt-z); mup=(z-1)*(tt-z)+z*(tt-z)-z*(z-1); mupp=-6*z+2*(tt+1)
    L=AA/z+BB-3*z; Lp=-AA/z**2-3; q=q1*(z-zz)
    F,Fp,Fpp=f(z),fp(z),fpp(z)
    C=mup*F*q-2*mu*F*q1-mu*Fp*q; D=q*(Fp*L+F*Lp)+0.5*q1*(mup*Fp+mu*Fpp)
    return C.max()<=0 and D.max()<=0 and F[-1]<0 and fp(np.array([zz]))[0]<=0
res={}
for p in (1,2,3,4,6):
  for dl in (0.01,0.1,0.5):
    ok=0
    for P in pts:
        tt,zz=P[0],P[1]; c=(1+dl)/(tt-zz)**p
        ok+=test(P,lambda z:1-c*(z-zz)**p,lambda z:-c*p*(z-zz)**(p-1),lambda z:-c*p*(p-1)*(z-zz)**max(p-2,0) if p>=2 else 0*z)
    res[(p,dl)]=ok
print(res)
