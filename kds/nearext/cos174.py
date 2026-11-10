import numpy as np, time, warnings; warnings.filterwarnings('ignore')
exec(open('phys149.py').read().split("rng=np.random.default_rng(1)")[0])
src=open('phys149.py').read()
def pdq(r0,r1,r2,m,frac):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0]); w=w2+frac*(w1-w2); om=m*w
    kap=[1j*Xi*(om*(x*x+a2)-a*m)/Dp(x) for x in rs]
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3)); X=[abs(kap[j]+kap[3]) for j in (0,1,2)]
    A=t*(X[0]**2-(X[1]+X[2])**2); C0=3/(L*(r1-r0)*(r2-r3)); c1=1+t+3*(1-al)*C0
    xi=a*om; l0=max(m,2); lamc=l0*(l0+1)-4-4*al-4*al*m*(1+al)-2*(1+al)**2*abs(xi)*(m+2)
    e1=2*(kap[1]+kap[3]); et=2*(kap[2]+kap[3]); q1v=(t-1)*e1.imag/2; qtv=-t*(t-1)*et.imag/2
    z0=1+(t-1)*q1v/(q1v-qtv)
    return t,z0,A,lamc*C0+c1,-q1v/(z0-1)
rng=np.random.default_rng(1); t0=time.time(); st={'n':0,'Ltneg':0,'q1sign+':0,'ok':0}; fails=[]
while time.time()-t0<45:
    r0=1-10**rng.uniform(-3.5,0)*0.999; r2=1+10**rng.uniform(-4,3.5); frac=rng.uniform(0.0001,0.99999); m=int(rng.integers(1,40))
    try: tt,zz,AA,BB,q1=pdq(r0,1.0,r2,m,frac)
    except Exception: continue
    st['n']+=1; Lt=AA/tt+BB-3*tt; st['q1sign+']+=q1>0
    if Lt>0: st['ok']+=1; continue
    st['Ltneg']+=1
    kap=1.05*(-Lt)*tt*(tt-1)/(2*q1*(tt-zz)**2)
    z=zz+(tt-zz)*np.linspace(1e-6,1-1e-6,2000); mu=z*(z-1)*(tt-z); mup=(z-1)*(tt-z)+z*(tt-z)-z*(z-1)
    Q=q1+kap*(tt-z); qt=(z-zz)*Q; qtp=Q-kap*(z-zz)
    Lp=-AA/z**2-3; corr=kap*(z-zz)**2*(2*q1+kap*(tt-z))/(z*(z-1)); 
    dcorr=np.gradient(corr,z); Ltp=Lp+dcorr
    C=mup*qt-2*mu*qtp; D=qt*Ltp-kap*mup
    good=(C.max()<=0) and (D.max()<=0)
    st['ok']+=good
    if not good: fails.append((C.max()>0,D.max()>0,q1>0))
print(st, fails[:5], 'Cfail',sum(f[0] for f in fails),'Dfail',sum(f[1] for f in fails))
