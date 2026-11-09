import numpy as np
exec(open('/tmp/win.py').read().split("worst=[]")[0])
def data(r0,r1,r2,m,frac,s=-2):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    L,a,al,w1=kds(r0,r1,r2); Xi=1+al; a2=a*a
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Om=[a/(x*x+a2) for x in rs]; Dabs=[abs(Dp(x)) for x in rs]
    w2=a*(1/Dabs[1]-1/Dabs[0])/(P[1]-P[0])
    w=w2+frac*(w1-w2); om=m*w
    kap=[1j*Xi*(om*(x*x+a2)-a*m)/Dp(x) for x in rs]
    t=(r1-r3)*(r2-r0)/((r1-r0)*(r2-r3))
    X=[abs(kap[j]+kap[3]) for j in (0,1,2)]
    A=t*(X[0]**2-(X[1]+X[2])**2)
    C0=3/(L*(r1-r0)*(r2-r3)); c1=1+t+3*(1-al)*C0
    lam=min(lam_min(m,s,max(al,1e-9),a*om) for s in (-2,2))
    B=lam*C0+c1
    e1=2*(kap[1]+kap[3]); et=2*(kap[2]+kap[3])
    q1v=(t-1)*e1.imag/2; qtv=-t*(t-1)*et.imag/2
    z0=1+(t-1)*q1v/(q1v-qtv)
    Lz=lambda z: A/z+B-3*z
    return dict(q1=-q1v/(z0-1) if True else 0,t=t,A=A,B=B,floor=5.5*t-1.25,lam=lam,z0=z0,Lt=Lz(t),L1=Lz(1),C0=C0)
for r0 in (0.9,0.99):
  for frac in (0.5,0.9,0.99):
    print(r0,frac,{k:round(float(v),3) for k,v in data(r0,1.0,10.0,2,frac).items()})
