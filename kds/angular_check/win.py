import numpy as np, sys
sys.argv=['x']
exec(open('/tmp/ang.py').read().split("print('check'")[0].replace("for m in ():","for m in []:"))
def kds(r0,r1,r2):
    r3=-(r0+r1+r2); rs=[r0,r1,r2,r3]
    e2=sum(rs[i]*rs[j] for i in range(4) for j in range(i+1,4)); e4=r0*r1*r2*r3
    # Lam from Vieta: e2=-(3/Lam)(1-al), e4=-3a^2/Lam, al=Lam a^2/3 -> al = -Lam^2 e4/9 ; e2 = -(3/Lam)(1+Lam^2 e4/9)
    # solve for Lam: e2*Lam = -3 - Lam^2 e4/3  -> (e4/3)Lam^2 + e2 Lam + 3 = 0
    A_,B_,C_=e4/3,e2,3.0; disc=B_*B_-4*A_*C_
    L=(-B_-np.sqrt(disc))/(2*A_); 
    if L<=0: L=(-B_+np.sqrt(disc))/(2*A_)
    a2=-L*e4/3; a=np.sqrt(a2); al=L*a2/3
    D=lambda r: -(L/3)*np.prod([r-x for x in rs])
    Dp=lambda r: -(L/3)*sum(np.prod([r-x for k,x in enumerate(rs) if k!=j]) for j in range(4))
    P=[(x*x+a2)/abs(Dp(x)) for x in rs]; Om=[a/(x*x+a2) for x in rs]
    w1=(P[1]*Om[1]+P[3]*Om[3])/(P[1]+P[3])
    return L,a,al,w1
worst=[]
for r0 in [0.3,0.6,0.85,0.95,0.99]:
  for r2 in [3,5,10,30,100]:
    L,a,al,w1=kds(r0,1.0,r2)
    for m in (1,2,3):
        xi=a*m*w1
        lm=min(lam_min(m,s,max(al,1e-9),xi) for s in (-2,2))
        worst.append((lm,m,r0,r2,round(al,4),round(xi,3)))
worst.sort(); print(worst[:8])
